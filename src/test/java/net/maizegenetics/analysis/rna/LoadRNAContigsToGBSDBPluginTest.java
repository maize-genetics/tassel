package net.maizegenetics.analysis.rna;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.tag.*;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.Utils;

import org.junit.Test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

import static org.junit.Assert.*;

/**
 * Tests {@link LoadRNAContigsToGBSDBPlugin} and {@link RNADeMultiplexProductionPlugin} on the
 * small-RNA reads in {@code dataFiles/RNA}.
 *
 * <p>The contigs are the most common small-RNA inserts in that FASTQ (the read up to the
 * Illumina adapter), chosen so that no two contigs share a word in either orientation.
 * Every read whose insert equals a contig must then be counted for that contig.</p>
 *
 * Created by edbuckler on 10/19/15.
 */
public class LoadRNAContigsToGBSDBPluginTest {

    private static final String RNA_DIR = GeneralConstants.DATA_DIR + "RNA/";
    private static final String FASTQ_FILE = RNA_DIR + "5452_7332_22978_HHM5MBGXX_150707_1_Leaves_1_R1.fastq";
    // Maps the FASTQ above to taxon A188, tissue Leaf. Its other rows name files that aren't in the archive.
    private static final String KEY_FILE = RNA_DIR + "Key_file_smallRNA_test_plugin.txt";
    private static final String TAXON = "A188";
    private static final String TISSUE = "Leaf";

    private static final Path TEMP_DIR = Paths.get(GeneralConstants.TEMP_DIR, "LoadRNAContigsToGBSDBPluginTest");
    private static final String ADAPTER_START = "AGATCGGA";
    private static final int WORD_LENGTH = 16;
    private static final int MIN_CONTIG_LENGTH = 20;
    private static final int NUM_CONTIGS = 20;

    @Test
    public void testPlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        List<String> contigs = selectContigs();
        Path contigFile = writeContigFasta("contigs.fa", contigs, true);
        String dbName = TEMP_DIR.resolve("loadContigs.db").toString();

        String returnedDB = new LoadRNAContigsToGBSDBPlugin()
                .deleteOldData(true)
                .contigFile(contigFile.toString())
                .outputDB(dbName)
                .runPlugin(null);
        assertEquals(dbName, returnedDB);

        TagDataSQLite tdw = new TagDataSQLite(dbName);
        try {
            Set<String> loaded = tdw.getTags().stream().map(Tag::sequence).collect(Collectors.toSet());
            assertEquals("The contig containing an N should be skipped", new HashSet<>(contigs), loaded);

            Map<Tag, String> names = tdw.getTagsNameMap();
            for (int i = 0; i < contigs.size(); i++) {
                Tag tag = TagBuilder.instance(contigs.get(i)).build();
                assertEquals("Each tag keeps its FASTA header as its name", contigHeader(i), names.get(tag));
            }
        } finally {
            tdw.close();
        }

        // Loading again with deleteOldData replaces the database instead of adding to it.
        new LoadRNAContigsToGBSDBPlugin()
                .deleteOldData(true)
                .contigFile(writeContigFasta("contigsSubset.fa", contigs.subList(0, 5), false).toString())
                .outputDB(dbName)
                .runPlugin(null);
        tdw = new TagDataSQLite(dbName);
        try {
            assertEquals(5, tdw.getTags().size());
        } finally {
            tdw.close();
        }
    }

    @Test
    public void testKmerMatch() throws Exception {
        List<String> contigs = selectContigs();
        String dbName = loadContigs("kmerMatch.db", contigs);

        TagDataSQLite tdw = new TagDataSQLite(dbName);
        Set<Tag> allContigs;
        try {
            allContigs = tdw.getTags();
        } finally {
            tdw.close();
        }

        for (FindMatchByWordHash.MatchType matchType : FindMatchByWordHash.MatchType.values()) {
            FindMatchByWordHash matcher = FindMatchByWordHash.getBuilder(allContigs)
                    .matchType(matchType)
                    .maxWordCopies(100)
                    .wordLength(WORD_LENGTH)
                    .build();
            assertEquals(contigs.size(), matcher.totalNumberOfTags());
            for (String contig : contigs) {
                FindMatchByWordHash.Match match = matcher.match(contig);
                assertTrue(matchType + " should find " + contig, match.isPresent());
                assertEquals(matchType + " matched the wrong contig", contig, match.sequence());
            }
        }

        FindMatchByWordHash matcher = FindMatchByWordHash.getBuilder(allContigs)
                .matchType(FindMatchByWordHash.MatchType.MODAL_KMER)
                .maxWordCopies(100)
                .wordLength(WORD_LENGTH)
                .build();
        for (String contig : contigs) {
            FindMatchByWordHash.Match reverse = matcher.match(BaseEncoder.getReverseComplement(contig));
            assertEquals("Reverse complement should match " + contig, contig, reverse.sequence());
            assertFalse("Reverse complement match should be flagged as reverse", reverse.direction());
            assertTrue(matcher.match(contig).direction());
        }
        assertTrue("A poly-A read shouldn't match any contig", matcher.match("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").isEmpty());
    }

    @Test
    public void testRNADeMultiplexProductionPlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        List<String> contigs = selectContigs();
        String dbName = loadContigs("deMultiplex.db", contigs);

        String returnedDB = new RNADeMultiplexProductionPlugin()
                .keyFile(KEY_FILE)
                .inputDir(RNA_DIR)
                .inputDB(dbName)
                .batchSize(4)
                .wordSizeMatching(WORD_LENGTH)
                .maxWordRepeats(100)
                .matchingType(FindMatchByWordHash.MatchType.MODAL_KMER)
                .runPlugin(null);
        assertEquals(dbName, returnedDB);

        Map<String, Integer> exactReads = countReadsByInsert();
        TagDataSQLite tdw = new TagDataSQLite(dbName);
        try {
            Set<String> taxa = new HashSet<>();
            tdw.getTaxaList().forEach(taxon -> taxa.add(taxon.getName()));
            assertEquals(new HashSet<>(Arrays.asList("A188", "F7", "NC250", "NC348", "B68")), taxa);
            assertEquals(new HashSet<>(Arrays.asList("Leaf", "Groot")), tdw.getAllTissue());

            TableReport report = tdw.getAllCountsForTaxonTissue(new Taxon(TAXON), TISSUE);
            assertNotNull(report);
            assertEquals("One row per contig", contigs.size(), report.getRowCount());
            Map<String, Integer> countBySequence = new HashMap<>();
            long totalCount = 0;
            for (long row = 0; row < report.getRowCount(); row++) {
                String sequence = report.getValueAt(row, 0).toString();
                int count = ((Number) report.getValueAt(row, 2)).intValue();
                assertEquals(contigHeader(contigs.indexOf(sequence)), report.getValueAt(row, 1));
                countBySequence.put(sequence, count);
                totalCount += count;
            }
            for (String contig : contigs) {
                int expectedAtLeast = exactReads.get(contig);
                assertTrue("Every read equal to " + contig + " (" + expectedAtLeast + ") should be counted, got "
                        + countBySequence.get(contig), countBySequence.get(contig) >= expectedAtLeast);
            }
            assertTrue("Can't count more reads than the FASTQ holds", totalCount <= countReads());

            String topContig = contigs.get(0);
            Phenotype perTag = tdw.getAllCountsForTagTissue(TagBuilder.instance(topContig).build(), TISSUE);
            assertEquals("Only " + TAXON + " has a FASTQ in the archive", 1, perTag.numberOfObservations());
            assertEquals(TAXON, perTag.taxa().get(0).getName());
            assertEquals(countBySequence.get(topContig).floatValue(), ((Number) perTag.value(0, 1)).floatValue(), 0f);

            TableReport noReads = tdw.getAllCountsForTaxonTissue(new Taxon("F7"), TISSUE);
            assertEquals("F7's files aren't in the archive, so it has no counts", 0, noReads.getRowCount());
        } finally {
            tdw.close();
        }
    }

    private static String loadContigs(String dbFileName, List<String> contigs) throws IOException {
        String dbName = TEMP_DIR.resolve(dbFileName).toString();
        new LoadRNAContigsToGBSDBPlugin()
                .deleteOldData(true)
                .contigFile(writeContigFasta(dbFileName + ".fa", contigs, false).toString())
                .outputDB(dbName)
                .runPlugin(null);
        return dbName;
    }

    /**
     * The most common inserts of at least {@link #MIN_CONTIG_LENGTH} bases, skipping any that
     * share a word with a contig already chosen, so each read can match only one contig.
     */
    private static List<String> selectContigs() throws IOException {
        Map<String, Integer> counts = countReadsByInsert();
        List<String> candidates = counts.keySet().stream()
                .filter(insert -> insert.length() >= MIN_CONTIG_LENGTH && insert.matches("[ACGT]+")
                        && !insert.matches(".*(A{" + WORD_LENGTH + "}|C{" + WORD_LENGTH + "}|G{" + WORD_LENGTH + "}|T{" + WORD_LENGTH + "}).*"))
                .sorted(Comparator.comparing((String insert) -> -counts.get(insert)).thenComparing(insert -> insert))
                .collect(Collectors.toList());
        List<String> contigs = new ArrayList<>();
        Set<String> usedWords = new HashSet<>();
        for (String candidate : candidates) {
            Set<String> words = words(candidate);
            if (!Collections.disjoint(words, usedWords)) continue;
            usedWords.addAll(words);
            contigs.add(candidate);
            if (contigs.size() == NUM_CONTIGS) break;
        }
        assertEquals("The archive FASTQ should yield enough distinct inserts", NUM_CONTIGS, contigs.size());
        return contigs;
    }

    private static Set<String> words(String sequence) {
        Set<String> words = new HashSet<>();
        String reverse = BaseEncoder.getReverseComplement(sequence);
        for (int i = 0; i + WORD_LENGTH <= sequence.length(); i++) {
            words.add(sequence.substring(i, i + WORD_LENGTH));
            words.add(reverse.substring(i, i + WORD_LENGTH));
        }
        return words;
    }

    /** Count of reads by insert: the read up to the adapter. Reads without an adapter are left out. */
    private static Map<String, Integer> countReadsByInsert() throws IOException {
        Map<String, Integer> counts = new HashMap<>();
        for (String read : readSequences()) {
            int adapter = read.indexOf(ADAPTER_START);
            if (adapter <= 0) continue;
            counts.merge(read.substring(0, adapter), 1, Integer::sum);
        }
        return counts;
    }

    private static int countReads() throws IOException {
        return readSequences().size();
    }

    private static List<String> readSequences() throws IOException {
        List<String> sequences = new ArrayList<>();
        try (BufferedReader br = Utils.getBufferedReader(FASTQ_FILE)) {
            while (br.readLine() != null) {
                sequences.add(br.readLine());
                br.readLine();
                br.readLine();
            }
        }
        return sequences;
    }

    private static Path writeContigFasta(String fileName, List<String> contigs, boolean addContigWithN) throws IOException {
        Files.createDirectories(TEMP_DIR);
        Path fasta = TEMP_DIR.resolve(fileName);
        try (BufferedWriter bw = Utils.getBufferedWriter(fasta.toString())) {
            for (int i = 0; i < contigs.size(); i++) {
                String contig = contigs.get(i);
                bw.write(contigHeader(i) + "\n");
                // Wrap the sequence so the plugin has to join lines.
                bw.write(contig.substring(0, 10) + "\n" + contig.substring(10) + "\n");
            }
            if (addContigWithN) {
                bw.write(">contigWithN\nACGTACGTACGTNACGTACGTACGT\n");
            }
        }
        return fasta;
    }

    private static String contigHeader(int index) {
        return ">contig" + index;
    }

}
