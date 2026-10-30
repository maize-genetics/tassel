package net.maizegenetics.analysis.gbs.repgen;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.BaseEncoder;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Small, deterministic rAmpSeq data set for the RepGen plugin tests.
 * <p>
 * The reference is 4 kb from each of chromosomes 9 and 10 of
 * {@code dataFiles/GBS/ZmB73_RefGen_v2_chr9_10_1st20MB.fasta}, with one primer pair spliced in at three
 * 200 bp amplicons per chromosome (forward primer at the amplicon start, reverse complement of the reverse
 * primer at its end). Six taxa carry either the reference allele or a SNP allele at each amplicon, so the
 * loader should keep exactly two tags per amplicon. Each taxon gets a FASTQ file of 160 bp forward reads,
 * and some files also hold reads the loader must discard: a sequencing error seen once, an allele private
 * to one taxon, a low-count variant, a read shorter than the minimum kmer length and a read containing N.
 */
final class RepGenSimData {

    static final String GBS_REFERENCE = GeneralConstants.DATA_DIR + "GBS/ZmB73_RefGen_v2_chr9_10_1st20MB.fasta";
    static final String FORWARD_PRIMER = "GCACAAGTTGTCCTGCTTCC";
    static final String REVERSE_PRIMER = "ACGTCGATCTAGGGGGTCTC";
    static final String[] CHROMOSOMES = {"9", "10"};
    static final int REGION_START = 1_000_000;
    static final int CHROMOSOME_LENGTH = 4000;
    /** 1-based amplicon starts, used on both chromosomes. */
    static final int[] AMPLICON_STARTS = {501, 1701, 2901};
    static final int AMPLICON_LENGTH = 200;
    static final int READ_LENGTH = 160;
    static final int KMER_LENGTH = 150;
    static final int SNP_OFFSET = 75;
    static final int QUALITY_BASE = 33;
    static final String[] TAXA = {"B73", "CML247", "Ky21", "M37W", "Mo17", "Oh43"};
    /** Taxa carrying the SNP allele, by amplicon (chromosome 9 amplicons first). */
    private static final int[][] ALT_TAXA = {{0, 1, 2}, {3, 4, 5}, {0, 2, 4}, {1, 3, 5}, {0, 1, 4, 5}, {2, 3}};

    final Path dir;
    final Path referenceFile;
    final Path fastqDir;
    final Path keyFile;
    final Map<String, String> referenceByChromosome = new LinkedHashMap<>();
    final List<Locus> loci = new ArrayList<>();
    /** Tags the loader must keep, keyed by their 150 bp sequence. */
    final Map<String, ExpectedTag> expectedTags = new LinkedHashMap<>();
    /** Present in two taxa at depth 2: kept only when minKmerCount is 4 or less. */
    final ExpectedTag lowCountTag;
    /** Sequences the loader must drop whatever the count settings. */
    final Set<String> discardedSequences = new HashSet<>();

    static final class Locus {
        final String chromosome;
        final int start;
        final String refAmplicon;
        final String altAmplicon;
        final Set<Integer> altTaxa;

        private Locus(String chromosome, int start, String refAmplicon, Set<Integer> altTaxa) {
            this.chromosome = chromosome;
            this.start = start;
            this.refAmplicon = refAmplicon;
            this.altAmplicon = mutate(refAmplicon, SNP_OFFSET);
            this.altTaxa = altTaxa;
        }

        String tag(boolean alt) {
            return (alt ? altAmplicon : refAmplicon).substring(0, KMER_LENGTH);
        }
    }

    static final class ExpectedTag {
        final String sequence;
        final Locus locus;
        final boolean alt;
        final int[] depths = new int[TAXA.length];
        /** Quality strings of the reads, trimmed to the kmer length. */
        final List<String> qualities = new ArrayList<>();

        private ExpectedTag(String sequence, Locus locus, boolean alt) {
            this.sequence = sequence;
            this.locus = locus;
            this.alt = alt;
        }

        int totalDepth() {
            return Arrays.stream(depths).sum();
        }

        /** Per-base average quality, truncated the same way the loader averages integer scores. */
        String averageQuality() {
            StringBuilder sb = new StringBuilder();
            for (int pos = 0; pos < KMER_LENGTH; pos++) {
                int sum = 0;
                for (String quality : qualities) sum += quality.charAt(pos) - QUALITY_BASE;
                sb.append((char) (sum / qualities.size() + QUALITY_BASE));
            }
            return sb.toString();
        }
    }

    private final Map<Integer, List<String>> readsByTaxon = new HashMap<>();
    private int readCounter = 0;

    private RepGenSimData(Path dir) throws IOException {
        this.dir = dir;
        this.referenceFile = dir.resolve("reference.fa");
        this.fastqDir = dir.resolve("fastq");
        this.keyFile = dir.resolve("repGen_key.txt");

        Map<String, String> regions = readRegions();
        int locusIndex = 0;
        for (String chromosome : CHROMOSOMES) {
            StringBuilder reference = new StringBuilder(regions.get(chromosome));
            for (int start : AMPLICON_STARTS) {
                reference.replace(start - 1, start - 1 + FORWARD_PRIMER.length(), FORWARD_PRIMER);
                int reverseStart = start - 1 + AMPLICON_LENGTH - REVERSE_PRIMER.length();
                reference.replace(reverseStart, reverseStart + REVERSE_PRIMER.length(), BaseEncoder.getReverseComplement(REVERSE_PRIMER));
                String amplicon = reference.substring(start - 1, start - 1 + AMPLICON_LENGTH);
                Set<Integer> altTaxa = Arrays.stream(ALT_TAXA[locusIndex++]).boxed().collect(Collectors.toSet());
                loci.add(new Locus(chromosome, start, amplicon, altTaxa));
            }
            referenceByChromosome.put(chromosome, reference.toString());
        }

        for (int locusIdx = 0; locusIdx < loci.size(); locusIdx++) {
            Locus locus = loci.get(locusIdx);
            for (boolean alt : new boolean[]{false, true}) {
                expectedTags.put(locus.tag(alt), new ExpectedTag(locus.tag(alt), locus, alt));
            }
            for (int taxon = 0; taxon < TAXA.length; taxon++) {
                boolean alt = locus.altTaxa.contains(taxon);
                int depth = 3 + (taxon + locusIdx) % 3;
                String read = (alt ? locus.altAmplicon : locus.refAmplicon).substring(0, READ_LENGTH);
                ExpectedTag expected = expectedTags.get(locus.tag(alt));
                for (int r = 0; r < depth; r++) {
                    String quality = quality(read.length(), r);
                    addRead(taxon, read, quality);
                    expected.depths[taxon]++;
                    expected.qualities.add(quality.substring(0, KMER_LENGTH));
                }
            }
        }

        String sequencingError = mutate(loci.get(0).refAmplicon, 30).substring(0, READ_LENGTH);
        addRead(0, sequencingError, quality(READ_LENGTH, 0));
        discardedSequences.add(sequencingError.substring(0, KMER_LENGTH));

        String privateAllele = mutate(loci.get(2).refAmplicon, 100).substring(0, READ_LENGTH);
        for (int r = 0; r < 4; r++) addRead(5, privateAllele, quality(READ_LENGTH, r));
        discardedSequences.add(privateAllele.substring(0, KMER_LENGTH));

        String lowCount = mutate(loci.get(3).refAmplicon, 120).substring(0, READ_LENGTH);
        lowCountTag = new ExpectedTag(lowCount.substring(0, KMER_LENGTH), loci.get(3), false);
        for (int taxon : new int[]{1, 2}) {
            for (int r = 0; r < 2; r++) {
                String quality = quality(READ_LENGTH, r);
                addRead(taxon, lowCount, quality);
                lowCountTag.depths[taxon]++;
                lowCountTag.qualities.add(quality.substring(0, KMER_LENGTH));
            }
        }

        String shortRead = loci.get(4).refAmplicon.substring(0, 100);
        for (int r = 0; r < 5; r++) addRead(4, shortRead, quality(shortRead.length(), r));
        discardedSequences.add(shortRead);

        StringBuilder withN = new StringBuilder(loci.get(5).refAmplicon.substring(0, READ_LENGTH));
        withN.setCharAt(50, 'N');
        for (int r = 0; r < 5; r++) addRead(3, withN.toString(), quality(READ_LENGTH, r));
        discardedSequences.add(withN.substring(0, KMER_LENGTH));
    }

    /** Writes the data set under {@code tempDir/RepGen/<label>}, replacing anything already there. */
    static RepGenSimData createUnder(String label) throws IOException {
        Path dir = Paths.get(GeneralConstants.TEMP_DIR, "RepGen", label);
        deleteRecursively(dir);
        RepGenSimData data = new RepGenSimData(dir);
        data.write();
        return data;
    }

    String referenceSequence(String chromosome, int start, int end) {
        return referenceByChromosome.get(chromosome).substring(start - 1, end);
    }

    private void addRead(int taxon, String sequence, String quality) {
        readCounter++;
        String header = "@M01032:387:000000000-ANP68:1:1101:" + (10000 + readCounter) + ":" + (2000 + taxon) + " 1:N:0:1";
        readsByTaxon.computeIfAbsent(taxon, t -> new ArrayList<>())
                .add(header + "\n" + sequence + "\n+\n" + quality + "\n");
    }

    /** Q40 everywhere except the first base (Q20 or Q30 alternating by read) and the bases past the kmer (Q2). */
    private static String quality(int length, int readIndex) {
        StringBuilder sb = new StringBuilder();
        for (int pos = 0; pos < length; pos++) {
            if (pos == 0) sb.append(readIndex % 2 == 0 ? '5' : '?');
            else if (pos >= KMER_LENGTH) sb.append('#');
            else sb.append('I');
        }
        return sb.toString();
    }

    private static String mutate(String sequence, int position) {
        char base = sequence.charAt(position);
        char replacement = "ACGT".charAt(("ACGT".indexOf(base) + 2) % 4);
        return sequence.substring(0, position) + replacement + sequence.substring(position + 1);
    }

    static String fastqName(int taxon) {
        return TAXA[taxon] + "_ANP68_1_Locus9_R1.fastq";
    }

    private void write() throws IOException {
        Files.createDirectories(fastqDir);
        try (BufferedWriter bw = Files.newBufferedWriter(referenceFile)) {
            for (Map.Entry<String, String> entry : referenceByChromosome.entrySet()) {
                bw.write(">" + entry.getKey() + "\n");
                String sequence = entry.getValue();
                for (int i = 0; i < sequence.length(); i += 60) {
                    bw.write(sequence, i, Math.min(60, sequence.length() - i));
                    bw.write("\n");
                }
            }
        }
        List<String> key = new ArrayList<>();
        key.add("FileName\tFullSampleName");
        for (int taxon = 0; taxon < TAXA.length; taxon++) {
            key.add(fastqName(taxon) + "\t" + TAXA[taxon]);
            Files.write(fastqDir.resolve(fastqName(taxon)), String.join("", readsByTaxon.get(taxon)).getBytes());
        }
        Files.write(keyFile, key);
        // Not in the key file, so the loader must ignore it.
        String unkeyed = loci.get(0).refAmplicon.substring(0, READ_LENGTH);
        StringBuilder reads = new StringBuilder();
        for (int r = 0; r < 10; r++) {
            reads.append("@M01032:387:000000000-ANP68:1:1101:9").append(r).append(":1 1:N:0:1\n")
                    .append(mutate(unkeyed, 10)).append("\n+\n").append(quality(READ_LENGTH, r)).append("\n");
        }
        Files.write(fastqDir.resolve("Undetermined_ANP68_1_R1.fastq"), reads.toString().getBytes());
        discardedSequences.add(mutate(unkeyed, 10).substring(0, KMER_LENGTH));
    }

    /** {@value #CHROMOSOME_LENGTH} bp from {@value #REGION_START} on each chromosome of the GBS test reference. */
    private static Map<String, String> readRegions() throws IOException {
        Map<String, StringBuilder> regions = new HashMap<>();
        try (BufferedReader br = Files.newBufferedReader(Paths.get(GBS_REFERENCE))) {
            String chromosome = null;
            int position = 0;
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    chromosome = line.substring(1).trim().replaceFirst("^chr", "");
                    position = 0;
                    continue;
                }
                int lineStart = position;
                position += line.length();
                if (position <= REGION_START || lineStart >= REGION_START + CHROMOSOME_LENGTH) continue;
                int from = Math.max(0, REGION_START - lineStart);
                int to = Math.min(line.length(), REGION_START + CHROMOSOME_LENGTH - lineStart);
                regions.computeIfAbsent(chromosome, c -> new StringBuilder()).append(line, from, to);
            }
        }
        Map<String, String> result = new HashMap<>();
        for (String chromosome : CHROMOSOMES) {
            String region = regions.get(chromosome).toString().toUpperCase();
            if (region.length() != CHROMOSOME_LENGTH || !region.matches("[ACGT]+")) {
                throw new IllegalStateException("Expected " + CHROMOSOME_LENGTH + " unambiguous bases on chromosome " + chromosome);
            }
            result.put(chromosome, region);
        }
        return result;
    }

    static void deleteRecursively(Path dir) throws IOException {
        if (!Files.exists(dir)) return;
        try (Stream<Path> paths = Files.walk(dir)) {
            for (Path path : paths.sorted(Comparator.reverseOrder()).collect(Collectors.toList())) {
                Files.delete(path);
            }
        }
    }

}
