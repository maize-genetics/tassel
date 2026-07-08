package net.maizegenetics.analysis.gbs.v2;

import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagData;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.util.Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

/**
 * Deterministic generator of a tiny, self-contained GBSv2 data set for unit tests.
 *
 * <p>Historically the GBSv2 tests depended on multi-megabyte FASTQ/reference/SAM fixtures
 * that had to be downloaded separately, plus an external aligner (bowtie2/BWA). That made the
 * tests un-runnable in CI and they were excluded from the build. This class produces everything
 * the pipeline needs from a fixed random seed, so the tests can build their own {@code GBSv2.db}
 * and assert against known properties (tag counts, injected SNP counts) rather than golden hashes.</p>
 *
 * <p>What it writes (all under a caller-supplied, git-ignored temp directory):</p>
 * <ul>
 *   <li>A reference FASTA with chromosomes {@code 9} and {@code 10}. Each locus embeds an ApeKI
 *       {@code GCWGC} cut site ({@code G} + a tag beginning with the {@code CAGC} cut remnant).</li>
 *   <li>A tab-delimited key file ({@code Flowcell/Lane/Barcode/FullSampleName/.../LibraryPrepID/Enzyme})
 *       compatible with both GBSv2 ({@code TaxaListIOUtils}) and legacy GBSv1 ({@code ParseBarcodeRead},
 *       which requires a non-empty column&nbsp;H).</li>
 *   <li>A gzipped FASTQ of reads shaped as {@code barcode + CAGC + genomic-kmer}, replicated across
 *       taxa above the min-kmer-count/replication thresholds, with single-base variants injected at
 *       SNP loci so that {@code DiscoverySNPCallerPluginV2} calls real SNPs.</li>
 *   <li>An aligner-free SAM produced by exact-matching each exported tag ({@code @tagSeq=} header)
 *       back to its known reference position, so every SAM tag already exists in the DB.</li>
 * </ul>
 *
 * <p>The public metadata fields ({@link #expectedDistinctTags}, {@link #numLoci},
 * {@link #numSnpLoci}, {@link #taxa}) let tests make property-based assertions.</p>
 *
 * @author generated for TAS GBS test rehabilitation
 */
public class GBSSimData {

    public static final String FLOWCELL = "SIMFC";
    public static final String LANE = "1";
    public static final String ENZYME = "ApeKI";
    public static final String CUT_REMNANT = "CAGC";
    public static final int TAG_LENGTH = 64;
    public static final int GENOMIC_LENGTH = TAG_LENGTH - CUT_REMNANT.length(); // 60
    /** Matches the {@code minKmerCount} used by the tests. */
    public static final int MIN_KMER_COUNT = 5;
    /** Depth (reads) per taxon for each tag. */
    public static final int DEPTH_PER_TAXON = 3;
    /** 0-based index within the tag where the injected SNP is placed. */
    public static final int SNP_TAG_INDEX = 40;
    /** Number of taxa carrying each allele at a locus. */
    private static final int TAXA_PER_ALLELE = 4;

    private static final String[] FORBIDDEN = {"GCAGC", "GCTGC", "GCAGAGAT", "GCTGAGAT"};
    private static final char[] BASES = {'A', 'C', 'G', 'T'};

    // Eight distinct fixed-length (6bp) barcodes; none a prefix of another.
    private static final String[] BARCODES = {
            "AACCGA", "AAGGTC", "ACACGT", "ACGTAC",
            "AGACCT", "AGGTCA", "ATCGGA", "ATGCAC"
    };

    /** Metadata for a single simulated tag. */
    public static class TagInfo {
        public final String sequence;
        public final String chrom;
        public final int cutPosition; // 1-based
        public final boolean variant; // true for the alternate allele tag at a SNP locus
        public final int[] taxaIndices; // taxa (by position in {@link #taxa}) carrying this tag

        TagInfo(String sequence, String chrom, int cutPosition, boolean variant, int[] taxaIndices) {
            this.sequence = sequence;
            this.chrom = chrom;
            this.cutPosition = cutPosition;
            this.variant = variant;
            this.taxaIndices = taxaIndices;
        }
    }

    public final Path baseDir;
    public final Path referenceFasta;
    public final Path keyFile;
    public final Path fastqDir;
    public final Path fastqFile;
    public final Path dbFile;
    public final Path exportFasta;
    public final Path samFile;

    public final List<String> taxa = new ArrayList<>();
    public final List<String> barcodes = new ArrayList<>();
    public final List<TagInfo> tagInfos = new ArrayList<>();

    public int expectedDistinctTags;
    public int numLoci;
    public int numSnpLoci;

    private final Random rng = new Random(20240521L);
    private final Set<String> usedTagSequences = new HashSet<>();
    private final Map<String, TagInfo> seqToTagInfo = new LinkedHashMap<>();

    // Reference sequence builders keyed by chromosome name.
    private final Map<String, StringBuilder> refByChrom = new LinkedHashMap<>();

    public GBSSimData(Path baseDir) {
        this.baseDir = baseDir;
        this.referenceFasta = baseDir.resolve("sim_reference.fa");
        this.keyFile = baseDir.resolve("sim_key.txt");
        this.fastqDir = baseDir.resolve("fastq");
        this.fastqFile = fastqDir.resolve(FLOWCELL + "_" + LANE + "_fastq.gz");
        this.dbFile = baseDir.resolve("GBSv2.db");
        this.exportFasta = baseDir.resolve("tagsForAlign.fa.gz");
        this.samFile = baseDir.resolve("tagsForAlign.sam");
    }

    /**
     * Generate all inputs (reference FASTA, key file, FASTQ) and populate metadata.
     * SNP loci come first per chromosome, followed by monomorphic loci.
     */
    public GBSSimData writeInputs() throws IOException {
        Files.createDirectories(baseDir);
        Files.createDirectories(fastqDir);
        for (int i = 0; i < BARCODES.length; i++) {
            taxa.add("SIMTAXON" + i);
            barcodes.add(BARCODES[i]);
        }

        String[] chroms = {"9", "10"};
        int snpLociPerChrom = 2;
        int monoLociPerChrom = 2;
        List<int[][]> readPlan = new ArrayList<>(); // per tag: nothing needed, handled inline

        List<String[]> fastqRecords = new ArrayList<>();

        for (String chrom : chroms) {
            refByChrom.put(chrom, new StringBuilder());
            for (int locus = 0; locus < snpLociPerChrom + monoLociPerChrom; locus++) {
                boolean snpLocus = locus < snpLociPerChrom;
                addLocus(chrom, snpLocus, fastqRecords);
            }
        }

        numLoci = (snpLociPerChrom + monoLociPerChrom) * chroms.length;
        numSnpLoci = snpLociPerChrom * chroms.length;
        expectedDistinctTags = tagInfos.size();

        writeReference();
        writeKeyFile();
        writeFastq(fastqRecords);
        return this;
    }

    private void addLocus(String chrom, boolean snpLocus, List<String[]> fastqRecords) {
        StringBuilder ref = refByChrom.get(chrom);
        // leading filler
        appendRandom(ref, 100);
        // preceding G forms the ApeKI GCWGC cut site with the tag's leading CAGC remnant
        ref.append('G');
        int cutPos = ref.length() + 1; // 1-based position where the tag begins
        String refTag = newTag();
        ref.append(refTag);
        // trailing filler
        appendRandom(ref, 35);

        int[] firstAllele = new int[TAXA_PER_ALLELE];
        for (int i = 0; i < TAXA_PER_ALLELE; i++) firstAllele[i] = i;
        registerTag(refTag, chrom, cutPos, false, firstAllele, fastqRecords);

        if (snpLocus) {
            String variantTag = makeVariant(refTag);
            int[] secondAllele = new int[TAXA_PER_ALLELE];
            for (int i = 0; i < TAXA_PER_ALLELE; i++) secondAllele[i] = TAXA_PER_ALLELE + i;
            registerTag(variantTag, chrom, cutPos, true, secondAllele, fastqRecords);
        }
    }

    private void registerTag(String tagSeq, String chrom, int cutPos, boolean variant,
                             int[] taxaIndices, List<String[]> fastqRecords) {
        TagInfo info = new TagInfo(tagSeq, chrom, cutPos, variant, taxaIndices);
        tagInfos.add(info);
        seqToTagInfo.put(tagSeq, info);
        // emit reads: for each taxon, DEPTH_PER_TAXON reads = barcode + tag
        for (int taxonIdx : taxaIndices) {
            String read = barcodes.get(taxonIdx) + tagSeq;
            for (int d = 0; d < DEPTH_PER_TAXON; d++) {
                fastqRecords.add(new String[]{read});
            }
        }
    }

    /** Create a new, globally-unique, ApeKI-safe tag of the form CAGC + 60 genomic bases. */
    private String newTag() {
        while (true) {
            StringBuilder sb = new StringBuilder(CUT_REMNANT);
            for (int i = 0; i < GENOMIC_LENGTH; i++) sb.append(randomBase());
            String tag = sb.toString();
            if (isTagSafe(tag) && usedTagSequences.add(tag)) {
                return tag;
            }
        }
    }

    /** Build a single-base variant of a reference tag at {@link #SNP_TAG_INDEX}. */
    private String makeVariant(String refTag) {
        while (true) {
            char orig = refTag.charAt(SNP_TAG_INDEX);
            char repl;
            do {
                repl = randomBase();
            } while (repl == orig);
            String variant = refTag.substring(0, SNP_TAG_INDEX) + repl + refTag.substring(SNP_TAG_INDEX + 1);
            if (isTagSafe(variant) && usedTagSequences.add(variant)) {
                return variant;
            }
        }
    }

    /**
     * A tag is "safe" when GBSSeqToTagDBPlugin will store it at full length: it must not begin with
     * an overlapping ApeKI cut site and must not contain any likely-read-end string (which would
     * truncate the tag in {@code removeSecondCutSiteIndexOf}).
     */
    private boolean isTagSafe(String tag) {
        if (tag.length() != TAG_LENGTH) return false;
        if (!tag.startsWith(CUT_REMNANT)) return false;
        if (tag.startsWith("CAGCTGC") || tag.startsWith("CTGCAGC")) return false;
        for (String bad : FORBIDDEN) {
            if (tag.indexOf(bad) >= 0) return false;
        }
        return true;
    }

    private void appendRandom(StringBuilder sb, int count) {
        for (int i = 0; i < count; i++) sb.append(randomBase());
    }

    private char randomBase() {
        return BASES[rng.nextInt(BASES.length)];
    }

    private void writeReference() throws IOException {
        try (BufferedWriter bw = Utils.getBufferedWriter(referenceFasta.toString())) {
            for (Map.Entry<String, StringBuilder> e : refByChrom.entrySet()) {
                bw.write(">" + e.getKey() + "\n");
                String seq = e.getValue().toString();
                // wrap at 60 chars for readability
                for (int i = 0; i < seq.length(); i += 60) {
                    bw.write(seq.substring(i, Math.min(seq.length(), i + 60)));
                    bw.write("\n");
                }
            }
        }
    }

    private void writeKeyFile() throws IOException {
        try (BufferedWriter bw = Utils.getBufferedWriter(keyFile.toString())) {
            // Column H (index 7) must be a non-empty LibraryPrepID for legacy ParseBarcodeRead.
            bw.write("Flowcell\tLane\tBarcode\tFullSampleName\tPlateName\tRow\tColumn\tLibraryPrepID\tEnzyme\n");
            for (int i = 0; i < taxa.size(); i++) {
                bw.write(FLOWCELL + "\t" + LANE + "\t" + barcodes.get(i) + "\t" + taxa.get(i)
                        + "\tPlate1\tA\t" + (i + 1) + "\t" + (100 + i) + "\t" + ENZYME + "\n");
            }
        }
    }

    private void writeFastq(List<String[]> fastqRecords) throws IOException {
        try (BufferedWriter bw = Utils.getBufferedWriter(fastqFile.toString())) {
            int readNum = 0;
            for (String[] rec : fastqRecords) {
                String read = rec[0];
                readNum++;
                // Header has >=5 colon-delimited fields so quality base is read as Phred+33.
                bw.write("@" + FLOWCELL + ":" + LANE + ":1:1000:" + readNum + ":1:1\n");
                bw.write(read + "\n");
                bw.write("+\n");
                StringBuilder qual = new StringBuilder();
                for (int i = 0; i < read.length(); i++) qual.append('I'); // Phred+33 => Q40
                bw.write(qual.toString() + "\n");
            }
        }
    }

    /** Run {@code GBSSeqToTagDBPlugin} to build the tag/taxa DB from the simulated inputs. */
    public GBSSimData buildTagDB() {
        new GBSSeqToTagDBPlugin()
                .enzyme(ENZYME)
                .inputDirectory(fastqDir.toString())
                .outputDatabaseFile(dbFile.toString())
                .keyFile(keyFile.toString())
                .kmerLength(TAG_LENGTH)
                .minKmerCount(MIN_KMER_COUNT)
                .minimumQualityScore(20)
                .deleteOldData(true)
                .performFunction(null);
        return this;
    }

    /**
     * Export tags to FASTQ ({@code @tagSeq=} headers) then synthesize a SAM by matching each exported
     * tag to its known reference position. Guarantees every SAM tag already exists in the DB.
     */
    public GBSSimData exportAndAlign() throws IOException {
        new TagExportToFastqPlugin()
                .inputDB(dbFile.toString())
                .outputFile(exportFasta.toString())
                .performFunction(null);
        writeSamFromExportedTags();
        return this;
    }

    private void writeSamFromExportedTags() throws IOException {
        try (BufferedReader br = Utils.getBufferedReader(exportFasta.toString());
             BufferedWriter bw = Utils.getBufferedWriter(samFile.toString())) {
            String line;
            while ((line = br.readLine()) != null) {
                if (!line.startsWith("@tagSeq=")) continue;
                String seq = line.substring("@tagSeq=".length()).trim();
                TagInfo info = seqToTagInfo.get(seq);
                if (info == null) {
                    // Unexpected tag - skip so we don't poison the DB import.
                    continue;
                }
                // name flag chr pos mapq cigar rnext pnext tlen seq qual [optional]
                bw.write("tagSeq=" + seq + "\t0\t" + info.chrom + "\t" + info.cutPosition
                        + "\t60\t" + TAG_LENGTH + "M\t*\t0\t0\t" + seq + "\t*"
                        + "\tNM:i:0\tMD:Z:" + TAG_LENGTH + "\tAS:i:128\n");
            }
        }
    }

    /** Import the synthesized SAM into the DB (minMAPQ 2 to mirror the real pipeline). */
    public GBSSimData importSam() {
        new SAMToGBSdbPlugin()
                .gBSDBFile(dbFile.toString())
                .sAMInputFile(samFile.toString())
                .minMAPQ(2)
                .deleteOldData(true)
                .performFunction(null);
        return this;
    }

    /** Convenience: build the DB and load alignments through the synthesized SAM. */
    public GBSSimData buildDatabaseThroughSam() throws IOException {
        buildTagDB();
        exportAndAlign();
        importSam();
        return this;
    }

    /** Run SNP discovery on the simulated data, optionally supplying the simulated reference. */
    public GBSSimData runDiscovery(boolean withReference) {
        // Ensure a clean, deterministic reference state regardless of prior tests in this JVM.
        DiscoverySNPCallerPluginV2 plugin = new DiscoverySNPCallerPluginV2()
                .inputDB(dbFile.toString())
                .minMinorAlleleFreq(0.1)
                .startChromosome(new Chromosome("9"))
                .endChromosome(new Chromosome("10"))
                .deleteOldData(true);
        if (withReference) {
            plugin.referenceGenomeFile(referenceFasta.toString());
        } else {
            DiscoverySNPCallerPluginV2.myRefSequence = null;
        }
        plugin.performFunction(null);
        return this;
    }

    /** Number of distinct tags currently stored in the DB. */
    public int dbTagCount() {
        TagData tagData = new TagDataSQLite(dbFile.toString());
        int size = tagData.getTags().size();
        try {
            ((TagDataSQLite) tagData).close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return size;
    }

    /** A known simulated tag sequence (the first reference-allele tag). */
    public String aKnownTagSequence() {
        return tagInfos.get(0).sequence;
    }

    /** A sequence guaranteed not to be a stored tag. */
    public String anUnknownTagSequence() {
        StringBuilder sb = new StringBuilder("CAGC");
        for (int i = 0; i < GENOMIC_LENGTH; i++) sb.append('A');
        String s = sb.toString();
        // extraordinarily unlikely to collide with any generated tag, but guard anyway
        return usedTagSequences.contains(s) ? s.substring(0, TAG_LENGTH - 1) + "T" : s;
    }

    /**
     * Convenience factory that creates a fresh simulated data set under
     * {@code tempDir/GBS/sim/<label>/} and writes all inputs.
     */
    public static GBSSimData createUnder(String label) throws IOException {
        Path base = Paths.get("tempDir", "GBS", "sim", label);
        if (Files.exists(base)) {
            deleteRecursively(base);
        }
        GBSSimData sim = new GBSSimData(base);
        sim.writeInputs();
        return sim;
    }

    private static void deleteRecursively(Path path) throws IOException {
        if (!Files.exists(path)) return;
        Files.walk(path)
                .sorted((a, b) -> b.getNameCount() - a.getNameCount())
                .forEach(p -> {
                    try {
                        Files.deleteIfExists(p);
                    } catch (IOException e) {
                        // best effort
                    }
                });
    }
}
