/*
 * GBSv1SimData
 */
package net.maizegenetics.analysis.gbs;

import net.maizegenetics.analysis.gbs.v2.GBSSimData;
import net.maizegenetics.util.Utils;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/**
 * Deterministic, aligner-free builder of legacy GBSv1 pipeline inputs.
 *
 * <p>The GBSv1 pipeline tests historically depended on a multi-hundred-megabyte raw FASTQ download
 * (not part of the {@code tassel_test_data} release) plus byte-exact {@code .h5}/{@code .topm} golden
 * fixtures and MD5 hashes. This helper mirrors the GBSv2 rehabilitation approach in
 * {@link GBSSimData}: it reuses the same deterministic reference/key/FASTQ that {@code GBSSimData}
 * generates from a fixed seed, and drives the real v1 plugins to reconstruct the intermediate
 * artifacts (master TagCounts, TOPM, TBT) so the tests can assert on pipeline properties instead of
 * golden hashes.</p>
 *
 * <p>Reconstructed chain (all under a git-ignored temp directory, no external aligner):</p>
 * <ul>
 *   <li>{@link #buildMasterTagCounts()} - {@code FastqToTagCountPlugin} then
 *       {@code MergeMultipleTagCountPlugin} produce the master {@code .cnt} tag list.</li>
 *   <li>{@link #buildTopm()} - a SAM is synthesized directly from the known tag/reference positions
 *       (every tag matches its own reference locus, exactly as a perfect aligner would report it),
 *       then {@code SAMConverterPlugin} converts it to a TOPM. The SAM query names carry the
 *       {@code length=..count=..} form that {@code SAMConverterPlugin} parses for the tag length.</li>
 *   <li>{@link #buildTbt()} - {@code SeqToTBTHDF5Plugin} builds a taxa-oriented TBT HDF5.</li>
 *   <li>{@link #pivotTbt()} - {@code ModifyTBTHDF5Plugin} transposes it into a tag-oriented TBT.</li>
 *   <li>{@link #runDiscovery()} - {@code DiscoverySNPCallerPlugin} calls SNPs into an output TOPM.</li>
 * </ul>
 */
public class GBSv1SimData {

    public static final String ENZYME = GBSSimData.ENZYME; // ApeKI
    public static final int START_CHR = 9;
    public static final int END_CHR = 10;

    public final GBSSimData sim;
    public final Path v1Dir;
    public final Path tagCountsDir;
    public final Path masterTagCounts;
    public final Path samFile;
    public final Path topmFile;
    public final Path tbtFile;
    public final Path pivotedTbtFile;
    public final Path discoveryTopmFile;
    /**
     * A copy of the simulated key file with column D renamed to "DNASample". {@code GBSSimData}
     * labels that column "FullSampleName" (accepted by GBSv2 / positional {@code ParseBarcodeRead}),
     * but the legacy {@code ProductionSNPCallerPlugin} strictly requires the header "DNASample" or
     * "Sample". The name ends in {@code _key.txt} as {@code ProductionPipelineMain} expects.
     */
    public final Path keyFileV1;

    public GBSv1SimData(GBSSimData sim) throws IOException {
        this.sim = sim;
        this.v1Dir = sim.baseDir.resolve("v1");
        Files.createDirectories(v1Dir);
        this.tagCountsDir = v1Dir.resolve("tagCounts");
        this.masterTagCounts = v1Dir.resolve("master.cnt");
        this.samFile = v1Dir.resolve("tagsForAlign.sam");
        this.topmFile = v1Dir.resolve("master.topm");
        this.tbtFile = v1Dir.resolve("tbt.h5");
        this.pivotedTbtFile = v1Dir.resolve("tbt_pivoted.h5");
        this.discoveryTopmFile = v1Dir.resolve("discovery.topm");
        this.keyFileV1 = v1Dir.resolve("prod_key.txt");
        writeV1KeyFile();
    }

    private void writeV1KeyFile() throws IOException {
        List<String> lines = Files.readAllLines(sim.keyFile);
        try (BufferedWriter bw = Utils.getBufferedWriter(keyFileV1.toString())) {
            for (int i = 0; i < lines.size(); i++) {
                if (i == 0) {
                    String[] header = lines.get(i).split("\t", -1);
                    if (header.length > 3) {
                        header[3] = "DNASample";
                    }
                    bw.write(String.join("\t", header));
                } else {
                    bw.write(lines.get(i));
                }
                bw.write("\n");
            }
        }
    }

    /** Create a fresh simulated data set under {@code tempDir/GBS/sim/<label>/} plus its v1 scaffolding. */
    public static GBSv1SimData createUnder(String label) throws IOException {
        return new GBSv1SimData(GBSSimData.createUnder(label));
    }

    /** Run {@code FastqToTagCountPlugin} + {@code MergeMultipleTagCountPlugin} to build the master {@code .cnt}. */
    public GBSv1SimData buildMasterTagCounts() throws IOException {
        Files.createDirectories(tagCountsDir);

        FastqToTagCountPlugin toTagCount = new FastqToTagCountPlugin();
        toTagCount.setParameters(new String[]{
                "-i", sim.fastqDir.toString(),
                "-o", tagCountsDir.toString(),
                "-k", sim.keyFile.toString(),
                "-e", ENZYME,
                "-s", "150000000",
                "-c", "1"
        });
        toTagCount.performFunction(null);

        MergeMultipleTagCountPlugin merge = new MergeMultipleTagCountPlugin();
        merge.setParameters(new String[]{
                "-i", tagCountsDir.toString(),
                "-o", masterTagCounts.toString(),
                "-c", "1"
        });
        merge.performFunction(null);
        return this;
    }

    /**
     * Synthesize a perfect-alignment SAM straight from the known tag/reference positions, then convert
     * it to a TOPM with {@code SAMConverterPlugin}. No external aligner or golden SAM is required.
     */
    public GBSv1SimData buildTopm() throws IOException {
        writeSam();
        SAMConverterPlugin converter = new SAMConverterPlugin();
        converter.setParameters(new String[]{
                "-i", samFile.toString(),
                "-o", topmFile.toString()
        });
        converter.performFunction(null);
        return this;
    }

    private void writeSam() throws IOException {
        try (BufferedWriter bw = Utils.getBufferedWriter(samFile.toString())) {
            for (GBSSimData.TagInfo info : sim.tagInfos) {
                String seq = info.sequence;
                int len = seq.length();
                // Query name mimics the TagCountToFastqPlugin header (length=..count=..) that
                // SAMConverterPlugin parses to recover the tag length. BWA-style optional fields
                // (X0/NM) mark a single perfect hit so the tag is recorded as uniquely aligned.
                bw.write("length=" + len + "count=1"
                        + "\t0"
                        + "\t" + info.chrom
                        + "\t" + info.cutPosition
                        + "\t60"
                        + "\t" + len + "M"
                        + "\t*\t0\t0"
                        + "\t" + seq
                        + "\t*"
                        + "\tX0:i:1\tNM:i:0\n");
            }
        }
    }

    /** Run {@code SeqToTBTHDF5Plugin} to build a taxa-oriented TBT HDF5 from the master {@code .cnt}. */
    public GBSv1SimData buildTbt() {
        SeqToTBTHDF5Plugin plugin = new SeqToTBTHDF5Plugin();
        plugin.setParameters(new String[]{
                "-i", sim.fastqDir.toString(),
                "-k", sim.keyFile.toString(),
                "-e", ENZYME,
                "-o", tbtFile.toString(),
                "-s", "150000000",
                "-L", v1Dir.resolve("seqToTBT.log").toString(),
                "-t", masterTagCounts.toString()
        });
        plugin.performFunction(null);
        return this;
    }

    /** Run {@code ModifyTBTHDF5Plugin -p} to transpose the taxa-oriented TBT into a tag-oriented TBT. */
    public GBSv1SimData pivotTbt() {
        ModifyTBTHDF5Plugin plugin = new ModifyTBTHDF5Plugin();
        plugin.setParameters(new String[]{
                "-o", tbtFile.toString(),
                "-p", pivotedTbtFile.toString()
        });
        plugin.performFunction(null);
        return this;
    }

    /** Run {@code DiscoverySNPCallerPlugin} on the pivoted TBT + TOPM, writing an output TOPM with variants. */
    public GBSv1SimData runDiscovery() {
        DiscoverySNPCallerPlugin plugin = new DiscoverySNPCallerPlugin();
        plugin.setParameters(new String[]{
                "-i", pivotedTbtFile.toString(),
                "-m", topmFile.toString(),
                "-o", discoveryTopmFile.toString(),
                "-ref", sim.referenceFasta.toString(),
                "-sC", String.valueOf(START_CHR),
                "-eC", String.valueOf(END_CHR),
                "-mnMAF", "0.01",
                "-mnMAC", "1",
                "-mnLCov", "0.1",
                "-mnF", "-2.0"
        });
        plugin.performFunction(null);
        return this;
    }

    /** Convenience: build the master tag counts and the TOPM (the inputs shared by most tests). */
    public GBSv1SimData buildTagCountsAndTopm() throws IOException {
        buildMasterTagCounts();
        buildTopm();
        return this;
    }
}
