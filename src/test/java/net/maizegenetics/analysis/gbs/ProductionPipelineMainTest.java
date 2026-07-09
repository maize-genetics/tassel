/*
 * ProductionPipelineMainTest
 */
package net.maizegenetics.analysis.gbs;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.util.Utils;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

/**
 * Property-based rehabilitation of the legacy GBSv1 {@code ProductionPipelineMain} test.
 *
 * <p>The original test drove the full production pipeline against a downloaded raw-FASTQ data set and
 * asserted a single hardcoded MD5 of a {@code .hmp.txt.gz} output. That output format no longer
 * matches the current {@code ProductionSNPCallerPlugin} (which writes HDF5 genotypes), and the raw
 * FASTQ is not part of the test-data release. This rehabilitation self-generates all inputs with
 * {@link GBSv1SimData}, drives {@code ProductionPipelineMain} entirely offline (a local run/config
 * file pair, SMTP neutralized to localhost with short timeouts), and asserts on the produced HDF5
 * genotypes instead of a brittle checksum.</p>
 *
 * @author Dallas Kroon (original), rehabilitated for self-generated data
 */
public class ProductionPipelineMainTest {

    @Test
    public void testProductionPipelineMain() throws Exception {
        GBSv1SimData sim = GBSv1SimData.createUnder("ProdPipeMain");
        sim.buildMasterTagCounts();
        sim.buildTopm();
        sim.buildTbt();
        sim.pivotTbt();
        sim.runDiscovery();

        Path pipelineDir = sim.v1Dir.resolve("pipeline");
        Path runDir = pipelineDir.resolve("run");
        Path outputDir = pipelineDir.resolve("out");
        Path archiveDir = pipelineDir.resolve("archive");
        Files.createDirectories(runDir);
        Files.createDirectories(outputDir);
        Files.createDirectories(archiveDir);

        // A .run file describing this pipeline run, all pointing at self-generated inputs.
        Path runFile = runDir.resolve("temp.run");
        try (BufferedWriter bw = Utils.getBufferedWriter(runFile.toString())) {
            bw.write("inputFolder=" + sim.sim.fastqDir + "\n");
            bw.write("enzyme=" + GBSv1SimData.ENZYME + "\n");
            bw.write("topmFile=" + sim.discoveryTopmFile + "\n");
            bw.write("outputFolder=" + outputDir + "\n");
            bw.write("keyFile=" + sim.keyFileV1 + "\n");
        }

        // Application configuration. emailHost is neutralized to localhost (no live send).
        Path propsFile = pipelineDir.resolve("production_pipeline.properties");
        try (BufferedWriter bw = Utils.getBufferedWriter(propsFile.toString())) {
            bw.write("runFileSuffix=.run\n");
            bw.write("emailHost=localhost\n");
            bw.write("emailAddress=test@example.com\n");
            bw.write("runDirectory=" + runDir + "\n");
            bw.write("archiveDirectory=" + archiveDir + "\n");
            bw.write("haplosDirectory=" + pipelineDir + "\n");
        }

        // Keep any SMTP attempt from blocking the test.
        System.setProperty("mail.smtp.connectiontimeout", "2000");
        System.setProperty("mail.smtp.timeout", "2000");
        System.setProperty("mail.smtp.writetimeout", "2000");

        // ProductionPipelineMain permanently redirects System.out/err to a log file; save and restore.
        PrintStream originalOut = System.out;
        PrintStream originalErr = System.err;
        try {
            boolean runCheckSum = false, runImputation = false;
            new ProductionPipelineMain(propsFile.toString(), runCheckSum, runImputation, null);
        } finally {
            System.setOut(originalOut);
            System.setErr(originalErr);
        }

        // Output HDF5 is named after the key file ("prod_key.txt" -> "prod.h5").
        Path outFile = outputDir.resolve("prod.h5");
        assertTrue("Production pipeline did not produce " + outFile, Files.exists(outFile));

        GenotypeTable genos = ImportUtils.readGuessFormat(outFile.toString());
        assertNotNull("Failed to read back HDF5 genotypes", genos);
        assertEquals("HDF5 taxa count should match the simulated taxa",
                sim.sim.taxa.size(), genos.numberOfTaxa());
        assertTrue("HDF5 genotypes should contain at least one site", genos.numberOfSites() > 0);
    }
}
