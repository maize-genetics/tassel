package net.maizegenetics.analysis.gbs.v2;

import net.maizegenetics.util.LoggingUtils;

import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.Assert.*;

/**
 * Tests {@link SNPQualityProfilerPlugin} against a self-generated {@link GBSSimData} database that
 * has been run through the full discovery pipeline. Asserts that quality statistics are produced
 * for the (nonzero) set of simulated SNPs.
 *
 * @author Ed Buckler (original), rewritten for simulated data
 */
public class SNPQualityProfilerPluginTest {

    @Test
    public void testSNPQualityProfilerPlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        System.out.println("Running SNPQualityProfilerPlugin");

        GBSSimData sim = GBSSimData.createUnder("SNPQuality");
        sim.buildDatabaseThroughSam();
        sim.runDiscovery(true);

        Path statFile = sim.baseDir.resolve("snpQualityProfilerOutput.txt");
        new SNPQualityProfilerPlugin()
                .dBFile(sim.dbFile.toString())
                .taxaListName("ALL")
                .statFile(statFile.toString())
                .performFunction(null);

        assertTrue("Stat file should have been written", Files.exists(statFile));
        List<String> lines = Files.readAllLines(statFile);
        long dataRows = lines.stream().skip(1).filter(l -> !l.isEmpty()).count();
        System.out.println("SNPQualityProfilerPlugin produced " + dataRows + " stat rows");
        assertTrue("Should produce quality stats for at least one SNP", dataRows > 0);
    }
}
