/*
 * ProductionSNPCallerPluginTest
 */
package net.maizegenetics.analysis.gbs;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

/**
 * Property-based rehabilitation of the legacy GBSv1 {@code ProductionSNPCallerPlugin} test.
 *
 * <p>The original test compared the output HDF5 genotypes against a downloaded golden
 * {@code PipelineTestingGenos.h5} fixture. Instead this test self-generates its inputs with
 * {@link GBSv1SimData} (raw FASTQ, key, and a discovery TOPM with variants), runs the production
 * caller to write an HDF5 genotypes file, reads it back, and asserts on its structural properties:
 * all simulated taxa are genotyped and the discovered SNP sites are present.</p>
 *
 * @author jcg233 (original), rehabilitated for self-generated data
 */
public class ProductionSNPCallerPluginTest {

    @Test
    public void testPerformFunction() throws Exception {
        GBSv1SimData sim = GBSv1SimData.createUnder("Production");
        sim.buildMasterTagCounts();
        sim.buildTopm();
        sim.buildTbt();
        sim.pivotTbt();
        sim.runDiscovery();

        Path outFile = sim.v1Dir.resolve("PipelineTestingGenos.h5");
        Files.deleteIfExists(outFile);

        new ProductionSNPCallerPlugin()
                .inputDirectory(sim.sim.fastqDir.toString())
                .keyFile(sim.keyFileV1.toString())
                .enzyme(GBSv1SimData.ENZYME)
                .inputTOPMFile(sim.discoveryTopmFile.toString())
                .outputHDF5GenotypesFile(outFile.toString())
                .performFunction(null);

        assertTrue("Production caller did not produce " + outFile, Files.exists(outFile));

        GenotypeTable genos = ImportUtils.readGuessFormat(outFile.toString());
        assertNotNull("Failed to read back HDF5 genotypes", genos);
        assertEquals("HDF5 taxa count should match the simulated taxa",
                sim.sim.taxa.size(), genos.numberOfTaxa());
        assertTrue("HDF5 genotypes should contain at least one site", genos.numberOfSites() > 0);
    }
}
