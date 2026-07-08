package net.maizegenetics.analysis.gbs.v2;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.util.LoggingUtils;
import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.Assert.*;

/**
 * Property-based test for {@link ProductionSNPCallerPluginV2} that builds its own tiny GBSv2
 * database from {@link GBSSimData} (deterministic simulated FASTQ/key/reference) rather than
 * relying on a pre-built {@code GBSv2.db} or downloaded fixtures.
 *
 * <p>The production caller is run twice against the same self-generated database, once writing a
 * VCF and once writing an HDF5 ({@code .h5}) genotypes file. The HDF5 write path is the point of
 * this test: it exercises the jhdf5 native library (bundled in {@code cisd:jhdf5} for the CI/dev
 * platforms) without any external aligner or download. Both outputs are read back and asserted to
 * agree on site/taxon counts instead of comparing golden hashes.</p>
 */
public class ProductionSNPCallerPluginV2Test {

    @Test
    public void testProcessData() throws Exception {
        LoggingUtils.setupDebugLogging();

        GBSSimData sim = GBSSimData.createUnder("ProductionV2");
        sim.buildDatabaseThroughSam();   // tag DB + synthesized SAM alignments
        sim.runDiscovery(true);          // populate SNP positions (required by the production caller)

        Path vcfOut = sim.baseDir.resolve("prod.vcf");
        Path h5Out = sim.baseDir.resolve("prod.h5");

        for (Path out : new Path[]{vcfOut, h5Out}) {
            Files.deleteIfExists(out);
            new ProductionSNPCallerPluginV2()
                    .enzyme(GBSSimData.ENZYME)
                    .inputDirectory(sim.fastqDir.toString())
                    .inputGBSDatabase(sim.dbFile.toString())
                    .keyFile(sim.keyFile.toString())
                    .outputGenotypesFile(out.toString())
                    .kmerLength(GBSSimData.TAG_LENGTH)
                    .minimumQualityScore(0)
                    .performFunction(null);
            assertTrue("Production caller did not produce " + out, Files.exists(out));
        }

        GenotypeTable vcf = ImportUtils.readGuessFormat(vcfOut.toString());
        GenotypeTable h5 = ImportUtils.readGuessFormat(h5Out.toString());

        assertNotNull("Failed to read back VCF genotypes", vcf);
        assertNotNull("Failed to read back HDF5 genotypes", h5);

        assertTrue("HDF5 genotypes should contain at least one site", h5.numberOfSites() > 0);
        assertEquals("VCF and HDF5 outputs should call the same number of sites",
                vcf.numberOfSites(), h5.numberOfSites());
        assertEquals("VCF and HDF5 outputs should contain the same number of taxa",
                vcf.numberOfTaxa(), h5.numberOfTaxa());
        assertEquals("HDF5 taxa count should match the simulated taxa",
                sim.taxa.size(), h5.numberOfTaxa());
    }
}
