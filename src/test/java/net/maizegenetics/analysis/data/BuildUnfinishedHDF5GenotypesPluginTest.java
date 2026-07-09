/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.analysis.data;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.Utils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Verifies that {@link BuildUnfinishedHDF5GenotypesPlugin} finalizes an "unfinished"
 * HDF5 genotype file (one closed with {@code closeUnfinished()} rather than {@code build()}).
 *
 * <p>Historically the input was produced by running the legacy GBSv1
 * {@code ProductionSNPCallerPlugin} with {@code keepGenotypesOpen(true)}, which required a
 * multi-hundred-MB raw FASTQ input directory that is not part of the test-data release. The
 * plugin under test does not depend on the v1 pipeline at all: it simply reopens an unfinished
 * HDF5 builder and calls {@code build()}. So the unfinished input is now generated directly with
 * {@link GenotypeTableBuilder#getTaxaIncremental} + {@code closeUnfinished()} from a small,
 * already-present genotype fixture, and assertions are derived from that fixture rather than
 * hard-coded golden counts.</p>
 *
 * @author jcg233
 */
public class BuildUnfinishedHDF5GenotypesPluginTest {

    String tempUnfinishedFile, tempFinishedFile;
    private int expectedSites, expectedTaxa;

    public BuildUnfinishedHDF5GenotypesPluginTest() {
    }

    @Before
    public void setUp() {
        // mkdirs and delete any previous files
        String date = "_" + new SimpleDateFormat("yyyyMMdd").format(new Date());
        File outdir = new File(GBSConstants.GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_DIR);
        outdir.mkdirs();
        tempUnfinishedFile = Utils.getFilename(GBSConstants.GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_UNFINISHED_FILE)+".h5";
        tempFinishedFile = Utils.getFilename(GBSConstants.GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_FINISHED_FILE)+".h5";
        tempFinishedFile = (tempFinishedFile.replace("__DATE__", date));
        if(outdir.listFiles()!=null) {
            for (File file : outdir.listFiles()) {
                if(file.getName().equals(tempUnfinishedFile) || file.getName().equals(tempFinishedFile)) {
                    file.delete();
                }
            }
        }
        tempUnfinishedFile = GBSConstants.GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_UNFINISHED_FILE;
        tempFinishedFile = GBSConstants.GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_FINISHED_FILE;
        tempFinishedFile = (tempFinishedFile.replace("__DATE__", date));

        // Delete any stale files so the incremental builder and the plugin's copy start clean.
        new File(tempUnfinishedFile).delete();
        new File(tempFinishedFile).delete();

        // Build an "unfinished" HDF5 genotype file directly from a present genotype fixture,
        // mirroring what the v1 caller's keepGenotypesOpen(true) used to produce.
        GenotypeTable source = ImportUtils.readGuessFormat(GBSConstants.GBS_EXPECTED_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_OUT_FILE);
        expectedSites = source.numberOfSites();
        expectedTaxa = source.numberOfTaxa();
        GenotypeTableBuilder gtb = GenotypeTableBuilder.getTaxaIncremental(source.positions(), tempUnfinishedFile);
        for (int i = 0; i < source.numberOfTaxa(); i++) {
            gtb.addTaxon(source.taxa().get(i), source.genotypeAllSites(i));
        }
        gtb.closeUnfinished();

        // Mirror the "-ko" (keep-genotypes-open) state the plugin expects: the taxa module is
        // locked (so it can be read) while the genotype module stays open for finalization.
        IHDF5Writer writer = HDF5Factory.open(tempUnfinishedFile);
        HDF5Utils.lockHDF5TaxaModule(writer);
        writer.close();
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of processData method, of class BuildUnfinishedHDF5GenotypesPlugin.
     */
    @Test
    public void testProcessData() {
        BuildUnfinishedHDF5GenotypesPlugin plugin = new BuildUnfinishedHDF5GenotypesPlugin()
            .inputFile(tempUnfinishedFile)
            .outputFile(tempFinishedFile)
            .dataSetName("ReallyBigData__DATE__")
            .dataSetDescription("This really big data set from __DATE__ has __SNPS__ SNPs and __TAXA__ samples in it!")
        ;
        plugin.performFunction(null);

        IHDF5Reader h5Reader = HDF5Factory.openForReading(tempFinishedFile);
        assertTrue("Finshed file is not TaxaLocked", HDF5Utils.isTaxaLocked(h5Reader));
        assertTrue("Finished file is not GenotypeLocked", HDF5Utils.isHDF5GenotypeLocked(h5Reader));
        h5Reader.close();

        GenotypeTable FinishedGenos = ImportUtils.readGuessFormat(tempFinishedFile);
        assertEquals("Finished file has an unexpected number of sites", expectedSites, FinishedGenos.numberOfSites());
        assertEquals("Finished file has an unexpected number of taxa", expectedTaxa, FinishedGenos.numberOfTaxa());
        String date = new SimpleDateFormat("yyyyMMdd").format(new Date());
        String expectedName = "ReallyBigData_" + date;
        assertEquals("The DataSetName differs from expected",expectedName,FinishedGenos.annotations().getTextAnnotation(GenotypeTable.ANNOTATION_DATA_SET_NAME)[0]);
        String expectedDescrip = "This really big data set from "+date+" has "+expectedSites+" SNPs and "+expectedTaxa+" samples in it!";
        assertEquals("The DataSetDescription differs from expected",expectedDescrip,FinishedGenos.annotations().getTextAnnotation(GenotypeTable.ANNOTATION_DATA_SET_DESCRIPTION)[0]);
    }

}
