/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.analysis.data;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;
import net.maizegenetics.analysis.gbs.ProductionSNPCallerPlugin;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.Utils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author jcg233
 */
public class BuildUnfinishedHDF5GenotypesPluginTest {

    String tempUnfinishedFile, tempFinishedFile;
    
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

        // run production SNP caller with the -ko option to generate an appropriate input file
        ProductionSNPCallerPlugin plugin = new ProductionSNPCallerPlugin()
            .inputDirectory(GBSConstants.GBS_INPUT_DIR)
            .keyFile(GBSConstants.GBS_TESTING_KEY_FILE)
            .enzyme("ApeKI")
            .inputTOPMFile(GBSConstants.GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE)
            .outputHDF5GenotypesFile(tempUnfinishedFile)
            .keepGenotypesOpen(true)
        ;
        plugin.performFunction(null);
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
        String date = new SimpleDateFormat("yyyyMMdd").format(new Date());
        String expectedName = "ReallyBigData_" + date;
        assertEquals("The DataSetName differs from expected",expectedName,FinishedGenos.annotations().getTextAnnotation(GenotypeTable.ANNOTATION_DATA_SET_NAME)[0]);
        String expectedDescrip = "This really big data set from "+date+" has 14754 SNPs and 192 samples in it!";
        assertEquals("The DataSetDescription differs from expected",expectedDescrip,FinishedGenos.annotations().getTextAnnotation(GenotypeTable.ANNOTATION_DATA_SET_DESCRIPTION)[0]);
    }
    
}
