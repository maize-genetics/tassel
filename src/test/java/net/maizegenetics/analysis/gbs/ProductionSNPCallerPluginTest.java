/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.analysis.gbs;

import java.io.File;
import java.util.Iterator;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.AlignmentTestingUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author jcg233
 */
public class ProductionSNPCallerPluginTest {

    GenotypeTable ExpectedGenos;
    
    public ProductionSNPCallerPluginTest() {
    }
    
    @Before
    public void setUp() {
        ExpectedGenos = ImportUtils.readGuessFormat(GBSConstants.GBS_EXPECTED_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_OUT_FILE);
    }

    /**
     * Test of performFunction method, of class ProductionSNPCallerPlugin.
     */
    @Test
    public void testPerformFunction() {
        if (!(new File(GBSConstants.GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_DIR)).mkdirs()) {
            throw new IllegalStateException("ProductionSNPCallerPluginTest: testPerformFunction: Can't create output directory: " + GBSConstants.GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_DIR);
        }
        ProductionSNPCallerPlugin plugin = new ProductionSNPCallerPlugin()
                .inputDirectory(GBSConstants.GBS_INPUT_DIR)
                .keyFile(GBSConstants.GBS_TESTING_KEY_FILE)
                .enzyme("ApeKI")
                .inputTOPMFile(GBSConstants.GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE)
                .outputHDF5GenotypesFile(GBSConstants.GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_OUT_FILE);
        plugin.performFunction(null);
        
        // open actual HDF5Genos and compare to expected
        GenotypeTable ActualGenos = ImportUtils.readGuessFormat(GBSConstants.GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_OUT_FILE);
        AlignmentTestingUtils.alignmentsEqual(ExpectedGenos, ActualGenos, true);
        
        // compare the Postions (not really needed)
        PositionList ExpectedPosits = ExpectedGenos.positions();
        PositionList ActualPosits = ActualGenos.positions();
        assertEquals("Expected and Actual PositionLists are not the same size", ExpectedPosits.size(), ActualPosits.size());
        Iterator<Position> ExpIter = ExpectedPosits.iterator();
        Iterator<Position> ActIter = ActualPosits.iterator();
        while (ExpIter.hasNext()) {
            Position expPos = ExpIter.next();
            Position actPos = ActIter.next();
            assertEquals("Expected and Actual Positions are not the same for postion "+expPos.getSNPID()+" vs "+actPos.getSNPID(), expPos.compareTo(actPos), 0);
        }
        
        // compare the depths, allowing taxa to be in a different sort order
        int[] t1to2 = new int[ActualGenos.numberOfTaxa()];
        for (int t1 = 0; t1 < ActualGenos.numberOfTaxa(); t1++) {
            t1to2[t1] = ActualGenos.taxa().indexOf(ExpectedGenos.taxa().get(t1));
        }
        for (int t = 0; t < ExpectedGenos.numberOfTaxa(); t++) {
            for (int s = 0; s < ExpectedGenos.numberOfSites(); s++) {
                int[] expDepths = ExpectedGenos.depthForAlleles(t, s);
                int[] actDepths = ActualGenos.depthForAlleles(t1to2[t], s);
                assertEquals(
                    "Number of alleles in depth array differs for taxon "+ExpectedGenos.taxaName(t)+" and site "+ExpectedGenos.siteName(s),expDepths.length,actDepths.length
                );
                for (int a = 0; a < expDepths.length; a++) {
                    assertEquals("Depth at allele "+a+" differs for taxon "+ExpectedGenos.taxaName(t)+" and site "+ExpectedGenos.siteName(s),expDepths[a],actDepths[a]);
                }
            }
        }
    }
}
