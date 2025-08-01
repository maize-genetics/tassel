/*
 * SeqToTBTHDF5PluginTest
 */
package net.maizegenetics.analysis.gbs;

import java.io.File;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.tag.TBTTestUtils;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TaxaGroups;
import org.junit.Test;

/**
 *
 * @author terry
 */
public class SeqToTBTHDF5PluginTest {

    public SeqToTBTHDF5PluginTest() {
    }

    /**
     * Test of performFunction method, of class SeqToTBTHDF5Plugin.
     */
    @Test
    public void testPerformFunction() {

        if (!(new File(GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_DIR)).mkdirs()) {
            throw new IllegalStateException("SeqToTBTHDF5PluginTest: testPerformFunction: Can't create output directory: " + GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_DIR);
        }

        String[] inputArgs = new String[]{
            "-k", GBSConstants.GBS_TESTING_KEY_FILE,
            "-e", "ApeKI",
            "-i", GBSConstants.GBS_INPUT_DIR,
            "-o", GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_FILE,
            "-s", "100000000",
            "-L", GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_DIR + "SeqToTBTHDF5Plugin.log",
            "-t", GBSConstants.GBS_EXPECTED_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_FILE
        };

        SeqToTBTHDF5Plugin plugin = new SeqToTBTHDF5Plugin();
        plugin.setParameters(inputArgs);
        plugin.performFunction(null);

        TagsByTaxaByteHDF5TaxaGroups expectedTBT = new TagsByTaxaByteHDF5TaxaGroups(GBSConstants.GBS_EXPECTED_SEQ_TO_TBT_HDF5_PLUGIN_FILE);
        TagsByTaxaByteHDF5TaxaGroups actualTBT = new TagsByTaxaByteHDF5TaxaGroups(GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_FILE);

        System.out.println("SeqToTBTHDF5PluginTest: Comparing TBT Expected: " + GBSConstants.GBS_EXPECTED_SEQ_TO_TBT_HDF5_PLUGIN_FILE + "  Actual: " + GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_FILE);
        TBTTestUtils.compareTBTs(expectedTBT, actualTBT);
    }
}