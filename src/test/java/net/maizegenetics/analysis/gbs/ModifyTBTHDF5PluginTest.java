/*
 * ModifyTBTHDF5PluginTest
 */
package net.maizegenetics.analysis.gbs;

import java.io.File;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.tag.TBTTestUtils;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TaxaGroups;
import org.junit.Test;

/**
 *
 * @author terry
 */
public class ModifyTBTHDF5PluginTest {

    public ModifyTBTHDF5PluginTest() {
    }

    /**
     * Test of -p option (Pivot / Transpose) of class ModifyTBTHDF5Plugin.
     */
    @Test
    public void testPerformFunctionPOption() {

        (new File(GBSConstants.GBS_TEMP_MODIFY_TBT_HDF5_PLUGIN_DIR)).mkdirs();

        String[] inputArgs = new String[]{
            "-o", GBSConstants.GBS_EXPECTED_SEQ_TO_TBT_HDF5_PLUGIN_FILE,
            "-p", GBSConstants.GBS_TEMP_MODIFY_TBT_HDF5_PLUGIN_PIVOTED_FILE
        };

        ModifyTBTHDF5Plugin plugin = new ModifyTBTHDF5Plugin();
        plugin.setParameters(inputArgs);
        plugin.performFunction(null);

        System.out.println("ModifyTBTHDF5PluginTest: Comparing TBT Expected: " + GBSConstants.GBS_EXPECTED_SEQ_TO_TBT_HDF5_PLUGIN_FILE + "  Actual: " + GBSConstants.GBS_TEMP_MODIFY_TBT_HDF5_PLUGIN_PIVOTED_FILE);

        TagsByTaxa expectedTBT = new TagsByTaxaByteHDF5TaxaGroups(GBSConstants.GBS_EXPECTED_SEQ_TO_TBT_HDF5_PLUGIN_FILE);
        TagsByTaxa actualTBT = new TagsByTaxaByteHDF5TagGroups(GBSConstants.GBS_TEMP_MODIFY_TBT_HDF5_PLUGIN_PIVOTED_FILE);

        TBTTestUtils.compareTBTs(expectedTBT, actualTBT);
    }

    /**
     * Test of -i option (Merge TBT Files) of class ModifyTBTHDF5Plugin.
     */
    public void testPerformFunctionIOption() {

        (new File(GBSConstants.GBS_TEMP_MODIFY_TBT_HDF5_PLUGIN_DIR)).mkdirs();

        String[] inputArgs = new String[]{
            "-o", GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_FILE,
            "-i", GBSConstants.GBS_EXPECTED_SEQ_TO_TBT_HDF5_PLUGIN__MIRROR_FILE
        };

        ModifyTBTHDF5Plugin plugin = new ModifyTBTHDF5Plugin();
        plugin.setParameters(inputArgs);
        plugin.performFunction(null);

        System.out.println("ModifyTBTHDF5PluginTest: Comparing TBT Expected: " + GBSConstants.GBS_EXPECTED_MODIFY_TBT_HDF5_PLUGIN_MERGE_FILE + "  Actual: " + GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_FILE);

        TagsByTaxa expectedTBT = new TagsByTaxaByteHDF5TaxaGroups(GBSConstants.GBS_EXPECTED_MODIFY_TBT_HDF5_PLUGIN_MERGE_FILE);
        TagsByTaxa actualTBT = new TagsByTaxaByteHDF5TaxaGroups(GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_FILE);

        TBTTestUtils.compareTBTs(expectedTBT, actualTBT);
    }

    /**
     * Test of -c option (Merge taxa with same LibraryPrepID) of class
     * ModifyTBTHDF5Plugin.
     */
    public void testPerformFunctionCOption() {

        (new File(GBSConstants.GBS_TEMP_MODIFY_TBT_HDF5_PLUGIN_DIR)).mkdirs();

        String[] inputArgs = new String[]{
            "-o", GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_FILE,
            "-c"
        };

        ModifyTBTHDF5Plugin plugin = new ModifyTBTHDF5Plugin();
        plugin.setParameters(inputArgs);
        plugin.performFunction(null);

        System.out.println("ModifyTBTHDF5PluginTest: Comparing TBT Expected: " + GBSConstants.GBS_EXPECTED_MODIFY_TBT_HDF5_PLUGIN_MERGE_TAXA_FILE + "  Actual: " + GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_FILE);

        TagsByTaxa expectedTBT = new TagsByTaxaByteHDF5TaxaGroups(GBSConstants.GBS_EXPECTED_MODIFY_TBT_HDF5_PLUGIN_MERGE_TAXA_FILE);
        TagsByTaxa actualTBT = new TagsByTaxaByteHDF5TaxaGroups(GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_FILE);

        TBTTestUtils.compareTBTs(expectedTBT, actualTBT);
    }
}
