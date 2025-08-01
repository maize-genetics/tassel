/*
 * DiscoverySNPCallerPluginTest
 */
package net.maizegenetics.analysis.gbs;

import java.io.File;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.util.CheckSum;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

/**
 *
 * @author terry
 */
public class DiscoverySNPCallerPluginTest {

    public DiscoverySNPCallerPluginTest() {
    }

    /**
     * Test of performFunction method, of class DiscoverySNPCallerPlugin.
     */
    @Test
    public void testPerformFunction() {

        if (!(new File(GBSConstants.GBS_TEMP_DISCOVERY_SNP_CALLER_PLUGIN_DIR)).mkdirs()) {
            throw new IllegalStateException("DiscoverySNPCallerPluginTest: testPerformFunction: Can't create output directory: " + GBSConstants.GBS_TEMP_DISCOVERY_SNP_CALLER_PLUGIN_DIR);
        }

        String[] inputArgs = new String[]{
            "-i", GBSConstants.GBS_EXPECTED_MODIFY_TBT_HDF5_PLUGIN_PIVOTED_FILE,
            "-m", GBSConstants.GBS_EXPECTED_SAM_CONVERTER_PLUGIN_FILE,
            "-o", GBSConstants.GBS_TEMP_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE,
            "-mnF", "0.8",
            "-mnMAF", "0.02",
            "-mnMAC", "100000",
            "-ref", GBSConstants.GBS_REFERENCE_GENOME,
            "-sC", "9",
            "-eC", "10"
        };

        DiscoverySNPCallerPlugin plugin = new DiscoverySNPCallerPlugin();
        plugin.setParameters(inputArgs);
        plugin.performFunction(null);

        String expectedMD5 = CheckSum.getMD5Checksum(GBSConstants.GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE);
        String actualMD5 = CheckSum.getMD5Checksum(GBSConstants.GBS_TEMP_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE);

        System.out.println("Expected: " + GBSConstants.GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE + ":\n  " + expectedMD5);
        System.out.println("Actual: " + GBSConstants.GBS_TEMP_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE + ":\n  " + actualMD5);

        assertEquals("DiscoverySNPCallerPluginTest Result MD5: " + GBSConstants.GBS_TEMP_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE, expectedMD5, actualMD5);
    }
}