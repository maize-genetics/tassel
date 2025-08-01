/*
 * SAMConverterPluginTest
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
public class SAMConverterPluginTest {

    public SAMConverterPluginTest() {
    }

    /**
     * Test of performFunction method, of class SAMConverterPlugin.
     */
    @Test
    public void testPerformFunction() {

        if (!(new File(GBSConstants.GBS_TEMP_SAM_CONVERTER_PLUGIN_DIR)).mkdirs()) {
            throw new IllegalStateException("SAMConverterPluginTest: testPerformFunction: Can't create output directory: " + GBSConstants.GBS_TEMP_SAM_CONVERTER_PLUGIN_DIR);
        }

        String[] inputArgs = new String[]{
            "-i", GBSConstants.GBS_EXPECTED_BOWTIE_SAM_FILE,
            "-o", GBSConstants.GBS_TEMP_SAM_CONVERTER_PLUGIN_FILE
        };

        SAMConverterPlugin plugin = new SAMConverterPlugin();
        plugin.setParameters(inputArgs);
        plugin.performFunction(null);

        String actualMD5 = CheckSum.getMD5Checksum(GBSConstants.GBS_TEMP_SAM_CONVERTER_PLUGIN_FILE);
        String expectedMD5 = CheckSum.getMD5Checksum(GBSConstants.GBS_EXPECTED_SAM_CONVERTER_PLUGIN_FILE);

        System.out.println("Expected: " + GBSConstants.GBS_EXPECTED_SAM_CONVERTER_PLUGIN_FILE + ": " + expectedMD5);
        System.out.println("Actual: " + GBSConstants.GBS_TEMP_SAM_CONVERTER_PLUGIN_FILE + ": " + actualMD5);

        assertEquals("SAMConverterPluginTest Result MD5: " + GBSConstants.GBS_TEMP_SAM_CONVERTER_PLUGIN_FILE, expectedMD5, actualMD5);

    }
}