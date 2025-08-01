/*
 * TagCountToFastqPluginTest
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
public class TagCountToFastqPluginTest {

    public TagCountToFastqPluginTest() {
    }

    /**
     * Test of performFunction method, of class TagExportToFastqPlugin.
     */
    @Test
    public void testPerformFunction() {

        if (!(new File(GBSConstants.GBS_TEMP_TAG_COUNT_TO_FASTQ_PLUGIN_DIR)).mkdirs()) {
            throw new IllegalStateException("TagCountToFastqPluginTest: testPerformFunction: Can't create output directory: " + GBSConstants.GBS_TEMP_TAG_COUNT_TO_FASTQ_PLUGIN_DIR);
        }

        String[] inputArgs = new String[]{
            "-i", GBSConstants.GBS_EXPECTED_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_FILE,
            "-o", GBSConstants.GBS_TEMP_TAG_COUNT_TO_FASTQ_PLUGIN_FILE,
            "-c", "5"
        };

        TagCountToFastqPlugin plugin = new TagCountToFastqPlugin();
        plugin.setParameters(inputArgs);
        plugin.performFunction(null);

        String actualMD5 = CheckSum.getMD5Checksum(GBSConstants.GBS_TEMP_TAG_COUNT_TO_FASTQ_PLUGIN_FILE);
        String expectedMD5 = CheckSum.getMD5Checksum(GBSConstants.GBS_EXPECTED_TAG_COUNT_TO_FASTQ_PLUGIN_FILE);

        System.out.println("Expected: " + GBSConstants.GBS_EXPECTED_TAG_COUNT_TO_FASTQ_PLUGIN_FILE + ": " + expectedMD5);
        System.out.println("Actual: " + GBSConstants.GBS_TEMP_TAG_COUNT_TO_FASTQ_PLUGIN_FILE + ": " + actualMD5);

        assertEquals("TagCountToFastqPluginTest Result MD5: " + GBSConstants.GBS_TEMP_TAG_COUNT_TO_FASTQ_PLUGIN_FILE, expectedMD5, actualMD5);

    }
}