/*
 * MergeMultipleTagCountPluginTest
 */
package net.maizegenetics.analysis.gbs;

import java.io.File;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagCountsTestUtils;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.util.CheckSum;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

/**
 *
 * @author terry
 */
public class MergeMultipleTagCountPluginTest {

    public MergeMultipleTagCountPluginTest() {
    }

    /**
     * Test of performFunction method, of class MergeMultipleTagCountPlugin.
     */
    @Test
    public void testPerformFunction() {

        if (!(new File(GBSConstants.GBS_TEMP_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_DIR)).mkdirs()) {
            throw new IllegalStateException("MergeMultipleTagCountPluginTest: testPerformFunction: Can't create output directory: " + GBSConstants.GBS_TEMP_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_DIR);
        }

        String[] inputArgs = new String[]{
            "-i", GBSConstants.GBS_EXPECTED_FASTQ_TO_TAG_COUNT_PLUGIN_DIR,
            "-o", GBSConstants.GBS_TEMP_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_FILE,
            "-c", "5"
        };

        MergeMultipleTagCountPlugin plugin = new MergeMultipleTagCountPlugin();
        plugin.setParameters(inputArgs);
        plugin.performFunction(null);


        System.out.println("Verifying TagCounts: " + GBSConstants.GBS_TEMP_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_FILE);

        TagCounts tc = new TagCounts(GBSConstants.GBS_TEMP_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_FILE, TagsByTaxa.FilePacking.Byte);

        TagCountsTestUtils.sanityCheck(tc, -1);

        String actualMD5 = CheckSum.getMD5Checksum(GBSConstants.GBS_TEMP_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_FILE);
        String expectedMD5 = CheckSum.getMD5Checksum(GBSConstants.GBS_EXPECTED_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_FILE);

        System.out.println("Expected: " + GBSConstants.GBS_EXPECTED_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_FILE + ": " + expectedMD5);
        System.out.println("Actual: " + GBSConstants.GBS_TEMP_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_FILE + ": " + actualMD5);

        assertEquals("TagCountToFastqPluginTest Result MD5: " + GBSConstants.GBS_TEMP_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_FILE, expectedMD5, actualMD5);
    }
}