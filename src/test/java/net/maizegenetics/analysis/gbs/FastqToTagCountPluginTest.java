/*
 * FastqToTagCountPluginTest
 */
package net.maizegenetics.analysis.gbs;

import java.io.File;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagCountsTestUtils;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.util.CheckSum;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;
import org.junit.Test;
import static org.junit.Assert.assertEquals;

/**
 *
 * @author terry
 */
public class FastqToTagCountPluginTest {

    public FastqToTagCountPluginTest() {
    }

    /**
     * Test of performFunction method, of class FastqToTagCountPlugin.
     */
    @Test
    public void testPerformFunction() {

        if (!(new File(GBSConstants.GBS_TEMP_FASTQ_TO_TAG_COUNT_PLUGIN_DIR)).mkdirs()) {
            throw new IllegalStateException("FastqToTagCountPluginTest: testPerformFunction: Can't create output directory: " + GBSConstants.GBS_TEMP_FASTQ_TO_TAG_COUNT_PLUGIN_DIR);
        }

        String[] inputArgs = new String[]{
            "-i", GBSConstants.GBS_INPUT_DIR,
            "-o", GBSConstants.GBS_TEMP_FASTQ_TO_TAG_COUNT_PLUGIN_DIR,
            "-k", GBSConstants.GBS_TESTING_KEY_FILE,
            "-e", "ApeKI",
            "-s", "150000000",
            "-c", "1"
        };

        FastqToTagCountPlugin plugin = new FastqToTagCountPlugin();
        plugin.setParameters(inputArgs);
        plugin.performFunction(null);

        File inputDirectory = new File(GBSConstants.GBS_INPUT_DIR);
        File[] fastqFiles = DirectoryCrawler.listFiles("(?i).*\\.fq$|.*\\.fq\\.gz$|.*\\.fastq$|.*_fastq\\.txt$|.*_fastq\\.gz$|.*_fastq\\.txt\\.gz$|.*_sequence\\.txt$|.*_sequence\\.txt\\.gz$", inputDirectory.getAbsolutePath());

        for (int i = 0; i < fastqFiles.length; i++) {

            String tagCountFilename = fastqFiles[i].getName().replaceAll("(?i)\\.fq$|\\.fq\\.gz$|\\.fastq$|_fastq\\.txt$|_fastq\\.gz$|_fastq\\.txt\\.gz$|_sequence\\.txt$|_sequence\\.txt\\.gz$", ".cnt");
            String expectedFilename = GBSConstants.GBS_EXPECTED_FASTQ_TO_TAG_COUNT_PLUGIN_DIR + tagCountFilename;
            tagCountFilename = GBSConstants.GBS_TEMP_FASTQ_TO_TAG_COUNT_PLUGIN_DIR + tagCountFilename;
            System.out.println("Verifying TagCounts: " + tagCountFilename);

            int numLines = Utils.getNumberLines(fastqFiles[i].getAbsolutePath());
            assertEquals("Input Fastq File should have Multiple of 4 Lines: " + numLines, 0, (numLines % 4));
            int numRawSeqs = numLines / 4;

            TagCounts tc = new TagCounts(tagCountFilename, TagsByTaxa.FilePacking.Byte);

            TagCountsTestUtils.sanityCheck(tc, numRawSeqs);

            String actualMD5 = CheckSum.getMD5Checksum(tagCountFilename);
            String expectedMD5 = CheckSum.getMD5Checksum(expectedFilename);

            System.out.println("Expected: " + expectedFilename + ": " + expectedMD5);
            System.out.println("Actual: " + tagCountFilename + ": " + actualMD5);

            assertEquals("TagCountToFastqPluginTest Result MD5: " + tagCountFilename, expectedMD5, actualMD5);
        }

    }
}