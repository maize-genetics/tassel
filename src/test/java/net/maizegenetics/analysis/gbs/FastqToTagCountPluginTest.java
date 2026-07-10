/*
 * FastqToTagCountPluginTest
 */
package net.maizegenetics.analysis.gbs;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagCountsTestUtils;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Property-based rehabilitation of the legacy GBSv1 {@code FastqToTagCountPlugin} test.
 *
 * <p>The original test compared each produced {@code .cnt} against a downloaded byte-exact golden
 * fixture (MD5). Instead this test self-generates a tiny deterministic data set with
 * {@link GBSv1SimData} (which reuses the aligner-free {@code GBSSimData} FASTQ/key), runs the plugin,
 * and asserts on structural properties of the resulting TagCounts: the file is produced, passes the
 * shared {@link TagCountsTestUtils#sanityCheck} (no duplicate/padding violations, positive counts,
 * accumulated count within the raw-read bound), and contains a non-empty tag list.</p>
 *
 * @author terry (original), rehabilitated for self-generated data
 */
public class FastqToTagCountPluginTest {

    @Test
    public void testPerformFunction() throws Exception {
        GBSv1SimData sim = GBSv1SimData.createUnder("FastqToTagCount");
        Path outputDir = sim.v1Dir.resolve("fastqToTagCount");
        Files.createDirectories(outputDir);

        FastqToTagCountPlugin plugin = new FastqToTagCountPlugin();
        plugin.setParameters(new String[]{
                "-i", sim.sim.fastqDir.toString(),
                "-o", outputDir.toString(),
                "-k", sim.sim.keyFile.toString(),
                "-e", GBSv1SimData.ENZYME,
                "-s", "150000000",
                "-c", "1"
        });
        plugin.performFunction(null);

        String fastqRegex = "(?i).*\\.fq$|.*\\.fq\\.gz$|.*\\.fastq$|.*_fastq\\.txt$|.*_fastq\\.gz$|.*_fastq\\.txt\\.gz$|.*_sequence\\.txt$|.*_sequence\\.txt\\.gz$";
        File[] fastqFiles = DirectoryCrawler.listFiles(fastqRegex, sim.sim.fastqDir.toString());
        assertTrue("Expected at least one simulated FASTQ file", fastqFiles.length > 0);

        int totalTags = 0;
        for (File fastq : fastqFiles) {
            String cntName = fastq.getName().replaceAll("(?i)\\.fq$|\\.fq\\.gz$|\\.fastq$|_fastq\\.txt$|_fastq\\.gz$|_fastq\\.txt\\.gz$|_sequence\\.txt$|_sequence\\.txt\\.gz$", ".cnt");
            File cntFile = new File(outputDir.toFile(), cntName);
            assertTrue("Expected TagCounts output file: " + cntFile, cntFile.exists());

            int numLines = Utils.getNumberLines(fastq.getAbsolutePath());
            assertEquals("Input FASTQ should have a multiple of 4 lines: " + numLines, 0, numLines % 4);
            int numRawSeqs = numLines / 4;

            TagCounts tc = new TagCounts(cntFile.toString(), FilePacking.Byte);
            TagCountsTestUtils.sanityCheck(tc, numRawSeqs);
            totalTags += tc.getTagCount();
        }

        assertTrue("FastqToTagCountPlugin should produce a non-empty tag list", totalTags > 0);
    }
}
