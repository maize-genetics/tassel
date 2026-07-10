/*
 * MergeMultipleTagCountPluginTest
 */
package net.maizegenetics.analysis.gbs;

import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagCountsTestUtils;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Property-based rehabilitation of the legacy GBSv1 {@code MergeMultipleTagCountPlugin} test.
 *
 * <p>The original test MD5-compared the merged master {@code .cnt} against a downloaded golden
 * fixture. Instead this test self-generates a tiny deterministic data set with {@link GBSv1SimData}
 * ({@code FastqToTagCountPlugin} then {@code MergeMultipleTagCountPlugin}) and asserts that the merged
 * master tag list passes the shared {@link TagCountsTestUtils#sanityCheck} and contains exactly the
 * number of distinct tags the simulator injected.</p>
 *
 * @author terry (original), rehabilitated for self-generated data
 */
public class MergeMultipleTagCountPluginTest {

    @Test
    public void testPerformFunction() throws Exception {
        GBSv1SimData sim = GBSv1SimData.createUnder("MergeMultipleTagCount");
        sim.buildMasterTagCounts();

        TagCounts masterTags = new TagCounts(sim.masterTagCounts.toString(), FilePacking.Byte);
        TagCountsTestUtils.sanityCheck(masterTags, -1);

        assertTrue("Merged master tag list should be non-empty", masterTags.getTagCount() > 0);
        assertEquals("Merged master tag list should contain every distinct simulated tag",
                sim.sim.expectedDistinctTags, masterTags.getTagCount());
    }
}
