package net.maizegenetics.analysis.gbs.v2;

import java.util.Map;
import java.util.Set;

import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TagDataSQLite;

import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

/**
 * Property-based rehabilitation of the legacy TAS-865 GBSv1-vs-GBSv2 tag comparison.
 *
 * <p>The original class contained no JUnit assertions: every method merely ran the v2 pipeline against
 * large downloaded fixtures and wrote diff/text files for a human to inspect. This rewrite preserves
 * the class's intent (verify the biology of the tags the GBSv2 pipeline stores) with real assertions
 * on a deterministic, self-generated data set from {@link GBSSimData}: the tag DB contains exactly the
 * simulated distinct tags, each stored tag is one of the known simulated sequences, and every
 * simulated tag is present with a positive depth.</p>
 */
public class GBSv2BiologyCompareTest {

    @Test
    public void compareTagsTest() throws Exception {
        GBSSimData sim = GBSSimData.createUnder("BiologyCompare");
        sim.buildTagDB();

        TagDataSQLite tdw = new TagDataSQLite(sim.dbFile.toString());
        try {
            Set<Tag> dbTags = tdw.getTags();
            assertEquals("GBSv2 DB should store exactly the distinct simulated tags",
                    sim.expectedDistinctTags, dbTags.size());

            for (GBSSimData.TagInfo info : sim.tagInfos) {
                Tag simTag = TagBuilder.instance(info.sequence).build();
                assertTrue("GBSv2 DB should contain simulated tag " + info.sequence,
                        dbTags.contains(simTag));
            }

            Map<Tag, Integer> tagsWithDepth = tdw.getTagsWithDepth(0);
            for (GBSSimData.TagInfo info : sim.tagInfos) {
                Tag simTag = TagBuilder.instance(info.sequence).build();
                Integer depth = tagsWithDepth.get(simTag);
                assertNotNull("Simulated tag should have a stored depth: " + info.sequence, depth);
                assertTrue("Simulated tag depth should be positive: " + info.sequence, depth > 0);
            }
        } finally {
            tdw.close();
        }
    }
}
