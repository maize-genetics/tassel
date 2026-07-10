/*
 * SAMConverterPluginTest
 */
package net.maizegenetics.analysis.gbs;

import java.util.HashSet;
import java.util.Set;

import net.maizegenetics.analysis.gbs.v2.GBSSimData;
import net.maizegenetics.dna.map.TagsOnPhysicalMap;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Property-based rehabilitation of the legacy GBSv1 {@code SAMConverterPlugin} test.
 *
 * <p>The original test needed a pre-made bowtie SAM and MD5-compared the resulting {@code .topm}
 * against a downloaded golden fixture. Instead this test lets {@link GBSv1SimData} synthesize a
 * perfect-alignment SAM directly from the known tag/reference positions (no external aligner) and run
 * {@code SAMConverterPlugin}. It then asserts every tag was uniquely aligned to its known
 * chromosome/position, dropping the byte-exact hash comparison.</p>
 *
 * @author terry (original), rehabilitated for self-generated data
 */
public class SAMConverterPluginTest {

    @Test
    public void testPerformFunction() throws Exception {
        GBSv1SimData sim = GBSv1SimData.createUnder("SAMConverter");
        sim.buildTagCountsAndTopm();

        Set<String> knownPositions = new HashSet<>();
        for (GBSSimData.TagInfo info : sim.sim.tagInfos) {
            knownPositions.add(info.chrom + ":" + info.cutPosition);
        }

        TagsOnPhysicalMap topm = new TagsOnPhysicalMap(sim.topmFile.toString(), true);
        assertEquals("TOPM should contain one entry per synthesized SAM tag",
                sim.sim.tagInfos.size(), topm.getSize());

        for (int tag = 0; tag < topm.getSize(); tag++) {
            int chromosome = topm.getChromosome(tag);
            assertTrue("Every simulated tag should be uniquely aligned (have a chromosome)",
                    chromosome != Integer.MIN_VALUE);
            int startPosition = topm.getStartPosition(tag);
            String key = chromosome + ":" + startPosition;
            assertTrue("Aligned position " + key + " should be a known simulated cut site",
                    knownPositions.contains(key));
        }
    }
}
