/*
 * ModifyTBTHDF5PluginTest
 */
package net.maizegenetics.analysis.gbs;

import net.maizegenetics.dna.tag.TBTTestUtils;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TaxaGroups;
import org.junit.Test;

/**
 * Property-based rehabilitation of the legacy GBSv1 {@code ModifyTBTHDF5Plugin} test.
 *
 * <p>The pivot ({@code -p}) option transposes a taxa-oriented TBT into a tag-oriented TBT; the two
 * orientations must hold identical data. Rather than comparing to a downloaded golden
 * {@code TBT_Pivoted.h5} fixture, this test builds a TBT from self-generated data, pivots it, and
 * asserts the round-trip is loss-free with {@link TBTTestUtils#compareTBTs}.</p>
 *
 * @author ed (original), rehabilitated for self-generated data
 */
public class ModifyTBTHDF5PluginTest {

    @Test
    public void testPerformFunctionPOption() throws Exception {
        GBSv1SimData sim = GBSv1SimData.createUnder("ModifyTBT");
        sim.buildMasterTagCounts();
        sim.buildTbt();
        sim.pivotTbt();

        TagsByTaxaByteHDF5TaxaGroups expectedTBT = new TagsByTaxaByteHDF5TaxaGroups(sim.tbtFile.toString());
        TagsByTaxaByteHDF5TagGroups actualTBT = new TagsByTaxaByteHDF5TagGroups(sim.pivotedTbtFile.toString());
        try {
            TBTTestUtils.compareTBTs(expectedTBT, actualTBT);
        } finally {
            expectedTBT.closeWriter();
            actualTBT.closeWriter();
        }
    }
}
