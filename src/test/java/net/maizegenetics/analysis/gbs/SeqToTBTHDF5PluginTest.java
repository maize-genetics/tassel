/*
 * SeqToTBTHDF5PluginTest
 */
package net.maizegenetics.analysis.gbs;

import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TaxaGroups;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Property-based rehabilitation of the legacy GBSv1 {@code SeqToTBTHDF5Plugin} test.
 *
 * <p>Instead of comparing against a downloaded golden {@code TBT_from_Raw_Seq.h5} fixture, this test
 * self-generates a tiny deterministic data set with {@link GBSv1SimData}, runs the plugin to build a
 * TagsByTaxa HDF5, and asserts on its structural properties: every simulated taxon is present, every
 * master tag is represented, and reads were actually matched into the TBT.</p>
 *
 * @author terry (original), rehabilitated for self-generated data
 */
public class SeqToTBTHDF5PluginTest {

    @Test
    public void testPerformFunction() throws Exception {
        GBSv1SimData sim = GBSv1SimData.createUnder("SeqToTBT");
        sim.buildMasterTagCounts();
        sim.buildTbt();

        TagCounts masterTags = new TagCounts(sim.masterTagCounts.toString(), FilePacking.Byte);
        int expectedTags = masterTags.getTagCount();
        assertTrue("Master tag list should be non-empty", expectedTags > 0);

        TagsByTaxaByteHDF5TaxaGroups tbt = new TagsByTaxaByteHDF5TaxaGroups(sim.tbtFile.toString());
        try {
            assertEquals("TBT taxa count should match the simulated taxa",
                    sim.sim.taxa.size(), tbt.getTaxaCount());
            assertEquals("TBT tag count should match the master tag list",
                    expectedTags, tbt.getTagCount());

            long totalReads = 0;
            for (int tag = 0; tag < tbt.getTagCount(); tag++) {
                for (int taxon = 0; taxon < tbt.getTaxaCount(); taxon++) {
                    totalReads += tbt.getReadCountForTagTaxon(tag, taxon);
                }
            }
            assertTrue("The TBT should contain matched reads", totalReads > 0);
        } finally {
            tbt.closeWriter();
        }
    }
}
