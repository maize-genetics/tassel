/*
 * TBTTestUtils
 */
package net.maizegenetics.dna.tag;

import net.maizegenetics.dna.BaseEncoder;
import static org.junit.Assert.assertEquals;

/**
 *
 * @author terry
 */
public class TBTTestUtils {

    private TBTTestUtils() {
        // utility
    }

    public static void compareTBTs(TagsByTaxa expected, TagsByTaxa actual) {

        assertEquals("Tag Count: ", expected.getTagCount(), actual.getTagCount());
        assertEquals("Tag Size in Long: ", expected.getTagSizeInLong(), actual.getTagSizeInLong());
        assertEquals("Taxa Count: ", expected.getTaxaCount(), actual.getTaxaCount());

        int numTags = expected.getTagCount();
        int numTaxa = expected.getTaxaCount();

        for (int i = 0; i < numTags; i++) {

            assertEquals("Tag Length at Index: " + i, expected.getTagLength(i), actual.getTagLength(i));

            long[] expectedTag = expected.getTag(i);
            String expectedTagStr = BaseEncoder.getSequenceFromLong(expectedTag);
            long[] actualTag = actual.getTag(i);
            String actualTagStr = BaseEncoder.getSequenceFromLong(actualTag);
            assertEquals("Tag at Index: " + i, expectedTagStr, actualTagStr);

            assertEquals("Taxa Read Bits at Index: " + i, expected.getTaxaReadBitsForTag(i), actual.getTaxaReadBitsForTag(i));

        }

        for (int i = 0; i < numTaxa; i++) {
            assertEquals("Taxa Name at Index: " + i, expected.getTaxaName(i), actual.getTaxaName(i));   
        }

        for (int r = 0; r < numTags; r++) {
            byte[] expectedReadCounts = expected.getTaxaReadCountsForTag(r);
            byte[] actualReadCounts = actual.getTaxaReadCountsForTag(r);
            for (int t = 0; t < numTaxa; t++) {
                assertEquals("Read Count at Tag/Taxon Index: " + r + "  Taxa Index: " + t, expected.getReadCountForTagTaxon(r, t), actual.getReadCountForTagTaxon(r, t));
                assertEquals("Read Count at Tag Index: " + r + "  Taxa Index: " + t, expectedReadCounts[t], actualReadCounts[t]);
            }
        }
    }
}
