/*
 * TagCountsTestUtils
 */
package net.maizegenetics.dna.tag;

import java.util.ArrayList;
import java.util.List;
import net.maizegenetics.dna.BaseEncoder;
import static org.junit.Assert.assertEquals;

/**
 *
 * @author terry
 */
public class TagCountsTestUtils {

    private TagCountsTestUtils() {
        // utility
    }

    public static void sanityCheck(TagCounts tagCounts, int initNumRawReads) {

        int numTags = tagCounts.getTagCount();
        int maxTagLength = tagCounts.getTagSizeInLong() * 32;

        List<String> tagList = new ArrayList<String>();
        int runningNumTags = 0;
        for (int i = 0; i < numTags; i++) {

            long[] currentTag = tagCounts.getTag(i);
            String currentTagStr = BaseEncoder.getSequenceFromLong(currentTag);

            assertEquals("Duplicate Tag Found at Index: " + i, false, tagList.contains(currentTagStr));
            tagList.add(currentTagStr);

            for (int j = tagCounts.getTagLength(i); j < maxTagLength; j++) {
                assertEquals("Tag: " + currentTagStr + " Length: " + tagCounts.getTagLength(i) + " Doesn't have correct number of padded As", 'A', currentTagStr.charAt(j));
            }

            int currentCount = tagCounts.getReadCount(i);
            assertEquals("Tag at Index: " + i + " has Count less than 1: " + currentCount, false, (currentCount <= 0));
            runningNumTags += currentCount;
        }

        System.out.println("Initial Num Raw Read: " + initNumRawReads + "  Number Tags: " + numTags + " Accumulated Count: " + runningNumTags);
        if (initNumRawReads != -1) {
            assertEquals("Initial Num Raw Read: " + initNumRawReads + " is Less Than Accumulated Count: " + runningNumTags, true, (runningNumTags <= initNumRawReads));
        }

    }
}
