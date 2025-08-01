/*
 * UTagPairFinder
 */
package net.maizegenetics.analysis.gbs;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.dna.tag.TagCounts;

/**
 *
 * @author Fei Lu
 */
public class UTagPairFinder {

    TagCounts tc;
    int tagSizeInLong;
    long[] lookup;
    int[] indexOfTc;

    public UTagPairFinder(TagCounts tc) {
        this.tc = tc;
        tagSizeInLong = tc.getTagSizeInLong();
        initialize();
    }

    public void initialize() {
        lookup = new long[tc.getSize() * tagSizeInLong];
        indexOfTc = new int[tc.getSize() * tagSizeInLong];
        for (int i = 0; i < tc.getSize(); i++) {
            for (int j = 0; j < tagSizeInLong; j++) {
                int dex = tagSizeInLong * i + j;
                lookup[dex] = tc.getTags()[j][i];
                indexOfTc[dex] = i;
            }
        }
        GenericSorting.quickSort(0, lookup.length, comp, swapper);
        this.reduceDuplicates();
        System.out.println("Initialization is done");
    }

    public ArrayList<Integer> findOneMismatch(long[] queryLongSeq) {
        ArrayList<Integer> hitIndex = new ArrayList();
        for (int i = 0; i < queryLongSeq.length; i++) {
            int hit = Arrays.binarySearch(lookup, queryLongSeq[i]);
            if (hit < 0) {
                continue;
            }
            while ((hit > 0) && (queryLongSeq[i] == lookup[hit - 1])) {
                hit--;
            }
            while ((hit < lookup.length) && (lookup[hit] == queryLongSeq[i])) {
                int count = 0;
                for (int j = 0; j < tagSizeInLong; j++) {
                    count += getMismatchInLong(queryLongSeq[j], tc.getTags()[j][indexOfTc[hit]]);
                }
                if (count == 1) {
                    hitIndex.add(indexOfTc[hit]);
                }
                hit++;
            }
        }
        return hitIndex;
    }

    public byte getMismatchInLong(long longSeq1, long longSeq2) {
        long mask = 3;
        long diff = longSeq1 ^ longSeq2;
        byte count = 0;
        for (int i = 0; i < 32; i++) {
            if ((diff & mask) > 0) {
                count++;
            }
            diff = diff >> 2;
        }
        return count;
    }

    private void reduceDuplicates() {
        int start = 0, end = -1, duplicated = 0;
        long currHap = lookup[0];
        for (int i = 0; i < lookup.length; i++) {
            if (lookup[i] == currHap) {
                end = i;
            } else {
                if (((end - start) > 1000)) {
                    //System.out.println(BaseEncoder.getSequenceFromInt(currHap)+" "+(end-start));
                    for (int j = start; j <= end; j++) {
                        lookup[j] = Long.MAX_VALUE;
                        duplicated++;
                    }
                }
                currHap = lookup[i];
                start = end = i;
            }
        }
        GenericSorting.quickSort(0, lookup.length, comp, swapper);
        long[] newlookup = new long[lookup.length - duplicated];
        int[] newindexOfRbt = new int[lookup.length - duplicated];
        System.arraycopy(lookup, 0, newlookup, 0, lookup.length - duplicated);
        System.arraycopy(indexOfTc, 0, newindexOfRbt, 0, indexOfTc.length - duplicated);
        System.out.println("Old Lookup Size:" + lookup.length + "  new size:" + newlookup.length);
        lookup = newlookup;
        indexOfTc = newindexOfRbt;
    }
    Swapper swapper = new Swapper() {
        public void swap(int a, int b) {
            long tl;
            int ti, tb;
            tl = lookup[a];
            lookup[a] = lookup[b];
            lookup[b] = tl;
            ti = indexOfTc[a];
            indexOfTc[a] = indexOfTc[b];
            indexOfTc[b] = ti;
        }
    };

    IntComparator comp = new IntComparator() {
        public int compare(int a, int b) {
            if (lookup[a] < lookup[b]) {
                return -1;
            } else if (lookup[a] > lookup[b]) {
                return 1;
            } else {
                return 0;
            }
        }
    };
}
