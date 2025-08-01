/*
 * UNetworkFilter
 */
package net.maizegenetics.analysis.gbs;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.BaseEncoder;

/**
 *
 * @author Fei Lu
 */
public class UNetworkFilter {

    TagAlignment[] tps;

    public UNetworkFilter(String mergedTagCountFileS, double etr, String tagPairFileS) {
        TagCounts tc = new TagCounts(mergedTagCountFileS, FilePacking.Byte);
        this.getTagPair(tc, etr);
        this.reciprocate();
        this.writeTagPair(tc, tagPairFileS);
    }

    public final void getTagPair(TagCounts tc, double etr) {
        UTagPairFinder tf = new UTagPairFinder(tc);
        System.out.println("Sequencing error tolerance rate is " + String.valueOf(etr));
        ArrayList<TagAlignment> tpList = new ArrayList();
        int cnt1 = 0, cnt2 = 0;
        for (int i = 0; i < tc.getSize(); i++) {
            long[] queryLongSeq = new long[2];
            queryLongSeq[0] = tc.getTags()[0][i];
            queryLongSeq[1] = tc.getTags()[1][i];
            ArrayList<Integer> hitIndex = tf.findOneMismatch(queryLongSeq);
            if (hitIndex.isEmpty()) {
                continue;
            }
            Integer[] hitIndexArray = hitIndex.toArray(new Integer[hitIndex.size()]);
            if (hitIndexArray.length == 1) {
                TagAlignment cl = new TagAlignment(i, hitIndexArray[0]);
                tpList.add(cl);
                cnt1++;
            } else {
                int maxCountHitIndex = this.errorFilter(tc, etr, i, hitIndexArray);
                if (maxCountHitIndex == -1) {
                    continue;
                }
                TagAlignment cl = new TagAlignment(i, maxCountHitIndex);
                tpList.add(cl);
                cnt2++;
            }
        }
        tps = tpList.toArray(new TagAlignment[tpList.size()]);
        System.out.println(tps.length + " TagAlignments found. " + cnt1 + " are size = 1, " + cnt2 + " are error curated from size > 1");
    }

    public int errorFilter(TagCounts rc, double errorRate, int queryIndex, Integer[] hitIndexArray) {// each other count < erroRate, errorrate is 0.01 or 0.02
        int[] counts = new int[hitIndexArray.length + 1];
        int queryCount = rc.getReadCount(queryIndex);
        counts[0] = queryCount;
        for (int i = 0; i < hitIndexArray.length; i++) {
            counts[i + 1] = rc.getReadCount(hitIndexArray[i]);
        }
        Arrays.sort(counts);
        int maxIndex = -1;
        int hit = Arrays.binarySearch(counts, queryCount);
        if (hit == (counts.length - 1) || hit == (counts.length - 2)) {
            if (queryCount == counts[counts.length - 3]) {
                return -1;
            }
            for (int i = 0; i < counts.length - 2; i++) {
                if ((double) counts[i] / (double) (counts[i] + queryCount) > errorRate) {
                    return -1;
                }
            }
            int max = 0;
            for (int i = 0; i < hitIndexArray.length; i++) {
                int cnt = rc.getReadCount(hitIndexArray[i]);
                if (cnt > max) {
                    max = cnt;
                    maxIndex = hitIndexArray[i];
                }
            }

        }
        return maxIndex;
    }

    public void reciprocate() {
        for (int i = 0; i < tps.length; i++) {
            tps[i].swapQueryHit();
        }
        Arrays.sort(tps);
        ArrayList<TagAlignment> tpList = new ArrayList();
        TagAlignment currentTagPair = tps[0];
        int count = 1;
        for (int i = 1; i < tps.length; i++) {
            if (tps[i].equals(currentTagPair)) {
                count++;
            } else {
                if (count == 2) {
                    tpList.add(currentTagPair);
                }
                currentTagPair = tps[i];
                count = 1;
            }
        }
        if (count == 2) {
            tpList.add(currentTagPair);
        }
        TagAlignment[] tpArray = tpList.toArray(new TagAlignment[tpList.size()]);
        tps = tpArray;
        System.out.println(tps.length + " reciprocal TagPairs found");
    }

    public void writeTagPair(TagCounts tc, String outfileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
            dos.writeInt(tps.length * 2);
            dos.writeInt(tc.getTagSizeInLong());
            for (int i = 0; i < tps.length; i++) {
                tps[i].writeTagPair(tc, dos, 2 * i);
            }
            dos.flush();
            dos.close();
            System.out.println("TagPair file is written to " + new File(outfileS).getAbsolutePath());
        } catch (Exception e) {
            System.out.println("Error occurred while writing TagPair file " + outfileS + " " + e.toString());
        }
    }

    public void screenPrintTps() {
        int total = 1000;
        for (int i = 0; i < total; i++) {
            System.out.println(tps[i].queryIndex + "\t" + tps[i].hitIndex);
        }
    }

    public void screenPrint(TagCounts rc) {
        int total = 10;
        if (tps.length < total) {
            total = tps.length;
        }
        for (int i = 0; i < total; i++) {
            long[] lseq = new long[2];
            lseq[0] = rc.getTags()[0][tps[i].queryIndex];
            lseq[1] = rc.getTags()[1][tps[i].queryIndex];
            System.out.println(BaseEncoder.getSequenceFromLong(lseq));
            lseq[0] = rc.getTags()[0][tps[i].hitIndex];
            lseq[1] = rc.getTags()[1][tps[i].hitIndex];
            System.out.println(BaseEncoder.getSequenceFromLong(lseq));
            System.out.println();
        }
    }

    class TagAlignment implements Comparable<TagAlignment> {

        int queryIndex;
        int hitIndex;

        TagAlignment(int queryIndex, int hitIndex) {
            this.queryIndex = queryIndex;
            this.hitIndex = hitIndex;
        }

        void swapQueryHit() {
            if (queryIndex > hitIndex) {
                int mid = queryIndex;
                queryIndex = hitIndex;
                hitIndex = mid;
            }
        }

        void writeTagPair(TagCounts tc, DataOutputStream dos, int order) {
            try {
                dos.writeLong(tc.getTags()[0][queryIndex]);
                dos.writeLong(tc.getTags()[1][queryIndex]);
                dos.writeByte(tc.getTagLength(queryIndex));
                dos.writeInt(order);
                dos.writeLong(tc.getTags()[0][hitIndex]);
                dos.writeLong(tc.getTags()[1][hitIndex]);
                dos.writeByte(tc.getTagLength(hitIndex));
                dos.writeInt(order + 1);
            } catch (Exception e) {
                System.out.println(e.toString());
            }
        }

        public boolean equals(TagAlignment o) {
            if (queryIndex == o.queryIndex && hitIndex == o.hitIndex) {
                return true;
            }
            return false;
        }

        public int compareTo(TagAlignment o) {
            if (queryIndex < o.queryIndex) {
                return -1;
            } else if (queryIndex > o.queryIndex) {
                return 1;
            } else {
                if (hitIndex < o.hitIndex) {
                    return -1;
                } else if (hitIndex > o.hitIndex) {
                    return 1;
                } else {
                    return 0;
                }
            }
        }

    }
}
