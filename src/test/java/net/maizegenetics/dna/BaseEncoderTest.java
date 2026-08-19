package net.maizegenetics.dna;

import static org.junit.Assert.*;

import java.util.Random;

import org.junit.Test;

import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

public class BaseEncoderTest {

    @Test
    public void testGetLongFromSeq() {
        String[] testSequence = new String[] {"ACGTTGGCAATCATGGTAACCTGCTACTAGCT",
                "GTACGTTAGGGTACACCGTACACCCATGCAGT",
                "GTACG-TAGGGTACACCGTACACCCATGC-GT",
                "GTACGT+AGGGTACACCGTACACCCATGCAGT",
                "GTACGTTASGGTACACCGTACACCCATGCAGT"};
        
        int n = testSequence.length;
        for (int i = 0; i < n; i++) {
            assertEquals(testSequence[i].length(), 32);
            //convert to bytes
            byte[] byteSeq = NucleotideAlignmentConstants.convertHaplotypeStringToAlleleByteArray(testSequence[i]);
            assertEquals(byteSeq.length, 32);
            if (i < 2) { //only ACGT
                long val = BaseEncoder.getLongFromSeq(testSequence[i]);
                assertFalse(val == -1L);
                String stringVal = BaseEncoder.getSequenceFromLong(val);
                assertEquals(testSequence[i], stringVal);
                val = BaseEncoder.getLongSeqFromByteArray(byteSeq);
                assertFalse(val == -1L);
                stringVal = BaseEncoder.getSequenceFromLong(val);
                assertEquals(testSequence[i], stringVal);
            } else { // contains other characters
                long val = BaseEncoder.getLongFromSeq(testSequence[i]);
                assertEquals(-1L, val);
                val = BaseEncoder.getLongSeqFromByteArray(byteSeq);
                assertEquals(-1L, val);
            }
        }

    }

    // ---- seqDifferences: reference copies of the ORIGINAL bit-loop implementations, kept here to
    // prove the new bit-parallel (collapse-then-popcount) versions in BaseEncoder are result-identical.

    private static byte refSeqDifferences(long seq1, long seq2) {
        long mask = 3;
        byte cnt = 0;
        long diff = seq1 ^ seq2;
        for (int x = 0; x < 32; x++) {
            if ((diff & mask) > 0) cnt++;
            diff = diff >> 2;
        }
        return cnt;
    }

    private static byte refSeqDifferences(long seq1, long seq2, int maxDivergence) {
        long mask = 3;
        byte cnt = 0;
        long diff = seq1 ^ seq2;
        for (int x = 0; x < 32 && cnt <= maxDivergence; x++) {
            if ((diff & mask) > 0) cnt++;
            diff = diff >> 2;
        }
        if (cnt > maxDivergence) cnt = (byte) 32;
        return cnt;
    }

    private static byte refSeqDifferencesForSubset(long seq1, long seq2, int lengthOfComp, int maxDivergence) {
        long mask = 3;
        byte cnt = 0;
        long diff = seq1 ^ seq2;
        diff = diff >> (2 * (32 - lengthOfComp));
        for (int x = 0; x < lengthOfComp && cnt < maxDivergence; x++) {
            if ((diff & mask) > 0) cnt++;
            diff = diff >> 2;
        }
        return cnt;
    }

    @Test
    public void testSeqDifferencesEquivalence() {
        Random r = new Random(42);
        for (int i = 0; i < 50_000; i++) {
            long a = r.nextLong();   // full 64-bit range: exercises the sign bit (base 31, bits 62-63)
            long b = r.nextLong();
            assertEquals("2-arg", refSeqDifferences(a, b), BaseEncoder.seqDifferences(a, b));
            for (int maxDiv = 0; maxDiv <= 33; maxDiv += 3) {
                assertEquals("3-arg maxDiv=" + maxDiv,
                        refSeqDifferences(a, b, maxDiv), BaseEncoder.seqDifferences(a, b, maxDiv));
            }
            for (int len = 0; len <= 32; len += 4) {
                for (int maxDiv = 0; maxDiv <= 33; maxDiv += 8) {
                    assertEquals("subset len=" + len + " maxDiv=" + maxDiv,
                            refSeqDifferencesForSubset(a, b, len, maxDiv),
                            BaseEncoder.seqDifferencesForSubset(a, b, len, maxDiv));
                }
            }
        }
    }

    @Test
    public void testSeqDifferencesEdgeCases() {
        long mixed = BaseEncoder.getLongFromSeq("ACGTACGTACGTACGTACGTACGTACGTACGT");
        long allA  = BaseEncoder.getLongFromSeq("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        long allT  = BaseEncoder.getLongFromSeq("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        long base0T = BaseEncoder.getLongFromSeq("TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); // only base 0 differs from allA

        assertEquals(0, BaseEncoder.seqDifferences(mixed, mixed));         // identical
        assertEquals(32, BaseEncoder.seqDifferences(allA, allT));          // every base differs
        assertEquals(1, BaseEncoder.seqDifferences(allA, base0T));         // A<->T is 1, NOT 2 (anti-bitCount)
        assertTrue("base0=T sets the sign bit", base0T < 0);              // difference is in the sign-bit base
        assertEquals(32, BaseEncoder.seqDifferences(allA, allT, 5));       // over threshold -> chunkSize
        assertEquals(0, BaseEncoder.seqDifferences(mixed, mixed, 5));
        assertEquals(1, BaseEncoder.seqDifferencesForSubset(allA, base0T, 1, 5));  // first 1 base differs
        assertEquals(1, BaseEncoder.seqDifferencesForSubset(allA, base0T, 32, 5)); // all 32, one differs
        assertEquals(0, BaseEncoder.seqDifferencesForSubset(allA, allT, 0, 5));    // lengthOfComp 0 -> 0
    }

}
