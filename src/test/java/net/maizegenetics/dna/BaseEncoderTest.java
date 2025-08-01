package net.maizegenetics.dna;

import static org.junit.Assert.*;

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

}
