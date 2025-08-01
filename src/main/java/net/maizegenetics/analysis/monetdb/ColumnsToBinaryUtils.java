/**
 * 
 */
package net.maizegenetics.analysis.monetdb;

import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

/**
 * This class contains utility functions to aid in creating monetdb binaries
 * 
 * @author lcj34
 *
 */
public class ColumnsToBinaryUtils {
    //NOTE: assumption is that a 1-based position was passed in
    // THis is true when we pass position from existing db file vs
    // positions from refGenomeFile
    public static String refAlleleFromChromPos(byte[] refChromBytes, int allelePos) {
        // The reference genome is stored in packed bytes. WHen chromosomeSequence()
        // is called, it unpacks the bytes and returns requested byte sequence
        // Test if we can get just 1
        
        // We  get the chromosome into bytes, then send in just this
        // chromosome and get a byte at a particular position.  Is is 0 or 1 based?
        // it looks like it fixes it, so can send in 1 based - see the GEnomeSequenceBUidler.test file
        // Acutally, chromosomeSequence() fixeds it.  But our bytes returned are
        // in a 0 based array.  IE, we translated to 0 before grabbing from the reference,
        // then stored the result in a 0 based array.
        byte[] oneAllele = new byte[1];
        // I am leaving the assumption that we passed in a 1-based position
        oneAllele[0] = refChromBytes[allelePos-1];
        return NucleotideAlignmentConstants.nucleotideBytetoString(oneAllele);        
    }
    
    // A copy of this lives in net.maizegenetics.dna.BaseEncoder
    // That copy does not contain the +/- (in/del) values
    public static int getIntFromSeq(String seq) {
        int chunkSize = 16;
        int seqLength = seq.length();
        int v = 0;
        for (int i = 0; i < seqLength; i++) {
            switch (seq.charAt(i)) {
                case 'A':
                case 'a':
                    v = v << 2;
                    break;
                case 'C':
                case 'c':
                    v = (v << 2) + (byte) 1;
                    break;
                case 'G':
                case 'g':
                    v = (v << 2) + (byte) 2;
                    break;
                case 'T':
                case 't':
                    v = (v << 2) + (byte) 3;
                    break;
                case '+':
                    v = (v << 2) + (byte) 4;
                    break;
                case '-':
                    v = (v << 2) + (byte) 5;
                    break;
                default:
                    return -1;
            }
        }
        if (seqLength == chunkSize) {
            return v;
        }
        if (seqLength > chunkSize) {
            return -1;
        }
        // LCJ - try for 0's (A's) in front
        //v = (v << (2 * (chunkSize - seqLength))); //if shorter fill with AAAA (which is 0000)
        return v;
    }
}
