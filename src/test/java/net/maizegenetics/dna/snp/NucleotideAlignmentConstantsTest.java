package net.maizegenetics.dna.snp;

import junit.framework.Assert;
import org.junit.Test;

import java.util.Map;

/**
 * Defines basic tests that there is integrity in round tripping of IUPAC codes
 *
 * @author Ed Buckler
 */
public class NucleotideAlignmentConstantsTest {

    @Test
    public void testRoundTripOfIUPACCodes() throws Exception {
        for (Map.Entry<Byte,String> byteStringEntry : NucleotideAlignmentConstants.NUCLEOTIDE_IUPAC_HASH.entrySet()) {
            byte d=NucleotideAlignmentConstants.getNucleotideDiploidByte(byteStringEntry.getValue().charAt(0));
            String s=NucleotideAlignmentConstants.getNucleotideIUPAC(d);
            Assert.assertEquals(byteStringEntry.getValue(),s);
            char c=NucleotideAlignmentConstants.getNucleotideIUPACChar(d);
            Assert.assertEquals(byteStringEntry.getValue().charAt(0),c);
        }

        for (int i=Byte.MIN_VALUE; i<=Byte.MAX_VALUE; i++) {
            byte bi=(byte)i;
            System.out.println(bi);
            Assert.assertEquals(NucleotideAlignmentConstants.getNucleotideIUPAC(bi).charAt(0),
                    NucleotideAlignmentConstants.getNucleotideIUPACChar(bi));
        }
    }

    @Test
    public void testComplementIUPACCodes() throws Exception {
        for (String iupac : NucleotideAlignmentConstants.NUCLEOTIDE_IUPAC_HASH.values()) {
            byte b=NucleotideAlignmentConstants.getNucleotideDiploidByte(iupac.charAt(0));
            System.out.println(b);
            System.out.println(NucleotideAlignmentConstants.getNucleotideDiploidComplement(b));
            char c1=iupac.charAt(0);
            char compC1=NucleotideAlignmentConstants.getNucleotideDiploidIUPACComplement(iupac.charAt(0));
            Assert.assertEquals(compC1,complement(c1));
            System.out.printf("%c %c %c %n",c1,compC1,complement(c1));
        }
    }

    public static char complement(char geno) {
        char comp = 'X';
        switch (geno) {
            case 'A':  comp = 'T';  break;
            case 'C':  comp = 'G';  break;
            case 'G':  comp = 'C';  break;
            case 'T':  comp = 'A';  break;
            case 'K':  comp = 'M';  break;
            case 'M':  comp = 'K';  break;
            case 'R':  comp = 'Y';  break;
            case 'S':  comp = 'S';  break;
            case 'W':  comp = 'W';  break;
            case 'Y':  comp = 'R';  break;
            case '-':  comp = '-';  break;  // both strands have the deletion
            case '+':  comp = '+';  break;  // both strands have the insertion
            case '0':  comp = '0';  break;
            case 'N':  comp = 'N';  break;
            case 'Z':  comp = 'Z';  break;
            default:   comp = 'N';  break;
        }
        return comp;
    }
}
