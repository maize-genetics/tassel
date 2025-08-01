package net.maizegenetics.analysis.distance;

import net.maizegenetics.dna.snp.GenotypeTable;
import org.junit.Test;

import static org.junit.Assert.*;

public class IBSDistanceMatrixTest {

    @Test
    public void testComputeHetDistances() throws Exception {
        System.out.println("IBSDistanceMatrixTest.testComputeHetDistances");
        byte[] seq1 = {0x00, 0x11, 0x22};
        byte[] seq2 = {0x00, 0x11, 0x22};
        assertEquals(0.0, IBSDistanceMatrix.computeHetDistances(seq1, seq2, 1)[0], 0.001);  //identical

        seq2 = new byte[]{0x00, 0x22, 0x22};
        assertEquals(0.33333, IBSDistanceMatrix.computeHetDistances(seq1, seq2, 1)[0], 0.001);  //1 diff

        seq2 = new byte[]{0x00, 0x13, 0x22};
        assertEquals(0.16666, IBSDistanceMatrix.computeHetDistances(seq1, seq2, 1)[0], 0.001);  //1 het diff

        seq2 = new byte[]{0x01, 0x13, 0x12};
        assertEquals(0.5, IBSDistanceMatrix.computeHetDistances(seq1, seq2, 1)[0], 0.001);  //3 het diff

        seq1 = new byte[]{0x10, 0x31, 0x12};
        seq2 = new byte[]{0x01, 0x13, 0x12};
        assertEquals(0.5, IBSDistanceMatrix.computeHetDistances(seq1, seq2, 1)[0], 0.001);  //all het to het
        assertEquals(3, IBSDistanceMatrix.computeHetDistances(seq1, seq2, 1)[1], 0.001);  //test missing site not counted

        seq1 = new byte[]{GenotypeTable.UNKNOWN_DIPLOID_ALLELE, 0x11, 0x22};
        seq2 = new byte[]{0x01, 0x11, 0x11};
        assertEquals(0.5, IBSDistanceMatrix.computeHetDistances(seq1, seq2, 1)[0], 0.001);  //test missing site not counted
        assertEquals(2, IBSDistanceMatrix.computeHetDistances(seq1, seq2, 1)[1], 0.001);  //test missing site not counted

        seq1 = new byte[]{GenotypeTable.UNKNOWN_DIPLOID_ALLELE, 0x11, 0x22};
        seq2 = new byte[]{0x01, 0x11, 0x11};
        assertEquals(Double.NaN, IBSDistanceMatrix.computeHetDistances(seq1, seq2, 3)[0], 0.001);  //test NaN thrown if too few sites
        assertEquals(2, IBSDistanceMatrix.computeHetDistances(seq1, seq2, 3)[1], 0.001);  //test missing site not counted

    }
}
