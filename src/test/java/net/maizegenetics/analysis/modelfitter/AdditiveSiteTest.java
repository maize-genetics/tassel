package net.maizegenetics.analysis.modelfitter;

import static org.junit.Assert.*;

import java.util.HashMap;

import net.maizegenetics.analysis.modelfitter.AdditiveSite.CRITERION;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

import org.junit.Ignore;
import org.junit.Test;

public class AdditiveSiteTest {
    @Ignore
    @Test
    public void testGenotypeAdditiveSite() {
        byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
        byte A = NucleotideAlignmentConstants.getNucleotideAlleleByte("A");
        byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
        byte AC = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
        byte CA = NucleotideAlignmentConstants.getNucleotideDiploidByte("CA");
        byte N = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        byte[] genotype = new byte[]{AA, AA, CC, CC, AC, CA, AA, CC, AC, AA, AA, AA, AA, N, AA, CC, N};
        
        double majorFreq = 0.633333;
        double meanCov = 2 * majorFreq;
        
        GenotypeAdditiveSite gas = new GenotypeAdditiveSite(1, "C1", 1, "test", CRITERION.pval, genotype, A, majorFreq);
        double[] cov = gas.getCovariate();
        double[] expectedCov = new double[genotype.length];
        HashMap<Byte, Double> covValue = new HashMap<>();
        covValue.put(AA, 2 - meanCov);
        covValue.put(CC, 0 - meanCov);
        covValue.put(AC, 1 - meanCov);
        covValue.put(CA, 1 - meanCov);
        covValue.put(N, 0.0);
        for (int i = 0; i < genotype.length; i++) expectedCov[i] = covValue.get(genotype[i]);
        
        assertEquals(cov.length, expectedCov.length);
        assertArrayEquals(expectedCov, cov, 1e-10);
        
    }

}
