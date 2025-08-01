package net.maizegenetics.analysis.imputation;

import static org.junit.Assert.*;

import java.util.Random;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;

import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

public class RandomGenotypeImputationTest {
	private static GenotypeTable myUnimputedGenotype;
	private static int ntaxa = 1000;
	private static Random ran = new Random();
	
	@BeforeClass
	public static void setUpClass() {
		TaxaListBuilder taxaBuilder = new TaxaListBuilder();
		for (int t = 0; t < ntaxa; t++) {
			taxaBuilder.add(new Taxon("t" + t));
		}
		
		GenotypeTableBuilder myBuilder = GenotypeTableBuilder.getSiteIncremental(taxaBuilder.build());
		
		//simulate some sites with missing data
		Chromosome mychr = new Chromosome("1");
		myBuilder.addSite(new GeneralPosition.Builder(mychr, 1).build(), simulateGenotype(0.2,0.2,0.4));
		myBuilder.addSite(new GeneralPosition.Builder(mychr, 2).build(), simulateGenotype(0.4,0.2,0.1));
		myBuilder.addSite(new GeneralPosition.Builder(mychr, 3).build(), simulateGenotype(0.1,0.1,0.2));
		myBuilder.addSite(new GeneralPosition.Builder(mychr, 4).build(), simulateGenotype(0.4,0.2,0.2));
		myBuilder.addSite(new GeneralPosition.Builder(mychr, 5).build(), simulateGenotype(0.4,0.3,0.1));
		
		myUnimputedGenotype = myBuilder.build();
	}
	
	public static byte[] simulateGenotype(double phom1, double phom2, double phet) {
		byte hom1 = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
		byte hom2 = NucleotideAlignmentConstants.getNucleotideDiploidByte("C");
		byte het = NucleotideAlignmentConstants.getNucleotideDiploidByte("M");
		byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
		
		double phom12 = phom1 + phom2;
		double pnotmissing = phom12 + phet;
		
		byte[] geno = new byte[ntaxa];
		for (int t = 0; t < ntaxa; t++) {
			double ranval = ran.nextDouble();
			if (ranval < phom1) geno[t] = hom1;
			else if (ranval < phom12) geno[t] = hom2;
			else if (ranval < pnotmissing) geno[t] = het;
			else geno[t] = NN;
		}
		return geno;
	}
	
	@Test
	public void test() {
		RandomGenotypeImputationPlugin rgip = new RandomGenotypeImputationPlugin(null, false);
		DataSet input = new DataSet(new Datum("name", myUnimputedGenotype,"nothing"), null);
		DataSet output = rgip.processData(input);
		GenotypeTable myImputedGenotype = (GenotypeTable) output.getData(0).getData();
		
		for (int i = 0; i < 5; i++) {
			int imputedNotMissing = myImputedGenotype.totalNonMissingForSite(i);
			assertEquals("Failed none missing test.", ntaxa, imputedNotMissing);
			Object[][] genoFreq1 = myUnimputedGenotype.genosSortedByFrequency(i);
			Object[][] genoFreq2 = myImputedGenotype.genosSortedByFrequency(i);
			
			int[] count1 = new int[3];
			int sum1 = 0;
			String[] genostr1 = new String[3];
			int ndx = 0;
			for (int j = 0; j < 4; j++) {
				if (!genoFreq1[0][j].equals("N")) {
					count1[ndx] = ((Integer) genoFreq1[1][j]).intValue();
					sum1 += count1[ndx];
					genostr1[ndx] = (String) genoFreq1[0][j];
					ndx++;
				}
			}
			
			int[] count2 = new int[3];
			int sum2 = 0;
			String[] genostr2 = new String[3];
			for (int j = 0; j < 3; j++) {
				count2[j] = ((Integer) genoFreq2[1][j]).intValue();
				sum2 += count2[j];
				genostr2[j] = (String) genoFreq2[0][j];
			}
			
			for (int g1 = 0; g1 < 3; g1++) {
				int g2 = 0;
				while (!genostr1[g1].equals(genostr2[g2])) g2++;
				double prg1 = ((double) count1[g1]) / (double) sum1;
				double prg2 = ((double) count2[g2]) / (double) sum2;
				assertEquals(String.format("Failed proportions equal for snp %d, genotype = %s", i, genostr1[g1]), prg1, prg2, 0.05);
//				System.out.printf("%s: unimputed p = %1.4f, imputed p = %1.4f\n", genostr1[g1], prg1, prg2);
			}
		}
	}

}
