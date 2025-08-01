package net.maizegenetics.analysis.filter;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import org.junit.BeforeClass;
import org.junit.Test;

public class FilterSiteNamePluginTest {
	static GenotypeTable genoTable;
	static DataSet myInputData;
	static final String covariatePrefix = "Site_Covariates_";
	static final String factorPrefix = "Site_Factors_";
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		genoTable = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_CHR_9_10_FILENAME);
		myInputData = new DataSet(new Datum("name", genoTable, "comment"), null);
	}
	

	@Test
	public void testIndicesToCovariates() {
		int[] siteNumbersToTest = new int[]{10, 20, 40};
		FilterSiteNamePlugin fsnp = new FilterSiteNamePlugin(null, false);
		fsnp.setSitesToCovariatesIndex(siteNumbersToTest);
		DataSet result = fsnp.performFunction(myInputData);
		Phenotype covariatePheno = (Phenotype) result.getData(0).getData();
		
		//test covariates
		for (int i = 0; i < siteNumbersToTest.length; i++) {
			int site = siteNumbersToTest[i];
			NumericAttribute myCovariate = (NumericAttribute) covariatePheno.attribute(i + 1);
			assertArrayEquals(covariatesFromSite(site), myCovariate.floatValues(), .001f);
		}
		
	}
	
	@Test
	public void testIndicesToFactors() {
		int[] siteNumbersToTest = new int[]{10, 20, 40};
		FilterSiteNamePlugin fsnp = new FilterSiteNamePlugin(null, false);
		fsnp.setSitesToFactorsIndex(siteNumbersToTest);
		DataSet result = fsnp.performFunction(myInputData);
		Phenotype factorPheno = (Phenotype) result.getData(0).getData();

		//test factors
		int ntaxa = genoTable.numberOfTaxa();
		for (int i = 0; i < siteNumbersToTest.length; i++) {
			int site = siteNumbersToTest[i];
			CategoricalAttribute myFactor = (CategoricalAttribute) factorPheno.attribute(i + 1);
			Object[][] freq = genoTable.genosSortedByFrequency(site);
			List<String> genoList = Arrays.stream(freq[0]).map(o -> (String) o).filter(str -> !str.equals("N")).collect(Collectors.toList());
			
			for (int t = 0; t < ntaxa; t++) {
				byte geno = genoTable.genotype(t, site);
				if (geno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) assertTrue(genoList.contains(myFactor.label(t)));
				else {
					assertEquals(NucleotideAlignmentConstants.getNucleotideIUPAC(geno), myFactor.label(t));
				}
			}
		}
		
	}
	
	@Test
	public void testNamesToCovariates() {
		int[] siteNumbersToTest = new int[]{10, 20, 40};
		String[] siteNames = Arrays.stream(siteNumbersToTest).mapToObj(s -> genoTable.siteName(s)).toArray(String[]::new);
		
		FilterSiteNamePlugin fsnp = new FilterSiteNamePlugin(null, false);
		fsnp.setSiteNamesToCovariates(siteNames);
		DataSet result = fsnp.performFunction(myInputData);
		Phenotype covariatePheno = (Phenotype) result.getData(0).getData();
		
		//test covariates
		for (int i = 0; i < siteNumbersToTest.length; i++) {
			int site = siteNumbersToTest[i];
			NumericAttribute myCovariate = (NumericAttribute) covariatePheno.attribute(i + 1);
			assertArrayEquals(covariatesFromSite(site), myCovariate.floatValues(), .001f);
		}
		
	}

	@Test
	public void testNamesToFactors() {
		int[] siteNumbersToTest = new int[]{10, 20, 40};
		String[] siteNames = Arrays.stream(siteNumbersToTest).mapToObj(s -> genoTable.siteName(s)).toArray(String[]::new);
		FilterSiteNamePlugin fsnp = new FilterSiteNamePlugin(null, false);
		fsnp.setSiteNamesToFactors(siteNames);
		DataSet result = fsnp.performFunction(myInputData);
		Phenotype factorPheno = (Phenotype) result.getData(0).getData();

		//test factors
		int ntaxa = genoTable.numberOfTaxa();
		for (int i = 0; i < siteNumbersToTest.length; i++) {
			int site = siteNumbersToTest[i];
			CategoricalAttribute myFactor = (CategoricalAttribute) factorPheno.attribute(i + 1);
			Object[][] freq = genoTable.genosSortedByFrequency(site);
			List<String> genoList = Arrays.stream(freq[0]).map(o -> (String) o).filter(str -> !str.equals("N")).collect(Collectors.toList());
			
			for (int t = 0; t < ntaxa; t++) {
				byte geno = genoTable.genotype(t, site);
				if (geno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) assertTrue(genoList.contains(myFactor.label(t)));
				else {
					assertEquals(NucleotideAlignmentConstants.getNucleotideIUPAC(geno), myFactor.label(t));
				}
			}
		}
		
	}
	private float[] covariatesFromSite(int site) {
		final float sitemean = (float) (2 * genoTable.majorAlleleFrequency(site));
		final int nsites = genoTable.numberOfSites();
		final int ntaxa = genoTable.numberOfTaxa();
		final byte major = genoTable.majorAllele(site);
		
		float[] result = new float[ntaxa];
		for (int t = 0; t < ntaxa; t++) {
			byte geno = genoTable.genotype(t, site);
			if (geno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) result[t] = sitemean;
			else {
				byte[] genoArray = GenotypeTableUtils.getDiploidValues(geno);
				float score = 0;
				if (genoArray[0] == major) score++;
				if (genoArray[1] == major) score++;
				result[t] = score;
			}
		}
		return result;
	}


}
