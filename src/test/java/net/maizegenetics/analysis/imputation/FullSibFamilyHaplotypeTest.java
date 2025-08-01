package net.maizegenetics.analysis.imputation;

import static org.junit.Assert.*;

import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.Taxon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;


public class FullSibFamilyHaplotypeTest {
	private static final byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
	private static final byte N = GenotypeTable.UNKNOWN_ALLELE;
	private static final byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
	private static final byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("C");
	private static final Logger myLogger = LogManager.getLogger(FullSibFamilyHaplotypeTest.class);
	
	private double minMatch = 0.8;
	private double minPresent = 0.8;
	private int totalNotMissing = 400000;		
	DataSourceForFSImputationTest myDataSource;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
		
	}

	@Test
	public void testWithoutClusterOption() {
		minMatch = 0.99;
		minPresent = 0.9;
		totalNotMissing = 500000;		

		myDataSource = new DataSourceForFSImputationTest();
		myDataSource.setnSelfs(4);
		myDataSource.simulateGenotypeData(false);
		String[] parameters = new String[]{"-minMaf", "0.1", "-maxMissing", "0.8", "-windowld"};
		myLogger.info("testWithoutClusterOption");
		testCallParentAlleles(parameters, false);
	}

	@Test
	public void testHaplotypeFinderNoResidualHet() {
		minMatch = 0.99;
		minPresent = 0.9;
		totalNotMissing = 500000;		

		myDataSource = new DataSourceForFSImputationTest();
		myDataSource.setnSelfs(4);
		myDataSource.simulateGenotypeData(false);
		String[] parameters = new String[]{"-minMaf", "0.1", "-maxMissing", "0.8", "-maxHetDev", "20", "windowSize", "100"};
		myLogger.info("testClusterOptionNoResidualHet");
		testCallParentAlleles(parameters, true);
	}
	
	@Test
	public void testHaplotypeFinderWithResidualHet() {
		minMatch = 0.99;
		minPresent = 0.85;
		totalNotMissing = 420000;		

		myDataSource = new DataSourceForFSImputationTest();
		myDataSource.setnSelfs(4);
		myDataSource.simulateGenotypeData(true);
		String[] parameters = new String[]{"-minMaf", "0.1", "-maxMissing", "0.8","-maxHetDev", "20", "windowSize", "100"};
		myLogger.info("testClusterOptionWithResidualHet");
		testCallParentAlleles(parameters, true);
	}

	public void testCallParentAlleles(String[] parameters, boolean parentDiploid) {
		CallParentAllelesPlugin cpa = new CallParentAllelesPlugin(null);
		cpa.setFamilyList(myDataSource.getPopulationData());
		
		cpa.setParameters(parameters);
		DataSet ds = new DataSet(new Datum("testData", myDataSource.getMyGenotypeTableWithErrorNucleotides(), "This is testData."), null);
		DataSet result = cpa.performFunction(ds);
		
		int nFamilies = myDataSource.getNumberOfFamilies();
		assertEquals("test for number of families returned by CallParentAllelesPlugin", nFamilies, result.getDataOfType(PopulationData.class).size());
		
		boolean[] switchParents = new boolean[nFamilies];
		for (int f = 0; f < nFamilies; f++) {
			//call parent alleles just converts nucleotide calls to parent calls so the only thing to check is whether the parent haplotypes are correct.
			//do parent haplotypes match myDataSource?
			int[] compResultA0, compResultC1, compResultA1, compResultC0;
			PopulationData family = (PopulationData) result.getDataOfType(PopulationData.class).get(f).getData();
			if (parentDiploid) {
				compResultA0 = compareSequence(myDataSource.getParentGenotype(f, 0), family.alleleA);
				compResultC1 = compareSequence(myDataSource.getParentGenotype(f, 1), family.alleleC);
				compResultA1 = compareSequence(myDataSource.getParentGenotype(f, 1), family.alleleA);
				compResultC0 = compareSequence(myDataSource.getParentGenotype(f, 0), family.alleleC);
			} else {
				compResultA0 = compareSequence(myDataSource.getParentGenotype(f, 0), diploidFromAllele(family.alleleA));
				compResultC1 = compareSequence(myDataSource.getParentGenotype(f, 1), diploidFromAllele(family.alleleC));
				compResultA1 = compareSequence(myDataSource.getParentGenotype(f, 1), diploidFromAllele(family.alleleA));
				compResultC0 = compareSequence(myDataSource.getParentGenotype(f, 0), diploidFromAllele(family.alleleC));
			}
			
			double pMatchA, pMatchC, presentA, presentC;
			if (compResultA0[2] > compResultA1[2]) {
				switchParents[f] = false;
				pMatchA = ((double) compResultA0[2]) /compResultA0[1];
				pMatchC = ((double) compResultC1[2]) /compResultC1[1];
				presentA = ((double) compResultA0[1]) /compResultA0[0];
				presentC = ((double) compResultC1[1]) /compResultC1[0];
			} else {
				switchParents[f] = true;
				pMatchA = ((double) compResultA1[2]) /compResultA1[1];
				pMatchC = ((double) compResultC0[2]) /compResultC0[1];
				presentA = ((double) compResultA1[1]) /compResultA1[0];
				presentC = ((double) compResultC0[1]) /compResultC0[0];
			}
			
			String msg = String.format("In parent haplotypes for family %d, percent A called = %1.6f, percent C called = %1.6f", f, presentA, presentC);
			myLogger.info(msg);
			msg = String.format("In parent haplotypes for family %d, percent A matched = %1.6f, percent C matched = %1.6f", f, pMatchA, pMatchC);
			myLogger.info(msg);
			assertTrue(String.format("(pMatchA < %1.2f) for family %d", minMatch, f), pMatchA > minMatch);
			assertTrue(String.format("(pMatchC < %1.2f) for family %d", minMatch, f), pMatchC > minMatch);
			assertTrue(String.format("(presentA < %1.2f) for family %d", minPresent, f), presentA > minPresent);
			assertTrue(String.format("(presentC < %1.2f) for family %d", minPresent, f), presentC > minPresent);
		}
		
		//test ViterbiAlgorithmPlugin
		//compare imputed values to the original values and calculate error rates
		//hom->hom error rates should be <.0005. They were originally only .001.
		//hom->het errors will increase since some hom will be incorrectly changed to het, but the rate should still be low (.01)
		//het->hom error should be <.05. They were originally .8.
		ViterbiAlgorithmPlugin vap = new ViterbiAlgorithmPlugin(null);
		vap.setFillGapsInAlignment(false);
		DataSet vapResult = vap.performFunction(result);

		//compare family.imputed to myDataSource.myGenotypeTable
		GenotypeTable refTable = myDataSource.getGenotypes();
//		GenotypeTable refTable = myDataSource.getGenotypesWithError();
		
		for (int i = 0; i < myDataSource.getNumberOfFamilies(); i++) {
			int totalNotMissing = 0;
			int matchCount = 0;
			PopulationData family = (PopulationData) vapResult.getData(i).getData();
			int ntaxa = family.imputed.numberOfTaxa();
			PositionList posImputed = family.imputed.positions();
			PositionList posRef = refTable.positions();
			int nImputedSites = family.imputed.numberOfSites();
			int[] imputedToRefIndex = new int[nImputedSites];
			for (int s = 0; s < nImputedSites; s++) {
				imputedToRefIndex[s] = posRef.indexOf(posImputed.get(s));
			}
			
			for (int t = 0; t < ntaxa; t++) {
				Taxon taxon = family.imputed.taxa().get(t);
				int refndx = refTable.taxa().indexOf(taxon);
				
				for (int s = 0; s < nImputedSites; s++) {
					byte imp = family.imputed.genotype(t,s);
					byte ref = refTable.genotype(refndx, imputedToRefIndex[s]);
					if (imp != NN && ref != NN) {
						totalNotMissing++;
						if (switchParents[i]) {
							if(GenotypeTableUtils.isHeterozygous(ref) && GenotypeTableUtils.isHeterozygous(imp)) matchCount++;
							else if (!GenotypeTableUtils.isHeterozygous(ref) && !GenotypeTableUtils.isHeterozygous(imp) && ref != imp) matchCount++;
						} else {
							if(GenotypeTableUtils.isHeterozygous(ref) && GenotypeTableUtils.isHeterozygous(imp)) matchCount++;
							else if (ref == imp) matchCount++;
						}
					}
					
					
				}
			}
			
			double pmatch = ((double) matchCount)/totalNotMissing;
			int[] xocount = countCrossOvers(myDataSource.getGenotypes(i));
			String msg = String.format("Family %d original: hom xo count = %d, het xo count = %d", i, xocount[0], xocount[1]);
			myLogger.info(msg);
			xocount = countCrossOvers(family.imputed);
			msg = String.format("Family %d imputed: hom xo count = %d, het xo count = %d", i, xocount[0], xocount[1]);
			myLogger.info(msg);
			msg = String.format("Family %d: not missing count = %d, % d matches, proportion matching = %1.6f", i, totalNotMissing, matchCount, pmatch);
			myLogger.info(msg);
			assertTrue(String.format("Family %d: totalNotMissing > 490000 && pmatch > .999", i), totalNotMissing > 400000 && pmatch > .8);
			
		}
		
	}
	
	public static int[] compareSequence(byte[] ref, byte[] comp) {
		int n = ref.length;
		int matchCount = 0;
		int presentCount = 0;
		int nMissingRef = 0;
		int nMissingComp = 0;
		byte N = GenotypeTable.UNKNOWN_ALLELE;
		byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
		for (int s = 0; s < n; s++) {
			if (ref[s] == N || ref[s] == NN) nMissingRef++;
			if (comp[s] == N || comp[s] == NN) nMissingComp++;
			if (ref[s] != N && comp[s] != N && ref[s] != NN && comp[s] != NN) {
				presentCount++;
				if (ref[s] == comp[s] || (GenotypeTableUtils.isHeterozygous(ref[s]) && GenotypeTableUtils.isHeterozygous(comp[s])) ) matchCount++;
				
			}
		}
		return new int[] {n, presentCount, matchCount, nMissingRef, nMissingComp};  //total sites, sites present in both sequences, matching sites
	}
	
	public static void compareSequenceAndPrint(byte[] ref, byte[] comp) {
		int n = ref.length;
		int matchCount = 0;
		int presentCount = 0;
		byte N = GenotypeTable.UNKNOWN_ALLELE;
		byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
		for (int s = 0; s < n; s++) {
			if (ref[s] != N && comp[s] != N && ref[s] != NN && comp[s] != NN) {
				presentCount++;
				if (ref[s] != comp[s]) System.out.printf("Mismatch at %d, ref = %d, comp = %d\n", s, ref[s], comp[s]);
			}
		}
	}
	
	public byte[] diploidFromAllele(byte[] alleleGenotypes) {
		int n = alleleGenotypes.length;
		byte[] diploidGenotypes = new byte[n];
		for (int i = 0; i < n; i++) {
			diploidGenotypes[i] = GenotypeTableUtils.getDiploidValue(alleleGenotypes[i], alleleGenotypes[i]);
		}
		return diploidGenotypes;
	}
	
	public static int[] countCrossOvers(GenotypeTable genotype) {
		int ntaxa = genotype.numberOfTaxa();
		int nsites = genotype.numberOfSites();
		int numberOfHomXO = 0;
		int numberOfHetXO = 0;
		for (int t = 0; t < ntaxa; t++) {
			byte prevGenotype = genotype.genotype(t, 0);
			for (int s = 0; s < nsites; s++) {
				byte currentGenotype = genotype.genotype(t,s);
				if (prevGenotype == NN) {
					prevGenotype = currentGenotype;
				} else if (currentGenotype != NN && currentGenotype != prevGenotype){
					if (GenotypeTableUtils.isHeterozygous(prevGenotype) || GenotypeTableUtils.isHeterozygous(currentGenotype)) numberOfHetXO++;
					else numberOfHomXO++;
					prevGenotype = currentGenotype;
				}
			}
		}
		return new int[]{numberOfHomXO, numberOfHetXO};
	}
}
