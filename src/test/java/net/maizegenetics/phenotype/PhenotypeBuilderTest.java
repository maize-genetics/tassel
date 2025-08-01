package net.maizegenetics.phenotype;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

import org.junit.*;

public class PhenotypeBuilderTest {
	Phenotype pheno1;
	
	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		List<Taxon> taxaList = new ArrayList<Taxon>();
		taxaList.add(new Taxon("t1"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t3"));
		taxaList.add(new Taxon("t1"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t3"));
		
		TaxaAttribute ta = new TaxaAttribute(taxaList);
		CategoricalAttribute factor = new CategoricalAttribute("factor", new String[]{"A","A","A","B","B","B"});
		NumericAttribute trait = new NumericAttribute("trait1", new float[]{1,2,3,4,5,6}, new OpenBitSet(6));
		
		List<PhenotypeAttribute> attributes = new ArrayList<PhenotypeAttribute>();
		attributes.add(ta);
		attributes.add(factor);
		attributes.add(trait);
		List<ATTRIBUTE_TYPE> types = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		types.add(ATTRIBUTE_TYPE.taxa);
		types.add(ATTRIBUTE_TYPE.factor);
		types.add(ATTRIBUTE_TYPE.data);
		
		pheno1 = new CorePhenotype(attributes, types, "pheno1");
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void basicImportTest() {
		String testfileName = GeneralConstants.DATA_DIR + "PhenotypeTests/basic_phenotype_test.txt";
		PhenotypeBuilder builder = new PhenotypeBuilder();
		Phenotype testPhenotype = builder.fromFile(testfileName).build().get(0);
		float floatDelta = 1e-12f;
		
		assertEquals(5, testPhenotype.numberOfAttributes());
		assertEquals(7, testPhenotype.numberOfObservations());
		assertEquals(1, testPhenotype.numberOfAttributesOfType(ATTRIBUTE_TYPE.taxa));
		assertEquals(1, testPhenotype.numberOfAttributesOfType(ATTRIBUTE_TYPE.factor));
		assertEquals(1, testPhenotype.numberOfAttributesOfType(ATTRIBUTE_TYPE.covariate));
		assertEquals(2, testPhenotype.numberOfAttributesOfType(ATTRIBUTE_TYPE.data));
		
		assertEquals("Taxa", testPhenotype.attributeName(0));
		assertEquals("Factor1", testPhenotype.attributeName(1));
		assertEquals("Trait1", testPhenotype.attributeName(2));
		assertEquals("Trait2", testPhenotype.attributeName(3));
		assertEquals("Trait3", testPhenotype.attributeName(4));
		
		assertEquals(ATTRIBUTE_TYPE.taxa, testPhenotype.attributeType(0));
		assertEquals(ATTRIBUTE_TYPE.factor, testPhenotype.attributeType(1));
		assertEquals(ATTRIBUTE_TYPE.data, testPhenotype.attributeType(2));
		
		assertEquals(new Taxon("t1"), testPhenotype.value(0, 0));
		assertEquals(new Taxon("t4"), testPhenotype.value(3, 0));
		assertArrayEquals(new float[]{1.1f,2.1f,3.1f,Float.NaN,5.1f,6.1f,7.1f}, (float[]) testPhenotype.attribute(2).allValues(), floatDelta);
		assertArrayEquals(new float[]{1.2f,2.2f,3.2f,4.2f,5.2f,6.2f,7.2f}, ((NumericAttribute) testPhenotype.attribute(3)).floatValues(), floatDelta);
		
		//CategoricalAttribute tests
		int[] factorIndices = testPhenotype.attributeIndicesOfType(ATTRIBUTE_TYPE.factor);
		CategoricalAttribute myFactor = (CategoricalAttribute) testPhenotype.attribute(factorIndices[0]);
		assertEquals(0, myFactor.intValue(0));
		assertArrayEquals(new int[]{0,0,0,0,1,1,1}, myFactor.allIntValues());
		assertEquals("B", myFactor.attributeLabelForIndex(1));
		assertEquals(1, myFactor.indexForAttrLabel("B"));
		assertEquals("A", myFactor.label(1));
		assertArrayEquals(new String[]{"A","A","A","A","B","B","B"}, myFactor.allLabels());
		assertArrayEquals(new String[]{"A", "B"}, myFactor.labelList().toArray(new String[2]));
		
		PhenotypeAttribute myPheno = testPhenotype.attribute(factorIndices[0]);
		assertArrayEquals(new String[]{"A","A","A","A","B","B","B"}, (String[]) myPheno.allValues());
		assertEquals(ATTRIBUTE_TYPE.factor, myPheno.getCompatibleTypes().get(0));
		assertFalse(myPheno.isMissing(1));
		assertEquals("Factor1", myPheno.name());
		assertEquals(7, myPheno.size());
		PhenotypeAttribute mySubset = myPheno.subset(new int[]{0,4}, null);
		assertTrue(mySubset instanceof CategoricalAttribute);
		assertEquals(2, mySubset.size());
		assertArrayEquals(new String[]{"A","B"}, (String[]) mySubset.allValues());
		
		//TaxaAttribute tests
		TaxaAttribute myTaxa = (TaxaAttribute) testPhenotype.attribute(0);
		TaxaList tlist = new TaxaListBuilder().addAll(new String[]{"t1","t2","t3","t4"}).build();
		Taxon[] myTaxaArray = new Taxon[]{tlist.get(0), tlist.get(1), tlist.get(2), tlist.get(3), tlist.get(0), tlist.get(1), tlist.get(2)};
		assertArrayEquals(myTaxaArray, myTaxa.allTaxa());
		assertEquals("t2", myTaxa.taxon(1).getName());
		
		myPheno = testPhenotype.attribute(0);
		assertArrayEquals(myTaxaArray, (Taxon[]) myPheno.allValues());
		assertEquals(ATTRIBUTE_TYPE.taxa, myPheno.getCompatibleTypes().get(0));
		assertFalse(myPheno.isMissing(1));
		assertEquals("Taxa", myPheno.name());
		assertEquals(7, myPheno.size());
		mySubset = myPheno.subset(new int[]{0,4}, null);
		assertTrue(mySubset instanceof TaxaAttribute);
		assertEquals(2, mySubset.size());
		assertArrayEquals(new Taxon[]{tlist.get(0), tlist.get(0)}, (Taxon[]) mySubset.allValues());
		
		//NumericAttribute tests
		NumericAttribute myNumbers = (NumericAttribute) testPhenotype.attribute(2);
		assertEquals(6.1f, myNumbers.floatValue(5), 1e-12f);
		assertArrayEquals(new float[]{1.1f,2.1f,3.1f,Float.NaN,5.1f,6.1f,7.1f}, myNumbers.floatValues(), floatDelta);

		myPheno = testPhenotype.attribute(2);
		assertArrayEquals(new float[]{1.1f,2.1f,3.1f,Float.NaN,5.1f,6.1f,7.1f}, (float[]) myPheno.allValues(), floatDelta);
		assertTrue(myPheno.getCompatibleTypes().contains(ATTRIBUTE_TYPE.data));
		assertTrue(myPheno.getCompatibleTypes().contains(ATTRIBUTE_TYPE.covariate));
		
		assertTrue(myPheno.isMissing(3));
		OpenBitSet expectedMissing = new OpenBitSet(7);
		expectedMissing.set(3);
		assertEquals(expectedMissing, myPheno.missing());
		assertEquals("Trait1", myPheno.name());
		assertEquals(7, myPheno.size());
		mySubset = myPheno.subset(new int[]{0,4}, null);
		assertTrue(mySubset instanceof NumericAttribute);
		assertEquals(2, mySubset.size());
		assertArrayEquals(new float[]{1.1f, 5.1f}, (float[]) mySubset.allValues(), floatDelta);

	}
	
	@Test
	public void version4FormatImportTest() {
		float floatDelta = 1e-12f;
		String testfileName = GeneralConstants.DATA_DIR + "PhenotypeTests/version4format_test.txt";
		Phenotype testPhenotype = new PhenotypeBuilder().fromFile(testfileName).build().get(0);

		assertEquals(6, testPhenotype.numberOfObservations());
		assertEquals(4, testPhenotype.numberOfAttributes());
		TaxaAttribute myTaxaAttribute = testPhenotype.taxaAttribute();
		String[] taxaNames = new String[6];
		for (int i = 0; i < 6; i++) taxaNames[i] = myTaxaAttribute.taxon(i).getName();
		assertArrayEquals(new String[]{"taxon1", "taxon2", "taxon3", "taxon1", "taxon2", "taxon3"}, taxaNames);
		String[] levels = (String[]) testPhenotype.attribute(testPhenotype.attributeIndexForName("loc")).allValues();
		assertArrayEquals(new String[]{"A", "A","A","B","B","B"}, levels);
		float[] values = (float[]) testPhenotype.attribute(testPhenotype.attributeIndexForName("trait1")).allValues();
		assertArrayEquals(new float[]{11,12,13,14,15,16}, values, floatDelta);
		values = (float[]) testPhenotype.attribute(testPhenotype.attributeIndexForName("trait2")).allValues();
		assertArrayEquals(new float[]{21,22,23,24,25,26}, values, floatDelta);
		
		String covfileName = GeneralConstants.DATA_DIR + "PhenotypeTests/version4format_covariate.txt";
		Phenotype testCovPhenotype = new PhenotypeBuilder().fromFile(covfileName).build().get(0);
		
		Phenotype mergedPhenotype = new PhenotypeBuilder().fromPhenotype(testPhenotype).fromPhenotype(testCovPhenotype).intersectJoin().build().get(0);
		assertEquals(6, mergedPhenotype.numberOfObservations());
		assertEquals(5, mergedPhenotype.numberOfAttributes());
		values = (float[]) mergedPhenotype.attribute(mergedPhenotype.attributeIndexForName("Q1")).allValues();
		assertArrayEquals(new float[]{0.1f,0.1f,0.2f,0.2f,0.3f,0.3f}, values, floatDelta);
		myTaxaAttribute = mergedPhenotype.taxaAttribute();
		taxaNames = new String[6];
		for (int i = 0; i < 6; i++) taxaNames[i] = myTaxaAttribute.taxon(i).getName();
		assertArrayEquals(new String[]{"taxon1", "taxon1", "taxon2", "taxon2", "taxon3", "taxon3"}, taxaNames);
		levels = (String[]) mergedPhenotype.attribute(mergedPhenotype.attributeIndexForName("loc")).allValues();
		assertArrayEquals(new String[]{"A", "B","A","B","A","B"}, levels);

	}

	@Test
	public void mergeTest1() {
		
		ArrayList<Taxon> taxaList = new ArrayList<Taxon>();
		taxaList.add(new Taxon("t1"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t3"));
		
		TaxaAttribute ta = new TaxaAttribute(taxaList);
		CategoricalAttribute factor = new CategoricalAttribute("factor", new String[]{"C","C","C"});
		NumericAttribute trait = new NumericAttribute("trait1", new float[]{11,12,13}, new OpenBitSet(3));
		
		ArrayList<PhenotypeAttribute> attributes = new ArrayList<PhenotypeAttribute>();
		attributes.add(ta);
		attributes.add(factor);
		attributes.add(trait);
		ArrayList<Phenotype.ATTRIBUTE_TYPE> types = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		types.add(ATTRIBUTE_TYPE.taxa);
		types.add(ATTRIBUTE_TYPE.factor);
		types.add(ATTRIBUTE_TYPE.data);
		
		Phenotype pheno2 = new CorePhenotype(attributes, types, "pheno2");
		List<Phenotype> phenoList = new ArrayList<Phenotype>();
		phenoList.add(pheno1);
		phenoList.add(pheno2);
		Phenotype pheno12 = new PhenotypeBuilder().fromPhenotypeList(phenoList).intersectJoin().build().get(0);
		
		//test the result
		assertEquals(9, pheno12.numberOfObservations());
		assertEquals(3, pheno12.numberOfAttributes());
		Taxon[] taxa = pheno12.taxaAttribute().allTaxa();
		String[] taxaNames = new String[9];
		for (int i = 0; i < 9; i++) taxaNames[i] = taxa[i].getName();
		Arrays.sort(taxaNames);
		assertArrayEquals(new String[]{"t1", "t1", "t1", "t2", "t2", "t2", "t3", "t3", "t3"}, taxaNames);
		String[] factorLevels = (String[]) pheno12.attribute(pheno12.attributeIndexForName("factor")).allValues();
		assertArrayEquals(new String[]{"A","B","C","A","B","C","A","B","C"}, factorLevels);
		float[] traitvalues = (float[]) pheno12.attribute(pheno12.attributeIndexForName("trait1")).allValues();
		assertArrayEquals(new float[]{1,4,11,2,5,12,3,6,13}, traitvalues, 1e-12f);
		
		//test the separate function
		int factorIndex = pheno12.attributeIndexForName("factor");
		List<PhenotypeAttribute> factorlist = new ArrayList<PhenotypeAttribute>();
		factorlist.add(pheno12.attribute(factorIndex));
		List<Phenotype> separatePhenotypes = new PhenotypeBuilder().fromPhenotype(pheno12).separateOn(factorlist).build();
		assertEquals(separatePhenotypes.size(), 3);
		for (int i = 0; i < 3; i++) {
			Phenotype testPheno = separatePhenotypes.get(i);
			assertEquals(3, testPheno.numberOfObservations());
			assertEquals(2, testPheno.numberOfAttributes());
			taxa = testPheno.taxaAttribute().allTaxa();
			taxaNames = new String[3];
			for (int j = 0; j < 3; j++) taxaNames[j] = taxa[j].getName();
			Arrays.sort(taxaNames);
			assertArrayEquals(new String[]{"t1", "t2", "t3"}, taxaNames);
			if (i == 0) {
				traitvalues = (float[]) testPheno.attribute(testPheno.attributeIndexForName("trait1_factor.A")).allValues();
				assertArrayEquals(new float[]{1,2,3}, traitvalues, 1e-12f);
			}
			else if (i == 1) {
				traitvalues = (float[]) testPheno.attribute(testPheno.attributeIndexForName("trait1_factor.B")).allValues();
				assertArrayEquals(new float[]{4,5,6}, traitvalues, 1e-12f);
			}
			else if (i == 2) {
				traitvalues = (float[]) testPheno.attribute(testPheno.attributeIndexForName("trait1_factor.C")).allValues();
				assertArrayEquals(new float[]{11,12,13}, traitvalues, 1e-12f);
			}
		}
		
		//remerge the separated phenotypes
		Phenotype remerged = new PhenotypeBuilder().fromPhenotypeList(separatePhenotypes).intersectJoin().build().get(0);
		assertEquals(4, remerged.numberOfAttributes());
		assertEquals(3, remerged.numberOfObservations());
		assertEquals("trait1_factor.A", remerged.attribute(1).name());
		assertEquals("trait1_factor.B", remerged.attribute(2).name());
		assertEquals("trait1_factor.C", remerged.attribute(3).name());
	}

	@Test
	public void mergeTest1WithMultipleObservations() {
		
		ArrayList<Taxon> taxaList = new ArrayList<Taxon>();
		taxaList.add(new Taxon("t1"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t3"));
		
		TaxaAttribute ta = new TaxaAttribute(taxaList);
		CategoricalAttribute factor = new CategoricalAttribute("factor", new String[]{"C","C","C","C","C"});
		NumericAttribute trait = new NumericAttribute("trait1", new float[]{11,12,13,14,15}, new OpenBitSet(5));
		
		ArrayList<PhenotypeAttribute> attributes = new ArrayList<PhenotypeAttribute>();
		attributes.add(ta);
		attributes.add(factor);
		attributes.add(trait);
		ArrayList<Phenotype.ATTRIBUTE_TYPE> types = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		types.add(ATTRIBUTE_TYPE.taxa);
		types.add(ATTRIBUTE_TYPE.factor);
		types.add(ATTRIBUTE_TYPE.data);
		
		Phenotype pheno2 = new CorePhenotype(attributes, types, "pheno2");
		List<Phenotype> phenoList = new ArrayList<Phenotype>();
		phenoList.add(pheno1);
		phenoList.add(pheno2);
		Phenotype pheno12 = new PhenotypeBuilder().fromPhenotypeList(phenoList).intersectJoin().build().get(0);
		
		//test the result
		assertEquals(11, pheno12.numberOfObservations());
		assertEquals(3, pheno12.numberOfAttributes());
		Taxon[] taxa = pheno12.taxaAttribute().allTaxa();
		String[] taxaNames = new String[11];
		for (int i = 0; i < 11; i++) taxaNames[i] = taxa[i].getName();
		Arrays.sort(taxaNames);
		assertArrayEquals(new String[]{"t1", "t1", "t1", "t2", "t2", "t2", "t2", "t2", "t3", "t3", "t3"}, taxaNames);
		String[] factorLevels = (String[]) pheno12.attribute(pheno12.attributeIndexForName("factor")).allValues();
		assertArrayEquals(new String[]{"A","B","C","A","B","C","C","C","A","B","C"}, factorLevels);
		float[] traitvalues = (float[]) pheno12.attribute(pheno12.attributeIndexForName("trait1")).allValues();
		assertArrayEquals(new float[]{1,4,11,2,5,12,13,14,3,6,15}, traitvalues, 1e-12f);
		
		//test the separate function
		int factorIndex = pheno12.attributeIndexForName("factor");
		List<PhenotypeAttribute> factorlist = new ArrayList<PhenotypeAttribute>();
		factorlist.add(pheno12.attribute(factorIndex));
		List<Phenotype> separatePhenotypes = new PhenotypeBuilder().fromPhenotype(pheno12).separateOn(factorlist).build();
		assertEquals(separatePhenotypes.size(), 3);
		for (int i = 0; i < 3; i++) {
			Phenotype testPheno = separatePhenotypes.get(i);
			if (i == 2) {
				assertEquals(5, testPheno.numberOfObservations());
				assertEquals(2, testPheno.numberOfAttributes());
				taxa = testPheno.taxaAttribute().allTaxa();
				taxaNames = new String[5];
				for (int j = 0; j < 5; j++) taxaNames[j] = taxa[j].getName();
				Arrays.sort(taxaNames);
				assertArrayEquals(new String[]{"t1", "t2", "t2", "t2", "t3"}, taxaNames);
			}
			else {
				assertEquals(3, testPheno.numberOfObservations());
				assertEquals(2, testPheno.numberOfAttributes());
				taxa = testPheno.taxaAttribute().allTaxa();
				taxaNames = new String[3];
				for (int j = 0; j < 3; j++) taxaNames[j] = taxa[j].getName();
				Arrays.sort(taxaNames);
				assertArrayEquals(new String[]{"t1", "t2", "t3"}, taxaNames);
			}
			if (i == 0) {
				traitvalues = (float[]) testPheno.attribute(testPheno.attributeIndexForName("trait1_factor.A")).allValues();
				assertArrayEquals(new float[]{1,2,3}, traitvalues, 1e-12f);
			}
			else if (i == 1) {
				traitvalues = (float[]) testPheno.attribute(testPheno.attributeIndexForName("trait1_factor.B")).allValues();
				assertArrayEquals(new float[]{4,5,6}, traitvalues, 1e-12f);
			}
			else if (i == 2) {
				traitvalues = (float[]) testPheno.attribute(testPheno.attributeIndexForName("trait1_factor.C")).allValues();
				assertArrayEquals(new float[]{11,12,13,14,15}, traitvalues, 1e-12f);
			}
		}
		
		//remerge the separated phenotypes
		Phenotype remerged = new PhenotypeBuilder().fromPhenotypeList(separatePhenotypes).intersectJoin().build().get(0);
		assertEquals(4, remerged.numberOfAttributes());
		assertEquals(5, remerged.numberOfObservations());
		assertEquals("trait1_factor.A", remerged.attribute(1).name());
		assertEquals("trait1_factor.B", remerged.attribute(2).name());
		assertEquals("trait1_factor.C", remerged.attribute(3).name());
	}

	@Test
	public void mergeTest2() {
		List<Taxon> taxaList = new ArrayList<Taxon>();
		taxaList.add(new Taxon("t1"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t3"));
		
		TaxaAttribute ta = new TaxaAttribute(taxaList);
		NumericAttribute trait = new NumericAttribute("trait2", new float[]{11,12,13}, new OpenBitSet(3));
		
		ArrayList<PhenotypeAttribute> attributes = new ArrayList<PhenotypeAttribute>();
		attributes.add(ta);
		attributes.add(trait);
		ArrayList<Phenotype.ATTRIBUTE_TYPE> types = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		types.add(ATTRIBUTE_TYPE.taxa);
		types.add(ATTRIBUTE_TYPE.data);
		
		Phenotype pheno2 = new CorePhenotype(attributes, types, "pheno2");
		Phenotype pheno12 = new PhenotypeBuilder().fromPhenotype(pheno1).fromPhenotype(pheno2).intersectJoin().build().get(0);
		
		//test the result
		assertEquals(6, pheno12.numberOfObservations());
		assertEquals(4, pheno12.numberOfAttributes());
		Taxon[] taxa = pheno12.taxaAttribute().allTaxa();
		String[] taxaNames = new String[6];
		for (int i = 0; i < 6; i++) taxaNames[i] = taxa[i].getName();
		Arrays.sort(taxaNames);
		assertArrayEquals(new String[]{"t1", "t1", "t2", "t2", "t3", "t3"}, taxaNames);
		String[] factorLevels = (String[]) pheno12.attribute(pheno12.attributeIndexForName("factor")).allValues();
		assertArrayEquals(new String[]{"A","B","A","B","A","B"}, factorLevels);
		float[] traitvalues = (float[]) pheno12.attribute(pheno12.attributeIndexForName("trait1")).allValues();
		assertArrayEquals(new float[]{1,4,2,5,3,6}, traitvalues, 1e-12f);
		traitvalues = (float[]) pheno12.attribute(pheno12.attributeIndexForName("trait2")).allValues();
		assertArrayEquals(new float[]{11,11,12,12,13,13}, traitvalues, 1e-12f);
	}
	
	@Test
	public void mergeTest3() {
		List<Taxon> taxaList = new ArrayList<Taxon>();
		taxaList.add(new Taxon("t1"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t3"));
		taxaList.add(new Taxon("t1"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t3"));
		
		TaxaAttribute ta = new TaxaAttribute(taxaList);
		CategoricalAttribute factor = new CategoricalAttribute("factor", new String[]{"A","A","A","B","B","B"});
		NumericAttribute trait = new NumericAttribute("trait1", new float[]{7,8,9,10,11,12}, new OpenBitSet(6));
		
		List<PhenotypeAttribute> attributes = new ArrayList<PhenotypeAttribute>();
		attributes.add(ta);
		attributes.add(factor);
		attributes.add(trait);
		List<ATTRIBUTE_TYPE> types = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		types.add(ATTRIBUTE_TYPE.taxa);
		types.add(ATTRIBUTE_TYPE.factor);
		types.add(ATTRIBUTE_TYPE.data);
		Phenotype pheno2 = new CorePhenotype(attributes, types, "pheno2");
		
		try {
			Phenotype pheno12 = new PhenotypeBuilder().fromPhenotype(pheno1).fromPhenotype(pheno2).intersectJoin().build().get(0);
			assertTrue("PhenotypeBuilder should have thrown an error but did not", false);
		} catch (Exception e) {
			assertEquals("Data sets will not be joined because both phenotypes have values for trait1", e.getMessage());
		}
	}
	
	@Test
	public void mergeTest4() {
		List<Taxon> taxaList = new ArrayList<Taxon>();
		taxaList.add(new Taxon("t1"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t3"));
		taxaList.add(new Taxon("t1"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t3"));
		
		TaxaAttribute ta = new TaxaAttribute(taxaList);
		CategoricalAttribute factor = new CategoricalAttribute("factor", new String[]{"A","A","A","B","B","B"});
		NumericAttribute trait = new NumericAttribute("trait2", new float[]{7,8,9,10,11,12}, new OpenBitSet(6));
		
		List<PhenotypeAttribute> attributes = new ArrayList<PhenotypeAttribute>();
		attributes.add(ta);
		attributes.add(factor);
		attributes.add(trait);
		List<ATTRIBUTE_TYPE> types = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		types.add(ATTRIBUTE_TYPE.taxa);
		types.add(ATTRIBUTE_TYPE.factor);
		types.add(ATTRIBUTE_TYPE.data);
		Phenotype pheno2 = new CorePhenotype(attributes, types, "pheno2");
		Phenotype pheno12 = new PhenotypeBuilder().fromPhenotype(pheno1).fromPhenotype(pheno2).intersectJoin().build().get(0);
		
		float[] data = (float[]) pheno12.attribute(pheno12.attributeIndexForName("trait1")).allValues();
		assertArrayEquals(new float[]{1,4,2,5,3,6}, data, 1e-12f);
		data = (float[]) pheno12.attribute(pheno12.attributeIndexForName("trait2")).allValues();
		assertArrayEquals(new float[]{7,10,8,11,9,12}, data, 1e-12f);
		Taxon[] taxa = pheno12.taxaAttribute().allTaxa();
		Collections.sort(taxaList);
		assertArrayEquals(taxaList.toArray(new Taxon[0]), taxa);
		String[] levels = (String[]) pheno12.attribute(pheno12.attributeIndexForName("factor")).allValues();
		assertArrayEquals(new String[]{"A","B","A","B","A","B"}, levels);
	}
	
	@Test
	public void mergeTest4WithMultipleObservations() {
		List<Taxon> taxaList = new ArrayList<Taxon>();
		taxaList.add(new Taxon("t1"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t3"));
		taxaList.add(new Taxon("t1"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t3"));
		taxaList.add(new Taxon("t3"));
		
		TaxaAttribute ta = new TaxaAttribute(taxaList);
		CategoricalAttribute factor = new CategoricalAttribute("factor", new String[]{"A","A","A","A","B","B","B","B"});
		NumericAttribute trait = new NumericAttribute("trait2", new float[]{7,8,9,10,11,12,13,14}, new OpenBitSet(8));
		
		List<PhenotypeAttribute> attributes = new ArrayList<PhenotypeAttribute>();
		attributes.add(ta);
		attributes.add(factor);
		attributes.add(trait);
		List<ATTRIBUTE_TYPE> types = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		types.add(ATTRIBUTE_TYPE.taxa);
		types.add(ATTRIBUTE_TYPE.factor);
		types.add(ATTRIBUTE_TYPE.data);
		Phenotype pheno2 = new CorePhenotype(attributes, types, "pheno2");
		Phenotype pheno12 = new PhenotypeBuilder().fromPhenotype(pheno1).fromPhenotype(pheno2).intersectJoin().build().get(0);
		
		float[] data = (float[]) pheno12.attribute(pheno12.attributeIndexForName("trait1")).allValues();
		assertArrayEquals(new float[]{1,4,2,2,5,3,6,6}, data, 1e-12f);
		data = (float[]) pheno12.attribute(pheno12.attributeIndexForName("trait2")).allValues();
		assertArrayEquals(new float[]{7,11,8,9,12,10,13,14}, data, 1e-12f);
		Taxon[] taxa = pheno12.taxaAttribute().allTaxa();
		Collections.sort(taxaList);
		assertArrayEquals(taxaList.toArray(new Taxon[0]), taxa);
		String[] levels = (String[]) pheno12.attribute(pheno12.attributeIndexForName("factor")).allValues();
		assertArrayEquals(new String[]{"A","B","A","A","B","A","B","B"}, levels);
	}

	@Test
	public void filterTaxaTest() {
		List<Taxon> taxaToRemove = new ArrayList<Taxon>();
		taxaToRemove.add(new Taxon("t2"));
		List<Taxon> taxaToKeep = new ArrayList<Taxon>();
		taxaToKeep.add(new Taxon("t1"));
		taxaToKeep.add(new Taxon("t3"));
		
		Phenotype filtered = new PhenotypeBuilder().fromPhenotype(pheno1).removeTaxa(taxaToRemove).build().get(0);
		testTaxaFilter(filtered);
		filtered = new PhenotypeBuilder().fromPhenotype(pheno1).keepTaxa(taxaToKeep).build().get(0);
		testTaxaFilter(filtered);
		
		taxaToKeep.remove(0);
		filtered = new PhenotypeBuilder().fromPhenotype(filtered).keepTaxa(taxaToKeep).build().get(0);

		assertEquals(2, filtered.numberOfObservations());
		assertEquals(3, filtered.numberOfAttributes());
		Taxon[] taxa = filtered.taxaAttribute().allTaxa();
		String[] taxaNames = new String[2];
		for (int i = 0; i < 2; i++) taxaNames[i] = taxa[i].getName();
		assertArrayEquals(new String[]{"t3", "t3"}, taxaNames);
		String[] factorLevels = (String[]) filtered.attribute(filtered.attributeIndexForName("factor")).allValues();
		assertArrayEquals(new String[]{"A","B"}, factorLevels);
		float[] traitvalues = (float[]) filtered.attribute(filtered.attributeIndexForName("trait1")).allValues();
		assertArrayEquals(new float[]{3,6}, traitvalues, 1e-12f);

	}
	
	@Test
	public void attributeFilterTest() {
		List<Taxon> taxaList = new ArrayList<Taxon>();
		taxaList.add(new Taxon("t1"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t3"));
		taxaList.add(new Taxon("t1"));
		taxaList.add(new Taxon("t2"));
		taxaList.add(new Taxon("t3"));
		
		TaxaAttribute ta = new TaxaAttribute(taxaList);
		CategoricalAttribute factor = new CategoricalAttribute("factor", new String[]{"A","A","A","B","B","B"});
		NumericAttribute trait = new NumericAttribute("trait2", new float[]{7,8,9,10,11,12}, new OpenBitSet(6));
		
		List<PhenotypeAttribute> attributes = new ArrayList<PhenotypeAttribute>();
		attributes.add(ta);
		attributes.add(factor);
		attributes.add(trait);
		List<ATTRIBUTE_TYPE> types = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		types.add(ATTRIBUTE_TYPE.taxa);
		types.add(ATTRIBUTE_TYPE.factor);
		types.add(ATTRIBUTE_TYPE.data);
		Phenotype pheno2 = new CorePhenotype(attributes, types, "pheno2");
		Phenotype pheno12 = new PhenotypeBuilder().fromPhenotype(pheno1).fromPhenotype(pheno2).intersectJoin().build().get(0);
		
		int[] index = new int[]{0,1,3};
		Phenotype pheno3 = new PhenotypeBuilder().fromPhenotype(pheno12).keepAttributes(index).build().get(0);
		testAttributeFilter(pheno3);
		
		ArrayList<PhenotypeAttribute> attributesToKeep = new ArrayList<>();
		attributesToKeep.add(pheno12.attribute(pheno12.attributeIndexForName(TaxaAttribute.DEFAULT_NAME)));
		attributesToKeep.add(pheno12.attribute(pheno12.attributeIndexForName("factor")));
		attributesToKeep.add(pheno12.attribute(pheno12.attributeIndexForName("trait2")));
		pheno3 = new PhenotypeBuilder().fromPhenotype(pheno12).keepAttributes(attributesToKeep).build().get(0);
		testAttributeFilter(pheno3);
		
		
	}

	@Test
	public void testPhenotypeImportErrors() {
		//Are there the same number of taxa names and data types?
		try {
			File tmpFile = File.createTempFile("phenotest", "txt");
			PrintWriter pw = new PrintWriter(tmpFile);
			pw.println("<Phenotype>");
			pw.println("taxa\tcovariate\tcovariate\tdata");
			pw.println("Taxa\tq1\tq2\ttrt1\ttrt2");
			pw.println("T1\t1\t2\t3.0\t4.0");
			pw.println("T2\t5\t6\t7.0\t8.0");
			pw.close();

			new PhenotypeBuilder().fromFile(tmpFile.toString()).build();
			Assert.fail("Exception not thrown for number of datatypes != number of phenotype names");
		} catch (IOException ioe) {
			throw new IllegalArgumentException("IO Error");
		} catch (IllegalArgumentException e) {
			String err = e.getMessage();
			System.out.println(err+"\n");
			assertTrue(err.endsWith("The number of data types must equal the number of phenotype names."));
		}

		//Do all rows have the correct number of data values?
		try {
			File tmpFile = File.createTempFile("phenotest", "txt");
			PrintWriter pw = new PrintWriter(tmpFile);
			pw.println("<Phenotype>");
			pw.println("taxa\tcovariate\tcovariate\tdata\tdata");
			pw.println("Taxa\tq1\tq2\ttrt1\ttrt2");
			pw.println("T1\t1\t2\t3.0\t4.0");
			pw.println("T2\t5\t6\t7.0\t8.0\t5.0");
			pw.close();

			new PhenotypeBuilder().fromFile(tmpFile.toString()).build();
			Assert.fail("Exception not thrown for number of data values != number of data types");
		} catch (IOException ioe) {
			throw new IllegalArgumentException("IO Error");
		} catch (IllegalArgumentException e) {
			String err = e.getMessage();
			System.out.println(err+"\n");
			assertTrue(err.endsWith("data types in the header."));
		}

		//Is there a taxa column?
		try {
			File tmpFile = File.createTempFile("phenotest", "txt");
			PrintWriter pw = new PrintWriter(tmpFile);
			pw.println("<Phenotype>");
			pw.println("factor\tcovariate\tcovariate\tdata\tdata");
			pw.println("Taxa\tq1\tq2\ttrt1\ttrt2");
			pw.println("T1\t1\t2\t3.0\t4.0");
			pw.println("T2\t5\t6\t7.0\t8.0");
			pw.close();

			new PhenotypeBuilder().fromFile(tmpFile.toString()).build();
			Assert.fail("Exception not thrown for data types do not include taxa");
		} catch (IOException ioe) {
			throw new IllegalArgumentException("IO Error");
		} catch (IllegalArgumentException e) {
			String err = e.getMessage();
			System.out.println(err+"\n");
			assertTrue(err.endsWith("data types must contain one taxa column"));
		}

	}

	private void testTaxaFilter(Phenotype filtered) {
		assertEquals(4, filtered.numberOfObservations());
		assertEquals(3, filtered.numberOfAttributes());
		Taxon[] taxa = filtered.taxaAttribute().allTaxa();
		String[] taxaNames = new String[4];
		for (int i = 0; i < 4; i++) taxaNames[i] = taxa[i].getName();
		assertArrayEquals(new String[]{"t1", "t3", "t1", "t3"}, taxaNames);
		String[] factorLevels = (String[]) filtered.attribute(filtered.attributeIndexForName("factor")).allValues();
		assertArrayEquals(new String[]{"A","A","B","B"}, factorLevels);
		float[] traitvalues = (float[]) filtered.attribute(filtered.attributeIndexForName("trait1")).allValues();
		assertArrayEquals(new float[]{1,3,4,6}, traitvalues, 1e-12f);
	}
	
	private void testAttributeFilter(Phenotype filtered) {
		assertEquals(6, filtered.numberOfObservations());
		assertEquals(3, filtered.numberOfAttributes());
		Taxon[] taxa = filtered.taxaAttribute().allTaxa();
		String[] taxaNames = new String[6];
		for (int i = 0; i < 6; i++) taxaNames[i] = taxa[i].getName();
		assertArrayEquals(new String[]{"t1", "t1", "t2", "t2", "t3", "t3"}, taxaNames);
		String[] factorLevels = (String[]) filtered.attribute(filtered.attributeIndexForName("factor")).allValues();
		assertArrayEquals(new String[]{"A","B","A","B","A","B"}, factorLevels);
		float[] traitvalues = (float[]) filtered.attribute(filtered.attributeIndexForName("trait2")).allValues();
		assertArrayEquals(new float[]{7,10,8,11,9,12}, traitvalues, 1e-12f);
	}
}
