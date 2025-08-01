package net.maizegenetics.analysis.association;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

import net.maizegenetics.analysis.numericaltransform.NumericalGenotypePlugin;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.GenotypeTable.GENOTYPE_TABLE_COMPONENT;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.GenotypePhenotypeBuilder;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.TableReport;

import org.junit.BeforeClass;
import org.junit.Test;

public class FastMultithreadedAssociationPluginTest {
    static GenotypeTable myGenotype;
    static GenotypeTable myNumericGenotype;
    static String phenotypeFile = "dataFiles/CandidateTests/mdp_traits_nomissing.txt";
    
    @BeforeClass
    public static void setUpBeforeClass() throws Exception {
        String genotypeFilename = "dataFiles/CandidateTests/mdp_genotype_chr9_10_nomissing.hmp.txt";
        myGenotype = ImportUtils.readFromHapmap(genotypeFilename);
        NumericalGenotypePlugin ngp = new NumericalGenotypePlugin(null, false);
        DataSet numgeno = ngp.processData(new DataSet(new Datum("name", myGenotype, ""), null));
        myNumericGenotype = (GenotypeTable) numgeno.getData(0).getData();

    }

    @Test
    public void testAdditive() {
        double maxP = 0.05;
        Phenotype initialPhenotype = new PhenotypeBuilder().fromFile(phenotypeFile).build().get(0);
        Phenotype oneTrait = new PhenotypeBuilder().fromPhenotype(initialPhenotype).keepAttributes(new int[]{0,1}).build().get(0);
        GenotypePhenotype myGenoPheno = new GenotypePhenotypeBuilder().genotype(myGenotype).phenotype(oneTrait).intersect().build();
        
        FastMultithreadedAssociationPlugin fastmap = new FastMultithreadedAssociationPlugin(null, false);
        fastmap.maxp(maxP);
        fastmap.saveAsFile(false);
        fastmap.maxThreads(3);
        DataSet inputData = new DataSet(new Datum("myGenoPheno", myGenoPheno, "comment"), null);
        DataSet resultData = fastmap.processData(inputData);

        TableReport resultTable = (TableReport) resultData.getData(0).getData();
        double[] expectedPvals = createValidationSetForAdditiveOnly(oneTrait, maxP);
        int nObsVals = (int) resultTable.getRowCount();
        double[] observedValues = new double[nObsVals];
        for (int i = 0; i < nObsVals; i++) {
            observedValues[i] = ((Double) resultTable.getValueAt(i, 6)).doubleValue();
        }
        Arrays.sort(observedValues);
        int nvals = Math.min(nObsVals, expectedPvals.length);

        for (int i = 0; i < nvals; i++) {
            String errMsg = String.format("Error at line %d testing array values:", i);
            double delta = expectedPvals[i] / 1e5;
            assertEquals(errMsg, expectedPvals[i], observedValues[i], delta);
        }
    }

    private double[] createValidationSetForAdditiveOnly(Phenotype myPheno, double maxP) {
        GenotypePhenotype myGenoPheno = new GenotypePhenotypeBuilder().genotype(myNumericGenotype).phenotype(myPheno).build();
        
        ReferenceProbabilityFELM rplm =
                new ReferenceProbabilityFELM(new Datum("genopheno", myGenoPheno, "genopheno"), null);

        rplm.maxP(maxP);

        rplm.solve();

        TableReport glmReport = rplm.siteReport();
        int nrows = (int) glmReport.getRowCount();
        double[] expectedPvals = IntStream.range(0, nrows)
                .mapToDouble(i -> ((Double) glmReport.getValueAt(i, 5)).doubleValue())
                .toArray();
        Arrays.sort(expectedPvals);
        return expectedPvals;
    }

    @Test
    public void testAdditiveWithCovariate() {
        double maxP = 0.1;
        Phenotype initialPhenotype = new PhenotypeBuilder().fromFile(phenotypeFile).build().get(0);
        List<PhenotypeAttribute> myAttr = new ArrayList<>();
        List<ATTRIBUTE_TYPE> myTypes = new ArrayList<>();
        for (int i = 0; i < 3; i++) myAttr.add(initialPhenotype.attribute(i));
        myTypes.add(ATTRIBUTE_TYPE.taxa);
        myTypes.add(ATTRIBUTE_TYPE.data);
        myTypes.add(ATTRIBUTE_TYPE.covariate);
        
        Phenotype oneTraitCovar = new PhenotypeBuilder().fromPhenotype(initialPhenotype).keepAttributes(new int[]{0,1}).build().get(0);
        GenotypePhenotype myGenoPheno = new GenotypePhenotypeBuilder().genotype(myGenotype).phenotype(oneTraitCovar).intersect().build();
        
        FastMultithreadedAssociationPlugin fastmap = new FastMultithreadedAssociationPlugin(null, false);
        fastmap.maxp(maxP);
        fastmap.saveAsFile(false);
        fastmap.maxThreads(3);

        DataSet inputData = new DataSet(new Datum("myGenoPheno", myGenoPheno, ""), null);

        DataSet resultData = fastmap.processData(inputData);

        TableReport resultTable = (TableReport) resultData.getData(0).getData();
        double[] expectedPvals = createValidationSetForAdditiveOnlyWithCovariate(oneTraitCovar, maxP);
        
        int nObsVals = (int) resultTable.getRowCount();
        double[] observedValues = new double[nObsVals];
        for (int i = 0; i < nObsVals; i++) {
            observedValues[i] = ((Double) resultTable.getValueAt(i, 6)).doubleValue();
        }
        Arrays.sort(observedValues);
        int nvals = Math.min(nObsVals, expectedPvals.length);

        for (int i = 0; i < nvals; i++) {
            String errMsg = String.format("Error at line %d testing array values:", i);
            double delta = expectedPvals[i] / 1e5;
            assertEquals(errMsg, expectedPvals[i], observedValues[i], delta);
        }

    }

    private double[] createValidationSetForAdditiveOnlyWithCovariate(Phenotype myPheno, double maxP) {
        
        GenotypePhenotype genoPheno = new GenotypePhenotypeBuilder().genotype(myNumericGenotype).phenotype(myPheno).intersect().build();
        ReferenceProbabilityFELM rplm =
                new ReferenceProbabilityFELM(new Datum("genopheno", genoPheno, "genopheno"), null);
        rplm.maxP(maxP);

        rplm.solve();

        TableReport glmReport = rplm.siteReport();
        int nrows = (int) glmReport.getRowCount();
        double[] expectedPvals = IntStream.range(0, nrows)
                .mapToDouble(i -> ((Double) glmReport.getValueAt(i, 5)).doubleValue())
                .toArray();
        Arrays.sort(expectedPvals);
        return expectedPvals;
    }
    
    @Test
    public void timingTest() {
        int ntraits = 10;
        double maxP = 0.01;
        
        //output file
//        String outfile = "/Volumes/Macintosh HD 2/temp/fast_assoc_test.txt";
        String filename = "dataFiles/CandidateTests/mdp_traits_nomissing.txt";
        Phenotype initialPhenotype = new PhenotypeBuilder().fromFile(filename).build().get(0);
        PhenotypeAttribute myTaxaAttr = initialPhenotype.taxaAttribute();
        NumericAttribute numericAttr =
                (NumericAttribute) initialPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data).get(0);
        List<PhenotypeAttribute> attrList = new ArrayList<>();
        List<ATTRIBUTE_TYPE> typeList = new ArrayList<>();
        attrList.add(initialPhenotype.taxaAttribute());
        typeList.add(ATTRIBUTE_TYPE.taxa);
        attrList.add(initialPhenotype.attribute(2));
        typeList.add(ATTRIBUTE_TYPE.covariate);

        float[] values = numericAttr.floatValues();
        for (int attr = 0; attr < ntraits; attr++) {
            String name = "trait" + attr;
            attrList.add(new NumericAttribute(name, values));
            typeList.add(ATTRIBUTE_TYPE.data);
        }

        Phenotype myPhenotype =
                new PhenotypeBuilder().fromAttributeList(attrList, typeList).build().get(0);
        GenotypePhenotype myGenoPheno = new GenotypePhenotypeBuilder().genotype(myGenotype)
        .phenotype(myPhenotype).intersect().build();
        
        FastMultithreadedAssociationPlugin fastmap = new FastMultithreadedAssociationPlugin(null, false);
        fastmap.maxp(maxP);
//        fastmap.saveAsFile(true);
//        fastmap.reportFilename(outfile);
        fastmap.maxThreads(4);
        DataSet inputData = new DataSet(new Datum("myGenoPheno", myGenoPheno, ""), null);
        int nsites = myGenotype.numberOfSites();
        int ntaxa = myGenotype.numberOfTaxa();
        long start = System.currentTimeMillis();
        
        DataSet ResultData = fastmap.processData(inputData);
        System.out.printf("%d traits analyzed in %d ms for %d sites and %d taxa.\n", ntraits, System.currentTimeMillis() - start, nsites, ntaxa);
    }

}
