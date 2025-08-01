package net.maizegenetics.analysis.association;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

import net.maizegenetics.analysis.numericaltransform.ImputationPlugin;
import net.maizegenetics.analysis.numericaltransform.NumericalGenotypePlugin;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable.GENOTYPE_TABLE_COMPONENT;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.GenotypePhenotypeBuilder;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.TableReport;

import org.junit.BeforeClass;
import org.junit.Test;

public class EqtlAssociationPluginTest {
    private static GenotypePhenotype myGenoPheno;
    private static GenotypePhenotype myGenoPhenoCov;
    private static GenotypePhenotype myGenoPhenoUnion;
    private static GenotypeTable myNumericGenotype;
    private static double maxP = .01;
    private static double[] expectedPvals;

    @BeforeClass
    public static void setUpBeforeClass() throws Exception {
        String genotypeFilename = "dataFiles/CandidateTests/mdp_genotype_chr9_10_nomissing.hmp.txt";
        GenotypeTable myGenotype = ImportUtils.readFromHapmap(genotypeFilename);
        NumericalGenotypePlugin ngp = new NumericalGenotypePlugin(null, false);
        DataSet numgeno = ngp.processData(new DataSet(new Datum("name", myGenotype, ""), null));
        myNumericGenotype = (GenotypeTable) numgeno.getData(0).getData();
        String filename = "dataFiles/CandidateTests/mdp_traits_nomissing.txt";
        Phenotype myPhenotype = new PhenotypeBuilder().fromFile(filename).build().get(0);
        myGenoPheno =
                new GenotypePhenotypeBuilder().genotype(myNumericGenotype)
                        .phenotype(myPhenotype).intersect().build();
        myGenoPhenoUnion =
                new GenotypePhenotypeBuilder().genotype(myNumericGenotype)
                        .phenotype(myPhenotype).union().build();
        Map<PhenotypeAttribute, ATTRIBUTE_TYPE> attrMap = new HashMap<>();
        attrMap.put(myPhenotype.attribute(2), ATTRIBUTE_TYPE.covariate);
        Phenotype phenoWithCovar = new PhenotypeBuilder().fromPhenotype(myPhenotype)
                .changeAttributeType(attrMap)
                .keepAttributes(new int[] { 0, 1, 2 })
                .build().get(0);
        myGenoPhenoCov =
                new GenotypePhenotypeBuilder().genotype(myNumericGenotype)
                        .phenotype(phenoWithCovar).intersect().build();

    }

    @Test
    public void testAdditive() {
        maxP = 0.01;
        EqtlAssociationPlugin eap = new EqtlAssociationPlugin(null, false);
        eap.additiveOnlyModel(true);
        eap.genotypeComponent(GENOTYPE_TABLE_COMPONENT.ReferenceProbability);
        eap.maxPValue(maxP);
        DataSet inputData = new DataSet(new Datum("myGenoPheno", myGenoPheno, ""), null);

        DataSet resultData = eap.processData(inputData);

        TableReport resultTable = (TableReport) resultData.getData(0).getData();
        double[] expectedPvals = createValidationSetForAdditiveOnly();
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

        DataSet unionData = new DataSet(new Datum("myGenoPhenoUnion", myGenoPhenoUnion, ""), null);
        DataSet unionResult = eap.processData(unionData);
        TableReport unionTable = (TableReport) unionResult.getData(0).getData();
        for (int i = 0; i < nObsVals; i++) {
            observedValues[i] = ((Double) unionTable.getValueAt(i, 6)).doubleValue();
        }
        Arrays.sort(observedValues);

        expectedPvals = createValidationSetForAdditiveOnly();
        nObsVals = (int) resultTable.getRowCount();
        nvals = Math.min(nObsVals, expectedPvals.length);

        //for now do not test. results will generally not be equal to GLM because missing values are handled differently
        //the above just tests that the union join does not throw an error
    }

    private double[] createValidationSetForAdditiveOnly() {

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
    public void testAdditivePlusDominance() {
        maxP = 0.05;
        EqtlAssociationPlugin eap = new EqtlAssociationPlugin(null, false);
        eap.additiveOnlyModel(false);
        eap.genotypeComponent(GENOTYPE_TABLE_COMPONENT.Genotype);
        eap.maxPValue(maxP);
        DataSet inputData = new DataSet(new Datum("myGenoPheno", myGenoPheno, ""), null);
        DataSet resultData = eap.processData(inputData);
        TableReport resultTable = (TableReport) resultData.getData(0).getData();
        double[] expectedPvals = createValidationSetForAD();
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
            //debug
            if (Math.abs(expectedPvals[i] - observedValues[i]) > delta) {
                System.out.println();
            }
            assertEquals(errMsg, expectedPvals[i], observedValues[i], delta);
        }
    }

    private double[] createValidationSetForAD() {
        DiscreteSitesFELM dslm =
                new DiscreteSitesFELM(new Datum("genopheno", myGenoPheno, "genopheno"), null);
        dslm.maxP(maxP);
        dslm.solve();
        TableReport glmReport = dslm.siteReport();
        int nrows = (int) glmReport.getRowCount();
        double[] expectedPvals = IntStream.range(0, nrows)
                .mapToDouble(i -> ((Double) glmReport.getValueAt(i, 5)).doubleValue())
                .toArray();
        Arrays.sort(expectedPvals);
        return expectedPvals;
    }

    @Test
    public void testAdditiveWithCovariate() {
        maxP = 0.05;
        EqtlAssociationPlugin eap = new EqtlAssociationPlugin(null, false);
        eap.additiveOnlyModel(true);
        eap.genotypeComponent(GENOTYPE_TABLE_COMPONENT.ReferenceProbability);
        eap.maxPValue(maxP);
        DataSet inputData = new DataSet(new Datum("myGenoPheno", myGenoPhenoCov, ""), null);

        DataSet resultData = eap.processData(inputData);

        TableReport resultTable = (TableReport) resultData.getData(0).getData();
        double[] expectedPvals = createValidationSetForAdditiveOnlyWithCovariate();
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

    private double[] createValidationSetForAdditiveOnlyWithCovariate() {
        ReferenceProbabilityFELM rplm =
                new ReferenceProbabilityFELM(new Datum("genopheno", myGenoPhenoCov, "genopheno"), null);
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
    public void testAdditivePlusDominanceWithCovariate() {
        maxP = 0.05;
        EqtlAssociationPlugin eap = new EqtlAssociationPlugin(null, false);
        eap.additiveOnlyModel(false);
        eap.genotypeComponent(GENOTYPE_TABLE_COMPONENT.Genotype);
        eap.maxPValue(maxP);
        DataSet inputData = new DataSet(new Datum("myGenoPheno", myGenoPhenoCov, ""), null);
        DataSet resultData = eap.processData(inputData);
        TableReport resultTable = (TableReport) resultData.getData(0).getData();
        double[] expectedPvals = createValidationSetForADWithCovar();
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

    private double[] createValidationSetForADWithCovar() {
        DiscreteSitesFELM dslm =
                new DiscreteSitesFELM(new Datum("genopheno", myGenoPhenoCov, "genopheno"), null);
        dslm.maxP(maxP);
        dslm.solve();
        TableReport glmReport = dslm.siteReport();
        int nrows = (int) glmReport.getRowCount();
        double[] expectedPvals = IntStream.range(0, nrows)
                .mapToDouble(i -> ((Double) glmReport.getValueAt(i, 5)).doubleValue())
                .toArray();
        Arrays.sort(expectedPvals);
        return expectedPvals;
    }

    @Test
    public void timingTestAdditiveWithCovariate() {
        int ntraits = 10;
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
        GenotypePhenotype myGenoPheno = new GenotypePhenotypeBuilder().genotype(myNumericGenotype)
        .phenotype(myPhenotype).intersect().build();
        
        maxP = 0.01;
        EqtlAssociationPlugin eap = new EqtlAssociationPlugin(null, false);
        eap.additiveOnlyModel(true);
        eap.genotypeComponent(GENOTYPE_TABLE_COMPONENT.ReferenceProbability);
        eap.maxPValue(maxP);
        DataSet inputData = new DataSet(new Datum("myGenoPheno", myGenoPheno, ""), null);
        
        int nsites = myNumericGenotype.numberOfSites();
        int ntaxa = myNumericGenotype.numberOfTaxa();
        long start = System.currentTimeMillis();
        
        DataSet resultData = eap.processData(inputData);
        System.out.printf("%d traits analyzed in %d ms for %d sites and %d taxa.\n", ntraits, System.currentTimeMillis() - start, nsites, ntaxa);
        TableReport resultTable = (TableReport) resultData.getData(0).getData();
        
        
    }

}
