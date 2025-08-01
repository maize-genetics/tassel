package net.maizegenetics.analysis.association;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javax.swing.ImageIcon;

import net.maizegenetics.stats.linearmodels.*;
import org.apache.commons.math3.distribution.FDistribution;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable.GENOTYPE_TABLE_COMPONENT;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.stats.linearmodels.SolveByOrthogonalizing.Marker;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportBuilder;

import org.apache.commons.math3.exception.OutOfRangeException;

public class EqtlAssociationPlugin extends AbstractPlugin {
    private GENOTYPE_TABLE_COMPONENT[] GENOTYPE_COMP = new GENOTYPE_TABLE_COMPONENT[] {
            GENOTYPE_TABLE_COMPONENT.Genotype, GENOTYPE_TABLE_COMPONENT.ReferenceProbability,
            GENOTYPE_TABLE_COMPONENT.AlleleProbability };

    private Datum myDatum;
    private GenotypePhenotype myGenoPheno;
    private GenotypeTable myGenotype;
    private Phenotype myPhenotype;
    private TableReportBuilder myReportBuilder;
    private SolveByOrthogonalizing orthogonalSolver;
    private List<String> phenotypeNames;
    private double minR2[];
    private FDistribution[] Fdist;
    private int numberOfObservations;

    //plugin parameter definitions
    private PluginParameter<Double> maxp =
            new PluginParameter.Builder<>("MaxPValue", .001, Double.class)
                    .guiName("MaxPValue")
                    .description("The maximum p-value that will be output by the analysis.")
                    .build();
    private PluginParameter<Boolean> addOnly =
            new PluginParameter.Builder<>("addOnly", false, Boolean.class)
                    .description("Should an additive only model be fit? If true, an additive model will be fit. If false, an additive + dominance model will be fit. Default = false.")
                    .guiName("Additive Only Model")
                    .build();
    private PluginParameter<GENOTYPE_TABLE_COMPONENT> myGenotypeTable =
            new PluginParameter.Builder<>("genotypeComponent", GENOTYPE_TABLE_COMPONENT.Genotype, GENOTYPE_TABLE_COMPONENT.class)
                    .genotypeTable()
                    .range(GENOTYPE_COMP)
                    .description("If the genotype table contains more than one type of genotype data, choose the type to use for the analysis.")
                    .build();
    private PluginParameter<Boolean> saveAsFile =
            new PluginParameter.Builder<>("writeToFile", false, Boolean.class)
                    .description("Should the results be saved to a file rather than stored in memory? It true, the results will be written to a file as each SNP is analyzed in order to reduce memory requirements"
                            + "and the results will NOT be saved to the data tree. Default = false.")
                    .guiName("Write to file")
                    .build();
    private PluginParameter<String> reportFilename =
            new PluginParameter.Builder<>("outputFile", null, String.class)
                    .outFile()
                    .dependentOnParameter(saveAsFile)
                    .description("The name of the file to which these results will be saved.")
                    .guiName("Output File")
                    .build();

    public EqtlAssociationPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet processData(DataSet input) {
        long start = System.currentTimeMillis();
        List<Datum> datumList = input.getDataOfType(GenotypePhenotype.class);
        if (datumList.size() != 1)
            throw new IllegalArgumentException("Exactly one data set with genotypes and phenotypes must be selected.");
        myDatum = datumList.get(0);
        myGenoPheno = (GenotypePhenotype) myDatum.getData();
        myGenotype = myGenoPheno.genotypeTable();
        myPhenotype = myGenoPheno.phenotype();
        numberOfObservations = myPhenotype.numberOfObservations();
        testMissingDataInTheBaseModel();
        initializeOutput();
        initializeOrthogonalizer();
        final int nsites = myGenotype.numberOfSites();
        Fdist = new FDistribution[2];
        Fdist[0] = new FDistribution(1, numberOfObservations - 1 - orthogonalSolver.baseDf());
        Fdist[1] = new FDistribution(2, numberOfObservations - 2 - orthogonalSolver.baseDf());
        calculateR2Fromp();
        System.out.printf("starting site loop at %d.\n", System.currentTimeMillis() - start);
        start = System.currentTimeMillis();
        if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype) {
            if (addOnly.value()) {
                IntStream.range(0, nsites)
                        .mapToObj(s -> orthogonalSolver.solveForR(myGenotype.positions().get(s), additiveSite(s)))
                        .forEach(this::updateOutputWithPvalues);
            } else {
                IntStream.range(0, nsites)
                        .mapToObj(s -> {
                            List<double[]> covars = additiveDominanceSite(s);
                            return orthogonalSolver.solveForR(myGenotype.positions().get(s), covars.get(0), covars.get(1));
                        })
                        .forEach(this::updateOutputWithPvalues);
            }
        } else if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability) {
            //    		IntStream.range(0, nsites)
            //    			.mapToObj(s-> orthogonalSolver.solveForR(myGenotype.positions().get(s), referenceProbabilitiesForSite(s)))
            //    			.forEach(this::updateOutputWithPvalues);
            for (int s = 0; s < nsites; s++) {
                Marker markerTest =
                        orthogonalSolver.solveForR(myGenotype.positions().get(s), referenceProbabilitiesForSite(s));
                updateOutputWithPvalues(markerTest);
            }

        } else if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability) {
            throw new UnsupportedOperationException("Fast association analysis of allele probabilities is not supported.");
            //TODO implement
        }
        System.out.printf("site loop finished after %d ms\n", System.currentTimeMillis() - start);
        if (!saveAsFile.value()) {
            String name = "FastAssociation_" + myDatum.getName();
            String comment = "Fast association output";
            Datum outDatum = new Datum(name, myReportBuilder.build(), comment);
            return new DataSet(outDatum, this);
        } else {
            myReportBuilder.build();
        }

        return null;
    }

    private void initializeOutput() {
        //output is a TableReport with p-value; site position information: chr, position, id; trait name
        //add separate values for additive test and dominant test later
        String[] columnNames = new String[] { AssociationConstants.STATS_HEADER_TRAIT,
                AssociationConstants.STATS_HEADER_MARKER,
                AssociationConstants.STATS_HEADER_CHR,
                AssociationConstants.STATS_HEADER_POSITION,
                "df",
                "r2",
                AssociationConstants.STATS_HEADER_P_VALUE };
        String name = "EqtlReport_" + myDatum.getName();
        if (saveAsFile.value())
            myReportBuilder =
                    TableReportBuilder.getInstance(name, columnNames, reportFilename.value());
        else
            myReportBuilder = TableReportBuilder.getInstance(name, columnNames);
    }

    private void updateOutputWithPvalues(SolveByOrthogonalizing.Marker markerResult) {
        double[] rvalues = markerResult.vector1();
        int npheno = rvalues.length;
        Position pos = markerResult.position();

        //debug
        if (markerResult.df > 0) {
            final double minrsq = minR2[markerResult.df - 1];
            int errdf = numberOfObservations - orthogonalSolver.baseDf() - markerResult.df;

            IntStream.range(0, npheno).sequential().filter(i -> rvalues[i] >= minrsq)
                    .forEach(i -> addToReport(new Object[] { phenotypeNames.get(i), pos.getSNPID(),
                            pos.getChromosome().getName(), pos.getPosition(),
                            markerResult.degreesOfFreedom(), rvalues[i],
                            pvalue(rvalues[i], markerResult.df, errdf) }));
        }
    }

    private double pvalue(double rvalue, int markerdf, int errordf) {
        double F = rvalue / (1 - rvalue) * errordf / markerdf;
        double p;
        try {
            p = LinearModelUtils.Ftest(F, markerdf, errordf);
        } catch (Exception e) {
            p = Double.NaN;
        }
        return p;

    }

    private void addToReport(Object[] row) {
        myReportBuilder.add(row);
    }

    private void initializeOrthogonalizer() {
        List<PhenotypeAttribute> phenotypeList =
                myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data);
        List<PhenotypeAttribute> covariateList =
                myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.covariate);
        List<PhenotypeAttribute> factorList =
                myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor);

        //build the model, no mean necessary because it will not be used
        List<ModelEffect> baseModel = new ArrayList<>();
        for (PhenotypeAttribute pa : factorList) {
            CategoricalAttribute ca = (CategoricalAttribute) pa;
            baseModel.add(new FactorModelEffect(ca.allIntValues(), true, ca.name()));
        }
        for (PhenotypeAttribute pa : covariateList) {
            NumericAttribute na = (NumericAttribute) pa;
            CovariateModelEffect cme =
                    new CovariateModelEffect(AssociationUtils.convertFloatArrayToDouble(na.floatValues()), na.name());
            baseModel.add(cme);
        }

        List<double[]> dataList = phenotypeList.stream()
                .map(pa -> (float[]) pa.allValues())
                .map(a -> AssociationUtils.convertFloatArrayToDouble(a))
                .collect(Collectors.toList());

        phenotypeNames =
                phenotypeList.stream().map(PhenotypeAttribute::name).collect(Collectors.toList());

        orthogonalSolver = SolveByOrthogonalizing.getInstanceFromModel(baseModel, dataList);
    }

    private void testMissingDataInTheBaseModel() {
        for (PhenotypeAttribute attr : myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor)) {
            if (attr.missing().cardinality() > 0) {
                String msg = "There is missing data in the factor " + attr.name();
                throw new IllegalArgumentException(msg);
            }
        }
        for (PhenotypeAttribute attr : myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.covariate)) {
            if (attr.missing().cardinality() > 0) {
                String msg = "There is missing data in the covariate " + attr.name();
                throw new IllegalArgumentException(msg);
            }
        }
        for (PhenotypeAttribute attr : myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data)) {
            if (attr.missing().cardinality() > 0) {
                String msg = "There is missing data in the phenotype " + attr.name();
                throw new IllegalArgumentException(msg);
            }
        }
    }

    private double[] referenceProbabilitiesForSite(int site) {
        double[] probs =
                AssociationUtils.convertFloatArrayToDouble(myGenoPheno.referenceProb(site));
        return replaceNansWithMean(probs);
    }

    private double[] additiveSite(int site) {
        //codes genotypes as homozygous major = 1, homozygous minor = -1, het = 0
        //set missing values to the mean
        int nobs = myGenoPheno.phenotype().numberOfObservations();
        byte major = myGenotype.majorAllele(site);
        byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;

        double[] code = IntStream.range(0, nobs).mapToDouble(t -> {
            byte geno = myGenoPheno.genotype(t, site);
            if (geno == NN)
                return Double.NaN;
            byte[] alleles = GenotypeTableUtils.getDiploidValues(geno);
            double val = 0;
            if (alleles[0] == major)
                val += 0.5;
            if (alleles[1] == major)
                val += 0.5;
            return val;
        }).toArray();

        return replaceNansWithMean(code);
    }

    private double[] replaceNansWithMean(double[] array) {
        //replaces the Nans in array with the mean and returns array as a convenience
        int n = array.length;
        double sum = 0;
        double count = 0;
        for (int i = 0; i < n; i++)
            if (!Double.isNaN(array[i])) {
                sum += array[i];
                count++;
            }
        double mean = sum / count;
        for (int i = 0; i < n; i++)
            if (Double.isNaN(array[i]))
                array[i] = mean;
        return array;
    }

    private List<double[]> additiveDominanceSite(int site) {
        //the first double[] is the additive site
        //the second double[] equals 1 for heterozygous sites, is 0 otherwise
        int ntaxa = myGenoPheno.numberOfObservations();
        List<double[]> result = new ArrayList<>();
        result.add(additiveSite(site));
        double[] dom = new double[ntaxa];
        for (int t = 0; t < ntaxa; t++)
            if (myGenoPheno.isHeterozygous(t, site))
                dom[t] = 1;
        result.add(dom);
        return result;
    }

    private void calculateR2Fromp() {
        //returns the value of R^2 corresponding to the value of F, f for which P(F>f) = alpha
        minR2 = new double[2];
        double p = 1 - maxp.value();
        int basedf = orthogonalSolver.baseDf();
        try {
            double F = Fdist[0].inverseCumulativeProbability(p);
            minR2[0] = F / (numberOfObservations - 1 - basedf + F);
        } catch (OutOfRangeException e) {
            e.printStackTrace();
            minR2[0] = Double.NaN;
        }

        try {
            double F = Fdist[1].inverseCumulativeProbability(p);
            minR2[1] = 2 * F / (numberOfObservations - 2 - basedf + 2 * F);
        } catch (OutOfRangeException e) {
            e.printStackTrace();
            minR2[1] = Double.NaN;
        }
    }

    //abstract plugin methods that need to be overridden
    @Override
    public ImageIcon getIcon() {
        URL imageURL =
                EqtlAssociationPlugin.class.getResource("/net/maizegenetics/analysis/images/speed.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Fast Association";
    }

    @Override
    public String getToolTipText() {
        return "Use a fixed effect linear model to test variants quickly.";
    }

    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(EqtlAssociationPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    public TableReport runPlugin(DataSet input) {
        return (TableReport) performFunction(input).getData(0).getData();
    }

    /**
     * The maximum p-value that will be output by the analysis.
     *
     * @return MaxPValue
     */
    public Double maxPValue() {
        return maxp.value();
    }

    /**
     * Set MaxPValue. The maximum p-value that will be output
     * by the analysis.
     *
     * @param value MaxPValue
     *
     * @return this plugin
     */
    public EqtlAssociationPlugin maxPValue(Double value) {
        maxp = new PluginParameter<>(maxp, value);
        return this;
    }

    /**
     * Should an additive only model be fit? If true, an additive
     * model will be fit. If false, an additive + dominance
     * model will be fit. Default = false.
     *
     * @return Additive Only Model
     */
    public Boolean additiveOnlyModel() {
        return addOnly.value();
    }

    /**
     * Set Additive Only Model. Should an additive only model
     * be fit? If true, an additive model will be fit. If
     * false, an additive + dominance model will be fit. Default
     * = false.
     *
     * @param value Additive Only Model
     *
     * @return this plugin
     */
    public EqtlAssociationPlugin additiveOnlyModel(Boolean value) {
        addOnly = new PluginParameter<>(addOnly, value);
        return this;
    }

    /**
     * If the genotype table contains more than one type of
     * genotype data, choose the type to use for the analysis.
     *
     * @return Genotype Component
     */
    public GENOTYPE_TABLE_COMPONENT genotypeComponent() {
        return myGenotypeTable.value();
    }

    /**
     * Set Genotype Component. If the genotype table contains
     * more than one type of genotype data, choose the type
     * to use for the analysis.
     *
     * @param value Genotype Component
     *
     * @return this plugin
     */
    public EqtlAssociationPlugin genotypeComponent(GENOTYPE_TABLE_COMPONENT value) {
        myGenotypeTable = new PluginParameter<>(myGenotypeTable, value);
        return this;
    }

    /**
     * Should the results be saved to a file rather than stored
     * in memory? It true, the results will be written to
     * a file as each SNP is analyzed in order to reduce memory
     * requirementsand the results will NOT be saved to the
     * data tree. Default = false.
     *
     * @return Write to file
     */
    public Boolean writeToFile() {
        return saveAsFile.value();
    }

    /**
     * Set Write to file. Should the results be saved to a
     * file rather than stored in memory? It true, the results
     * will be written to a file as each SNP is analyzed in
     * order to reduce memory requirementsand the results
     * will NOT be saved to the data tree. Default = false.
     *
     * @param value Write to file
     *
     * @return this plugin
     */
    public EqtlAssociationPlugin writeToFile(Boolean value) {
        saveAsFile = new PluginParameter<>(saveAsFile, value);
        return this;
    }

    /**
     * The name of the file to which these results will be
     * saved.
     *
     * @return Output File
     */
    public String outputFile() {
        return reportFilename.value();
    }

    /**
     * Set Output File. The name of the file to which these
     * results will be saved.
     *
     * @param value Output File
     *
     * @return this plugin
     */
    public EqtlAssociationPlugin outputFile(String value) {
        reportFilename = new PluginParameter<>(reportFilename, value);
        return this;
    }

    @Override
    public String getCitation() {
        String citation =
                "Shabalin, AA. (2012) Matrix eQTL: ultra fast eQTL analysis via large matrix operations. "
                        + "Bioinformatics 28:1353-1358";
        return citation;
    }
}
