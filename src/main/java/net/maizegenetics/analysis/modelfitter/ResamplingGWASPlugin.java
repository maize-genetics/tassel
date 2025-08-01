package net.maizegenetics.analysis.modelfitter;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javax.swing.ImageIcon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.analysis.association.FixedEffectLMPlugin;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.stats.linearmodels.BasicShuffler;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportBuilder;

public class ResamplingGWASPlugin extends AbstractPlugin {
    //input parameters :
    //enterlimit = the p-value cutoff for a term to enter the model
    //select next term using residuals
    //number of iterations
    //resample method (with/without replacement)
    //resample proportion
    //maximum number of terms in the model

    //restrictions:
    //there must be no missing data in the Phenotype
    //the input data set must contain one GenotypeTable and one Phenotype
    //at most one factor is allowed (expected to be family or pop)

    private Random randomGen = new Random();
    int numberOfFactors;

    private static Logger myLogger = LogManager.getLogger(ResamplingGWASPlugin.class);

    private PluginParameter<Double> enterLimit =
            new PluginParameter.Builder<>("enterLimit", 1e-8, Double.class)
                    .description("A new term entering the model must have a p-value equal to or less than the enter limit.")
                    .guiName("Enter Limit")
                    .build();
    private PluginParameter<Integer> maxModelTerms =
            new PluginParameter.Builder<>("maxterms", 100, Integer.class)
                    .description("The maximum number of variants that will be fit. If the chromosome residuals are being fit, the maximum number of variants fit per chromosome.")
                    .guiName("Max terms")
                    .build();

    private PluginParameter<Boolean> useResiduals =
            new PluginParameter.Builder<>("residuals", false, Boolean.class)
                    .description("Should new terms be tested using residuals from the prior model? The analysis runs faster using this option.")
                    .guiName("Use residuals")
                    .build();

    private PluginParameter<Integer> numberOfIterations =
            new PluginParameter.Builder<>("nIterations", 100, Integer.class)
                    .description("The number of times the data should be resampled.")
                    .guiName("Number of Iterations")
                    .build();

    private PluginParameter<Double> resampleProportion =
            new PluginParameter.Builder<>("resampleProportion", 0.8, Double.class)
                    .description("The size of the sample is resample proportion times the number of observations in the complete data. For bootstrap, set this value to 1 and with replacement to true.")
                    .guiName("")
                    .build();

    private PluginParameter<Boolean> withReplacement =
            new PluginParameter.Builder<>("replacement", false, Boolean.class)
                    .description("Should the sample be formed by sampling with replacement?  For bootstrap, set resample proportion to 1 and this value to true.")
                    .guiName("With Replacement")
                    .build();

    private PluginParameter<Boolean> useSerialFile =
            new PluginParameter.Builder<>("useSitefile", false, Boolean.class)
                    .description("Use an additive site file as the source of genotypes.")
                    .guiName("Use Site File")
                    .build();

    private PluginParameter<String> serialFilename =
            new PluginParameter.Builder<>("sitefile", null, String.class)
                    .description("The name of the file containing genotypes stored in site objects.")
                    .guiName("Site Filename")
                    .dependentOnParameter(useSerialFile)
                    .build();

    private PluginParameter<Integer> maxThreads =
            new PluginParameter.Builder<>("maxThreads", -1, Integer.class)
                    .description("The maximum number of threads to available to the plugin. If maxThreads = -1, then the value is set to the number of available processors - 2.")
                    .guiName("Maximum Thread Number")
                    .build();

    public ResamplingGWASPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        Phenotype pheno;
        //check that there is only one GenotypePhenotype
        if (useSerialFile.value()) {
            //check that there is only one Phenotype
            List<Datum> datumList = input.getDataOfType(Phenotype.class);
            if (datumList.size() != 1) {
                throw new IllegalArgumentException("Exactly one joined Phenotype dataset must be supplied as input to Resample GWAS.");
            }
            pheno = (Phenotype) datumList.get(0).getData();
        } else {
            //check that there is only one GenotypePhenotype
            List<Datum> datumList = input.getDataOfType(GenotypePhenotype.class);
            if (datumList.size() != 1) {
                throw new IllegalArgumentException("Exactly one joined Genotype-Phenotype dataset must be supplied as input to Resample GWAS.");
            }
            pheno = ((GenotypePhenotype) datumList.get(0).getData()).phenotype();
        }

        //if the Phenotype has more than one factor level, throw an error
        numberOfFactors = pheno.numberOfAttributesOfType(ATTRIBUTE_TYPE.factor);
        if (numberOfFactors > 1)
            throw new IllegalArgumentException("Phenotype supplied to Resample GWAS can have at most one factor");

        //if there is any missing data in the phenotype throw an error
        boolean anyMissing = pheno.attributeStream().anyMatch(pa -> pa.missing().cardinality() > 0);
        if (anyMissing)
            throw new IllegalArgumentException("No missing phenotype data allow as input to Resample GWAS");
    }

    @Override
    public DataSet processData(DataSet input) {
        long mainStart = System.nanoTime();
        Phenotype pheno;
        String dataname;
        List<Datum> datumList;
        if (useSerialFile.value()) {
            datumList = input.getDataOfType(Phenotype.class);
            pheno = (Phenotype) datumList.get(0).getData();
        } else {
            datumList = input.getDataOfType(GenotypePhenotype.class);
            pheno = ((GenotypePhenotype) datumList.get(0).getData()).phenotype();
        }
        dataname = datumList.get(0).getName();

        //create the model fitter
        ForwardRegression modelfitter;
        if (useResiduals.value()) {
            //not implemented
            modelfitter = null;
        } else {
            if (useSerialFile.value()) {
                long start = System.nanoTime();
                if (maxThreads.value() > 0) {
                    modelfitter = new AdditiveModelForwardRegression(serialFilename.value(), pheno, maxThreads.value());
                } else {
                    modelfitter = new AdditiveModelForwardRegression(serialFilename.value(), pheno);
                }
                myLogger.debug(String.format("Serialized sites loaded in %d ms", (System.nanoTime() - start) / 1000000));
                pheno = modelfitter.phenotype();
            } else {
                GenotypePhenotype myGenoPheno = (GenotypePhenotype) datumList.get(0).getData();
                if (maxThreads.value() > 0) {
                    modelfitter = new AdditiveModelForwardRegression(myGenoPheno, maxThreads.value());
                } else {
                    modelfitter = new AdditiveModelForwardRegression(myGenoPheno);
                }
            }
        }
        
        int numberOfTraits = pheno.numberOfAttributesOfType(ATTRIBUTE_TYPE.data);
        List<int[]> factorLevelList = new ArrayList<>();

        //create a factor level list for creating subsamples
        //if there are nf factors, the list will contain nf int arrays, each containing the observation numbers for one factor level
        if (numberOfFactors == 1) {
            CategoricalAttribute myFactor =
                    (CategoricalAttribute) pheno.attributeListOfType(ATTRIBUTE_TYPE.factor).get(0);
            int[] levels = myFactor.allIntValues();
            factorLevelList = new ArrayList<>();
            Map<Integer, List<Integer>> factorGroups =
                    IntStream.range(0, levels.length).boxed().collect(Collectors.groupingBy(i -> new Integer(levels[i])));
            for (int i = 0; i < factorGroups.size(); i++) {
                List<Integer> factorMembers = factorGroups.get(i);
                factorLevelList.add(factorMembers.stream().mapToInt(I -> I.intValue()).toArray());
            }
        } else {  //number of factors = 0
            factorLevelList.add(IntStream.range(0, pheno.numberOfObservations()).toArray());
        }

        //report builder column labels
        String[] columns =
                new String[] { "trait", "interation", "step", "SnpID", "Chr", "Pos", "p-value",
                        "-log10p" };
        TableReportBuilder reportBuilder =
                TableReportBuilder.getInstance("Resample terms_" + dataname, columns);

        //for each trait
        long start = System.nanoTime();
        List<PhenotypeAttribute> dataAttributes = pheno.attributeListOfType(ATTRIBUTE_TYPE.data);
        for (int ph = 0; ph < numberOfTraits; ph++) {
            modelfitter.resetModel(ph, enterLimit.value(), maxModelTerms.value());
            myLogger.debug(String.format("Analyzing phenotype %d, %s.", ph, dataAttributes.get(ph).name()));

            //for each iteration
            for (int iter = 0; iter < numberOfIterations.value(); iter++) {
                myLogger.debug(String.format("phenotype %d, iteration %d", ph, iter));

                //create a random subsample
                int[] subsample = randomSample(factorLevelList);

                //generate a model
                modelfitter.fitModelForSubsample(subsample, iter);

            }

            //add the model terms and p-values to the result TableReportBuilder
            for (Object[] row : modelfitter.fittedModel())
                reportBuilder.add(row);

        }
        myLogger.debug(String.format("Resample GWAS model fit in %d ms.", (System.nanoTime() - start) / 1000000));
        //return a DataSet that contains the TableReportBuilder
        myLogger.debug(String.format("Elapse time for plugin = %d ms.", (System.nanoTime() - mainStart) / 1000000));
        Datum theResult = new Datum("name", reportBuilder.build(), "comment");
        return new DataSet(theResult, this);
    }

    private int[] randomSample(List<int[]> factorLevelList) {
        //draw samples within factor levels
        double rp = resampleProportion.value().doubleValue();
        if (withReplacement.value()) {
            return factorLevelList.stream()
                    .flatMapToInt(iarray -> randomGen.ints((int) Math.round(iarray.length * rp), 0, iarray.length).map(i -> iarray[i]))
                    .toArray();
        } else {
            return factorLevelList.stream()
                    .flatMapToInt(iarray -> {
                        int[] copy = Arrays.copyOf(iarray, iarray.length);
                        BasicShuffler.shuffle(copy);
                        return Arrays.stream(copy).limit((int) Math.round(iarray.length * rp));
                    })
                    .toArray();
        }
    }

    @Override
    public String pluginDescription() {
        return "ResamplingGWASPlugin uses forward regression to fit a multiple SNP model to each of a number of samples drawn from the phenotype data. ";
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL =
                FixedEffectLMPlugin.class.getResource("/net/maizegenetics/analysis/images/resample.png");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Resample GWAS";
    }

    @Override
    public String getToolTipText() {
        return "GWAS with resampling";
    }

    /**
     * Initializes the random number generator with a seed. Primarily used for testing to generate reproducible results.
     * @param seed      the seed for the Random object
     */
    public void setRandomSeed(int seed) {
        randomGen = new Random(seed);
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(ResamplingGWASPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    public TableReport runPlugin(DataSet input) {
        return (TableReport) performFunction(input).getData(0).getData();
    }

    /**
     * A new term entering the model must have a p-value equal
     * to or less than the enter limit. (Default = 1e-8)
     *
     * @return Enter Limit
     */
    public Double enterLimit() {
        return enterLimit.value();
    }

    /**
     * Set Enter Limit. A new term entering the model must
     * have a p-value equal to or less than the enter limit.
     * (Default = 1e-8)
     *
     * @param value Enter Limit
     *
     * @return this plugin
     */
    public ResamplingGWASPlugin enterLimit(Double value) {
        enterLimit = new PluginParameter<>(enterLimit, value);
        return this;
    }

    /**
     * Should new terms be tested using residuals from the
     * prior model? The analysis runs faster using this option.
     * (Default = true)
     *
     * @return Use residuals
     */
    public Boolean useResiduals() {
        return useResiduals.value();
    }

    /**
     * Set Use residuals. Should new terms be tested using
     * residuals from the prior model? The analysis runs faster
     * using this option. (Default = true)
     *
     * @param value Use residuals
     *
     * @return this plugin
     */
    public ResamplingGWASPlugin useResiduals(Boolean value) {
        useResiduals = new PluginParameter<>(useResiduals, value);
        return this;
    }

    /**
     * The number of times the data should be resampled. (Default
     * = 100)
     *
     * @return Number of Iterations
     */
    public Integer numberOfIterations() {
        return numberOfIterations.value();
    }

    /**
     * Set Number of Iterations. The number of times the data
     * should be resampled. (Default = 100)
     *
     * @param value Number of Iterations
     *
     * @return this plugin
     */
    public ResamplingGWASPlugin numberOfIterations(Integer value) {
        numberOfIterations = new PluginParameter<>(numberOfIterations, value);
        return this;
    }

    /**
     * The size of the sample is resample proportion times
     * the number of observations in the complete data. For
     * bootstrap, set this value to 1 and with replacement
     * to true. (Default = 0.8)
     *
     * @return Resample Proportion
     */
    public Double resampleProportion() {
        return resampleProportion.value();
    }

    /**
     * Set Resample Proportion. The size of the sample is
     * resample proportion times the number of observations
     * in the complete data. For bootstrap, set this value
     * to 1 and with replacement to true. (Default = 0.8)
     *
     * @param value Resample Proportion
     *
     * @return this plugin
     */
    public ResamplingGWASPlugin resampleProportion(Double value) {
        resampleProportion = new PluginParameter<>(resampleProportion, value);
        return this;
    }

    /**
     * Should the sample be formed by sampling with replacement?
     * (Default = false)
     *
     * @return With Replacement
     */
    public Boolean withReplacement() {
        return withReplacement.value();
    }

    /**
     * Set With Replacement. Should the sample be formed by
     * sampling with replacement? (Default = false)
     *
     * @param value With Replacement
     *
     * @return this plugin
     */
    public ResamplingGWASPlugin withReplacement(Boolean value) {
        withReplacement = new PluginParameter<>(withReplacement, value);
        return this;
    }

}
