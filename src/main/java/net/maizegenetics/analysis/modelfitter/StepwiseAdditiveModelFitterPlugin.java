package net.maizegenetics.analysis.modelfitter;

import java.awt.Frame;
import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import javax.swing.ImageIcon;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


import net.maizegenetics.analysis.modelfitter.AdditiveSite.CRITERION;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable.GENOTYPE_TABLE_COMPONENT;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportUtils;

public class StepwiseAdditiveModelFitterPlugin extends AbstractPlugin {
    //needs to set the following parameters in StepwiseAdditiveModelFitter:
    //    private int numberOfPermutations = 1000;
    //    private double permutationAlpha = 0.05;
    //    private double enterLimit = 1e-5;
    //    private double exitLimit = 2e-5;
    //    private boolean useReferenceProbability = true;
    //    private boolean isNested = true;
    //    private String nestingEffectName = "family";
    //    private AdditiveSite.CRITERION modelSelectionCriterion = AdditiveSite.CRITERION.pval;
    //    private int maxSitesInModel = 1000;
    private static Logger myLogger = LogManager.getLogger(StepwiseAdditiveModelFitterPlugin.class);
    private List<String> myFactorNameList;
    private GenotypePhenotype myGenoPheno;
    private String datasetName;

    private PluginParameter<AdditiveSite.CRITERION> modelCriterion =
            new PluginParameter.Builder<>("criterion", AdditiveSite.CRITERION.pval, AdditiveSite.CRITERION.class)
                    .range(AdditiveSite.CRITERION.values())
                    .guiName("Model selection criterion")
                    .description("The model selection criterion used to determine which terms enter the model and how many. Value must be one of pval, bic, mbic, or aic")
                    .build();

    private PluginParameter<Boolean> useResiduals =
            new PluginParameter.Builder<>("useResidual", false, Boolean.class)
                    .description("Should the each new term be tested using residuals from the previous model instead of fitting a complete model each time?")
                    .guiName("Fit using residuals")
                    .build();

    private PluginParameter<Boolean> usePermutations =
            new PluginParameter.Builder<>("usePerm", true, Boolean.class)
                    .description("Should permutations be used to set the enter and exit limits for stepwise regression? A permutation test will be used to determine the enter limit. The exit limit will be set to 2 times the enter limit.")
                    .guiName("Use permutations")
                    .build();

    private PluginParameter<Integer> numberOfPermutations =
            new PluginParameter.Builder<>("nPerm", 1000, Integer.class)
                    .description("The number of permutations used to determine the enter limit.")
                    .guiName("Number of permutations")
                    .dependentOnParameter(usePermutations)
                    .build();

    private PluginParameter<Double> permutationAlpha =
            new PluginParameter.Builder<>("permAlpha", 0.05, Double.class)
                    .description("Type I errors will be controlled at this level.")
                    .guiName("Alpha for permutations")
                    .dependentOnParameter(usePermutations)
                    .build();

    private PluginParameter<Double> enterLimit =
            new PluginParameter.Builder<>("enterLimit", 1e-5, Double.class)
                    .description("When p-value is the model selection criteria, model fitting will stop when the next term chosen has a p-value greater than the enterLimit. This value will be over-ridden the permutation test, if used.")
                    .guiName("enterLimit")
                    .dependentOnParameter(usePermutations, false)
                    .build();

    private PluginParameter<Double> exitLimit =
            new PluginParameter.Builder<>("exitLimit", 2e-5, Double.class)
                    .description("During the backward step of model fitting if p-value has been chosen as the selection criterion, if the term in model with the highest p-value has a p-value > exitLimit, it will be removed from the model.")
                    .guiName("exitLimit")
                    .dependentOnParameter(usePermutations, false)
                    .build();

    private PluginParameter<Integer> maxTerms =
            new PluginParameter.Builder<>("maxterms", 1000, Integer.class)
                    .description("The maximum number of SNPs/markers that will be fit. If the model selection criterion is met first, then the fitting process will stop at that point.")
                    .guiName("Max SNPs/markers")
                    .build();

    private PluginParameter<Boolean> isNested =
            new PluginParameter.Builder<>("isNested", true, Boolean.class)
                    .description("Should SNPs/markers be nested within a factor, such as family?")
                    .guiName("")
                    .build();

    private PluginParameter<String> nestingFactor =
            new PluginParameter.Builder<>("nestFactor", null, String.class)
                    .guiName("Nesting factor")
                    .description("Nest markers within this factor. This parameter cannot be set from the command line. Instead, the first factor in the data set will be used.")
                    .dependentOnParameter(isNested)
                    .objectListSingleSelect()
                    .build();

    private PluginParameter<GenotypeTable.GENOTYPE_TABLE_COMPONENT> myGenotypeTable =
            new PluginParameter.Builder<>("genotypeComponent", GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT.class)
                    .genotypeTable()
                    .range(GenotypeTable.GENOTYPE_TABLE_COMPONENT.values())
                    .description("If the genotype table contains more than one type of genotype data, choose the type to use for the analysis.")
                    .build();

    private PluginParameter<Boolean> createAnova =
            new PluginParameter.Builder<>("anova", true, Boolean.class)
                    .description("Create pre- and post-scan anova reports.")
                    .guiName("Create anova reports")
                    .build();

    private PluginParameter<Boolean> createEffects =
            new PluginParameter.Builder<>("effects", true, Boolean.class)
                    .description("Create a report of marker effects based on the scan results.")
                    .guiName("Create effects report")
                    .build();

    private PluginParameter<Boolean> createEffectsPrescan =
            new PluginParameter.Builder<>("effectsPrescan", false, Boolean.class)
                    .description("Create a report of marker effects based on the results pre-scan.")
                    .guiName("Create prescan effects")
                    .build();

    private PluginParameter<Boolean> createStep =
            new PluginParameter.Builder<>("step", true, Boolean.class)
                    .description("Create a report of the which markers enter and leave the model as it is being fit.")
                    .guiName("Create step report")
                    .build();

    private PluginParameter<Boolean> createResiduals =
            new PluginParameter.Builder<>("residuals", false, Boolean.class)
                    .description("Create a phenotype dataset of model residuals for each chromosome. For each chromosome, the residuals will be calculated from a model with all terms EXCEPT the markers on that chromosome.")
                    .guiName("Create residuals")
                    .build();

    private PluginParameter<Boolean> writeFiles =
            new PluginParameter.Builder<>("saveToFile", false, Boolean.class)
                    .description("Should the requested output be written to files?")
                    .guiName("Write to files")
                    .build();

    private PluginParameter<String> outputName =
            new PluginParameter.Builder<>("savePath", "", String.class)
                    .description("The base file path for the save files. Each file saved will add a descriptive name to the base name.")
                    .guiName("Base file path")
                    .outFile()
                    .dependentOnParameter(writeFiles)
                    .build();

    private PluginParameter<Boolean> fitAddDom =
            new PluginParameter.Builder<>("addDom", false, Boolean.class)
                    .description("Should the plugin fit an additive plus dominance model? If false, an additive only model will be fit. Note this option does not implement a nested model. The isNested parameter will be ignored.")
                    .guiName("Fit AddDom Model")
                    .build();

    private PluginParameter<Integer> minHets =
            new PluginParameter.Builder<>("minHets", 20, Integer.class)
                    .dependentOnParameter(fitAddDom)
                    .description("The minimum number of individuals that are heterozygous at a site for the dominance term to be included. If the number of hets is less than minHets at a site, an additive only model will be fit for that site.")
                    .guiName("Minimum number of heterozygotes")
                    .build();

    private PluginParameter<Boolean> imputedDominance =
            new PluginParameter.Builder<>("imputeDom", true, Boolean.class)
                    .description("Should the plugin impute a value for dominance if data is missing. If false, the dominance score in the design matrix will be 0 for missing genotypes.")
                    .guiName("Impute Dominance Score")
                    .build();

    public StepwiseAdditiveModelFitterPlugin() {
        super(null, false);
    }

    public StepwiseAdditiveModelFitterPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public String pluginDescription() {
        return "A stepwise model fitter designed to select variants tested for association with a large " +
                "number of nucleotide variants. It can be multi-threaded to handle large number of variants efficiently. " +
                "It can fit an additive only model or an additive + dominance model.";
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        DoubleMatrixFactory.setDefault(DoubleMatrixFactory.FactoryType.ejml);
        //input data should be a single GenotypePhenotype
        List<Datum> datumList = input.getDataOfType(GenotypePhenotype.class);
        if (datumList.size() != 1)
            throw new IllegalArgumentException("Choose exactly one dataset that has combined genotype and phenotype data.");
        myGenoPheno = (GenotypePhenotype) datumList.get(0).getData();
        datasetName = datumList.get(0).getName();

        myFactorNameList =
                myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.factor).stream()
                        .map(pa -> pa.name())
                        .collect(Collectors.toList());

        if (myFactorNameList.isEmpty()) myFactorNameList.add("None");
        nestingFactor = new PluginParameter<>(nestingFactor, myFactorNameList);
    }

    @Override
    protected void postProcessParameters() {
        if (isNested.value()) {
            if (myFactorNameList.size() > 1 && nestingFactor.value().isEmpty()) throw new IllegalArgumentException("Is Nested was checked but no nesting factor was selected."); 
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        StepwiseAdditiveModelFitter stamFitter;
        if (fitAddDom.value()) {
            stamFitter = new StepwiseAddDomModelFitter(myGenoPheno, datasetName);
            //initialize constants
            AddPlusDomModelEffect.MIN_HETS = minHets.value();
            AddPlusDomModelEffect.IMPUTE_DOM = imputedDominance.value();

        } else {
            stamFitter = new StepwiseAdditiveModelFitter(myGenoPheno, datasetName);
        }

        //set parameters
        if (usePermutations.value()) {
            stamFitter.numberOfPermutations(numberOfPermutations.value().intValue());
            stamFitter.permutationAlpha(permutationAlpha.value());
        } else {
            stamFitter.enterLimit(enterLimit.value());
            stamFitter.exitLimit(exitLimit.value());
        }

        if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability)
            stamFitter.useReferenceProbability(true);
        else if (myGenoPheno.genotypeTable().hasGenotype() == false
                && myGenoPheno.genotypeTable().hasReferenceProbablity() == true)
            stamFitter.useReferenceProbability(true);
        else
            stamFitter.useReferenceProbability(false);

        if (isNested.value()) {
            if (isInteractive()) {
                String nestingList = nestingFactor.value();
                if (nestingList.isEmpty()) {
                    if (myFactorNameList.get(0).equals("None")) {
                        stamFitter.isNested(false);
                    } else {
                        stamFitter.isNested(true);
                        stamFitter.nestingEffectName(nestingList);
                    }
                } else {
                    String factorName = nestingList;
                    if (factorName.equals("None"))
                        stamFitter.isNested(false);
                    else {
                        stamFitter.isNested(true);
                        stamFitter.nestingEffectName(factorName);
                    }
                }
            } else {
                List<PhenotypeAttribute> factorList = myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.factor);
                if (factorList.isEmpty()) stamFitter.isNested(false);
                else {
                    stamFitter.isNested(true);
                    stamFitter.nestingEffectName(factorList.get(0).name());
                }
            }

        } else {
            stamFitter.isNested(false);
        }

        stamFitter.enterLimit(enterLimit.value());
        stamFitter.exitLimit(exitLimit.value());
        stamFitter.maxSitesInModel(maxTerms.value());
        stamFitter.modelSelectionCriterion(modelCriterion.value());
        stamFitter.createAnovaReport(createAnova.value());
        stamFitter.createPostScanEffectsReport(createEffects.value());
        stamFitter.createPreScanEffectsReport(createEffectsPrescan.value());
        stamFitter.createResidualsByChr(createResiduals.value());
        stamFitter.createStepReport(createStep.value());

        long start = System.nanoTime();
        stamFitter.runAnalysis();
        myLogger.debug(String.format("ran analysis in %d ms.",  (System.nanoTime() - start)/1000000));
        
        List<Datum> resultData = new ArrayList<>();
        if (createStep.value()) {
            String name = "Steps_" + datasetName;
            String comment = "Stepwise regression results:\nModel fitting steps\n";
            TableReport output = stamFitter.getSteps();
            resultData.add(new Datum(name, output, comment));

            if (writeFiles.value()) {
                String filename = outputName.value() + "_steps.txt";
                TableReportUtils.saveDelimitedTableReport(output, new File(filename));
            }
        }
        if (createAnova.value()) {
            TableReport output1 = stamFitter.getAnovaReport();
            String name = "Prescan_anova_" + datasetName;
            String comment = "Stepwise regression results:\nAnova for the initial model prior to rescanning\n";
            resultData.add(new Datum(name, output1, comment));

            TableReport output2 = stamFitter.getAnovaReportWithCI();
            name = "Anova_" + datasetName;
            comment = "Stepwise regression results:\nAnova for the final model with support intervals\n";
            resultData.add(new Datum(name, output2, comment));

            if (writeFiles.value()) {
                String filename = outputName.value() + "_prescan_anova.txt";
                TableReportUtils.saveDelimitedTableReport(output1, new File(filename));
                filename = outputName.value() + "_anova.txt";
                TableReportUtils.saveDelimitedTableReport(output2, new File(filename));
            }
        }
        if (createEffectsPrescan.value()) {
            String name = "Prescan_effects_" + datasetName;
            String comment = "Stepwise regression results:\nMarker effects estimated from the initial model before rescanning\n";
            TableReport output = stamFitter.getMarkerEffectReport();
            resultData.add(new Datum(name, output, comment));

            if (writeFiles.value()) {
                String filename = outputName.value() + "_prescan_effects.txt";
                TableReportUtils.saveDelimitedTableReport(output, new File(filename));
            }
        }
        if (createEffects.value()) {
            String name = "Effects_regression_" + datasetName;
            String comment = "Stepwise regression results:\nMarker effects estimated from the final model\n";
            TableReport output = stamFitter.getMarkerEffectReportWithCI();
            resultData.add(new Datum(name, output, comment));

            if (writeFiles.value()) {
                String filename = outputName.value() + "_effects.txt";
                TableReportUtils.saveDelimitedTableReport(output, new File(filename));
            }
        }
        if (createResiduals.value()) {
            List<Phenotype> phenoList = stamFitter.getResidualPhenotypesByChromosome();
            phenoList.stream().forEach(pheno -> {
                String traitname = pheno.attributeName(pheno.numberOfAttributes() - 1);
                String name = "Residuals_" + traitname + "_" + datasetName;
                String comment = "Stepwise regression results:\nResiduals for " + traitname + "\n";
                resultData.add(new Datum(name, pheno, comment));
            });

            if (writeFiles.value()) {
                phenoList.stream().forEach(pheno -> {
                    String traitname = pheno.attributeName(pheno.numberOfAttributes() - 1);
                    String filename = outputName.value() + "_" + traitname + ".txt";
                    PhenotypeUtils.write(pheno, filename);
                });
            }
        }

        return new DataSet(resultData, this);
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = StepwiseOLSModelFitterPlugin.class.getResource("stepwise.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Stepwise-Multithread";
    }

    @Override
    public String getToolTipText() {
        return "Fit a model using stepwise forward-backward regression.";
    }

//      The following getters and setters were auto-generated.
//      Please use this method to re-generate.

//      public static void main(String[] args) {
//          GeneratePluginCode.generate(StepwiseAdditiveModelFitterPlugin.class);
//      }

     /**
      * The model selection criterion used to determine which
      * terms enter the model and how many. Value must be one
      * of pval, bic, mbic, or aic
      *
      * @return Model selection criterion
      */
     public CRITERION modelCriterion() {
         return modelCriterion.value();
     }

     /**
      * Set Model selection criterion. The model selection
      * criterion used to determine which terms enter the model
      * and how many. Value must be one of pval, bic, mbic,
      * or aic
      *
      * @param value Model selection criterion
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin modelCriterion(CRITERION value) {
         modelCriterion = new PluginParameter<>(modelCriterion, value);
         return this;
     }

     /**
      * Should the each new term be tested using residuals
      * from the previous model instead of fitting a complete
      * model each time?
      *
      * @return Fit using residuals
      */
     public Boolean useResiduals() {
         return useResiduals.value();
     }

     /**
      * Set Fit using residuals. Should the each new term be
      * tested using residuals from the previous model instead
      * of fitting a complete model each time?
      *
      * @param value Fit using residuals
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin useResiduals(Boolean value) {
         useResiduals = new PluginParameter<>(useResiduals, value);
         return this;
     }

     /**
      * Should permutations be used to set the enter and exit
      * limits for stepwise regression? A permutation test
      * will be used to determine the enter limit. The exit
      * limit will be set to 2 times the enter limit.
      *
      * @return Use permutations
      */
     public Boolean usePermutations() {
         return usePermutations.value();
     }

     /**
      * Set Use permutations. Should permutations be used to
      * set the enter and exit limits for stepwise regression?
      * A permutation test will be used to determine the enter
      * limit. The exit limit will be set to 2 times the enter
      * limit.
      *
      * @param value Use permutations
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin usePermutations(Boolean value) {
         usePermutations = new PluginParameter<>(usePermutations, value);
         return this;
     }

     /**
      * The number of permutations used to determine the enter
      * limit.
      *
      * @return Number of permutations
      */
     public Integer numberOfPermutations() {
         return numberOfPermutations.value();
     }

     /**
      * Set Number of permutations. The number of permutations
      * used to determine the enter limit.
      *
      * @param value Number of permutations
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin numberOfPermutations(Integer value) {
         numberOfPermutations = new PluginParameter<>(numberOfPermutations, value);
         return this;
     }

     /**
      * Type I errors will be controlled at this level.
      *
      * @return Alpha for permutations
      */
     public Double permutationAlpha() {
         return permutationAlpha.value();
     }

     /**
      * Set Alpha for permutations. Type I errors will be controlled
      * at this level.
      *
      * @param value Alpha for permutations
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin permutationAlpha(Double value) {
         permutationAlpha = new PluginParameter<>(permutationAlpha, value);
         return this;
     }

     /**
      * When p-value is the model selection criteria, model
      * fitting will stop when the next term chosen has a p-value
      * greater than the enterLimit. This value will be over-ridden
      * the permutation test, if used.
      *
      * @return enterLimit
      */
     public Double enterLimit() {
         return enterLimit.value();
     }

     /**
      * Set enterLimit. When p-value is the model selection
      * criteria, model fitting will stop when the next term
      * chosen has a p-value greater than the enterLimit. This
      * value will be over-ridden the permutation test, if
      * used.
      *
      * @param value enterLimit
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin enterLimit(Double value) {
         enterLimit = new PluginParameter<>(enterLimit, value);
         return this;
     }

     /**
      * During the backward step of model fitting if p-value
      * has been chosen as the selection criterion, if the
      * term in model with the highest p-value has a p-value
      * > exitLimit, it will be removed from the model.
      *
      * @return exitLimit
      */
     public Double exitLimit() {
         return exitLimit.value();
     }

     /**
      * Set exitLimit. During the backward step of model fitting
      * if p-value has been chosen as the selection criterion,
      * if the term in model with the highest p-value has a
      * p-value > exitLimit, it will be removed from the model.
      *
      * @param value exitLimit
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin exitLimit(Double value) {
         exitLimit = new PluginParameter<>(exitLimit, value);
         return this;
     }

     /**
      * The maximum number of SNPs/markers that will be fit.
      * If the model selection criterion is met first, then
      * the fitting process will stop at that point.
      *
      * @return Max SNPs/markers
      */
     public Integer maxTerms() {
         return maxTerms.value();
     }

     /**
      * Set Max SNPs/markers. The maximum number of SNPs/markers
      * that will be fit. If the model selection criterion
      * is met first, then the fitting process will stop at
      * that point.
      *
      * @param value Max SNPs/markers
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin maxTerms(Integer value) {
         maxTerms = new PluginParameter<>(maxTerms, value);
         return this;
     }

     /**
      * Should SNPs/markers be nested within a factor, such
      * as family?
      *
      * @return Is Nested
      */
     public Boolean isNested() {
         return isNested.value();
     }

     /**
      * Set Is Nested. Should SNPs/markers be nested within
      * a factor, such as family?
      *
      * @param value Is Nested
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin isNested(Boolean value) {
         isNested = new PluginParameter<>(isNested, value);
         return this;
     }

     /**
      * Nest markers within this factor. This parameter cannot
      * be set from the command line. Instead, the first factor
      * in the data set will be used.
      *
      * @return Nesting factor
      */
     public String nestingFactor() {
         return nestingFactor.value();
     }

     /**
      * Set Nesting factor. Nest markers within this factor.
      * This parameter cannot be set from the command line.
      * Instead, the first factor in the data set will be used.
      *
      * @param value Nesting factor
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin nestingFactor(List<String> value) {
         nestingFactor = new PluginParameter<>(nestingFactor, value);
         return this;
     }

     /**
      * If the genotype table contains more than one type of
      * genotype data, choose the type to use for the analysis.
      *
      * @return Genotype Component
      */
     public GENOTYPE_TABLE_COMPONENT genotypeTable() {
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
     public StepwiseAdditiveModelFitterPlugin genotypeTable(GENOTYPE_TABLE_COMPONENT value) {
         myGenotypeTable = new PluginParameter<>(myGenotypeTable, value);
         return this;
     }

     /**
      * Create pre- and post-scan anova reports.
      *
      * @return Create anova reports
      */
     public Boolean createAnova() {
         return createAnova.value();
     }

     /**
      * Set Create anova reports. Create pre- and post-scan
      * anova reports.
      *
      * @param value Create anova reports
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin createAnova(Boolean value) {
         createAnova = new PluginParameter<>(createAnova, value);
         return this;
     }

     /**
      * Create a report of marker effects based on the scan
      * results.
      *
      * @return Create effects report
      */
     public Boolean createEffects() {
         return createEffects.value();
     }

     /**
      * Set Create effects report. Create a report of marker
      * effects based on the scan results.
      *
      * @param value Create effects report
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin createEffects(Boolean value) {
         createEffects = new PluginParameter<>(createEffects, value);
         return this;
     }

     /**
      * Create a report of marker effects based on the results
      * pre-scan.
      *
      * @return Create prescan effects
      */
     public Boolean createEffectsPrescan() {
         return createEffectsPrescan.value();
     }

     /**
      * Set Create prescan effects. Create a report of marker
      * effects based on the results pre-scan.
      *
      * @param value Create prescan effects
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin createEffectsPrescan(Boolean value) {
         createEffectsPrescan = new PluginParameter<>(createEffectsPrescan, value);
         return this;
     }

     /**
      * Create a report of the which markers enter and leave
      * the model as it is being fit.
      *
      * @return Create step report
      */
     public Boolean createStep() {
         return createStep.value();
     }

     /**
      * Set Create step report. Create a report of the which
      * markers enter and leave the model as it is being fit.
      *
      * @param value Create step report
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin createStep(Boolean value) {
         createStep = new PluginParameter<>(createStep, value);
         return this;
     }

     /**
      * Create a phenotype dataset of model residuals for each
      * chromosome. For each chromosome, the residuals will
      * be calculated from a model with all terms EXCEPT the
      * markers on that chromosome.
      *
      * @return Create residuals
      */
     public Boolean createResiduals() {
         return createResiduals.value();
     }

     /**
      * Set Create residuals. Create a phenotype dataset of
      * model residuals for each chromosome. For each chromosome,
      * the residuals will be calculated from a model with
      * all terms EXCEPT the markers on that chromosome.
      *
      * @param value Create residuals
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin createResiduals(Boolean value) {
         createResiduals = new PluginParameter<>(createResiduals, value);
         return this;
     }

     /**
      * Should the requested output be written to files?
      *
      * @return Write to files
      */
     public Boolean writeFiles() {
         return writeFiles.value();
     }

     /**
      * Set Write to files. Should the requested output be
      * written to files?
      *
      * @param value Write to files
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin writeFiles(Boolean value) {
         writeFiles = new PluginParameter<>(writeFiles, value);
         return this;
     }

     /**
      * The base file path for the save files. Each file saved
      * will add a descriptive name to the base name.
      *
      * @return Base file path
      */
     public String outputName() {
         return outputName.value();
     }

     /**
      * Set Base file path. The base file path for the save
      * files. Each file saved will add a descriptive name
      * to the base name.
      *
      * @param value Base file path
      *
      * @return this plugin
      */
     public StepwiseAdditiveModelFitterPlugin outputName(String value) {
         outputName = new PluginParameter<>(outputName, value);
         return this;
     }

    /**
     * Should the plugin fit an additive plus dominance model?
     * If false, an additive only model will be fit. Note
     * this option does not implement a nested model. The
     * isNested parameter will be ignored.
     *
     * @return Fit AddDom Model
     */
    public Boolean fitAddDom() {
        return fitAddDom.value();
    }

    /**
     * Set Fit AddDom Model. Should the plugin fit an additive
     * plus dominance model? If false, an additive only model
     * will be fit. Note this option does not implement a
     * nested model. The isNested parameter will be ignored.
     *
     * @param value Fit AddDom Model
     *
     * @return this plugin
     */
    public StepwiseAdditiveModelFitterPlugin fitAddDom(Boolean value) {
        fitAddDom = new PluginParameter<>(fitAddDom, value);
        return this;
    }

    /**
     * The minimum number of individuals that are heterozygous
     * at a site for the dominance term to be included. If
     * the number of hets is less than minHets at a site,
     * an additive only model will be fit for that site.
     *
     * @return Minimum number of heterozygotes
     */
    public Integer minHets() {
        return minHets.value();
    }

    /**
     * Set Minimum number of heterozygotes. The minimum number
     * of individuals that are heterozygous at a site for
     * the dominance term to be included. If the number of
     * hets is less than minHets at a site, an additive only
     * model will be fit for that site.
     *
     * @param value Minimum number of heterozygotes
     *
     * @return this plugin
     */
    public StepwiseAdditiveModelFitterPlugin minHets(Integer value) {
        minHets = new PluginParameter<>(minHets, value);
        return this;
    }

    /**
     * Should the plugin impute a value for dominance if data
     * is missing. If false, the dominance score in the design
     * matrix will be 0 for missing genotypes.
     *
     * @return Impute Dominance Score
     */
    public Boolean imputedDominance() {
        return imputedDominance.value();
    }

    /**
     * Set Impute Dominance Score. Should the plugin impute
     * a value for dominance if data is missing. If false,
     * the dominance score in the design matrix will be 0
     * for missing genotypes.
     *
     * @param value Impute Dominance Score
     *
     * @return this plugin
     */
    public StepwiseAdditiveModelFitterPlugin imputedDominance(Boolean value) {
        imputedDominance = new PluginParameter<>(imputedDominance, value);
        return this;
    }

}
