package net.maizegenetics.analysis.modelfitter;

import com.google.common.collect.Range;

import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.util.TableReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;

import java.awt.*;
import java.net.URL;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Stepwise Ordinary Least Squares model fitter
 *
 * @author Alex Lipka
 * @author Peter Bradbury
 */
public class StepwiseOLSModelFitterPlugin extends AbstractPlugin {

    private List<String> myFactorList;
    private GenotypePhenotype myGenoPheno;
    private String datasetName;
    private final String NONE = "None";

    private PluginParameter<StepwiseOLSModelFitter.MODEL_TYPE> modelType =
            new PluginParameter.Builder<>("modelType", StepwiseOLSModelFitter.MODEL_TYPE.pvalue, StepwiseOLSModelFitter.MODEL_TYPE.class)
                    .range(StepwiseOLSModelFitter.MODEL_TYPE.values())
                    .guiName("Model type")
                    .description("The model selection criteria used to determine which terms enter the model and how many. Value must be one of pvalue, bic, mbic, or aic")
                    .build();
    private PluginParameter<Double> enterlimit =
            new PluginParameter.Builder<>("enter", 0.01, Double.class)
                    .range(Range.closed(0.0, 1.0))
                    .guiName("Entry limit")
                    .description("The enter limit or maximum p-value for which a term can enter the model")
                    .build();
    private PluginParameter<Double> exitlimit =
            new PluginParameter.Builder<>("exit", 0.01, Double.class)
                    .range(Range.closed(0.0, 1.0))
                    .guiName("Exit limit")
                    .description("A term exits the model on a backward step if its p-value is greater than this value")
                    .build();
    private PluginParameter<Integer> maxNumberOfMarkers =
            new PluginParameter.Builder<>("maxMarkers", 20, Integer.class)
                    .range(Range.closed(0, 10000))
                    .guiName("Maximum markers")
                    .description("The maximum number of markers that will be fit, if the enter limit is not reached first")
                    .build();
    private PluginParameter<Boolean> nestMarkers =
            new PluginParameter.Builder<>("nestMarkers", false, Boolean.class)
                    .guiName("Nest markers")
                    .description("Should markers be nested within a model factor")
                    .build();
    private PluginParameter<String> nestingFactor =
            new PluginParameter.Builder<>("nestFactor", null, String.class)
                    .guiName("Nesting factor")
                    .description("Nest markers within this factor.")
                    .dependentOnParameter(nestMarkers)
                    .objectListSingleSelect()
                    .build();
    private PluginParameter<Integer> numberOfPermutations =
            new PluginParameter.Builder<>("nperm", 0, Integer.class)
                    .range(Range.closed(0, 100000))
                    .guiName("Number of permutations")
                    .description("Number of permutations for the model to determine an empirical alpha")
                    .build();
//    private PluginParameter<Boolean> chromosomeResiduals =
//            new PluginParameter.Builder<>("chrResidual", false, Boolean.class)
//                    .guiName("Output chromosome residuals")
//                    .description("Should a dataset of chromosome residuals be created for each phenotype? The output datasets will include all factors but no covariates from the original phenotype data.")
//                    .build();
//    private PluginParameter<Boolean> residualsAsFile =
//            new PluginParameter.Builder<>("resAsFile", false, Boolean.class)
//                    .guiName("Save residuals as file?")
//                    .description("Should the chromosome residuals to be saved to separate files rather than stored in memory?")
//                    .dependentOnParameter(chromosomeResiduals)
//                    .build();
//    private PluginParameter<String> residualFilebase =
//            new PluginParameter.Builder<>("resFilename", null, String.class)
//                    .guiName("Residual file name")
//                    .description("The base name for the residual files. _chrname.txt will be appended to each file.")
//                    .dependentOnParameter(residualsAsFile)
//                    .outFile()
//                    .build();

    private static final Logger myLogger = LogManager.getLogger(StepwiseOLSModelFitterPlugin.class);
    private double alpha = 0.05;

    //TODO need to change this to a list

    public StepwiseOLSModelFitterPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public StepwiseOLSModelFitterPlugin() {
        super(null, false);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> datasets = input.getDataOfType(new Class[] { GenotypePhenotype.class });
        if (datasets.size() < 1) {
            String msg = "Error in performFunction: No appropriate dataset selected.";
            throw new IllegalArgumentException(msg);
        }

        //only analyze one dataset
        if (datasets.size() > 1) {
            String msg = "Multiple datasets selected. Only one dataset is allowed.";
            throw new IllegalArgumentException(msg);
        }

        myGenoPheno = (GenotypePhenotype) datasets.get(0).getData();
        datasetName = datasets.get(0).getName();
        myFactorList = myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.factor).stream()
                .map(pa -> pa.name())
                .collect(Collectors.toList());

        if (myFactorList.isEmpty()) {
            List<String> noneList = Arrays.asList(NONE);
            nestingFactor = new PluginParameter(nestingFactor, noneList);
        } else {
            nestingFactor = new PluginParameter(nestingFactor, myFactorList);
        }
    }

    @Override
    protected void postProcessParameters() {
        if (nestMarkers.value() && nestingFactor.value().isEmpty()) {
            if (myFactorList.size() == 1) {
                nestingFactor(myFactorList.get(0));
            } else if (myFactorList.size() > 1) {
                throw new IllegalArgumentException("Nest markers was checked (set to true), but a single factor was not selected to nest markers within. This must be corrected before the analysis will run.");
            }
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        StepwiseOLSModelFitter modelFitter = new StepwiseOLSModelFitter(myGenoPheno, datasetName);
        modelFitter.setEnterlimit(enterlimit.value());
        modelFitter.setExitlimit(exitlimit.value());
        modelFitter.setMaxNumberOfMarkers(maxNumberOfMarkers.value());
        modelFitter.setNested(nestMarkers.value());
        if (nestMarkers.value()) {
            int ndx = myGenoPheno.phenotype().attributeIndexForName(nestingFactor.value());
            if (ndx < 0)
                modelFitter.setNested(false);
            else
                modelFitter.setNestingEffectIndex(ndx);
        }

        modelFitter.setModelType(modelType.value());
        modelFitter.setNumberOfPermutations(numberOfPermutations.value());
        modelFitter.setAlpha(alpha);

        modelFitter.runAnalysis();

        TableReport trResults = modelFitter.getAnovaReport();
        TableReport trEffects = modelFitter.getMarkerEffectReport();
        TableReport trResultsAfterCIScan = modelFitter.getAnovaReportWithCI();
        TableReport trEffectsAfterCIScan = modelFitter.getMarkerEffectReportWithCI();
        TableReport trPermPvalues = modelFitter.getPermutationReport();

        LinkedList<Datum> datumList = new LinkedList<Datum>();
        if (trResults != null)
            datumList.add(new Datum("ANOVA_stepwise_" + datasetName, trResults, "comments"));
        if (trEffects != null)
            datumList.add(new Datum("Marker_estimates_" + datasetName, trEffects, "comments"));
        if (trResultsAfterCIScan != null)
            datumList.add(new Datum("ANOVA_stepwise_After_CI_Scan" + datasetName, trResultsAfterCIScan, "comments"));
        if (trEffectsAfterCIScan != null)
            datumList.add(new Datum("Marker_estimates_After_CI_Scan" + datasetName, trEffectsAfterCIScan, "comments"));
        if (trPermPvalues != null)
            datumList.add(new Datum("Permuted_Pvalues" + datasetName, trPermPvalues, "comments"));

        DataSet myResult = new DataSet(datumList, this);
        fireDataSetReturned(myResult);

        return myResult;
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
        return "Stepwise";
    }

    @Override
    public String getToolTipText() {
        return "Fit multiple markers in a single model (experimental).";
    }

    public String getCitation() {
        String citation = "Written in 2013 by Peter Bradbury and Alex Lipka";
        return citation;
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(StepwiseOLSModelFitterPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public TableReport runPlugin(DataSet input) {
        return (TableReport) performFunction(input).getData(0).getData();
    }

    /**
     * The model selection criteria used to determine which
     * terms enter the model and how many. Value must be one
     * of pvalue, bic, mbic, or aic
     *
     * @return Model type
     */
    public StepwiseOLSModelFitter.MODEL_TYPE modelType() {
        return modelType.value();
    }

    /**
     * Set Model type. The model selection criteria used to
     * determine which terms enter the model and how many.
     * Value must be one of pvalue, bic, mbic, or aic
     *
     * @param value Model type
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin modelType(StepwiseOLSModelFitter.MODEL_TYPE value) {
        modelType = new PluginParameter<>(modelType, value);
        return this;
    }

    /**
     * The enter limit or maximum p-value for which a term
     * can enter the model
     *
     * @return Entry limit
     */
    public Double enterlimit() {
        return enterlimit.value();
    }

    /**
     * Set Entry limit. The enter limit or maximum p-value
     * for which a term can enter the model
     *
     * @param value Entry limit
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin enterlimit(Double value) {
        enterlimit = new PluginParameter<>(enterlimit, value);
        return this;
    }

    /**
     * A term exits the model on a backward step if its p-value
     * is greater than this value
     *
     * @return Exit limit
     */
    public Double exitlimit() {
        return exitlimit.value();
    }

    /**
     * Set Exit limit. A term exits the model on a backward
     * step if its p-value is greater than this value
     *
     * @param value Exit limit
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin exitlimit(Double value) {
        exitlimit = new PluginParameter<>(exitlimit, value);
        return this;
    }

    /**
     * The maximum number of markers that will be fit, if
     * the enter limit is not reached first
     *
     * @return Maximum markers
     */
    public Integer maxNumberOfMarkers() {
        return maxNumberOfMarkers.value();
    }

    /**
     * Set Maximum markers. The maximum number of markers
     * that will be fit, if the enter limit is not reached
     * first
     *
     * @param value Maximum markers
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin maxNumberOfMarkers(Integer value) {
        maxNumberOfMarkers = new PluginParameter<>(maxNumberOfMarkers, value);
        return this;
    }

    /**
     * Should markers be nested within a model factor
     *
     * @return Nest markers
     */
    public Boolean nestMarkers() {
        return nestMarkers.value();
    }

    /**
     * Set Nest markers. Should markers be nested within a
     * model factor
     *
     * @param value Nest markers
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin nestMarkers(Boolean value) {
        nestMarkers = new PluginParameter<>(nestMarkers, value);
        return this;
    }

    /**
     * Nest markers within this factor.
     *
     * @return Nesting factor
     */
    public String nestingFactor() {
        return nestingFactor.value();
    }

    /**
     * Set Nesting factor. Nest markers within this factor.
     *
     * @param value Nesting factor
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin nestingFactor(String value) {
        nestingFactor = new PluginParameter<>(nestingFactor, value);
        return this;
    }

    /**
     * Number of permutations for the model to determine an
     * empirical alpha
     *
     * @return Number of permutations
     */
    public Integer numberOfPermutations() {
        return numberOfPermutations.value();
    }

    /**
     * Set Number of permutations. Number of permutations
     * for the model to determine an empirical alpha
     *
     * @param value Number of permutations
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin numberOfPermutations(Integer value) {
        numberOfPermutations = new PluginParameter<>(numberOfPermutations, value);
        return this;
    }

}
