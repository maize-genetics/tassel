package net.maizegenetics.analysis.modelfitter;

import java.util.List;

import net.maizegenetics.phenotype.Phenotype;

public interface ForwardRegression {

    /**
     * @return  a list of terms in the model. Each item in the list is an an array whose first element is a Position. The second element an Integer equal to the site index of the term.
     * The third element is a Double equal to the p-value for the term.
     */
    List<Object[]> fittedModel();

    /**
     * @param subSample         an array of the sample indices to included in the subsample. Indices can be repeated, for example for bootstrap samples.
     * @param iteration         the iteration or the sequential index of models fit
     */
    void fitModelForSubsample(int[] subSample, int iteration);

    void fitModel();

    public Phenotype phenotype();

    /**
     * This method resets the fitted model to the base model and sets the index of the trait to be fit. 
     * @param traitIndex        index of the trait to be fit. The index is the index within the list of data attributes.
     * @param enterLimit        the maximum p-value allowed for a term entering the model
     * @param maxVariants       the maximum number of variants that will be fit in a model
     */
    void resetModel(int phenotypeIndex, double enterLimit, int maxVariants);

    public static String[] columnLabels() {
        return new String[] { "trait", "iteration", "step", "SnpID", "Chr", "Pos", "p-value",
                "-log10p" };
    }
}