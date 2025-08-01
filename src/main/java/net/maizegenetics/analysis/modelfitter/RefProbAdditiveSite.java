package net.maizegenetics.analysis.modelfitter;

import java.util.Arrays;
import java.util.List;

/**
 * An AdditiveSite that takes a float[] valued covariate as an argument. The covariate should equal the probability that
 * a given allele (major or minor, for example) would be selected at random from that site. It also stores the value of a
 * statistic determined by selectionCriteria which results from fitting a linear model.
 */
public class RefProbAdditiveSite extends AbstractAdditiveSite {

    private static final long serialVersionUID = 2040665024409852166L;
    private int ntaxa;
    private float[] cov;
    private int[] taxaIndex = null;

    public RefProbAdditiveSite(int site, String chr, int pos, String id,
            CRITERION selectionCriteria, float[] covariate) {
        super(site, chr, pos, id, selectionCriteria);
        cov = covariate;
        ntaxa = cov.length;
    }

    @Override
    public double[] getCovariate() {
        double[] dcov = new double[ntaxa];
        if (taxaIndex == null) {
            for (int i = 0; i < ntaxa; i++)
                dcov[i] = cov[i];
            return dcov;
        } else {
            for (int i = 0; i < ntaxa; i++)
                dcov[i] = cov[taxaIndex[i]];
            return dcov;
        }
    }

    @Override
    public double[] getCovariate(int[] subset) {
        if (taxaIndex == null)
            return Arrays.stream(subset).mapToDouble(i -> cov[i]).toArray();
        return Arrays.stream(subset).mapToDouble(i -> cov[taxaIndex[i]]).toArray();
    }

    @Override
    public void reindexTaxa(int[] taxaIndex, List<Integer> uniqueTaxa) {
        this.taxaIndex = taxaIndex;
    }

    @Override
    public AdditiveSite copy() {
        return new RefProbAdditiveSite(siteIndex, chrName, position, name, selectionCriterion, Arrays.copyOf(cov, cov.length));
    }
}
