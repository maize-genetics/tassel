package net.maizegenetics.analysis.modelfitter;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.Arrays;

/**
 * A ModelEffect that takes an AdditiveSite as an argument and uses the additiveCovariate
 * from that object to compute a covariate representing dominance. In general, the AdditiveSite will be supplied as the id,
 * although that is not a requirement. The dominance values are calculated as 1 - abs(x - 1) where x is the addivitive value.
 * The values of x are expected to be 0, 1, 2, where 0 and 1 are homozygous genotypes and 1 is the het.
 *
 */
public class AddPlusDomModelEffect implements ModelEffect {

    private AdditiveSite id;
    private final double[] addCovariate;
    private final double[] domCovariate;
    private final CovariateModelEffect addModelEffect;
    private final CovariateModelEffect domModelEffect;
    private final DoubleMatrix X;
    private final double delta = 1e-8;
    public static double MIN_HETS = 20;
    public static boolean IMPUTE_DOM = true;

    public AddPlusDomModelEffect(AdditiveSite id, AdditiveSite addSite) {
        this(id, addSite.getCovariate());
    }

    public AddPlusDomModelEffect(AdditiveSite id, double[] additiveCovariate) {
        this.id = id;
        addCovariate = additiveCovariate;
        double lower = 1 - delta;
        double upper = 1 + delta;

        if (IMPUTE_DOM) domCovariate = Arrays.stream(addCovariate).map(add -> 1.0 - Math.abs(add - 1)).toArray();
        else domCovariate = Arrays.stream(addCovariate).map(add -> (add > lower && add < upper) ? 1.0 : 0.0).toArray();
        double domSum = Arrays.stream(domCovariate).sum();

        addModelEffect = new CovariateModelEffect(addCovariate, id);
        if (domSum >= MIN_HETS) {
            domModelEffect = new CovariateModelEffect(domCovariate, id);
            X = addModelEffect.getX().concatenate(domModelEffect.getX(), false);
        } else {
            domModelEffect = null;
            X = addModelEffect.getX();
        }
    }

    @Override
    public Object getID() {
        return id;
    }

    @Override
    public void setID(Object id) {
        this.id = (AdditiveSite) id;
    }

    @Override
    public int getSize() {
        return addCovariate.length;
    }

    @Override
    public DoubleMatrix getX() {
        return X;
    }

    @Override
    public DoubleMatrix getXtX() {
        return getX().crossproduct();
    }

    @Override
    public DoubleMatrix getXty(double[] y) {
        return getX().crossproduct(DoubleMatrixFactory.DEFAULT.make(y.length, 1, y));
    }

    @Override
    public DoubleMatrix getyhat(DoubleMatrix beta) {
        return getX().crossproduct(beta);
    }

    @Override
    public DoubleMatrix getyhat(double[] beta) {
        return getyhat(DoubleMatrixFactory.DEFAULT.make(1, beta.length, beta));
    }

    @Override
    public int[] getLevelCounts() {
        if (domModelEffect == null) {
            return new int[]{addCovariate.length};
        } else return new int[]{addCovariate.length, domCovariate.length};
    }

    @Override
    public int getNumberOfLevels() {
        if (domModelEffect == null) {
            return 1;
        } else return 2;

    }

    @Override
    public int getEffectSize() {
        return getNumberOfLevels();
    }

    @Override
    public ModelEffect getCopy() {
        return new AddPlusDomModelEffect(id.copy(), Arrays.copyOf(addCovariate, addCovariate.length));
    }

    @Override
    public ModelEffect getSubSample(int[] sample) {
        double[] subsample = new double[sample.length];
        for (int i = 0; i < sample.length; i++) {
            subsample[i] = addCovariate[sample[i]];
        }
        return new AddPlusDomModelEffect(id, subsample);
    }
}
