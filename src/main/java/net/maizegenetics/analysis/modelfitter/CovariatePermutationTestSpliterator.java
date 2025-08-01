package net.maizegenetics.analysis.modelfitter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import org.apache.commons.math3.distribution.FDistribution;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;

/**
 * A Spliterator that uses the covariate from an AdditiveSite to test it against a baseModel and a list of permuted
 * data sets using a fixed effect linear model. This class is used as part of stepwise model fitting to process
 * a list of AdditiveSites. The permutation test is used to determine the entry and exit limits to be used for model fitting.
 */
public class CovariatePermutationTestSpliterator implements Spliterator<double[]> {
    protected List<double[]> myPermutedData;
    protected List<AdditiveSite> mySites;
    protected List<ModelEffect> myBaseModel;
    protected int origin;
    protected final int end;

    public CovariatePermutationTestSpliterator(List<double[]> permutedData,
            List<AdditiveSite> siteList, List<ModelEffect> baseModel) {
        myPermutedData = permutedData;
        mySites = siteList;
        myBaseModel = baseModel;
        origin = 0;
        end = siteList.size();
        int numberOfEffects = baseModel.size();
        DoubleMatrix[][] components = new DoubleMatrix[1][numberOfEffects];
        for (int i = 0; i < numberOfEffects; i++) {
            components[0][i] = myBaseModel.get(i).getX();
        }

    }

    @Override
    public boolean tryAdvance(Consumer<? super double[]> action) {
        if (origin == end)
            return false;
        AdditiveSite as = mySites.get(origin);
        List<ModelEffect> myModel = new ArrayList<ModelEffect>(myBaseModel);
        ModelEffect me;
        CovariateModelEffect cme = new CovariateModelEffect(as.getCovariate());
        myModel.add(cme);

        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, myPermutedData.get(0));
        double dfError = sflm.getResidualSSdf()[1];
        DoubleMatrix G = sflm.getInverseOfXtX();
        FDistribution fdist = new FDistribution(1, dfError);
        double[] pvals =
                myPermutedData.stream().map(d -> DoubleMatrixFactory.DEFAULT.make(d.length, 1, d))
                        .mapToDouble(y -> {
                            double[] yarray = y.to1DArray();
                            int nbase = myBaseModel.size();
                            DoubleMatrix[][] xtyMatrices = new DoubleMatrix[nbase + 1][1];
                            for (int i = 0; i < nbase; i++) xtyMatrices[i][0] = myBaseModel.get(i).getXty(yarray);
                            xtyMatrices[nbase][0] = cme.getXty(yarray);
                            DoubleMatrix Xty = DoubleMatrixFactory.DEFAULT.compose(xtyMatrices);
                            DoubleMatrix beta = G.mult(Xty);
                            double ssTotal = y.crossproduct().get(0, 0);
                            double ssModel = Xty.crossproduct(beta).get(0, 0);
                            double ssError = ssTotal - ssModel;
                            double kb =
                                    beta.get(beta.numberOfRows() - 1, beta.numberOfColumns() - 1);
                            double kgk = G.get(G.numberOfRows() - 1, G.numberOfColumns() - 1);
                            double F = kb * kb / kgk / ssError * dfError;
                            double p = 1 - fdist.cumulativeProbability(F);
                            return p;
                        }).toArray();

        action.accept(pvals);
        origin++;
        return true;
    }

    @Override
    public Spliterator<double[]> trySplit() {
        int numberRemaining = end - origin;
        if (numberRemaining < 50)
            return null;
        int mid = origin + numberRemaining / 2;
        List<AdditiveSite> splitSublist = mySites.subList(origin, mid);
        origin = mid;
        List<double[]> permutedDataCopy =
                myPermutedData.stream().map(d -> Arrays.copyOf(d, d.length)).collect(Collectors.toList());

        return new CovariatePermutationTestSpliterator(permutedDataCopy, splitSublist, myBaseModel);
    }

    @Override
    public long estimateSize() {
        return end - origin;
    }

    @Override
    public int characteristics() {
        return Spliterator.IMMUTABLE + Spliterator.NONNULL + Spliterator.SIZED
                + Spliterator.SUBSIZED;
    }
}
