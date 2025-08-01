package net.maizegenetics.analysis.modelfitter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.NestedCovariateModelEffect;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;

import org.apache.commons.math3.distribution.FDistribution;

public class NestedCovariatePermutationTestSpliterator implements Spliterator<double[]> {
    private List<double[]> myPermutedData;
    private List<AdditiveSite> mySites;
    private List<ModelEffect> myBaseModel;
    private FactorModelEffect myOuter;
    private int origin;
    private final int end;

    public NestedCovariatePermutationTestSpliterator(List<double[]> permutedData,
            List<AdditiveSite> siteList, List<ModelEffect> baseModel, ModelEffect outerEffect) {
        myPermutedData = permutedData;
        mySites = siteList;
        myBaseModel = baseModel;
        if (!(outerEffect instanceof FactorModelEffect))
            throw new IllegalArgumentException(String.format("The outer effect, %s, is not a factor and cannot be used for nesting.", outerEffect.getID().toString()));
        myOuter = (FactorModelEffect) outerEffect;
        origin = 0;
        end = siteList.size();
        int numberOfEffects = baseModel.size();
        DoubleMatrix[][] components = new DoubleMatrix[1][numberOfEffects];
        for (int i = 0; i < numberOfEffects; i++) {
            components[0][i] = myBaseModel.get(i).getX();
        }
//        baseX = DoubleMatrixFactory.DEFAULT.compose(components);

    }

    @Override
    public boolean tryAdvance(Consumer<? super double[]> action) {
        if (origin == end)
            return false;
        AdditiveSite as = mySites.get(origin);
        List<ModelEffect> myModel = new ArrayList<ModelEffect>(myBaseModel);
        CovariateModelEffect cme = new CovariateModelEffect(as.getCovariate());
        NestedCovariateModelEffect ncme = new NestedCovariateModelEffect(cme, myOuter);
        myModel.add(ncme);

        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, myPermutedData.get(0));
        //Q = (KB)'inv(KGK')(KB)
        //F = Q/rK/msError, where rK = dfCovar
        double dfError = sflm.getResidualSSdf()[1];
        double dfCovar = sflm.getIncrementalSSdf(myBaseModel.size())[1];
        DoubleMatrix G = sflm.getInverseOfXtX();
        int numberOfOuterLevels = myOuter.getNumberOfLevels();
        int ncolBase = myBaseModel.stream().mapToInt(me -> me.getEffectSize()).sum();
        
        int[] betaSelection = IntStream.range(ncolBase, ncolBase + numberOfOuterLevels).toArray();
        DoubleMatrix invKGK = G.getSelection(betaSelection, betaSelection).inverse();
        FDistribution fdist = new FDistribution(dfCovar, dfError);

        double[] pvals =
                myPermutedData.stream().map(d -> DoubleMatrixFactory.DEFAULT.make(d.length, 1, d))
                        .mapToDouble(y -> {
                            double[] yarray = y.to1DArray();
                            int nbase = myBaseModel.size();
                            DoubleMatrix[][] xtyMatrices = new DoubleMatrix[nbase + 1][1];
                            for (int i = 0; i < nbase; i++) xtyMatrices[i][0] = myBaseModel.get(i).getXty(yarray);
                            xtyMatrices[nbase][0] = ncme.getXty(yarray);
                            DoubleMatrix Xty = DoubleMatrixFactory.DEFAULT.compose(xtyMatrices);
                            DoubleMatrix beta = G.mult(Xty);
                            DoubleMatrix Kbeta = beta.getSelection(betaSelection, null);
                            double Q = Kbeta.crossproduct(invKGK.mult(Kbeta)).get(0, 0);
                            double totalSS = y.crossproduct().get(0, 0);
                            double modelSS = Xty.crossproduct(beta).get(0, 0);
                            double residualSS = totalSS - modelSS;
                            double errorMS = residualSS / dfError;
                            double F = Q / dfCovar /errorMS;
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
        myBaseModel.get(0).getCopy();
        List<ModelEffect> baseModelCopy = myBaseModel.stream().map(me -> me.getCopy()).collect(Collectors.toList());
        return new NestedCovariatePermutationTestSpliterator(permutedDataCopy, splitSublist, baseModelCopy, myOuter.getCopy());
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
