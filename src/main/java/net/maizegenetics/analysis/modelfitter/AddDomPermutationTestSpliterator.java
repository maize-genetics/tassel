package net.maizegenetics.analysis.modelfitter;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import org.apache.commons.math3.distribution.FDistribution;

import java.util.ArrayList;
import java.util.List;
import java.util.Spliterator;
import java.util.function.Consumer;

/**
 * This class calculates p-values for a test of a single SNP with a list of permuted data sets.
 * myPermutedData is a list of double[] arrays, each of which is a data set to be tested.
 * This class extends a parent class which tests an additive only model to test an additive + dominance model instead.
 */
public class AddDomPermutationTestSpliterator extends CovariatePermutationTestSpliterator {

    public AddDomPermutationTestSpliterator(List<double[]> permutedData,
                                               List<AdditiveSite> siteList, List<ModelEffect> baseModel) {
        super(permutedData, siteList, baseModel);
    }

    @Override
    public boolean tryAdvance(Consumer<? super double[]> action) {
        if (origin == end)
            return false;
        AdditiveSite as = mySites.get(origin);
        List<ModelEffect> myModel = new ArrayList<ModelEffect>(myBaseModel);
        ModelEffect effectToAdd = new AddPlusDomModelEffect(as,as);
        myModel.add(effectToAdd);

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
                            xtyMatrices[nbase][0] = effectToAdd.getXty(yarray);
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
}
