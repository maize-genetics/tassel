package net.maizegenetics.analysis.modelfitter;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import org.apache.commons.math3.distribution.FDistribution;

import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;

public class AdditiveResidualForwardRegression extends AbstractForwardRegression {
    double[] residuals;
    
    public AdditiveResidualForwardRegression(GenotypePhenotype data) {
        super(data);
    }

    @Override
    public void fitModel() {
        int maxModelSize = myModel.size() + maxVariants;
        int step = 0;
        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, y);
        residuals = sflm.getResiduals().to1DArray();
        
        while (step < maxModelSize && forwardStepParallel(true, step)) {
           step++;
        }
    }
    
    @Override
    public void fitModelForSubsample(int[] subSample, int iteration) {
        //create myModel from myBaseModel for this subsample
        myModel = myBaseModel.stream().map(me -> me.getSubSample(subSample)).collect(Collectors.toList());
        double[] ySubSample = Arrays.stream(subSample).mapToDouble(i -> y[i]).toArray();
        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, ySubSample);
        residuals = sflm.getResiduals().to1DArray();
        
        int maxModelSize = myModel.size() + maxVariants;
        int step = 0;
        
        while (step < maxModelSize && forwardStepParallel(subSample, true, iteration, step)) {
            step++;
        }
        
    }

    private boolean forwardStepParallel(boolean doParallel, int step) {
        int nsamples = residuals.length;
        ModelEffect meanMe = new FactorModelEffect(new int[nsamples], false, "mean");
        List<ModelEffect> residualModel = Arrays.asList(meanMe);

        AdditiveSite bestSite = StreamSupport.stream(new ForwardStepAdditiveSpliterator(siteList, residualModel, residuals), doParallel)
                .max((a,b) -> a.compareTo(b)).get();
        
        ModelEffect siteEffect = new CovariateModelEffect(bestSite.getCovariate());
        residualModel.add(siteEffect);
        SweepFastLinearModel sflm = new SweepFastLinearModel(residualModel, residuals);
        double[] errorSSdf = sflm.getResidualSSdf();
        double[] siteSSdf = sflm.getIncrementalSSdf(myModel.size() - 1);
        double F,p;
        if (siteSSdf[1] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY || errorSSdf[0] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY) {
            F = Double.NaN;
            p = Double.NaN;
        } else {
            F = siteSSdf[0] / siteSSdf[1] / errorSSdf[0] * errorSSdf[1];
            p = 1 - (new FDistribution(siteSSdf[1], errorSSdf[1]).cumulativeProbability(F));
        }
        
        if (!Double.isNaN(p) && p <= enterLimit) {
            addVariant(bestSite, p, 0, step);
            residuals = sflm.getResiduals().to1DArray();
            return true;
        }
        else return false;
    }
    
    private boolean forwardStepParallel(int[] subset, boolean doParallel, int iteration, int step) {
        int nsamples = residuals.length;
        ModelEffect meanMe = new FactorModelEffect(new int[nsamples], false, "mean");
        List<ModelEffect> residualModel = Arrays.asList(meanMe);

        AdditiveSite bestSite = StreamSupport.stream(new ForwardStepSubsettingAdditiveSpliterator(siteList, residualModel, residuals, subset), doParallel)
                .max((a,b) -> a.compareTo(b)).get();

        ModelEffect siteEffect = new CovariateModelEffect(bestSite.getCovariate(subset));
        residualModel.add(siteEffect);
        SweepFastLinearModel sflm = new SweepFastLinearModel(residualModel, residuals);
        double[] errorSSdf = sflm.getResidualSSdf();
        double[] siteSSdf = sflm.getIncrementalSSdf(myModel.size() - 1);
        double F,p;
        if (siteSSdf[1] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY || errorSSdf[0] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY) {
            F = Double.NaN;
            p = Double.NaN;
        } else {
            F = siteSSdf[0] / siteSSdf[1] / errorSSdf[0] * errorSSdf[1];
            p = 1 - (new FDistribution(siteSSdf[1], errorSSdf[1]).cumulativeProbability(F));
        }
        
        if (!Double.isNaN(p) && p <= enterLimit) {
            //columns in myFittedVariants: "trait","SnpID","Chr","Pos", "p-value", "-log10p"
            addVariant(bestSite, p, iteration, step);
            residuals = sflm.getResiduals().to1DArray();
            return true;
        }
        else return false;
        
    }

}
