package net.maizegenetics.analysis.modelfitter;

import java.util.Arrays;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import org.apache.commons.math3.distribution.FDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;

public class AdditiveModelForwardRegression extends AbstractForwardRegression {
    private static Logger myLogger = LogManager.getLogger(AdditiveModelForwardRegression.class);
    double highestSS;
    int bestSite;
    int maxThreads;
    ForkJoinPool myPool = null;
    
    public AdditiveModelForwardRegression(GenotypePhenotype data, int maxThreads) {
        super(data);
        this.maxThreads = maxThreads;
        myPool = new ForkJoinPool(maxThreads);
    }

    public AdditiveModelForwardRegression(String serialFilename, Phenotype pheno, int maxThreads) {
        super(serialFilename, pheno);
        this.maxThreads = maxThreads;
        myPool = new ForkJoinPool(maxThreads);
    }
    
    public AdditiveModelForwardRegression(GenotypePhenotype data) {
        super(data);
        maxThreads = -1;
    }

    public AdditiveModelForwardRegression(String serialFilename, Phenotype pheno) {
        super(serialFilename, pheno);
        maxThreads = -1;
    }


    @Override
    public void fitModel() {
        int maxModelSize = myModel.size() + maxVariants;
        int step = 0;
        if (maxThreads > 0) {
            while (forwardStepParallelUseMaxThreads(true, step++) && myModel.size() < maxModelSize)
                ;
        } else {
            while (forwardStepParallel(true, step++) && myModel.size() < maxModelSize)
                ;
        }
        
    }

    @Override
    public void fitModelForSubsample(int[] subSample, int iteration) {

        //create myModel from myBaseModel for this subsample
        myModel =
                myBaseModel.stream().map(me -> me.getSubSample(subSample)).collect(Collectors.toList());
        double[] original = y;
        y = Arrays.stream(subSample).mapToDouble(i -> original[i]).toArray();

        int maxModelSize = myModel.size() + maxVariants;
        int step = 0;
        if (maxThreads > 0) {
            while (forwardStepParallelUseMaxThreads(subSample, true, iteration, step++)
                    && myModel.size() < maxModelSize)
                ;
        } else {
            while (forwardStepParallel(subSample, true, iteration, step++)
                    && myModel.size() < maxModelSize)
                ;
        }

        y = original;
    }

    private boolean forwardStepParallel(boolean doParallel, int step) {

        AdditiveSite bestSite =
                StreamSupport.stream(new ForwardStepAdditiveSpliterator(siteList, myModel, y), doParallel)
                        .max((a, b) -> a.compareTo(b)).get();

        ModelEffect siteEffect = new CovariateModelEffect(bestSite.getCovariate());
        myModel.add(siteEffect);
        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, y);
        double[] errorSSdf = sflm.getResidualSSdf();
        double[] siteSSdf = sflm.getIncrementalSSdf(myModel.size() - 1);
        double F, p;
        if (siteSSdf[1] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY
                || errorSSdf[0] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY) {
            F = Double.NaN;
            p = Double.NaN;
        } else {
            F = siteSSdf[0] / siteSSdf[1] / errorSSdf[0] * errorSSdf[1];
            p = 1 - (new FDistribution(siteSSdf[1], errorSSdf[1]).cumulativeProbability(F));
        }

        if (!Double.isNaN(p) && p <= enterLimit) {
            addVariant(bestSite, p, 0, step);
            return true;
        }
        else
            return false;
    }

    private boolean forwardStepParallelUseMaxThreads(boolean doParallel, int step) {
        
        AdditiveSite bestSite;
        try {
            bestSite = myPool.submit(() -> StreamSupport.stream(new ForwardStepAdditiveSpliterator(siteList, myModel, y), doParallel)
                    .max((a, b) -> a.compareTo(b)).get()).get();
            ModelEffect siteEffect = new CovariateModelEffect(bestSite.getCovariate());
            myModel.add(siteEffect);
            SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, y);
            double[] errorSSdf = sflm.getResidualSSdf();
            double[] siteSSdf = sflm.getIncrementalSSdf(myModel.size() - 1);
            double F, p;
            if (siteSSdf[1] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY
                    || errorSSdf[0] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY) {
                F = Double.NaN;
                p = Double.NaN;
            } else {
                F = siteSSdf[0] / siteSSdf[1] / errorSSdf[0] * errorSSdf[1];
                p = 1 - (new FDistribution(siteSSdf[1], errorSSdf[1]).cumulativeProbability(F));
            }

            if (!Double.isNaN(p) && p <= enterLimit) {
                addVariant(bestSite, p, 0, step);
                return true;
            }
            else
                return false;
            
        } catch (InterruptedException | ExecutionException e) {
            myLogger.error("Thread execution failed in forwardStepParallelUseMaxThreads");
            e.printStackTrace();
            return false;
        }
    }

    private boolean forwardStepParallel(int[] subset, boolean doParallel, int iteration, int step) {

        AdditiveSite bestSite =
                StreamSupport.stream(new ForwardStepSubsettingAdditiveSpliterator(siteList, myModel, y, subset), doParallel)
                        .max((a, b) -> a.compareTo(b)).get();

        ModelEffect siteEffect = new CovariateModelEffect(bestSite.getCovariate(subset));
        myModel.add(siteEffect);
        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, y);
        double[] errorSSdf = sflm.getResidualSSdf();
        double[] siteSSdf = sflm.getIncrementalSSdf(myModel.size() - 1);
        double F, p;
        if (siteSSdf[1] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY
                || errorSSdf[0] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY) {
            F = Double.NaN;
            p = Double.NaN;
        } else {
            F = siteSSdf[0] / siteSSdf[1] / errorSSdf[0] * errorSSdf[1];
            p = 1 - (new FDistribution(siteSSdf[1], errorSSdf[1]).cumulativeProbability(F));
        }

        if (!Double.isNaN(p) && p <= enterLimit) {
            //columns in myFittedVariants: "trait","SnpID","Chr","Pos", "p-value", "-log10p"
            addVariant(bestSite, p, iteration, step);
            return true;
        }
        else
            return false;

    }

    private boolean forwardStepParallelUseMaxThreads(int[] subset, boolean doParallel, int iteration, int step) {
        
        AdditiveSite bestSite;
        try {
            bestSite = myPool.submit(() -> StreamSupport.stream(new ForwardStepSubsettingAdditiveSpliterator(siteList, myModel, y, subset), doParallel)
                    .max((a, b) -> a.compareTo(b)).get()).get();
            ModelEffect siteEffect = new CovariateModelEffect(bestSite.getCovariate(subset));
            myModel.add(siteEffect);
            SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, y);
            double[] errorSSdf = sflm.getResidualSSdf();
            double[] siteSSdf = sflm.getIncrementalSSdf(myModel.size() - 1);
            double F, p;
            if (siteSSdf[1] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY
                    || errorSSdf[0] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY) {
                F = Double.NaN;
                p = Double.NaN;
            } else {
                F = siteSSdf[0] / siteSSdf[1] / errorSSdf[0] * errorSSdf[1];
                p = 1 - (new FDistribution(siteSSdf[1], errorSSdf[1]).cumulativeProbability(F));
            }

            if (!Double.isNaN(p) && p <= enterLimit) {
                //columns in myFittedVariants: "trait","SnpID","Chr","Pos", "p-value", "-log10p"
                addVariant(bestSite, p, iteration, step);
                return true;
            }
            else
                return false;
            
        } catch (InterruptedException | ExecutionException e) {
            myLogger.error("Thread execution failed in forwardStepParallelUseMaxThreads");
            e.printStackTrace();
            return false;
        }

    }

    @Override
    protected void finalize() throws Throwable {
        if (myPool != null) myPool.shutdown();
        super.finalize();
    }

}
