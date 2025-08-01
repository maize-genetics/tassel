package net.maizegenetics.analysis.modelfitter;

import java.util.Arrays;
import java.util.List;
import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.PartitionedLinearModel;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;

/**
 * This Spliterator is used to fit a baseModel plus an additive effect covariate to a data set y. It takes an AdditiveSite
 * as input, calculates the test statistic specified by selectionCriterion, and sets its value in the AdditiveSite.
 */
public class ForwardStepAdditiveSpliterator implements Spliterator<AdditiveSite> {
    protected final PartitionedLinearModel plm;
    protected List<AdditiveSite> mySites;
    protected final List<ModelEffect> baseModel;
    protected final double[] y;
    protected int origin;
    protected final int end;
    protected final double baseModeldf;
    protected final int nobs;
    protected int nsites;

    public ForwardStepAdditiveSpliterator(List<AdditiveSite> siteList, List<ModelEffect> baseModel,
            double[] y) {
        SweepFastLinearModel sflm = new SweepFastLinearModel(baseModel, y);
        baseModeldf = sflm.getFullModelSSdf()[1];
        nobs = y.length;
        plm = new PartitionedLinearModel(baseModel, sflm);
        mySites = siteList;
        this.baseModel = baseModel;
        this.y = y;
        origin = 0;
        end = siteList.size();
        nsites = siteList.size();
    }

    protected ForwardStepAdditiveSpliterator(List<AdditiveSite> siteList,
            List<ModelEffect> baseModel, double[] y, int numberOfSites) {
        this(siteList, baseModel, y);
        nsites = numberOfSites;
    }

    @Override
    public boolean tryAdvance(Consumer<? super AdditiveSite> action) {
        if (origin == end)
            return false;
        AdditiveSite as = mySites.get(origin);

        switch (as.selectionCriterion()) {
            case pval:
                as.criterionValue(plm.testNewModelEffect(as.getCovariate()));
                break;
            case aic:
                plm.testNewModelEffect(as.getCovariate());
                double rss = plm.getErrorSS();
                as.criterionValue(nobs * Math.log(rss / nobs) + 2 * (baseModeldf + 1));
                break;
            case bic:
                plm.testNewModelEffect(as.getCovariate());
                rss = plm.getErrorSS();
                as.criterionValue(nobs * Math.log(rss / nobs) + Math.log(nobs) * (baseModeldf + 1));
                break;
            case mbic:
                plm.testNewModelEffect(as.getCovariate());
                rss = plm.getErrorSS();
                as.criterionValue(nobs * Math.log(rss / nobs) + Math.log(nobs) * (baseModeldf + 1) + 2
                        * (baseModeldf + 1) * Math.log(nsites / 2.2 - 1));
                break;
        }

        action.accept(as);
        origin++;
        return true;
    }

    @Override
    public Spliterator<AdditiveSite> trySplit() {
        int numberRemaining = end - origin;
        if (numberRemaining < 50)
            return null;
        int mid = origin + numberRemaining / 2;
        List<AdditiveSite> splitSublist = mySites.subList(origin, mid);
        origin = mid;
        double[] yCopy = Arrays.copyOf(y, y.length);
        List<ModelEffect> baseModelCopy =
                baseModel.stream().map(me -> me.getCopy()).collect(Collectors.toList());
        return new ForwardStepAdditiveSpliterator(splitSublist, baseModelCopy, yCopy, nsites);
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
