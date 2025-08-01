package net.maizegenetics.analysis.modelfitter;

import java.util.Arrays;
import java.util.List;
import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import net.maizegenetics.stats.linearmodels.ModelEffect;

public class ForwardStepSubsettingAdditiveSpliterator extends ForwardStepAdditiveSpliterator {
    private final int[] subset;

    public ForwardStepSubsettingAdditiveSpliterator(List<AdditiveSite> siteList,
            List<ModelEffect> baseModel, double[] y, int[] subset) {
        super(siteList, baseModel, y);
        this.subset = subset;
    }

    private ForwardStepSubsettingAdditiveSpliterator(List<AdditiveSite> siteList,
            List<ModelEffect> baseModel, double[] y, int numberOfSites, int[] subset) {
        super(siteList, baseModel, y, numberOfSites);
        this.subset = subset;
    }

    @Override
    public boolean tryAdvance(Consumer<? super AdditiveSite> action) {
        if (origin == end)
            return false;
        AdditiveSite as = mySites.get(origin);

        switch (as.selectionCriterion()) {
        case pval:
            as.criterionValue(plm.testNewModelEffect(as.getCovariate(subset)));
            break;
        case aic:
            plm.testNewModelEffect(as.getCovariate(subset));
            double rss = plm.getErrorSS();
            as.criterionValue(nobs * Math.log(rss / nobs) + 2 * (baseModeldf + 1));
            break;
        case bic:
            plm.testNewModelEffect(as.getCovariate(subset));
            rss = plm.getErrorSS();
            as.criterionValue(nobs * Math.log(rss / nobs) + Math.log(nobs) * (baseModeldf + 1));
            break;
        case mbic:
            plm.testNewModelEffect(as.getCovariate(subset));
            rss = plm.getErrorSS();
            as.criterionValue(nobs * Math.log(rss / nobs) + Math.log(nobs) * (baseModeldf + 1) + 2
                    * (baseModeldf + 1) * Math.log(nsites / 2.2 - 1));
            break;
        //        case pval:
        //            double modelss = plm.testNewModelEffect(as.getCovariate(subset));
        //            plm.setModelSS(modelss);
        //            as.criterionValue(plm.getp());
        //            break;
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
        int[] subsetCopy = Arrays.copyOf(subset, subset.length);
        return new ForwardStepSubsettingAdditiveSpliterator(splitSublist, baseModelCopy, yCopy, nsites, subsetCopy);
    }

}
