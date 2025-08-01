package net.maizegenetics.analysis.modelfitter;

import java.util.Arrays;
import java.util.List;
import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.NestedCovariateModelEffect;

public class ForwardStepNestedAdditiveSpliterator extends ForwardStepAdditiveSpliterator {
    protected FactorModelEffect nestingFactor;

    public ForwardStepNestedAdditiveSpliterator(List<AdditiveSite> siteList,
            List<ModelEffect> baseModel, double[] y, FactorModelEffect nestingFactor) {
        super(siteList, baseModel, y);
        this.nestingFactor = nestingFactor;
    }

    public ForwardStepNestedAdditiveSpliterator(List<AdditiveSite> siteList,
            List<ModelEffect> baseModel, double[] y, int numberOfSites,
            FactorModelEffect nestingFactor) {
        super(siteList, baseModel, y, numberOfSites);
        this.nestingFactor = nestingFactor;
    }

    @Override
    public boolean tryAdvance(Consumer<? super AdditiveSite> action) {
        if (origin == end)
            return false;
        AdditiveSite as = mySites.get(origin);
        ModelEffect me = new NestedCovariateModelEffect(as.getCovariate(), nestingFactor);
        switch (as.selectionCriterion()) {
        case pval:
            plm.testNewModelEffect(me);
            as.criterionValue(plm.getModelSS());
            break;
        case aic:
            plm.testNewModelEffect(me);
            double rss = plm.getErrorSS();
            as.criterionValue(nobs * Math.log(rss / nobs) + 2 * (baseModeldf + 1));
            break;
        case bic:
            plm.testNewModelEffect(me);
            rss = plm.getErrorSS();
            as.criterionValue(nobs * Math.log(rss / nobs) + Math.log(nobs) * (baseModeldf + 1));
            break;
        case mbic:
            plm.testNewModelEffect(me);
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
        return new ForwardStepNestedAdditiveSpliterator(splitSublist, baseModelCopy, yCopy, nsites, nestingFactor);
    }

}
