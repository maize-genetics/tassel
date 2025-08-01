package net.maizegenetics.analysis.modelfitter;

import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.stats.linearmodels.*;
import net.maizegenetics.util.TableReportBuilder;
import org.apache.commons.math3.distribution.FDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.*;
import java.util.concurrent.atomic.LongAdder;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Fits an additive plus dominance model in a stepwise fashion. See StepwiseAdditiveModelFitter for a desription of
 * outputs available.
 */
public class StepwiseAddDomModelFitter extends StepwiseAdditiveModelFitter {

    private static Logger myLogger = LogManager.getLogger(StepwiseAddDomModelFitter.class);

    /**
     *
     * @param genopheno     a GenotypePhenotype object
     * @param datasetName   a name for the genopheno
     * @throws IllegalArgumentException if any phenotype data is missing
     */
    public StepwiseAddDomModelFitter(GenotypePhenotype genopheno, String datasetName) {

        super(genopheno, datasetName);
        markerEffectReportBuilder =
                TableReportBuilder.getInstance("Marker Effects", new String[] { "Trait", "SiteID", "Chr", "Position",
                        "Additive", "Dominance" });
        markerEffectCIReportBuilder =
                TableReportBuilder.getInstance("Marker Effects", new String[] { "Trait", "SiteID", "Chr", "Position",
                        "Additive", "Dominance" });
    }

    @Override
    public void runAnalysis() {
        //load the markers into the appropriate additive site list
        if (useReferenceProbability) {
            mySites =
                    IntStream.range(0, myGenotype.numberOfSites())
                            .mapToObj(s -> {
                                int ntaxa = myPhenotype.numberOfObservations();
                                float[] cov = myGenoPheno.referenceProb(s);
                                return new RefProbAdditiveSite(s, myGenotype.chromosomeName(s), myGenotype.chromosomalPosition(s), myGenotype.siteName(s), modelSelectionCriterion, cov);
                            })
                            .collect(Collectors.toList());
        } else {  // use genotype
            mySites =
                    IntStream.range(0, myGenotype.numberOfSites())
                            .mapToObj(s -> new GenotypeAdditiveSite(s, myGenotype.chromosomeName(s), myGenotype.chromosomalPosition(s), myGenotype.siteName(s),
                                    modelSelectionCriterion, myGenoPheno.genotypeAllTaxa(s), myGenotype.majorAllele(s), myGenotype.majorAlleleFrequency(s)))
                            .collect(Collectors.toList());
        }

        //for each phenotype:
        for (PhenotypeAttribute phenoAttr : dataAttributeList) {
            currentTraitName = phenoAttr.name();

            //build the base model
            List<ModelEffect> myBaseModel = baseModel(phenoAttr);
            myModel = new ArrayList<>(myBaseModel);
            numberOfBaseEffects = myModel.size();

            //call fitModel()
            fitModel();

            //add to reports
            if (createAnovaReport)
                addToAnovaReport(Optional.empty());
            if (createPreScanEffectsReport)
                addToMarkerEffectReport(false);

            //call scanFindCI()
            long start = System.nanoTime();
            List<int[]> intervalList = scanToFindCI();
            myLogger.info(String.format("Rescan in %d ms", (System.nanoTime() - start) / 1000000));

            //created a new scanned model
            myModel = new ArrayList<>(myBaseModel);
            for (int[] interval : intervalList) {
                AdditiveSite as = mySites.get(interval[0]);
                //changed for add-dom model
                myModel.add(new AddPlusDomModelEffect(as,as));
            }
            mySweepFast = new SweepFastLinearModel(myModel, y);

            //add to reports
            if (createAnovaReport)
                addToAnovaReport(Optional.of(intervalList));
            if (createPostScanEffectsReport)
                addToMarkerEffectReport(true);

        }
    }

    @Override
    protected double forwardStep(double prevCriterionValue) {
        //do this in parallel
        //create a stream returning AdditiveSites that have an ordering; select the max
        //criteria can be one of SS, pvalue, aic, bic, mbic (handled by ForwardStepAdditiveSpliterator)

        Spliterator<AdditiveSite> siteEvaluator;
        siteEvaluator = new ForwardStepAddDomSpliterator(mySites, myModel, y);
        LongAdder counter = new LongAdder();
        Optional<AdditiveSite> bestSite =
                StreamSupport.stream(siteEvaluator, true).peek(s -> counter.increment()).max((a, b) -> a.compareTo(b));
        System.out.println(counter.longValue() + " sites evaluated.");
        if (!bestSite.isPresent())
            return Double.NaN;

        ModelEffect nextEffect;
        nextEffect = new AddPlusDomModelEffect(bestSite.get(), bestSite.get());

        myModel.add(nextEffect);
        mySweepFast = new SweepFastLinearModel(myModel, y);
        double[] siteSSdf = mySweepFast.getIncrementalSSdf(myModel.size() - 1);
        double[] errorSSdf = mySweepFast.getResidualSSdf();
        double F, p;
        F = siteSSdf[0] / siteSSdf[1] / errorSSdf[0] * errorSSdf[1];
        p = LinearModelUtils.Ftest(F, siteSSdf[1], errorSSdf[1]);

        boolean addToModel = false;
        double criterionValue = Double.NaN;
        switch (modelSelectionCriterion) {
            case pval:
                criterionValue = p;
                if (p < enterLimit)
                    addToModel = true;
                break;
            case aic:
                criterionValue = aic(errorSSdf[0], y.length, mySweepFast.getFullModelSSdf()[0]);
                if (criterionValue < prevCriterionValue)
                    addToModel = true;
                break;
            case bic:
                criterionValue = bic(errorSSdf[0], y.length, mySweepFast.getFullModelSSdf()[0]);
                if (criterionValue < prevCriterionValue)
                    addToModel = true;
                break;
            case mbic:
                criterionValue =
                        mbic(errorSSdf[0], y.length, mySweepFast.getFullModelSSdf()[0], mySites.size());
                if (criterionValue < prevCriterionValue)
                    addToModel = true;
                break;

        }

        if (addToModel) {
            addToStepsReport(bestSite.get().siteNumber(), mySweepFast, "add", siteSSdf, errorSSdf, F, p);
            return criterionValue;
        }

        addToStepsReport(bestSite.get().siteNumber(), mySweepFast, "stop", siteSSdf, errorSSdf, F, p);
        myModel.remove(myModel.size() - 1);
        mySweepFast = new SweepFastLinearModel(myModel, y);
        return Double.NaN;

    }

    @Override
    protected List<int[]> scanToFindCI() {
        //define an IntFunction that finds interval endpoints
        //the interval is bounded by the first points that when added to the model result in the marginal p of the test site <= alpha
        Function<ModelEffect, int[]> intervalFinder = me -> {
            //scan steps:
            //1. find interval end points
            //2. determine if any point in the interval gives a better model fit (ssmodel) than the original
            //3. if no, return support interval
            //4. if yes, replace the original with that point and rescan then return support interval

            AdditiveSite scanSite = (AdditiveSite) me.getID();
            myLogger.info(String.format("Scanning site %d, %s, pos = %d", scanSite.siteNumber(), myGenotype.chromosome(scanSite.siteNumber()), myGenotype.chromosomalPosition(scanSite.siteNumber())));
            int[] support = findCI(me, myModel);
            List<ModelEffect> baseModel = new ArrayList<>(myModel);
            baseModel.remove(me);
            AdditiveSite bestSite = bestTerm(baseModel, support);
            if (!bestSite.equals(scanSite)) {
                ModelEffect bestEffect;
                bestEffect = new AddPlusDomModelEffect(bestSite, bestSite);
                baseModel.add(bestEffect);
                support = findCI(bestEffect, baseModel);
            }
            return support;
        };

        return myModel.stream().skip(numberOfBaseEffects).parallel().map(intervalFinder).collect(Collectors.toList());
    }

    @Override
    protected double testAddedTerm(int testedTerm, AdditiveSite addedTerm, List<ModelEffect> theModel) {
        List<ModelEffect> testingModel = new ArrayList<>(theModel);

        //changed for add-dom model
        AddPlusDomModelEffect apdme = new AddPlusDomModelEffect(addedTerm, addedTerm);
        testingModel.add(apdme);

        SweepFastLinearModel sflm = new SweepFastLinearModel(testingModel, y);
        sflm.getResidualSSdf();
        double[] residualSSdf = sflm.getResidualSSdf();
        double[] marginalSSdf = sflm.getMarginalSSdf(testedTerm);
        double F = marginalSSdf[0] / marginalSSdf[1] / residualSSdf[0] * residualSSdf[1];

        //debug
        double prob = 1;
        try {
            prob -= (new FDistribution(marginalSSdf[1], residualSSdf[1]).cumulativeProbability(F));
        } catch(Exception e) {
            //do nothing
        }
        return prob;
    }

    @Override
    protected AdditiveSite bestTerm(List<ModelEffect> baseModel, int[] interval) {
        List<AdditiveSite> intervalList = mySites.subList(interval[1], interval[2]);
        PartitionedLinearModel plm =
                new PartitionedLinearModel(baseModel, new SweepFastLinearModel(baseModel, y));
        return intervalList.stream()
                .map(s -> {
                    plm.testNewModelEffect(new AddPlusDomModelEffect(s,s));
                    s.criterionValue(plm.getp());
                    return s;
                })
                //change from >= to <= for add-dom model
                .reduce((a, b) -> a.criterionValue() <= b.criterionValue() ? a : b)
                .get();

    }

    @Override
    public void runPermutationTest() {
        //parallel version of permutation test
        int enterLimitIndex = (int) (permutationAlpha * numberOfPermutations);  //index of percentile to be used for the enter limit

        //create the permutedData
        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, y);
        double[] yhat = sflm.getPredictedValues().to1DArray();
        double[] residuals = sflm.getResiduals().to1DArray();

        BasicShuffler.shuffle(residuals);
        List<double[]> permutedData = Stream.iterate(residuals, BasicShuffler.shuffleDouble())
                .limit(numberOfPermutations)
                .map(a -> {
                    double[] permutedValues = Arrays.copyOf(a, a.length);
                    for (int i = 0; i < a.length; i++)
                        permutedValues[i] += yhat[i];
                    return permutedValues;
                })
                .collect(Collectors.toList());

        //find the minimum p values for each site
        double[] maxP = new double[numberOfPermutations];
        Arrays.fill(maxP, 1.0);
        double[] minP;
        List<double[]> plist = new ArrayList<>();

        minP =
                StreamSupport.stream(new AddDomPermutationTestSpliterator(permutedData, mySites, myModel), true).reduce(maxP, (a, b) -> {
                    int n = a.length;
                    for (int i = 0; i < n; i++) {
                        if (a[i] > b[i])
                            a[i] = b[i];
                    }
                    return a;
                });

        Arrays.sort(minP);
        enterLimit = minP[enterLimitIndex];
        exitLimit = 2 * enterLimit;

        myLogger.info(String.format("Additive + Dominance Permutation results for %s: enterLimit = %1.5e, exitLimit = %1.5e\n", currentTraitName, enterLimit, exitLimit));

        //add values to permutation report : "Trait","p-value"
        Arrays.stream(minP).forEach(d -> permutationReportBuilder.add(new Object[] {
                currentTraitName, new Double(d) }));

    }

    @Override
    protected void addToMarkerEffectReport(boolean CI) {
        //header: "Trait", "SiteID", "Chr", "Position","Additive", "Dominance"
        double[] beta = mySweepFast.getBeta();
        int numberOfEffects = myModel.size();
        int numberOfMarkerEffects = numberOfEffects - numberOfBaseEffects;
        int baseDf = myModel.stream().limit(numberOfBaseEffects).mapToInt(me -> me.getEffectSize()).sum();

        int betaCount = baseDf;
        for (int me = 0; me <  numberOfMarkerEffects; me++) {
            Object[] row = new Object[6];
            int col = 0;
            row[col++] = currentTraitName;
            AddPlusDomModelEffect adModelEffect = (AddPlusDomModelEffect) myModel.get(numberOfBaseEffects + me);
            AdditiveSite mySite = (AdditiveSite) adModelEffect.getID();

            row[col++] = myGenotype.siteName(mySite.siteNumber());
            row[col++] = myGenotype.positions().chromosomeName(mySite.siteNumber());
            row[col++] = myGenotype.positions().get(mySite.siteNumber()).getPosition();
            row[col++] = new Double(beta[betaCount++]);
            if (adModelEffect.getEffectSize() == 2) row[col++] = new Double(beta[betaCount++]);
            else row[col++] = "NA";
            if (CI)
                markerEffectCIReportBuilder.add(row);
            else
                markerEffectReportBuilder.add(row);
        }

    }
}
