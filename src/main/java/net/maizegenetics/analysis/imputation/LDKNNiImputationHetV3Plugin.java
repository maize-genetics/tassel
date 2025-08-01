package net.maizegenetics.analysis.imputation;

import com.google.common.collect.*;
import com.google.common.primitives.Bytes;
import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.analysis.popgen.LDResult;
import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.Tuple;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static net.maizegenetics.dna.snp.GenotypeTable.RARE_DIPLOID_ALLELE;
import static net.maizegenetics.dna.snp.GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
import static net.maizegenetics.dna.snp.GenotypeTableUtils.*;

/**
 * This imputation algorithm uses LD to identify good predictors for each SNP,
 * and then uses the high LD SNPs to identify K- Nearest Neighbors.
 * The genotype is called with a weighted mode of the KNNs.
 *
 * @author Daniel Money (Developer and developed the algorithm)
 * @author Ed Buckler (assisted in conversion to TASSEL)
 */
public class LDKNNiImputationHetV3Plugin extends AbstractPlugin {
    private static int minAlleleDivisorForLDMin=3;  //Sets a threshold for minimum of minor alleles in the LD test

    private PluginParameter<Integer> highLDSSites = new PluginParameter.Builder<>("highLDSSites", 30, Integer.class)
            .range(Range.closed(2, 2000))
            .guiName("High LD Sites")
            .description("Number of sites in high LD to use in imputation")
            .build();

    private PluginParameter<Integer> knnTaxa = new PluginParameter.Builder<>("knnTaxa", 10, Integer.class)
            .range(Range.closed(2, 200))
            .guiName("Number of nearest neighbors")
            .description("Number of neighbors to use in imputation")
            .build();

    private PluginParameter<Integer> maxDistance = new PluginParameter.Builder<>("maxLDDistance", -1, Integer.class)
            .guiName("Max distance between site to find LD")
            .description("Maximum physical distance between sites to search for LD (-1 for no distance cutoff - unlinked chromosomes will be tested)")
            .build();

    private PluginParameter<Double> maxDistanceFromNN = new PluginParameter.Builder<>("maxDistanceFromNN", 0.1, Double.class)
            .guiName("Maximum distance from Nearest Neighbor")
            .description("Maximum distance from Nearest Neighbor")
            .build();

    private PluginParameter<Double> duplicateHetsThreshold = new PluginParameter.Builder<>("hetThreshold", 0.05, Double.class)
            .guiName("HeterozygousThreshold")
            .description("Threshold for defining heterozygous sites")
            .build();

    private PluginParameter<Double> automaticMajorMAF = new PluginParameter.Builder<>("autoMajorMAF", 0.01, Double.class)
            .guiName("Automatic MajorGenotype if MAF")
            .description("Set to Major genotype if no imputation result and MAF is below threshold")
            .build();

    private PluginParameter<Double> minCoverageForDonors = new PluginParameter.Builder<>("minCoverageForDonors", 0.5, Double.class)
            .guiName("Minimum coverage for donors and LD")
            .description("Minimum coverage for donor genotype and LD calculation")
            .build();

    private PluginParameter<Double> minCallBestGenoRatio = new PluginParameter.Builder<>("minCallBestGenoRatio", 10.0, Double.class)
            .guiName("Minimum support ratio for best genotype")
            .description("Minimum ratio between best and second best genotype to make a call")
            .build();

    private PluginParameter<Integer> maxCores = new PluginParameter.Builder<>("maxCores", 8, Integer.class)
            .guiName("Maximum number of cores for processing")
            .description("Maximum number of cores to be used for processing")
            .build();




    private static final Logger myLogger = LogManager.getLogger(LDKNNiImputationHetV3Plugin.class);

    public LDKNNiImputationHetV3Plugin() {
        super(null, false);
    }

    public LDKNNiImputationHetV3Plugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }


    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
        if (alignInList.size() != 1) {
            throw new IllegalArgumentException("LDKNNiImputationPlugin: preProcessParameters: Please select one Genotype Table.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        // Load in the genotype table
        Datum genoDatum=input.getDataOfType(GenotypeTable.class).get(0);
        GenotypeTable genotypeTable = (GenotypeTable)genoDatum.getData();
        System.out.println("genotypeTable = " + genotypeTable.toString());

        StatsOnSites statsOnSites=new StatsOnSites();

        final int numberSites=genotypeTable.numberOfSites();
        final int numberTaxa=genotypeTable.numberOfTaxa();
        final double[] taxaCoverage=IntStream.range(0,numberTaxa)
                .mapToDouble(taxaIndex -> (double)genotypeTable.totalNonMissingForTaxon(taxaIndex)/(double)numberSites)
                .toArray();
        final double[] siteCoverage=IntStream.range(0, numberSites)
                .mapToDouble(siteIndex -> (double) genotypeTable.totalNonMissingForSite(siteIndex) / (double) numberTaxa)
                .toArray();
        //site heterozygosity weighted depth of coverage.  Low coverage will always appear homozygous.
        final double[] weightHets=IntStream.range(0, numberSites).parallel().mapToDouble(siteIndex -> {
            double weightedHets=0, weightHomozygous=0;
            for (int taxaIndex = 0; taxaIndex < genotypeTable.numberOfTaxa(); taxaIndex++) {
                byte g=genotypeTable.genotype(taxaIndex,siteIndex);
                if(g==UNKNOWN_DIPLOID_ALLELE) continue;
                if(isHeterozygous(g)) {weightedHets+=taxaCoverage[taxaIndex];}
                else weightHomozygous+=taxaCoverage[taxaIndex];
            }
            return weightedHets/(weightedHets+weightHomozygous);
            }).toArray();

        GenotypeTableBuilder incSiteBuilder = GenotypeTableBuilder.getSiteIncremental(genotypeTable.taxa());
        int[][] chrDivisions=FILLINFindHaplotypesPlugin.divideChromosome(genotypeTable,2000,true);
        //Start imputing site by site
        long time=System.nanoTime();
        LongAdder sitesDone=new LongAdder();

        IntStream.range(0,maxCores()).parallel().forEach(core -> {
            IntStream.range(0, chrDivisions.length).
                filter(chrDivisionIndex -> chrDivisionIndex % maxCores() == core)
                .forEach(chrDivisionIndex ->
                {   Map<Position,Integer> highLDSiteCnt=new TreeMap<>();
                    IntStream.range(chrDivisions[chrDivisionIndex][0],chrDivisions[chrDivisionIndex][1]).forEach(posIndex -> {
                        PositionList pl = getHighLDPositionList(genotypeTable, posIndex, highLDSSites(), weightHets, siteCoverage);
                        pl.stream().forEach(p -> highLDSiteCnt.merge(p, 1, (oldValue, newValue) -> oldValue + newValue));
                    });
                    PositionList positionList=highLDSiteCnt.entrySet().stream().filter(entry -> entry.getValue()>100).map(entry -> entry.getKey()).collect(PositionList.collectReorder());
                    GenotypeTable ldGenoTable = GenotypeTableBuilder.getGenotypeCopyInstance(FilterGenotypeTable.getInstance(genotypeTable, positionList));
                    double[][] ldDM= new double[ldGenoTable.numberOfTaxa()][ldGenoTable.numberOfTaxa()];
                    for (int i = 0; i < ldGenoTable.numberOfTaxa(); i++) {
                        for (int j = 0; j < i; j++) {
                            ldDM[i][j]=ldDM[j][i]=dist(ldGenoTable.genotypeAllSites(i), ldGenoTable.genotypeAllSites(j),10)[0];
                        }
                        ldDM[i][i]=(double) ldGenoTable.totalNonMissingForTaxon(i) / (double) ldGenoTable.numberOfSites();  //set the diagonal to coverage

                    }
                    List<Multimap<Double, Integer>> closeTaxa = IntStream.range(0, ldGenoTable.numberOfTaxa())
                            .mapToObj(taxonIndex -> getClosestNonMissingTaxa(taxonIndex, ldDM, knnTaxa()))
                            .collect(Collectors.toList());
                    for (int taxon = 0; taxon < ldGenoTable.numberOfTaxa(); taxon++) {
                        Multimap<Double, Integer> closeGenotypes1 = getClosestNonMissingTaxa(taxon, ldDM, knnTaxa());
                    }
                    IntStream.range(chrDivisions[chrDivisionIndex][0],chrDivisions[chrDivisionIndex][1])
                            .forEach(posIndex -> {
                                sitesDone.increment();
                                Position position = genotypeTable.positions().get(posIndex);

                                double maf = genotypeTable.minorAlleleFrequency(posIndex);
                                byte majorAllele = genotypeTable.majorAllele(posIndex);
                                byte minorAllele = genotypeTable.minorAllele(posIndex);
                                if (genotypeTable.isPolymorphic(posIndex) == false) {
                                    System.out.println(Arrays.toString(genotypeTable.genotypeAllTaxa(posIndex)));
                                }
                                StatsOnOneSite statsOnOneSite = new StatsOnOneSite(genotypeTable.minorAlleleFrequency(posIndex),
                                        weightHets[posIndex], genotypeTable.isPolymorphic(posIndex), positionList.size(), majorAllele, minorAllele);

                                byte[] currGenos = genotypeTable.genotypeAllTaxa(posIndex);
                                byte[] impGenos = new byte[currGenos.length];

                                for (int taxon = 0; taxon < currGenos.length; taxon++) {
                                    Multimap<Double, Byte> closeGenotypes = getClosestNonMissingCalls(currGenos,
                                            closeTaxa.get(taxon));
                                    impGenos[taxon] = (closeGenotypes.isEmpty()) ? UNKNOWN_DIPLOID_ALLELE : impute(closeGenotypes, highLDSSites());
                                    if (impGenos[taxon] == UNKNOWN_DIPLOID_ALLELE && maf < automaticMajorMAF()) {
                                        impGenos[taxon] = getDiploidValue(majorAllele, majorAllele);  //set to major genotype for rare allele
                                    }
//                                    System.out.println(closeGenotypes.size()+":"+impGenos[taxon]);  //problem closeGenotypes of large size are not imputing
                                    statsOnOneSite.updateStats(currGenos[taxon], impGenos[taxon]);

                                }
                                byte[] resolvedGeno = resolveGenotypes(currGenos, impGenos, statsOnOneSite.hetFreq() > duplicateHetsThreshold());

                                GeneralPosition.Builder gpb = new GeneralPosition.Builder(position)
                                        .addAnno("ImpHomoAccuracy", statsOnOneSite.homozygousAcc())
                                        .addAnno("ImpMinorAccuracy", statsOnOneSite.minorAcc());
                                if (weightHets[posIndex] > duplicateHetsThreshold()) gpb.addAnno("DUP");
                                incSiteBuilder.addSite(gpb.build(), resolvedGeno);

                                statsOnSites.addStats(statsOnOneSite.hetFreq(), statsOnOneSite.getCnts());

                                if ((posIndex + 1) % 100 == 0) {
                                    fireProgress(33 + ((int) (66 * sitesDone.longValue()) / genotypeTable.numberOfSites()));
                                    System.out.println(sitesDone.longValue() + ": ms/site" + ((System.nanoTime() - time) / 1_000_000) / sitesDone.longValue());
                                    statsOnSites.printToStdOut();
                                    System.out.println(reportingParameters() + statsOnSites.homozygousAcc(2) + "\t" + statsOnSites.recallPowerOfHomozgyous(2));
                                }
                            });
                });
        });
        System.out.println("Final:"+sitesDone.longValue() + ": ms/site" + ((double)(System.nanoTime() - time) / 1_000_000D) / (double)sitesDone.longValue());
        statsOnSites.printToStdOut();
        System.out.println("Stats\thighLDSSites\tknnTaxa\tmaxDistance\tmaxDistanceFromNN\tduplicateHetsThreshold\tHomozygousAcc\tHomozygousPower");
        System.out.printf("Final\t%s\t%.4g\t%.4g\n",reportingParameters(),statsOnSites.homozygousAcc(2),statsOnSites.recallPowerOfHomozgyous(2));
        GenotypeTable impGenotypeTable = incSiteBuilder.build();

        return new DataSet(new Datum(genoDatum.getName()+"_KNNimp",impGenotypeTable,"Imputed genotypes by KNN imputation"),this);
    }

    private byte[] resolveGenotypes(byte[] currGenos, byte[] impGenos, boolean isDuplicated) {
        byte[] resolvedGeno=new byte[currGenos.length];
        for (int taxaIndex = 0; taxaIndex < resolvedGeno.length; taxaIndex++) {
            if(currGenos[taxaIndex]==impGenos[taxaIndex]) {
                resolvedGeno[taxaIndex]=currGenos[taxaIndex];
            } else if(currGenos[taxaIndex]==UNKNOWN_DIPLOID_ALLELE) {
                resolvedGeno[taxaIndex]=impGenos[taxaIndex];
            } else if(isDuplicated) {
                resolvedGeno[taxaIndex]=currGenos[taxaIndex];
            } else {
                //todo one could consider dealing with the het disagreement and making them homozygous
                resolvedGeno[taxaIndex]=currGenos[taxaIndex];
            }
        }
        return resolvedGeno;
    }

    private String reportingParameters() {
        StringBuilder sb=new StringBuilder();
        sb.append(highLDSSites()+"\t");
        sb.append(knnTaxa()+"\t");
        sb.append(maxDistance()+"\t");
        sb.append(maxDistanceFromNN()+"\t");
        sb.append(duplicateHetsThreshold()+"\t");
        return sb.toString();
    }


//    /**
//     * Create a multimap of distances between the target taxon with the closest observed genotypes.  Distance is calculated
//     * between the target taxon and all other taxa for the ldGenoTable (subset of high LD sites).
//     * @param inputTaxon Taxon being imputed
//     * @param genotypeTable Master genotype table
//     * @param ldGenoTable Subset genotype table with only the high LD sites
//     * @param targetPosition The position in the master genotype table being imputed.
//     * @param numberOfTaxa number of genotypes to retain
//     * @return Map of distance to genotype call for the target position
//     */
//    private Multimap<Double, Byte> getClosestNonMissingTaxa(Taxon inputTaxon, GenotypeTable genotypeTable, GenotypeTable ldGenoTable,
//                                                            Position targetPosition, double[] inputCoverage, int numberOfTaxa) {
//        final int targetPosIdx = genotypeTable.positions().indexOf(targetPosition);
//        final int inputTaxonIdx = genotypeTable.taxa().indexOf(inputTaxon);
//        byte[] inputTaxonGenotypes = ldGenoTable.genotypeAllSites(inputTaxonIdx);
//
//
//        MinMaxPriorityQueue<Tuple<Double, Byte>> topTaxa = IntStream.range(0, genotypeTable.numberOfTaxa())
//                .filter(closeTaxonIdx -> closeTaxonIdx != inputTaxonIdx)  //do not test itself
//                .filter(closeTaxonIdx -> inputCoverage[closeTaxonIdx] > minCoverageForDonors())
//                .filter(closeTaxonIdx -> inputCoverage[closeTaxonIdx] * inputCoverage[inputTaxonIdx] * (double) ldGenoTable.numberOfSites() > 10)  //skip tests with
//                .filter(closeTaxonIdx -> genotypeTable.genotype(closeTaxonIdx, targetPosIdx) != GenotypeTable.UNKNOWN_DIPLOID_ALLELE)
//                .mapToObj(closeTaxonIdx -> new Tuple<>(dist(inputTaxonGenotypes, ldGenoTable.genotypeAllSites(closeTaxonIdx), 10)[0], genotypeTable.genotype(closeTaxonIdx, targetPosIdx)))  //calculate the distance
//                .filter(distanceTaxon -> !Double.isNaN(distanceTaxon.x))  //skip is too few sites (<10 results in NaN)
//                .filter(distanceTaxon -> distanceTaxon.x <= maxDistanceFromNN())  //skip greater than max distance
//                .collect(Collectors.toCollection(() -> MinMaxPriorityQueue.maximumSize(numberOfTaxa).create()));
//        final Multimap<Double, Byte> distGenoMap = ArrayListMultimap.create();
//        topTaxa.stream().forEach(distGeno -> distGenoMap.put(distGeno.x, distGeno.y));
//        return distGenoMap;
//    }

    /**
     * Create a multimap of distances between the target taxon with the closest observed genotypes.  Distance is calculated
     * between the target taxon and all other taxa for the ldGenoTable (subset of high LD sites).
     * @param inputTaxonIdx Taxon being imputed
     * @param numberOfTaxa number of genotypes to retain
     * @return Map of distance to genotype call for the target position
     */
    private Multimap<Double, Integer> getClosestNonMissingTaxa(int inputTaxonIdx, double[][] distCoverageMatrix, int numberOfTaxa) {
        System.exit(0);
//        MinMaxPriorityQueue<Tuple<Double, Integer>> topTaxa = IntStream.range(0, distCoverageMatrix.length)
//                .filter(closeTaxonIdx -> closeTaxonIdx != inputTaxonIdx)  //do not test itself
//                .filter(closeTaxonIdx -> distCoverageMatrix[closeTaxonIdx][closeTaxonIdx] > minCoverageForDonors())
//                .mapToObj(closeTaxonIdx -> new Tuple<>(distCoverageMatrix[inputTaxonIdx][closeTaxonIdx], closeTaxonIdx))  //lookup the distance
//                .filter(distanceTaxon -> !Double.isNaN(distanceTaxon.x))  //skip is too few sites (<10 results in NaN)
//                .filter(distanceTaxon -> distanceTaxon.x <= maxDistanceFromNN())  //skip greater than max distance
//                .collect(Collectors.toCollection(() -> MinMaxPriorityQueue.maximumSize(numberOfTaxa).create()));
//        double[] distances=new double[topTaxa.size()];
//        int[] taxaIndices=new int[topTaxa.size()];
//        AtomicInteger ai=new AtomicInteger();
//        topTaxa.stream().forEach(distGeno -> {
//            distances[ai.get()]=distGeno.x;
//            taxaIndices[ai.get()]= distGeno.y;
//        }
//
//        Tuple<double[],int[]> distGenoArray=new Tuple<>(distances,taxaIndices);  //todo
//        final Multimap<Double, Integer> distGenoMap = ArrayListMultimap.create();
//        topTaxa.stream().forEach(distGeno -> distGenoMap.put(distGeno.x, distGeno.y));
//        return distGenoMap;
        return null;
    }

    /**
     * Create a multimap of distances between the target taxon with the closest observed genotypes.  Distance is calculated
     * between the target taxon and all other taxa for the ldGenoTable (subset of high LD sites).
     * @return Map of distance to genotype call for the target position
     */
    private Multimap<Double, Byte> getClosestNonMissingCalls(byte[] genosForSite, Multimap<Double, Integer> closestTaxa) {
        final Multimap<Double, Byte> distGenoMap = ArrayListMultimap.create(closestTaxa.size(),10);
        closestTaxa.entries().stream()
                .forEach(entry -> distGenoMap.put(entry.getKey(), genosForSite[entry.getValue()]));
        return distGenoMap;
    }

    /**
     * Calculates high LD sites with the target site.  Key issues is that we use high coverage taxa , and
     * we ensure
     *
     */
    private PositionList getHighLDPositionList(GenotypeTable genotypeTable, int posIndex, int numberOfSNPs,
                                               double[] hetFreq, double[] coverage) {
        MinMaxPriorityQueue<LDResult> highestLD = MinMaxPriorityQueue.orderedBy(LDResult.byR2Ordering.reverse()).maximumSize(numberOfSNPs).create();
        int minorAlleleCnt = (int) genotypeTable.allelePresenceForAllTaxa(posIndex, WHICH_ALLELE.Minor).cardinality();
        minorAlleleCnt = minorAlleleCnt / minAlleleDivisorForLDMin;  //TODO consider whether this is needed
        if(minorAlleleCnt<2) minorAlleleCnt=2;

        for (int site2 = 0; site2 < genotypeTable.numberOfSites(); site2++) {
            if (posIndex == site2) continue;
            if(coverage[site2]<minCoverageForDonors()) continue;
            if (maxDistance() > -1 && Math.abs(genotypeTable.chromosomalPosition(posIndex) - genotypeTable.chromosomalPosition(site2)) > maxDistance()) {
                continue;
            }
            if (hetFreq[site2] > duplicateHetsThreshold()) continue;
            LDResult ld = LinkageDisequilibrium.calculateBitLDForHaplotype(20, minorAlleleCnt,genotypeTable, posIndex, site2);
            if (Double.isNaN(ld.r2())) {
                continue;
            }
            highestLD.add(ld);
        }
        List<Position> positionList = new ArrayList<>();
        for (LDResult result : highestLD) {
            positionList.add(genotypeTable.positions().get(result.site2()));
        }
        return PositionListBuilder.getInstance(positionList);
    }

    /**
     * Imputes to the most common genotype weighted by distance
     * @param distGeno Multimap of the closest genotypes and their distance
     * @param useLDSites Number of high LD sites used.
     * @return The imputed genotype
     */
    private byte impute(Multimap<Double, Byte> distGeno, int useLDSites) {
        // useLDSites  is used to scale distance so is similar to DMs original implementation.
        // Seems to have at most a small effect on accuracy.  Could be removed?

        // Create an array to store the weighted counts of each genotype
        double[] weightedCount = new double[256];

        // For each distance to genotype / genotype pair update the weighted counts
        distGeno.entries().forEach(entry -> {
            // +128 is because bytes have values from -128..127 but we want 0..255 for array indexes
            if(entry.getValue()!=UNKNOWN_DIPLOID_ALLELE) weightedCount[entry.getValue() + 128] += 1.0 / (1.0 + useLDSites * entry.getKey());
        });

        // Find the best genotype - the one with the maximum rate
        int bestGeno = 0;
        double bestWeightedCount = weightedCount[0];
        double secondBest = bestWeightedCount;
        for (int i = 1; i < 256; i++) {
            if (weightedCount[i] > bestWeightedCount) {
                secondBest=bestWeightedCount;
                bestWeightedCount = weightedCount[i];
                bestGeno = i;
            }
        }
        if(bestWeightedCount<(minCallBestGenoRatio() *secondBest)) return UNKNOWN_DIPLOID_ALLELE;

        //Return the best genotype.  -128 is for the same reason we added it above.
        return (byte) (bestGeno - 128);
    }

    @Override
    public String getCitation() {
        return "Daniel Money, Kyle Gardner, Heidi Schwaninger, Gan-Yuan Zhong, Sean Myles. (In Review) " +
                " LinkImpute: fast and accurate genotype imputation for non-model organisms";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "LD KNNi Imputation";
    }

    @Override
    public String getToolTipText() {
        return "LD KNNi Imputation";
    }


    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    public static void main(String[] args) {
         GeneratePluginCode.generate(LDKNNiImputationHetV3Plugin.class);
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Number of sites in high LD to use in imputation
     *
     * @return High LD Sites
     */
    public Integer highLDSSites() {
        return highLDSSites.value();
    }

    /**
     * Set High LD Sites. Number of sites in high LD to use
     * in imputation
     *
     * @param value High LD Sites
     *
     * @return this plugin
     */
    public LDKNNiImputationHetV3Plugin highLDSSites(Integer value) {
        highLDSSites = new PluginParameter<>(highLDSSites, value);
        return this;
    }

    /**
     * Number of neighbors to use in imputation
     *
     * @return Number of nearest neighbors
     */
    public Integer knnTaxa() {
        return knnTaxa.value();
    }

    /**
     * Set Number of nearest neighbors. Number of neighbors
     * to use in imputation
     *
     * @param value Number of nearest neighbors
     *
     * @return this plugin
     */
    public LDKNNiImputationHetV3Plugin knnTaxa(Integer value) {
        knnTaxa = new PluginParameter<>(knnTaxa, value);
        return this;
    }

    /**
     * Maximum physical distance between sites to search for
     * LD (-1 for no distance cutoff - unlinked chromosomes
     * will be tested)
     *
     * @return Max distance between site to find LD
     */
    public Integer maxDistance() {
        return maxDistance.value();
    }

    /**
     * Set Max distance between site to find LD. Maximum physical
     * distance between sites to search for LD (-1 for no
     * distance cutoff - unlinked chromosomes will be tested)
     *
     * @param value Max distance between site to find LD
     *
     * @return this plugin
     */
    public LDKNNiImputationHetV3Plugin maxDistance(Integer value) {
        maxDistance = new PluginParameter<>(maxDistance, value);
        return this;
    }

    /**
     * Maximum distance from Nearest Neighbor
     *
     * @return Maximum distance from Nearest Neighbor
     */
    public Double maxDistanceFromNN() {
        return maxDistanceFromNN.value();
    }

    /**
     * Set Maximum distance from Nearest Neighbor. Maximum
     * distance from Nearest Neighbor
     *
     * @param value Maximum distance from Nearest Neighbor
     *
     * @return this plugin
     */
    public LDKNNiImputationHetV3Plugin maxDistanceFromNN(Double value) {
        maxDistanceFromNN = new PluginParameter<>(maxDistanceFromNN, value);
        return this;
    }

    /**
     * Threshold for defining heterozygous sites
     *
     * @return HeterozygousThreshold
     */
    public Double duplicateHetsThreshold() {
        return duplicateHetsThreshold.value();
    }

    /**
     * Set HeterozygousThreshold. Threshold for defining heterozygous
     * sites
     *
     * @param value HeterozygousThreshold
     *
     * @return this plugin
     */
    public LDKNNiImputationHetV3Plugin duplicateHetsThreshold(Double value) {
        duplicateHetsThreshold = new PluginParameter<>(duplicateHetsThreshold, value);
        return this;
    }

    /**
     * Set to Major genotype if no imputation result and MAF
     * is below threshold
     *
     * @return Automatic MajorGenotype if MAF
     */
    public Double automaticMajorMAF() {
        return automaticMajorMAF.value();
    }

    /**
     * Set Automatic MajorGenotype if MAF. Set to Major genotype
     * if no imputation result and MAF is below threshold
     *
     * @param value Automatic MajorGenotype if MAF
     *
     * @return this plugin
     */
    public LDKNNiImputationHetV3Plugin automaticMajorMAF(Double value) {
        automaticMajorMAF = new PluginParameter<>(automaticMajorMAF, value);
        return this;
    }

    /**
     * Minimum coverage for donor genotype and LD calculation
     *
     * @return Minimum coverage for donors and LD
     */
    public Double minCoverageForDonors() {
        return minCoverageForDonors.value();
    }

    /**
     * Set Minimum coverage for donors and LD. Minimum coverage
     * for donor genotype and LD calculation
     *
     * @param value Minimum coverage for donors and LD
     *
     * @return this plugin
     */
    public LDKNNiImputationHetV3Plugin minCoverageForDonors(Double value) {
        minCoverageForDonors = new PluginParameter<>(minCoverageForDonors, value);
        return this;
    }

    /**
     * Minimum ratio between best and second best genotype
     * to make a call
     *
     * @return Minimum support ratio for best genotype
     */
    public Double minCallBestGenoRatio() {
        return minCallBestGenoRatio.value();
    }

    /**
     * Set Minimum support ratio for best genotype. Minimum
     * ratio between best and second best genotype to make
     * a call
     *
     * @param value Minimum support ratio for best genotype
     *
     * @return this plugin
     */
    public LDKNNiImputationHetV3Plugin minCallBestGenoRatio(Double value) {
        minCallBestGenoRatio = new PluginParameter<>(minCallBestGenoRatio, value);
        return this;
    }

    /**
     * Maximum number of cores to be used for processing
     *
     * @return Maximum number of cores for processing
     */
    public Integer maxCores() {
        return maxCores.value();
    }

    /**
     * Set Maximum number of cores for processing. Maximum
     * number of cores to be used for processing
     *
     * @param value Maximum number of cores for processing
     *
     * @return this plugin
     */
    public LDKNNiImputationHetV3Plugin maxCores(Integer value) {
        maxCores = new PluginParameter<>(maxCores, value);
        return this;
    }


    /**
    Alternative to current IBS distance measure
    AA <> AA = 0
    Aa <> Aa = 0 distance (normal IBS distance this is 0.5)
    AA <> aa = 1 distance
     */
    //TODO Terry revisit where this should go
    //TAS-787
    public static double[] dist(byte[] b1, byte[] b2, int min) {
        int distance = 0;
        int count = 0;

        for (int i = 0; i < b1.length; i++) {
            //Make sure hets are represented the same way
            byte p1 = getUnphasedSortedDiploidValue(b1[i]);
            byte p2 = getUnphasedSortedDiploidValue(b2[i]);
            //Ignore the case where either genotype is unknown
            if ((p1 != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) && (p2 != GenotypeTable.UNKNOWN_DIPLOID_ALLELE)) {
                // count counts how many snps we've actually used to scale the
                // distance with since some snps will be unknown
                count++;
                //If the genotypes are unknown then we need to increase the distance
                if (p1 != p2) {
                    //If either genotype is a het we must be either AA <> Aa or aa <> Aa so add one
                    if (GenotypeTableUtils.isHeterozygous(p1) ||
                            GenotypeTableUtils.isHeterozygous(p2)) {
                        distance += 1;
                    //else we must be aa <> AA so add two
                    } else {
                        distance += 2;
                    }
                }
            }
        }

        // If we haven't found enough snps to include in the calculation return NaN
        if (count < min) {
            return new double[]{Double.NaN,count};
        }
        //Else return the scaled distance
        else {
            return new double[]{((double) distance / (double) (2 * count)),count};
        }
    }

    class StatsOnSites {
        private int[][][] totalCnts=new int[3][5][5];
        private int siteCnt=0;
        private int[] classCnt=new int[3];

        private StatsOnSites() {
        }

        synchronized void addStats(double hetFreq, int[][] cnts) {
            int classIndex;
            siteCnt++;
            double hetAcc=(double)cnts[1][1]/(double)(cnts[1][0]+cnts[1][1]+cnts[1][2]);
            if(hetFreq>duplicateHetsThreshold()) {
                classIndex = (hetAcc>0.5)?0:1;
            } else {
                classIndex=2;
            }
            classCnt[classIndex]++;
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 5; j++) {
                    totalCnts[classIndex][i][j]+=cnts[i][j];
                }
            }
        }

        double homozygousAcc(int siteClass) {
            return (double)(totalCnts[siteClass][0][0]+totalCnts[siteClass][2][2])/
                    (double)(totalCnts[siteClass][0][0]+totalCnts[siteClass][0][1]+totalCnts[siteClass][0][2]+
                    totalCnts[siteClass][2][0]+totalCnts[siteClass][2][1]+totalCnts[siteClass][2][2]);
        }

        double recallPowerOfHomozgyous(int siteClass) {
            return (double)(totalCnts[siteClass][0][0]+totalCnts[siteClass][2][2])/
                    (double)(IntStream.of(totalCnts[siteClass][0]).sum()+IntStream.of(totalCnts[siteClass][2]).sum());
        }

        void printToStdOut() {
            for (int siteClass = 0; siteClass < 3; siteClass++) {
                System.out.println("SiteClass"+siteClass+" Count:"+classCnt[siteClass]);
                for (int i = 0; i < totalCnts[0].length; i++) {
                    for (int j = 0; j < totalCnts[0][0].length; j++) {
                        System.out.print(totalCnts[siteClass][i][j]+"\t");
                    }
                    System.out.println();
                }
                System.out.printf("Homozygous Acc:%.4g %n",homozygousAcc(siteClass));
            }
            System.out.println("StatsOnSites siteCnt = " + siteCnt);
        }
    }

    class StatsOnOneSite {

        private int[][] cnts =new int[6][6];
        private final double maf, hetFreq;
        private final boolean isPolymorphic;
        private final int ldSiteCnt;
        Map<Byte,Integer> genotypeToIndexMap=new HashMap<>();
        List<Byte> headers;

        private StatsOnOneSite(double maf, double hetFreq, boolean isPolymorphic, int ldSiteCnt, byte majorAllele, byte minorAllele) {
            this.maf=maf;
            this.hetFreq=hetFreq;
            this.isPolymorphic = isPolymorphic;
            this.ldSiteCnt=ldSiteCnt;
            if(isPolymorphic == false) {
                System.out.println(majorAllele+":"+minorAllele);
            }
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(majorAllele, majorAllele), 0);
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(majorAllele,minorAllele),1);
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(minorAllele,majorAllele),1);
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(minorAllele,minorAllele),2);
            genotypeToIndexMap.put(GenotypeTable.UNKNOWN_DIPLOID_ALLELE,3);
            genotypeToIndexMap.put(NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE,4);
            headers= Bytes.asList(
                    GenotypeTableUtils.getDiploidValue(majorAllele, majorAllele),
                    GenotypeTableUtils.getDiploidValue(majorAllele, minorAllele),
                    GenotypeTableUtils.getDiploidValue(minorAllele, minorAllele),
                    GenotypeTable.UNKNOWN_DIPLOID_ALLELE,
                    NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE,
                    RARE_DIPLOID_ALLELE
            );
        }

        void updateStats(byte origGeno, byte impGeno) {
            int originalIndex= genotypeToIndexMap.getOrDefault(origGeno,5);
            int impIndex= genotypeToIndexMap.getOrDefault(impGeno,5);
            cnts[originalIndex][impIndex]++;
        }

        double minorAcc() {
            return (double)(cnts[2][2])/(double)(cnts[2][0]+cnts[2][1]+cnts[2][2]+0.5);
        }

        double homozygousAcc() {
            return (double)(cnts[0][0]+cnts[2][2])/(double)(cnts[0][0]+cnts[0][1]+cnts[0][2]+cnts[2][0]+cnts[2][1]+cnts[2][2]+0.5);
        }

        double hetAcc() {
            return (double)cnts[1][1]/(double)(cnts[1][0]+cnts[1][1]+cnts[1][2]+0.5);
        }

        double hetFreq() {
            return hetFreq;
        }

        double gapCallRatioOnMajor() {
            return (double)cnts[0][4]/(double)IntStream.of(cnts[0]).sum();
        }

        double gapCallRatioOnMissing() {
            return (double)cnts[3][4]/(double)IntStream.of(cnts[3]).sum();
        }

        double recallPowerOfHomozgyous() {
            return (double)(cnts[0][0]+cnts[2][2])/
                    (double)(IntStream.of(cnts[0]).sum()+IntStream.of(cnts[2]).sum());
        }

        int[][] getCnts() {
            return cnts;
        }

        void printToStdOut() {
            System.out.printf("isPolymorphic:%s MAF:%g LDSites:%d %n", isPolymorphic, maf, ldSiteCnt);
            for (int i = 0; i < cnts[0].length; i++) {
                byte label=headers.get(i);
                System.out.print(NucleotideAlignmentConstants.getNucleotideIUPAC(label)+"\t");
                for (int j = 0; j < cnts[0].length; j++) {
                    System.out.print(cnts[i][j]+"\t");
                }
                System.out.println();
            }
        }
    }
}
