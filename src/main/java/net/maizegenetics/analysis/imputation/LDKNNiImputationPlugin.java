package net.maizegenetics.analysis.imputation;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.MinMaxPriorityQueue;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import net.maizegenetics.analysis.popgen.LDResult;
import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.stats.statistics.FisherExact;
import net.maizegenetics.util.BitSet;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.IntStream;

import static net.maizegenetics.dna.snp.GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
import static net.maizegenetics.dna.snp.GenotypeTableUtils.getDiploidValue;
import static net.maizegenetics.dna.snp.GenotypeTableUtils.getUnphasedSortedDiploidValue;
import net.maizegenetics.plugindef.GeneratedGuiBoilerplate;

/**
 * This imputation algorithm uses LD to identify good predictors for each SNP,
 * and then uses the high LD SNPs to identify K- Nearest Neighbors.
 * The genotype is called with a weighted mode of the KNNs.
 *
 * @author Daniel Money (Developer and developed the algorithm)
 * @author Ed Buckler (assisted in conversion to TASSEL)
 */
public class LDKNNiImputationPlugin extends AbstractPlugin {

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

    private PluginParameter<Integer> maxDistance = new PluginParameter.Builder<>("maxLDDistance", 10_000_000, Integer.class)
            .guiName("Max distance between site to find LD")
            .description("Maximum physical distance between sites to search for LD (-1 for no distance cutoff - unlinked chromosomes will be tested)")
            .build();


    private static final Logger myLogger = LogManager.getLogger(LDKNNiImputationPlugin.class);

    public LDKNNiImputationPlugin() {
        super(null, false);
    }

    public LDKNNiImputationPlugin(Frame parentFrame, boolean isInteractive) {
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

        // Create a multimap of the SNPs in highest LD with each SNP
/*Debatable on what to calc*/ //        Multimap<Position, Position> highLDMap = getHighLDMap(GenotypeTableBuilder.getHomozygousInstance(genotypeTable), highLDSSites());
        Multimap<Position, Position> highLDMap = getHighLDMap(genotypeTable, highLDSSites());
        System.out.println("LD calculated");

        GenotypeTableBuilder incSiteBuilder = GenotypeTableBuilder.getSiteIncremental(genotypeTable.taxa());
        //Start imputing site by site
        long time=System.nanoTime();
        LongAdder sites1Kdone=new LongAdder();
        IntStream.range(0, genotypeTable.numberOfSites()).parallel().forEach(posIndex ->
        {
            Position position = genotypeTable.positions().get(posIndex);
            PositionList positionList = PositionListBuilder.getInstance(new ArrayList<>(highLDMap.get(position)));
            byte[] currGenos = genotypeTable.genotypeAllTaxa(posIndex);
            byte[] newGenos = new byte[currGenos.length];
            //set monomorphic sites to the major allele
            if (!genotypeTable.isPolymorphic(posIndex)) {
                byte monomorphicGenotype = getDiploidValue(genotypeTable.majorAllele(posIndex), genotypeTable.majorAllele(posIndex));
                for (int i = 0; i < newGenos.length; i++) {
                    newGenos[i] = (currGenos[i] == UNKNOWN_DIPLOID_ALLELE) ? monomorphicGenotype : currGenos[i];
                }
            } else {
                //create a site specific in memory filtered alignment to use
                final GenotypeTable ldGenoTable = GenotypeTableBuilder.getGenotypeCopyInstance(FilterGenotypeTable.getInstance(genotypeTable, positionList));
                final int numberSites = ldGenoTable.numberOfSites();
                double[] taxaCoverage = IntStream.range(0, ldGenoTable.numberOfTaxa()).sequential()  //used determine when insufficient overlap is likely
                        .mapToDouble(t -> (double) ldGenoTable.totalNonMissingForTaxon(t) / (double) numberSites)
                        .toArray();
                final byte[][] ldGenotypes = new byte[ldGenoTable.numberOfTaxa()][];
                for (int t = 0; t < ldGenotypes.length; t++) {
                    ldGenotypes[t] = ldGenoTable.genotypeAllSites(t);
                }
                for (int taxon = 0; taxon < currGenos.length; taxon++) {
                    newGenos[taxon] = currGenos[taxon];
                    if (currGenos[taxon] == UNKNOWN_DIPLOID_ALLELE) {  //starting imputing
                        Multimap<Double, Byte> closeGenotypes = getClosestNonMissingTaxa(taxon, currGenos, ldGenotypes,
                                taxaCoverage, knnTaxa());
                        if (closeGenotypes.isEmpty()) {  //this is empty when two few high LD sites are shared between the target taxon and all others
                            newGenos[taxon] = UNKNOWN_DIPLOID_ALLELE;
                        } else {
                            newGenos[taxon] = impute(closeGenotypes, highLDSSites());  //set to weight mode genotype
                        }
                    }

                }//);
            }
            incSiteBuilder.addSite(position, newGenos);
            if ((posIndex + 1) % 100 == 0) {
                sites1Kdone.add(100);
                fireProgress(33 + ((int) (66 * sites1Kdone.longValue()) / genotypeTable.numberOfSites()));
                System.out.println(sites1Kdone.longValue() + ":" + ((System.nanoTime() - time) / 1_000_000) / sites1Kdone.longValue());
            }
        });
        GenotypeTable impGenotypeTable = incSiteBuilder.build();

        return new DataSet(new Datum(genoDatum.getName()+"_KNNimp",impGenotypeTable,"Imputed genotypes by KNN imputation"),this);
    }

    /**
     * Create a multimap of distances between the target taxon with the closest observed genotypes.  Distance is calculated
     * between the target taxon and all other taxa over the high LD sites.
     * @param inputTaxonIdx Index of the taxon being imputed
     * @param targetGenotypes Genotypes of every taxon at the site being imputed
     * @param ldGenotypes Genotypes of every taxon at the high LD sites, indexed by taxon
     * @param inputCoverage Fraction of the high LD sites scored in each taxon
     * @param numberOfTaxa number of genotypes to retain
     * @return Map of distance to genotype call for the target position
     */
    private static Multimap<Double, Byte> getClosestNonMissingTaxa(int inputTaxonIdx, byte[] targetGenotypes, byte[][] ldGenotypes,
                                                                   double[] inputCoverage, int numberOfTaxa) {
        final byte[] inputTaxonGenotypes = ldGenotypes[inputTaxonIdx];
        final int numberOfLDSites = inputTaxonGenotypes.length;

        //The numberOfTaxa smallest (distance, genotype) pairs, in ascending order
        final double[] topDistances = new double[numberOfTaxa];
        final byte[] topGenotypes = new byte[numberOfTaxa];
        int size = 0;
        for (int closeTaxonIdx = 0; closeTaxonIdx < targetGenotypes.length; closeTaxonIdx++) {
            if (closeTaxonIdx == inputTaxonIdx) continue;  //do not test itself
            //skip tests with too little overlap; written as !(x > 10) so a NaN coverage is skipped too
            if (!(inputCoverage[closeTaxonIdx] * inputCoverage[inputTaxonIdx] * (double) numberOfLDSites > 10)) continue;
            final byte genotype = targetGenotypes[closeTaxonIdx];
            if (genotype == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) continue;  //ignore taxa with the genotype not scored
            final double distance = dist(inputTaxonGenotypes, ldGenotypes[closeTaxonIdx], 10)[0];
            if (Double.isNaN(distance)) continue;  //skip is too few sites (<10 results in NaN)
            if (size == numberOfTaxa
                    && compareDistanceGenotype(distance, genotype, topDistances[size - 1], topGenotypes[size - 1]) >= 0) {
                continue;
            }
            int i = (size < numberOfTaxa) ? size++ : size - 1;
            while (i > 0 && compareDistanceGenotype(distance, genotype, topDistances[i - 1], topGenotypes[i - 1]) < 0) {
                topDistances[i] = topDistances[i - 1];
                topGenotypes[i] = topGenotypes[i - 1];
                i--;
            }
            topDistances[i] = distance;
            topGenotypes[i] = genotype;
        }
        final Multimap<Double, Byte> distGenoMap = ArrayListMultimap.create();
        for (int i = 0; i < size; i++) {
            distGenoMap.put(topDistances[i], topGenotypes[i]);
        }
        return distGenoMap;
    }

    private static int compareDistanceGenotype(double distance1, byte genotype1, double distance2, byte genotype2) {
        int result = Double.compare(distance1, distance2);
        return (result != 0) ? result : Byte.compare(genotype1, genotype2);
    }

    private Multimap<Position, Position> getHighLDMap(GenotypeTable genotypeTable, int numberOfSNPs) {
        final int numberOfSites = genotypeTable.numberOfSites();
        // ArrayListMultimap is not thread-safe, so the parallel loop fills one slot per site
        // and the multimap is assembled afterwards.
        @SuppressWarnings("unchecked")
        final List<Position>[] highLDSites = new List[numberOfSites];
        final BitSet[] majorPresence = new BitSet[numberOfSites];
        final BitSet[] minorPresence = new BitSet[numberOfSites];
        for (int site = 0; site < numberOfSites; site++) {
            majorPresence[site] = genotypeTable.allelePresenceForAllTaxa(site, WHICH_ALLELE.Major);
            minorPresence[site] = genotypeTable.allelePresenceForAllTaxa(site, WHICH_ALLELE.Minor);
        }
        final FisherExact fisherExact = FisherExact.getInstance((2 * genotypeTable.numberOfTaxa()) + 10);
        LongAdder sites1Kdone=new LongAdder();
        IntStream.range(0, genotypeTable.numberOfSites()).parallel()
                .forEach(posIndex ->
                {
                    MinMaxPriorityQueue<LDResult> highestLD = MinMaxPriorityQueue.orderedBy(LDResult.byR2Ordering.reverse())
                            .maximumSize(numberOfSNPs).create();
                    for (int site2 = 0; site2 < numberOfSites; site2++) {
                        if (posIndex == site2) {
                            continue;
                        }
                        if(maxDistance()>-1 && Math.abs(genotypeTable.chromosomalPosition(posIndex)-genotypeTable.chromosomalPosition(site2))>maxDistance()) {
                            continue;
                        }
                        LDResult ld = LinkageDisequilibrium.getLDForSitePair(majorPresence[posIndex], minorPresence[posIndex],
                                majorPresence[site2], minorPresence[site2], 2, 20, -1.0f, fisherExact, posIndex, site2);
                        if (Double.isNaN(ld.r2())) {
                            continue;
                        }
                        highestLD.add(ld);
                    }
                    List<Position> positionList = new ArrayList<>();
                    for (LDResult result : highestLD) {
                        positionList.add(genotypeTable.positions().get(result.site2()));
                    }
                    highLDSites[posIndex] = positionList;
                    if((posIndex+1)%1000==0) {
                        sites1Kdone.add(1000);
                        fireProgress((int)(33*sites1Kdone.longValue())/numberOfSites);
                        System.out.println(sites1Kdone.longValue());
                    }
                });

        Multimap<Position, Position> highLDMap = ArrayListMultimap.create();
        for (int posIndex = 0; posIndex < numberOfSites; posIndex++) {
            highLDMap.putAll(genotypeTable.positions().get(posIndex), highLDSites[posIndex]);
        }
        return highLDMap;
    }

    /**
     * Imputes to the most common genotype weighted by distance
     * @param distGeno Multimap of the closest genotypes and their distance
     * @param useLDSites Number of high LD sites used.
     * @return The imputed genotype
     */
    static byte impute(Multimap<Double, Byte> distGeno, int useLDSites) {
        // useLDSites  is used to scale distance so is similar to DMs original implementation.
        // Seems to have at most a small effect on accuracy.  Could be removed?

        // Create an array to store the weighted counts of each genotype
        double[] weightedCount = new double[256];

        // For each distance to genotype / genotype pair update the weighted counts
        distGeno.entries().forEach(entry -> {
            // +128 is because bytes have values from -128..127 but we want 0..255 for array indexes
            weightedCount[entry.getValue() + 128] += 1.0 / (1.0 + useLDSites * entry.getKey());
        });

        // Find the best genotype - the one with the maximum rate
        int bestGeno = 0;
        double bestWeightedCount = weightedCount[0];
        for (int i = 1; i < 256; i++) {
            if (weightedCount[i] > bestWeightedCount) {
                bestWeightedCount = weightedCount[i];
                bestGeno = i;
            }
        }

        //Return the best genotype.  -128 is for the same reason we added it above.
        return (byte) (bestGeno - 128);
    }

    /*
    Alternative approach to determine the best genotype
    This is an allele based approaches.  Provides a weighting of each allele, and Heterozygous is called currently
    if the second allele has the count more than 1/2 the most common.

     Provides roughly the same results as the genotype approach.  An ensemble is best but lower call rate.
     Both approaches have value.
     */
//    public static byte impute(Multimap<Double, Byte> distGeno, int useLDSites) {
//        // Create an array to store the weighted counts
//        TByteDoubleMap alleleWeight=new TByteDoubleHashMap(16);
//
//        // For each distance to genotype / genotype pair update the weighted counts
//        distGeno.entries().forEach(entry -> {
//            byte[] alleles= GenotypeTableUtils.getDiploidValues(entry.getValue());
//            double weight= 1.0 / (1.0 + useLDSites * entry.getKey());
//            alleleWeight.adjustOrPutValue(alleles[0],weight,weight);
//            alleleWeight.adjustOrPutValue(alleles[1],weight,weight);
//        });
//
//        double[] weights= alleleWeight.values();
//        Arrays.sort(weights);
//
//        if(weights.length==1) {
//            return GenotypeTableUtils.getUnphasedDiploidValue(alleleWeight.keys()[0], alleleWeight.keys()[0]);
//        } else if(weights[weights.length-1]>(2*weights[weights.length-2])) {
//            alleleWeight.retainEntries((a,w) -> w>=weights[weights.length-1]);
//            byte[] alleles=alleleWeight.keys();
//            return GenotypeTableUtils.getUnphasedDiploidValue(alleles[0], alleles[0]);
//        } else {
//            alleleWeight.retainEntries((a,w) -> w>=weights[weights.length-2]);
//            byte[] alleles=alleleWeight.keys();
//            return GenotypeTableUtils.getUnphasedDiploidValue(alleles[0], alleles[1]);
//        }
//    }

    @Override
    public String getCitation() {
        return "Daniel Money, Kyle Gardner, Heidi Schwaninger, Gan-Yuan Zhong, Sean Myles. (In Review) " +
                " LinkImpute: fast and accurate genotype imputation for non-model organisms";
    }

    @Override
    @GeneratedGuiBoilerplate
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    @GeneratedGuiBoilerplate
    public String getButtonName() {
        return "LD KNNi Imputation";
    }

    @Override
    @GeneratedGuiBoilerplate
    public String getToolTipText() {
        return "LD KNNi Imputation";
    }


    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    public static void main(String[] args) {
         GeneratePluginCode.generate(LDKNNiImputationPlugin.class);
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    @GeneratedGuiBoilerplate
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Maximum number of sites in high LD to use in imputation
     *
     * @return High LD Sites
     */
    @GeneratedGuiBoilerplate
    public Integer highLDSSites() {
        return highLDSSites.value();
    }

    /**
     * Set High LD Sites. Maximum number of sites in high
     * LD to use in imputation
     *
     * @param value High LD Sites
     *
     * @return this plugin
     */
    @GeneratedGuiBoilerplate
    public LDKNNiImputationPlugin highLDSSites(Integer value) {
        highLDSSites = new PluginParameter<>(highLDSSites, value);
        return this;
    }

    /**
     * Maximum number of neighbours to use in imputation
     *
     * @return Number of nearest neighbors
     */
    @GeneratedGuiBoilerplate
    public Integer knnTaxa() {
        return knnTaxa.value();
    }

    /**
     * Set Number of nearest neighbors. Maximum number of
     * neighbours to use in imputation
     *
     * @param value Number of nearest neighbors
     *
     * @return this plugin
     */
    @GeneratedGuiBoilerplate
    public LDKNNiImputationPlugin knnTaxa(Integer value) {
        knnTaxa = new PluginParameter<>(knnTaxa, value);
        return this;
    }

    /**
     * Maximum physical distance between sites to look for
     * LD (-1 for no distance cutoff - unlinked chromosomes
     * will be tested)
     *
     * @return Max distance between site to find LD
     */
    @GeneratedGuiBoilerplate
    public Integer maxDistance() {
        return maxDistance.value();
    }

    /**
     * Set Max distance between site to find LD. Maximum physical
     * distance between sites to look for LD (-1 for no distance
     * cutoff - unlinked chromosomes will be tested)
     *
     * @param value Max distance between site to find LD
     *
     * @return this plugin
     */
    @GeneratedGuiBoilerplate
    public LDKNNiImputationPlugin maxDistance(Integer value) {
        maxDistance = new PluginParameter<>(maxDistance, value);
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
}
