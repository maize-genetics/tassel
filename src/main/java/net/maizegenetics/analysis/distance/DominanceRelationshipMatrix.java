/*
 *  DominanceRelationshipMatrix
 * 
 *  Created on Oct 31, 2015
 */
package net.maizegenetics.analysis.distance;

import java.util.Arrays;
import java.util.Optional;
import java.util.Spliterator;
import static java.util.Spliterator.IMMUTABLE;
import java.util.function.Consumer;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.DistanceMatrixBuilder;
import net.maizegenetics.util.GeneralAnnotationStorage;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Tuple;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Compute Dominance Relationship Matrix for all pairs of taxa. Missing sites
 * are ignored. http://www.genetics.org/content/198/4/1759.abstract
 *
 * @author Terry Casstevens
 */
public class DominanceRelationshipMatrix {

    private static final Logger myLogger = LogManager.getLogger(DominanceRelationshipMatrix.class);

    private static final int DEFAULT_MAX_ALLELES = 6;
    private static final KinshipPlugin.ALGORITHM_VARIATION DEFAULT_ALGORITHM_VARIATION = KinshipPlugin.ALGORITHM_VARIATION.Observed_Allele_Freq;

    private DominanceRelationshipMatrix() {
        // utility
    }

    /**
     * Compute Dominance Relationship Matrix for all pairs of taxa. Missing
     * sites are ignored.
     *
     * @param genotype Genotype Table used to compute dominance relationship
     *
     * @return Dominance Relationship Matrix
     */
    public static DistanceMatrix getInstance(GenotypeTable genotype) {
        return getInstance(genotype, DEFAULT_MAX_ALLELES, DEFAULT_ALGORITHM_VARIATION, null);
    }

    /**
     * Compute Dominance Relationship Matrix for all pairs of taxa. Missing
     * sites are ignored.
     *
     * @param genotype Genotype Table used to compute dominance relationship
     * @param maxAlleles
     * @param variation
     * @param listener progress listener
     *
     * @return Dominance Relationship Matrix
     */
    public static DistanceMatrix getInstance(GenotypeTable genotype, int maxAlleles, KinshipPlugin.ALGORITHM_VARIATION variation, ProgressListener listener) {
        return computeDominanceRelationships(genotype, maxAlleles, variation, listener);
    }

    private static DistanceMatrix computeDominanceRelationships(GenotypeTable genotype, int maxAlleles, KinshipPlugin.ALGORITHM_VARIATION variation, ProgressListener listener) {

        if ((maxAlleles < 2) || (maxAlleles > 6)) {
            throw new IllegalArgumentException("DominanceRelationshipMatrix: computeDominanceRelationships: max alleles must be between 2 and 6 inclusive.");
        }

        if ((variation != KinshipPlugin.ALGORITHM_VARIATION.Observed_Allele_Freq) && (variation != KinshipPlugin.ALGORITHM_VARIATION.Proportion_Heterozygous)) {
            throw new IllegalArgumentException("DominanceRelationshipMatrix: computeDominanceRelationships: variation must be: " + KinshipPlugin.ALGORITHM_VARIATION.Observed_Allele_Freq + " or " + KinshipPlugin.ALGORITHM_VARIATION.Proportion_Heterozygous);
        }

        int numTaxa = genotype.numberOfTaxa();
        long time = System.currentTimeMillis();

        //
        // Sets up parellel stream to divide up sites for processing.
        // Also reduces the distance sums and sum of frequencies into one instance.
        //
        Optional<CountersDistances> optional = stream(genotype, maxAlleles, variation, listener).reduce((CountersDistances t, CountersDistances u) -> {
            t.addAll(u);
            return t;
        });

        if (!optional.isPresent()) {
            return null;
        }
        CountersDistances counters = optional.get();
        double sumpk = counters.mySumOfVariances;
        float[] distances = counters.myDistances;

        //
        // This does the final division of the frequency sum into
        // the distance sums.
        //
        GeneralAnnotationStorage.Builder annotations = GeneralAnnotationStorage.getBuilder();
        annotations.addAnnotation(DistanceMatrixBuilder.MATRIX_TYPE, KinshipPlugin.KINSHIP_METHOD.Dominance_Centered_IBS.toString());
        annotations.addAnnotation(DistanceMatrixBuilder.MATRIX_ALGORITHM_VARIATION, variation.toString());

        DistanceMatrixBuilder builder = DistanceMatrixBuilder.getInstance(genotype.taxa());
        builder.annotation(annotations.build());

        int index = 0;
        for (int t = 0; t < numTaxa; t++) {
            for (int i = 0, n = numTaxa - t; i < n; i++) {
                builder.set(t, t + i, distances[index] / sumpk);
                index++;
            }
        }

        myLogger.info("DominanceRelationshipMatrix: computeDominanceRelationships time: " + (System.currentTimeMillis() - time) / 1000 + " seconds");
        return builder.build();

    }

    private static void fireProgress(int percent, ProgressListener listener) {
        if (listener != null) {
            if (percent > 100) {
                percent = 100;
            }
            listener.progress(percent, null);
        }
    }

    //
    // Each CPU thread (process) creates an instance of this class
    // to acculate terms of the Dominance equation. These are
    // combined with addAll() to result in one instance at the end.
    //
    private static class CountersDistances {

        private double mySumOfVariances = 0.0;
        private final float[] myDistances;
        private final int myNumTaxa;

        public CountersDistances(int numTaxa) {
            myNumTaxa = numTaxa;
            myDistances = new float[myNumTaxa * (myNumTaxa + 1) / 2];
        }

        public void addAll(CountersDistances counters) {
            float[] otherDistances = counters.myDistances;
            for (int t = 0, n = myDistances.length; t < n; t++) {
                myDistances[t] += otherDistances[t];
            }
            mySumOfVariances += counters.mySumOfVariances;
        }

    }

    //
    // This pre-calculates the number of occurances of the major
    // for all possible diploid major values.  Numbers 0 through 7
    // represent A, C, G, T, -, +, N respectively.  First three bits
    // codes the major.  Remaining six bits codes the diploid
    // major values. The stored counts are encodings.  Value 7 (bits 111) means
    // it's not a comparable combination because either major major
    // is unknown or the diploid major value is unknown.
    // Code 1 (bits 001) is heterozygous.
    // Code 2 (bits 010) is homozygous.
    //
    private static final byte[] PRECALCULATED_COUNTS = new byte[512];

    static {
        for (int allele = 0; allele < 8; allele++) {
            for (int a = 0; a < 8; a++) {
                for (int b = 0; b < 8; b++) {
                    int temp = (allele << 6) | (a << 3) | b;
                    if ((allele == 7) || ((a == 7) && (b == 7))) {
                        PRECALCULATED_COUNTS[temp] = 7;
                    } else if (((allele == a) || (allele == b)) && (a != b)) {
                        PRECALCULATED_COUNTS[temp] = 1;
                    } else {
                        PRECALCULATED_COUNTS[temp] = 2;
                    }
                }
            }
        }
    }

    private static final int NUM_CORES_TO_USE = Runtime.getRuntime().availableProcessors() - 1;

    //
    // Used to report progress.  This is not thread-safe but
    // works well enough for this purpose.
    //
    private static int myNumSitesProcessed = 0;

    //
    // Creates stream from DominanceSiteSpliterator and Genotype Table
    //
    private static Stream<CountersDistances> stream(GenotypeTable genotypes, int maxAlleles, KinshipPlugin.ALGORITHM_VARIATION variation, ProgressListener listener) {
        myNumSitesProcessed = 0;
        return StreamSupport.stream(new DominanceSiteSpliterator(genotypes, 0, genotypes.numberOfSites(), maxAlleles, variation, listener), true);
    }

    //
    // Spliterator that splits the sites
    //
    static class DominanceSiteSpliterator implements Spliterator<CountersDistances> {

        private int myCurrentSite;
        private final int myFence;
        private final GenotypeTable myGenotypes;
        private final int myNumTaxa;
        private final int myNumSites;
        private final ProgressListener myProgressListener;
        private int myMinSitesToProcess;
        private final int myMaxAlleles;
        private final KinshipPlugin.ALGORITHM_VARIATION myVariation;

        DominanceSiteSpliterator(GenotypeTable genotypes, int currentIndex, int fence, int maxAlleles, KinshipPlugin.ALGORITHM_VARIATION variation, ProgressListener listener) {
            myGenotypes = genotypes;
            myNumTaxa = myGenotypes.numberOfTaxa();
            myNumSites = myGenotypes.numberOfSites();
            myCurrentSite = currentIndex;
            myFence = fence;
            myMaxAlleles = maxAlleles;
            myVariation = variation;
            myProgressListener = listener;
            myMinSitesToProcess = myNumSites / NUM_CORES_TO_USE;
            if (myMinSitesToProcess == 0) {
                myMinSitesToProcess = myNumSites;
            }
        }

        @Override
        public void forEachRemaining(Consumer<? super CountersDistances> action) {

            int numSitesProcessed = myFence - myCurrentSite;
            CountersDistances result = new CountersDistances(myNumTaxa);
            float[] distances = result.myDistances;
            double[] sumOfVariances = new double[1];

            float[] answer1 = new float[32768];
            float[] answer2 = new float[32768];
            float[] answer3 = new float[32768];

            for (; myCurrentSite < myFence;) {

                //
                // This keeps track of number of sites processed.  The blocks
                // of sites may contain entries for minor major, 2nd minor
                // major, etc.
                int[] realSites = new int[1];

                //
                // Pre-calculates possible terms and gets counts for
                // three blocks for five sites.
                //
                Tuple<short[], float[]> firstBlock = getBlockOfSites(myCurrentSite, sumOfVariances, realSites);
                float[] possibleTerms = firstBlock.y;
                short[] dominanceEffect1 = firstBlock.x;

                Tuple<short[], float[]> secondBlock = getBlockOfSites(myCurrentSite + realSites[0], sumOfVariances, realSites);
                float[] possibleTerms2 = secondBlock.y;
                short[] dominanceEffect2 = secondBlock.x;

                Tuple<short[], float[]> thirdBlock = getBlockOfSites(myCurrentSite + realSites[0], sumOfVariances, realSites);
                float[] possibleTerms3 = thirdBlock.y;
                short[] dominanceEffect3 = thirdBlock.x;

                myCurrentSite += realSites[0];

                //
                // Using possible terms, calculates all possible answers
                // for each site block.
                //
                for (int i = 0; i < 32768; i++) {
                    answer1[i] = possibleTerms[(i & 0x7000) >>> 12] + possibleTerms[((i & 0xE00) >>> 9) | 0x8] + possibleTerms[((i & 0x1C0) >>> 6) | 0x10] + possibleTerms[((i & 0x38) >>> 3) | 0x18] + possibleTerms[(i & 0x7) | 0x20];
                    answer2[i] = possibleTerms2[(i & 0x7000) >>> 12] + possibleTerms2[((i & 0xE00) >>> 9) | 0x8] + possibleTerms2[((i & 0x1C0) >>> 6) | 0x10] + possibleTerms2[((i & 0x38) >>> 3) | 0x18] + possibleTerms2[(i & 0x7) | 0x20];
                    answer3[i] = possibleTerms3[(i & 0x7000) >>> 12] + possibleTerms3[((i & 0xE00) >>> 9) | 0x8] + possibleTerms3[((i & 0x1C0) >>> 6) | 0x10] + possibleTerms3[((i & 0x38) >>> 3) | 0x18] + possibleTerms3[(i & 0x7) | 0x20];
                }

                //
                // Iterates through all pair-wise combinations of taxa adding
                // distance comparisons and site counts.
                //
                int index = 0;
                for (int firstTaxa = 0; firstTaxa < myNumTaxa; firstTaxa++) {
                    //
                    // Can skip inter-loop if all fifteen sites for first
                    // taxon is Unknown diploid major values
                    //
                    if ((dominanceEffect1[firstTaxa] != 0x7FFF) || (dominanceEffect2[firstTaxa] != 0x7FFF) || (dominanceEffect3[firstTaxa] != 0x7FFF)) {
                        for (int secondTaxa = firstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                            //
                            // Combine first taxon's dominance effect with
                            // second taxon's dominance effect to
                            // create index into pre-calculated answers
                            //
                            distances[index] += answer1[dominanceEffect1[firstTaxa] | dominanceEffect1[secondTaxa]] + answer2[dominanceEffect2[firstTaxa] | dominanceEffect2[secondTaxa]] + answer3[dominanceEffect3[firstTaxa] | dominanceEffect3[secondTaxa]];
                            index++;
                        }
                    } else {
                        index += myNumTaxa - firstTaxa;
                    }
                }
            }

            result.mySumOfVariances = sumOfVariances[0];
            action.accept(result);
            myNumSitesProcessed += numSitesProcessed;
            fireProgress((int) ((double) myNumSitesProcessed / (double) myNumSites * 100.0), myProgressListener);
        }

        private static final int NUM_SITES_PER_BLOCK = 5;

        private Tuple<short[], float[]> getBlockOfSites(int currentSite, double[] sumOfVariances, int[] realSites) {

            int currentSiteNum = 0;

            //
            // This hold possible terms for the dominance calculation.
            // First three bits
            // identifies relative site (0, 1, 2, 3, 4).  Remaining three bits
            // whether heterozygous, homozygous, or unknown.
            //
            float[] possibleTerms = new float[40];

            //
            // This holds count of major for each taxa.
            // Each short holds dominance effect (heterozygous, missing,
            // and homozygous) for all NUM_SITES_PER_BLOCK sites
            // at given taxon.  The encodings are stored in three
            // bits each.
            //
            short[] dominanceEffect = new short[myNumTaxa];

            //
            // This initializes the counts to 0x7FFF.  That means
            // diploid allele values for the four sites are Unknown.
            //
            Arrays.fill(dominanceEffect, (short) 0x7FFF);

            while ((currentSiteNum < NUM_SITES_PER_BLOCK) && (currentSite < myFence)) {

                byte[] genotypes = myGenotypes.genotypeAllTaxa(currentSite);
                int[][] alleles = AlleleFreqCache.allelesSortedByFrequencyNucleotide(genotypes);
                int numAlleles = Math.min(alleles[0].length - 1, myMaxAlleles - 1);

                if (numAlleles + currentSiteNum <= NUM_SITES_PER_BLOCK) {

                    //
                    // Calculates total number of haploid alleles that
                    // are not missing.
                    //
                    int totalAlleleCount = 0;
                    for (int i = 0; i < alleles[1].length; i++) {
                        totalAlleleCount += alleles[1][i];
                    }

                    for (int a = 0; a < numAlleles; a++) {

                        byte allele = (byte) alleles[0][a];
                        float standardizedTerm = 0.0f;
                        if (myVariation == KinshipPlugin.ALGORITHM_VARIATION.Observed_Allele_Freq) {
                            float alleleFreq = (float) alleles[1][a] / (float) totalAlleleCount;
                            standardizedTerm = 2.0f * alleleFreq * (1.0f - alleleFreq);
                        } else if (myVariation == KinshipPlugin.ALGORITHM_VARIATION.Proportion_Heterozygous) {
                            standardizedTerm = (float) AlleleFreqCache.proportionHeterozygous(genotypes);
                        }
                        sumOfVariances[0] += standardizedTerm * (1.0 - standardizedTerm);

                        //
                        // Temporarily stores component terms of equation
                        //
                        float[] term = new float[2];

                        //
                        // If allele is Unknown, the entire
                        // site is skipped.
                        //
                        if (allele != GenotypeTable.UNKNOWN_ALLELE) {

                            term[0] = 1.0f - standardizedTerm;
                            term[1] = -standardizedTerm;

                            //
                            // Pre-calculates all possible terms of the summation
                            // for this current (pseudo-) site.
                            // Counts (0,0; 0,1; 0,2; 1,1; 1,2; 2,2)
                            //
                            int siteNumIncrement = currentSiteNum * 8;
                            possibleTerms[siteNumIncrement + 1] = term[0] * term[0];
                            possibleTerms[siteNumIncrement + 3] = term[0] * term[1];
                            possibleTerms[siteNumIncrement + 2] = term[1] * term[1];

                            //
                            // Records allele counts for current site in
                            // three bits.
                            //
                            int temp = (allele & 0x7) << 6;
                            int shift = (NUM_SITES_PER_BLOCK - currentSiteNum - 1) * 3;
                            int mask = ~(0x7 << shift) & 0x7FFF;
                            for (int i = 0; i < myNumTaxa; i++) {
                                dominanceEffect[i] = (short) (dominanceEffect[i] & (mask | PRECALCULATED_COUNTS[temp | ((genotypes[i] & 0x70) >>> 1) | (genotypes[i] & 0x7)] << shift));
                            }

                        }

                        currentSiteNum++;
                    }
                } else {
                    return new Tuple<>(dominanceEffect, possibleTerms);
                }

                currentSite++;
                realSites[0]++;
            }

            return new Tuple<>(dominanceEffect, possibleTerms);

        }

        @Override
        public boolean tryAdvance(Consumer<? super CountersDistances> action) {
            if (myCurrentSite < myFence) {
                forEachRemaining(action);
                return true;
            } else {
                return false;
            }
        }

        @Override
        /**
         * Splits sites
         */
        public Spliterator<CountersDistances> trySplit() {
            int lo = myCurrentSite;
            int mid = lo + myMinSitesToProcess;
            if (mid < myFence) {
                myCurrentSite = mid;
                return new DominanceSiteSpliterator(myGenotypes, lo, mid, myMaxAlleles, myVariation, myProgressListener);
            } else {
                return null;
            }
        }

        @Override
        public long estimateSize() {
            return (long) (myFence - myCurrentSite);
        }

        @Override
        public int characteristics() {
            return IMMUTABLE;
        }
    }

}
