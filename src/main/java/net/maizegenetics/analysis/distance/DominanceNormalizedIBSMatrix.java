/*
 *  DominanceNormalizedIBSMatrix
 * 
 *  Created on Nov 26, 2015
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
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.DistanceMatrixBuilder;
import net.maizegenetics.util.GeneralAnnotationStorage;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Tuple;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class DominanceNormalizedIBSMatrix {

    private static final Logger myLogger = LogManager.getLogger(DominanceNormalizedIBSMatrix.class);

    private DominanceNormalizedIBSMatrix() {
        // utility
    }

    /**
     * http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4375616/pdf/main.pdf p. 378
     *
     * @param genotype Genotype Table used to compute matrix
     *
     * @return Dominance Normalized IBS Matrix
     */
    public static DistanceMatrix getInstance(GenotypeTable genotype) {
        return getInstance(genotype, null);
    }

    /**
     * Same as other getInstance() but reports progress.
     *
     * @param genotype Genotype Table used to compute matrix
     * @param listener Progress listener
     *
     * @return Dominance Normalized IBS Matrix
     */
    public static DistanceMatrix getInstance(GenotypeTable genotype, ProgressListener listener) {
        return computeDominanceNormalizedIBSMatrix(genotype, listener);
    }

    private static DistanceMatrix computeDominanceNormalizedIBSMatrix(GenotypeTable genotype, ProgressListener listener) {

        int numTaxa = genotype.numberOfTaxa();
        long time = System.currentTimeMillis();

        //
        // Sets up parellel stream to divide up sites for processing.
        // Also reduces the distance sums and site counters into one instance.
        //
        Optional<CountersDistances> optional = stream(genotype, listener).reduce((CountersDistances t, CountersDistances u) -> {
            t.addAll(u);
            return t;
        });

        if (!optional.isPresent()) {
            return null;
        }
        CountersDistances counters = optional.get();
        int[] counts = counters.myCounters;
        float[] distances = counters.myDistances;

        //
        // This does the final division of the site counts into
        // the distance sums.
        //
        GeneralAnnotationStorage.Builder annotations = GeneralAnnotationStorage.getBuilder();
        annotations.addAnnotation(DistanceMatrixBuilder.MATRIX_TYPE, KinshipPlugin.KINSHIP_METHOD.Dominance_Normalized_IBS.toString());
        DistanceMatrixBuilder builder = DistanceMatrixBuilder.getInstance(genotype.taxa());
        builder.annotation(annotations.build());
        int index = 0;
        for (int t = 0; t < numTaxa; t++) {
            for (int i = 0, n = numTaxa - t; i < n; i++) {
                builder.set(t, t + i, distances[index] / (double) counts[index]);
                builder.setCount(t, i, counts[index]);
                index++;
            }
        }

        myLogger.info("DominanceNormalizedIBSMatrix: computeDominanceNormalizedIBSMatrix time: " + (System.currentTimeMillis() - time) / 1000 + " seconds");
        return builder.build();

    }

    protected static void fireProgress(int percent, ProgressListener listener) {
        if (listener != null) {
            if (percent > 100) {
                percent = 100;
            }
            listener.progress(percent, null);
        }
    }

    //
    // Each CPU thread (process) creates an instance of this class
    // to acculate terms of the equation and the number of
    // sites involved for each pair-wise calculation.  These are
    // combined with addAll() to result in one instance at the end.
    //
    private static class CountersDistances {

        private final int[] myCounters;
        private final float[] myDistances;
        private final int myNumTaxa;

        public CountersDistances(int numTaxa) {
            myNumTaxa = numTaxa;
            myCounters = new int[myNumTaxa * (myNumTaxa + 1) / 2];
            myDistances = new float[myNumTaxa * (myNumTaxa + 1) / 2];
        }

        public void addAll(CountersDistances counters) {
            float[] otherDistances = counters.myDistances;
            for (int t = 0, n = myCounters.length; t < n; t++) {
                myDistances[t] += otherDistances[t];
            }
            otherDistances = null;
            int[] otherCounters = counters.myCounters;
            for (int t = 0, n = myCounters.length; t < n; t++) {
                myCounters[t] += otherCounters[t];
            }
        }

    }

    //
    // This pre-calculates the state of the major / minor allele
    // for all possible diploid allele values.  Numbers 0 through 7
    // represent A, C, G, T, -, +, N respectively.  First three bits
    // codes the major allele.  Remaining six bits codes the diploid
    // allele values. The stored counts are encodings.  Value 7 (bits 111) means
    // it's not a comparable combination because either major allele
    // is unknown or the diploid allele value is unknown.
    // Code 2 (bits 010) is heterzygous major / minor.
    // Code 4 (bits 100) is homozygous major or homozygous minor.
    //
    private static final byte[] PRECALCULATED_COUNTS = new byte[512];

    static {
        for (int major = 0; major < 8; major++) {
            for (int a = 0; a < 8; a++) {
                for (int b = 0; b < 8; b++) {
                    int temp = (major << 6) | (a << 3) | b;
                    if ((major == 7) || ((a == 7) && (b == 7))) {
                        PRECALCULATED_COUNTS[temp] = 7;
                    } else if (a == major) {
                        if (b == major) {
                            PRECALCULATED_COUNTS[temp] = 4;
                        } else {
                            PRECALCULATED_COUNTS[temp] = 2;
                        }
                    } else if (b == major) {
                        PRECALCULATED_COUNTS[temp] = 2;
                    } else {
                        PRECALCULATED_COUNTS[temp] = 4;
                    }
                }
            }
        }
    }

    //
    // This pre-calculates the number of sites involved in a pair-wise
    // comparison.  Counts are the number of sites involved in the
    // calculation (up to 5 sites).
    // Count value of 7 is coded when diploid allele value is
    // GenotypeTable.UNKNOWN_DIPLOID_ALLELE.  Any pair-wise comparison when
    // either taxa has GenotypeTable.UNKNOWN_DIPLOID_ALLELE at a given site,
    // is not involved in the calulation. The index of this array represents
    // every bitwise OR combination of genotype class (2, 4) and UNKNOWN (7)
    // for five consecutive sites.  Each three bits encodes two genotype classes.
    // Those three bits times five sites equals 32768 combinations.
    // Code 010 - both genotype class heterozygous
    // Code 110 - one heterozygous, one homozygous
    // Code 100 - both genotype class homozygous
    //
    private static final byte[] INCREMENT = new byte[32768];

    static {
        for (int a = 1; a < 8; a++) {
            int temp = a << 12;
            for (int b = 1; b < 8; b++) {
                int temp2 = b << 9;
                for (int c = 1; c < 8; c++) {
                    int temp3 = c << 6;
                    for (int d = 1; d < 8; d++) {
                        int temp4 = d << 3;
                        for (int e = 1; e < 8; e++) {
                            int incrementIndex = temp | temp2 | temp3 | temp4 | e;
                            if (a != 7) {
                                INCREMENT[incrementIndex]++;
                            }
                            if (b != 7) {
                                INCREMENT[incrementIndex]++;
                            }
                            if (c != 7) {
                                INCREMENT[incrementIndex]++;
                            }
                            if (d != 7) {
                                INCREMENT[incrementIndex]++;
                            }
                            if (e != 7) {
                                INCREMENT[incrementIndex]++;
                            }
                        }
                    }
                }
            }
        }
    }

    private static final int NUM_CORES_TO_USE = TasselPrefs.getMaxThreads();

    //
    // Used to report progress.  This is not thread-safe but
    // works well enough for this purpose.
    //
    private static int myNumSitesProcessed = 0;

    //
    // Creates stream from DominanceNormalizedIBSSiteSpliterator and Genotype Table
    //
    private static Stream<CountersDistances> stream(GenotypeTable genotypes, ProgressListener listener) {
        myNumSitesProcessed = 0;
        return StreamSupport.stream(new DominanceNormalizedIBSSiteSpliterator(genotypes, 0, genotypes.numberOfSites(), listener), true);
    }

    //
    // Spliterator that splits the sites into halves each time for
    // processing.
    //
    static class DominanceNormalizedIBSSiteSpliterator implements Spliterator<CountersDistances> {

        private int myCurrentSite;
        private final int myFence;
        private final GenotypeTable myGenotypes;
        private final int myNumTaxa;
        private final int myNumSites;
        private final ProgressListener myProgressListener;
        private final int myMinSitesToProcess;
        private final int myNumSitesPerBlockForProgressReporting;

        DominanceNormalizedIBSSiteSpliterator(GenotypeTable genotypes, int currentIndex, int fence, ProgressListener listener) {
            myGenotypes = genotypes;
            myNumTaxa = myGenotypes.numberOfTaxa();
            myNumSites = myGenotypes.numberOfSites();
            myCurrentSite = currentIndex;
            myFence = fence;
            myProgressListener = listener;
            myMinSitesToProcess = Math.max(myNumSites / NUM_CORES_TO_USE, 1000);
            myNumSitesPerBlockForProgressReporting = (myFence - myCurrentSite) / 10;
        }

        @Override
        public void forEachRemaining(Consumer<? super CountersDistances> action) {

            CountersDistances result = new CountersDistances(myNumTaxa);
            int[] counts = result.myCounters;
            float[] distances = result.myDistances;;

            float[] answer1 = new float[32768];
            float[] answer2 = new float[32768];
            float[] answer3 = new float[32768];

            for (; myCurrentSite < myFence;) {

                int currentBlockFence = Math.min(myCurrentSite + myNumSitesPerBlockForProgressReporting, myFence);

                int numSitesProcessed = currentBlockFence - myCurrentSite;

                for (; myCurrentSite < currentBlockFence;) {

                    int[] numSites = new int[1];

                    //
                    // Pre-calculates possible terms and gets counts for
                    // three blocks for five sites.
                    //
                    Tuple<short[], float[]> firstBlock = getBlockOfSites(myCurrentSite, numSites);
                    float[] possibleTerms = firstBlock.y;
                    short[] genotypeClass1 = firstBlock.x;

                    Tuple<short[], float[]> secondBlock = getBlockOfSites(myCurrentSite + numSites[0], numSites);
                    float[] possibleTerms2 = secondBlock.y;
                    short[] genotypeClass2 = secondBlock.x;

                    Tuple<short[], float[]> thirdBlock = getBlockOfSites(myCurrentSite + numSites[0], numSites);
                    float[] possibleTerms3 = thirdBlock.y;
                    short[] genotypeClass3 = thirdBlock.x;

                    myCurrentSite += numSites[0];

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
                        // taxon is Unknown diploid allele values
                        //
                        if ((genotypeClass1[firstTaxa] != 0x7FFF) || (genotypeClass2[firstTaxa] != 0x7FFF) || (genotypeClass3[firstTaxa] != 0x7FFF)) {
                            for (int secondTaxa = firstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                                //
                                // Combine first taxon's major allele counts with
                                // second taxon's major allele counts to
                                // create index into pre-calculated answers
                                // and site counts.
                                //
                                distances[index] += answer1[genotypeClass1[firstTaxa] | genotypeClass1[secondTaxa]] + answer2[genotypeClass2[firstTaxa] | genotypeClass2[secondTaxa]] + answer3[genotypeClass3[firstTaxa] | genotypeClass3[secondTaxa]];
                                counts[index] += INCREMENT[genotypeClass1[firstTaxa] | genotypeClass1[secondTaxa]] + INCREMENT[genotypeClass2[firstTaxa] | genotypeClass2[secondTaxa]] + INCREMENT[genotypeClass3[firstTaxa] | genotypeClass3[secondTaxa]];
                                index++;
                            }
                        } else {
                            index += myNumTaxa - firstTaxa;
                        }
                    }
                }

                myNumSitesProcessed += numSitesProcessed;
                fireProgress((int) ((double) myNumSitesProcessed / (double) myNumSites * 100.0), myProgressListener);

            }

            action.accept(result);
        }

        private static final int NUM_SITES_PER_BLOCK = 5;

        private Tuple<short[], float[]> getBlockOfSites(int currentSite, int[] numSites) {

            int currentSiteNum = 0;

            //
            // This hold possible terms for the summation given
            // site's major allele frequency.  First two bits
            // identifies relative site (0, 1, 2, 3, 4).  Remaining three bits
            // the major allele counts encoding.
            //
            float[] possibleTerms = new float[40];

            //
            // This holds genotype class for each taxa.
            // Each short holds class (4, 2, 7) for all sites
            // at given taxon.  The classes are stored in four bits each.
            // This leaves the two higher bits for each empty for shifting.
            //
            short[] genotypeClasses = new short[myNumTaxa];

            //
            // This initializes the counts to 0x7FFF.  That means
            // diploid allele values for the sites are Unknown.
            //
            Arrays.fill(genotypeClasses, (short) 0x7FFF);

            while ((currentSiteNum < NUM_SITES_PER_BLOCK) && (currentSite < myFence)) {

                byte[] genotypes = myGenotypes.genotypeAllTaxa(currentSite);
                int[][] alleleCounts = AlleleFreqCache.allelesSortedByFrequencyNucleotide(genotypes);
                byte major = AlleleFreqCache.majorAllele(alleleCounts);
                float minorFreq = 1.0f - (float) AlleleFreqCache.majorAlleleFrequency(alleleCounts);
                float majorFreqSqrTimes2 = minorFreq * minorFreq * 2.0f;
                float denominatorTerm = minorFreq * minorFreq * 4.0f * (1.0f - minorFreq) * (1.0f - minorFreq);

                //
                // If major allele is Unknown or major allele frequency
                // equals 1.0 (resulting in denominator 0.0), the entire
                // site is skipped.
                //
                if ((major != GenotypeTable.UNKNOWN_ALLELE) && (denominatorTerm != 0.0)) {

                    float hetTerm = 1.0f - majorFreqSqrTimes2;
                    float homoTerm = 0.0f - majorFreqSqrTimes2;

                    //
                    // Pre-calculates all possible terms of the summation
                    // for this current site.  Genotype Classes (het:het; homo:homo; het:homo)
                    //
                    int siteNumIncrement = currentSiteNum * 8;
                    possibleTerms[siteNumIncrement + 2] = hetTerm * hetTerm / denominatorTerm;
                    possibleTerms[siteNumIncrement + 4] = homoTerm * homoTerm / denominatorTerm;
                    possibleTerms[siteNumIncrement + 6] = hetTerm * homoTerm / denominatorTerm;

                    //
                    // Records genotype class for current site in
                    // three bits.
                    //
                    int temp = (major & 0x7) << 6;
                    int shift = (NUM_SITES_PER_BLOCK - currentSiteNum - 1) * 3;
                    int mask = ~(0x7 << shift) & 0x7FFF;
                    for (int i = 0; i < myNumTaxa; i++) {
                        genotypeClasses[i] = (short) (genotypeClasses[i] & (mask | PRECALCULATED_COUNTS[temp | ((genotypes[i] & 0x70) >>> 1) | (genotypes[i] & 0x7)] << shift));
                    }

                    currentSiteNum++;
                }

                currentSite++;
                numSites[0]++;
            }

            return new Tuple<>(genotypeClasses, possibleTerms);

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
                return new DominanceNormalizedIBSSiteSpliterator(myGenotypes, lo, mid, myProgressListener);
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
