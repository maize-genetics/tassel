/*
 *  IBSDistanceMatrix3Alleles
 * 
 *  Created on Jul 10, 2015
 */
package net.maizegenetics.analysis.distance;

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
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class IBSDistanceMatrix3Alleles {

    private static final Logger myLogger = LogManager.getLogger(IBSDistanceMatrix3Alleles.class);

    private IBSDistanceMatrix3Alleles() {
        // utility
    }

    public static DistanceMatrix getInstance(GenotypeTable genotype) {
        return getInstance(genotype, 0, false, null);
    }

    public static DistanceMatrix getInstance(GenotypeTable genotype, ProgressListener listener) {
        return getInstance(genotype, 0, false, listener);
    }

    public static DistanceMatrix getInstance(GenotypeTable genotype, int minSiteComp, boolean trueIBS, ProgressListener listener) {
        return computeHetBitDistances(genotype, listener, trueIBS, minSiteComp);
    }

    private static DistanceMatrix computeHetBitDistances(GenotypeTable genotype, ProgressListener listener, boolean isTrueIBS, int minSitesComp) {

        int numSeqs = genotype.numberOfTaxa();
        double avgTotalSites = 0.0;
        long time = System.currentTimeMillis();

        Counters temp = new Counters(numSeqs);
        stream(genotype, listener).forEach((long[] t) -> {
            temp.add(t);
        });

        int[][] counters = temp.myCounters;

        DistanceMatrixBuilder builder = DistanceMatrixBuilder.getInstance(genotype.taxa());
        long count = 0;
        for (int i = 0; i < numSeqs; i++) {
            int index = 0;
            for (int j = i; j < numSeqs; j++) {
                if (j == i && !isTrueIBS) {
                    builder.set(i, i, 0.0);
                    index += 3;
                } else {
                    int sameCount = counters[i][index++];
                    int diffCount = counters[i][index++];
                    int hetCount = counters[i][index++];
                    long sites = sameCount + diffCount - hetCount;
                    double identity = ((double) (sameCount) - 0.5 * hetCount) / (double) (sites);
                    double dist = 1 - identity;

                    if (sites < minSitesComp) {
                        dist = Double.NaN;
                    }
                    builder.set(i, j, dist);
                    avgTotalSites += sites;  //this assumes not hets
                    count++;
                }
            }
        }

        avgTotalSites /= (double) count;
        
        GeneralAnnotationStorage.Builder annotations = GeneralAnnotationStorage.getBuilder();
        annotations.addAnnotation(DistanceMatrixBuilder.MATRIX_TYPE, DistanceMatrixBuilder.IBS_DISTANCE_MATRIX_TYPE);
        annotations.addAnnotation(DistanceMatrixBuilder.IBS_DISTANCE_MATRIX_NUM_ALLELES, "3");
        annotations.addAnnotation(DistanceMatrixBuilder.IBS_DISTANCE_MATRIX_TRUE_IBS, String.valueOf(isTrueIBS));
        annotations.addAnnotation(DistanceMatrixBuilder.IBS_DISTANCE_MATRIX_AVE_TOTAL_SITES, String.valueOf(avgTotalSites));
        
        builder.annotation(annotations.build());
        
        myLogger.info("IBSDistanceMatrix3Alleles: computeHetBitDistances time = " + (System.currentTimeMillis() - time) / 1000 + " seconds");
        
        return builder.build();

    }

    public static double[] computeHetDistances(byte[] first, byte[] second, int minSitesComp) {
        return null;
    }

    protected static void fireProgress(int percent, ProgressListener listener) {
        if (listener != null) {
            listener.progress(percent, null);
        }

    }

    //
    // Each CPU thread (process) creates an instance of this class
    // to acculate counters of the IBS Distance Matrix. The add()
    // method parses out the three counts from each long that's
    // coming from the stream. These are
    // combined with addAll() to result in one instance at the end.
    // Each three consecutive int holds the same, different, and het
    // count for a pair-wise comparison.
    //
    private static class Counters {

        private final int[][] myCounters;
        private final int myNumTaxa;

        public Counters(int numTaxa) {
            myNumTaxa = numTaxa;
            myCounters = new int[myNumTaxa][];
            for (int i = 0; i < myNumTaxa; i++) {
                myCounters[i] = new int[(myNumTaxa - i) * 3];
            }
        }

        public synchronized void add(long[] values) {
            int index = 0;
            for (int i = 0; i < myNumTaxa; i++) {
                for (int j = 0; j < myCounters[i].length; j += 3) {
                    myCounters[i][j] += (int) (values[index] & 0x1FFFFFl);
                    myCounters[i][j + 1] += (int) ((values[index] >>> 21) & 0x1FFFFFl);
                    myCounters[i][j + 2] += (int) ((values[index] >>> 42) & 0x1FFFFFl);
                    index++;
                }
            }
        }

        public void addAll(Counters counters) {
            int[][] other = counters.myCounters;
            for (int t = 0; t < myNumTaxa; t++) {
                for (int i = 0, n = myCounters[t].length; i < n; i++) {
                    myCounters[t][i] += other[t][i];
                }
            }
        }

    }

    //
    // These constants named whether a pair-wise comparison is
    // SAME_DIFFERENT.  The value has a 1 in the appropriate
    // 3 x 20 bits depending whether same, different, or het
    //
    private static final long TRUE_TRUE_LONG = 0x40000200001l;
    private static final long TRUE_FALSE_LONG = 0x1l;
    private static final long FALSE_TRUE_LONG = 0x200000l;
    private static final long FALSE_FALSE_LONG = 0x0l;

    //
    // This precalculates the counts for every combination
    // of three sites.
    //
    private static long[] PRECALCULATED_COUNTS = null;

    static {

        long[] possibleTerms = new long[32];
        possibleTerms[22] = TRUE_FALSE_LONG;
        possibleTerms[20] = TRUE_TRUE_LONG;
        possibleTerms[18] = TRUE_TRUE_LONG;
        possibleTerms[6] = FALSE_TRUE_LONG;
        possibleTerms[2] = FALSE_TRUE_LONG;
        possibleTerms[0] = FALSE_FALSE_LONG;
        possibleTerms[21] = TRUE_TRUE_LONG;
        possibleTerms[17] = TRUE_TRUE_LONG;
        possibleTerms[4] = TRUE_TRUE_LONG;
        possibleTerms[1] = TRUE_TRUE_LONG;
        possibleTerms[5] = FALSE_TRUE_LONG;
        possibleTerms[19] = TRUE_TRUE_LONG;
        possibleTerms[3] = TRUE_TRUE_LONG;
        possibleTerms[14] = TRUE_FALSE_LONG;
        possibleTerms[10] = TRUE_TRUE_LONG;
        possibleTerms[11] = TRUE_TRUE_LONG;
        possibleTerms[7] = TRUE_FALSE_LONG;

        PRECALCULATED_COUNTS = new long[23255];

        for (int i = 0; i < 23255; i++) {
            int firstCode = i & 0x1f;
            int secondCode = (i >>> 5) & 0x1f;
            int thirdCode = (i >>> 10) & 0x1f;
            PRECALCULATED_COUNTS[i] = possibleTerms[firstCode] + possibleTerms[secondCode] + possibleTerms[thirdCode];
        }

    }

    //
    // This defines the codes for each possible state at a given
    // site and taxon.
    //
    private static final byte[] PRECALCULATED_ENCODINGS = new byte[8];

    static {
        // 22, 21, 19, 14, 11, 7, 0
        PRECALCULATED_ENCODINGS[1] = 0x16; // Major
        PRECALCULATED_ENCODINGS[3] = 0x15; // Major and Minor
        PRECALCULATED_ENCODINGS[5] = 0x13; // Major and Second Minor
        PRECALCULATED_ENCODINGS[2] = 0xE; // Minor
        PRECALCULATED_ENCODINGS[6] = 0xB; // Minor and Second Minor
        PRECALCULATED_ENCODINGS[4] = 0x7; // Second Minor
        PRECALCULATED_ENCODINGS[0] = 0x0; // Unknown
    }

    private static final int NUM_CORES_TO_USE = TasselPrefs.getMaxThreads();

    //
    // Used to report progress.  This is not thread-safe but
    // works well enough for this purpose.
    //
    private static int myNumSitesProcessed = 0;

    private static final int MAX_NUMBER_20_BITS = 0xFFFFF;

    //
    // Creates stream from IBSSiteSpliterator and Genotype Table
    //
    private static Stream<long[]> stream(GenotypeTable genotypes, ProgressListener listener) {
        myNumSitesProcessed = 0;
        return StreamSupport.stream(new IBSSiteSpliterator(genotypes, 0, genotypes.numberOfSites(), listener), true);
    }

    //
    // Spliterator that splits the sites into halves each time for
    // processing.
    //
    static class IBSSiteSpliterator implements Spliterator<long[]> {

        private int myCurrentSite;
        private final int myFence;
        private final GenotypeTable myGenotypes;
        private final int myNumTaxa;
        private final int myNumSites;
        private final ProgressListener myProgressListener;
        private int myMinSitesToProcess;

        IBSSiteSpliterator(GenotypeTable genotypes, int currentIndex, int fence, ProgressListener listener) {
            myGenotypes = genotypes;
            myNumTaxa = myGenotypes.numberOfTaxa();
            myNumSites = myGenotypes.numberOfSites();
            myCurrentSite = currentIndex;
            myFence = fence;
            myProgressListener = listener;
            myMinSitesToProcess = myNumSites / NUM_CORES_TO_USE;
            if (myMinSitesToProcess == 0) {
                myMinSitesToProcess = myNumSites;
            }
        }

        @Override
        public void forEachRemaining(Consumer<? super long[]> action) {

            int numSitesProcessed = myFence - myCurrentSite;

            //
            // This prevents overrunning the max number that can
            // be held in 20 bits of the long.
            //
            for (; myCurrentSite < myFence;) {

                int currentBlockFence = Math.min(myCurrentSite + MAX_NUMBER_20_BITS, myFence);
                long[] counts = new long[myNumTaxa * (myNumTaxa + 1) / 2];

                for (; myCurrentSite < currentBlockFence;) {

                    int[] numSites = new int[1];

                    //
                    // Gets encodings for several blocks of sites.
                    //
                    short[] encodings1 = getBlockOfSites(myCurrentSite, numSites, currentBlockFence);

                    short[] encodings2 = getBlockOfSites(myCurrentSite + numSites[0], numSites, currentBlockFence);

                    short[] encodings3 = getBlockOfSites(myCurrentSite + numSites[0], numSites, currentBlockFence);

                    short[] encodings4 = getBlockOfSites(myCurrentSite + numSites[0], numSites, currentBlockFence);

                    short[] encodings5 = getBlockOfSites(myCurrentSite + numSites[0], numSites, currentBlockFence);

                    short[] encodings6 = getBlockOfSites(myCurrentSite + numSites[0], numSites, currentBlockFence);

                    myCurrentSite += numSites[0];

                    //
                    // Iterates through all pair-wise combinations of taxa
                    //
                    int index = 0;
                    for (int firstTaxa = 0; firstTaxa < myNumTaxa; firstTaxa++) {
                        //
                        // Can skip inter-loop if all sites for first
                        // taxon is Unknown diploid allele values
                        //
                        if ((encodings1[firstTaxa] != 0x0) || (encodings2[firstTaxa] != 0x0) || (encodings3[firstTaxa] != 0x0)
                                || (encodings4[firstTaxa] != 0x0) || (encodings5[firstTaxa] != 0x0) || (encodings6[firstTaxa] != 0x0)) {
                            for (int secondTaxa = firstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                                //
                                // Combine first taxon's encoding with
                                // second taxon's encoding to
                                // create index into pre-calculated counts
                                //
                                counts[index] += PRECALCULATED_COUNTS[encodings1[firstTaxa] & encodings1[secondTaxa]]
                                        + PRECALCULATED_COUNTS[encodings2[firstTaxa] & encodings2[secondTaxa]]
                                        + PRECALCULATED_COUNTS[encodings3[firstTaxa] & encodings3[secondTaxa]]
                                        + PRECALCULATED_COUNTS[encodings4[firstTaxa] & encodings4[secondTaxa]]
                                        + PRECALCULATED_COUNTS[encodings5[firstTaxa] & encodings5[secondTaxa]]
                                        + PRECALCULATED_COUNTS[encodings6[firstTaxa] & encodings6[secondTaxa]];
                                index++;
                            }
                        } else {
                            index += myNumTaxa - firstTaxa;
                        }
                    }
                }

                action.accept(counts);
            }
            myNumSitesProcessed += numSitesProcessed;
            fireProgress((int) ((double) myNumSitesProcessed / (double) myNumSites * 100.0), myProgressListener);
        }

        private static final int NUM_SITES_PER_BLOCK = 3;

        private short[] getBlockOfSites(int currentSite, int[] numSites, int currentBlockFence) {

            int currentSiteNum = 0;

            //
            // This holds the encoding for every taxa.  Each
            // short has encodings in each 5 bits.
            //
            short[] encodings = new short[myNumTaxa];

            while ((currentSiteNum < NUM_SITES_PER_BLOCK) && (currentSite < currentBlockFence)) {

                byte[] genotype = myGenotypes.genotypeAllTaxa(currentSite);
                int[][] alleles = AlleleFreqCache.allelesSortedByFrequencyNucleotide(genotype);
                int numAlleles = alleles[0].length;

                //
                // If whole site is Unknown, then skip the site.
                //
                if (numAlleles != 0) {

                    //
                    // Records presence of major, minor, and second minor alleles
                    // for current site in 5 bits.
                    //
                    for (int i = 0; i < myNumTaxa; i++) {
                        byte first = (byte) (genotype[i] & 0xf);
                        byte second = (byte) (genotype[i] >>> 4 & 0xf);
                        int allelePresent = 0;
                        if ((alleles[0][0] == first) || (alleles[0][0] == second)) {
                            allelePresent = 0x1;
                        }
                        if (numAlleles >= 2) {
                            if ((alleles[0][1] == first) || (alleles[0][1] == second)) {
                                allelePresent |= 0x2;
                            }
                            if (numAlleles >= 3) {
                                if ((alleles[0][2] == first) || (alleles[0][2] == second)) {
                                    allelePresent |= 0x4;
                                }
                            }
                        }

                        encodings[i] = (short) (encodings[i] << 5 | PRECALCULATED_ENCODINGS[allelePresent]);
                    }

                    currentSiteNum++;
                }

                currentSite++;
                numSites[0]++;
            }

            return encodings;

        }

        @Override
        public boolean tryAdvance(Consumer<? super long[]> action) {
            if (myCurrentSite < myFence) {

                long[] counts = new long[myNumTaxa * (myNumTaxa + 1) / 2];

                int[] numSites = new int[1];

                short[] encodings1 = getBlockOfSites(myCurrentSite, numSites, myFence);

                short[] encodings2 = getBlockOfSites(myCurrentSite + numSites[0], numSites, myFence);

                short[] encodings3 = getBlockOfSites(myCurrentSite + numSites[0], numSites, myFence);

                myCurrentSite += numSites[0];

                //
                // Iterates through all pair-wise combinations of taxa
                //
                int index = 0;
                for (int firstTaxa = 0; firstTaxa < myNumTaxa; firstTaxa++) {
                    //
                    // Can skip inter-loop if all sites for first
                    // taxon is Unknown diploid allele values
                    //
                    if ((encodings1[firstTaxa] != 0x0) || (encodings2[firstTaxa] != 0x0) || (encodings3[firstTaxa] != 0x0)) {
                        for (int secondTaxa = firstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                            //
                            // Combine first taxon's encoding with
                            // second taxon's encoding to
                            // create index into pre-calculated counts
                            //
                            counts[index] += PRECALCULATED_COUNTS[encodings1[firstTaxa] & encodings1[secondTaxa]]
                                    + PRECALCULATED_COUNTS[encodings2[firstTaxa] & encodings2[secondTaxa]]
                                    + PRECALCULATED_COUNTS[encodings3[firstTaxa] & encodings3[secondTaxa]];
                            index++;
                        }
                    } else {
                        index += myNumTaxa - firstTaxa;
                    }
                }

                action.accept(counts);

                return true;
            } else {
                return false;
            }
        }

        @Override
        /**
         * Splits sites
         */
        public Spliterator<long[]> trySplit() {
            int lo = myCurrentSite;
            int mid = lo + myMinSitesToProcess;
            if (mid < myFence) {
                myCurrentSite = mid;
                return new IBSSiteSpliterator(myGenotypes, lo, mid, myProgressListener);
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
