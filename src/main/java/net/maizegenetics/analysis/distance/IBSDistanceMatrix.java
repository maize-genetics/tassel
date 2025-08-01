/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and 
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General 
 * public license. 
 *
 */
package net.maizegenetics.analysis.distance;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.ProgressListener;

import java.util.stream.IntStream;

import static net.maizegenetics.dna.WHICH_ALLELE.Major;
import static net.maizegenetics.dna.WHICH_ALLELE.Minor;
import static net.maizegenetics.dna.WHICH_ALLELE.Minor2;
import net.maizegenetics.util.Tuple;

/**
 * This class calculates an identity by state matrix. It is scaled so only
 * non-missing comparison are used. It conducts bit level calculations of IBS
 * for genotypes. Only the two most common alleles are used in the distance
 * calculations.
 * <p>
 * Please note that when heterozygous genotypes are used, Het to Het distance is
 * 0.5 NOT 0.0. The default along the identity diagonal is 0 (isTrueIBS =
 * false), but changing isTrueIBS = true will calculate the identity.
 * <p>
 * The distance estimates become wildly inaccurate when too few sites are used
 * to calculate distance. The minSiteComp parameter can be used to control the
 * minimum number of sites used for a calculation. If there are insufficient
 * sites in the estimate, then Double.NaN is returned.
 *
 * @author Ed Buckler
 * @version 1.0
 */
public class IBSDistanceMatrix {

    private IBSDistanceMatrix() {
        // utility
    }

    /**
     * Compute observed distances for all taxa. Missing sites are ignored.
     *
     * @param theAlignment Alignment used to computed distances
     */
    public static DistanceMatrix getInstance(GenotypeTable theAlignment) {
        return getInstance(theAlignment, null);
    }

    /**
     * Compute observed distances for all taxa. Missing sites are ignored.
     *
     * @param theAlignment Alignment used to computed distances
     * @param listener Listener to track progress in calculations
     */
    public static DistanceMatrix getInstance(GenotypeTable theAlignment, ProgressListener listener) {
        return getInstance(theAlignment, 0, listener);
    }

    /**
     * Compute observed distances for all taxa. Missing sites are ignored.
     *
     * @param theAlignment Alignment used to computed distances
     * @param minSiteComp Minimum number of sites needed to estimate distance
     * @param listener Listener to track progress in calculations
     */
    public static DistanceMatrix getInstance(GenotypeTable theAlignment, int minSiteComp, ProgressListener listener) {
        return getInstance(theAlignment, minSiteComp, false, listener, true);
    }

    /**
     * Compute observed distances for all taxa. Missing sites are ignored.
     *
     * @param theAlignment Alignment used to computed distances
     * @param minSiteComp Minimum number of sites needed to estimate distance
     * @param trueIBS estimate diagonal distance based IBS (default = false,
     * i=i=0.0)
     * @param listener Listener to track progress in calculations
     * @param useThirdState
     */
    public static DistanceMatrix getInstance(GenotypeTable theAlignment, int minSiteComp, boolean trueIBS, ProgressListener listener, boolean useThirdState) {
        if (useThirdState) {
            return IBSDistanceMatrix3Alleles.getInstance(theAlignment, minSiteComp, trueIBS, listener);
        } else {
            return IBSDistanceMatrix2Alleles.getInstance(theAlignment, minSiteComp, trueIBS, listener);
        }
    }

    public static double[] computeHetDistances(byte[] first, byte[] second, int minSitesComp) {
        setupHetBitDistancesIncrementors();
        long counts = 0;
        for (int i = 0; i < first.length; i++) {
            int key1 = (first[i] & 0x70) >>> 1 | (first[i] & 0x7);
            int key2 = (second[i] & 0x70) >>> 1 | (second[i] & 0x7);
            counts += INCREMENT_FUNCTIONS[key1][key2];
        }

        int sameCount = (int) (counts & 0x1FFFFFl);
        int diffCount = (int) ((counts >>> 21) & 0x1FFFFFl);
        int hetCount = (int) ((counts >>> 42) & 0x1FFFFFl);
        long sites = sameCount + diffCount - hetCount;
        double identity = ((double) (sameCount) - 0.5 * hetCount) / (double) (sites);
        double dist = 1 - identity;

        if (sites < minSitesComp) {
            dist = Double.NaN;
        }

        return new double[]{dist, sites};
    }

    /**
     * This is a cleanest, fastest and most accurate way to calculate distance.
     */
    private static Tuple<double[][], Double> computeHetBitDistances(boolean useThirdState, ProgressListener listener, GenotypeTable theTBA, boolean isTrueIBS, int minSitesComp) {
        int numSeqs = theTBA.numberOfTaxa();
        double avgTotalSites = 0;
        //LongAdder count = new LongAdder();
        //note this distance object is modified by a parallel stream, but each element is only touched once
        double[][] distance = new double[numSeqs][numSeqs];
        long numberOfTests = numSeqs * (numSeqs - 1) / 2;
        long time = System.currentTimeMillis();
        IntStream.range(0, numSeqs).parallel().forEach(i -> {
            long[] iMj = theTBA.allelePresenceForAllSites(i, Major).getBits();
            long[] iMn = theTBA.allelePresenceForAllSites(i, Minor).getBits();
            long[] iMn2 = null;
            if (useThirdState) {
                iMn2 = theTBA.allelePresenceForAllSites(i, Minor2).getBits();
            }
            for (int j = i; j < numSeqs; j++) {
                if (j == i && !isTrueIBS) {
                    distance[i][i] = 0;
                } else {
                    long[] jMj = theTBA.allelePresenceForAllSites(j, Major).getBits();
                    long[] jMn = theTBA.allelePresenceForAllSites(j, Minor).getBits();
                    double[] result;
                    if (useThirdState) {
                        long[] jMn2 = theTBA.allelePresenceForAllSites(j, Minor2).getBits();
                        result = computeHetBitDistancesThirdState(iMj, iMn, iMn2, jMj, jMn, jMn2, minSitesComp);
                    } else {
                        result = computeHetBitDistances(iMj, iMn, jMj, jMn, minSitesComp);
                    }
                    distance[i][j] = distance[j][i] = result[0];
                    //avgTotalSites += result[1];  //this assumes not hets
                    //count.increment();
                }
            }
            //fireProgress((int) ((double) (i + 1) / (double) numSeqs * 100.0), listener);
        });
        //avgTotalSites /= (double) count.longValue();
        System.out.println("computeHetBitDistances time = " + (System.currentTimeMillis() - time) / 1000 + " seconds");
        return new Tuple<>(distance, avgTotalSites);
    }

    /**
     * Compute distance for a pair of taxa.
     *
     * @param theTBA input alignment
     * @param taxon1 index of taxon 1
     * @param taxon2 index of taxon 2
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(GenotypeTable theTBA, int taxon1, int taxon2) {
        return computeHetBitDistances(theTBA, taxon1, taxon2, 0);
    }

    /**
     * Compute distance for a pair of taxa.
     *
     * @param theTBA input alignment
     * @param taxon1 index of taxon 1
     * @param taxon2 index of taxon 2
     * @param minSitesCompared Minimum number of sites needed to estimate
     * distance
     *
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(GenotypeTable theTBA, int taxon1, int taxon2, int minSitesCompared) {
        long[] iMj = theTBA.allelePresenceForAllSites(taxon1, Major).getBits();
        long[] iMn = theTBA.allelePresenceForAllSites(taxon1, Minor).getBits();
        long[] iMn2 = theTBA.allelePresenceForAllSites(taxon1, Minor2).getBits();
        long[] jMj = theTBA.allelePresenceForAllSites(taxon2, Major).getBits();
        long[] jMn = theTBA.allelePresenceForAllSites(taxon2, Minor).getBits();
        long[] jMn2 = theTBA.allelePresenceForAllSites(taxon2, Minor2).getBits();
        return computeHetBitDistancesThirdState(iMj, iMn, iMn2, jMj, jMn, jMn2, minSitesCompared, 0, iMj.length - 1);
    }

    /**
     * Compute distance for a pair of taxa. Optimized for calculations sites
     * within a certain range of underlying word (64 sites chunks) in the TBit
     * array
     *
     * @param theTBA input alignment
     * @param taxon1 index of taxon 1
     * @param taxon2 index of taxon 2
     * @param minSitesCompared Minimum number of sites needed to estimate
     * distance
     * @param firstWord starting word for calculating distance
     * site=(firstWord*64)
     * @param lastWord ending word for calculating distance inclusive
     * site=(lastWord*64+63)
     * @param maskBadSet Optional mask for sites (those set to 1 are kept)
     *
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(GenotypeTable theTBA, int taxon1, int taxon2,
            int minSitesCompared, int firstWord, int lastWord, BitSet maskBadSet) {
        long[] iMj = theTBA.allelePresenceForAllSites(taxon1, Major).getBits();
        long[] iMn = theTBA.allelePresenceForAllSites(taxon1, Minor).getBits();
        long[] iMn2 = theTBA.allelePresenceForAllSites(taxon1, Minor2).getBits();
        if (maskBadSet != null) {
            long[] maskBad = maskBadSet.getBits();
            for (int i = 0; i < iMj.length; i++) {
                iMj[i] = iMj[i] & maskBad[i];
            }
            for (int i = 0; i < iMn.length; i++) {
                iMn[i] = iMn[i] & maskBad[i];
            }
            for (int i = 0; i < iMn2.length; i++) {
                iMn2[i] = iMn2[i] & maskBad[i];
            }
        }
        long[] jMj = theTBA.allelePresenceForAllSites(taxon2, Major).getBits();
        long[] jMn = theTBA.allelePresenceForAllSites(taxon2, Minor).getBits();
        long[] jMn2 = theTBA.allelePresenceForAllSites(taxon2, Minor2).getBits();
        return computeHetBitDistancesThirdState(iMj, iMn, iMn2, jMj, jMn, jMn2, minSitesCompared, firstWord, lastWord);
    }

    /**
     * Calculation of distance using the bit vector of major and minor alleles.
     *
     * @param iMj Vector of major alleles for taxon i
     * @param iMn Vector of minor alleles for taxon i
     * @param jMj Vector of major alleles for taxon j
     * @param jMn Vector of minor alleles for taxon j
     * @param minSitesCompared Minimum number of sites needed to estimate
     * distance
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(long[] iMj, long[] iMn, long[] jMj, long[] jMn, int minSitesCompared) {
        return computeHetBitDistances(iMj, iMn, jMj, jMn, minSitesCompared, 0, iMj.length - 1);
    }

    /**
     * Calculation of distance using the bit vector of the first three alleles.
     *
     * @param iMj Vector of major alleles for taxon i
     * @param iMn Vector of minor alleles for taxon i
     * @param jMj Vector of major alleles for taxon j
     * @param jMn Vector of minor alleles for taxon j
     * @param minSitesCompared Minimum number of sites needed to estimate
     * distance
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistancesThirdState(long[] iMj, long[] iMn, long[] iMn2, long[] jMj, long[] jMn, long[] jMn2,
            int minSitesCompared) {
        return computeHetBitDistancesThirdState(iMj, iMn, iMn2, jMj, jMn, jMn2, minSitesCompared, 0, iMj.length - 1);
    }

    /**
     * Calculation of distance using the bit vector of major and minor alleles.
     *
     * @param iMj Vector of major alleles for taxon i
     * @param iMn Vector of minor alleles for taxon i
     * @param jMj Vector of major alleles for taxon j
     * @param jMn Vector of minor alleles for taxon j
     * @param minSitesCompared Minimum number of sites needed to estimate
     * distance
     * @param firstWord first world for calculating distance
     * @param lastWord last word for calculating distance inclusive
     * site=(endWord*64+63)
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(long[] iMj, long[] iMn, long[] jMj, long[] jMn,
            int minSitesCompared, int firstWord, int lastWord) {
        int sameCnt = 0, diffCnt = 0, hetCnt = 0;
        for (int x = firstWord; x <= lastWord; x++) {
            long same = (iMj[x] & jMj[x]) | (iMn[x] & jMn[x]);
            long diff = (iMj[x] & jMn[x]) | (iMn[x] & jMj[x]);
            long hets = same & diff;
            sameCnt += BitUtil.pop(same);
            diffCnt += BitUtil.pop(diff);
            hetCnt += BitUtil.pop(hets);
        }
        int sites = sameCnt + diffCnt - hetCnt;
        double identity = ((double) (sameCnt) - 0.5 * hetCnt) / (double) (sites);
        double dist = 1 - identity;
        if (sites > minSitesCompared) {
            return new double[]{dist, (double) sites};
        } else {
            return new double[]{Double.NaN, (double) sites};
        }
    }

    /**
     * Calculation of distance using the bit vectors of the first three alleles.
     *
     * @param iMj Vector of major alleles for taxon i
     * @param iMn Vector of minor alleles for taxon i
     * @param iMn2 Vector of second minor alleles for taxon i
     * @param jMj Vector of major alleles for taxon j
     * @param jMn Vector of minor alleles for taxon j
     * @param jMn2 Vector of second minor alleles for taxon j
     * @param minSitesCompared Minimum number of sites needed to estimate
     * distance
     * @param firstWord first world for calculating distance
     * @param lastWord last word for calculating distance inclusive
     * site=(endWord*64+63)
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistancesThirdState(long[] iMj, long[] iMn, long[] iMn2, long[] jMj, long[] jMn, long[] jMn2,
            int minSitesCompared, int firstWord, int lastWord) {
        int sameCnt = 0, diffCnt = 0, hetCnt = 0;
        for (int x = firstWord; x <= lastWord; x++) {
            long same = (iMj[x] & jMj[x]) | (iMn[x] & jMn[x]) | (iMn2[x] & jMn2[x]);
            long diff = (iMj[x] & jMn[x]) | (iMn[x] & jMj[x]) | (iMj[x] & jMn2[x])
                    | (iMn2[x] & jMj[x]) | (iMn[x] & jMn2[x]) | (iMn2[x] & jMn[x]);
            long hets = same & diff;
            sameCnt += BitUtil.pop(same);
            diffCnt += BitUtil.pop(diff);
            hetCnt += BitUtil.pop(hets);
        }
        int sites = sameCnt + diffCnt - hetCnt;
        double identity = ((double) (sameCnt) - 0.5 * hetCnt) / (double) (sites);
        double dist = 1 - identity;
        if (sites > minSitesCompared) {
            return new double[]{dist, (double) sites};
        } else {
            return new double[]{Double.NaN, (double) sites};
        }
    }

    private static final long TRUE_TRUE_LONG = 0x40000200001l;
    private static final long TRUE_FALSE_LONG = 0x1l;
    private static final long FALSE_TRUE_LONG = 0x200000l;
    private static final long FALSE_FALSE_LONG = 0x0l;
    private static long[][] INCREMENT_FUNCTIONS = null;

    private static void setupHetBitDistancesIncrementors() {

        if (INCREMENT_FUNCTIONS != null) {
            return;
        }

        INCREMENT_FUNCTIONS = new long[64][64];

        int[] possibleValues = new int[]{0, 1, 2, 3, 4, 5, 15};

        for (int a = 0; a < 7; a++) {
            for (int b = 0; b < 7; b++) {
                for (int c = 0; c < 7; c++) {
                    for (int d = 0; d < 7; d++) {

                        int key1 = (possibleValues[a] & 0x7) << 3 | (possibleValues[b] & 0x7);
                        int key2 = (possibleValues[c] & 0x7) << 3 | (possibleValues[d] & 0x7);

                        byte[] target = new byte[2];
                        target[0] = (byte) possibleValues[a];
                        target[1] = (byte) possibleValues[b];

                        byte[] match = new byte[2];
                        match[0] = (byte) possibleValues[c];
                        match[1] = (byte) possibleValues[d];

                        int[] result = new int[3];
                        result[0] = 0;
                        result[1] = 0;

                        if (target[0] != GenotypeTable.UNKNOWN_ALLELE) {
                            if ((target[0] == match[0]) || (target[0] == match[1])) {
                                result[0] = 1;
                            }
                            if ((match[0] != GenotypeTable.UNKNOWN_ALLELE) && (match[0] != target[0])) {
                                result[1] = 1;
                            } else if ((match[1] != GenotypeTable.UNKNOWN_ALLELE) && (match[1] != target[0])) {
                                result[1] = 1;
                            }
                        }
                        if ((result[0] == 1) && (result[1] == 1)) {
                            INCREMENT_FUNCTIONS[key1][key2] = TRUE_TRUE_LONG;
                            continue;
                        }
                        if (target[1] != GenotypeTable.UNKNOWN_ALLELE) {
                            if ((target[1] == match[0]) || (target[1] == match[1])) {
                                result[0] = 1;
                            }
                            if ((match[0] != GenotypeTable.UNKNOWN_ALLELE) && (match[0] != target[1])) {
                                result[1] = 1;
                            } else if ((match[1] != GenotypeTable.UNKNOWN_ALLELE) && (match[1] != target[1])) {
                                result[1] = 1;
                            }
                        }

                        if ((result[0] == 1) && (result[1] == 1)) {
                            INCREMENT_FUNCTIONS[key1][key2] = TRUE_TRUE_LONG;
                            continue;
                        }

                        if (result[0] == 1) {
                            INCREMENT_FUNCTIONS[key1][key2] = TRUE_FALSE_LONG;
                        } else if (result[1] == 1) {
                            INCREMENT_FUNCTIONS[key1][key2] = FALSE_TRUE_LONG;
                        } else {
                            INCREMENT_FUNCTIONS[key1][key2] = FALSE_FALSE_LONG;
                        }

                    }
                }
            }
        }

    }

}
