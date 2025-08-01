package net.maizegenetics.analysis.imputation;

import com.google.common.collect.*;
import com.google.common.primitives.Ints;
import net.maizegenetics.analysis.popgen.DonorHypoth;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import static net.maizegenetics.dna.snp.GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
import static net.maizegenetics.dna.WHICH_ALLELE.Major;
import static net.maizegenetics.dna.WHICH_ALLELE.Minor;
import static net.maizegenetics.dna.snp.GenotypeTableUtils.isHeterozygous;

/**
 * Basic utility functions to support imputation by blocks.
 *
 * @author Ed Buckler
 */
public class FILLINImputationUtils {

    /**
     * Counts union and intersection of major and minor alleles between the target genotype and potential
     * donor genotypes.  The counts are done by 64 sites block.  These counts can be quickly used to estimate
     * distance for the set of blocks.
     * @param modBitsOfTarget major and minor presence bits for target genotype (must be aligned same as donor)
     * @param donorAlign genotypeTable with potential donor genotypes
     * @return array with [donor index][sites, same count, diff count, het count index][block index]
     */
    public static byte[][][] calcAllelePresenceCountsBtwTargetAndDonors(BitSet[] modBitsOfTarget, GenotypeTable donorAlign) {
        int blocks=modBitsOfTarget[0].getNumWords();
        byte[][][] allDist=new byte[donorAlign.numberOfTaxa()][4][blocks];
        long[] iMj=modBitsOfTarget[0].getBits();
        long[] iMn=modBitsOfTarget[1].getBits();
        for (int donor1 = 0; donor1 < allDist.length; donor1++) {
            long[] jMj=donorAlign.allelePresenceForAllSites(donor1, Major).getBits();
            long[] jMn=donorAlign.allelePresenceForAllSites(donor1, Minor).getBits();
            for (int i = 0; i <blocks; i++) {
                long same = (iMj[i] & jMj[i]) | (iMn[i] & jMn[i]);
                long diff = (iMj[i] & jMn[i]) | (iMn[i] & jMj[i]);
                long hets = same & diff;
                int sameCnt = BitUtil.pop(same);
                int diffCnt = BitUtil.pop(diff);
                int hetCnt = BitUtil.pop(hets);
                int sites = sameCnt + diffCnt - hetCnt;
                allDist[donor1][2][i]=(byte)diffCnt;
                allDist[donor1][3][i]=(byte)hetCnt;
                allDist[donor1][0][i]=(byte)sites;
                allDist[donor1][1][i]=(byte)sameCnt;
            }
        }
        return allDist;
    }

    /**
     *Simple algorithm that tests every possible haplotype as a homozygous donor to minimize
     * the number of unmatched informative alleles.  Currently, there is little tie
     * breaking, longer matches are favored by add 0.5 errors to the error rate calculation within
     * DonorHypoth.  Distance calculation are made between first and last block. This code is relatively
     * fast compared to calculating distance as it works with the precomputed distances.
     * @param targetTaxon index of target taxon only used to annotated DonorHypoth
     * @param firstBlock  index of first 64 site block
     * @param lastBlock  inclusive index of last 64 site block
     * @param focusBlock  index of the focus block (only used for annotation of DonorHypoth)
     * @param donor1indices  array of donor indices to tests
     * @param targetToDonorDistances precomputed block distances
     * @param minTestSites minimum number of comparable sites to be included the analysis
     * @param maxDonorHypotheses maximum number of donors to retain
     * @return sorted list of homozygous donors
     */
    public static DonorHypoth[] findHomozygousDonorHypoth(int targetTaxon, int firstBlock, int lastBlock, int focusBlock,
                                                          int[] donor1indices, byte[][][] targetToDonorDistances,
                                                          int minTestSites, int maxDonorHypotheses) {
        MinMaxPriorityQueue<DonorHypoth> bestDonors=MinMaxPriorityQueue.orderedBy(DonorHypoth.byErrorRateOrdering)
                .maximumSize(maxDonorHypotheses).create();
        for (int d1 : donor1indices) {
            int sameCnt = 0, diffCnt = 0, hetCnt = 0;
            for (int i = firstBlock; i <=lastBlock; i++) {
                sameCnt+=targetToDonorDistances[d1][1][i];
                diffCnt+=targetToDonorDistances[d1][2][i];
                hetCnt+=targetToDonorDistances[d1][3][i];
            }
            int testSites= sameCnt + diffCnt - hetCnt;
            if(testSites<minTestSites) continue;
            int totalMendelianErrors=diffCnt-(hetCnt/2);
            DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d1, firstBlock,
                    focusBlock, lastBlock, testSites, totalMendelianErrors);
            bestDonors.add(theDH);
        }
        DonorHypoth[] result=bestDonors.toArray(new DonorHypoth[0]);
        Arrays.sort(result,DonorHypoth.byErrorRateOrdering);  //Ques keep the top values, but not ordered
        return result;
    }

    /**
     * Produces a sort list of most prevalent donors across the donorAlignment.
     * It looks at all focus blocks and weights blocks by rank
     * @param allDH homozygous donor hypotheses across entire region
     * @param maxHypotheses min number of focus blocks that the hypotheses needs to show up in
     * @return array of indices of the best donors
     */
    public static int[] mostFrequentDonorsAcrossFocusBlocks(DonorHypoth[][] allDH, int maxHypotheses) {
        Map<Integer,Integer> bd=new HashMap<>();
        for (int i = 0; i < allDH.length; i++) {
            int rank=allDH[i].length;
            for (DonorHypoth dh: allDH[i]) {
                if (dh==null) continue;
                if (bd.containsKey(dh.donor1Taxon)) {  //counts the frequency of best hit
                    bd.put(dh.donor1Taxon, bd.get(dh.donor1Taxon)+ rank);
                } else {
                    bd.put(dh.donor1Taxon, rank);
                }
                rank--;
                //break;//this line allows only the best hypothesis to go into set (original method)
            }
        }
        TreeMultimap<Integer,Integer> rankings= TreeMultimap.create(Ordering.natural().reverse(),Ordering.natural());
        for (Map.Entry<Integer, Integer> dhs : bd.entrySet()) {
            rankings.put(dhs.getValue(),dhs.getKey());
        }
        int resultSize=(rankings.size()<maxHypotheses)?rankings.size():maxHypotheses;
        int[] result=new int[resultSize];
        Iterator<Integer> iDH=rankings.values().iterator();
        for (int i=0; i<resultSize; i++) {
            result[i]=iDH.next();
        }
        return result;
    }

    /**
     * Produces a sort list of most prevalent donorHypotheses across the donorAlign.
     * It looks at all focus blocks and weight blocks by rank
     * @param targetToDonorDistances precomputed block distances
     * @param minTestSites minimum number of comparable sites to be included the analysis
     * @param maxDonorHypotheses maximum number of donor hypotheses to retain
     * @return array of indices of the best donors
     */
    public static int[] bestDonorsAcrossEntireRegion(byte[][][] targetToDonorDistances,
                                                     int minTestSites, int maxDonorHypotheses) {
        int[] donor1indices=fillInc(0,targetToDonorDistances.length-1);
        DonorHypoth[] bestDH=findHomozygousDonorHypoth(-1,0,targetToDonorDistances[0][0].length-1,0,donor1indices,
                targetToDonorDistances,minTestSites,maxDonorHypotheses);
        int resultSize=(bestDH.length<maxDonorHypotheses)?bestDH.length:maxDonorHypotheses;
        int[] result=new int[resultSize];
        int i=0;
        for (DonorHypoth donorHypoth : bestDH) {
            result[i++]=donorHypoth.donor1Taxon;
        }
        return result;
    }

    public static int sumOf(int... integers) {
        int total = 0;
        for (int i = 0; i < integers.length; total += integers[i++]);
        return total;
    }

    public static int sumOf(byte... integers) {
        int total = 0;
        for (int i = 0; i < integers.length; total += integers[i++]);
        return total;
    }

    public static int[] fillInc(int first, int last) {
        int[] total = new int[last-first+1];
        for (int i = 0; i < total.length; i++) total[i]=first++;
        return total;
    }



    /**
     *
     * Simple algorithm that does a one dimensional test of donors.   A single donor with a list of others donor
     * combinations to minimize the number of unmatched informative alleles.
     * @param targetTaxon index of target taxon only used to annotated DonorHypoth
     * @param mjT masked bitset for the major allele
     * @param mnT  masked bitset for the minor allele
     * @param firstBlock index of first 64 site block
     * @param lastBlock  inclusive index of last 64 site block
     * @param focusBlock index of the focus block (only used for annotation of DonorHypoth)
     * @param donorAlign genotypeTable with potential donor genotypes
     * @param d1 fixed donor
     * @param donor2Indices list of second potential donors
     * @param maxDonorHypotheses maximum number of donor hypotheses to retain
     * @param minTestSites minimum number of comparable sites to be included the analysis
     * @return  array of DonorHypoth sorted by error rate
     */
    public static DonorHypoth[] findHeterozygousDonorHypoth(int targetTaxon, long[] mjT, long[] mnT,
              int firstBlock, int lastBlock, int focusBlock, GenotypeTable donorAlign, int d1, int[] donor2Indices,
                                              int maxDonorHypotheses, int minTestSites) {
        MinMaxPriorityQueue<DonorHypoth> bestDonors=MinMaxPriorityQueue.orderedBy(DonorHypoth.byErrorRateOrdering)
                .maximumSize(maxDonorHypotheses).create();
        long[] mj1=donorAlign.allelePresenceForSitesBlock(d1, Major, firstBlock, lastBlock + 1);
        long[] mn1=donorAlign.allelePresenceForSitesBlock(d1, Minor, firstBlock, lastBlock + 1);
        for (int d2 : donor2Indices) {
            long[] mj2=donorAlign.allelePresenceForSitesBlock(d2, Major, firstBlock, lastBlock + 1);
            long[] mn2=donorAlign.allelePresenceForSitesBlock(d2, Minor, firstBlock, lastBlock + 1);
            int[] mendErr=mendelErrorComparison(mjT, mnT, mj1, mn1, mj2, mn2);
            if(mendErr[1]<minTestSites) continue;
            DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d2, firstBlock,
                    focusBlock, lastBlock, mendErr[1], mendErr[0]);
            bestDonors.add(theDH);
        }
        DonorHypoth[] result=bestDonors.toArray(new DonorHypoth[0]);
        Arrays.sort(result,DonorHypoth.byErrorRateOrdering);  //Ques keep the top values, but not ordered
        return result;
    }

    /**
     *
     * Simple algorithm that does a two dimensional test of donors.   A list donor versus with a list of others donor
     * combinations to minimize the number of unmatched informative alleles.
     * @param targetTaxon index of target taxon only used to annotated DonorHypoth
     * @param mjT masked bitset for the major allele
     * @param mnT  masked bitset for the minor allele
     * @param firstBlock index of first 64 site block
     * @param lastBlock  inclusive index of last 64 site block
     * @param focusBlock index of the focus block (only used for annotation of DonorHypoth)
     * @param donorAlign genotypeTable with potential donor genotypes
     * @param donor1Indices fixed donor
     * @param donor2Indices list of second potential donors
     * @param maxDonorHypotheses maximum number of donor hypotheses to retain
     * @param minTestSites minimum number of comparable sites to be included the analysis
     * @return  array of DonorHypoth sorted by error rate
     */
    public static DonorHypoth[] findHeterozygousDonorHypoth(int targetTaxon, long[] mjT, long[] mnT,
                 int firstBlock, int lastBlock, int focusBlock, GenotypeTable donorAlign, int[] donor1Indices, int[] donor2Indices,
                                                            int maxDonorHypotheses, int minTestSites) {
        Multimap<Integer,Integer> tests=HashMultimap.create();
        for (int d1 : donor1Indices) {
            for (int d2 : donor2Indices) {
                if(d1<d2) {tests.put(d1,d2);}
                else {tests.put(d2,d1);}
            }
        }
        MinMaxPriorityQueue<DonorHypoth> bestDonors=MinMaxPriorityQueue.orderedBy(DonorHypoth.byErrorRateOrdering)
                .maximumSize(maxDonorHypotheses).create();
        for (int d1 : tests.keySet()) {
            int[] d2donors=Ints.toArray(tests.get(d1));
            DonorHypoth[] oneDimenHypoth=findHeterozygousDonorHypoth(targetTaxon, mjT, mnT, firstBlock, lastBlock,
                    focusBlock, donorAlign, d1, d2donors, maxDonorHypotheses, minTestSites);
            for (DonorHypoth donorHypoth : oneDimenHypoth) {bestDonors.add(donorHypoth);}
        }
        DonorHypoth[] result=bestDonors.toArray(new DonorHypoth[0]);
        Arrays.sort(result,DonorHypoth.byErrorRateOrdering);  //Ques keep the top values, but not ordered
        return result;
    }

    /**
     * Combines arrays of donorHypoth, sorts them, and returns the best limited by maxDonorHypotheses
     * @param maxDonorHypotheses maximum number of donor hypotheses to retain
     * @param dhs arrays of DonorHypoth[]
     * @return combined arrays of DonorHypoth (length<=maxDonorHypotheses) ordered by error rate
     */
    public static DonorHypoth[] combineDonorHypothArrays(int maxDonorHypotheses, DonorHypoth[] ... dhs) {
        MinMaxPriorityQueue<DonorHypoth> bestDonors=MinMaxPriorityQueue.orderedBy(DonorHypoth.byErrorRateOrdering)
                .maximumSize(maxDonorHypotheses).create();
        for (DonorHypoth[] aDHArray : dhs) {
            for (DonorHypoth donorHypoth : aDHArray) {
                bestDonors.add(donorHypoth);
            }
        }
        DonorHypoth[] result=bestDonors.toArray(new DonorHypoth[0]);
        Arrays.sort(result,DonorHypoth.byErrorRateOrdering);  //Ques keep the top values, but not ordered
        return result;
    }

    /**
     * Given a start 64 site block, it expands to the left and right until it hits
     * the minimum Minor Site count or MajorSiteCount in the target taxon
     * @param mnT - minor allele bit presence in a series of longs
     * @param focusBlock  index of the focus block
     * @param minMinorCnt minimum count to stop expanding for minor allele sites
     * @param minMajorCnt minimum count to stop expanding for major allele sites
     * @return arrays of blocks {startBlock, focusBlock, endBlock}
     */
    public static int[] getBlockWithMinMinorCount(long[] mjT, long[] mnT, int focusBlock, int minMinorCnt,
                                            int minMajorCnt) {
        int blocks=mjT.length;
        int majorCnt=Long.bitCount(mjT[focusBlock]);
        int minorCnt=Long.bitCount(mnT[focusBlock]);
        int endBlock=focusBlock, startBlock=focusBlock;
        while((minorCnt<minMinorCnt)&&(majorCnt<minMajorCnt)) {
            boolean preferMoveStart=(focusBlock-startBlock<endBlock-focusBlock)?true:false;
            if(startBlock==0) {preferMoveStart=false;}
            if(endBlock==blocks-1) {preferMoveStart=true;}
            if((startBlock==0)&&(endBlock==blocks-1)) break;
            if(preferMoveStart) {//expand start
                startBlock--;
                minorCnt+=Long.bitCount(mnT[startBlock]);
                majorCnt+=Long.bitCount(mjT[startBlock]);
            } else { //expand end
                endBlock++;
                minorCnt+=Long.bitCount(mnT[endBlock]);
                majorCnt+=Long.bitCount(mjT[startBlock]);
            }
        }
        int[] result={startBlock, focusBlock, endBlock};
        return result;
    }

    /**
     * Determines the number of sites in which the target (T) sequence cannot be explained by the genotypes of
     * either donor (1 & 2).  Only sites where the genotype for all taxa can be tested.
     * @param mjT major allele bits of target
     * @param mnT minor allele bits of target
     * @param mj1 major allele bits of donor 1
     * @param mn1 minor allele bits of donor 1
     * @param mj2 major allele bits of donor 2
     * @param mn2 minor allele bits of donor 2
     * @return array with [count of mendelian errors, total sites tested]
     */
    public static int[] mendelErrorComparison(long[] mjT, long[] mnT, long[] mj1, long[] mn1,
                                        long[] mj2, long[] mn2) {
        int mjUnmatched=0;
        int mnUnmatched=0;
        int testSites=0;
        for (int i = 0; i < mjT.length; i++) {
            long siteMask=(mjT[i]|mnT[i])&(mj1[i]|mn1[i])&(mj2[i]|mn2[i]);
            mjUnmatched+=Long.bitCount(siteMask&mjT[i]&(mjT[i]^mj1[i])&(mjT[i]^mj2[i]));
            mnUnmatched+=Long.bitCount(siteMask&mnT[i]&(mnT[i]^mn1[i])&(mnT[i]^mn2[i]));
            testSites+=Long.bitCount(siteMask);
        }
        int totalMendelianErrors=mjUnmatched+mnUnmatched;
        return (new int[] {totalMendelianErrors, testSites});
    }

    /**
     * Sums the number of unknown and heterozgyous sites in a byte genotype
     * @param a a byte genotype
     * @return sum of unknown and heterozgyous sites
     */
    public static int[] countUnknownAndHeterozygotes(byte[] a) {
        int cnt=0, cntHets=0;
        for (int i = 0; i < a.length; i++) {
            if(a[i]==UNKNOWN_DIPLOID_ALLELE) {cnt++;}
            else if(isHeterozygous(a[i])) {cntHets++;}
        }
        return new int[]{cnt,cntHets};
    }
}
