/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.analysis.popgen;

import com.google.common.collect.Ordering;

/**
 *
 * @author edbuckler
 */
public class DonorHypoth implements Comparable<DonorHypoth>{
    public int targetTaxon = -1;
    public int donor1Taxon = -1;
    public int donor2Taxon = -1;
    public int startBlock = -1;
    public int focusBlock = -1;
    public int endBlock = -1;
    public double pError = 1;
    public double pHeterozygous = -1;
    public double pHomoD1 = -1;
    public double pHomoD2 = -11;
    public int totalSites = 0;
    public int mendelianErrors = 0;
    public final int startSite;
    public final int endSite;
    public byte[] phasedResults=null;
    private double errorRate=Double.NaN;
    private int prevalenceScore=Integer.MIN_VALUE;
    private int myHash=Integer.MIN_VALUE;


    public static Ordering<DonorHypoth> byErrorRateOrdering = new Ordering<DonorHypoth>() {
        public int compare(DonorHypoth left, DonorHypoth right) {
            return Double.compare(left.getErrorRate(),right.getErrorRate());
        }
    };

    public static Ordering<DonorHypoth> byPrevalenceOrdering = new Ordering<DonorHypoth>() {
        public int compare(DonorHypoth left, DonorHypoth right) {
            return Double.compare(left.getErrorRate(),right.getErrorRate());
        }
    };


    public DonorHypoth(int targetTaxon, int donor1Taxon, int donor2Taxon, int startBlock, int focusBlock, int endBlock,
                       int totalSites, int mendelianErrors) {
        this(targetTaxon, donor1Taxon, donor2Taxon, startBlock, focusBlock, endBlock);
        this.totalSites = totalSites;
        this.mendelianErrors = mendelianErrors;
        errorRate= ((double)mendelianErrors+0.5) / (double) totalSites;
    }

    public DonorHypoth(int targetTaxon, int donor1Taxon, int donor2Taxon, int startBlock, int focusBlock, int endBlock) {
        this.targetTaxon = targetTaxon;
        if (donor1Taxon < donor2Taxon) {
            this.donor1Taxon = donor1Taxon;
            this.donor2Taxon = donor2Taxon;
        } else {
            this.donor1Taxon = donor2Taxon;
            this.donor2Taxon = donor1Taxon;
        }
        this.startBlock = startBlock;
        this.focusBlock = focusBlock;
        this.endBlock = endBlock;
        this.startSite=startBlock*64;
        this.endSite=(endBlock*64)+63;
    }
    
    public byte getPhaseForSite(int site) {
        if(phasedResults==null) return (byte)1;
        return phasedResults[site-startSite];
    }

    public double getErrorRate() {
        return errorRate;
    }
    
    public boolean isInbred() {
        return (donor1Taxon==donor2Taxon);
    }
    
    public int getFocusStartSite() {
        return focusBlock*64;
    }
    
    public int getFocusEndSite() {
        return (focusBlock*64)+63;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final DonorHypoth other = (DonorHypoth) obj;
        if (this.targetTaxon != other.targetTaxon) {
            return false;
        }
        if (this.donor1Taxon != other.donor1Taxon) {
            return false;
        }
        if (this.donor2Taxon != other.donor2Taxon) {
            return false;
        }
//        if (this.focusBlock != other.focusBlock) {
//            return false;
//        }
        return true;
    }

    @Override
    public int compareTo(DonorHypoth o) {
        if(o.equals(this)) {return 0;}
        DonorHypoth aDH=o;
        if(aDH.donor1Taxon<this.donor1Taxon) return 1;
        return -1;
        
    }

    @Override
    public int hashCode() {
        if(myHash!=Integer.MIN_VALUE) return myHash;
        int hash = 5;
        hash = 83 * hash + this.targetTaxon;
        hash = 83 * hash + this.donor1Taxon;
        hash = 83 * hash + this.donor2Taxon;
        hash = 83 * hash + this.focusBlock;
        myHash=hash;
        return myHash;
    }

    public String toString() {
        return String.format("TargetTx:%d Donor11Tx:%d Donor2Tx:%d FirstBlock:%d FocusBlock:%d LastBlock:%d TotalSites:%d MendelErrors:%d ",
                targetTaxon, donor1Taxon, donor2Taxon, startBlock, focusBlock, endBlock, totalSites, mendelianErrors);
    }
    
}
