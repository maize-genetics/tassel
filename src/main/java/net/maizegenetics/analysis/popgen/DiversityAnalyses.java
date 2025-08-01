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
package net.maizegenetics.analysis.popgen;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.util.AbstractTableReport;
import net.maizegenetics.util.TableReport;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.DistanceMatrixBuilder;

/**
 *This method calculated estimates of nucleotide diversity (pi, theta, etc).
 *
 * In order to provide high performance, it establishes one nucleotide SitePattern, and then
 * manipulates the weighting patterns for sliding windows and different types of sites.
 *
 *Total segregating sites needs to be adjusted for missing and gap data, as this will modify Theta and Tajima estimates
 *
 * @author Ed Buckler
 * @version 1.0
 */
public class DiversityAnalyses extends AbstractTableReport implements TableReport, Serializable {
    
    private static final int NUM_OF_COLUMNS = 14;

    /** the limits of analysis */
    int startSite, endSite;
    /** values for the sliding window*/
    int step, window;
    /** number of site in the subsetted SitePattern*/
    int currNumSites;
    /** whether to use a sliding window*/
    boolean slideWindow = false;
    /** the base annotated alignment, it should be unprocessed (raw alignment)*/
    GenotypeTable theAAlignment;
    /** the nature of the sites being evaluated, currently set to ALL POLYMORPHISMs, need to add specifics indels and SNPs.*/

    List<DiversityResults> diversityResultsVector = new ArrayList<>();
    PolymorphismDistribution thePolymorphismDistribution = null;

    /**
     * Constructor that uses a sliding window
     * @param aa an annotated alignment
     * @param slidingWindow if true a sliding window will be run
     *
     */
    public DiversityAnalyses(GenotypeTable aa, boolean slidingWindow, int start, int end, int window, int step,
            PolymorphismDistribution thePolymorphismDistribution) {
        this.startSite = start;
        this.endSite = end;
        this.step = step;
        this.window = window;
        this.slideWindow = slidingWindow;
        this.theAAlignment = aa;

        this.thePolymorphismDistribution = thePolymorphismDistribution;
        runAnalyses();
    }

    public DiversityAnalyses(GenotypeTable aa, boolean slidingWindow, int start, int end, int window, int step) {
        this(aa, slidingWindow, start, end, window, step, null);
    }

    /**
     * This will run the basic analyses and controls the sliding window analyses
     */
    private void runAnalyses() {
        if (!slideWindow) {
            runAnalysisForRegion(startSite, endSite);
            return;
        }
        //slide the window
        for (int i = startSite; i < (endSite - window); i += step) {
            runAnalysisForRegion(i, i + window - 1);
        }
    }

    /**
     * This will determine what analyses are to be run and run them
     */
    private void runAnalysisForRegion(int start, int end) {
        Chromosome locus = theAAlignment.chromosome(start);
        String chromosome = locus.getName();
        double startChrPosition = theAAlignment.chromosomalPosition(start);
        double endChrPosition = theAAlignment.chromosomalPosition(end);
        GenotypeTable theFilteredAlignment = FilterGenotypeTable.getInstance(theAAlignment, start, end);
        DistanceMatrix adm = IBSDistanceMatrix.getInstance(theFilteredAlignment);
        diversityResultsVector.add(evaluate(theFilteredAlignment, adm, start, end, chromosome, startChrPosition, endChrPosition));
        if (thePolymorphismDistribution != null) {
            thePolymorphismDistribution.addDistribution("ALL" + "s" + start + "-e" + end, theFilteredAlignment, true);
        }
    }

    DiversityResults evaluate(GenotypeTable theAlignment, DistanceMatrix dm,
            int start, int end, String chromosome, double startChrPosition, double endChrPosition) {
        int sites = end - start + 1;
        DiversityResults theDiversityResults = new DiversityResults(start, end, chromosome, startChrPosition, endChrPosition);
        if (dm == null) {
            theDiversityResults.totalSites = 0;
            theDiversityResults.pipbp = Double.NaN;
            theDiversityResults.thetapbp = Double.NaN;
            theDiversityResults.segregatingSites = 0;
            return theDiversityResults;
        }
        double pipbp = dm.meanDistance();
        int segSites = countSegregatingSites(theAlignment);
        int taxa = theAlignment.numberOfTaxa();
        theDiversityResults.pipbp = pipbp;
        theDiversityResults.avgSiteCoverage = dm.annotations().getQuantAnnotation(DistanceMatrixBuilder.IBS_DISTANCE_MATRIX_AVE_TOTAL_SITES)[0];
        //theDiversityResults.avgSiteCoverage = dm.getAverageTotalSites();
        theDiversityResults.totalSites = sites;
        theDiversityResults.segregatingSites = segSites;
        theDiversityResults.thetapbp = estimateThetaPerbp(segSites, sites, theDiversityResults.avgSiteCoverage, taxa);

        // theDiversityResults.theta=estimateTheta(segSites,sites,theDiversityResults.avgSiteCoverage, theAlignment.numberOfTaxa());
        theDiversityResults.tajimaD = estimateTajimaD(segSites, theDiversityResults.totalSites, theDiversityResults.avgSiteCoverage,
                theAlignment.numberOfTaxa(), theDiversityResults.pipbp, theDiversityResults.thetapbp);
        return theDiversityResults;
    }

    public static double estimateTheta(int segSites, int totalSites, double averageSiteCoverage, int taxa) {
        double a = 0.0;
        double n = taxa * averageSiteCoverage / totalSites;
        for (double i = 1; i < n; i += 1.0) {
            a += 1 / i;
        }
        return segSites / a;
    }

    public static double estimateThetaPerbp(int segSites, int totalSites, double averageSiteCoverage, int taxa) {
        return estimateTheta(segSites, totalSites, averageSiteCoverage, taxa) / totalSites;
    }

    public static double estimatePi(int totalSites, double avgPairwiseDivergence, double averageSiteCoverage) {
        //this function is only needed to clearly show what we are doing with missing sites
        return totalSites * avgPairwiseDivergence / averageSiteCoverage;
    }

    public static double estimatePiPerbp(int totalSites, double avgPairwiseDivergence, double averageSiteCoverage) {
        //this function is only needed to clearly show what we are doing with missing sites
        return estimatePi(totalSites, avgPairwiseDivergence, averageSiteCoverage) / averageSiteCoverage;
    }

    public static double estimateTajimaD(int segSites, double totalSites, double averageSiteCoverage, double taxa, double pipbp, double thetapbp) {
        //this is the Tajima 1989 formula, but I don't know how it will deal with the
        //differential missing data
        double a1 = 0.0, a2 = 0.0;
        double n = taxa * averageSiteCoverage / totalSites;
        for (double i = 1; i < n; i += 1.0) {
            a1 += (1 / i);
            a2 += (1 / (i * i));
        }
        double b1 = (n + 1) / (3 * (n - 1));
        double b2 = (2 * (n * n + n + 3)) / (9 * n * (n - 1));
        double c1 = b1 - (1 / a1);
        double c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1));
        double e1 = c1 / a1;
        double e2 = c2 / (a1 * a1 + a2);
        double D = (pipbp - thetapbp) / (Math.sqrt((e1 * segSites) + (e2 * segSites * (segSites - 1))) / totalSites);
        return D;
    }

    int countSegregatingSites(GenotypeTable theAlignment) {
        int total = 0;
        if (theAlignment.isAllPolymorphic()) {
            return theAlignment.numberOfSites();
        }
        for (int i = 0; i < theAlignment.numberOfSites(); i++) {
            if (theAlignment.isPolymorphic(i)) {
                total++;
            }
        }
        return total;
    }

    //Implementation of TableReport Interface
    public Object[] getTableColumnNames() {
        String[] basicLabels = {"Site_Type", "Chromosome", "StartChrPosition", "EndChrPosition", "StartSite", "EndSite", "MidSite",
            "SiteCount", "AvgSiteCount", "SegSites", "PiPerBP",
            "ThetaPerBP", "Haplotypes", "TajimaD"};
        return basicLabels;
    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    @Override
    public Object[] getRow(long row) {

        Object[] data;
        java.text.NumberFormat nf = new java.text.DecimalFormat();
        nf.setMaximumFractionDigits(5);
        java.text.NumberFormat nf2 = new java.text.DecimalFormat();
        nf2.setMaximumFractionDigits(1);
        data = new String[NUM_OF_COLUMNS];
        DiversityResults theDiversityResults = diversityResultsVector.get((int) row);
        int labelOffset = 0;
        data[labelOffset++] = "ALL";
        data[labelOffset++] = "" + theDiversityResults.chromosome;
        data[labelOffset++] = "" + nf2.format(theDiversityResults.startChrPosition);
        data[labelOffset++] = "" + nf2.format(theDiversityResults.endChrPosition);
        data[labelOffset++] = "" + theDiversityResults.startSite;
        data[labelOffset++] = "" + theDiversityResults.endSite;
        data[labelOffset++] = "" + (theDiversityResults.startSite + theDiversityResults.endSite) / 2;
        data[labelOffset++] = "" + theDiversityResults.totalSites;
        data[labelOffset++] = "" + nf.format(theDiversityResults.avgSiteCoverage);
        data[labelOffset++] = "" + theDiversityResults.segregatingSites;
        data[labelOffset++] = "" + nf.format(theDiversityResults.pipbp);
        data[labelOffset++] = "" + nf.format(theDiversityResults.thetapbp);
        data[labelOffset++] = "NotAvail";//+theDiversityResults.haplotypes;
        data[labelOffset++] = "" + nf.format(theDiversityResults.tajimaD);
        return data;

    }

    @Override
    public String getTableTitle() {
        return "Diversity estimates";
    }

    @Override
    public String toString() {
        if (diversityResultsVector.size() == 0) {
            return "Needs to be run";
        }
        StringBuffer cs = new StringBuffer();
        Object[] header = this.getTableColumnNames();
        for (int i = 0; i < header.length; i++) {
            cs.append(header[i]);
        }
        cs.append("\n");
        for (int i = 0; i < getRowCount(); i++) {
            Object[] data = getRow(i);
            for (int j = 0; j < getColumnCount(); j++) {
                cs.append(data[j]);
            }
            cs.append("\n");
        }
        return cs.toString();
    }

    @Override
    public long getRowCount() {
        return diversityResultsVector.size();
    }

    @Override
    public long getElementCount() {
        return getRowCount() * getColumnCount();
    }

    @Override
    public int getColumnCount() {
        return NUM_OF_COLUMNS;
    }

    void testCode() {
        int n = 10, sites = 41, segSites = 16;
        double pi = 3.88;
        double theta = DiversityAnalyses.estimateTheta(segSites, sites, sites, n);
        double pipbp = pi / sites;
        double td = DiversityAnalyses.estimateTajimaD((int) segSites, sites, sites, n, pipbp, theta);

    }
}

class DiversityResults implements Serializable {

    protected double pipbp, thetapbp, totalSites, avgSiteCoverage, tajimaD, startChrPosition, endChrPosition;
    protected int startSite, endSite, haplotypes, segregatingSites;
    private int index;
    protected String chromosome;

    public DiversityResults(int start, int end, String chromosome,
            double startChrPosition, double endChrPosition) {
        this.startSite = start;
        this.endSite = end;
        this.chromosome = chromosome;
        this.startChrPosition = startChrPosition;
        this.endChrPosition = endChrPosition;
    }

    public boolean equals(Object anObject) {
        DiversityResults x = (DiversityResults) anObject;
        if (x.index == this.index) {
            return true;
        } else {
            return false;
        }
    }

    public void setIndex(int theIndex) {
        this.index = theIndex;
    }


    public String toString() {
        String result = "Diversity Results for " + "ALL" + "\n";
        result += "Pi =" + pipbp + "\n";
        result += "Theta =" + thetapbp + "\n";
        result += "Segregrating Sites =" + segregatingSites + "\n";
        result += "Total Sites =" + totalSites + "\n";
        result += "Start Site =" + startSite + "\n";
        result += "End Site =" + endSite + "\n";
        return result;
    }
}
