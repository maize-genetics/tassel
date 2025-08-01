/*
 * FilterSubsetPlugin.java
 *
 * Created on February 28, 2014 by Jason Wallace
 *
 */
package net.maizegenetics.analysis.filter;

import java.awt.Frame;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import javax.swing.*;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Jason Wallace
 */
public class FilterSubsetPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(FilterSubsetPlugin.class);
    private boolean myIsRandom = true;
    private double mySiteSubset = -1;
    private int[] mySitesToKeep = null;
    private double myTaxaSubset = -1;
    private int[] myTaxaToKeep = null;
    private TaxaList myTaxaList = null;

    /**
     * Creates a new instance of FilterSiteNamePlugin
     */
    public FilterSubsetPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        try {

            List inputData = input.getDataOfType(GenotypeTable.class);
            if (inputData.size() != 1) {
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection. Please select a single alignment.");
                } else {
                    myLogger.error("performFunction: Please input a single alignment.");
                }
                return null;
            }

            Datum td = processDatum((Datum) inputData.get(0), isInteractive());
            if (td == null) {
                return null;
            }

            DataSet output = new DataSet(td, this);

            fireDataSetReturned(new PluginEvent(output, FilterSubsetPlugin.class));

            return output;

        } finally {
            fireProgress(100);
        }
    }

    private Datum processDatum(Datum inDatum, boolean isInteractive) {

        final GenotypeTable alignment = (GenotypeTable) inDatum.getData();

        if (isInteractive) {
            myLogger.error("Interactive random subset not implemented yet.");
        }

        //Check that have something valid to filter with
        if (mySiteSubset <= 0 && myTaxaSubset <= 0) {
            myLogger.error("Must provide a positive value (decimal or integer) to use for subsetting");
            return null;
        }

        GenotypeTable result = alignment;
        //Filter on sites if necessary
        if (mySiteSubset > 0 && mySiteSubset < alignment.numberOfSites()) {
            determineSitesToKeep(alignment);
            result = FilterGenotypeTable.getInstance(result, mySitesToKeep);
        }
        //Filter on taxa if necessary
        if (myTaxaSubset > 0 && myTaxaSubset < alignment.numberOfTaxa()) {
            determineTaxaToKeep(alignment);
            result = FilterGenotypeTable.getInstance(result, myTaxaList);
        }

        //Return result
        String theName = inDatum.getName() + "_" + result.numberOfSites() + "_Sites";
        String theComment = "Subset of " + result.numberOfSites() + " from " + alignment.numberOfSites() + " Sites\n" + inDatum.getComment();
        return new Datum(theName, result, theComment);

    }

    private void determineSitesToKeep(GenotypeTable gt) {
        mySitesToKeep = getRandomSubset(mySiteSubset, gt.numberOfSites(), "sites"); //mySitesToKeep = sitesToKeep;
    }

    private void determineTaxaToKeep(GenotypeTable gt) {
        myTaxaToKeep = getRandomSubset(myTaxaSubset, gt.numberOfTaxa(), "taxa"); //mySitesToKeep = sitesToKeep;
        String[] taxa = new String[myTaxaToKeep.length];
        for (int i = 0; i < myTaxaToKeep.length; i++) {
            taxa[i] = gt.taxaName(myTaxaToKeep[i]);
        }
        TaxaListBuilder builder = new TaxaListBuilder();
        builder.addAll(taxa);
        myTaxaList = builder.build();
    }

    private int[] getRandomSubset(double subset_val, int max, String type) {
        //Process subset value. Values >=1 are made into the final count of sites/taxa. Values <1 are assumed to be a fraction
        int n;
        if (subset_val >= 1) {
            n = (int) Math.round(subset_val);
        } else {
            n = (int) Math.round(subset_val * max);
            if (n < 1) {
                myLogger.warn("Given values would result in <1 " + type + "; forcing minimum of 1");
                n = 1;
            }
        }

        //Check to see if n is more than half the size of max. If so, flag so that chosen sites are _excluded_, not included
        boolean flip = false;
        if (n > (double) (max) / 2.0) {
            flip = true;
            n = max - n;
        }

        //Generate indices to keep. Kept in HashMap for convenience while generating random numbers.
        HashMap<Integer, Boolean> chosen = new HashMap<>();
        if (myIsRandom) {	//Random subset (default)
            while (chosen.size() < n) {
                int index = (int) (Math.random() * max);
                if (index >= max) {	//Don't let go above maximum value
                    index = max - 1;
                }
                chosen.put(index, true);
            }
        } else {	//Stepwise subset
            double step = ((double) max) / n;
            double current = 0;
            while (current < max && chosen.size() < n) {
                chosen.put((int) Math.floor(current), true);
                current += step;
            }
        }

        //If specified subset was over half the available size, flip to _exclude_ the chosen sites
        if (flip) {
            HashMap<Integer, Boolean> flipped = new HashMap();
            for (int i = 0; i < max; i++) {
                if (!chosen.containsKey(i)) {
                    flipped.put(i, true);
                }
            }
            chosen = flipped;
            n = max - n;	//Reverse previous change to value of n (when first set "flip" to true
        }

        //Copy over chosen Integers to int[] array and then sort.
        int[] result = new int[chosen.size()];
        if (result.length != n) {
            myLogger.error("Output set for " + type + " has length " + result.length + " instead of specified length " + n);
        }
        int i = 0;
        for (Integer c : chosen.keySet()) {
            result[i] = c;
            i++;
        }
        Arrays.sort(result);

        return result;
    }

    public void setSiteSubset(double subset_val) {
        mySiteSubset = subset_val;
    }

    public void setTaxaSubset(double subset_val) {
        myTaxaSubset = subset_val;
    }

    public void setIsRandom(boolean isRandom) {
        myIsRandom = isRandom;
    }

    public double siteSubset() {
        return mySiteSubset;
    }

    public double taxaSubset() {
        return myTaxaSubset;
    }

    public int[] sitesToKeep() {
        return mySitesToKeep;
    }

    public int[] taxaToKeep() {
        return myTaxaToKeep;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
