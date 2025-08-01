/*
 * FilterTaxaPropertiesPlugin.java
 *
 * Created on May 9, 2013
 *
 */
package net.maizegenetics.analysis.filter;

import com.google.common.collect.Range;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.plugindef.PluginParameter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.List;

/**
 *
 * @author Terry Casstevens
 */
public class FilterTaxaPropertiesPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(FilterTaxaPropertiesPlugin.class);

    private PluginParameter<Double> myMinNotMissing = new PluginParameter.Builder<Double>("minNotMissing", TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING_DEFAULT, Double.class)
            .guiName("Min Proportion of Sites Present").range(Range.closed(0.0, 1.0)).build();
    private PluginParameter<Double> myMinHeterozygous = new PluginParameter.Builder<Double>("minHeterozygous", TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MIN_HET_DEFAULT, Double.class)
            .guiName("Min Heterozygous Proportion").range(Range.closed(0.0, 1.0)).build();
    private PluginParameter<Double> myMaxHeterozygous = new PluginParameter.Builder<Double>("maxHeterozygous", TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MAX_HET_DEFAULT, Double.class)
            .guiName("Max Heterozygous Proportion").range(Range.closed(0.0, 1.0)).build();

    /**
     * Creates a new instance of FilterTaxaPropertiesPlugin
     */
    public FilterTaxaPropertiesPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {

        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);

        if (alignInList.size() != 1) {
            throw new IllegalArgumentException("FilterTaxaPropertiesPlugin: preProcessParameters: Must select one Genotype Table.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        List<Datum> datumList = input.getDataOfType(GenotypeTable.class);
        Datum inDatum = datumList.get(0);
        GenotypeTable aa = (GenotypeTable) inDatum.getData();

        GenotypeTable result = getFilteredAlignment(aa);

        StringBuilder builder = new StringBuilder();
        builder.append("Filter Alignment by Taxa Properties...\n");
        builder.append("   Min. Proportion of Sites Present: ");
        builder.append(minProportionOfSitesPresent());
        builder.append("\n");
        builder.append("   Heterozygous Proportion: ");
        builder.append(minHeterozygousProportion());
        builder.append(" - ");
        builder.append(maxHeterozygousProportion());
        builder.append("\n");
        String theComment = builder.toString();

        if (result == aa) {
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), "The Alignment is Unchanged.");
            } else {
                myLogger.warn("The Alignment is Unchanged: " + inDatum.getName());
            }
        }

        if (result.numberOfTaxa() != 0) {
            String theName = inDatum.getName() + "_" + result.numberOfTaxa() + "Taxa";
            myLogger.info("Resulting Number Sequences: " + result.numberOfTaxa());
            return new DataSet(new Datum(theName, result, theComment), this);
        } else {
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), "No remaining Taxa given filter parameters.");
            } else {
                myLogger.warn("No remaining Taxa given filter parameters.");

            }
            return null;
        }
    }

    private GenotypeTable getFilteredAlignment(GenotypeTable alignment) {
        int numSites = alignment.numberOfSites();
        int numTaxa = alignment.numberOfTaxa();
        TaxaList ids = alignment.taxa();

        TaxaListBuilder keepTaxaList = new TaxaListBuilder();
        for (int t = 0; t < numTaxa; t++) {

            progress((int) ((double) t / (double) numTaxa * 100.0), null);

            if (minProportionOfSitesPresent() != 0.0) {
                int totalNotMissing = alignment.totalNonMissingForTaxon(t);
                double percentNotMissing = (double) totalNotMissing / (double) numSites;
                if (percentNotMissing < minProportionOfSitesPresent()) {
                    continue;
                }
            }

            if ((minHeterozygousProportion() != 0.0) || (maxHeterozygousProportion() != 1.0)) {
                int numHeterozygous = alignment.heterozygousCountForTaxon(t);
                int totalSitesNotMissing = alignment.totalNonMissingForTaxon(t);
                double percentHets = (double) numHeterozygous / (double) totalSitesNotMissing;
                if ((percentHets < minHeterozygousProportion()) || (percentHets > maxHeterozygousProportion())) {
                    continue;
                }
            }

            keepTaxaList.add(ids.get(t));
        }
        return FilterGenotypeTable.getInstance(alignment, keepTaxaList.build(), false);
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(FilterTaxaPropertiesPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Min Proportion of Sites Present
     *
     * @return Min Proportion of Sites Present
     */
    public Double minProportionOfSitesPresent() {
        return myMinNotMissing.value();
    }

    /**
     * Set Min Proportion of Sites Present. Min Proportion of Sites Present
     *
     * @param value Min Proportion of Sites Present
     *
     * @return this plugin
     */
    public FilterTaxaPropertiesPlugin minProportionOfSitesPresent(Double value) {
        myMinNotMissing = new PluginParameter<>(myMinNotMissing, value);
        return this;
    }

    /**
     * Min Heterozygous Proportion
     *
     * @return Min Heterozygous Proportion
     */
    public Double minHeterozygousProportion() {
        return myMinHeterozygous.value();
    }

    /**
     * Set Min Heterozygous Proportion. Min Heterozygous Proportion
     *
     * @param value Min Heterozygous Proportion
     *
     * @return this plugin
     */
    public FilterTaxaPropertiesPlugin minHeterozygousProportion(Double value) {
        myMinHeterozygous = new PluginParameter<>(myMinHeterozygous, value);
        return this;
    }

    /**
     * Max Heterozygous Proportion
     *
     * @return Max Heterozygous Proportion
     */
    public Double maxHeterozygousProportion() {
        return myMaxHeterozygous.value();
    }

    /**
     * Set Max Heterozygous Proportion. Max Heterozygous Proportion
     *
     * @param value Max Heterozygous Proportion
     *
     * @return this plugin
     */
    public FilterTaxaPropertiesPlugin maxHeterozygousProportion(Double value) {
        myMaxHeterozygous = new PluginParameter<>(myMaxHeterozygous, value);
        return this;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = FilterTaxaPropertiesPlugin.class.getResource("/net/maizegenetics/analysis/images/Filter_horizontal.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Taxa";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Filter Alignment Based Taxa Properties";
    }
}
