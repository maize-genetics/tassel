/*
 *  RemoveNaNFromDistanceMatrixPlugin
 *
 *  Created on Oct 5, 2015
 */
package net.maizegenetics.analysis.distance;

import com.google.common.collect.Range;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.DistanceMatrixBuilder;
import net.maizegenetics.taxa.distance.DistanceMatrixWithCounts;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Terry Casstevens
 */
public class RemoveNaNFromDistanceMatrixPlugin extends AbstractPlugin {

    private PluginParameter<Double> myMaxPercentNaN = new PluginParameter.Builder<>("maxPercentNaN", 0.0, Double.class)
            .description("Maximum percent of NaN allowed")
            .guiName("Maximum Percent NaN")
            .range(Range.closed(0.0, 100.0))
            .build();

    public RemoveNaNFromDistanceMatrixPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        List<Datum> temp = input.getDataOfType(DistanceMatrix.class);
        if (temp.size() != 1) {
            throw new IllegalArgumentException("RemoveNaNFromDistanceMatrixPlugin: processData: Must select one Distance Matrix");
        }

        DistanceMatrix origMatrix = (DistanceMatrix) temp.get(0).getData();
        String origName = temp.get(0).getName();

        TaxaList taxa = origMatrix.getTaxaList();
        int numTaxa = taxa.numberOfTaxa();
        BitSet[] whichNaN = new BitSet[numTaxa];

        for (int t = 0; t < numTaxa; t++) {
            whichNaN[t] = new OpenBitSet(numTaxa);
        }

        for (int x = 0; x < numTaxa; x++) {
            for (int y = x; y < numTaxa; y++) {
                if (Double.isNaN(origMatrix.getDistance(x, y))) {
                    whichNaN[x].fastSet(y);
                    whichNaN[y].fastSet(x);
                }
            }
        }

        long maxNaNToKeep = (long) Math.floor((double) numTaxa * maxPercentNaN());

        List<Integer> taxaToRemove = new ArrayList<>();
        boolean notFinished = true;
        while (notFinished) {
            long highestCount = 0;
            int highestTaxon = -1;
            for (int t = 0; t < numTaxa; t++) {
                long currentCount = whichNaN[t].cardinality();
                if (currentCount > highestCount) {
                    highestCount = currentCount;
                    highestTaxon = t;
                }
            }
            if (highestCount > maxNaNToKeep) {
                taxaToRemove.add(highestTaxon);
                for (int t = 0; t < numTaxa; t++) {
                    whichNaN[t].fastClear(highestTaxon);
                    whichNaN[highestTaxon].fastClear(t);
                }
            } else {
                notFinished = false;
            }
        }

        if (!taxaToRemove.isEmpty()) {

            int newNumTaxa = numTaxa - taxaToRemove.size();
            int[] origIndex = new int[newNumTaxa];
            TaxaListBuilder builder = new TaxaListBuilder();
            int count = 0;
            for (int t = 0; t < numTaxa; t++) {
                if (!taxaToRemove.contains(t)) {
                    builder.add(taxa.get(t));
                    origIndex[count++] = t;
                }
            }
            TaxaList newTaxaList = builder.build();

            DistanceMatrixBuilder result = DistanceMatrixBuilder.getInstance(newTaxaList);
            result.annotation(origMatrix.annotations());
            if (origMatrix instanceof DistanceMatrixWithCounts) {
                for (int x = 0; x < newNumTaxa; x++) {
                    for (int y = 0; y <= x; y++) {
                        result.set(x, y, origMatrix.getDistance(origIndex[x], origIndex[y]));
                        result.setCount(x, y, ((DistanceMatrixWithCounts) origMatrix).getCount(origIndex[x], origIndex[y]));
                    }
                }
            } else {
                for (int x = 0; x < newNumTaxa; x++) {
                    for (int y = 0; y <= x; y++) {
                        result.set(x, y, origMatrix.getDistance(origIndex[x], origIndex[y]));
                    }
                }
            }

            return new DataSet(new Datum(origName + " with " + maxPercentNaN()+ " percent NaN", result.build(), null), this);

        } else {
            return input;
        }

    }

    public static DistanceMatrix runPlugin(DistanceMatrix matrix) {
        return (DistanceMatrix) (new RemoveNaNFromDistanceMatrixPlugin(null, false).processData(DataSet.getDataSet(matrix))).getData(0).getData();
    }

    public static DataSet runPlugin(DataSet input) {
        return new RemoveNaNFromDistanceMatrixPlugin(null, false).processData(input);
    }

    /**
     * Maximum percent of NaN allowed
     *
     * @return Maximum Percent NaN
     */
    public Double maxPercentNaN() {
        return myMaxPercentNaN.value();
    }

    /**
     * Set Maximum Percent NaN. Maximum percent of NaN allowed
     *
     * @param value Maximum Percent NaN
     *
     * @return this plugin
     */
    public RemoveNaNFromDistanceMatrixPlugin maxPercentNaN(Double value) {
        myMaxPercentNaN = new PluginParameter<>(myMaxPercentNaN, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Remove NaN From Distance Matrix";
    }

    @Override
    public String getToolTipText() {
        return "Remove NaN From Distance Matrix";
    }

    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/RemoveNaN/RemoveNaN";
    }

}
