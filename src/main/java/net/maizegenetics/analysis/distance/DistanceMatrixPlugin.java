/*
 * DistanceMatrixPlugin
 */
package net.maizegenetics.analysis.distance;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.ProgressListener;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author Terry Casstevens
 */
public class DistanceMatrixPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(DistanceMatrixPlugin.class);

    public DistanceMatrixPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
        if (alignInList.size() < 1) {
            throw new IllegalArgumentException("DistanceMatrixPlugin: performFunction: Please select a genotype table.");
        }

        List<DataSet> result = new ArrayList<>();
        Iterator<Datum> itr = alignInList.iterator();
        while (itr.hasNext()) {
            Datum current = itr.next();
            GenotypeTable aa = (GenotypeTable) current.getData();
            DistanceMatrix adm = getDistanceMatrix(aa, this);
            DataSet tds = new DataSet(new Datum("Matrix:" + current.getName(), adm, "Distance Matrix"), this);
            result.add(tds);
        }

        return DataSet.getDataSet(result, this);

    }

    public static DistanceMatrix getDistanceMatrix(GenotypeTable alignment) {
        return IBSDistanceMatrix3Alleles.getInstance(alignment);
    }

    public static DistanceMatrix getDistanceMatrix(GenotypeTable alignment, ProgressListener listener) {
        return IBSDistanceMatrix3Alleles.getInstance(alignment, listener);
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    @Override
    public ImageIcon getIcon() {
        URL imageURL = DistanceMatrixPlugin.class.getResource("/net/maizegenetics/analysis/images/DistanceMatrix.gif");
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
    @Override
    public String getButtonName() {
        return "Distance Matrix";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Create a distance matrix";
    }
}
