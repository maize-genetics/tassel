/*
 * HetsToUnknownPlugin.java
 *
 * Created on March 13, 2014
 *
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.net.URL;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.plugindef.Datum;

import javax.swing.*;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class HetsToUnknownPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(HetsToUnknownPlugin.class);

    public HetsToUnknownPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    /**
     * Returns Homozygous version of input Genotype Table.
     *
     * @param input Genotype Table.
     */
    public DataSet performFunction(DataSet input) {

        try {

            List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);

            if (alignInList.size() != 1) {
                String gpMessage = "Invalid selection.  Please select one genotype alignment.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), gpMessage);
                } else {
                    myLogger.error(gpMessage);
                }
                return null;
            }

            Datum current = alignInList.get(0);
            GenotypeTable genotypeTable = (GenotypeTable) current.getData();
            String name = current.getName();

            Datum outputDatum = new Datum(name + "_Homozygous", GenotypeTableBuilder.getHomozygousInstance(genotypeTable), "Homozygous " + name);
            DataSet output = new DataSet(outputDatum, this);
            fireDataSetReturned(new PluginEvent(output, HetsToUnknownPlugin.class));

            return output;

        } finally {
            fireProgress(100);
        }

    }

    public String getToolTipText() {
        return "Change Heterozygous to Unknown";
    }

    public ImageIcon getIcon() {
        URL imageURL = HetsToUnknownPlugin.class.getResource("/net/maizegenetics/analysis/images/homozygous.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    public String getButtonName() {
        return "Homozygous Genotype";
    }

}
