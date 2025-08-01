/*
 *  GetPositionListPlugin
 * 
 *  Created on Jul 7, 2014
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;

import java.net.URL;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.Datum;

import javax.swing.*;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class GetPositionListPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(GetPositionListPlugin.class);

    public GetPositionListPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        try {

            List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);

            if (alignInList.size() != 1) {
                throw new IllegalArgumentException("Invalid selection.  Please select one genotype table.");
            }

            Datum current = alignInList.get(0);
            GenotypeTable genotypeTable = (GenotypeTable) current.getData();
            String name = current.getName();

            Datum outputDatum = new Datum(name + "_PositionList", genotypeTable.positions(), "Position List from " + name);
            DataSet output = new DataSet(outputDatum, this);

            return output;

        } finally {
            fireProgress(100);
        }

    }

    @Override
    public String getToolTipText() {
        return "Get Position List";
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = GetPositionListPlugin.class.getResource("/net/maizegenetics/analysis/images/lists.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Get Position List";
    }

}
