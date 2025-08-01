package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.List;

/**
 * @author Terry Casstevens Created June 01, 2017
 */
public class IndelsToUnknownPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(HetsToUnknownPlugin.class);

    public IndelsToUnknownPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);

        if (alignInList.size() != 1) {
            throw new IllegalArgumentException("IndelsToUnknownPlugin: processData: Please select one genotype table.");
        }

        Datum current = alignInList.get(0);
        GenotypeTable genotypeTable = (GenotypeTable) current.getData();
        String name = current.getName();

        Datum result = new Datum(name + "_NoIndels", GenotypeTableBuilder.getInstanceMaskIndels(genotypeTable), "Indels changed to Unknown " + name);
        return new DataSet(result, this);

    }

    @Override
    public String getToolTipText() {
        return "Change Indels to Unknown";
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = HetsToUnknownPlugin.class.getResource("/net/maizegenetics/analysis/images/homozygous.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    public String getButtonName() {
        return "Change Indels to Unknown";
    }

}

