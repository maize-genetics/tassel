/*
 * GetTaxaListPlugin.java
 *
 * Created on June 16, 2014
 *
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Terry Casstevens
 */
public class GetTaxaListPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(GetTaxaListPlugin.class);

    public GetTaxaListPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet processData(DataSet input) {

        List<Datum> result = new ArrayList<>();

        List<Datum> datumList = input.getDataOfType(GenotypeTable.class);
        for (Datum current : datumList) {
            String name = current.getName();
            GenotypeTable genotypeTable = (GenotypeTable) current.getData();
            result.add(new Datum(name + "_TaxaList", genotypeTable.taxa(), "Taxa List from " + name));
        }

        datumList = input.getDataOfType(DistanceMatrix.class);
        for (Datum current : datumList) {
            String name = current.getName();
            DistanceMatrix matrix = (DistanceMatrix) current.getData();
            result.add(new Datum(name + "_TaxaList", matrix.getTaxaList(), "Taxa List from " + name));
        }

        datumList = input.getDataOfType(Phenotype.class);
        for (Datum current : datumList) {
            String name = current.getName();
            Phenotype phenotype = (Phenotype) current.getData();
            result.add(new Datum(name + "_TaxaList", phenotype.taxa(), "Taxa List from " + name));
        }

        if (result.isEmpty()) {
            throw new IllegalArgumentException("GetTaxaListPlugin: processData: nothing is selected that has a taxa list.");
        } else {
            return new DataSet(result, this);
        }

    }

    public String getToolTipText() {
        return "Get Taxa List";
    }

    public ImageIcon getIcon() {
        URL imageURL = GetTaxaListPlugin.class.getResource("/net/maizegenetics/analysis/images/lists.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    public String getButtonName() {
        return "Get Taxa List";
    }

}
