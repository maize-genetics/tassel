/*
 *  SortTaxaAlphabeticallyPlugin
 * 
 *  Created on Mar 2, 2017
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.net.URL;
import java.util.List;
import javax.swing.ImageIcon;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

/**
 *
 * @author Terry Casstevens
 */
public class SortTaxaAlphabeticallyPlugin extends AbstractPlugin {

    public SortTaxaAlphabeticallyPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        if (input == null) {
            throw new IllegalArgumentException("SortTaxaAlphabeticallyPlugin: processData: Please select one Genotype Table.");
        }
        List<Datum> genotypeTables = input.getDataOfType(GenotypeTable.class);
        if (genotypeTables.size() != 1) {
            throw new IllegalArgumentException("SortTaxaAlphabeticallyPlugin: processData: Please select one Genotype Table.");
        }
        GenotypeTable genotypes = (GenotypeTable) genotypeTables.get(0).getData();
        GenotypeTable result = FilterGenotypeTable.getInstanceSortTaxaAlphabetically(genotypes);

        if (result == genotypes) {
            DialogUtils.showWarning("Taxa already sorted alphabetically", getParentFrame());
            return input;
        } else {
            return new DataSet(new Datum(genotypeTables.get(0).getName() + "_Sorted_Taxa", result, null), this);
        }

    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = SortTaxaAlphabeticallyPlugin.class.getResource("/net/maizegenetics/analysis/images/sort.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Sort Taxa Alphabetically";
    }

    @Override
    public String getToolTipText() {
        return "Sort Taxa Alphabetically";
    }

}
