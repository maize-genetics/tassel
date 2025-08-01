package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.io.BufferedReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ImageIcon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.DifferenceGenotypeCallTable;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Utils;

public class SubtractGenotypesPlugin extends AbstractPlugin {
	private static final Logger myLogger = LogManager.getLogger(SubtractGenotypesPlugin.class);

    private PluginParameter<String> myHybridFile = new PluginParameter.Builder<>("hybridFile", null, String.class)
            .description("Three column tab-delimited file defining parent crosses. Columns are hybrid, firstParent, secondParent. "
            		+ "The second parent will be imputed from the hybrid by subtracting the first parent, which should be homozygous.")
            .required(true)
            .inFile()
            .build();

    public SubtractGenotypesPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> data = input.getDataOfType(GenotypeTable.class);
        if (data.size() != 1) {
            throw new IllegalArgumentException("SubtractGenotypesPlugin: preProcessParameters: must input 1 GenotypeTable.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        Datum data = input.getDataOfType(GenotypeTable.class).get(0);
        GenotypeTable genotypeTable = (GenotypeTable) data.getData();
        TaxaList origTaxa = genotypeTable.taxa();

        TaxaListBuilder taxaBuilder = new TaxaListBuilder();
        try (BufferedReader reader = Utils.getBufferedReader(myHybridFile.value())) {
            List<Integer> hybrids = new ArrayList<>();
            List<Integer> parents = new ArrayList<>();
            int lineNum = 1;
            String line = reader.readLine();
            while (line != null) {
                String[] tokens = line.trim().split("\t");
                if (tokens.length != 3) {
                    throw new IllegalArgumentException("SubtractGenotypePlugin: processData: a hybrid and two parents per line in file: "
                            + myHybridFile.value() + ". Problem on line: " + lineNum);
                }
                int firstIndex = origTaxa.indexOf(tokens[0]);
                int secondIndex = origTaxa.indexOf(tokens[1]);
                if (firstIndex == -1) {
                    myLogger.warn("processData: line: " + lineNum + " taxon name: " + tokens[0] + " not in the genotype.");
                } else if (secondIndex == -1) {
                    myLogger.warn("processData: line: " + lineNum + " taxon name: " + tokens[1] + " not in the genotype.");
                } else {
                	hybrids.add(firstIndex);
                	parents.add(secondIndex);
                    taxaBuilder.add(new Taxon(tokens[2].trim()));
                }
                line = reader.readLine();
                lineNum++;
            }

            GenotypeTable result = GenotypeTableBuilder.getInstance(new DifferenceGenotypeCallTable(genotypeTable.genotypeMatrix(), hybrids, parents), genotypeTable.positions(), taxaBuilder.build());
            return new DataSet(new Datum(data.getName() + "_ParentsBySubtraction", result, null), this);

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("SubtractGenotypesPlugin: processData: problem reading hybrid file: " + myHybridFile.value());
        }

    }

    /**
     * @return
     */
    public String hybridFile() {
        return myHybridFile.value();
    }

    /**
     * Sets the name of the hybrid file. 
     * @param filename
     * @return this plugin
     */
    public SubtractGenotypesPlugin hybridFile(String filename) {
        myHybridFile = new PluginParameter<>(myHybridFile, filename);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = SubtractGenotypesPlugin.class.getResource("/net/maizegenetics/analysis/images/hybrid.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Subtract Genotypes";
    }

    @Override
    public String getToolTipText() {
        return "Subtract Parent from Hybrid";
    }

}
