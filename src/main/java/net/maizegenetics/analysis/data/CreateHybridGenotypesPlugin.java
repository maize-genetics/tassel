/*
 *  CreateHybridGenotypesPlugin
 * 
 *  Created on May 17, 2016
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.HybridGenotypeCallTable;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.dna.snp.score.ReferenceProbabilityBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Terry Casstevens
 *
 * This plugin creates a hybrid genotype table from a reference genotype table and a hybrid file
 * that defines the parents. The hybrid genotype is a combination of two homozyous values or unknown
 * if either is heterozygous. For the numeric option, the numeric genotype values are averaged.
 */
public class CreateHybridGenotypesPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(CreateHybridGenotypesPlugin.class);

    private PluginParameter<String> myHybridFile = new PluginParameter.Builder<>("hybridFile", null, String.class)
            .description("Two column tab-delimited file defining parent crosses.")
            .required(true)
            .inFile()
            .build();

    private PluginParameter<String> myHybridChar = new PluginParameter.Builder<>("hybridChar", "/", String.class)
            .description("String used to combine taxa names to create hybrid name.")
            .build();

    private PluginParameter<Boolean> myNumeric = new PluginParameter.Builder<>("numericGenotype", false, Boolean.class)
            .description("Whether to create a numeric hybrid genotype.")
            .build();

    public CreateHybridGenotypesPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> data = input.getDataOfType(GenotypeTable.class);
        if (data.size() != 1) {
            throw new IllegalArgumentException("CreateHybridGenotypesPlugin: preProcessParameters: must input 1 GenotypeTable.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        Datum data = input.getDataOfType(GenotypeTable.class).get(0);
        GenotypeTable genotypeTable = (GenotypeTable) data.getData();
        TaxaList origTaxa = genotypeTable.taxa();

        TaxaListBuilder taxa = new TaxaListBuilder();
        try (BufferedReader reader = Utils.getBufferedReader(hybridFile())) {
            List<Integer> firstParents = new ArrayList<>();
            List<Integer> secondParents = new ArrayList<>();
            int lineNum = 1;
            String line = reader.readLine();
            while (line != null) {
                String[] tokens = line.trim().split("\t");
                if (tokens.length != 2) {
                    throw new IllegalArgumentException("CreateHybridGenotypePlugin: processData: Must have two parents per line in file: "
                            + hybridFile() + ". Problem on line: " + lineNum);
                }
                int firstIndex = origTaxa.indexOf(tokens[0]);
                int secondIndex = origTaxa.indexOf(tokens[1]);
                if (firstIndex == -1) {
                    myLogger.warn("processData: line: " + lineNum + " taxon name: " + tokens[0] + " not in the genotype.");
                } else if (secondIndex == -1) {
                    myLogger.warn("processData: line: " + lineNum + " taxon name: " + tokens[1] + " not in the genotype.");
                } else {
                    firstParents.add(firstIndex);
                    secondParents.add(secondIndex);
                    taxa.add(new Taxon(tokens[0] + hybridChar() + tokens[1]));
                }
                line = reader.readLine();
                lineNum++;
            }

            if (numeric()) {
                return new DataSet(new Datum(data.getName() + "_Numeric_Hybrids", numericHybridGenotype(genotypeTable, firstParents, secondParents, taxa.build()), null), this);
            } else {
                GenotypeTable result = GenotypeTableBuilder.getInstance(new HybridGenotypeCallTable(genotypeTable.genotypeMatrix(), firstParents, secondParents), genotypeTable.positions(), taxa.build());
                return new DataSet(new Datum(data.getName() + "_Hybrids", result, null), this);
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("CreateHybridGenotypePlugin: processData: problem reading hybrid file: " + hybridFile());
        }

    }

    /**
     * Create a numeric hybrid genotype. The hybrid genotype is the average of the
     * two parent's numeric genotype values.
     */
    private static GenotypeTable numericHybridGenotype(GenotypeTable table, List<Integer> firstParents, List<Integer> secondParents, TaxaList taxa) {

        int nsites = table.numberOfSites();
        int ntaxa = taxa.numberOfTaxa();

        ReferenceProbability original = table.referenceProbability();

        // build new ReferenceProbability
        ReferenceProbabilityBuilder refBuilder = ReferenceProbabilityBuilder.getInstance(ntaxa, nsites, table.taxa());
        for (int t = 0; t < ntaxa; t++) {
            float[] hybrid = new float[nsites];
            for (int s = 0; s < nsites; s++) {
                hybrid[s] = (original.value(firstParents.get(t), s) + original.value(secondParents.get(t), s)) / 2.0f;
            }
            refBuilder.addTaxon(t, hybrid);
        }

        // build new GenotypeTable
        return GenotypeTableBuilder.getInstance(null, table.positions(),
                taxa, null, null, refBuilder.build(), null,
                null);

    }

    /**
     * Hybrid File
     *
     * @return Hybrid File
     */
    public String hybridFile() {
        return myHybridFile.value();
    }

    /**
     * Set Hybrid File. Hybrid File
     *
     * @param value Hybrid File
     *
     * @return this plugin
     */
    public CreateHybridGenotypesPlugin hybridFile(String value) {
        myHybridFile = new PluginParameter<>(myHybridFile, value);
        return this;
    }

    /**
     * String used to combine taxa names to create hybrid name.
     *
     * @return Cross Char
     */
    public String hybridChar() {
        return myHybridChar.value();
    }

    /**
     * Set Cross Char. String used to combine taxa names to create hybrid name.
     *
     * @param value Cross Char
     *
     * @return this plugin
     */
    public CreateHybridGenotypesPlugin hybridChar(String value) {
        myHybridChar = new PluginParameter<>(myHybridChar, value);
        return this;
    }

    /**
     * Whether to create a numeric hybrid genotype.
     *
     * @return Numeric Genotype
     */
    public Boolean numeric() {
        return myNumeric.value();
    }

    /**
     * Set Numeric Genotype. Whether to create a numeric hybrid
     * genotype.
     *
     * @param value Numeric Genotype
     *
     * @return this plugin
     */
    public CreateHybridGenotypesPlugin numeric(Boolean value) {
        myNumeric = new PluginParameter<>(myNumeric, value);
        return this;
    }

    @Override
    public String getToolTipText() {
        return "Create Hybrid Genotypes";
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = CreateHybridGenotypesPlugin.class.getResource("/net/maizegenetics/analysis/images/hybrid.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Create Hybrid Genotypes";
    }

    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/CreateHybridGenotypes/CreateHybridGenotypes";
    }

}
