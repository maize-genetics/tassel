/*
 * MergeRenameDeleteTaxaPlugin
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.genotypecall.BasicGenotypeMergeRule;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.List;

/**
 * @author jcg233
 */
public class MergeRenameDeleteTaxaPlugin extends net.maizegenetics.plugindef.AbstractPlugin {

    @Override
    public String pluginDescription() {
        return "Rename taxa, merging those with the same new name in the taxa rename key file. "
                + "Taxa with the new name \"delete\" will be removed.";
    }

    private PluginParameter<String> taxaRenameKey = new PluginParameter.Builder<>("renameKey", null, String.class)
            .required(true)
            .inFile()
            .guiName("Taxa Rename Key")
            .description("Tab-delimited file with original and new taxa names. Taxa with the same new name will be merged. "
                    + "Taxa with the new name \"delete\" will be removed. Any other columns (and the header line) are ignored.")
            .build();
    private PluginParameter<String> outputHDF5Genotypes = new PluginParameter.Builder<>("o", null, String.class)
            .guiName("Output HDF Genotypes")
            .required(true)
            .outFile()
            .description("Output HDF5 genotypes file").build();
    private PluginParameter<Double> avgSeqErrorRate = new PluginParameter.Builder<>("eR", 0.01, Double.class)
            .guiName("Avg Seq Error Rate")
            .description("Average sequencing error rate per base (used to decide between heterozygous and homozygous calls when merging taxa)")
            .build();
    private PluginParameter<Boolean> noDepthOutput = new PluginParameter.Builder<>("ndo", false, Boolean.class)
            .guiName("No Depth Output")
            .description("No depth output: do not write depths to the output HDF5 genotypes file")
            .build();
    private PluginParameter<String> dataSetName = new PluginParameter.Builder<>("name", null, String.class)
            .guiName("Data set name")
            .required(false)
            .description("(Optional) Short data set name to be added as an root level annotation under \"dataSetName\"")
            .build();
    private PluginParameter<String> dataSetDescription = new PluginParameter.Builder<>("desc", null, String.class)
            .guiName("Data set description")
            .required(false)
            .description("(Optional) Short data set description to be added as an root level annotation under \"dataSetDescription\"")
            .build();

    public MergeRenameDeleteTaxaPlugin() {
        super(null, false);
    }

    public MergeRenameDeleteTaxaPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    private static final Logger myLogger = LogManager.getLogger(MergeRenameDeleteTaxaPlugin.class);
    String dataSetDescrip, date;
    String errorMessage;
    private GenotypeTable inputGenotypes = null;
    private String inputGenosName = null;
    private TreeMap<String, TreeSet<String>> newNameToOldNames = new TreeMap();
    private BasicGenotypeMergeRule genoMergeRule = null;
    private GenotypeTableBuilder genos = null; //output genotype table
    private TaxaList taxaList = null;

    @Override
    protected void preProcessParameters(DataSet input) {
        if (input == null) {
            throw new IllegalArgumentException("MergeRenameDeleteTaxaPlugin: preProcessParameters: Please select one Genotype Table.");
        }
        List<Datum> genotypeTables = input.getDataOfType(GenotypeTable.class);

        inputGenosName = genotypeTables.get(0).getName();

        myLogger.info("\n" + pluginDescription() + "\n");

        myLogger.info("Input genotype name: " + inputGenosName);

        if (genotypeTables.size() == 1) {
            inputGenotypes
                    = (GenotypeTable) genotypeTables.get(0).getData();
            taxaList = inputGenotypes.taxa();
        } else {
            throw new IllegalArgumentException("MergeRenameDeleteTaxaPlugin: preProcessParameters: Please select one Genotype Table.");
        }
        date = new SimpleDateFormat("yyyyMMdd").format(new Date());
        String outfile = outputHDF5Genotypes();
        outputHDF5Genotypes(outfile.replace("__DATE__", "_" + date));
    }

    public DataSet processData(DataSet input) {
        readTaxaRenameKey();
        setUpHDF5GenotypeTableBuilder();
        int nTaxaAdded = addRenamedMergedTaxa();
        if (dataSetName() != null) {
            genos.dataSetName(parseDataSetName(dataSetName()));
        }
        if (dataSetDescription() != null) {
            genos.dataSetDescription(parseDataSetDescription(dataSetDescription(), nTaxaAdded));
        }
        genos.build();
        myLogger.info("\n\nFinished creating new HDF5 genotpye file with merged and renamed taxa.\n\n");
        return null;
    }

    private void readTaxaRenameKey() {
        myLogger.info("\nReading the taxaRenameKey file:\n   " + taxaRenameKey() + "\n");
        BufferedReader taxaRenameKeyReader = Utils.getBufferedReader(taxaRenameKey());
        String line;
        int nLinesRead = 0;
        try {
            while ((line = taxaRenameKeyReader.readLine()) != null) {
                nLinesRead++;
                if (nLinesRead == 1) {
                    continue;  // skip the header
                }
                String[] values = line.split("\t", -1);
                String oldName = values[0];
                String newName = values[1];
                if (!newNameToOldNames.containsKey(newName)) {
                    newNameToOldNames.put(newName, new TreeSet());
                }
                newNameToOldNames.get(newName).add(oldName);
            }
        } catch (IOException e) {
            System.err.println("\n\nProblem reading the taxaRenameKey file (" + taxaRenameKey() + "):\n\t" + e);
            System.exit(1);
        }
        myLogger.info("\nFinished reading the taxaRenameKey file (nTaxa=" + (nLinesRead - 1) + ")\n");
    }


    private void setUpHDF5GenotypeTableBuilder() {
        genoMergeRule = new BasicGenotypeMergeRule(avgSeqErrorRate());
        File hdf5File = new File(outputHDF5Genotypes());
        if (hdf5File.exists()) {
            errorMessage = "\nERROR: the output HDF5 genotypes file:\n   " + outputHDF5Genotypes() + "\n already exists\n\n";
            myLogger.error(errorMessage);
            throw new IllegalStateException(errorMessage);
        } else {
            myLogger.info("\nInitializing the output HDF5 file:\n   " + outputHDF5Genotypes() + "\n\n");
            genos = GenotypeTableBuilder.getTaxaIncrementalWithMerging(outputHDF5Genotypes(), inputGenotypes.positions(), genoMergeRule);
        }
    }

    private int addRenamedMergedTaxa() {
        int nTaxaAdded = 0;
        int nTaxaDeleted = 0;
        for (Map.Entry<String, TreeSet<String>> newNameAndOldNames : newNameToOldNames.entrySet()) {
            String newName = newNameAndOldNames.getKey();
            int nTaxaToMergeOrDelete = newNameAndOldNames.getValue().size();
            if (newName.equals("delete") || newName.equals("remove")) {
                nTaxaDeleted += nTaxaToMergeOrDelete;
                continue;
            }
            StringBuilder oldNames = new StringBuilder();
            Taxon.Builder TaxonBuilder = null;
            int[][] alleleDepths = new int[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][inputGenotypes.positions().numberOfSites()];
            byte[] taxonGenos = null;
            for (String oldName : newNameAndOldNames.getValue()) {
                oldNames.append(oldName + ",");
                Taxon oldTaxon = taxaList.get(taxaList.indexOf(oldName));
                if (TaxonBuilder == null) {
                    TaxonBuilder = new Taxon.Builder(oldTaxon);  // adds all the annotations from the first oldName taxon
                    TaxonBuilder = TaxonBuilder.name(newName);
                } else {
                    GeneralAnnotation oldAnnos = oldTaxon.getAnnotation();
                    for (Map.Entry<String, String> oldAnno : oldAnnos.getAllAnnotationEntries()) {
                        TaxonBuilder = TaxonBuilder.addAnno(oldAnno.getKey(), oldAnno.getValue());
                    }
                }
                TaxonBuilder = TaxonBuilder.addAnno("OldName", oldName);

                if (inputGenotypes.hasDepth()) {
                    for (int site = 0; site < inputGenotypes.positions().numberOfSites(); site++) {
                        int[] alleleDepthsAtSite = inputGenotypes.depthForAlleles(taxaList.indexOf(oldName), site);
                        for (int allele = 0; allele < NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES; allele++) {
                            alleleDepths[allele][site] += alleleDepthsAtSite[allele];
                        }
                    }
                }
                if (nTaxaToMergeOrDelete == 1) {
                    taxonGenos = inputGenotypes.genotypeAllSites(taxaList.indexOf(oldName));
                }
            }
            if (nTaxaToMergeOrDelete > 1) {
                if (inputGenotypes.hasDepth()) {
                    taxonGenos = resolveGenosForTaxon(alleleDepths);
                } else {
                    throw new IllegalStateException("\n\nERROR: Merging genotypes across replicate taxa is not allowed when there is no depth\n\n");
                }
            }
            if (noDepthOutput()) {
                genos.addTaxon(TaxonBuilder.build(), taxonGenos, null);
            } else {
                genos.addTaxon(TaxonBuilder.build(), alleleDepths, taxonGenos);
            }
            myLogger.info("  ...finished calling/adding genotypes for " + newName + "   OldName(s):" + oldNames.deleteCharAt(oldNames.length() - 1).toString());
            nTaxaAdded++;
        }
        myLogger.info("\nFinished adding genotypes for " + nTaxaAdded + " taxa.  nTaxaDeleted=" + nTaxaDeleted);
        return nTaxaAdded;
    }

    private byte[] resolveGenosForTaxon(int[][] depthsForTaxon) {
        int nAlleles = depthsForTaxon.length;
        int[] depthsAtSite = new int[nAlleles];
        int nSites = depthsForTaxon[0].length;
        byte[] genos = new byte[nSites];
        for (int site = 0; site < nSites; site++) {
            for (int allele = 0; allele < nAlleles; allele++) {
                depthsAtSite[allele] = depthsForTaxon[allele][site];
            }
            genos[site] = genoMergeRule.callBasedOnDepth(depthsAtSite);
        }
        return genos;
    }

    private String parseDataSetName(String dataSetName) {
        return dataSetName.replace("__DATE__", "_" + date);
    }

    private String parseDataSetDescription(String dataSetDescrip, int nTaxa) {
        int nSNPs = inputGenotypes.numberOfSites();
        return dataSetDescrip.replace("__SNPS__", "" + nSNPs).replace("__TAXA__", "" + nTaxa).replace("__DATE__", date);
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = MergeRenameDeleteTaxaPlugin.class
                .getResource("/net/maizegenetics/analysis/images/lowDepthToMissing.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "MergeRenameTaxa";
    }

    @Override
    public String getToolTipText() {
        return "Rename taxa, merging those with the same new name";
    }

    /**
     * Convenience method to run plugin with input and output GenotypeTable objects
     * (rather than DataSets)
     *
     * @param inputGenos Input GenotypeTable.
     *
     * @return GenotypeTable where genotypes with depth below the specified minimum are set to missing.
     */
    public void runPlugin(GenotypeTable inputGenos) {
        DataSet input = new DataSet(new Datum("inputGenotypes", inputGenos, null), null);
        runPlugin(input);
    }

    /**
     * Convenience method to run plugin.
     */
    public void runPlugin(DataSet input) {
        performFunction(input);
    }

    /**
     * Tab-delimited file with original and new taxa names.
     * Taxa with the same new name will be merged. Taxa with
     * the new name "delete" will be removed. Any other columns
     * (and the header line) are ignored.
     *
     * @return Taxa Rename Key
     */
    public String taxaRenameKey() {
        return taxaRenameKey.value();
    }

    /**
     * Set Taxa Rename Key. Tab-delimited file with original
     * and new taxa names. Taxa with the same new name will
     * be merged. Taxa with the new name "delete" will be
     * removed. Any other columns (and the header line) are
     * ignored.
     *
     * @param value Taxa Rename Key
     *
     * @return this plugin
     */
    public MergeRenameDeleteTaxaPlugin taxaRenameKey(String value) {
        taxaRenameKey = new PluginParameter<>(taxaRenameKey, value);
        return this;
    }

    /**
     * Output genotypes file
     *
     * @return Output Genotypes
     */
    public String outputHDF5Genotypes() {
        return outputHDF5Genotypes.value();
    }

    /**
     * Set Output Genotypes. Output genotypes file
     *
     * @param value Output Genotypes
     *
     * @return this plugin
     */
    public MergeRenameDeleteTaxaPlugin outputHDF5Genotypes(String value) {
        outputHDF5Genotypes = new PluginParameter<>(outputHDF5Genotypes, value);
        return this;
    }

    /**
     * Average sequencing error rate per base (used to decide
     * between heterozygous and homozygous calls when merging
     * taxa)
     *
     * @return Avg Seq Error Rate
     */
    public Double avgSeqErrorRate() {
        return avgSeqErrorRate.value();
    }

    /**
     * Set Avg Seq Error Rate. Average sequencing error rate
     * per base (used to decide between heterozygous and homozygous
     * calls when merging taxa)
     *
     * @param value Avg Seq Error Rate
     *
     * @return this plugin
     */
    public MergeRenameDeleteTaxaPlugin avgSeqErrorRate(Double value) {
        avgSeqErrorRate = new PluginParameter<>(avgSeqErrorRate, value);
        return this;
    }

    /**
     * No depth output: do not write depths to the output
     * genotypes file (applies only to hdf5 or VCF)
     *
     * @return No Depth Output
     */
    public Boolean noDepthOutput() {
        return noDepthOutput.value();
    }

    /**
     * Set No Depth Output. No depth output: do not write
     * depths to the output genotypes file (applies only to
     * hdf5 or VCF)
     *
     * @param value No Depth Output
     *
     * @return this plugin
     */
    public MergeRenameDeleteTaxaPlugin noDepthOutput(Boolean value) {
        noDepthOutput = new PluginParameter<>(noDepthOutput, value);
        return this;
    }

    /**
     * (Optional) Short data set name to be added as an root
     * level annotation under "/DataSetName"
     *
     * @return Data set name
     */
    public String dataSetName() {
        return dataSetName.value();
    }

    /**
     * Set Data set name. (Optional) Short data set name to
     * be added as an root level annotation under "dataSetName"
     *
     * @param value Data set name
     *
     * @return this plugin
     */
    public MergeRenameDeleteTaxaPlugin dataSetName(String value) {
        dataSetName = new PluginParameter<>(dataSetName, value);
        return this;
    }

    /**
     * (Optional) Short data set description to be added as
     * an root level annotation under "dataSetDescription"
     *
     * @return Data set description
     */
    public String dataSetDescription() {
        return dataSetDescription.value();
    }

    /**
     * Set Data set description. (Optional) Short data set
     * description to be added as an root level annotation
     * under "/DataSetDescription"
     *
     * @param value Data set description
     *
     * @return this plugin
     */
    public MergeRenameDeleteTaxaPlugin dataSetDescription(String value) {
        dataSetDescription = new PluginParameter<>(dataSetDescription, value);
        return this;
    }

}
