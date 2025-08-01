package net.maizegenetics.analysis.data;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.map.TOPMInterface;
import net.maizegenetics.dna.map.TOPMUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Created by jgw87 on 5/28/14.
 * This plugin is meant to take an HDF5 file and output a report of various summary statistics or information.
 *  INclude
 */
//TODO: Add support for chromosomes, positions, etc.
public class HDF5SummaryPlugin extends AbstractPlugin {

    private Logger logger = LogManager.getLogger(HDF5SummaryPlugin.class);
    private enum H5FileType {GENOTYPE, TOPM, TBT, UNKNOWN}
    private H5FileType myFileType = H5FileType.UNKNOWN;

    //Various handles to be used by the class functions
    GenotypeTable genos = null;
    TagsByTaxa tbt = null;
    TOPMInterface topm = null;
    IHDF5Reader h5reader = null;
    BufferedWriter outputWriter = null;


    private PluginParameter<String> inputFile
            = new PluginParameter.Builder<>("input", null, String.class)
            .description("TASSEL HDF5 file to get summary data from")
            .guiName("HDF5 file")
            .required(true)
            .inFile()
            .build();
    private PluginParameter<String> outputFile
            = new PluginParameter.Builder<>("output", null, String.class)
            .description("File to write summary data to")
            .guiName("Output file")
            .required(true)
            .outFile()
            .build();
    private PluginParameter<Boolean> taxaCount
            = new PluginParameter.Builder<>("taxaCount", false, Boolean.class)
            .description("Output number of taxa in file")
            .guiName("Output count of taxa?")
            .build();
    private PluginParameter<Boolean> taxaNames
            = new PluginParameter.Builder<>("taxaNames", false, Boolean.class)
            .description("Output names of all taxa in file")
            .guiName("Output taxa names")
            .build();
    private PluginParameter<Boolean> siteCount
            = new PluginParameter.Builder<>("siteCount", false, Boolean.class)
            .description("Output number of sites in file")
            .guiName("Output site count")
            .build();
    private PluginParameter<Boolean> siteNames
            = new PluginParameter.Builder<>("siteNames", false, Boolean.class)
            .description("Output names of all sites in file")
            .guiName("Output site names")
            .build();
    private PluginParameter<Boolean> tagCount
            = new PluginParameter.Builder<>("tagCount", false, Boolean.class)
            .description("Output number of sequence tags in file")
            .guiName("Output tag count")
            .build();
    private PluginParameter<Boolean> tagSeqs
            = new PluginParameter.Builder<>("tagNames", false, Boolean.class)
            .description("Output sequence of all tags in file")
            .guiName("Output tag sequences")
            .build();
    private PluginParameter<Boolean> hasDepth
            = new PluginParameter.Builder<>("hasDepth", false, Boolean.class)
            .description("Output whether file contains read depth information")
            .guiName("Output if has read depth")
            .build();
    private PluginParameter<Boolean> printAll
            = new PluginParameter.Builder<>("all", false, Boolean.class)
            .description("Output all available information (overrides other options)")
            .guiName("Output all available information")
            .build();
    private PluginParameter<Boolean> rawData
            = new PluginParameter.Builder<>("rawData", false, Boolean.class)
            .description("Include only the raw output (no metadata or descriptions)")
            .guiName("Output only raw data")
            .build();

    public HDF5SummaryPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }


    @Override
    public DataSet processData(DataSet input) {

        h5reader = HDF5Factory.openForReading(inputFile());
        determineFileType();
        setUpDataStructures();

        //Check for each option in succession and output if needed
        //The check for printAll() is here rather than in postProcessParameters to make it easier to remember for when adding new parameters
        try {
            outputWriter = Utils.getBufferedWriter(outputFile());

            if(! rawData()){
                writeFileData();
            }

            //First batch is the one-line summary stats
            if (siteCount() || printAll()) {
                writeSiteCount();
            }
            if (taxaCount() || printAll()) {
                writeTaxaCount();
            }
            if ( tagCount() || printAll()){
                writeTagCount();
            }
            if (hasDepth() || printAll()) {
                writeHasDepth();
            }

            //Next batch is the larger datasets
            if (siteNames() || printAll()) {
                writeSiteNames();
            }
            if (taxaNames() || printAll()) {
                writeTaxaNames();
            }
            if (tagSeqs() || printAll()){
                writeTagSequences();
            }

            outputWriter.close();
        } catch (IOException e) {
            logger.error("Error writing report to " + outputFile() + ":\n" + e.getStackTrace());
        }

        return null;
    }

    /**
     * Determine the file type by checking for various fields
     */
    private void determineFileType(){
        if(HDF5Utils.doesGenotypeModuleExist(h5reader)){
            myFileType =  H5FileType.GENOTYPE;
        }
        else if(HDF5Utils.doTagsByTaxaExist(h5reader)){
            myFileType = H5FileType.TBT;
        }
        else if(HDF5Utils.doTagsExist(h5reader) && !HDF5Utils.doTagsByTaxaExist(h5reader)){ // Extra check for lack of tags-by-taxa is in case this statement ever gets reordered
            myFileType = H5FileType.TOPM;
        } else{
            myFileType = H5FileType.UNKNOWN;
        }
    }

    private void setUpDataStructures(){
        switch (myFileType){
            case GENOTYPE:
                genos = ImportUtils.readGuessFormat(inputFile());
                break;
            case TBT:
                //tbt = TagsByTaxaHDF5Builder.openForReading(inputFile());
                break;
            case TOPM:
                topm = TOPMUtils.readTOPM(inputFile());
                break;
            default:
                //do nothing
            }
    }

    //Methods to write the different selected data
    private void writeFileData() throws IOException {
        File infile = new File(inputFile());
        String path = infile.getCanonicalFile().getAbsolutePath();
        String output ="### Summary data for TASSEL HDF5 file " + infile.getName() + " (" + path + ")###\n";

        switch(myFileType){
            case GENOTYPE:
                output += "File type:\tGenotypeTable";
                if(HDF5Utils.isTASSEL4HDF5Format(h5reader)){
                    output += " (TASSEL 4 formatted)";
                }
                output += "\n";
                break;
            case TOPM:
                output += "File type:\tTagsOnPhysicalMap\n";
                break;
            case TBT:
                output += "File type:\tTagsByTaxa\n";
                break;
            default:
                output += "File type:\tUnknown\n";
        }

        outputWriter.append(output);
    }

    private void writeSiteCount() throws IOException {
        //Get the number of sites only if a genotype file
        String nsites;
        switch(myFileType){
            case GENOTYPE:
                nsites = "" + genos.numberOfSites();
                break;
            default:
                nsites = "n/a";
        }

        String output;
        if(rawData()) {
            output = nsites + "\n";
        }else{
            output = "numberOfSites:\t" + nsites + "\n";
        }
        outputWriter.append(output);
    }

    private void writeTaxaCount()  throws IOException {
        String ntaxa;
        switch(myFileType) {
            case GENOTYPE: // may need to invoke getHDF5GenotypeTaxaCount
            case TBT:
                ntaxa = "" + HDF5Utils.getHDF5TaxaNumTaxa(h5reader);
                break;
            case TOPM:
            default:
                ntaxa = "n/a";
        }

        String output;
        if(rawData()) {
            output =  ntaxa + "\n";
        }else{
            output = "numberOfTaxa:\t" + ntaxa + "\n";
        }
        outputWriter.append(output);
    }

    private void writeTagCount()  throws IOException {
        String ntags;
        switch(myFileType) {
            case TBT:
            case TOPM:
                ntags = "" + HDF5Utils.getHDF5TagCount(h5reader);
                break;
            default:
                ntags = "n/a";
        }

        String output;
        if(rawData()) {
            output =  ntags + "\n";
        }else{
            output = "numberOfTags:\t" + ntags + "\n";
        }
        outputWriter.append(output);
    }

    private void writeHasDepth()  throws IOException {
        String isDepth;
        switch(myFileType) {
            case GENOTYPE:
                isDepth = "" + HDF5Utils.doesGenotypeDepthExist(h5reader);
                break;
            case TOPM:
                isDepth = "n/a";
                break;
            case TBT:
                isDepth = "true";
                break;
            default:
                isDepth = "unknown";
        }
        String output;
        if(rawData()) {
            output = isDepth + "\n";
        }else{
            output = "hasDepth:\t" +  isDepth + "\n";
        }
        outputWriter.append(output);
    }

    private void writeSiteNames()  throws IOException {
        StringBuilder sites = new StringBuilder();
        if(!rawData()) {
            sites.append("###Site Names###\n");
        }

        switch(myFileType){
            case GENOTYPE:
                for(int i=0; i<genos.numberOfSites(); i++){
                    sites.append(genos.siteName(i) + "\n");
                }
                break;
            default:
                sites.append("(not applicable for this file type)\n");
        }

        outputWriter.append(sites.toString());
    }

    private void writeTaxaNames()  throws IOException {
        StringBuilder taxa = new StringBuilder();
        if(!rawData()) {
            taxa.append("###Taxa Names###\n");
        }

        switch(myFileType){
            case GENOTYPE:  // Original code commented out to test if same funciton works on TBT and Genotype
                /*for(int i=0; i<genos.numberOfTaxa(); i++){
                    taxa.append(genos.taxaName(i) + "\n");
                }
                break;*/
            case TBT:
                List<String> myTaxaList = HDF5Utils.getAllTaxaNames(h5reader);
                for(String t: myTaxaList){
                    taxa.append(t + "\n");
                }
                break;
            default:
                taxa.append("(not applicable for this file type)\n");
        }

        outputWriter.append(taxa.toString());
    }

    private void writeTagSequences()  throws IOException {
        StringBuilder tags = new StringBuilder();
        if(!rawData()) {
            tags.append("###Tag Sequences###\n");
        }

        switch(myFileType){
            case TBT:
            case TOPM:
                long[][] myTags = HDF5Utils.getTags(h5reader);
                for (int i = 0; i < myTags.length; i++) {
                    tags.append(BaseEncoder.getSequenceFromLong(myTags[i]) + "\n");
                }
                break;
            default:
                tags.append("(not applicable for this file type)\n");
        }

        outputWriter.append(tags.toString());
    }

    //Overridden methods from AbstractPlugin
    /*@Override
    protected void postProcessParameters() {
        if (!(taxaCount() || taxaNames() || siteCount() || siteNames() || hasDepth())) {
            throw new IllegalArgumentException("\n\nMust select at least one option to output.\n\n");
        }
    }*/

    @Override
    public String pluginDescription(){
        return "This plugin takes a TASSEL-generated HDF5 file and prints out a set of summary data depending on which " +
                "command-line flags are passed to it. It is meant to allow quick retrieval of certain basic data (taxa names, " +
                "site count, etc.) that are not easily available without operations that could take considerable time " +
                "(e.g., a genotype summary report). This plugin currently supports Genotype, TOPM (tags on physical map), and" +
                "TBT (tags by taxa) files. (Note that the implementation of TOPM and TBT formats in HDF5 is still ongoing, " +
                "so this plugin may not work on them.)\n";
    }

    //Parameter set functions
    public HDF5SummaryPlugin inputFile(String filename) {
        inputFile = new PluginParameter<>(inputFile, filename);
        return this;
    }

    public HDF5SummaryPlugin outputFile(String filename) {
        outputFile = new PluginParameter<>(outputFile, filename);
        return this;
    }

    public HDF5SummaryPlugin printAll(Boolean value) {
        printAll = new PluginParameter<>(printAll, value);
        return this;
    }

    public HDF5SummaryPlugin taxaCount(Boolean value) {
        taxaCount = new PluginParameter<>(taxaCount, value);
        return this;
    }

    public HDF5SummaryPlugin taxaNames(Boolean value) {
        taxaNames = new PluginParameter<>(taxaNames, value);
        return this;
    }

    public HDF5SummaryPlugin siteCount(Boolean value) {
        siteCount = new PluginParameter<>(siteCount, value);
        return this;
    }

    public HDF5SummaryPlugin siteNames(Boolean value) {
        siteNames = new PluginParameter<>(siteNames, value);
        return this;
    }

    public HDF5SummaryPlugin hasDepth(Boolean value) {
        hasDepth = new PluginParameter<>(hasDepth, value);
        return this;
    }

    public HDF5SummaryPlugin rawData(Boolean value) {
        rawData = new PluginParameter<>(rawData, value);
        return this;
    }

    //Parameter get functions
    public String inputFile() {
        return inputFile.value();
    }

    public String outputFile() {
        return outputFile.value();
    }

    public Boolean printAll() {
        return printAll.value();
    }

    public Boolean taxaCount() {
        return taxaCount.value();
    }

    public Boolean taxaNames() {
        return taxaNames.value();
    }

    public Boolean siteCount() {
        return siteCount.value();
    }

    public Boolean siteNames() {
        return siteNames.value();
    }

    public Boolean tagCount() {
        return tagCount.value();
    }

    public Boolean tagSeqs() {
        return tagSeqs.value();
    }

    public Boolean hasDepth() {
        return hasDepth.value();
    }

    public Boolean rawData() {
        return rawData.value();
    }

    //GUI-required methods
    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "HDF5 Summary";
    }

    @Override
    public String getToolTipText() {
        return "HDF5 Summary";
    }

}
