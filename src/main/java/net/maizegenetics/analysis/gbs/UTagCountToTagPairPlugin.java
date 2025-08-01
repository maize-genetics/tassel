/*
 * UTagCountToTagPairPlugin
 */
package net.maizegenetics.analysis.gbs;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.File;

/**
 *
 * @author Fei Lu
 */
public class UTagCountToTagPairPlugin extends AbstractPlugin {

    private static final Logger logger = LogManager.getLogger(UTagCountToTagPairPlugin.class);

    private PluginParameter<Double> errorTolerance
            = new PluginParameter.Builder<>("errorTolerance", 0.03, Double.class)
            .description("What fraction of errors to tolerate when filtering by UNEAK")
            .required(false)
            .guiName("Error tolerance")
            .build();
    private PluginParameter<String> infile
            = new PluginParameter.Builder<>("inputFile", null, String.class)
            .required(true)
            .inFile()
            .guiName("Input file")
            .description("Input file of merged tag counts")
            .build();
    private PluginParameter<String> outfile
            = new PluginParameter.Builder<>("outputFile", null, String.class)
            .required(true)
            .outFile()
            .guiName("Output file")
            .description("Output file of matched tag pairs")
            .build();


    public UTagCountToTagPairPlugin() {
        super(null, false);
    }

    public UTagCountToTagPairPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }


    @Override
    public DataSet processData(DataSet input) {
        String mergedTagCountOfAllS, tagPairS;
        mergedTagCountOfAllS = new File(inputFile()).getAbsolutePath();
        tagPairS = new File(outputFile()).getAbsolutePath();
        UNetworkFilter unf = new UNetworkFilter(mergedTagCountOfAllS, errorTolerance(), tagPairS);
        return null;
    }

    @Override
    public String pluginDescription(){
        return "This plugin takes a set of merged tag counts from the TASSEL GBS pipeline and converts it into a set " +
                "of tag pairs according to the UNEAK filter (see citation). This plugin is intended for GBS with organisms " +
                "that lack a reference genome.";
    }

    @Override
    public String getCitation(){
        return "Lu F, Lipka AE, Elshire RJ, Glaubitz JC, Cherney JH, Casler MD, Buckler ES, Costich DE. (2013)"
         + " Switchgrass genomic diversity, ploidy and evolution: novel insights from a network-based SNP discovery protocol. PLoS Genetics 9(1):e1003215.";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Tag Counts to Tag Pairs";
    }

    @Override
    public String getToolTipText() {
        return "Deprecated: Reference Pipeline is Better";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(UTagCountToTagPairPlugin.class);
    // }

    /**
     * What fraction of errors to tolerate when filtering
     * by UNEAK
     *
     * @return Error tolerance
     */
    public Double errorTolerance() {
        return errorTolerance.value();
    }

    /**
     * Set Error tolerance. What fraction of errors to tolerate
     * when filtering by UNEAK
     *
     * @param value Error tolerance
     *
     * @return this plugin
     */
    public UTagCountToTagPairPlugin errorTolerance(Double value) {
        errorTolerance = new PluginParameter<>(errorTolerance, value);
        return this;
    }

    /**
     * Input file of merged tag counts
     *
     * @return Input file
     */
    public String inputFile() {
        return infile.value();
    }

    /**
     * Set Input file. Input file of merged tag counts
     *
     * @param value Input file
     *
     * @return this plugin
     */
    public UTagCountToTagPairPlugin inputFile(String value) {
        infile = new PluginParameter<>(infile, value);
        return this;
    }

    /**
     * Output file of matched tag pairs
     *
     * @return Output file
     */
    public String outputFile() {
        return outfile.value();
    }

    /**
     * Set Output file. Output file of matched tag pairs
     *
     * @param value Output file
     *
     * @return this plugin
     */
    public UTagCountToTagPairPlugin outputFile(String value) {
        outfile = new PluginParameter<>(outfile, value);
        return this;
    }

}
