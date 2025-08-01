/*
 * BinaryToTextPlugin
 */
package net.maizegenetics.analysis.gbs;

import com.google.common.collect.Range;

import net.maizegenetics.dna.map.TagsOnPhysicalMap;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.tag.TagsByTaxaByte;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.util.Arrays;

/**
 *
 * @author Terry Casstevens
 */
public class BinaryToTextPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(BinaryToTextPlugin.class);

    public static enum FILE_TYPES {

        TOPM, TagCounts, TBTByte

    };

    private PluginParameter<String> myInputFile = new PluginParameter.Builder<String>("i", null, String.class).guiName("Input File").required(true).inFile().build();
    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<String>("o", null, String.class).guiName("Output File").required(true).outFile().build();
    private PluginParameter<FILE_TYPES> myFileType = new PluginParameter.Builder<FILE_TYPES>("t", FILE_TYPES.TOPM, FILE_TYPES.class).guiName("File Type").range(Range.encloseAll(Arrays.asList(FILE_TYPES.values()))).build();

    public BinaryToTextPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        switch (fileType()) {
            case TOPM:
                TagsOnPhysicalMap topm = new TagsOnPhysicalMap(inputFile(), true);
                topm.writeTextFile(new File(outputFile()));
                break;
            case TagCounts:
                TagCounts tc = new TagCounts(inputFile(), FilePacking.Byte);
                tc.writeTagCountFile(outputFile(), FilePacking.Text, 0);
                break;
            case TBTByte:
                TagsByTaxaByte tbtbyte = new TagsByTaxaByte(inputFile(), FilePacking.Byte);
                tbtbyte.writeDistFile(new File(outputFile()), FilePacking.Text, 0);
                break;
        }

        return null;

    }

    @Override
    public String pluginDescription() {
        return "Converts Binary GBS Files into Text format.";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(BinaryToTextPlugin.class);
    // }
    
    /**
     * Input File
     *
     * @return Input File
     */
    public String inputFile() {
        return myInputFile.value();
    }

    /**
     * Set Input File. Input File
     *
     * @param value Input File
     *
     * @return this plugin
     */
    public BinaryToTextPlugin inputFile(String value) {
        myInputFile = new PluginParameter<>(myInputFile, value);
        return this;
    }

    /**
     * Output File
     *
     * @return Output File
     */
    public String outputFile() {
        return myOutputFile.value();
    }

    /**
     * Set Output File. Output File
     *
     * @param value Output File
     *
     * @return this plugin
     */
    public BinaryToTextPlugin outputFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
        return this;
    }

    /**
     * File Type
     *
     * @return File Type
     */
    public FILE_TYPES fileType() {
        return myFileType.value();
    }

    /**
     * Set File Type. File Type
     *
     * @param value File Type
     *
     * @return this plugin
     */
    public BinaryToTextPlugin fileType(FILE_TYPES value) {
        myFileType = new PluginParameter<>(myFileType, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Binary to Text";
    }

    @Override
    public String getToolTipText() {
        return "Binary to Text";
    }
}
