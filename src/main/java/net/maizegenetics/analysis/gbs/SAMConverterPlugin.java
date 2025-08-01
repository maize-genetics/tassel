/*
 * SAMConverterPlugin
 */
package net.maizegenetics.analysis.gbs;

import java.awt.Frame;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.map.TagsOnPhysicalMap;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * This class can read in a CBSU TagMapFile into the gbs.TagsOnPhysicalMap data
 * structure.
 *
 * @author harriman
 *
 */
public final class SAMConverterPlugin extends AbstractPlugin {

    boolean cleanCutSites = true;
    private static final Logger myLogger = LogManager.getLogger(SAMConverterPlugin.class);

    private PluginParameter<String> myInputFile = new PluginParameter.Builder<String>("i", null, String.class).guiName("Input File").required(true).inFile()
            .description("Name of input file in SAM text format").build();
    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<String>("o", null, String.class).guiName("Output File").required(true).outFile()
            .description("Name of output file (Default: output.topm.bin)").build();
    private PluginParameter<Integer> myTagLengthInNumLongs = new PluginParameter.Builder<Integer>("l", 2, Integer.class).guiName("Tag Length")
            .description("tag length in integer multiples of 32 bases").build();
    private PluginParameter<Boolean> myTextOutputFormat = new PluginParameter.Builder<Boolean>("t", false, Boolean.class).guiName("Text Output Format")
            .description("Specifies text output format").build();

    public SAMConverterPlugin() {
        super(null, false);
    }

    public SAMConverterPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void postProcessParameters() {
        if ((outputFile() == null) || (outputFile().length() == 0)) {
            if (inputFile() != null) {
                outputFile(Utils.getDirectory(inputFile()) + File.separator + "output.topm.bin");
            }
        }
        if (textOutputFormat()) {
            if (outputFile() != null) {
                outputFile(outputFile().replace(".bin", ".txt"));
            }
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        TagsOnPhysicalMap topm = new TagsOnPhysicalMap();
        topm.readSAMFile(inputFile(), tagLength());
        topm.sort();
        try {
            if (textOutputFormat()) {
                topm.writeTextFile(new File(outputFile()));
            } else {
                topm.writeBinaryFile(new File(outputFile()));
            }
        } catch (Exception e) {
            System.out.println("Catch in writing binary topm file: " + e);
        }
        writeLogFile(topm);
        return null;
    }

    private void writeLogFile(TagsOnPhysicalMap topm) {
        try {
            DataOutputStream report = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile() + ".log"), 65536));
            int[] aligned = topm.mappedTags();
            int unique = 0, multi = 1;  // the indices of aligned
            int unaligned = topm.getTagCount() - aligned[unique] - aligned[multi];
            report.writeBytes(
                    "Input file: " + inputFile() + "\n"
                    + "Output file: " + outputFile() + "\n"
                    + "Total " + topm.getTagCount() + " tags\n\t"
                    + aligned[unique] + " were aligned to unique postions\n\t"
                    + aligned[multi] + " were aligned to multiple postions\n\t"
                    + unaligned + " could not be aligned.\n\n");
            int[] dist = topm.mappingDistribution();
            report.writeBytes("nPositions  nTags\n");
            for (int i = 0; i < dist.length; i++) {
                if (dist[i] > 0) {
                    if (i < 10) {
                        report.writeBytes(i + "           " + dist[i] + "\n");
                    } else if (i < 100) {
                        report.writeBytes(i + "          " + dist[i] + "\n");
                    } else if (i < 1000) {
                        report.writeBytes(i + "         " + dist[i] + "\n");
                    }
                }
            }
            report.close();
        } catch (Exception e) {
            myLogger.warn("Caught exception while writing log file: " + e);
        }
    }

    public String inputFile() {
        return myInputFile.value();
    }

    public SAMConverterPlugin inputFile(String value) {
        myInputFile = new PluginParameter<>(myInputFile, value);
        return this;
    }

    public String outputFile() {
        return myOutputFile.value();
    }

    public SAMConverterPlugin outputFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
        return this;
    }

    public boolean textOutputFormat() {
        return myTextOutputFormat.value();
    }

    public SAMConverterPlugin textOutputFormat(boolean value) {
        myTextOutputFormat = new PluginParameter<>(myTextOutputFormat, value);
        return this;
    }

    public int tagLength() {
        return myTagLengthInNumLongs.value();
    }

    public SAMConverterPlugin tagLength(int value) {
        myTagLengthInNumLongs = new PluginParameter<>(myTagLengthInNumLongs, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "SAM to TOPM Converter";
    }

    @Override
    public String getToolTipText() {
        return "SAM to TOPM Converter";
    }
}
