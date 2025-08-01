/*
 * Plugin for the UNEAK pipeline to convert the Tag Pair file into a (fake) TOPM file for reintegration with the normal GBS pipeline
 */
package net.maizegenetics.analysis.gbs;

import com.google.common.collect.Range;
import net.maizegenetics.dna.map.TagsOnPhysicalMap;
import net.maizegenetics.dna.tag.AbstractTags;
import net.maizegenetics.dna.tag.UTagPairs;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.ArgsEngine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.File;

/**
 *
 * @author Jason Wallace
 *
 * This plugin takes the TagPairs file produced by earlier steps in the UNEAK
 * pipeline and converts it to a TOPM alignment file that can be used by
 * DiscoverySNPCallerPlugin to call SNPs as if from a reference genome.
 *
 */
public class UTagPairToTOPMPlugin extends AbstractPlugin {

    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(UTagPairToTOPMPlugin.class);

    //Internal data for converting tag pairs
    private UTagPairs tp;   //Structure to read in tag pairs
    byte myStrand = 1, myMultimaps = 1, myDcoP = Byte.MIN_VALUE, myMapP = Byte.MIN_VALUE;   //Dummy values to feed to the TOPM structure
    //int myMaxMapping = 1, myMaxVariants = 8;    //HDF5-specific dummy values for output; not currently in use since DiscoverySNPCallerPlugin doesn't take HDF5-TOPM as input

    //Parameters for this plugin
    private PluginParameter<Integer> chrom
            = new PluginParameter.Builder<>("chrom", 1, Integer.class)
            .description("Chromosome to start numbering at")
            .guiName("Start chromosome")
            .build();
    private PluginParameter<Integer> distance
            = new PluginParameter.Builder<>("distance", 1000, Integer.class)
            .description("Distance to pad between each tag pair")
            .guiName("Pad distance")
            .units("base pairs")
            .range(Range.closed(1, Integer.MAX_VALUE))
            .build();
    private PluginParameter<String> infile
            = new PluginParameter.Builder<>("input", null, String.class)
            .description("Input file of matched tag pairs")
            .guiName("Input tag pair file")
            .required(true)
            .inFile()
            .build();
    private PluginParameter<String> textOutputFile
            = new PluginParameter.Builder<>("toText", null, String.class)
            .description("File to output TOPM in text format")
            .guiName("Output file (text)")
            .required(false)
            .outFile()
            .build();
    private PluginParameter<String> binaryOutputFile
            = new PluginParameter.Builder<>("toBinary", null, String.class)
            .description("File to output TOPM in binary format")
            .guiName("Output file (binary)")
            .required(false)
            .outFile()
            .build();


    public UTagPairToTOPMPlugin() {
        super(null, false);
    }

    public UTagPairToTOPMPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        //Make TOPM
        tp = new UTagPairs(input());
        TagsOnPhysicalMap myTopm = makeTopmFromTagPairs(tp);

        //Output to specified file(s)
        if (toText() != null) {
            logger.info("Outputting TOPM in text format to " + toText());
            myTopm.writeTextFile(new File(toText()));
        }
        if (toBinary() != null) {
            logger.info("Outputting TOPM in binary format to " + toBinary());
            myTopm.writeBinaryFile(new File(toBinary()));
        }

        return null;
    }

    /**
     * Convert the TagPair object to a TOPM object in memory, filling in things
     * with dummy data where needed
     */
    private TagsOnPhysicalMap makeTopmFromTagPairs(UTagPairs tp) {
        HelperTags tempTags = new HelperTags(tp.getTagNum(), tp.getTag(0).length);
        //Load in all the data for each tag pair; having trouble finding the function to load actual sequences
        for (int i = 0; i < tp.getTagNum(); i++) {
            tempTags.setTag(i, tp.getTag(i), tp.getTagLength(i));
        }
        TagsOnPhysicalMap myTopm = new TagsOnPhysicalMap(tempTags);
        long currPos = 1;
        int currChrom = startChromosome();
        for (int i = 0; i < tp.getTagNum(); i++) {

            myTopm.setChromoPosition(i, currChrom, myStrand, (int) currPos, (int) currPos + tp.getTagLength(i) - 1);
            myTopm.setDivergence(i, (byte) 0);

            //These may not be necessary; don't know
            myTopm.setDcoP(i, Byte.MIN_VALUE);  //May have to alter TOPM class to do this
            myTopm.setMapP(i, Byte.MIN_VALUE);
            myTopm.setMultimaps(i, (byte) 1);

            //Increment position after odd-numbered tags (so pairs are at the same position)
            if (i % 2 == 1) {
                currPos += padDistance();
            }
            //If over max interger value, increment chromosome and start over
            if (currPos >= Integer.MAX_VALUE ) {
                currChrom++;
                currPos = 1;
            }
        }
        return myTopm;
    }

    @Override
    protected void postProcessParameters(){
        //Test that at least one output file supplied
        if( (toBinary() == null || "".equals(toBinary())) && (toText() == null || "".equals(toText())) ){
            throw new IllegalArgumentException("\n\nMust specify at least one output file (text or binary).\n\n");
        }
        if(padDistance() < 100){
            logger.warn("Warning! Setting the pad distance to <100 base pairs risks overlapping tags pairs with each other (and thus calling false SNPs)");
        }
    }

    @Override
    public String pluginDescription(){
        return "This plugin takes a tag-pair file (from UTagCountToTagPairPlugin) and converts it into a TOPM (tags on physical map) " +
                "file for use in the GBS pipeline. The resulting chromosome coordinates and other data are just filler to comply with " +
                "the TOPM file specifications.";
    }

    @Override
    public ImageIcon getIcon() {
        return(null);
    }

    @Override
    public String getButtonName() {
        return("Tag Pairs to TOPM");
    }

    @Override
    public String getToolTipText() {
        return "Deprecated: Reference Pipeline is Better";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(UTagPairToTOPMPlugin.class);
    // }

    /**
     * Chromosome to start numbering at
     *
     * @return Start chromosome
     */
    public Integer startChromosome() {
        return chrom.value();
    }

    /**
     * Set Start chromosome. Chromosome to start numbering
     * at
     *
     * @param value Start chromosome
     *
     * @return this plugin
     */
    public UTagPairToTOPMPlugin startChromosome(Integer value) {
        chrom = new PluginParameter<>(chrom, value);
        return this;
    }

    /**
     * Distance to pad between each tag pair
     *
     * @return Pad distance
     */
    public Integer padDistance() {
        return distance.value();
    }

    /**
     * Set Pad distance. Distance to pad between each tag
     * pair
     *
     * @param value Pad distance
     *
     * @return this plugin
     */
    public UTagPairToTOPMPlugin padDistance(Integer value) {
        distance = new PluginParameter<>(distance, value);
        return this;
    }

    /**
     * Input file of matched tag pairs
     *
     * @return Input
     */
    public String input() {
        return infile.value();
    }

    /**
     * Set Input. Input file of matched tag pairs
     *
     * @param value Input
     *
     * @return this plugin
     */
    public UTagPairToTOPMPlugin input(String value) {
        infile = new PluginParameter<>(infile, value);
        return this;
    }

    /**
     * File to output TOPM in text format
     *
     * @return To Text
     */
    public String toText() {
        return textOutputFile.value();
    }

    /**
     * Set To Text. File to output TOPM in text format
     *
     * @param value To Text
     *
     * @return this plugin
     */
    public UTagPairToTOPMPlugin toText(String value) {
        textOutputFile = new PluginParameter<>(textOutputFile, value);
        return this;
    }

    /**
     * File to output TOPM in binary format
     *
     * @return To Binary
     */
    public String toBinary() {
        return binaryOutputFile.value();
    }

    /**
     * Set To Binary. File to output TOPM in binary format
     *
     * @param value To Binary
     *
     * @return this plugin
     */
    public UTagPairToTOPMPlugin toBinary(String value) {
        binaryOutputFile = new PluginParameter<>(binaryOutputFile, value);
        return this;
    }
}

/*
 * A small helper class that exists solely to pass tag info to TOPM
 */
class HelperTags extends AbstractTags {

    /* Inherited
     protected int tagLengthInLong;  //
     protected long[][] tags;  //Index in rows, actual tag components in columns
     protected byte[] tagLength;  // length of tag (number of bases)  // 1 byte
     */
    public HelperTags(int numtags, int myTagLengthInLong) {
        tagLengthInLong = myTagLengthInLong;
        tags = new long[tagLengthInLong][numtags];
        tagLength = new byte[numtags];
    }

    public void setTag(int index, long[] tagValue, byte myTagLength) {
        tagLength[index] = myTagLength;
        for (int i = 0; i < tagValue.length; i++) {
            tags[i][index] = tagValue[i];
        }
    }
}
