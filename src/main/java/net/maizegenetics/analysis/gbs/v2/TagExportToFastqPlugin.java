/*
 * TagCountToFastqPlugin
 */
package net.maizegenetics.analysis.gbs.v2;

import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagData;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.concurrent.atomic.LongAdder;

/**
 * Converts a TagCounts binary (*.cnt) file (presumably a master tag list) to a fastq file that can be used as input
 * for BWA or bowtie2 (and possibly additional aligners).  The same function can be performed with
 * MergeMultipleTagCountPlugin using the -t option and a single Master Tag List file in the input directory, but
 * having a separate plugin to do this reduces confusion and eliminates the risk of merging the master tag list back on
 * itself.
 *
 * @author Jeff Glaubitz
 * @author Ed Buckler
 */
public class TagExportToFastqPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(TagExportToFastqPlugin.class);

    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with tags and taxa distribution").build();
    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<String>("o", null, String.class).guiName("Output File").required(true).outFile()
            .description("Output fastq file to use as input for BWA or bowtie2").build();
    private PluginParameter<Integer> myMinCount = new PluginParameter.Builder<Integer>("c", 1, Integer.class).guiName("Min Count")
            .description("Minimum count of reads for a tag to be output").build();


    public TagExportToFastqPlugin() {
        super(null, false);
    }

    public TagExportToFastqPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public TagExportToFastqPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        try {
            BufferedWriter bw=Utils.getBufferedWriter(outputFile());
            TagData tagData=new TagDataSQLite(inputDB());
            LongAdder count=new LongAdder();
            tagData.getTagsWithDepth(minCount())
                    .forEach((tag,depth)->{
                        writeFASTQ(bw, tag); // the tag sequence becomes the header
                        count.increment();
                    });
            bw.close();
            ((TagDataSQLite)tagData).close();  //todo autocloseable should do this but it is not working.

            myLogger.info("Finished converting binary tag count file to fastq."
                    + "\nTotal number of tags written: " + count.longValue() + " (above minCount of " + minCount() + ")"
                    + "\nOuput fastq file: " + outputFile() + "\n\n");
        } catch (Exception e) {
            myLogger.info("Catch in reading TagCount file e=" + e);
            e.printStackTrace();
        }
        return null;
    }

    private void writeFASTQ(BufferedWriter outStream, Tag tag) {
        try {  //build a string first so that if needed this could be parallized in writing
            // With the change for TAS-722 we decided the sequence becomes the header. After aligning,
            // this allows the preservation of the original tag.
            StringBuilder sb=new StringBuilder("@tagSeq=" + tag.sequence() + "\n");
            sb.append(tag.sequence() + "\n+\n");    //Sequence and "+" symbol
            for (int i = 0; i < tag.seqLength(); i++) {
                sb.append("f");
            }           //Bogus quality string
            sb.append("\n");
            outStream.write(sb.toString());
        } catch (IOException e) {
            myLogger.info("Catch in writing TagCount file e=" + e);
            e.printStackTrace();
        }
    }


    @Override
    public String getToolTipText() {
        return "Export unique tags to fastQ file in a format readable by aligners";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Tag Export to Fastq";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(TagExportToFastqPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public String runPlugin(DataSet input) {
        return (String) performFunction(input).getData(0).getData();
    }

    /**
     * Input database file with tags and taxa distribution
     *
     * @return Input DB
     */
    public String inputDB() {
        return myDBFile.value();
    }

    /**
     * Set Input DB. Input database file with tags and taxa
     * distribution
     *
     * @param value Input DB
     *
     * @return this plugin
     */
    public TagExportToFastqPlugin inputDB(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }

    /**
     * Output fastq file to use as input for BWA or bowtie2
     *
     * @return Output File
     */
    public String outputFile() {
        return myOutputFile.value();
    }

    /**
     * Set Output File. Output fastq file to use as input
     * for BWA or bowtie2
     *
     * @param value Output File
     *
     * @return this plugin
     */
    public TagExportToFastqPlugin outputFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
        return this;
    }

    /**
     * Minimum count of reads for a tag to be output
     *
     * @return Min Count
     */
    public Integer minCount() {
        return myMinCount.value();
    }

    /**
     * Set Min Count. Minimum count of reads for a tag to
     * be output
     *
     * @param value Min Count
     *
     * @return this plugin
     */
    public TagExportToFastqPlugin minCount(Integer value) {
        myMinCount = new PluginParameter<>(myMinCount, value);
        return this;
    }
}
