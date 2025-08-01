/*
 * TagCountToFastqPlugin
 */
package net.maizegenetics.analysis.gbs;

import java.awt.Frame;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.IOException;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/*
 * Converts a TagCounts binary (*.cnt) file (presumably a master tag list) to a fastq file that can be used as input
 * for BWA or bowtie2 (and possibly additional aligners).  The same function can be performed with
 * MergeMultipleTagCountPlugin using the -t option and a single Master Tag List file in the input directory, but
 * having a separate plugin to do this reduces confusion and eliminates the risk of merging the master tag list back on
 * itself.
 *
 * @author glaubitz (jcg233) (modified from MergeMultipleTagCountPlugin)
 */
public class TagCountToFastqPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(TagCountToFastqPlugin.class);

    private PluginParameter<String> myInputFile = new PluginParameter.Builder<String>("i", null, String.class).guiName("Input File").required(true).inFile()
            .description("Input binary tag count (*.cnt) file").build();
    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<String>("o", null, String.class).guiName("Output File").required(true).outFile()
            .description("Output fastq file to use as input for BWA or bowtie2").build();
    private PluginParameter<Integer> myMinCount = new PluginParameter.Builder<Integer>("c", 1, Integer.class).guiName("Min Count")
            .description("Minimum count of reads for a tag to be output").build();

    private DataInputStream inStream;
    private DataOutputStream outStream;
    private int nTags, tagLengthInLong;
    private long[] tag = new long[2]; // [indexOfTagLong], emptyLong=0, fileFinish=Long.MAX
    private int tagCount; // tag count
    private byte tagLength;
    private int tagsRead = 0;

    public TagCountToFastqPlugin() {
        super(null, false);
    }

    public TagCountToFastqPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public TagCountToFastqPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        try {
            inStream = new DataInputStream(new BufferedInputStream(new FileInputStream(inputFile()), 655360));
            nTags = inStream.readInt();
            tagLengthInLong = inStream.readInt();
            myLogger.info("Opened the input file: " + inputFile() + "  nTags=" + nTags);
            //outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileName), 655360));
            outStream = Utils.getDataOutputStream(outputFile(), 655360);
            int tagsWritten = 0;
            while (inStream.available() != 0) {
                readNextTag();
                if (tagCount >= minCount()) {
                    writeFASTQ();
                    tagsWritten++;
                }
                if (tagsRead % 500000 == 1) {
                    System.out.printf("tagsRead=%d tagsWritten=%d %n", tagsRead, tagsWritten);
                    myLogger.info(BaseEncoder.getSequenceFromLong(tag));
                }
            }
            outStream.flush();
            outStream.close();
            myLogger.info("Finished converting binary tag count file to fastq."
                    + "\nTotal number of tags read: " + tagsRead
                    + "\nTotal number of tags written: " + tagsWritten + " (above minCount of " + minCount() + ")"
                    + "\nOuput fastq file: " + outputFile() + "\n\n");
        } catch (Exception e) {
            myLogger.info("Catch in reading TagCount file e=" + e);
            e.printStackTrace();
        }
        return null;
    }

    private void readNextTag() {
        try {
            for (int j = 0; j < tagLengthInLong; j++) {
                tag[j] = inStream.readLong();
            }
            tagLength = inStream.readByte();
            tagCount = inStream.readInt();
            tagsRead++;
        } catch (IOException eof) {
            try {
                myLogger.info("Finished reading input file.");
                inStream.close();
                inStream = null;
            } catch (IOException eof2) {
                myLogger.info("Catch closing" + eof2);
                inStream = null;
            }
        }
    }

    private void writeFASTQ() {
        try {
            outStream.writeBytes("@length=" + tagLength + "count=" + tagCount + "\n");   //Length & count header
            String tagSequence = BaseEncoder.getSequenceFromLong(tag);
            tagSequence = tagSequence.substring(0, tagLength);  //Remove any poly-A padding
            outStream.writeBytes(tagSequence + "\n+\n");    //Sequence and "+" symbol
            for (int i = 0; i < tagLength; i++) {
                outStream.writeBytes("f");
            }           //Bogus quality string
            outStream.writeBytes("\n");
        } catch (IOException e) {
            myLogger.info("Catch in writing TagCount file e=" + e);
            e.printStackTrace();
        }
    }

    public String inputFile() {
        return myInputFile.value();
    }

    public TagCountToFastqPlugin inputFile(String value) {
        myInputFile = new PluginParameter<>(myInputFile, value);
        return this;
    }

    public String outputFile() {
        return myOutputFile.value();
    }

    public TagCountToFastqPlugin outputFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
        return this;
    }

    public Integer minCount() {
        return myMinCount.value();
    }

    public TagCountToFastqPlugin minCount(Integer value) {
        myMinCount = new PluginParameter<>(myMinCount, value);
        return this;
    }

    @Override
    public String getToolTipText() {
        return "Tag Count to Fastq";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Tag Count to Fastq";
    }
}
