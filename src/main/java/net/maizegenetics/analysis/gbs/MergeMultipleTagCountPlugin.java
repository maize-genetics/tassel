/*
 * MergeMultipleTagCountPlugin
 */
package net.maizegenetics.analysis.gbs;

import java.awt.Frame;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.tag.AbstractTags;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.DirectoryCrawler;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/*
 * Implements an external mergesort to combine multiple tag-count files.
 * @author edbuckler
 */
public class MergeMultipleTagCountPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(MergeMultipleTagCountPlugin.class);

    PluginParameter<String> myInputDir = new PluginParameter.Builder<String>("i", null, String.class).guiName("Input Directory").required(true).inDir()
            .description("Input directory containing .cnt files.").build();
    PluginParameter<String> myOutputFile = new PluginParameter.Builder<String>("o", null, String.class).guiName("Output File").required(true).outFile()
            .description("Output file name.").build();
    PluginParameter<Integer> myMinCount = new PluginParameter.Builder<Integer>("c", 1, Integer.class).guiName("Min Count")
            .description("Minimum count of reads to be output.").build();
    PluginParameter<Boolean> myIsTextOutput = new PluginParameter.Builder<Boolean>("t", false, Boolean.class).guiName("Text Output")
            .description("Specifies that reads should be output in FASTQ text format.").build();

    private long[][] myCtags;  //array of current tags, [indexOfFile][indexOfTagLong], emptyLong=0, fileFinish=Long.MAX
    private int[] myCtagCnt; //array of tag counts [indexOfFile]
    private byte[] myCtagLength;  //array of tag length [indexOfFile]
    private int[] myChunkTagSizes;  // array of number of tags in each file  [indexOfFile]
    private int myTagLengthInLong = 2;
    private int myNumTagsRead = 0;
    private int myOutCnt = 0;
    private int myNumInputStreamsOpen = 0;
    private DataInputStream[] myInputStreams;
    private DataOutputStream myOutStream;

    public MergeMultipleTagCountPlugin() {
        super(null, false);
    }

    public MergeMultipleTagCountPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        String[] inputFileNames = DirectoryCrawler.listFileNames(".*\\.cnt", inputDirectory());
        if (inputFileNames == null || inputFileNames.length == 0) {
            throw new IllegalArgumentException("Couldn't find any files ending in \".cnt\" in the directory you specified: " + inputDirectory());
        }
        myLogger.info("Merging the following .cnt files...");
        for (String filename : inputFileNames) {
            myLogger.info(filename);
        }
        myLogger.info("...to \"" + outputFile() + "\".");

        mergeChunks(inputFileNames, outputFile(), minCount());
        return null;
    }

    public void mergeChunks(String[] chunkFileNames, String outputFileName, int minCount) {
        myInputStreams = new DataInputStream[chunkFileNames.length];
        myChunkTagSizes = new int[myInputStreams.length];
        myCtagCnt = new int[myInputStreams.length];
        myCtagLength = new byte[myInputStreams.length];
        myNumInputStreamsOpen = myInputStreams.length;
        try {
            for (int f = 0; f < myInputStreams.length; f++) {
                String infile = chunkFileNames[f];
                myInputStreams[f] = new DataInputStream(new BufferedInputStream(new FileInputStream(infile), 4000000));
                myChunkTagSizes[f] = myInputStreams[f].readInt();
                myTagLengthInLong = myInputStreams[f].readInt();
                myLogger.info("Opened :" + infile + " tags=" + myChunkTagSizes[f]);
            }
            myOutStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileName + ".fq"), 655360));
            myCtags = new long[myInputStreams.length][myTagLengthInLong];
            outStack("B:");
            int t = 0;
            while (myNumInputStreamsOpen > 0) {
                long[] minTag = updateCurrentTags();
                //     myLogger.info(BaseEncoder.getSequenceFromLong(minTag));
                //     outStack("U:"+t);
                if (textOutput()) {
                    writeFASTQ(minTag, minCount);
                } else {
                    writeTags(minTag, minCount);
                }
                //     outStack("A:"+t);
                if (t % 1000000 == 0) {
                    System.out.printf("t=%d tagsRead=%d outCnt=%d rwOpen=%d %n", t, myNumTagsRead, myOutCnt, myNumInputStreamsOpen);
                    outStack("A:" + t);
                    myLogger.info(BaseEncoder.getSequenceFromLong(minTag));
                }
                t++;
            }
            myOutStream.flush();
            myOutStream.close();
        } catch (Exception e) {
            myLogger.info("Catch in reading TagCount file e=" + e);
            e.printStackTrace();
        }
        //Binary files need a header
        if (!textOutput()) {
            prependHeader(outputFileName, myOutCnt, myTagLengthInLong);
        }
    }

    private void outStack(String prefix) {
        myLogger.info(prefix + ":");
        for (int f = 0; f < myInputStreams.length; f++) {
            myLogger.info(myCtags[f][0] + ":");
        }
        myLogger.info("");
    }

    private void writeTags(long[] writeTag, int minCount) {
        int count = 0;
        int tagLength = -1;

        //Loop through input buffers, compare current records, write smallest to output buffer
        for (int f = 0; f < myInputStreams.length; f++) {
            if (AbstractTags.compareTags(myCtags[f], writeTag) == 0) {
                count += myCtagCnt[f];
                tagLength = myCtagLength[f];
                for (int j = 0; j < myTagLengthInLong; j++) {
                    myCtags[f][j] = 0;
                }
            }
        }
        if (count >= minCount) {
            myOutCnt++;
        } else {
            return;
        }
        //write to file here.
        try {
            long test1;
            for (int i = 0; i < myTagLengthInLong; i++) {
                test1 = writeTag[i];
                myOutStream.writeLong(test1);
            }
            byte test2 = (byte) tagLength;
            int test3 = count;
            myOutStream.writeByte(test2);
            myOutStream.writeInt(test3);
            if (count == 0) {
                myLogger.info("");
            }
        } catch (IOException e) {
            myLogger.info("Catch in writing TagCount file e=" + e);
            e.printStackTrace();
        }
    }

    private void writeFASTQ(long[] writeTag, int minCount) {
        int count = 0;
        int tagLength = -1;
        String tagSequence = "";

        //Loop through input buffers, compare current records, write smallest to output buffer
        for (int f = 0; f < myInputStreams.length; f++) {
            if (AbstractTags.compareTags(myCtags[f], writeTag) == 0) {
                count += myCtagCnt[f];
                tagLength = myCtagLength[f];
                for (int j = 0; j < myTagLengthInLong; j++) {
                    myCtags[f][j] = 0;
                }
            }
        }
        if (count >= minCount) {
            myOutCnt++;
        } else {
            return;
        }
        //write to file here.
        try {
            myOutStream.writeBytes("@length=" + tagLength + "count=" + count + "\n");   //Length & count header
            tagSequence = BaseEncoder.getSequenceFromLong(writeTag);
            tagSequence = tagSequence.substring(0, tagLength);  //Remove any poly-A padding
            myOutStream.writeBytes(tagSequence + "\n+\n");    //Sequence and "+" symbol
            for (int i = 0; i < tagLength; i++) {
                myOutStream.writeBytes("f");
            }           //Bogus quality string
            myOutStream.writeBytes("\n");
        } catch (IOException e) {
            myLogger.info("Catch in writing TagCount file e=" + e);
            e.printStackTrace();
        }
    }

    private static void prependHeader(String fileName, int tagCount, int tagLengthInLong) {
        myLogger.info("Adding header to " + fileName + ".");
        File inputFile = new File(fileName + ".fq");
        File outputFile = new File(fileName);
        int tagsWritten = 0;
        try {
            DataInputStream input = new DataInputStream(new BufferedInputStream(new FileInputStream(inputFile), 4000000));
            DataOutputStream output = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile), 65536));

            output.writeInt(tagCount);
            output.writeInt(tagLengthInLong);

            for (int i = 0; i < tagCount; i++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    output.writeLong(input.readLong());
                }
                output.writeByte(input.readByte());
                output.writeInt(input.readInt());
                tagsWritten++;
                if (tagsWritten % 1000000 == 0) {
                    myLogger.info("Wrote " + tagsWritten + " records.");
                }
            }
            input.close();
            output.close();

            if (!inputFile.delete()) {
                myLogger.info("WARNING: Failure to delete file:\n\t" + inputFile.getCanonicalPath());  //Delete old file
            }
        } catch (Exception e) {
            myLogger.info("Caught exception while prepending read count to read count file: " + e);
            e.printStackTrace();
        }

    }

    private long[] updateCurrentTags() {
        long[] minTag = new long[myTagLengthInLong];
        minTag[0] = Long.MAX_VALUE;
        for (int f = 0; f < myInputStreams.length; f++) {
            if (myCtags[f][0] == 0) {
                readNextTag(f);
            }
            if (AbstractTags.compareTags(myCtags[f], minTag) < 0) {
                minTag = myCtags[f].clone();
            }
        }
        return minTag;
    }

    private void readNextTag(int f) {
        if (myInputStreams[f] == null) {
            return;
        }
        try {
            for (int j = 0; j < myTagLengthInLong; j++) {
                myCtags[f][j] = myInputStreams[f].readLong();
            }
            myCtagLength[f] = myInputStreams[f].readByte();
            myCtagCnt[f] = myInputStreams[f].readInt();
            myNumTagsRead++;
        } catch (IOException eof) {
            try {
                myLogger.info("Finished reading file " + f + ".");
                myInputStreams[f].close();
                myInputStreams[f] = null;
                for (int i = 0; i < myTagLengthInLong; i++) {

                    myCtags[f][i] = Long.MAX_VALUE;
                }
                myNumInputStreamsOpen--;
            } catch (IOException eof2) {
                myLogger.info("Catch closing" + eof2);
                myInputStreams[f] = null;
            }
        }
    }

    public String inputDirectory() {
        return myInputDir.value();
    }

    public MergeMultipleTagCountPlugin inputDirectory(String value) {
        myInputDir = new PluginParameter<>(myInputDir, value);
        return this;
    }

    public String outputFile() {
        return myOutputFile.value();
    }

    public MergeMultipleTagCountPlugin outputFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
        return this;
    }

    public Integer minCount() {
        return myMinCount.value();
    }

    public MergeMultipleTagCountPlugin minCount(Integer value) {
        myMinCount = new PluginParameter<>(myMinCount, value);
        return this;
    }

    public Boolean textOutput() {
        return myIsTextOutput.value();
    }

    public MergeMultipleTagCountPlugin textOutput(Boolean value) {
        myIsTextOutput = new PluginParameter<>(myIsTextOutput, value);
        return this;
    }

    @Override
    public String getToolTipText() {
        return "Merge Multiple Tag Count Files";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Merge Multiple Tag Count Files";
    }
}
