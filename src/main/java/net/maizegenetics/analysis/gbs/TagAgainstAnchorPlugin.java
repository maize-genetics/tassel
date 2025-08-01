package net.maizegenetics.analysis.gbs;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.MultiMemberGZIPInputStream;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.*;

/** 
 * Genetic mapping of GBS tags.
 * Steps:
 * 1. Using -pc and -t option to calculate the number of TBT chunks.
 * 2. When -pc is 0, run the program on cluster. Submit jobs to nodes by specifying chunkStartIndex and chunkEndIndex
 * @author Fei Lu
 */
public class TagAgainstAnchorPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(TagAgainstAnchorPlugin.class);
    int preCalculation = 0;
    String hapMapHDF5 = null;
    String tbtHDF5 = null;
    String blockFileS = null;
    String outfileS = null;
    double pThresh = 0.000001;
    int minCount = 20;
    int coreNum = -1;
    int chunkSize = 32768;
    int chunkStartIndex = 0;
    int chunkEndIndex = 1;

    public TagAgainstAnchorPlugin() {
        super(null, false);
    }

    public TagAgainstAnchorPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -pc pre-caculation for chunk number, value = 0/1. Default = 0; when 1, should only be used with -t option together\n"
                + " -g  input anchor map file, SimpleGenotypeSBit format\n"
                + " -t  input TBT file, TagsByTaxaByteHDF5TagGroup format\n"
                + " -b  input TagBlockPosition file, correspongding to tags in TBT. Used to block the marker coming from the tag to be mapped. Default = null\n"
                + " -o  output file of genetic mapping\n"
                + " -m  minimum count when tag appear in taxa, default = 20, too low number lacks statistical power\n"
                + " -c  coreNum, value = max/Integer. Default:max, which means using all cores in a node. When the coreNum is set less than total core number, which means using coreNum cores, each core runs 1 thread\n"
                + " -s  chunkSize, number of tags in a chunk. This determines the time usage in a node/computer. Default = 32768\n"
                + " -cs chunkStartIndex, start index of chunk\n"
                + " -ce chunkEndIndex, end index of chunk\n\n");
    }

    public DataSet performFunction(DataSet input) {
        if (preCalculation == 1) {
            TagAgainstAnchor.getChunkNum(tbtHDF5, chunkSize);
        }
        else {
            TagAgainstAnchor taa = new TagAgainstAnchor(hapMapHDF5, tbtHDF5, blockFileS, outfileS, pThresh, minCount, coreNum, chunkSize, chunkStartIndex, chunkEndIndex);
        }    
        return null;
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine == null) {
            engine = new ArgsEngine();
            engine.add("-pc", "--pre-calculation", true);
            engine.add("-g", "--anchor-file", true);
            engine.add("-t", "--TBT-file", true);
            engine.add("-b", "--TagBlock-file", true);
            engine.add("-o", "--output-file", true);
            engine.add("-m", "--min-tagInTaxa", true);
            engine.add("-c", "--core-num", true);
            engine.add("-s", "--chunk-size", true);
            engine.add("-cs", "--chunkStart-index", true);
            engine.add("-ce", "--chunkEnd-index", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-pc")) {
            preCalculation = Integer.valueOf(engine.getString("-pc"));
        } 

        if (engine.getBoolean("-g")) {
            hapMapHDF5 = engine.getString("-g");
        } 

        if (engine.getBoolean("-t")) {
            tbtHDF5 = engine.getString("-t");
        } 

        if (engine.getBoolean("-b")) {
            blockFileS = engine.getString("-b");
        }

        if (engine.getBoolean("-o")) {
            outfileS = engine.getString("-o");
        }

        if (engine.getBoolean("-m")) {
            minCount = Integer.valueOf(engine.getString("-m"));
        } 
        
        if (engine.getBoolean("-c")) {
            if (engine.getString("-c").substring(0, 1).matches("\\d")){
                int numOfProcessors = Runtime.getRuntime().availableProcessors();
                coreNum = Integer.valueOf(engine.getString("-c"));
                if (coreNum > numOfProcessors) coreNum = numOfProcessors;
                System.out.println("This node has " + numOfProcessors + " processors");
            }
        } 
        
        if (engine.getBoolean("-s")) {
            chunkSize = Integer.valueOf(engine.getString("-s"));
        }
        
        if (engine.getBoolean("-cs")) {
            chunkStartIndex = Integer.valueOf(engine.getString("-cs"));
        } 
        
        if (engine.getBoolean("-ce")) {
            chunkEndIndex = Integer.valueOf(engine.getString("-ce"));
        }
    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
