package net.maizegenetics.analysis.gbs.pana;

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
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TagGroups;

/** 
 * Split large TagsByTaxaByteHDF5TagGroup file into small sub TBTs. Designed to submit genetic mapping jobs in cluster
 * 
 * @author Fei Lu
 */
public class PanASplitTBTPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(PanASplitTBTPlugin.class);
    
    String tbtHDF5 = null;
    String outDirS = null;
    int chunkSize = 65536;

    public PanASplitTBTPlugin() {
        super(null, false);
    }

    public PanASplitTBTPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -i  input TagsByTaxa(TBT) file, TagsByTaxaByteHDF5TagGroups format\n"
                + " -s  chunkSize, number of tags in a sub TBT. This determines the mapping calculation time usage in a node/computer. Default = 65536\n"        
                + " -o  output directory of sub TBTs\n");
    }

    public DataSet performFunction(DataSet input) {
        TagsByTaxaByteHDF5TagGroups tbt = new TagsByTaxaByteHDF5TagGroups (tbtHDF5);
        int tagNum = tbt.getTagCount();
        int fileNum;
        int base = tagNum%chunkSize;
        if (base == 0) fileNum = tagNum/chunkSize;
        else fileNum = tagNum/chunkSize+1;
        System.out.println("TBT will be split into " + String.valueOf(fileNum) + " subTBTs");
        int big = 10;
        while (fileNum>big) big*=10;
        for (int i = 0; i < fileNum; i++) {
            int actualSize = chunkSize;
            if (i == fileNum-1 && base !=0) actualSize = base;
            int actualStart = i*chunkSize;
            int actualEnd = actualStart+actualSize;
            String fileName = "0"+String.valueOf(big+i).substring(1);
            fileName = "pivotTBT_"+fileName+"_"+String.valueOf(actualStart)+"_"+String.valueOf(actualSize)+".h5";
            fileName = new File(outDirS, fileName).getAbsolutePath();
            int[] selectIndex = new int[actualSize];
            for (int j = 0; j < actualSize; j++) {
                selectIndex[j] = actualStart+j;
            }
            tbt.writeDistFile(fileName, selectIndex);
            System.out.println("Completed "+String.valueOf((double)(i+1)/fileNum));
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
            engine.add("-i", "--input-TBT", true);
            engine.add("-s", "--chunk-size", true);
            engine.add("-o", "--output-dir", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-i")) {
            tbtHDF5 = engine.getString("-i");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (engine.getBoolean("-s")) {
            chunkSize = Integer.valueOf(engine.getString("-s"));
        } 
        else {
            System.out.println("Default chunkSize is "+String.valueOf(chunkSize));
            printUsage();
        }
        
        if (engine.getBoolean("-o")) {
            outDirS = engine.getString("-o");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
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
