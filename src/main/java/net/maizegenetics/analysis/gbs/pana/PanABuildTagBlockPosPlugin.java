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
import net.maizegenetics.analysis.gbs.TagBlockPosition;

/** 
 * Build blocked physical position of tags. The blocked physical position of tag is used to block the corresponding marker in genetic mapping if the tag is mapping to the marker coming from itself
 * Positions come from TOPM alignment hypothesis or the best position from machine learning prediction
 * 
 * @author Fei Lu
 */
public class PanABuildTagBlockPosPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(PanABuildTagBlockPosPlugin.class);
    
    String tbtHDF5 = null;
    String topmFileS = null;
    String tagBlockFileS = null;
    int topmVersion = 2;

    public PanABuildTagBlockPosPlugin() {
        super(null, false);
    }

    public PanABuildTagBlockPosPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -t  input TagsByTaxa(TBT) file, TagsByTaxaByteHDF5TagGroup format\n"
                + " -p  input TOPM file\n"        
                + " -v  TOPM version value. Binary file = 1; HDF5 file = 2\n"        
                + " -o  output tagBlockPosition file\n");
    }

    public DataSet performFunction(DataSet input) {
        TagBlockPosition tbp = new TagBlockPosition(tbtHDF5, topmFileS, topmVersion);
        tbp.writeTagBlockPosition(tagBlockFileS);   
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
            engine.add("-t", "--input-TBT", true);
            engine.add("-p", "--TOPM-file", true);
            engine.add("-v", "--TOPM-version", true);
            engine.add("-o", "--output-tagBlock", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-t")) {
            tbtHDF5 = engine.getString("-t");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-p")) {
            topmFileS = engine.getString("-p");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (engine.getBoolean("-v")) {
            topmVersion = Integer.valueOf(engine.getString("-v"));
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-o")) {
            tagBlockFileS = engine.getString("-o");
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
