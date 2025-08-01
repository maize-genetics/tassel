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
import net.maizegenetics.dna.map.TagGWASMap;

/** 
 * Convert mapping result to HDF5 {@link TagGWASMap} format
 * 
 * @author Fei Lu
 */
public class PanABuildTagGWASMapPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(PanABuildTagGWASMapPlugin.class);
    
    String mappingResultFileS = null;
    String tagCountFileS = null;
    String tagGWASMapFileS = null;

    public PanABuildTagGWASMapPlugin() {
        super(null, false);
    }

    public PanABuildTagGWASMapPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -i  tag GWAS mapping result file\n"
                + " -t  tagCount file\n"        
                + " -o  output file in tagGWASMap format\n");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        new TagGWASMap (this.mappingResultFileS, this.tagCountFileS, this.tagGWASMapFileS);
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
            engine.add("-i", "--mapping-result", true);
            engine.add("-t", "--tagCount-file", true);
            engine.add("-o", "--tagGWASMap-file", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-i")) {
            mappingResultFileS = engine.getString("-i");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (engine.getBoolean("-t")) {
            tagCountFileS = engine.getString("-t");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-o")) {
            tagGWASMapFileS = engine.getString("-o");
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
