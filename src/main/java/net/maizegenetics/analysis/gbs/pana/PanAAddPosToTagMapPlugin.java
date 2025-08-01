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
import net.maizegenetics.dna.map.TagsOnPhysicalMapV3;

/** 
 * Add alignment information from Bowtie2 to tagMap. Used to identify unique ref tags for model training
 * 
 * @author Fei Lu
 */
public class PanAAddPosToTagMapPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(PanAAddPosToTagMapPlugin.class);
    
    String tagMapFileS = null;
    String topmV3FileS = null;

    public PanAAddPosToTagMapPlugin() {
        super(null, false);
    }

    public PanAAddPosToTagMapPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -i  tagMap (e.g. tagGWASMap) file\n"
                + " -t  multiple position TOPM file\n");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        TagGWASMap tgm = new TagGWASMap(this.tagMapFileS);
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(this.topmV3FileS);
        tgm.addAlignment(topm);
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
            engine.add("-i", "--tagMap-file", true);
            engine.add("-t", "--topmV3-file", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-i")) {
            tagMapFileS = engine.getString("-i");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-t")) {
            topmV3FileS = engine.getString("-t");
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
