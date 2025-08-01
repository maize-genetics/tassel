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
import net.maizegenetics.analysis.gbs.AnnotateTOPM;
import net.maizegenetics.dna.map.TagGWASMap;
import net.maizegenetics.dna.map.TagsOnPhysicalMapV3;

/** 
 * Split large TagsByTaxaByteHDF5TagGroup file into small sub TBTs. Designed to submit genetic mapping jobs in cluster
 * 
 * @author Fei Lu
 */
public class PanASamToMultiPositionTOPMPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(PanASamToMultiPositionTOPMPlugin.class);
    
    String samFileS = null;
    String tagGWASMapFileS = null;
    String topmV3FileS = null;
    int maxNumAlignment = 2;

    public PanASamToMultiPositionTOPMPlugin() {
        super(null, false);
    }

    public PanASamToMultiPositionTOPMPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -i  bowtie2 alignemnt file in SAM format\n"
                + " -t  TagGWASMap file\n"        
                + " -o  output multiple position TOPM file\n");
    }

    public DataSet performFunction(DataSet input) {
        TagGWASMap tgm = new TagGWASMap(this.tagGWASMapFileS); 
        TagsOnPhysicalMapV3.createFile(tgm, topmV3FileS);
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmV3FileS);
        AnnotateTOPM anno = new AnnotateTOPM (topm);
        anno.annotateWithBowtie2(this.samFileS, this.maxNumAlignment);
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
            engine.add("-i", "--sam-file", true);
            engine.add("-t", "--tagGWASMap-file", true);
            engine.add("-o", "--topmV3-file", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-i")) {
            samFileS = engine.getString("-i");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (engine.getBoolean("-t")) {
            tagGWASMapFileS = engine.getString("-t");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-o")) {
            topmV3FileS = engine.getString("-o");
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
