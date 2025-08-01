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
import net.maizegenetics.analysis.gbs.SimpleGenotypeSBit;

/** 
 * Reformat HDF5 genotype to {@link SimpleGenotypeSBit} anchor.
 * 
 * @author Fei Lu
 */
public class PanAH5ToAnchorPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(PanAH5ToAnchorPlugin.class);
    
    String h5GenotypeFileS = null;
    String sBitGenotypeFileS = null;

    public PanAH5ToAnchorPlugin() {
        super(null, false);
    }

    public PanAH5ToAnchorPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -i  HDF5 format genotype file\n"
                + " -o  site bit genotype file in SimpleGenotypeSBit format\n");
    }

    public DataSet performFunction(DataSet input) {
        new SimpleGenotypeSBit(h5GenotypeFileS, sBitGenotypeFileS);
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
            engine.add("-o", "--output-dir", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-i")) {
            h5GenotypeFileS = engine.getString("-i");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        
        if (engine.getBoolean("-o")) {
            sBitGenotypeFileS = engine.getString("-o");
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
