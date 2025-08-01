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
 * Output sequence in tagMap (e.g. {@link TagGWASMap}) file in Fasta format, which is used in Bowtie2 alignment
 * 
 * @author Fei Lu
 */
public class PanATagMapToFastaPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(PanATagMapToFastaPlugin.class);
    
    String tagMapFileS = null;
    String fastaFileS = null;

    public PanATagMapToFastaPlugin() {
        super(null, false);
    }

    public PanATagMapToFastaPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -i  tagMap(e.g. TagGWASMap) file\n"     
                + " -o  output Fasta format sequence file of TagGWASMap\n");
    }

    public DataSet performFunction(DataSet input) {
        TagGWASMap tgm = new TagGWASMap(this.tagMapFileS);
        tgm.writeToFasta(fastaFileS);
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
            engine.add("-o", "--fasta-file", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-i")) {
            this.tagMapFileS = engine.getString("-i");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-o")) {
            this.fastaFileS = engine.getString("-o");
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
