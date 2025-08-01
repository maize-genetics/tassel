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
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.map.TagGWASMap;
import net.maizegenetics.dna.map.TagGWASMapInfo;

/** 
 * Filtering high quality mapped tag and write them into a file, based on the predicted distance
 * 
 * @author Fei Lu
 */
public class PanAFilteringTagMapPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(PanAFilteringTagMapPlugin.class);
    
    String tagMapFileS = null;
    String anchorFileS = null;
    int distanceCutoff;

    public PanAFilteringTagMapPlugin() {
        super(null, false);
    }

    public PanAFilteringTagMapPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -t  tagMap file\n"
                + " -a  output file\n"        
                + " -c  distance cutoff\n");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        TagGWASMap tgm = new TagGWASMap (this.tagMapFileS);
        double logCut = Math.log10(this.distanceCutoff);
        System.out.println("Distance cutoff is " + String.valueOf(this.distanceCutoff) + " bp");
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(this.anchorFileS), 65536);
            bw.write("Tag\tTagLength\tGChr\tGPos\tIfPAV\tPredictedDistance");
            bw.newLine();
            long[] t;
            TagGWASMapInfo info;
            int cnt = 0;
            for (int i = 0; i < tgm.getTagCount(); i++) {
                info = tgm.getTagGWASMapInfo(i);
                if (info.predictedDistance > logCut) continue;
                t = tgm.getTag(i);
                bw.write(BaseEncoder.getSequenceFromLong(t)+"\t"+String.valueOf(tgm.getTagLength(i))+"\t");
                bw.write(String.valueOf(info.gChr)+"\t"+String.valueOf(info.gPos)+"\t"+String.valueOf(Byte.MIN_VALUE)+"\t"+String.valueOf(info.predictedDistance));
                bw.newLine();
                cnt++;
                if (cnt%100000 == 0) System.out.println(String.valueOf(cnt+1)+" anchors are written");
            }
            bw.flush();
            bw.close();
            System.out.println(String.valueOf(cnt) + " (" + String.valueOf((double)cnt/tgm.getTagCount())+ ") anchors are written to " + this.anchorFileS);
        }
        catch (Exception e) {
            e.printStackTrace();
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
            engine.add("-t", "--tagMap-file", true);
            engine.add("-a", "--anchor-file", true);
            engine.add("-c", "--dictance-cutoff", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-t")) {
            this.tagMapFileS = engine.getString("-t");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (engine.getBoolean("-a")) {
            anchorFileS = engine.getString("-a");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-c")) {
            this.distanceCutoff = Integer.valueOf(engine.getString("-c"));
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
