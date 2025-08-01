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
import java.util.Arrays;

/** 
 * Split large TagsByTaxaByteHDF5TagGroup file into small sub TBTs. Designed to submit genetic mapping jobs in cluster
 * 
 * @author Fei Lu
 */
public class PanAMergeMappingResultPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(PanAMergeMappingResultPlugin.class);
    
    String subResultDirS = null;
    String mergedResultFileS = null;

    public PanAMergeMappingResultPlugin() {
        super(null, false);
    }

    public PanAMergeMappingResultPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -i  directory of mapping results of sub TBTs\n"
                + " -o  filename of merged mapping result\n");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        File[] files = new File (subResultDirS).listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.toLowerCase().endsWith("txt");
            }
        });
        Arrays.sort(files);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(this.mergedResultFileS), 65536);
            BufferedReader br = new BufferedReader (new FileReader(files[0].getAbsolutePath()), 65536);
            String header = br.readLine();
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < files.length; i++) {
                br = new BufferedReader (new FileReader(files[i].getAbsolutePath()), 65536);
                br.readLine();
                String temp;
                while ((temp = br.readLine()) != null) {
                    bw.write(temp);
                    bw.newLine();
                }
                System.out.println(files[i].getAbsolutePath() + " is merged");      
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
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
            engine.add("-o", "--output-dir", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-i")) {
            subResultDirS = engine.getString("-i");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-o")) {
            mergedResultFileS = engine.getString("-o");
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
