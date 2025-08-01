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
import java.util.concurrent.TimeUnit;
import net.maizegenetics.dna.map.TagGWASMap;

/** 
 * Make predictions based on trained machine learning model. Write predicted values into tagMap file
 * 
 * @author Fei Lu
 */
public class PanAPredictionPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(PanAPredictionPlugin.class);
    
    String tagMap = null;
    String wekaPath = null;
    String modelFileS = null;
    String boxcoxParemeterFileS = null;
    
    public PanAPredictionPlugin() {
        super(null, false);
    }

    public PanAPredictionPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -t  input tagMap (e.g. tagGWASMap) file\n"
                + " -m  trained machine learning model\n"
                + " -b  boxcox parameter file\n"
                + " -w  path of weka library\n");
    }

    public DataSet performFunction(DataSet input) {
        double[] lamdas = null;
        try {
            BufferedReader br = new BufferedReader (new FileReader(this.boxcoxParemeterFileS), 65536);
            br.readLine();
            String[] temp = br.readLine().split("\t");
            lamdas = new double[temp.length];
            for (int i = 0; i < temp.length; i++) lamdas[i] = Double.valueOf(temp[i]);
        }
        catch(Exception e) {
            e.printStackTrace();
        }
        TagGWASMap tgm = new TagGWASMap(this.tagMap);
        File arffFile = new File (new File(this.tagMap).getParent(), "block.arff");
        int blockNum = tgm.getBlockNum();
        for (int i = 0; i < blockNum; i++) {
            int startIndex = i*tgm.getBlockSize();
            int endIndex = startIndex+tgm.getBlockSize();
            if (endIndex > tgm.getTagCount()) endIndex = tgm.getTagCount();
            System.out.println("Start predicting block(Index) " + String.valueOf(i));
            this.generateARFFFileS(tgm, lamdas, startIndex, endIndex, arffFile);
            double[] predictedValue = this.mkPrediction(arffFile.getAbsolutePath(), startIndex, endIndex);
            for (int j = 0; j < predictedValue.length; j++) {
                tgm.getTagGWASMapInfo(j+startIndex).setPredictedDistance(predictedValue[j]);
            }
            tgm.writeBlock(i);
            System.out.println("Predicted distance of block(Index) "+ String.valueOf(i)+" is written");
            System.out.println("");
        }
        arffFile.delete();
        System.out.println("Prediction completed in " + this.tagMap);
        return null;
    }
    
    private double[] mkPrediction (String arffFileS, int startIndex, int endIndex) {
        String cmd = "java -Xms500m -Xmx5g -cp " + this.wekaPath.replace("\\", "/") + " weka.classifiers.rules.M5Rules -p 0 -T " + arffFileS.replace("\\", "/") + " -l " + this.modelFileS.replace("\\", "/");
        System.out.println(cmd);
        Runtime rt = Runtime.getRuntime();
        Process p;
        double[] predictedValue = new double[endIndex-startIndex];
        try {
            p = rt.exec(cmd);
            //p.waitFor();
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()), 65536);
            
            String temp;
            String[] tem;
            for (int i = 0; i < 5; i++) br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                if(temp.isEmpty()) continue;
                tem = temp.trim().split("\\s+");
                predictedValue[cnt] = Double.valueOf(tem[2]);
                cnt++;
            }
            if (cnt!=predictedValue.length) {
                System.out.println("Need to run weka prediction from command line");
                System.out.println(cmd);
                System.exit(1);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return predictedValue;
    }
    
    private void generateARFFFileS (TagGWASMap tgm, double[] lamdas, int startIndex, int endIndex, File arffFile) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(arffFile), 65536);
            bw.write("@relation predictionTag\n\n");
            String[] temp = "TagCount\tTagTaxaCount\tGBinomP\tLRatioSB\tLRatioMB\tGNumSigChr\tGNumSigSite\tGNumSigSiteBC\tGSigWidthBC\tGDist".split("\t");
            for (int i = 0; i < temp.length; i++) {
                bw.write("@attribute " + temp[i] + " numeric\n");
            }
            bw.write("\n@data\n");
            int cnt = 0;
            for (int i = 0; i < endIndex-startIndex; i++) {
                bw.write(tgm.getTagGWASMapInfo(i+startIndex).getBoxcoxAttributesStr(lamdas, ","));
                bw.newLine();
                if (cnt%100000 == 0) System.out.println(String.valueOf(cnt+1)+" transformed instances are written");
                cnt++;
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("Prediciton input file with tag index from " + String.valueOf(startIndex) + " to " + String.valueOf(endIndex) + " is written");
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
            engine.add("-m", "--ml-model", true);
            engine.add("-b", "--boxcox-file", true);
            engine.add("-w", "--weka-path", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-t")) {
            this.tagMap = engine.getString("-t");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (engine.getBoolean("-m")) {
            modelFileS = engine.getString("-m");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-b")) {
            this.boxcoxParemeterFileS = engine.getString("-b");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-w")) {
            this.wekaPath = engine.getString("-w");
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
