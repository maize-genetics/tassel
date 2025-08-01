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
import java.util.logging.Level;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.map.TagGWASMap;

/** 
 * Build training data set from tagMap, including boxcox transformation and converting to ARFF format
 * 
 * @author Fei Lu
 */
public class PanABuildTrainingSetPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(PanABuildTrainingSetPlugin.class);
    
    String tagMap = null;
    String trainingSetFileS = null;
    String rPath = null;
    String boxcoxParemeterFileS = null;
    int maxInstance = 30000;

    public PanABuildTrainingSetPlugin() {
        super(null, false);
    }

    public PanABuildTrainingSetPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -m  tagMap file\n"
                + " -t  training data set file\n"
                + " -i  max number of instances\n"        
                + " -r  R path\n"
                + " -b  boxcox parameter file\n");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        TagGWASMap tgm = new TagGWASMap(tagMap);
        this.writeOriginalTrainingSet(tgm);
        this.creatLamdaFile();
        this.transformTrainingSet(tgm);
        return null;
    }
    
    private void transformTrainingSet (TagGWASMap tgm) {
        double[] lamdas;
        File transformFile = new File(new File(this.trainingSetFileS).getParent(), "trans.arff");
        try {
            BufferedReader br = new BufferedReader(new FileReader(this.boxcoxParemeterFileS), 65536);
            br.readLine();
            String[] temp = br.readLine().split("\t");
            lamdas = new double[temp.length];
            for (int i = 0; i < lamdas.length; i++) {
                lamdas[i] = Double.valueOf(temp[i]);
            }
            br.close();
            br = new BufferedReader(new FileReader(this.trainingSetFileS), 65536);
            String header = br.readLine();
            br.close();
            BufferedWriter bw = new BufferedWriter(new FileWriter(transformFile), 65536);
            bw.write("@relation uniqueRefTag\n\n");
            temp = header.split("\t");
            for (int i = 0; i < temp.length; i++) {
                bw.write("@attribute " + temp[i] + " numeric\n");
            }
            bw.write("\n@data\n");
            int cnt = 0;
            for (int i = 0; i < tgm.getTagCount(); i++) {
                if (!tgm.getTagGWASMapInfo(i).isUniqueRef()) continue;
                bw.write(tgm.getTagGWASMapInfo(i).getBoxcoxAttributesStr(lamdas, ","));
                bw.newLine();
                if (cnt%100000 == 0) System.out.println(String.valueOf(cnt+1)+" transformed instances are written");
                cnt++;
                if (cnt == this.maxInstance) break;
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        new File(this.trainingSetFileS).delete();
        transformFile.renameTo(new File(this.trainingSetFileS));
    }
    
    private void creatLamdaFile () {
        File scriptFile = new File (new File (this.boxcoxParemeterFileS).getParent(), "box.r");
        String temp;
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(scriptFile), 65536);
            bw.write("library(MASS)\n");
            temp = this.trainingSetFileS;
            temp = temp.replace("\\", "/");
            bw.write("data <- read.table(\""+temp+"\", header=TRUE, sep=\"\t\")\n");
            bw.write("lamdas = matrix(nrow=1,ncol=ncol(data))\n");
            bw.write("for (i in 1:ncol(data)) {\n");
            bw.write("b=boxcox(data[,i]~1)\n");
            bw.write("ymax <- max(b$y, na.rm=T)\n");
            bw.write("lamda <- b$x[b$y==ymax]\n");
            bw.write("lamdas[,i]=lamda\n");
            bw.write("}\n");
            bw.write("colnames(lamdas)=colnames(data)\n");
            temp = this.boxcoxParemeterFileS;
            temp = temp.replace("\\", "/");
            bw.write("write.table(lamdas, \"" + temp + "\", sep=\"\\t\", col.names=T, row.names=F)\n");
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        String cmd = this.rPath + " " + scriptFile.getAbsolutePath();
        Runtime rt = Runtime.getRuntime();
        Process p;
        try {
            p = rt.exec(cmd);
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        scriptFile.delete();
        System.out.println("Boxcox paremeter files is generated");
    }
    
    private void writeOriginalTrainingSet (TagGWASMap tgm) {
        System.out.println("Start writing training set with original values");
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(this.trainingSetFileS), 65536);
            bw.write("TagCount\tTagTaxaCount\tGBinomP\tLRatioSB\tLRatioMB\tGNumSigChr\tGNumSigSite\tGNumSigSiteBC\tGSigWidthBC\tGDist");
            bw.newLine();
            int cnt = 0;
            for (int i = 0; i < tgm.getTagCount(); i++) {
                if (!tgm.getTagGWASMapInfo(i).isUniqueRef()) continue;
                bw.write(tgm.getTagGWASMapInfo(i).getAttributesStr("\t"));
                bw.newLine();
                if (cnt%100000 == 0) System.out.println(String.valueOf(cnt+1)+" instances are written");
                cnt++;
                if (cnt == this.maxInstance) break;
            }
            System.out.println(String.valueOf(cnt) + " instances in total");
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("Original training set written to " + this.trainingSetFileS);
    }
    
    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine == null) {
            engine = new ArgsEngine();
            engine.add("-m", "--tagMap-file", true);
            engine.add("-t", "--training-file", true);
            engine.add("-i", "--max-instance", true);
            engine.add("-r", "--r-path", true);
            engine.add("-b", "--boxcox-dir", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-m")) {
            tagMap = engine.getString("-m");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (engine.getBoolean("-t")) {
            trainingSetFileS = engine.getString("-t");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-i")) {
            this.maxInstance = Integer.valueOf(engine.getString("-i"));
        } 
        
        if (engine.getBoolean("-r")) {
            rPath = engine.getString("-r");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
         if (engine.getBoolean("-b")) {
            boxcoxParemeterFileS = engine.getString("-b");
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
