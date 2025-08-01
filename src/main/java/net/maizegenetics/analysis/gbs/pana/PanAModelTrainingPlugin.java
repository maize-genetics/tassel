package net.maizegenetics.analysis.gbs.pana;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

/**
 * Training data with M5Rules model, generate model file and training report files
 *
 * @author Fei Lu
 */
public class PanAModelTrainingPlugin extends AbstractPlugin {
    
    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(PanAModelTrainingPlugin.class);
    
    String trainingSetFileS = null;
    String wekaPath = null;
    String modelFileS = null;
    String reportDirS = null;
    
    String predictionFileS = "prediction.txt";
    String accuracyTableFileS = "accuracyTable.txt";
    public PanAModelTrainingPlugin() {
        super(null, false);
    }
    
    public PanAModelTrainingPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }
    
    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                        + " -t  training data set file, including unique reference tag\n"
                        + " -w  path of weka library\n"
                        + " -m  M5 model file\n"
                        + " -r  directory of training report\n");
    }
    
    @Override
    public DataSet performFunction(DataSet input) {
        this.trainModel();
        this.generateReport();
        return null;
    }
    
    private void generateReport () {
        int[] predictionCut = {10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000};
        int[] accuracyCut = {10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000};
        double[] logPreCut = new double[predictionCut.length];
        for (int i = 0; i < logPreCut.length; i++) logPreCut[i] = Math.log10(predictionCut[i]);
        double[] logAccCut = new double[accuracyCut.length];
        for (int i = 0; i < logAccCut.length; i++) logAccCut[i] = Math.log10(accuracyCut[i]);
        double[] actualValue = null;
        double[] predictionValue = null;
        try {
            BufferedReader br = new BufferedReader(new FileReader(this.predictionFileS), 65536);
            br.readLine();
            String temp;
            String[] tem;
            int cnt = 0;
            while ((temp = br.readLine()) != null) cnt++;
            actualValue = new double[cnt];
            predictionValue = new double[cnt];
            br = new BufferedReader(new FileReader(this.predictionFileS), 65536);
            br.readLine();
            for (int i = 0; i < cnt; i++) {
                tem = br.readLine().split("\t");
                actualValue[i] = Double.valueOf(tem[1]);
                predictionValue[i] = Double.valueOf(tem[2]);
            }
            br.close();
            BufferedWriter bw = new BufferedWriter (new FileWriter(this.accuracyTableFileS), 65536);
            bw.write("PredictionCutoff\tProportionRemain");
            for (int i = 0; i < accuracyCut.length; i++) {
                bw.write("\t<"+String.valueOf(accuracyCut[i]/1000)+ " kb");
            }
            bw.newLine();
            for (int i = 0; i < logPreCut.length; i++) {
                double[] ratio = new double[accuracyCut.length];
                cnt = 0;
                for (int j = 0; j < actualValue.length; j++) {
                    if (predictionValue[j] > logPreCut[i]) continue;
                    for (int k = 0; k < logAccCut.length; k++) {
                        if (actualValue[j] < logAccCut[k]) ratio[k]++;
                    }
                    cnt++;
                }
                bw.write(String.valueOf(predictionCut[i]/1000)+"Kb\t"+String.valueOf((double)cnt/actualValue.length));
                for (int j = 0; j < logAccCut.length; j++) {
                    bw.write("\t"+String.valueOf((double)ratio[j]/cnt));
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("Accuracy table is at " + this.accuracyTableFileS);
    }
    
    private void trainModel () {
        //only works on windows, > is not supported in Runtime exec
        //String cmd = "java -Xms500m -Xmx5g -cp " + this.wekaPath.replace("\\", "/") + " weka.classifiers.rules.M5Rules -p 0 -t " + this.trainingSetFileS.replace("\\", "/") + " -d " + this.modelFileS.replace("\\", "/") + " > " + this.predictionFileS.replace("\\", "/");
        String cmd = "java -Xms500m -Xmx5g -cp " + this.wekaPath.replace("\\", "/") + " weka.classifiers.rules.M5Rules -p 0 -t " + this.trainingSetFileS.replace("\\", "/") + " -d " + this.modelFileS.replace("\\", "/");
        System.out.println(cmd);
        System.out.println("\nStart model training");
        System.out.println("Be patient. It takes a while...\n");
        Runtime rt = Runtime.getRuntime();
        rt.gc();
        Process p;
        ArrayList<Double> obList = new ArrayList();
        ArrayList<Double> preList = new ArrayList();
        try {
            p = rt.exec(cmd);
            //p.waitFor();
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()), 655360);
            BufferedWriter bw = new BufferedWriter(new FileWriter(this.predictionFileS), 65536);
            bw.write("Instance\tAcutal\tPredicted\tError");
            bw.newLine();
            String temp;
            String[] tem;
            for (int i = 0; i < 5; i++) { 
                System.out.println(br.readLine());
            }
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
                if(temp.isEmpty()) continue;
                tem = temp.trim().split("\\s+");
                bw.write(String.valueOf(cnt)+"\t"+tem[1]+"\t"+tem[2]+"\t"+tem[3]);
                obList.add(Double.valueOf(tem[1]));
                preList.add(Double.valueOf(tem[2]));
                bw.newLine();
                cnt++;
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        double[] ob = new double[obList.size()];
        double[] pre = new double[preList.size()];
        for (int i = 0; i < ob.length; i++) {
            ob[i] = obList.get(i);
            pre[i] = preList.get(i);
        }
        double r = new PearsonsCorrelation().correlation(ob, pre);
        System.out.println("Pearson correlation coefficient (observation Vs prediction):\t" + String.valueOf(r));
        System.out.println("Training model is generated in " + this.modelFileS);
        System.out.println("Ten fold cross validation prediction at " + this.predictionFileS);
    }
    
    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine == null) {
            engine = new ArgsEngine();
            engine.add("-t", "--training-data", true);
            engine.add("-w", "--weka-path", true);
            engine.add("-m", "--output-model", true);
            engine.add("-r", "--report-dir", true);
            engine.parse(args);
        }
        
        if (engine.getBoolean("-t")) {
            this.trainingSetFileS = engine.getString("-t");
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
        
        
        if (engine.getBoolean("-m")) {
            modelFileS = engine.getString("-m");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-r")) {
            reportDirS = engine.getString("-r");
            File reportDir = new File(reportDirS);
            reportDir.mkdirs();
            this.predictionFileS = new File(reportDir.getAbsoluteFile(), predictionFileS).getAbsolutePath();
            this.accuracyTableFileS = new File(reportDir.getAbsoluteFile(), accuracyTableFileS).getAbsolutePath();
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
