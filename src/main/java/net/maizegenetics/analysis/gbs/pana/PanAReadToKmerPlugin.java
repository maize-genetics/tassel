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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeSet;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.tag.TagCountMutable;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.tag.UTagCountMutable;

/** 
 * Generate Kmers from Fastq/Qseq data
 * Do not use underscore in taxa name in keyfile
 * @author Fei Lu
 */
public class PanAReadToKmerPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = LogManager.getLogger(PanAReadToKmerPlugin.class);
    
    String polyA = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    
    //user definable parameters
    String rawSeqDirS = null;
    String keyFileS = null;
    int customTagLength = 0;
    String outDirS = null;
    int inputFormatValue;
    
    //calculatable parameters
    int tagLengthInLong;
    String[] lanes = null;
    HashMap<String, String> laneTaxaMap = null;

    public PanAReadToKmerPlugin() {
        super(null, false);
    }

    public PanAReadToKmerPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -i  input directory of Fastq or Qseq files\n"
                + " -f  input format value.  0 = Fastq. 1 = Qseq\n"
                + " -k  key file which links Fastq/Qseq file with samples\n"
                + " -l  customed tag length\n" 
                + " -o  output directory of tag count files\n");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        this.getLaneTaxaInfo();
        if (this.inputFormatValue == 0) {
            this.processFastq();
        }
        else if (this.inputFormatValue == 1) {
            this.processQseq();
        }
        else {
            
        }
        
        this.mergeTagCounts();
        return null;
    }
    
    private void mergeTagCounts () {
        File[] tcFiles = new File (outDirS).listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.toLowerCase().endsWith("cnt");
            }
        });
        Arrays.sort(tcFiles);
        TreeSet<String> taxaNameSet = new TreeSet();
        for (int i = 0; i < tcFiles.length; i++) {
            String[] temp = tcFiles[i].getName().substring(0, tcFiles[i].getName().length()-4).split("_");
            String theTaxon = temp[temp.length-1];
            taxaNameSet.add(theTaxon);
        }
        String[] taxaNames = taxaNameSet.toArray(new String[taxaNameSet.size()]);
        for (int i = 0; i < taxaNames.length; i++) {
            System.out.println("Merging tagCount files for taxa " + taxaNames[i]);
            ArrayList<File> fileList = new ArrayList();
            for (int j = 0; j < tcFiles.length; j++) {
                String[] temp = tcFiles[j].getName().substring(0, tcFiles[j].getName().length()-4).split("_");
                String theTaxon = temp[temp.length-1];
                if (taxaNames[i].equals(theTaxon)) fileList.add(tcFiles[j]);
            }
            if (fileList.size() == 1) {
                File dest = new File(fileList.get(0).getParent(), taxaNames[i]+".cnt");
                fileList.get(0).renameTo(dest);
                System.out.println("File is renamed to " + dest.getAbsolutePath());
            }
            else if (fileList.size() > 1){
                File dest = new File(fileList.get(0).getParent(), taxaNames[i]+".cnt");
                UTagCountMutable tcu = new UTagCountMutable(this.tagLengthInLong, 0);
                for (int j = 0; j < fileList.size(); j++) {
                    TagCounts tc = new TagCounts(fileList.get(j).getAbsolutePath(), FilePacking.Byte);
                    for (int k = 0; k < tc.getTagCount(); k++) {
                        long[] t = tc.getTag(k);
                        tcu.addReadCount(t, tc.getTagLength(k), tc.getReadCount(k));
                    }
                }
                tcu.toArray();
                tcu.collapseCounts();
                tcu.writeTagCountFile(dest.getAbsolutePath(), FilePacking.Byte, 1);
                for (int j = 0; j < fileList.size(); j++) {
                    fileList.get(j).delete();
                }
                System.out.println(String.valueOf(fileList.size()) + " files are merged into " + dest.getAbsolutePath());
            }
            System.out.println();
        }
    }
    
    private void processFastq () {
        BufferedReader br;
        for (int i = 0; i < lanes.length; i++) {
            File infile = new File(rawSeqDirS, lanes[i]);
            if (!infile.exists()) {
                System.out.println(infile.getAbsolutePath() + " doesn't exist. Program stops");
                System.exit(1);
            }
            String infileS = infile.getAbsolutePath();
            try {
                String outfileS = new File (this.outDirS,lanes[i]+"_"+laneTaxaMap.get(lanes[i])+".cnt").getAbsolutePath();
                if (infileS.endsWith(".gz")) {
                    br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(infile))));
                } 
                else {
                    br = new BufferedReader(new FileReader(infile), 65536);
                }
                System.out.println("Start reading " + infile.getAbsolutePath());
                String temp = br.readLine();
                temp = br.readLine();
                System.out.println("Sequence length is " + temp.length() + " bps in file " + infileS);
                int actualLength = temp.length();
                if (temp.length() < this.customTagLength) {
                    System.out.println("Custom tag length is longer than actual sequence length in file. Programs stops");
                    System.exit(0);
                }
                br.readLine();
                br.readLine();
                //int numKmersLine = actualLength-this.customTagLength+1;
                int numKmersLine = 1;
                
                UTagCountMutable tcu = new UTagCountMutable(this.tagLengthInLong, 0);
                String tag = null;
                long[] t;
                long tagCnt = 0;
                while ((temp = br.readLine()) != null) {
                    temp = br.readLine();
                    for (int j = 0; j < numKmersLine; j++) {
                        tag = temp.substring(j, j+this.customTagLength);
                        if (!this.isBadSequence(tag)) {
                            t = BaseEncoder.getLongArrayFromSeq(tag);
                            tcu.addReadCount(t, this.customTagLength, 1);
                            tagCnt++;
                        }
                    }
                    br.readLine();
                    br.readLine();
                }
                br.close();
                tcu.toArray();
                tcu.collapseCounts();
                tcu.writeTagCountFile(outfileS, TagsByTaxa.FilePacking.Byte, 1);
                System.out.println("Generated " + String.valueOf(tagCnt)+ " Kmers/tags from " + infile.getAbsolutePath());
                System.out.println();
            }
            catch (Exception e)     {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }
    
    private void processQseq () {
        BufferedReader br;
        for (int i = 0; i < lanes.length; i++) {
            File infile = new File(rawSeqDirS, lanes[i]);
            if (!infile.exists()) {
                System.out.println(infile.getAbsolutePath() + " doesn't exist. Program stops");
                System.exit(1);
            }
            String infileS = infile.getAbsolutePath();
            try {
                String outfileS = new File (this.outDirS,lanes[i]+"_"+laneTaxaMap.get(lanes[i])+".cnt").getAbsolutePath();
                if (infileS.endsWith(".gz")) {
                    br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(infile))));
                } 
                else {
                    br = new BufferedReader(new FileReader(infile), 65536);
                }
                System.out.println("Start reading " + infile.getAbsolutePath());
                String temp = br.readLine();
                String[] tem = temp.split("\t");
                System.out.println("Sequence length is " + tem[8].length() + " bps in file " + infileS);
                int actualLength = tem[8].length();
                if (tem[8].length() < this.customTagLength) {
                    System.out.println("Custom tag length is longer than actual sequence length in file. Programs stops");
                    System.exit(0);
                }
                //int numKmersLine = actualLength-this.customTagLength+1;
                int numKmersLine = 1;
                
                UTagCountMutable tcu = new UTagCountMutable(this.tagLengthInLong, 0);
                String tag = null;
                long[] t;
                long tagCnt = 0;
                while ((temp = br.readLine()) != null) {
                    for (int j = 0; j < numKmersLine; j++) {
                        tag = temp.substring(j, j+this.customTagLength);
                        if (!this.isBadSequence(tag)) {
                            t = BaseEncoder.getLongArrayFromSeq(tag);
                            tcu.addReadCount(t, this.customTagLength, 1);
                            tagCnt++;
                        }
                    }
                }
                br.close();
                tcu.toArray();
                tcu.collapseCounts();
                tcu.writeTagCountFile(outfileS, TagsByTaxa.FilePacking.Byte, 1);
                System.out.println("Generated " + String.valueOf(tagCnt)+ " Kmers/tags from " + infile.getAbsolutePath());
                System.out.println();
            }
            catch (Exception e)     {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }
    
    private boolean isBadSequence (String seq) {
        if (seq.indexOf("N") != -1) return true;
        if (seq.indexOf(".") != -1) return true;
        return false;
    }
    
    private void getLaneTaxaInfo () {
        laneTaxaMap = new HashMap();
        ArrayList<String> laneList =new ArrayList();
        try {
            BufferedReader br = new BufferedReader(new FileReader(keyFileS), 65536);
            String temp = br.readLine();
            while ((temp = br.readLine())!=null) {
                String[] tem = temp.split("\t");
                laneList.add(tem[0]);
                laneTaxaMap.put(tem[0], tem[1]);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        lanes = laneList.toArray(new String[laneList.size()]);
    }
    
    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine == null) {
            engine = new ArgsEngine();
            engine.add("-i", "--input-DirS", true);
            engine.add("-f", "--input-format", true);
            engine.add("-k", "--key-file", true);
            engine.add("-l", "--costum-tagLength", true);
            engine.add("-o", "--output-DirS", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-i")) {
            rawSeqDirS = engine.getString("-i");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-f")) {
            this.inputFormatValue = Integer.valueOf(engine.getString("-f"));
            if (this.inputFormatValue != 0 && this.inputFormatValue != 1) {
                printUsage();
                throw new IllegalArgumentException("\n\nPlease double input format value.\n\n");
            }
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-k")) {
            keyFileS = engine.getString("-k");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
            
        if (engine.getBoolean("-l")) {
            customTagLength = Integer.valueOf(engine.getString("-l"));
            int base = customTagLength%BaseEncoder.chunkSize;
            if (base == 0) {
                this.tagLengthInLong = customTagLength/BaseEncoder.chunkSize;
            }
            else {
                this.tagLengthInLong = customTagLength/BaseEncoder.chunkSize+1;
            }
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-o")) {
            outDirS = engine.getString("-o");
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
