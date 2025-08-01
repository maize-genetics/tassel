package net.maizegenetics.util;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.ImageIcon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.google.common.collect.ImmutableMap;

import net.maizegenetics.analysis.gbs.Barcode;
import net.maizegenetics.analysis.gbs.v2.BarcodeTrie;
import net.maizegenetics.analysis.gbs.v2.EnzymeList;
import net.maizegenetics.analysis.gbs.v2.GBSUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;

/**
 * This plugin takes a single fastq or a directory of fastq files
 * in the "old" format and converts them to the "modern" illumina format.  
 * 
 * The old format had multiple taxa combined in a single file, with
 * barcodes attaached to each read sequence.
 * 
 * The new fastQ format has the taxa (samples) each in individual files.
 * 
 * THis method takes the old fastq and for each read block, subsets
 * the sequence to start after the barcode (but including the initial cut site)
 * and end after the ending cut site.  The quality score is subsetted 
 * identically to the sequence so the 2 match.
 * 
 * Files are written to a user defined directory.
 * 
 * @author lcj34
 *
 */
public class ConvertOldFastqToModernFormatPlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(ConvertOldFastqToModernFormatPlugin.class);
    
    private PluginParameter<String> inputFile = new PluginParameter.Builder<>("i", null, String.class).guiName("Input Directory").required(true)
            .description("Input file or directory containing FASTQ files in text or gzipped text.\n"
                    + "     NOTE: Directory will be searched recursively").build();
    private PluginParameter<String> keyFile = new PluginParameter.Builder<>("k", null, String.class).guiName("Key File").required(true).inFile()
            .description("Key file listing barcodes distinguishing the samples").build();
    private PluginParameter<String> projectName = new PluginParameter.Builder<>("p", null, String.class).guiName("Project File").required(true)
            .description("Name for this project.  Project name becomes part of the newly created fastq file name").build();
    private PluginParameter<String> outputDir = new PluginParameter.Builder<>("o", null, String.class).guiName("Output Directory").required(true).outDir()
            .description("Output directory where new fastQ files will be written").build();
    private PluginParameter<String> myEnzyme = new PluginParameter.Builder<>("e", null, String.class).guiName("Enzyme").required(true)
            .description("Enzyme used to create the GBS library").build();
    
    public static final String inputFileGlob="glob:*{.fq,fq.gz,fastq,fastq.txt,fastq.gz,fastq.txt.gz,_sequence.txt,_sequence.txt.gz}";
    public static final String sampleNameField="FullSampleName";
    public static final String flowcellField="Flowcell";
    public static final String laneField="Lane";
    public static final String barcodeField="Barcode";
    public static final String tissueNameField = "Tissue";
    public static final String fileNameField = "FileName";
    
    static List<Path> infiles;
    static String[] likelyReadEndStrings;
    protected static int readEndCutSiteRemnantLength;
    
    public ConvertOldFastqToModernFormatPlugin() {
        super(null, false);
    }

    public ConvertOldFastqToModernFormatPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    @Override
    protected void postProcessParameters() {
        
        // create list of directory files
        File dirList = new File(inputFile());
        if (!(dirList.exists())) {
            throw new IllegalStateException("Input file or directory not found !!") ;
        } 
        if (dirList.isDirectory()) {
            infiles = DirectoryCrawler.listPaths(GBSUtils.inputFileGlob, Paths.get(inputFile.value()).toAbsolutePath());
            if(infiles.isEmpty()) {
                myLogger.warn("No files matching:"+ GBSUtils.inputFileGlob);
                throw new IllegalStateException("no .txt files found in input directory !!");
            }
        } else {
            infiles=new ArrayList<>();
            Path filePath= Paths.get(inputFile()).toAbsolutePath(); 
            infiles.add(filePath);
        }
        if (!myEnzyme.isEmpty()) {
            // Add likelyReadEnds for later processing
            EnzymeList.Enzyme enzyme = EnzymeList.defaultCache.getEnzyme(enzyme()); 
            likelyReadEndStrings = enzyme.likelyReadEnd; // for removeSecondCutSiteIndexOf()
            readEndCutSiteRemnantLength = enzyme.readEndCutSiteRemnantLength;
        }
    }
    
    @Override
    public  DataSet processData(DataSet input) {
        
        // Check if keyFile contains FullSampleName field.
        // If not, create masterTaxaList using
        // sample:Flowcell:lane:barcode concatenated columns as header
        
        myLogger.info("Checking Key format");
        boolean keyFormatIsNew = determineKeyFormat(keyFile());
        
        TaxaList masterTaxaList;
        if (keyFormatIsNew) {
            masterTaxaList= TaxaListIOUtils.readTaxaAnnotationFile(keyFile(), GBSUtils.sampleNameField, 
                    new HashMap<>(),true);
        } else {
            masterTaxaList = readOldTaxaAnnotationFile(keyFile(), new HashMap<>(),true);                   
        }
                 
        myLogger.info("Created masterTaxaList - start processing fastq files");

        int batchSize = 8;
        for (int idx = 0; idx < infiles.size(); idx+=batchSize) {
            int end = idx+batchSize;
            if (end > infiles.size()) end = infiles.size();
            ArrayList<Path> sub = new ArrayList<Path>();
            for (int j = idx; j < end; j++) sub.add(infiles.get(j));
            myLogger.info("\nStart processing batch " + String.valueOf(idx/batchSize+1));
            sub.parallelStream()
            .forEach(inputSeqFile -> {
                try {
                    //
                    Path keyPath= Paths.get(keyFile()).toAbsolutePath();
                    myLogger.info(" - processing file " + inputSeqFile);
                    Map<Taxon,Tuple<BufferedWriter,String>> taxaWriterMap = new HashMap<Taxon,Tuple<BufferedWriter,String>>();
                    processFastQFile(keyPath, inputSeqFile, enzyme(), masterTaxaList,taxaWriterMap,
                             outputDir(), projectName(), keyFormatIsNew);
                    // End of file - close the writers;
                    myLogger.info("Calling Shutdown to close the writers for " + inputSeqFile);
                    Shutdown(taxaWriterMap);

                } catch (StringIndexOutOfBoundsException oobe) {
                    return;
                }
            });
        }       

        return null;
    }
    
    // Determine the key file format.  IF headers do not include "FullSampleName",
    // it is an old format.
    public static boolean determineKeyFormat(String keyFile){
        
        try (BufferedReader br = Utils.getBufferedReader(keyFile)) {
            String header = br.readLine();
            if (header.contains(sampleNameField)) return true;
            else return false;
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
        return false;
    }
    
    public static ArrayList<Taxon> getLaneAnnotatedTaxaList(Path keyPath, Path fastQpath, boolean useNew) {
        String[] filenameField = fastQpath.getFileName().toString().split("_");
        ArrayList<Taxon> annoTL;
        if (filenameField.length == 3) {
            if (useNew) {
                annoTL = TaxaListIOUtils.readTaxaAnnotationFileAL(keyPath.toAbsolutePath().toString(), sampleNameField,
                        ImmutableMap.of(flowcellField, filenameField[0], laneField, filenameField[1])); 
            } else {
                annoTL = readOldFormatKeyFile(keyPath.toAbsolutePath().toString(), 
                        ImmutableMap.of(flowcellField, filenameField[0], laneField, filenameField[1]), true);

            }
        } else if (filenameField.length == 4) {
            if (useNew) {
                annoTL = TaxaListIOUtils.readTaxaAnnotationFileAL(keyPath.toAbsolutePath().toString(),sampleNameField,
                        ImmutableMap.of(flowcellField, filenameField[0], laneField, filenameField[2]));
            } else {
                annoTL = readOldFormatKeyFile(keyPath.toAbsolutePath().toString(),
                        ImmutableMap.of(flowcellField, filenameField[0], laneField, filenameField[2]),true);
            }
            
        }
        else if (filenameField.length == 5) {
            if (useNew) {
                annoTL = TaxaListIOUtils.readTaxaAnnotationFileAL(keyPath.toAbsolutePath().toString(),sampleNameField,
                        ImmutableMap.of(flowcellField, filenameField[1], laneField, filenameField[3]));
            } else {
                annoTL = readOldFormatKeyFile(keyPath.toAbsolutePath().toString(),
                        ImmutableMap.of(flowcellField, filenameField[1], laneField, filenameField[3]),true);
            }
            
        } else {
            myLogger.error("Error in parsing file name: " + fastQpath.toString());
            myLogger.error("   The filename does not contain either 3, 4, or 5 underscore-delimited values.");
            myLogger.error("   Expect: flowcell_lane_fastq.txt.gz OR flowcell_s_lane_fastq.txt.gz OR code_flowcell_s_lane_fastq.txt.gz");
            return null;
        }
        return annoTL;
    }
    
    public static TaxaList readOldTaxaAnnotationFile(String fileName,  Map<String, String> filters, boolean mergeSameNames) {    
        // create list
        ArrayList<Taxon> taxaAL = readOldFormatKeyFile( fileName,  filters, mergeSameNames);
        if (taxaAL == null) return null;
        TaxaListBuilder tlb = new TaxaListBuilder();
        taxaAL.stream().forEach(taxa -> {
                if (mergeSameNames) {
                        tlb.addOrMerge(taxa);
                } else {
                        tlb.add(taxa);
                }
        });
        return tlb.sortTaxaAlphabetically().build();    
    }
    // Call used for old format of files
    public static ArrayList<Taxon> readOldFormatKeyFile(String fileName, Map<String, 
            String> filters, boolean mergeSameNames) {
        try {
            BufferedReader fileIn = Utils.getBufferedReader(fileName, 1000000);
            fileIn.mark(1 << 16);
            String line = fileIn.readLine();
            ArrayList<Taxon> taxaAL = new ArrayList<Taxon>();

            // Name becomes concatenation of 4 fields:  Sample, Flowcell, Lane, Barcode - in that order
            int indexOfSample = 0;
            int indexOfFlowcell = 0;
            int indexOfLane = 0;
            int indexOfBarcode = 0;
            
            //parse headers
            List<String> headers = new ArrayList<>();
            List<Boolean> isQuant = new ArrayList<>();
            if (line.contains("Flowcell")) {
                int idx = 0;
                for (String header : line.split("\\t")) {
                    if (header.equalsIgnoreCase("Sample")) {
                        indexOfSample= idx;
                    }
                    if (header.equalsIgnoreCase("Flowcell")) {
                        indexOfFlowcell= idx;
                    }
                    if (header.equalsIgnoreCase("Lane")) {
                        indexOfLane= idx;
                    }
                    if (header.equalsIgnoreCase("Barcode")) {
                        indexOfBarcode= idx;
                    }
                    isQuant.add(header.startsWith("#") || header.startsWith("<#"));
                    headers.add(header.replace(">", "").replace("<", "").replace("#", ""));
                    idx++;
                }
            } else {
                fileIn.reset();
            }
            //parse taxa rows
            while ((line = fileIn.readLine()) != null) {
                String[] stokens = line.split("\\t");
                StringBuilder taxonSB = new StringBuilder();
                taxonSB.append(stokens[indexOfSample]).append(":");
                taxonSB.append(stokens[indexOfFlowcell]).append(":");
                taxonSB.append(stokens[indexOfLane]).append(":");
                taxonSB.append(stokens[indexOfBarcode]);
                String aTaxonName = taxonSB.toString(); 
                Taxon.Builder anID = new Taxon.Builder(aTaxonName);
                for (int idx = 0; idx < stokens.length; idx++) {

                    String[] cs = stokens[idx].split(";");
                    for (String ta : cs) {
                        if (ta == null || ta.isEmpty()) {
                            continue;
                        }
                        if (isQuant.get(idx)) {
                            if (ta.equals("NA")) {
                                anID.addAnno(headers.get(idx), Double.NaN);
                            } else {
                                anID.addAnno(headers.get(idx), Double.parseDouble(ta));
                            }
                        } else {
                            anID.addAnno(headers.get(idx), ta);
                        }
                    }
                }
                Taxon myTaxon = anID.build();
                if (TaxaListIOUtils.doesTaxonHaveAllAnnotations(myTaxon, filters)) {                  
                        taxaAL.add(myTaxon);
                }
            }
            // Sort alphabetically based on name.  This is to remain consistent
            // with the readTaxaAnnotationFile(), which alphabetizes the taxaList
            Collections.sort(taxaAL, new Comparator<Taxon>() {
                @Override
                public int compare(Taxon taxa1, Taxon taxa2) {
                    return taxa1.getName().compareTo(taxa2.getName());
                }
            });
            return taxaAL;
 
        } catch (Exception e) {
            System.err.println("Error in Reading Annotated Taxon File:" + fileName);
            e.printStackTrace();
        }
        return null;
    }
    
    // THis method will take the taxon/outputWriter map, a fastqfile, and an enzyme name.
    // From this data a barcode trie is created.  Then for each "read" in the fastq file,
    // the taxa is identified from the barcode.  The barcode is stripped, and the read
    // is written to the appropriate file based on the taxaFileMap
    public static void processFastQFile(Path keyPath, Path fastQPath, String enzymeName,
            TaxaList masterTaxaList,Map<Taxon,Tuple<BufferedWriter,String>> taxaFileMap,
            String outputDir, String projectName, boolean keyFormatIsNew){

        int allReads=0, goodBarcodedReads = 0, barCodeNotFound = 0, badEndReads = 0;
        int currentTotal = 0;
        int first20hasN = 0;
    
        ArrayList<Taxon> tl = getLaneAnnotatedTaxaList(keyPath, fastQPath, keyFormatIsNew);;
        
        // Creating a hashmap of taxon name to taxon object.  The file name we want to
        // write to has flowcell and lane in the name.  These annotations are stored with
        // the tl arraylist taxons, but NOT in the barcode taxons.  It is necessary to
        // take the taxon name from the barcode, find the taxon in the tl list, then
        // retrieve the flowcell and lane.
        Map<String,Taxon> nameToTaxon = new HashMap<String,Taxon>();
        for (Taxon taxon : tl) {
            String name = taxon.getName();
            nameToTaxon.put(name, taxon);
        }
        if (tl.size() == 0) return;
        BarcodeTrie barcodeTrie=GBSUtils.initializeBarcodeTrie(tl, masterTaxaList, EnzymeList.defaultCache.getEnzyme(enzymeName));

        String[] seqAndQual;
        BufferedReader br = Utils.getBufferedReader(fastQPath.toString(), 1 << 22);
        long time=System.nanoTime();
        try {
            while ((seqAndQual=readAllFastQBlock(br,allReads)) != null) {
                // Unlike GBS, we read all 4 lines, as all 4 lines will be written to new file
                allReads++;
                currentTotal++;
                //After quality score is read, decode barcode using the current sequence & quality  score
                Barcode barcode=barcodeTrie.longestPrefix(seqAndQual[1]);
                if(barcode==null) {
                    if (seqAndQual[1].contains("N")) first20hasN++;
                    barCodeNotFound++;
                    continue; // couldn't find taxa - skip this read
                }
                int barcodeLen = barcode.getBarLength(); 
                String adjustedSeq = seqAndQual[1].substring(barcodeLen);
                String adjustedQual = seqAndQual[3].substring(barcodeLen);
                
                Tuple<Integer,Integer> stringStartEnd = truncateToSecondCutSite(adjustedSeq,enzymeName); 
                if (stringStartEnd == null) {
                    badEndReads++;
                    continue;
                }
                adjustedSeq = adjustedSeq.substring(stringStartEnd.x, stringStartEnd.y);
                adjustedQual = adjustedQual.substring(stringStartEnd.x, stringStartEnd.y);
                
                goodBarcodedReads++;
                
                // write entries to file for this taxon
                Taxon readTaxon = barcode.getTaxon(); 
                Tuple<BufferedWriter,String> origTuple = taxaFileMap.get(readTaxon);
                if (origTuple == null) {
                    // add entry
                    Taxon taxonWithAnno = nameToTaxon.get(readTaxon.getName());
                    GeneralAnnotation annotation = taxonWithAnno.getAnnotation();
                    String outFile = outputDir + "/" + projectName + "_" + readTaxon.getName() + "_" +
                            annotation.getTextAnnotation(flowcellField)[0] + "_" + annotation.getTextAnnotation("Lane")[0] +
                            "_" + barcode.getBarcodeString() + ".fastq";
                    taxaFileMap.put(readTaxon,new Tuple<BufferedWriter,String>(null,outFile));
                    origTuple = taxaFileMap.get(readTaxon);
                }

                BufferedWriter bw = origTuple.x;
                if (bw == null) {
                    // First time we've written this file - create writer and store it                  
                    bw = Utils.getBufferedWriter(origTuple.y);
                    taxaFileMap.put(readTaxon, new Tuple<BufferedWriter,String>(bw,origTuple.y));
                }
                bw.write(seqAndQual[0]);
                bw.write("\n");
                bw.write(adjustedSeq);
                bw.write("\n");
                bw.write(seqAndQual[2]);
                bw.write("\n");               
                bw.write(adjustedQual);
                bw.write("\n");   
                if (currentTotal > 1000000) {
                    System.out.println("Reads processed so far: " + allReads);
                    currentTotal = 0;
                }
            } 
            myLogger.info("Summary for "+fastQPath.toString()+"\n"+
                    "Total number of reads in lane=" + allReads +"\n"+
                    "Total number of good barcoded reads=" + goodBarcodedReads+"\n"+
                    "Total number of bar codes not found=" + barCodeNotFound+"\n"+
                    "Total number where first 20 bps of seq has N =" + first20hasN+"\n"+
                    "Total number of read End not found=" + badEndReads+"\n"+
                    "Timing process (parsing, collapsing, and writing readBlock to file)."+"\n"+
                    "Process took " + (System.nanoTime() - time)/1e9 + " seconds.");
        } catch (Exception exc) {
            exc.printStackTrace();
        }
    }
    
    
    public static String[] readAllFastQBlock(BufferedReader bw, int currentRead) throws IOException {
        //consider converting this into a stream of String[]
        String[] result=new String[4];
        try{
            result[0]=bw.readLine();
            result[1] = bw.readLine();
            if(result[0]==null) {
                return null;
            }
            result[2]=bw.readLine();
            result[3] = bw.readLine();
 
            return result;
        } catch (IOException e) {
            e.printStackTrace();
            myLogger.error("Unable to correctly parse the fastQ read nead line: " + currentRead*4
                    + " from fastq file.  Your fastq file may have been corrupted.");
            return null;
        }
    }
 
    // This method determines the start and total length of the passed sequence
    // based on the enzyme and its likely read-end strings.
    // The start position may change based on handling of the overlapping
    // cut sites.  The end site should change based on likely read ends.
    // These values are used in calling method to adjust both the sequence
    // and the quality strings.
    private static Tuple<Integer,Integer> truncateToSecondCutSite(String seq,String enzyme) {
        
        int seqStart = 0;
        // handle overlapping cutsite for ApeKI enzyme
        if (enzyme.equalsIgnoreCase("ApeKI")) {
            if (seq.startsWith("CAGCTGC") || seq.startsWith("CTGCAGC")) {
                seq = seq.substring(3,seq.length());
                seqStart = 3;
            }             
        }
        int indexOfReadEnd = -1;
        // Starts at 20 so don't have overlapping start/end cutsites.
        String shortSeq = seq.substring(20); // substr from 20 to end of sequence
        for (String readEnd: likelyReadEndStrings){
            int indx = shortSeq.indexOf(readEnd);
            if (indx > 0 ) {
                if (indexOfReadEnd < 0 || indx < indexOfReadEnd) {
                    indexOfReadEnd = indx;
                } 
            }
        }
        
        if (indexOfReadEnd < 0) {
            return new Tuple<Integer,Integer>(seqStart,seq.length()-seqStart);
        } else {
            // Add back the 20 from the beginning that we chopped off.
            int tagLen = indexOfReadEnd + 20 + readEndCutSiteRemnantLength;
            return new Tuple<Integer,Integer>(seqStart,tagLen);
        }

    }
    private static void Shutdown(Map<Taxon,Tuple<BufferedWriter,String>> taxaWriterMap) {
        if (taxaWriterMap == null) return;
        try {
            for (Map.Entry<Taxon, Tuple<BufferedWriter,String>> entry: taxaWriterMap.entrySet()) {
                BufferedWriter bw = entry.getValue().x;
                if (bw != null) bw.close();
            }
                
        } catch (Exception e) {
            System.out.println("Problem with shutdown");
            e.printStackTrace();
        }
    }
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
//     public static void main(String[] args) {
//         //GeneratePluginCode.generate(ConvertOldFastqPlugin.class);
//     }
     
    @Override
    public ImageIcon getIcon() {       
        return null;
    }

    @Override
    public String getButtonName() {
        return ("ConvertOldFastqToModernFormat");
        
    }

    @Override
    public String getToolTipText() {       
        return ("Fastq files formatted with multiple taxa and bar codes are transformed into Illumina format with 1 taxa per file and no bar codes in the read sequence");
    }
    
    /**
     * Input file or directory containing FASTQ files in text
     * or gzipped text.
     *      NOTE: Directory will be searched recursively
     *
     * @return Input Directory
     */
    public String inputFile() {
        return inputFile.value();
    }

    /**
     * Set Input Directory. Input file or directory containing
     * FASTQ files in text or gzipped text.
     *      NOTE: Directory will be searched recursively
     *
     * @param value Input Directory
     *
     * @return this plugin
     */
    public ConvertOldFastqToModernFormatPlugin inputFile(String value) {
        inputFile = new PluginParameter<>(inputFile, value);
        return this;
    }

    /**
     * Key file listing barcodes distinguishing the samples
     *
     * @return Key File
     */
    public String keyFile() {
        return keyFile.value();
    }

    /**
     * Set Key File. Key file listing barcodes distinguishing
     * the samples
     *
     * @param value Key File
     *
     * @return this plugin
     */
    public ConvertOldFastqToModernFormatPlugin keyFile(String value) {
        keyFile = new PluginParameter<>(keyFile, value);
        return this;
    }

    /**
     * Name for this project.  Project name becomes part of
     * the newly created fastq file name
     *
     * @return Project File
     */
    public String projectName() {
        return projectName.value();
    }

    /**
     * Set Project File. Name for this project.  Project name
     * becomes part of the newly created fastq file name
     *
     * @param value Project File
     *
     * @return this plugin
     */
    public ConvertOldFastqToModernFormatPlugin projectName(String value) {
        projectName = new PluginParameter<>(projectName, value);
        return this;
    }

    /**
     * Output directory where new fastQ files will be written
     *
     * @return Output Directory
     */
    public String outputDir() {
        return outputDir.value();
    }

    /**
     * Set Output Directory. Output directory where new fastQ
     * files will be written
     *
     * @param value Output Directory
     *
     * @return this plugin
     */
    public ConvertOldFastqToModernFormatPlugin outputDir(String value) {
        outputDir = new PluginParameter<>(outputDir, value);
        return this;
    }
    
    /**
     * Enzyme used to create the GBS library
     *
     * @return Enzyme
     */
    public String enzyme() {
        return myEnzyme.value();
    }

    /**
     * Set Enzyme. Enzyme used to create the GBS library
     *
     * @param value Enzyme
     *
     * @return this plugin
     */
    public ConvertOldFastqToModernFormatPlugin enzyme(String value) {
        myEnzyme = new PluginParameter<>(myEnzyme, value);
        return this;
    }

}

