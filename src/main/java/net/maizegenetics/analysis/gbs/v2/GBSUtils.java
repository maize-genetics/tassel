/**
 * 
 */
package net.maizegenetics.analysis.gbs.v2;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.analysis.gbs.Barcode;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.Utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Ordering;

/**
 * This class contains methods and constants used by various classes in the GBSv2
 * pipeline.
 * @author lcj34
 *
 */
public class GBSUtils {

    private static final Logger myLogger = LogManager.getLogger(GBSUtils.class);
    public static final String inputFileGlob="glob:*{.fq,fq.gz,fastq,fastq.txt,fastq.gz,fastq.txt.gz,_sequence.txt,_sequence.txt.gz}";
    public static final String sampleNameField="FullSampleName";
    public static final String flowcellField="Flowcell";
    public static final String laneField="Lane";
    public static final String barcodeField="Barcode";
    public static final String tissueNameField = "Tissue";
    public static final String fileNameField = "FileName";
    
    private GBSUtils() {
    }
    
    /**
     * Method for reading FastQ four line structure, and returning a string array with [sequence, qualityScore]
     */
    public static String[] readFastQBlock(BufferedReader bw, int currentRead) throws IOException {
        //consider converting this into a stream of String[]
        String[] result=new String[2];
        try{
            bw.readLine();
            result[0]=bw.readLine();
            bw.readLine();
            result[1]=bw.readLine();
            if(result[0]==null) {
                return null;
            }
            return result;
        } catch (IOException e) {
            e.printStackTrace();
            myLogger.error("Unable to correctly parse the sequence and quality score near line: " + currentRead*4
                    + " from fastq file.  Your fastq file may have been corrupted.");
            return null;
        }
    }
 
    /**
     * Method for reading FastQ four line structure, and returning a string array with [sequence, qualityScore]
     */
    public static String[] readDeMultiPlexFastQBlock(BufferedReader bw, int currentRead) throws IOException {
        //consider converting this into a stream of String[]
        String[] result=new String[3];
        try{
            // Grab the barcode from the first line of the fastq sequence
            String barCode = bw.readLine();
            if (barCode == null) {
                return null;
            }
            int index = barCode.lastIndexOf(":");
            result[2] = barCode.substring(index+1);

            // sequence
            result[0]=bw.readLine();
            bw.readLine();//quality header - thrown away
            // Second entry is the quality score
            result[1]=bw.readLine();
            return result;
        } catch (IOException e) {
            e.printStackTrace();
            myLogger.error("Unable to correctly parse the sequence and quality score near line: " + currentRead*4
                    + " from fastq file.  Your fastq file may have been corrupted.");
            return null;
        }
    }
    /**
     * Method for reading FastQ four line structure, and returning a string array with [sequence, qualityScore]
     */
    public static int determineQualityScoreBase(Path fastqFile) throws IOException {
        try{BufferedReader bw = Utils.getBufferedReader(fastqFile.toString());
            int headerParts=bw.readLine().split(":").length;
            int base=(headerParts<5)?64:33;
            myLogger.info(fastqFile.toString()+": Quality score base:"+base);
            return base;
        } catch (IOException e) {
            e.printStackTrace();
            myLogger.error("Unable to correctly parse the quality score base from fastq file.  " +
                    "Your fastq file may have been corrupted.");
            return 0;
        }
    }
    
    /**
     * Returns an annotated taxaList based on a Keyfile for GBS
     * @param keyPath
     * @param fastQpath
     * @return
     */
    public static ArrayList<Taxon> getLaneAnnotatedTaxaList(Path keyPath, Path fastQpath) {
        String[] filenameField = fastQpath.getFileName().toString().split("_");
        ArrayList<Taxon> annoTL;
        if (filenameField.length == 3) {
            annoTL = TaxaListIOUtils.readTaxaAnnotationFileAL(keyPath.toAbsolutePath().toString(), sampleNameField,
                    ImmutableMap.of(flowcellField, filenameField[0], laneField, filenameField[1])); 
        } else if (filenameField.length == 4) {
            annoTL = TaxaListIOUtils.readTaxaAnnotationFileAL(keyPath.toAbsolutePath().toString(),sampleNameField,
                    ImmutableMap.of(flowcellField, filenameField[0], laneField, filenameField[2]));
        }
        else if (filenameField.length == 5) {
            annoTL = TaxaListIOUtils.readTaxaAnnotationFileAL(keyPath.toAbsolutePath().toString(),sampleNameField,
                    ImmutableMap.of(flowcellField, filenameField[1], laneField, filenameField[3]));
        } else {
            myLogger.error("Error in parsing file name: " + fastQpath.toString());
            myLogger.error("   The filename does not contain either 3, 4, or 5 underscore-delimited values.");
            myLogger.error("   Expect: flowcell_lane_fastq.txt.gz OR flowcell_s_lane_fastq.txt.gz OR code_flowcell_s_lane_fastq.txt.gz");
            return null;
        }
        return annoTL;
    }
    
    /**
     * Produces a trie for sorting the read
     * @param taxaList the taxaList of the current flowcell lanes that is annotated with barcode information
     * @param masterTaxaList  the mastertaxaList provides the taxaIndex
     * @param myEnzyme
     * @return Barcode trie for examining the prefixes
     */

    public static BarcodeTrie initializeBarcodeTrie(ArrayList<Taxon> taxaList, TaxaList masterTaxaList, 
             EnzymeList.Enzyme myEnzyme){
        BarcodeTrie aTrie=new BarcodeTrie();
        for (Taxon taxon : taxaList) {
            int masterIndex=masterTaxaList.indexOf(taxon.getName());
            GeneralAnnotation annotation = taxon.getAnnotation();
            String[] myTissues = annotation.getTextAnnotation("Tissue");
            // Tissue should be stored as annotation against the taxon
            Barcode theBC = null;
            theBC = new Barcode(annotation.getTextAnnotation(barcodeField)[0], myEnzyme.initialCutSiteRemnant(), taxon.getName(),
                    masterIndex,annotation.getTextAnnotation(flowcellField)[0],annotation.getTextAnnotation("Lane")[0]);
            aTrie.addBarcode(theBC);
        }
        return aTrie;
    }  
    
    // THis one sends in a  tissue list - is called from the RNASeq pipeline
    public static BarcodeTrie initializeBarcodeTrie(ArrayList<Taxon> taxaList, TaxaList masterTaxaList, 
            ArrayList masterTissueList, EnzymeList.Enzyme myEnzyme){
        BarcodeTrie aTrie=new BarcodeTrie();
        for (Taxon taxon : taxaList) {
            int masterIndex=masterTaxaList.indexOf(taxon.getName());
            GeneralAnnotation annotation = taxon.getAnnotation();
            String[] myTissues = annotation.getTextAnnotation("Tissue");
            // Tissue should be stored as annotation against the taxon
            Barcode theBC = null;
            if (myTissues.length > 0) {
                int masterTissueIndex = masterTissueList.indexOf(myTissues[0]);
                // keyfile had Tissue column: tissues were added to taxon annotations
                theBC = new Barcode(annotation.getTextAnnotation(barcodeField)[0], myEnzyme.initialCutSiteRemnant(), taxon.getName(),
                        masterIndex,myTissues[0], masterTissueIndex, 
                        annotation.getTextAnnotation(flowcellField)[0],annotation.getTextAnnotation("Lane")[0]);
            } else { // no tissue variables in the taxon annotations
                theBC = new Barcode(annotation.getTextAnnotation(barcodeField)[0], myEnzyme.initialCutSiteRemnant(), taxon.getName(),
                        masterIndex,annotation.getTextAnnotation(flowcellField)[0],annotation.getTextAnnotation("Lane")[0]);
            }
            aTrie.addBarcode(theBC);
        }
        return aTrie;
    }   
    /**
     * Produces a list of fastq files that are represented by the plugin's keyfile
     * @param directoryFiles:  List of all the files in the directory
     * @return filesToProcess:  List of only those files that should be processed
     */
    public static List<Path> culledFiles(List<Path>directoryFiles,Path keyFile ) {
        
        List<Path> filesToProcess = new ArrayList<Path>();
        // Get map  of flowcell/lanes from the key file
        String keyFileName = keyFile.toString();
        ListMultimap<String, String> keyFileValues = parseKeyfileIntoMap(keyFileName); 
        if (keyFileValues.isEmpty()) return filesToProcess; // no entries

        // for each file in the directory, check if the flowcell and lane are represented 
        // The directoryFile list is in alphabetical order.  It is quicker to run a non-parallel
        // stream and skip sorting than run with parallel and have to sort at the end (entries
        // in filesToProcess are not in alphabetical order when parallelStream is used). 
        // Alphabetical order is necessary to ensure consistency of tags removed by 
        // "removeTagsWithoutReplication" when multiple runs are performed.
        directoryFiles.stream()
        .forEach(directoryFile -> {             
                String[] filenameField = directoryFile.getFileName().toString().split("_");
            if (filenameField.length == 3) {
               if (keyFileValues.containsEntry(filenameField[0],filenameField[1])) {
                   filesToProcess.add(directoryFile);
               }
            } else if (filenameField.length == 4) {
                if (keyFileValues.containsEntry(filenameField[0],filenameField[2])) {
                   filesToProcess.add(directoryFile);
                }
            }
            else if (filenameField.length == 5) {
                if (keyFileValues.containsEntry(filenameField[1],filenameField[3])) {
                   filesToProcess.add(directoryFile);
                }
            }
            else {
                myLogger.error("Error in parsing file name: " + directoryFile.toString());
                myLogger.error("   The filename does not contain either 3, 4, or 5 underscore-delimited values.");
                myLogger.error("   Expect: flowcell_lane_fastq.txt.gz OR flowcell_s_lane_fastq.txt.gz OR code_flowcell_s_lane_fastq.txt.gz");
                myLogger.error("   " + directoryFile.toString() + " will not be processed.");                
            }
        });             
        return filesToProcess; 
    }
    
    /**
     * Parses a tab-delimited keyFile storing the flow cell and lane values into a multimap.
     * The flow cell is the key, which may have multiple associated lanes.
     *
     * @param
     * @return
     */
    public static ListMultimap<String, String> parseKeyfileIntoMap(String fileName) {
        if (fileName == null) {
            return null;
        }
        ImmutableListMultimap.Builder<String, String> mMap = new ImmutableListMultimap.Builder<String, String>()
                .orderKeysBy(Ordering.natural()).orderValuesBy(Ordering.natural());
        try {
            BufferedReader fileIn = Utils.getBufferedReader(fileName, 1000000);
            fileIn.mark(1 << 16);
            String line = fileIn.readLine();
            int indexOfFlowcell = 0, indexOfLane = 0;
            //parse headers
            if (line.contains(flowcellField)) {
                int idx = 0;
                for (String header : line.split("\\t")) {
                    if (header.equals(flowcellField)) {
                        indexOfFlowcell = idx;
                    }
                    if (header.equals(laneField)) {
                        indexOfLane = idx;
                    }
                    idx++;
                }
            } else {
                fileIn.reset();
            }
            // create list of flowcells and lanes
            while ((line = fileIn.readLine()) != null) {
                String[] myString = line.split("\\t");
                String myFlowCell = myString[indexOfFlowcell];
                String myLane = myString[indexOfLane];
                mMap.put(myFlowCell,myLane);
            }
        } catch (Exception e) {
            System.err.println("Error in Reading Parsing Key File:" + fileName);
            e.printStackTrace();
        }
        return mMap.build();
    }
}
