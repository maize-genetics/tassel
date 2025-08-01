package net.maizegenetics.analysis.gbs.v2;
/**
 * Tests from this file were used to compare results from GBSv2 pipeline 
 * to old GBS pipeline as part of TAS-865.
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TagDataWriter;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;

import org.junit.Test;

public class GBSv2BiologyCompareTest {
    static String depthDiffFile = GBSConstants.GBS_TEMP_DIR + "tagsInBothWithDepths.txt";
    static String depthDiffFileWithLen = GBSConstants.GBS_TEMP_DIR + "tagsInBothWithDepthsAndLen.txt";
    static String diffTagDepths = GBSConstants.GBS_TEMP_DIR + "differingTagDepths.txt";
    static String diffTagDepthsWithLen = GBSConstants.GBS_TEMP_DIR + "differingTagDepthsWithLen.txt";


    @Test
    public void compareTagsTest () throws Exception {
        System.out.println("Begin compareTagsTest\n");
        ArrayList<String> origGBSTags = getTagsFromOrigGBS();
        ArrayList<String> barcodes = getBarCodesFromKeyFile(GBSConstants.GBS_TESTING_KEY_FILE);

        System.out.println(Paths.get("Creating DataBase: " + GBSConstants.GBS_GBS2DB_FILE).toAbsolutePath().toString());
        try{
            Files.deleteIfExists(Paths.get(GBSConstants.GBS_GBS2DB_FILE));
        } catch (IOException e) {
            e.printStackTrace();
        }
        long time=System.nanoTime();
        System.out.println();

        new GBSSeqToTagDBPlugin()
        .enzyme("ApeKI")
        .inputDirectory(GBSConstants.GBS_INPUT_DIR)
        .outputDatabaseFile(GBSConstants.GBS_GBS2DB_FILE)
        .keyFile(GBSConstants.GBS_TESTING_KEY_FILE)
        .kmerLength(64)
        .batchSize(1)
        .minimumKmerLength(20) // added - default is 20
        .minKmerCount(5) // added - default = 10
        .minimumQualityScore(0) // lcj - was 20, kept lower to include more tags
        .performFunction(null);

        TagDataWriter tdw =  new TagDataSQLite(GBSConstants.GBS_GBS2DB_FILE);

        // here - print tags with total depth for comparison
        Set<Tag> gbsv2Tags = tdw.getTags();
        ArrayList<String> gbsv2TagSequence = new ArrayList<>();; 
        gbsv2Tags.stream()
        .forEach(tag -> {
            gbsv2TagSequence.add(tag.sequence());
        });
        Collections.sort(gbsv2TagSequence);
        System.out.println("CompareTagsTest: gbsv2Tags size: " + gbsv2Tags.size() + ", gbsv2TagSequence size: " + gbsv2TagSequence.size());

        List<String> tagsInOldNotNew = origGBSTags.stream()
                .filter(tagSequence -> {                            
                    // find tags in original GBS   
                    return (!gbsv2Tags.contains(TagBuilder.instance(tagSequence).build()));
                })
                .collect(Collectors.toList());

        List<Tag> tagsInNewNotOld = gbsv2Tags.stream()
                .filter(tag -> {
                    return (!origGBSTags.contains(tag.sequence()));
                })
                .collect(Collectors.toList());

        String outputFile = GBSConstants.GBS_TEMP_DIR +"tagsNotInGBSv2.txt";
        // write to a file I can look at
        BufferedWriter fileWriter = null;
        StringBuilder strB = new StringBuilder();    
        tagsInOldNotNew.stream().forEach(entry -> {
            strB.append(entry);
            strB.append("\n");
        });                
        try {  
            fileWriter = new BufferedWriter(new FileWriter(outputFile));
            fileWriter.write(strB.toString());
        }
        catch(IOException e) {
            System.out.println("Caught Exception in compareTagsTest writing first file");
            System.out.println(e);
        }
        fileWriter.close();

        // create 2nd file
        StringBuilder strB2 = new StringBuilder();
        outputFile = GBSConstants.GBS_TEMP_DIR + "tagsNotInOrigGBS.txt";
        BufferedWriter fileWriter2 = null;    
        tagsInNewNotOld.stream().forEach(tag -> {
            strB2.append(tag.sequence());
            strB2.append("\n");
        });                
        try {  
            fileWriter2 = new BufferedWriter(new FileWriter(outputFile));
            fileWriter2.write(strB2.toString());
        }
        catch(IOException e) {
            System.out.println("Caught Exception in compareTagsTest writing second file");
            System.out.println(e);
        }
        fileWriter2.close();

        StringBuilder strB3 = new StringBuilder();
        outputFile = GBSConstants.GBS_TEMP_DIR + "gbsv2TagsWithDepth.txt";
        BufferedWriter fileWriter3 = null;
        System.out.println("compareTagsTest: gbsv2TagSequence size: " + gbsv2TagSequence.size());
        Map<Tag,Integer> gbsv2TagsWithDepth = tdw.getTagsWithDepth(0);
        gbsv2TagsWithDepth.entrySet().stream().forEach(entry -> {
            strB3.append(entry.getKey().sequence());
            strB3.append("\t");
            strB3.append(entry.getKey().seqLength());
            strB3.append("\t");
            strB3.append(entry.getValue());
            strB3.append("\n");
        });                
        try {  
            fileWriter3 = new BufferedWriter(new FileWriter(outputFile));
            fileWriter3.write(strB3.toString());
        }
        catch(IOException e) {
            System.out.println("Caught Exception in compareTagsTest writing third file");
            System.out.println(e);
        }
        fileWriter3.close();
        
        System.out.println("\n\nFinished CompareTagsTest: tags in gbsV2: " + gbsv2Tags.size()
                + " tags in orig GBS: " + origGBSTags.size());
        System.out.println("    tags in old GBS not in v2: " + tagsInOldNotNew.size());
        System.out.println("    tags in v2 GBS not in old: " + tagsInNewNotOld.size());  
        System.out.println("\nNOTE: Reasons for differences in tag counts polyA additions in oldGBS "
                + "\n     and differences in when/how tags are removed from the db at the end");

    }

    @Test
    public void compareTagDepths() throws Exception {
        String oldGBSMergedFile = GBSConstants.GBS_EXPECTED_DIR + "MergeMultipleTagCountPlugin/Merged_Tag_CountsConverted.txt";
        String gbsv2WithDepthFile = GBSConstants.GBS_TEMP_DIR + "gbsv2TagsWithDepth.txt";
        Map<String, Tuple<Integer,Integer>> oldGBSTagsDepth =getTagDepthFromFile(oldGBSMergedFile);
        Map<String, Tuple<Integer,Integer>> v2TagsDepth = getTagDepthFromFile(gbsv2WithDepthFile);
        
        // create file with tag, olddepth, newdepth
        BufferedWriter fileWriter = null;
        StringBuilder strB = new StringBuilder(); 
        strB.append("Tag\tOldLength\tOldDepth\tV2Length\tV2Depth\n");
        v2TagsDepth.entrySet().stream().forEach(entry -> {
            String v2Tag = entry.getKey();
            if (oldGBSTagsDepth.containsKey(v2Tag)) {
                strB.append(v2Tag);
                strB.append("\t");
                strB.append(oldGBSTagsDepth.get(v2Tag).x);
                strB.append("\t");
                strB.append(oldGBSTagsDepth.get(v2Tag).y);
                strB.append("\t");
                strB.append(v2TagsDepth.get(v2Tag).x);
                strB.append("\t");
                strB.append(v2TagsDepth.get(v2Tag).y);
                strB.append("\n");
            }
        });
             
        try {  
            fileWriter = new BufferedWriter(new FileWriter(depthDiffFileWithLen));
            fileWriter.write(strB.toString());
        }
        catch(IOException e) {
            System.out.println("Caught Exception in compareTagsTest writing first file");
            System.out.println(e);
        }
        fileWriter.close();
        
        // create file that only has tag sequences with depths that are different.
        // file shows tag, depthOld, depthV2
        createDepthDiffFile(diffTagDepthsWithLen);
        System.out.println("\n\nTest compareTagDepths complete, file written to " + depthDiffFileWithLen);
        
    }
    
    @Test
    public void mapTaxonToTag() throws Exception {
        Map<String,String> taxonToBarcode = getBarCodeTaxonFromKeyFile(GBSConstants.GBS_TESTING_KEY_FILE);
        String oldGBSTaxonNames = GBSConstants.GBS_EXPECTED_DIR + "FastqToTagCountPlugin/fastqToTagCount2_taxonNames.txt";
        String outputFile = GBSConstants.GBS_TEMP_DIR + "oldGBSTaxonNameToBarcode_126x32.txt";
        
        List<String> origGBSTaxonNames = new ArrayList<>();
        
        try {
            BufferedReader fileIn = Utils.getBufferedReader(oldGBSTaxonNames, 1000000);
            // First value is Sequence string, second is length, third is count
            String line; // no header to skip
            while ((line = fileIn.readLine()) != null) {
                origGBSTaxonNames.add(line);
            }
        } catch (Exception e) {
            System.err.println("Error reading oldGBSTaxonNames file:" + oldGBSTaxonNames);
            e.printStackTrace();
        }       
        
        BufferedWriter fileWriter = null;
        StringBuilder strB = new StringBuilder();   
        
        for (String name:origGBSTaxonNames) {
            String barcode = taxonToBarcode.get(name);
            strB.append(name);
            strB.append("\t");
            strB.append(barcode);
            strB.append("\n");
        }
      
        try {  
            fileWriter = new BufferedWriter(new FileWriter(outputFile));
            fileWriter.write(strB.toString());
        }
        catch(IOException e) {
            System.out.println("Caught Exception in mapTaxonToTag writing  file");
            System.out.println(e);
        }
        System.out.println("mapTaxonToTag test complete: file written to " + outputFile);
        fileWriter.close();       
                
    }
    
    private static ArrayList<String> getTagsFromOrigGBS() {
        ArrayList<String> origGBSTags = new ArrayList<>();
        // the file below is the ExpectedREsults/GBS/Chr9_10-20000000/MergeMultipleTagCountPlugin/merged_Tag_Counts.cnt
        // file converted to tab-delimited file via BinaryToTextPlugin (in GBS folder)
        String fileName = GBSConstants.GBS_EXPECTED_DIR + "MergeMultipleTagCountPlugin/Merged_Tag_CountsConverted.txt";
        try {
            BufferedReader fileIn = Utils.getBufferedReader(fileName, 1000000);
            // First value is Sequence string, second is length, third is count
            String line = fileIn.readLine(); // read/skip header, which is number of lines
            while ((line = fileIn.readLine()) != null) {
                String[] tokens = line.split("\\t");
                if (tokens.length != 3) {
                    System.err.println("Error in SNP Position QualityScore file format:" + fileName);
                    System.err.println("Expecting tab-delimited file with 2 integer and 1 float value per row "
                            + " with header values CHROM POS QUALITYSCORE");
                }
                origGBSTags.add(tokens[0]);
            }
        } catch (Exception e) {
            System.err.println("Error in original gbs File:" + fileName);
            e.printStackTrace();
        }
        return origGBSTags;    
    }

    private static ArrayList<String> getBarCodesFromKeyFile(String keyfile) {
        ArrayList<String> barcodes = new ArrayList<>();
        try {
            BufferedReader fileIn = Utils.getBufferedReader(keyfile, 1000000);
            // First value is Sequence string, second is length, third is count
            String line = fileIn.readLine(); // read/skip header, which is number of lines
            while ((line = fileIn.readLine()) != null) {
                String[] tokens = line.split("\\t");
                //System.out.println("LCJ -adding  barcode: " + tokens[2]);
                barcodes.add(tokens[2]);
            }
        } catch (Exception e) {
            System.err.println("Error in key file:" + keyfile);
            e.printStackTrace();
        }
        return barcodes;    
    }
 
    private static Map<String, String> getBarCodeTaxonFromKeyFile(String keyFile) {
        Map<String, String> barCodeTaxonMap = new HashMap<String,String>();
        try {
            BufferedReader fileIn = Utils.getBufferedReader(keyFile,1000000);
            String line = fileIn.readLine(); // skip header
            while ((line = fileIn.readLine()) != null) {
                String[] tokens = line.split("\\t");
                barCodeTaxonMap.put(tokens[18], tokens[2]);
                       
            }
        } catch (Exception ex) {
            System.err.println("Error in key file:" + keyFile);
            ex.printStackTrace();
        }
        return barCodeTaxonMap;
    }
    private static Map<String,Tuple<Integer,Integer>> getTagDepthFromFile(String fileName) {
        Map<String,Tuple<Integer,Integer>> tagDepthMap = new HashMap<String,Tuple<Integer,Integer>>();
        try {
            BufferedReader fileIn = Utils.getBufferedReader(fileName, 1000000);
            // First value is Sequence string, second is length, third is count
            String line; 
            while ((line = fileIn.readLine()) != null) {
                String[] tokens = line.split("\\t"); 
                tagDepthMap.put(tokens[0],new Tuple<>(Integer.parseInt(tokens[1]),Integer.parseInt(tokens[2])));
            }
        } catch (Exception e) {
            System.err.println("getTagDepthFromFile: Error in original gbs File:" + fileName);
            e.printStackTrace();
        }
        return tagDepthMap;           
    }
    
    private static void createDepthDiffFile(String outfile) throws Exception{
        StringBuilder strB = new StringBuilder();
        try {
            BufferedReader fileIn = Utils.getBufferedReader(depthDiffFileWithLen, 1000000);
            // First value is Sequence string, second is length, third is count
            String line = fileIn.readLine();  // skip first line
            strB.append("Tag\tLenOrigGBS\tDepthOrigGBS\tLenV2\tDepthV2");
            while ((line = fileIn.readLine()) != null) {
                String[] tokens = line.split("\\t");
                if (!tokens[2].equals(tokens[4])) {
                    strB.append("\n");
                    strB.append(tokens[0]);
                    strB.append("\t");
                    strB.append(tokens[1]);
                    strB.append("\t");
                    strB.append(tokens[2]);
                    strB.append("\t");
                    strB.append(tokens[3]);
                    strB.append("\t");
                    strB.append(tokens[4]);
                }
            }
        } catch (Exception e) {
            System.err.println("createDepthDiffFile: Error in File:" + depthDiffFileWithLen);
            e.printStackTrace();
        }
        
        BufferedWriter fileWriter = null;
        try {  
            fileWriter = new BufferedWriter(new FileWriter(outfile));
            fileWriter.write(strB.toString());
        }
        catch(IOException e) {
            System.out.println("Caught Exception in writing depthDiffFile");
            System.out.println(e);
        }
        fileWriter.close();
    }  
    
    @Test
    public void verifyTag () throws Exception {
        // This test verifies the GetTagSequenceFromDBPlugin will return
        // a tab-delimited file.  
        String outFile = GBSConstants.GBS_TEMP_DIR + "tagInDB.txt";
        new GetTagSequenceFromDBPlugin()
        .tagSequence("CAGCAAAAAAAAACGAAGAAGAGATGTTCCTCTATATCTATACTATACTTAAAAACAGTGTCTT")
        .inputDB(GBSConstants.GBS_GBS2DB_FILE)
        .outputFile(outFile)
        .performFunction(null);
    }
}
