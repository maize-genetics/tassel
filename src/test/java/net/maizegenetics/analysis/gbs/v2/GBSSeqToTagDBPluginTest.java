/*
 * FastqToTagCountPluginTest
 */
package net.maizegenetics.analysis.gbs.v2;

import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.map.*;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TagData;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Tuple;

import org.junit.Test;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import java.util.Random;
import java.util.Set;

/**
 *
 * @author Ed Buckler
 */
public class GBSSeqToTagDBPluginTest {
    private final static String outTagFasta=GBSConstants.GBS_TEMP_DIR + "tagsForAlign.fa.gz";
    private final static String inTagSAM=GBSConstants.GBS_TEMP_DIR + "tagsForAlign910auto.sam";

    public GBSSeqToTagDBPluginTest() {
    }

    /**
     * Test of performFunction method, of class FastqToTagCountPlugin.
     */
    @Test
    public void testGBSSeqToTagDBPlugin() {
        LoggingUtils.setupDebugLogging();
        System.out.println(Paths.get(GBSConstants.GBS_GBS2DB_FILE).toAbsolutePath().toString());
//        try{
//            Files.deleteIfExists(Paths.get(GBSConstants.GBS_GBS2DB_FILE));
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
        long time=System.nanoTime();
        System.out.println();

        new GBSSeqToTagDBPlugin()
                .enzyme("ApeKI")
                .inputDirectory(GBSConstants.GBS_INPUT_DIR)
                //.inputDirectory("/Users/lcj34/notes_files/gbsv2/debug_problems/shuzhen_noGBSSeqData/fastq")
                //.outputDatabaseFile("/Users/lcj34/notes_files/gbsv2/debug_problems/shuzhen_noGBSSeqData/gbsv2.db")
                .outputDatabaseFile(GBSConstants.GBS_GBS2DB_FILE)
                //.outputDatabaseFile("/Volumes/Samsung_T1/gbsv2.db")// - can't write to external drive!
                .keyFile(GBSConstants.GBS_TESTING_KEY_FILE)
                //.keyFile("/Users/lcj34/notes_files/gbsv2/debug_problems/shuzhen_noGBSSeqData/Key1-5.txt")
                .kmerLength(64)
                //.minimumKmerLength(64) // LCJ - this is normally not here!  COmment it back out
                .minKmerCount(5)
                .minimumQualityScore(20)
                .deleteOldData(false)
                //.maximumMapMemoryInMb(5500)
                .performFunction(null);
        System.out.printf("TotalTime %g sec%n", (double) (System.nanoTime() - time) / 1e9);
        //String outName= GeneralConstants.TEMP_DIR+"TaxaTest.db";
        //TagsByTaxaHDF5Builder.create(, (Map<Tag,TaxaDistribution>)ds.getData(0).getData(),null);
    }
    
    @Test
    public void testKeepOldData() {
        // THis should be committed against tas 1120
        LoggingUtils.setupDebugLogging();
        System.out.println(Paths.get(GBSConstants.GBS_GBS2DB_FILE).toAbsolutePath().toString());
        // start with a clean db
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
                .minKmerCount(5)
                .minimumQualityScore(20)
                .deleteOldData(true)
                .performFunction(null);
        System.out.printf("TotalTime %g sec%n", (double) (System.nanoTime() - time) / 1e9);
        //String outName= GeneralConstants.TEMP_DIR+"TaxaTest.db";
        //TagsByTaxaHDF5Builder.create(, (Map<Tag,TaxaDistribution>)ds.getData(0).getData(),null);
        
        
     // Verify total number of tags:
        TagData tagData =new TagDataSQLite(GBSConstants.GBS_GBS2DB_FILE);
        Set<Tag> gbsv2Tags1 = tagData.getTags();
        try {
            ((TagDataSQLite) tagData).close();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
        System.out.println("\nRunning GBSSeqToTagDBPLugin second time, do NOT delete data: ... ");
        // LCJ - then re-run and verify the size of the tags in the table
        new GBSSeqToTagDBPlugin()
        .enzyme("ApeKI")
        .inputDirectory(GBSConstants.GBS_INPUT_DIR)
        .outputDatabaseFile(GBSConstants.GBS_GBS2DB_FILE)
        .keyFile(GBSConstants.GBS_TESTING_KEY_FILE)
        .kmerLength(64)
        .minKmerCount(5)
        .minimumQualityScore(20)
        .deleteOldData(false)
        .performFunction(null);
        System.out.printf("TotalTime %g sec%n", (double) (System.nanoTime() - time) / 1e9);
        
        Set<Tag> gbsv2Tags2 = tagData.getTags();
        // Number of tags should be the same
        // The actual test performed was to run the plugin once.
        // Then run again with delete-old-data = false.  In the plugin code,
        // after the existing tags/taxadist were pulled from the db, run
        // "removeTagsWithoutReplication(tagCntMap);", then return.  Verify
        // the number of tags in the DB matches the number of tags existing after
        // the first run.  This verified that depths were maintained.  A previous
        // bug created the new tag/taxadist map without correct depths, and
        // the "removeTagsWithoutReplication() method deleted them.
        // Those test lines are preserved in the source file, but are commented out.
        assertEquals(gbsv2Tags1.size(),gbsv2Tags2.size());

    }

    @Test
    public void testTagExportPlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        System.out.println("Running testTagExportPlugin");
        new TagExportToFastqPlugin()
                .inputDB(GBSConstants.GBS_GBS2DB_FILE)
                .outputFile(outTagFasta)
                .performFunction(null);

    }

    @Test
    public void testSAMImportPlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        System.out.println("Running testTagExportPlugin");
        new SAMToGBSdbPlugin()
                .gBSDBFile(GBSConstants.GBS_GBS2DB_FILE)
                //.sAMInputFile(GBSConstants.GBS_EXPECTED_BOWTIE_SAM_FILE)
                .sAMInputFile(inTagSAM)
                .performFunction(null);
    }

    @Test
    public void GBSSeqToTagDBPluginAppendTest() throws Exception {
        // Testing duplicate tagTagIDMap when GBSSEq run twice in a row (run this test twice in a row!)
        System.out.println(Paths.get(GBSConstants.GBS_GBS2DB_FILE).toAbsolutePath().toString());
        try{
            Files.deleteIfExists(Paths.get(GBSConstants.GBS_GBS2DB_FILE));
        } catch (IOException e) {
            e.printStackTrace();
        }
        long time=System.nanoTime();
        System.out.println("Running GBSSeqToTagDBPluginAppendTest");

        new GBSSeqToTagDBPlugin()
        .enzyme("ApeKI")
        .inputDirectory(GBSConstants.GBS_INPUT_DIR)
        .outputDatabaseFile(GBSConstants.GBS_GBS2DB_FILE)
        .keyFile(GBSConstants.GBS_TESTING_KEY_FILE)
        .kmerLength(64)
        .minKmerCount(5)
        .minimumQualityScore(20)
        .performFunction(null);

        // Grab existing data from db, append to tagCntMap
        TagData tdw =new TagDataSQLite(GBSConstants.GBS_GBS2DB_FILE);
        Set<Tag> firstRunTags = tdw.getTags();
        Map<Tag, TaxaDistribution> firstTDM = tdw.getAllTagsTaxaMap(); 
        ((TagDataSQLite)tdw).close();

        BufferedWriter fileWriter = null;
        try { 
            fileWriter = new BufferedWriter(new FileWriter(GBSConstants.GBS_TEMP_DIR + "dbTagsFirstRun.txt"));
            String headers = "sequence\ttotalDepth\tTaxaDist_ToString";
            fileWriter.write(headers);
            for (Map.Entry<Tag, TaxaDistribution> entry : firstTDM.entrySet()) {
                fileWriter.write("\n");
                Tag tag = entry.getKey();
                TaxaDistribution td = entry.getValue();
                fileWriter.write(tag.sequence());
                fileWriter.write("\t");
                String totalDepth = "totalDepth=" + td.totalDepth(); // This prints depth correctly
                //fileWriter.write(td.totalDepth()); // gave weird output in file
                fileWriter.write(totalDepth);
                fileWriter.write("\t");
                fileWriter.write(td.toString());
            }
            fileWriter.close();
        }catch(IOException e) {
            System.out.println(e);
        }
        fileWriter.close();

        System.out.println("Calling GBSSeqToTagDBPlugin 2nd time to append data.");
        // Run again, this time append tags.  This uses the same fastq files, so
        // we expect the same tags, but with increased total depth at all kept positions
        new GBSSeqToTagDBPlugin()
        .enzyme("ApeKI")
        .inputDirectory(GBSConstants.GBS_INPUT_DIR)
        .outputDatabaseFile(GBSConstants.GBS_GBS2DB_FILE)
        .keyFile(GBSConstants.GBS_TESTING_KEY_FILE)
        .kmerLength(64)
        .minKmerCount(5)
        .minimumQualityScore(20)
        .performFunction(null);

        // Grab existing data from db, store in map
        tdw =new TagDataSQLite(GBSConstants.GBS_GBS2DB_FILE);
        Set<Tag> secondRunTags = tdw.getTags();
        Map<Tag, TaxaDistribution> secondTDM = tdw.getAllTagsTaxaMap(); 
        ((TagDataSQLite)tdw).close();  //todo autocloseable should do this but it is not working.

        // Verify number of tags are the same.  BEcause we are using the
        // same files it will have the same tags.  The tags/taxa are the same,
        // the depths at each tag should be greater after the second run that appends.
        assertEquals(firstRunTags.size(),secondRunTags.size());

        // Write files to the database so tags can be compared between files
        // and differences analysed.
        try { 
            fileWriter = new BufferedWriter(new FileWriter(GBSConstants.GBS_TEMP_DIR + "dbTagsAfterAppendRun.txt"));
            String headers = "sequence\ttotalDepth\tTaxaDist_ToString";
            fileWriter.write(headers);
            for (Map.Entry<Tag, TaxaDistribution> entry : secondTDM.entrySet()) {
                fileWriter.write("\n");
                Tag tag = entry.getKey();
                TaxaDistribution td = entry.getValue();
                fileWriter.write(tag.sequence());
                fileWriter.write("\t");
                String totalDepth = "totalDepth=" + td.totalDepth(); 
                fileWriter.write(totalDepth);
                fileWriter.write("\t");
                fileWriter.write(td.toString());
            }
            fileWriter.close();
        }catch(IOException e) {
            System.out.println(e);
        }
        fileWriter.close();

        // LCJ - I found "short" tags (ie tags with length < maxTagL size) never
        // increased their tag depths value.  Need to understand why this is.  The
        // values were always the same before and after append, the counts were
        // often greater than minTagCount so that wasn't the issue.

        //      Set<Tag> firstTags = firstTDM.keySet();
        //        for (Map.Entry<Tag, TaxaDistribution> secondEntry : secondTDM.entrySet()) {                    
        //            Tag secondTag = secondEntry.getKey();
        //            if (firstTags.contains(secondTag)) {
        //                System.out.printf("Assert Tag %s depths2 %d greater than depths1 %d\n",
        //                        secondTag.sequence(), secondEntry.getValue().totalDepth(),firstTDM.get(secondTag).totalDepth());
        //                assertTrue(secondEntry.getValue().totalDepth() > firstTDM.get(secondTag).totalDepth());
        //            }                    
        //        }
        System.out.println("Finished GBSSeqToTagDBPluginAppendTest.");
    }
    @Test
    public void testRemoveSecondCutSiteIndexOf() {
        
        GBSSeqToTagDBPlugin GBSSeqToTagdb = new GBSSeqToTagDBPlugin();
        Class GBSSeqClass = GBSSeqToTagdb.getClass();
        Method removeSecondCutSiteIndexOf;
        try {
            Class[] args = new Class[]{ String.class,Integer.TYPE};
            removeSecondCutSiteIndexOf = GBSSeqClass.getDeclaredMethod("removeSecondCutSiteIndexOf", args);
            removeSecondCutSiteIndexOf.setAccessible(true);
            GBSSeqToTagdb.enzyme("ApeKI");
            GBSSeqToTagdb.likelyReadEndStrings = new String[]{"ATGCAT", "ATGCAAGAT"};
            GBSSeqToTagdb.readEndCutSiteRemnantLength = 5;
            String seq1 = "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATGCAAGATAGGAACAGCGCTAGGGGAATG";
            String seq2 = "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATGCATGATAGGAACAGCGCTAGGGGAATG";
            String expectedTag = "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATGCA";
            // Test string with second cut site = likelyReadEndString[1] 
            // ATGCAAGAT at position 41 in string
            Tag tag = (Tag) removeSecondCutSiteIndexOf.invoke(GBSSeqToTagdb, (Object)seq1,55);
          // System.out.println("\nLCJ - tag sequence from seq1 is: " + tag.sequence());
          // System.out.println("LCJ - expecated tag sequence :   " + "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATGCA");
            assertTrue(expectedTag.equals(tag.sequence()));
            
            // Test string with 2nd cut site = likelyReadEndString[0]
            // ATGCAT at position 41 in string
            tag = (Tag) removeSecondCutSiteIndexOf.invoke(GBSSeqToTagdb, (Object)seq2,55);
           // System.out.println("\nLCJ - tag sequence from seq2 is: " + tag.sequence());
           // System.out.println("LCJ - expecated tag sequence :   " + "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATGCA");
            assertTrue(expectedTag.equals(tag.sequence()));
            
            // Test sequence with NO second cut site appearing 
            expectedTag = "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATTCATGATAGGAAC";
            String seq3 = "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATTCATGATAGGAACAGCGCTAGGGGAATG";
            tag = (Tag) removeSecondCutSiteIndexOf.invoke(GBSSeqToTagdb, (Object)seq3,55);
            //tag = removeSecondCutSiteIndexOf(seq3,55);
            assertTrue(expectedTag.equals(tag.sequence()));
        } catch (NoSuchMethodException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (SecurityException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (InvocationTargetException ite) {
            System.out.println("LCJ - ite cause:");
            ite.getCause();
            System.out.println("LCJ - ite stackback");
            ite.printStackTrace();
        } catch (Exception exc)  {
            exc.printStackTrace();;
            
        }
    }

//    private Tag removeSecondCutSiteIndexOf(String seq, int preferredLength ) {
//        Tag tag = null;
//        String[] likelyReadEndStrings = new String[]{"ATGCAT", "ATGCAAGAT"};
//        int readEndCutSiteRemnantLength = 5;
//
//        // handle overlapping cutsite for ApeKI enzyme
////        if (enzyme().equalsIgnoreCase("ApeKI")) {
////            if (seq.startsWith("CAGCTGC") || seq.startsWith("CTGCAGC")) {
////                seq = seq.substring(3,seq.length());
////            }             
////        }
//        int indexOfReadEnd = -1;
//        String shortSeq = seq.substring(20);
//        System.out.println("LCJ - shortSeq is:" + shortSeq);
//        for (String readEnd: likelyReadEndStrings){
//            int indx = shortSeq.indexOf(readEnd);
//            System.out.println("LCJ - checking for readEnd:" + readEnd + " indx is:" + indx);
//            if (indx > 0 ) {
//                if (indexOfReadEnd < 0 || indx < indexOfReadEnd) {
//                    System.out.println("LCJ - indexOfReadEnd now set to: " + indx);
//                    indexOfReadEnd = indx;
//                } 
//            }
//        }
//
//        System.out.println("LCJ - indexOfReadEnd: " + indexOfReadEnd + " preferredLength: " + preferredLength);
//        if (indexOfReadEnd > 0 &&
//                (indexOfReadEnd + 20 + readEndCutSiteRemnantLength < preferredLength)) {
//            // trim tag to sequence up to & including the cut site
//            tag = TagBuilder.instance(seq.substring(0, indexOfReadEnd + 20 + readEndCutSiteRemnantLength)).build();
//           // System.out.println("\nLCJ - tag sequence from IF is: " + tag.sequence());
//           // System.out.println("LCJ - expecated tag sequence : " + "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATGCA");
// 
//        } else {
//            int seqEnd = (byte) Math.min(seq.length(), preferredLength);
//            tag = TagBuilder.instance(seq.substring(0,seqEnd)).build();
//           // System.out.println("\nLCJ - tag sequence from ELSE is: " + tag.sequence());
//            //System.out.println("LCJ - expecated tag sequence   : " + "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATGCA");
//        }
//        return tag;
//    }
    // full pipeline and biology evaluation helper methods moved
    // to EvaluateSNPCallQualityOfPipelineTest.java
}