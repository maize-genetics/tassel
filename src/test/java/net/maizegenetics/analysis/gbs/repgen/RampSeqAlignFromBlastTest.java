/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.Test;

import com.google.common.collect.Multimap;

import net.maizegenetics.analysis.gbs.v2.GBSUtils;
import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.LoggingUtils;

/**
 * @author lcj34
 *
 */
public class RampSeqAlignFromBlastTest {
    public static final String nextGenDir="/Users/lcj34/notes_files/repGen/";

     public static final String primerFile = nextGenDir + "primers.txt";
     public static final String referenceGenome = "/Volumes/Samsung_T1/machineLearningFiles/refGenomeFiles/Zea_mays.AGPv4.dna.genome.fa";

     public static final String repGenDB_FILE=nextGenDir + "repGen.db";
     public static final String forwardp = "GCACAAGTTGTCCTGCTTCC";
     public static final String reversep = "ACGTCGATCTAGGGGGTCTC";
     public static final String blastFile = "/Users/lcj34/notes_files/repgen/LDmetric_graphs/junit/blast_98per_3cols.txt";
     

     @Test
     public void testGetPEARAssemblies () {
         // This test executes the PEAR (Paired-End reAd mergeR).
         // It is assumed PEAR lives on the user's system.
         // user should adjust the forwardDir, reverseDir and outputDir's to
         // match their system.
          
         String forwardDir = "/Users/lcj34/notes_files/repgen/anp68_ije_R1_fastq";
         String reverseDir = "/Users/lcj34/notes_files/repgen/anp68_ije_R2_fastq";
         String outputDir = "/Users/lcj34/notes_files/repgen/pear_output/";
         
         // Get list of 5' (forward strand) files
         List<Path> forwardFiles= DirectoryCrawler.listPaths(GBSUtils.inputFileGlob, Paths.get(forwardDir).toAbsolutePath());
         if(forwardFiles.isEmpty()) {
             System.out.println("No forward files matching:"+ GBSUtils.inputFileGlob);
             return ;
         }
         
         if (forwardFiles.size() == 0) {
             System.out.println("testGetPEARAssemblies - found NO forward files to process.");
             System.out.println("Please verify your directory contains fastq files.");
             return ; // no files in this directory to process
         }
         Collections.sort(forwardFiles);
         System.out.println("Found " + forwardFiles.size() + " forward files to process");
         
         
         // Get list of 3' (reverse strand) fastq files
         List<Path> reverseFiles= DirectoryCrawler.listPaths(GBSUtils.inputFileGlob, Paths.get(reverseDir).toAbsolutePath());
         if(reverseFiles.isEmpty()) {
             System.out.println("No reverse files matching:"+ GBSUtils.inputFileGlob);
             return ;
         }
         
         if (reverseFiles.size() == 0) {
             System.out.println("testGetPEARAssemblies - found NO reverse files to process.");
             System.out.println("Please verify your reverse directory contains fastq files.");
             return ; // no files in this directory to process
         }
         System.out.println("Found " + reverseFiles.size() + "  reverse files to process");        
         Collections.sort(reverseFiles);
         
         // Need to be matching forward/reverse files. The assmption is the name is identical
         // on each except for R1 vs R2.  Sorting the list should result in the lists matching
         // such that traversing the lists together gives complementary files to send to PEAR.
         if (forwardFiles.size() != reverseFiles.size()) {
             System.out.println("ERROR - number of forward files " + forwardFiles.size() 
             + " does not match number of reverse files " + reverseFiles.size());          
             return;
         }
         
         long time=System.nanoTime();
         
         int fileCount = 0;
         int totalFilesProcessed = 0;
         for (int idx = 0; idx < forwardFiles.size(); idx++) {
             fileCount++;
             totalFilesProcessed++;
             Path forwardPath = forwardFiles.get(idx);
             Path reversePath = reverseFiles.get(idx);
             
             // Create name for output files.  the "pear_output" variable is the
             // prefix PEAR adds to each file.  For example:  If the prefix were "combined",
             // PEAR would create files called 
             //     combined.assembled.fastq
             //     combined.discarded.fastq
             //     combined.unassembled.forward.fastq
             //     combined.unassembled.reverse.fastq
             // It is the *.assembled.fastq files that will be loaded into RepGenLoadSeqToDBPlugin
             
             // THis code takes everything up to the last underscore.  Depending on your file,
             // you could take up to the index of the first period, or any other filtering you like.
             // This means the files look like:
             //   6970_1465_34282_ANP68_10357781_CI_7_Locus9_ACTGAGCG_CCTAGAGT.assembled.fastq
             // An alternative would be to do lastindexof ("_")
             String forwardName = forwardPath.getFileName().toString();
             String sampleName = forwardName.substring(0,forwardName.lastIndexOf("_")); // for writing the assembly
             
             String pear_output = outputDir + sampleName;
             
             String cmd="/usr/local/bin/pear  -j 8 -f " + forwardPath + " -r " + reversePath 
                     + " -o " + pear_output ;
             try {
                 Process start = Runtime.getRuntime().exec(cmd);
                 BufferedReader br = new BufferedReader(
                         new InputStreamReader(start.getErrorStream()));
                 String line = null;
                 while ((line = br.readLine()) != null)
                 {
                     System.out.println(line);
                 }
             } catch (Exception exc) {
                 exc.printStackTrace();
             }
             if (fileCount > 19) {
                 System.out.println("PEAR finished processing " + totalFilesProcessed);
                 fileCount = 0;
             }
         }
         System.out.println("Finished PEAR processing " + forwardFiles.size() + " files to directory " + outputDir);
         
         System.out.printf("\nTotalTime for testGetPEARAssemblies: %g sec%n", (double) (System.nanoTime() - time) / 1e9);
         
     }
     
     @Test
     public void testRampSeqAlignFromBlastPlugin() {
         LoggingUtils.setupDebugLogging();
 
         long time=System.nanoTime();
 

         System.out.println("testRampSeqAlignFromBlastPlugin, calling RampSeqAlignFromBlastPlugin  ....");
         // Run the analysis to load the aligner tables
         new RampSeqAlignFromBlastTags()
          .dBFile(repGenDB_FILE)
          .blastFile(blastFile)
          .refGenome(referenceGenome)
          .forwardp(forwardp)
          .reversep(reversep)
         .performFunction(null);
     }
        

     @Test
     public void testRetrievingAlignments() {
         // THis test assumes you have a populated DB.  
         RepGenDataWriter repGenData=new RepGenSQLite(repGenDB_FILE);
         Set<Tag> tagsToAlign = repGenData.getTags();
         List<Tag> tagList = new ArrayList<Tag>(tagsToAlign);
         List<Tag> shortTagList = new ArrayList<Tag>();
         
         System.out.println("LCJ - total number of tags in db: " + tagsToAlign.size());

         int numTags = 0;
         System.out.println("LCJ - first 1 tags from db:");
         while(numTags < 1) { // want to test with a small number
             shortTagList.add(tagList.get(numTags));
             numTags++;
         }
         
         long time=System.nanoTime();
         // Test getTagAlignmentsForTags - just 1 tag
         System.out.println("\ngetting tag alignments for 1 tag");
         Multimap<Tag,AlignmentInfo> oneTagALignments = repGenData.getTagAlignmentsForTags(shortTagList, 0);
         System.out.println("total alignments for tag: " + oneTagALignments.size());
         System.out.printf("TotalTime for getTagAlignmentsForTag: %g sec%n", (double) (System.nanoTime() - time) / 1e9);
         
         System.out.println("\ngetting reftag alignments for 1 tag");
         time=System.nanoTime();
         Multimap<Tag,AlignmentInfo> oneTagRefAlignments = repGenData.getRefAlignmentsForTags(shortTagList,0);
         System.out.printf("TotalTime for getRefAlignmentsForTag: %g sec%n", (double) (System.nanoTime() - time) / 1e9);
         System.out.println("Total ref alignments for tag - should be twice number of ref tags: " + oneTagRefAlignments.size());
         
         Set<RefTagData> refTags = repGenData.getRefTags();
         List<RefTagData> refTagList = new ArrayList<RefTagData>(refTags);
         List<RefTagData> shortRefTagList = new ArrayList<RefTagData>();
         
         int numRefTags = 0;
         System.out.println("LCJ - first 1 tags from db:");
         while(numRefTags < 1) { // want to test with a small number
             shortRefTagList.add(refTagList.get(numRefTags));
             numRefTags++;
         }
         
         Collection<AlignmentInfo> aiCollection = oneTagRefAlignments.get(shortTagList.get(0));
         List<AlignmentInfo> aiList = new ArrayList<AlignmentInfo> (aiCollection);
         AlignmentInfo ai = aiList.get(0);
         RefTagData refTagForTest = new RefTagData(ai.tag2(), ai.tag2chrom(), ai.tag2pos(), ai.ref_genome());
         
         System.out.println("Calling getTagAlignmentsForRefTag");
         time=System.nanoTime();
         Multimap<Tag,AlignmentInfo> tagAlignmentsForRef = repGenData.getTagAlignmentsForRefTag(refTagForTest,0);
         System.out.printf("TotalTime for getTagAlignmentsForRef: %g sec%n", (double) (System.nanoTime() - time) / 1e9);
         System.out.println("Total tag alignments for reftag : " + tagAlignmentsForRef.size());
         
         System.out.println("\nFinished junit !!");
     }
}