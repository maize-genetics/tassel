/*
 * FastqToTagCountPluginTest
 */
package net.maizegenetics.analysis.gbs.repgen;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import org.junit.Test;
import org.junit.Ignore;
import java.lang.reflect.Method;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;

//import jdk.nashorn.internal.ir.annotations.Ignore;
import net.maizegenetics.analysis.gbs.repgen.RepGenLoadSeqToDBPlugin.TagCountQualityScoreMap;
import net.maizegenetics.analysis.gbs.v2.DiscoverySNPCallerPluginV2;
import net.maizegenetics.analysis.gbs.v2.GBSUtils;
import net.maizegenetics.analysis.gbs.v2.SAMToGBSdbPlugin;
import net.maizegenetics.analysis.gbs.v2.TagExportToFastqPlugin;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;

/**
 * Testing of the pipeline using repeat GBS data
 * 
 * Changes made to use the RepGenSchema.  All calls to gbvs2 plugins
 * using db created in RepGenLoadSeqToDBPlugin will fail due to
 * calling TagSQLite vs RepGenSQLite.
 * 
 * @author Ed Buckler
 */
public class RepGenLoadSeqToDBPluginTest {
//    public static final String nextGenDir="/Users/edbuckler/SolexaAnal/";
//    public static final String repGenInputDir="/Users/edbuckler/SolexaAnal/rGBS/Maize3/";
//    public static final String repGenInputDir="/Volumes/Nextgen/rGBS/160429_M01032_0387_000000000-ANP68/";
    public static final String nextGenDir="/Users/lcj34/notes_files/repGen/";
   // public static final String repGenInputDir="/Users/lcj34/notes_files/repGen/repGen_files/";
    public static final String repGenInputDir="/Volumes/Samsung_T1/repgen/single_fastq/"; // was anp68_Locus9_fastq
    private final static String outTagFasta=repGenInputDir + "tagsForAlign.fa.gz";
    private final static String inTagSAM=repGenInputDir + "tagsForAlignauto.sam";

    //public static final String repGenDB_FILE=rGBSInputDir+"rGBSx.db";
    //public static final String repGenDB_FILE="/Users/edbuckler/SolexaAnal/repGen/rGBSx.db";
    public static final String repGenDB_FILE=nextGenDir + "repGen.db"; // remember - can't write to Samsung drive - SQL barfs!
    //public static final String repGen_TESTING_KEY_FILE=repGenInputDir+"repGen.key.txt";
    public static final String repGen_TESTING_KEY_FILE=nextGenDir +"repGen_keyFiles/anp68_locus9_key.txt";
    public static final String rGBS_VCF=repGenInputDir+"repGenOutput2.vcf";


    public RepGenLoadSeqToDBPluginTest() {
    }

    
    // This test written to create a key file when needed.
    // WIll there always be the same number of underscores in the file names ??
    // Can I always pull the FullSampleName from the 5th  tokenized place in the file?
    @Test
    public void testCreateKeyFileFromFastq() {
        LoggingUtils.setupDebugLogging();
        System.out.println("Begin testCreateKeyFileFromFastq");
        
        //String inputDir="/Volumes/Samsung_T1/repgen/160429_M01032_0387_000000000-ANP68/";
        String inputDir = "/Volumes/Samsung_T1/repgen/170111_M01032_0465_000000000-AYH7H/";
        //String outputFile = "/Users/lcj34/notes_files/repGen/repGen_keyfiles/anp68_locus9_key.txt";
        //String outputFile = "/Users/lcj34/notes_files/repGen/repGen_keyfiles/all_locus_key.txt"; // MAY BE WRONG - DE_2 sample !!
        String outputFile = "/Users/lcj34/notes_files/repGen/repGen_keyfiles/ije_key.txt";
        
        List<Path> directoryFiles= DirectoryCrawler.listPaths(GBSUtils.inputFileGlob, Paths.get(inputDir).toAbsolutePath());
        if(directoryFiles.isEmpty()) {
            System.out.println("No files matching:"+ GBSUtils.inputFileGlob);
            return;
        }
        BufferedWriter bw=Utils.getBufferedWriter(outputFile);
        try {
            bw.write("FileName\tFullSampleName\tGID\n");
            int GID = 3000; // picked at random
      
            for (Path inputSeqFile : directoryFiles) {
                String[] filenameField = inputSeqFile.getFileName().toString().split("_");
                // ***NOTE 0 change the fullSampleName index to match your input file !!
                String fullSampleName = filenameField[9]; // will it always be in this spot?
                String fileLine = inputSeqFile.getFileName().toString() + "\t" + fullSampleName 
                        + "\t" + Integer.toString(GID) + "\n";
                bw.write(fileLine);
                GID++;
            } 
            bw.close();
        } catch (IOException ioe) {
            System.out.println("LCJ - error processing files in testCreateKeyFileFromFastq");
            ioe.printStackTrace();
        } 
        System.out.println("Key File written to " + outputFile);
    }

    @Test
    public void testcalculateTagAveQS() throws Exception {
        String tagSeq1 = "GCACAAGTTGTCCTGC";
        String tagSeq2 = "ATCGAGCCTGGCCTGC";
        String qulSeq1 = "ABD;->DFABF-../9";
        String qulSeq2 = "GB-;->DFABF-../G";
        Tag tag1 = TagBuilder.instance(tagSeq1).build();
        Tag tag2 = TagBuilder.instance(tagSeq2).build();
        
        TagCountQualityScoreMap tagCntQSMap = new TagCountQualityScoreMap(50000000,0.95f, 128, 4);
        List<String> tagScores1 =  Collections.synchronizedList(new ArrayList<String>());
        tagScores1.add(qulSeq1); // add same one 3 times
        tagScores1.add(qulSeq1);
        tagScores1.add(qulSeq1);
        tagCntQSMap.put(tag1, tagScores1);
        
        List<String> tagScores2 =  Collections.synchronizedList(new ArrayList<String>());
        tagScores2.add(qulSeq2); // add same one 3 times
        tagScores2.add(qulSeq2);
        tagScores2.add(qulSeq2);
        tagCntQSMap.put(tag2, tagScores2);

        // First test:  Create a couple tags.
        // Create a score string for each of the tags.
        // Add the score string to each tag multiple times
        // test that the calculated average score value equals the original score value.
        // This shows we added/divided correctly and got original result back
                 
        
        Map<Tag,Tuple<Integer,String>> tagInfoMap= callCalculateTagAveQS(tagCntQSMap);
        // TODO then fix the print statements below to test values we expect
//        tagInfoMap.forEach((t,s) -> {
//            String tagSeq = t.sequence();
//            String tagScore = s.y;
//            int numInstances = s.x;
//            assertEquals(3,numInstances);
//            if (tagSeq.equals(tagSeq1)) {
//                assertTrue(qulSeq1.equals(tagScore));
//                System.out.println("testcalculateTagAveQS: tag1 score matches as expected");
//            } else if (tagSeq.equals(tagSeq2)) {
//                assertTrue(qulSeq2.equals(tagScore));
//                System.out.println("testcalculateTagAveQS: tag2 score matches as expected");
//            } else {
//                System.out.println("testcalculateTagAveQS: Tag from calculateTagAveQS does not match either tag we passed in !!");
//            }
//            
//        });

    }
    
    private Map<Tag,Tuple<Integer,String>> callCalculateTagAveQS(TagCountQualityScoreMap tagCntQSMap) throws Exception {
        RepGenLoadSeqToDBPlugin repGenLoadSeqToDBPlugin = new RepGenLoadSeqToDBPlugin();
        Class theClass = repGenLoadSeqToDBPlugin.getClass();
        
        // lcj - this is debug
        Method[] methods = theClass.getDeclaredMethods();
        for(int i = 0; i < methods.length; i++) {
           System.out.println("method = " + methods[i].toString());
        }
        // lcj end debug
        
        // Have not been able to get this call to work.  Still working on it.
        Class tagCntMap = Class.forName("net.maizegenetics.analysis.gbs.repgen.RepGenLoadSeqToDBPlugin$TagCountQualityScoreMap");
       // Constructor constructor = tagCntMap.getDeclaredConstructor(RepGenLoadSeqToDBPlugin.class);
        return null;

//        Object obj = constructor.newInstance(repGenLoadSeqToDBPlugin);
//          // this call isn't working ....
//        Method theMethod = theClass.getDeclaredMethod("calculateTagAveQS", new Class[] {obj.class, int.class} );
//        // need to call the method with proper arguments.
//        theMethod.setAccessible(true);
//
//        int qualityScoreBase = 33;
//        // Hard coding true - need to fix this
//        return (Map<Tag,Tuple<Integer,String>>)theMethod.invoke(repGenLoadSeqToDBPlugin,(Object)tagCntQSMap,qualityScoreBase);
    }

    
    /**
     * Test of performFunction method, of class RepGenLoadSeqToDBPlugin.
     */
    @Test
    public void testRepGenLoadSeqToDBPlugin() {
        LoggingUtils.setupDebugLogging();
        System.out.println(Paths.get(repGenDB_FILE).toAbsolutePath().toString());
        try{
            Files.deleteIfExists(Paths.get(repGenDB_FILE));
        } catch (IOException e) {
            e.printStackTrace();
        }
        long time=System.nanoTime();
        System.out.println();

        new RepGenLoadSeqToDBPlugin()
                .inputDirectory(repGenInputDir)
                .outputDatabaseFile(repGenDB_FILE)
                .keyFile(repGen_TESTING_KEY_FILE)
                .kmerLength(150)
                .minKmerCount(10)
                //.minimumQualityScore(20)  //todo needs to be fixed
                //.deleteOldData(true)
                //.maximumMapMemoryInMb(5500)
                .performFunction(null);
        System.out.printf("TotalTime %g sec%n", (double) (System.nanoTime() - time) / 1e9);
        //String outName= GeneralConstants.TEMP_DIR+"TaxaTest.db";
        //TagsByTaxaHDF5Builder.create(, (Map<Tag,TaxaDistribution>)ds.getData(0).getData(),null);
    }

    // NOTE:  These will fail until all plugins use the RepGenSchema.
    // Everyone plugin that is gbsv2 will fail.
    @Test
    public void testTagExportPlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        System.out.println("Running testTagExportPlugin");
        new TagExportToFastqPlugin()
                .inputDB(repGenDB_FILE)
                .outputFile(outTagFasta)
                .performFunction(null);

    }

    // THis method no longer needed for basic runs.  I have moved the
    // file created from bowtie into the GBS expected dir under
    //  bowtie2.1/tagsForAlign910auto.sam
    // This allows users without bowtie to run these tests
    @Ignore
    @Test
    public void runBowtie() throws Exception{
        String cmd=nextGenDir+"AGPv3/bowtie2-2.2.3/bowtie2 " +
                "-M 4 -x "+nextGenDir+"AGPv3/ZmAGP3_910 " +
                "-U " + Paths.get(outTagFasta).toAbsolutePath().toString() + " "+
                "-S " + Paths.get(inTagSAM).toAbsolutePath().toString();
        //System.out.println(cmd);
        Process start = Runtime.getRuntime().exec(cmd);
        BufferedReader r = new BufferedReader(
                new InputStreamReader(start.getErrorStream()));
        String line = null;
        while ((line = r.readLine()) != null)
        {
            System.out.println(line);
        }
    }

    @Test
    public void runBWA() throws Exception{
        String cmd=nextGenDir+"bwa/bwa aln " +
        nextGenDir+"AGPv3/wholegenome/Zea_mays.AGPv3.22.dna_sm.genome.fa.gz" +
                " " +Paths.get(outTagFasta).toAbsolutePath().toString()+
                " > "+Paths.get(inTagSAM).toAbsolutePath().toString();
        System.out.println(cmd);
        Process start = Runtime.getRuntime().exec(cmd);
        BufferedReader r = new BufferedReader(
                new InputStreamReader(start.getErrorStream()));
        String line = null;
        while ((line = r.readLine()) != null)
        {
            System.out.println(line);
        }
        //"/Users/edbuckler/SolexaAnal/bwa/bwa samse /Users/edbuckler/SolexaAnal/AGPv3/wholegenome/Zea_mays.AGPv3.22.dna_sm.genome.fa.gz /Users/edbuckler/SolexaAnal/rGBS/Maize3/tagsForAlignauto.sai /Users/edbuckler/SolexaAnal/rGBS/Maize3/tagsForAlign.fa.gz > /Users/edbuckler/SolexaAnal/rGBS/Maize3/tagsForAlignauto.sam"
    }

    @Test
    public void testSAMImportPlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        System.out.println("Running testTagExportPlugin");
        new SAMToGBSdbPlugin()
                .gBSDBFile(repGenDB_FILE)
                //.sAMInputFile(GBSConstants.GBS_EXPECTED_BOWTIE_SAM_FILE)
                .sAMInputFile(inTagSAM)
                .performFunction(null);
    }

    @Test
    public void runCaller() throws Exception{

        new DiscoverySNPCallerPluginV2()
                .inputDB(repGenDB_FILE)
                .minLocusCoverage(0.99)
                .maxTagsPerCutSite(3)
                .minMinorAlleleFreq(0.15)

                .deleteOldData(true)
                .startChromosome(new Chromosome("1"))
                .endChromosome(new Chromosome("10"))
                //.gapAlignmentThreshold(0.12) //  default is 1.0 (keep all alignments)
                .performFunction(null);

//        // Verify number of SNPs called that should not have been called, ie conserved sites
//        double typeIError=evaluateConservedSites();
//        // Find percentage of HM3 identified SNPs the pipeline also identified.
//        double power= evaluateHM3Sites();
//        System.out.printf("Power:%g TypeIError:%g RatioP/E:%g %n",power,typeIError,power/typeIError);
    }

    @Test
    public void testProcessData() throws Exception {
        //  This test assumes the database has already been populated.
        // The test will fail if db does not already exist.  Run the
        // DiscoveryPipeline test to create the database
        // (EvaluateSNPCallQualityOfPipelineTest.java)
        LoggingUtils.setupDebugLogging();
        System.out.println(Paths.get(repGenDB_FILE).toAbsolutePath().toString());

        // Delete output file if it exists.
        try{
            Files.deleteIfExists(Paths.get(rGBS_VCF));
        } catch (IOException e) {
            e.printStackTrace();
        }

        long time=System.nanoTime();
        System.out.println();

        new RGBSProductionSNPCallerPlugin()
                .inputDirectory(repGenInputDir)
                .inputGBSDatabase(repGenDB_FILE)
                .keyFile(repGen_TESTING_KEY_FILE)
                .outputGenotypesFile(rGBS_VCF)
                .kmerLength(150)
                .minimumQualityScore(0) // tag read quality score
                //.positionQualityScore(5.0) // snp position quality score
                //.maximumMapMemoryInMb(5500)
                .performFunction(null);
        System.out.printf("TotalTime for vcf: %g sec%n", (double) (System.nanoTime() - time) / 1e9);

        time=System.nanoTime();
        System.out.println();

    }
    
    // This is initial testing of method to create averaged quality score.
    // Need to get the reflection method working to call RepGenLoadSeqToDB:calculateTagAveQS to work.
    // Leaving this here to show there was some testing of calculateTagAveQS's algorithm 
    @Test
    public void testCreateQualityScore () {
        String qualSFFShort1 = "ABD;->DFABF-../9";
        String qualSFFShort2 = "BBD;->DEABF-../9";
        String qualSFFShort3 = "GB-;->DFABF-../G";
        int[] qualShortNum = new int[qualSFFShort1.length()];
        
 
        // Default qualiytScoreBase. In source, this is determined by fastQ files.
        int qualityScoreBase = 33;
        
        for (int idx = 0; idx < qualSFFShort1.length(); idx++) {
            qualShortNum[idx] = (qualSFFShort1.charAt(idx) - qualityScoreBase);
        }
        // print values and create qualityScoreSTring at same time
        StringBuilder sb = new StringBuilder();
        System.out.print("LCJ - qualShortNum values are: ");
        for (int idx = 0; idx < qualShortNum.length; idx++){
            char ch  = (char)(qualShortNum[idx]+ qualityScoreBase) ;
            System.out.print(" " + qualShortNum[idx]);
            sb.append(ch);
        }
        System.out.println();
        // translate back to char string
        System.out.println("LCJ - the recreated quality string is: " + sb.toString());
        System.out.println("LCJ - original qualSFFSHort is       : " + qualSFFShort1);
        
        // test adding quality scores, then dividing by number of scores,
        // and print it back out.
        // First try adding qualSFFShort1 3 times, then divide by 3,
        // show the number.  This should come out to be the same as the original.
        
        // qualShortNum array currently contains values from just 1 interation
        // add the others
        System.out.print("LCJ - tripled qualShortNum values are: ");
        for (int idx = 0; idx < qualSFFShort1.length(); idx++) {
            qualShortNum[idx] += (qualSFFShort1.charAt(idx) - qualityScoreBase);
            qualShortNum[idx] += (qualSFFShort1.charAt(idx) - qualityScoreBase);  
            System.out.print(" " + qualShortNum[idx]);
        }
        System.out.println();
        // divide each by 3, print the average.  In this case, the values will
        // divide evenly by 3 because all we did was triple them.  
        System.out.print("LCJ - averaged qualShortNum values are: ");
        sb.setLength(0);
        for (int idx = 0; idx < qualSFFShort1.length(); idx++) {
            qualShortNum[idx] = qualShortNum[idx]/3;
            char ch  = (char)(qualShortNum[idx]+ qualityScoreBase) ;
            sb.append(ch);
            System.out.print(" " + qualShortNum[idx]);
        }
        System.out.println("\n");
        // translate back to char string
        System.out.println("LCJ - the quality string after averaging: " + sb.toString());
        System.out.println("LCJ - original qualSFFSHort is          : " + qualSFFShort1);
        String sbToString = sb.toString();
        assertTrue(sbToString.equals(qualSFFShort1));
        
        // now add in different numbers:
        System.out.print("LCJ - added all 3 strings, qualShortNum values are: ");
        for (int idx = 0; idx < qualSFFShort1.length(); idx++) {
            qualShortNum[idx] += (qualSFFShort2.charAt(idx) - qualityScoreBase);
            qualShortNum[idx] += (qualSFFShort3.charAt(idx) - qualityScoreBase);  
            System.out.print(" " + qualShortNum[idx]);
        }
        System.out.println();
        // divide each by 3, print the average.  In this case, the values will
        // divide evenly by 3 because all we did was triple them.  IN real instance,
        // this won't always be true as you could add values of 1+1+3 and 5/3 is not even
        // I guess we just take the integer value, which is rounded one way.  SHOuld be fine.
        
        // also recreate quality score string
        System.out.print("LCJ - 2nd averaged qualShortNum values are: ");
        sb.setLength(0);
        for (int idx = 0; idx < qualSFFShort1.length(); idx++) {
            qualShortNum[idx] = qualShortNum[idx]/3;
            char ch  = (char)(qualShortNum[idx]+ qualityScoreBase) ;
            sb.append(ch);
            System.out.print(" " + qualShortNum[idx]);
        }
        System.out.println("\n");
        // translate back to char string
        System.out.println("LCJ - the quality string after 2nd averaging: " + sb.toString());
        System.out.println("LCJ - original qualSFFSHort is              : " + qualSFFShort1);             
    }
    
    @Test
    public void checkTagSequences() {
        System.out.println("LCJ - begin checkTagSequences ...");
        String mydb = "/Users/lcj34/notes_files/repgen/repGen_allLocus_justLoading.db";
        RepGenDataWriter repGenData=new RepGenSQLite(mydb);
        Set<Tag> tagsToAlign = repGenData.getTags();
        System.out.println("LCJ - number of tags in db: " + tagsToAlign.size());
        String outputFile = "/Users/lcj34/notes_files/repgen/ed_debug_out/allLocusTags_c100mintaxa5batch100.txt";
        
        
        // Write the tags to a file
        BufferedWriter bw = Utils.getBufferedWriter(outputFile);
        StringBuffer sb = new StringBuffer();
        int count = 0;
        try {
            for (Tag tag : tagsToAlign) {
                count++;
                sb.append(tag.sequence());
                sb.append("\n");
                if (count > 10000) {
                    bw.write(sb.toString());
                    count = 0;
                    sb.setLength(0);
                    System.out.println("LCJ - wrote next 10000 tags");
                }               
            }
            if (count > 0) {
                bw.write(sb.toString());
            }
            bw.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
        System.out.println("Test finished - tags written to: " + outputFile);
    }

}
