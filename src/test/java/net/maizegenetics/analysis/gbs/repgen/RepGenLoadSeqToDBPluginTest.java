/*
 * RepGenLoadSeqToDBPluginTest
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
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.Ignore;
import org.junit.Test;

import net.maizegenetics.analysis.gbs.repgen.RepGenLoadSeqToDBPlugin.TagCountQualityScoreMap;
import net.maizegenetics.analysis.gbs.v2.DiscoverySNPCallerPluginV2;
import net.maizegenetics.analysis.gbs.v2.GBSUtils;
import net.maizegenetics.analysis.gbs.v2.SAMToGBSdbPlugin;
import net.maizegenetics.analysis.gbs.v2.TagExportToFastqPlugin;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;

/**
 * Loads the {@link RepGenSimData} reads into a RepGen database and checks the tags, per-taxon depths,
 * average quality strings and the reads that must be filtered out.
 * <p>
 * The ignored tests below are the original developer runs on the ANP68 rAmpSeq data. They are kept for
 * reference until #43 decides whether the RepGen pipeline stays. Calls to GBSv2 plugins on a database made by
 * RepGenLoadSeqToDBPlugin fail, because those plugins open it as a TagDataSQLite rather than a RepGenSQLite.
 * 
 * @author Ed Buckler
 */
public class RepGenLoadSeqToDBPluginTest {

    private static final int MIN_KMER_COUNT = 5;

    private static final String nextGenDir = GeneralConstants.DATA_DIR + "RepGen/";
    private static final String repGenInputDir = nextGenDir + "single_fastq/";
    private static final String repGenOutputDir = GeneralConstants.TEMP_DIR + "RepGen/";
    private static final String outTagFasta = repGenOutputDir + "tagsForAlign.fa.gz";
    private static final String inTagSAM = repGenOutputDir + "tagsForAlignauto.sam";
    private static final String repGenDB_FILE = repGenOutputDir + "repGen.db";
    private static final String repGen_TESTING_KEY_FILE = nextGenDir + "repGen_keyFiles/anp68_locus9_key.txt";
    private static final String rGBS_VCF = repGenOutputDir + "repGenOutput2.vcf";

    @Test
    public void testcalculateTagAveQS() {
        String qulSeq1 = "ABD;->DFABF-../9";
        String qulSeq2 = "BBD;->DEABF-../9";
        String qulSeq3 = "GB-;->DFABF-../G";
        Tag tag1 = TagBuilder.instance("GCACAAGTTGTCCTGC").build();
        Tag tag2 = TagBuilder.instance("ATCGAGCCTGGCCTGC").build();

        TagCountQualityScoreMap tagCntQSMap = new TagCountQualityScoreMap(100, 0.95f, 128, 1);
        tagCntQSMap.put(tag1, Collections.synchronizedList(new ArrayList<>(List.of(qulSeq1, qulSeq1, qulSeq1))));
        tagCntQSMap.put(tag2, Collections.synchronizedList(new ArrayList<>(List.of(qulSeq1, qulSeq2, qulSeq3))));

        Map<Tag, Tuple<Integer, String>> tagInfoMap = RepGenLoadSeqToDBPlugin.calculateTagAveQS(tagCntQSMap, 33);

        assertEquals(2, tagInfoMap.size());
        assertEquals(3, tagInfoMap.get(tag1).x.intValue());
        assertEquals("Averaging three identical quality strings returns the original", qulSeq1, tagInfoMap.get(tag1).y);
        assertEquals(3, tagInfoMap.get(tag2).x.intValue());
        // Per-base mean of the Phred scores, truncated: e.g. the first base is (32 + 33 + 38) / 3 = 34, 'C'.
        assertEquals("CB<;->DEABF-../=", tagInfoMap.get(tag2).y);
    }

    /**
     * Test of performFunction method, of class RepGenLoadSeqToDBPlugin.
     */
    @Test
    public void testRepGenLoadSeqToDBPlugin() throws Exception {
        RepGenSimData data = RepGenSimData.createUnder("load");
        String db = loadDB(data, MIN_KMER_COUNT);

        Map<String, TaxaDistribution> tagTaxa;
        TaxaList taxaList;
        RepGenSQLite repGenData = new RepGenSQLite(db);
        try {
            taxaList = repGenData.getTaxaList();
            tagTaxa = repGenData.getAllTagsTaxaMap().entrySet().stream()
                    .collect(Collectors.toMap(e -> e.getKey().sequence(), Map.Entry::getValue));
        } finally {
            repGenData.close();
        }

        assertEquals(RepGenSimData.TAXA.length, taxaList.numberOfTaxa());
        assertEquals("Two alleles per amplicon; errors, private alleles, low counts, short reads and N reads are dropped",
                data.expectedTags.keySet(), tagTaxa.keySet());
        for (String discarded : data.discardedSequences) {
            assertFalse(tagTaxa.containsKey(discarded));
        }
        assertFalse("Total depth is below minKmerCount", tagTaxa.containsKey(data.lowCountTag.sequence));

        for (Map.Entry<String, RepGenSimData.ExpectedTag> entry : data.expectedTags.entrySet()) {
            String sequence = entry.getKey();
            assertEquals("Tags are truncated to the kmer length", RepGenSimData.KMER_LENGTH, sequence.length());
            int[] depths = tagTaxa.get(sequence).depths();
            for (int taxon = 0; taxon < RepGenSimData.TAXA.length; taxon++) {
                int index = taxaList.indexOf(RepGenSimData.TAXA[taxon]);
                assertEquals("Depth of " + RepGenSimData.TAXA[taxon], entry.getValue().depths[taxon], depths[index]);
            }
        }

        Map<String, Tuple<Integer, String>> stored = readTagQualities(db);
        assertEquals(data.expectedTags.size(), stored.size());
        for (Map.Entry<String, RepGenSimData.ExpectedTag> entry : data.expectedTags.entrySet()) {
            Tuple<Integer, String> instancesAndQuality = stored.get(entry.getKey());
            assertEquals(entry.getValue().totalDepth(), instancesAndQuality.x.intValue());
            assertEquals(entry.getValue().averageQuality(), instancesAndQuality.y);
        }
    }

    @Test
    public void testMinKmerCount() throws Exception {
        RepGenSimData data = RepGenSimData.createUnder("loadMinKmerCount");

        RepGenSQLite repGenData = new RepGenSQLite(loadDB(data, data.lowCountTag.totalDepth()));
        try {
            TaxaList taxaList = repGenData.getTaxaList();
            Map<Tag, TaxaDistribution> tagTaxa = repGenData.getAllTagsTaxaMap();
            assertEquals(data.expectedTags.size() + 1, tagTaxa.size());
            int[] depths = tagTaxa.get(TagBuilder.instance(data.lowCountTag.sequence).build()).depths();
            for (int taxon = 0; taxon < RepGenSimData.TAXA.length; taxon++) {
                assertEquals(data.lowCountTag.depths[taxon], depths[taxaList.indexOf(RepGenSimData.TAXA[taxon])]);
            }
        } finally {
            repGenData.close();
        }
    }

    static String loadDB(RepGenSimData data, int minKmerCount) {
        String db = data.dir.resolve("repGen.db").toString();
        new RepGenLoadSeqToDBPlugin()
                .inputDirectory(data.fastqDir.toString())
                .outputDatabaseFile(db)
                .keyFile(data.keyFile.toString())
                .kmerLength(RepGenSimData.KMER_LENGTH)
                .minKmerCount(minKmerCount)
                .maximumKmerNumber(10_000)
                .performFunction(null);
        return db;
    }

    /** Tag sequence to (number of instances, average quality string), read straight from the tag table. */
    private static Map<String, Tuple<Integer, String>> readTagQualities(String db) throws Exception {
        Map<String, Tuple<Integer, String>> result = new HashMap<>();
        try (Connection connection = DriverManager.getConnection("jdbc:sqlite:" + db);
             ResultSet rs = connection.createStatement()
                     .executeQuery("select sequence, seqlen, qualityScore, numTagInstances from tag where isReference = 0")) {
            while (rs.next()) {
                String sequence = TagBuilder.instance(rs.getBytes("sequence"), rs.getShort("seqlen")).build().sequence();
                result.put(sequence, new Tuple<>(rs.getInt("numTagInstances"), rs.getString("qualityScore")));
            }
        }
        return result;
    }

    // This test written to create a key file when needed.
    // WIll there always be the same number of underscores in the file names ??
    // Can I always pull the FullSampleName from the 5th  tokenized place in the file?
    @Ignore("#43: writes a key file for the AYH7H run's FASTQ files, which are not in tassel_test_data; "
            + "a utility that asserts nothing. RepGenSimData writes the key file the tests above use")
    @Test
    public void testCreateKeyFileFromFastq() {
        LoggingUtils.setupDebugLogging();
        System.out.println("Begin testCreateKeyFileFromFastq");
        
        String inputDir = nextGenDir + "170111_M01032_0465_000000000-AYH7H/";
        String outputFile = repGenOutputDir + "ije_key.txt";
        
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

    // NOTE:  These will fail until all plugins use the RepGenSchema.
    // Everyone plugin that is gbsv2 will fail.
    @Ignore("#43: TagExportToFastqPlugin is a GBSv2 plugin and can't read a RepGen database; this also needs "
            + "a database loaded from the ANP68 reads, which are not in tassel_test_data")
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
    @Ignore("#43: runs a bowtie2 binary from a developer's machine; its output is in the archive as bowtie2.1/tagsForAlign910auto.sam")
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

    @Ignore("#43: runs a bwa binary against the full AGPv3 genome on a developer's machine; neither is in "
            + "tassel_test_data, and it only prints bwa's stderr")
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
    }

    @Ignore("#43: SAMToGBSdbPlugin is a GBSv2 plugin and can't write to a RepGen database; this also needs "
            + "the SAM from runBowtie")
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

    @Ignore("#43: DiscoverySNPCallerPluginV2 is a GBSv2 plugin and can't read a RepGen database; this also needs "
            + "a database loaded from the ANP68 reads, which are not in tassel_test_data")
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

    @Ignore("#43: RGBSProductionSNPCallerPlugin opens the database as a GBSv2 TagDataSQLite and needs SNP "
            + "positions, which a RepGen database never has; this also needs the ANP68 reads, which are not in tassel_test_data")
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
    @Ignore("#43: hand-computes quality-score averaging without calling the plugin; RepGenQualityScoreUnitTest "
            + "and testcalculateTagAveQS above cover RepGenLoadSeqToDBPlugin.calculateTagAveQS")
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
    
    @Ignore("#43: needs a developer's RepGen database of all ANP68 loci, which is not in tassel_test_data; "
            + "writes the tag sequences to a file and asserts nothing")
    @Test
    public void checkTagSequences() {
        System.out.println("LCJ - begin checkTagSequences ...");
        String mydb = nextGenDir + "repGen_allLocus_justLoading.db";
        RepGenDataWriter repGenData=new RepGenSQLite(mydb);
        Set<Tag> tagsToAlign = repGenData.getTags();
        System.out.println("LCJ - number of tags in db: " + tagsToAlign.size());
        String outputFile = repGenOutputDir + "allLocusTags_c100mintaxa5batch100.txt";
        
        
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
