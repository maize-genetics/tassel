/**
 * 
 */
package net.maizegenetics.analysis.gbs.v2;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.Random;

import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.map.PositionListIOUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.tag.TagData;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.util.LoggingUtils;

/**
 * This class contains junit tests to verify the biology of the
 * GBSV2 Pipeline.
 * 
 * @author lcj34
 *
 */
public class EvaluateSNPCallQualityOfPipelineTest {
    private final static String outTagFasta=GBSConstants.GBS_TEMP_DIR + "tagsForAlign.fa.gz";
   // private final static String inTagSAM=GBSConstants.GBS_TEMP_DIR + "tagsForAlign910auto.sam"; // used if running bowtie during test
    private final static String inTagSAM=GBSConstants.GBS_EXPECTED_DIR + "bowtie2.1/tagsForAlign910auto.sam";
    private final static String conservedSitesFile= GeneralConstants.DATA_DIR + "CandidateTests/chr910SNPConserves.txt.gz";
    private final static String HM3SitesFile= GeneralConstants.DATA_DIR + "CandidateTests/c910_hmp31_q30_qualityPositions_1-20M.txt.gz";
    private final static Random random=new Random(0);
    public static final String REFERENCE_GENOME_EXPECTED_DIR = "dataFiles/CandidateTests/";
   // public static final String REFERENCE_GENOME_FILE = REFERENCE_GENOME_EXPECTED_DIR + "Zea_maysParsed.AGPv3.23.dna.genome.fa";
    public static final String REFERENCE_GENOME_FILE = REFERENCE_GENOME_EXPECTED_DIR + "ZmB73_RefGen_v3_chr9_10_1st25MB.fasta.gz";
    public static final String REFERENCE_GENOME_FILEv2 = GBSConstants.GBS_DATA_DIR + "ZmB73_RefGen_v2_chr9_10_1st20MB.fasta";


    // temp storage place when running bowtie2 below with "very-sensitive-local" - this matches V1 setting
    public static String samFileAGPv2_local = "dataFiles/CandidateTests/tagsForAlign910_agpv2_local.sam";
    @Test
    public void testBiologyOfDiscoveryTBT() throws Exception {
        //TODO implement this code

        // Moved from GBSSeqToTagDBPluginTest
        // Tag comparison to GBSv1 done with tests and manual
        // comparison via GBSv2BiologyCompareTest.java
        
        LoggingUtils.setupDebugLogging();
        System.out.println("Running testBiologyOfDiscoveryTBT");
        //run for minimum quality of 0
        //minimum number of tags
        //sum of the total depths
        //test for min tag length
        //test for max tag length

        //run for minimum quality of 30
        //minimum number of tags
        //sum of the total depths
        //test for max tag length

        //test for tags from known good SNPs, and total depth and distribution
        //perhaps choose 2 loci - one easy & one with lots of indels

    }

    @Test
    public void pipeLineWithVariableSitesTest() throws Exception {
        LoggingUtils.setupDebugLogging();
        runPipeline(0.01); // for variable sites
    }

    @Test
    public void pipelineIncludingInvariantSites() throws Exception {
        LoggingUtils.setupDebugLogging();
        runPipeline(0.0); // for invariant sites
    }

    public void runPipeline(double MAF) throws Exception{
        System.out.println(Paths.get(GBSConstants.GBS_GBS2DB_FILE).toAbsolutePath().toString());
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
        //.maximumMapMemoryInMb(5500)
        .performFunction(null);

        new TagExportToFastqPlugin()
        .inputDB(GBSConstants.GBS_GBS2DB_FILE)
        .outputFile(outTagFasta)
        .performFunction(null);

        // Bowtie no longer needs to be run for junits using the
        // short chrom 9-10 test files.  The bowtie2 output has been
        // stored in the expected results file defined in variable "inTagSAM"
        // This allows users without bowtie on their system to run this
        // junit pipeline without changes to the code.
        //runBowtie(); 
        new SAMToGBSdbPlugin()
        .gBSDBFile(GBSConstants.GBS_GBS2DB_FILE)
        .sAMInputFile(inTagSAM)
        .minMAPQ(2)
        .deleteOldData(true)
        .performFunction(null);	  

        new DiscoverySNPCallerPluginV2()
        .inputDB(GBSConstants.GBS_GBS2DB_FILE)
        .minMinorAlleleFreq(MAF)
        .startChromosome(new Chromosome("9"))
        .endChromosome(new Chromosome("10"))
        .referenceGenomeFile(REFERENCE_GENOME_FILE)
        .deleteOldData(true)
        //.gapAlignmentThreshold(0.12) //  default is 1.0 (keep all alignments)
        .performFunction(null);

        // Verify number of SNPs called that should not have been called, ie conserved sites
        double typeIError=evaluateConservedSites();
        // Find percentage of HM3 identified SNPs the pipeline also identified.
        double power= evaluateHM3Sites();
        System.out.printf("Power:%g TypeIError:%g RatioP/E:%g %n",power,typeIError,power/typeIError);
        
        // Verify total number of SNPs:
        TagData tagData =new TagDataSQLite(GBSConstants.GBS_GBS2DB_FILE);
        PositionList positsV2 = tagData.getSNPPositions();
        if (MAF == 0.01) {
            // actual number SNPs in tassel 5.2.24 is 9631
            assertTrue(positsV2.size() > 9500);
            assertTrue(positsV2.size() < 9700);
            // check conserved and HM3 Site data falls within expected range
            double pToERatio = power/typeIError;
         // actual value in 5.2.24 was 30.0552, 29.8754 when run without a referenece genome
            assertTrue(pToERatio > 28); 
        }
        if (MAF == 0) {
            // THe values are way high - §38717 SNPs ?
            // this is with a reference genome.  There are more SNPs
            // if we don't use the referenece gnome
            assertTrue(positsV2.size() > 638000);
            assertTrue(positsV2.size() < 640000);
            //  When MAF=0, ratio: (power of finding HM3 sites) to Error is quite poor
            // for the 5.2.24 load,
            //Power:0.0258292 TypeIError:0.159280 RatioP/E:0.162162 
        }       
        ((TagDataSQLite)tagData).close();
    }

    @Test
    public void agpv2WithVariableSitesTest() throws Exception {
        LoggingUtils.setupDebugLogging();
        runPipelineAGPv2(0.01); // for variable sites
    }

    @Test
    public void agpv2IncludingInvariantSites() throws Exception {
        LoggingUtils.setupDebugLogging();
        runPipelineAGPv2(0.0); // for invariant sites
    }
    public void runPipelineAGPv2(double MAF) throws Exception{
        // This comes from net.maizegenetics.analysis.gbs.ProductionSNPCallerPluginTest.java
        GenotypeTable ExpectedGenosV1 = ImportUtils.readGuessFormat(GBSConstants.GBS_EXPECTED_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_OUT_FILE);

        System.out.println(Paths.get(GBSConstants.GBS_GBS2DB_FILE).toAbsolutePath().toString());
        String samFileAGPv2 = "dataFiles/CandidateTests/tagsForAlign910_agpv2.sam";
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
        .minimumQualityScore(20) // test with 0 for V1 compatibility
        //.maximumMapMemoryInMb(5500)
        .performFunction(null);

        new TagExportToFastqPlugin()
        .inputDB(GBSConstants.GBS_GBS2DB_FILE)
        .outputFile(outTagFasta)
        .performFunction(null);
 
        //runBowtie(); // comment out
        new SAMToGBSdbPlugin()
        .gBSDBFile(GBSConstants.GBS_GBS2DB_FILE)
        .sAMInputFile(samFileAGPv2) // using stored .sam for AGPv2 created with -M 4
        //.sAMInputFile(samFileAGPv2_local) // using very-sensitive-local
        .deleteOldData(true)
        .minMAPQ(2)
        .performFunction(null); 
        
        new DiscoverySNPCallerPluginV2()
        .inputDB(GBSConstants.GBS_GBS2DB_FILE)
        .minMinorAlleleFreq(MAF)
        .startChromosome(new Chromosome("9"))
        .endChromosome(new Chromosome("10"))
        .referenceGenomeFile(REFERENCE_GENOME_FILEv2)
        .deleteOldData(true)
        //.gapAlignmentThreshold(0.12) //  default is 1.0 (keep all alignments)
        .performFunction(null);

        // Verify number of SNPs called that should not have been called, ie conserved sites
        double typeIError=evaluateConservedSites();
        // Find percentage of HM3 identified SNPs the pipeline also identified.
        double power= evaluateHM3Sites();
        System.out.printf("Power:%g TypeIError:%g RatioP/E:%g %n",power,typeIError,power/typeIError); 
        
        // Get position list
        TagData tagData =new TagDataSQLite(GBSConstants.GBS_GBS2DB_FILE);
        PositionList positsV2 = tagData.getSNPPositions();
        System.out.println("LCJ - number of SNPs in SNPPosition table: " + positsV2.size());
        ((TagDataSQLite)tagData).close();

        
        // compare the Postions 
        PositionList positsV1 = ExpectedGenosV1.positions();
        System.out.println("SNP Positions v1: " + positsV1.size() + ", SNP Positions v2: " + positsV2.size() );
        System.out.println("Intersection:"+positsV1.stream().filter(positsV2::contains).count());
        positsV1.stream()
                .filter(p1 -> !positsV2.contains(p1))
                .limit(100)
                .forEach(System.out::println);

       // assertEquals("Expected and Actual PositionLists are not the same size", PositsV1.size(), PositsV2.size());
        Iterator<Position> V1Iter = positsV1.iterator();
        Iterator<Position> V2Iter = positsV2.iterator();
        
        String inV2notV1 = "/Users/lcj34/notes_files/gbsv2/v2SNPCalling_junitTests/SNPSinV2notV1.txt";
        String inV1notV2 = "/Users/lcj34/notes_files/gbsv2/v2SNPCalling_junitTests/SNPSinV1notV2.txt";
        String inBoth = "/Users/lcj34/notes_files/gbsv2/v2SNPCalling_junitTests/inBoth.txt";

        BufferedWriter fileWriterBoth = null;
        BufferedWriter fileWriterV2Only = null;
        BufferedWriter fileWriterV1Only = null;
        try {
            fileWriterV2Only = new BufferedWriter(new FileWriter(inV2notV1));
            fileWriterV1Only = new BufferedWriter(new FileWriter(inV1notV2));
            fileWriterBoth = new BufferedWriter(new FileWriter(inBoth));
            while(V2Iter.hasNext()) {
                Position pos = V2Iter.next();
                Chromosome chrom = pos.getChromosome();
                String chromPos = chrom.getName() + ":" + Integer.toString(pos.getPosition());
                if (!(positsV1.contains(pos))) {
                    fileWriterV2Only.write(chromPos);
                    fileWriterV2Only.write("\n");
                } else {
                    
                    fileWriterBoth.write(chromPos);
                    fileWriterBoth.write("\n");
                }
            }
            
            while(V1Iter.hasNext()) {
                Position pos = V1Iter.next();
                Chromosome chrom = pos.getChromosome();
                String chromPos = chrom.getName() + ":" + Integer.toString(pos.getPosition());
                if (!(positsV2.contains(pos))) {
                    fileWriterV1Only.write(chromPos);
                    fileWriterV1Only.write("\n");
                }
            }   
            fileWriterV2Only.close();
            fileWriterV1Only.close();
            fileWriterBoth.close();
           
        }catch(IOException e) {
            System.out.println(e);
        }        
 
     }

    // THis method no longer needed for basic runs.  I have moved the
    // file created from bowtie into the GBS expected dir under 
    //  bowtie2.1/tagsForAlign910auto.sam
    // This allows users without bowtie to run these tests
    private static void runBowtie() throws Exception{

        // NOTE: when comparing to V1, we need to use bowtie with AGPv2, not AGPv3
//        String cmd="/Users/lcj34/development/bowtie2-2.2.3/bowtie2 " +
//               // " --very-sensitive-local -x /Users/lcj34/development/AGPv2/ZMAGP2_910 " +
//               "-M 4 -x /Users/lcj34/development/AGPv2/ZMAGP2_910 " +
          // THis one if AGPv3
        String cmd="/Users/lcj34/development/bowtie2-2.2.8/bowtie2 " +
               "-M 4 -x /Users/lcj34/development/AGPv3/ZmB73_AGPv3 " +
//        String cmd="/Users/edbuckler/NextGen/AGPv3/bowtie2-2.2.3/bowtie2 " +
//                  "-M 4 -x /Users/edbuckler/NextGen/AGPv3/ZmAGP3_910 " +
                "-U " + Paths.get(outTagFasta).toAbsolutePath().toString() + " "+
                "-S " + Paths.get(inTagSAM).toAbsolutePath().toString();
               // "-S " + Paths.get(samFileAGPv2_local).toAbsolutePath().toString(); // used for very-sensitive-local
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
    private static double evaluateConservedSites() throws Exception{
        PositionList conservedList=PositionListIOUtils.readSNPConserveFile(conservedSitesFile);
        System.out.println("conservedList.numberOfSites() = " + conservedList.numberOfSites());
        PositionListBuilder plb=new PositionListBuilder();
        conservedList.stream()
        .forEach(p -> plb.add(getRandomOffsetSite(conservedList,p)));
        PositionList shiftList=plb.build();
        TagData tagData =new TagDataSQLite(GBSConstants.GBS_GBS2DB_FILE);
        long overlap=tagData.getSNPPositions().stream().filter(p -> conservedList.contains(p)).count();
        System.out.println("overlap = " + overlap);
        long overlapShift=tagData.getSNPPositions().stream().filter(p -> shiftList.contains(p)).count();
        System.out.println("overlapShift = " + overlapShift);
        double typeIerror=(double)overlap/conservedList.numberOfSites();
        System.out.println("Conserved site Type I error = " + typeIerror);

        tagData.getSNPPositions().stream()
        .filter(conservedList::contains)
        .forEach(p -> System.out.println(p.toString()));
        ((TagDataSQLite)tagData).close();
        return typeIerror;
    }

    private static void evaluatePowerOfConservedSites() throws Exception{
        //todo evaluate conserved sites do not vary
        PositionList conservedList=PositionListIOUtils.readSNPConserveFile(conservedSitesFile);
        System.out.println("conservedList.numberOfSites() = " + conservedList.numberOfSites());
        PositionListBuilder plb=new PositionListBuilder();
        conservedList.stream()
        .forEach(p -> plb.add(getRandomOffsetSite(conservedList,p)));
        PositionList shiftList=plb.build();
        TagData tagData =new TagDataSQLite(GBSConstants.GBS_GBS2DB_FILE);
        long overlap=tagData.getSNPPositions().stream().filter(p -> conservedList.contains(p)).count();
        System.out.println("overlap = " + overlap);
        long overlapShift=tagData.getSNPPositions().stream().filter(p -> shiftList.contains(p)).count();
        System.out.println("overlapShift = " + overlapShift);
        double typeIerror=(double)overlap/conservedList.numberOfSites();

        tagData.getSNPPositions().stream()
        .filter(conservedList::contains)
        .forEach(p -> System.out.println(p.toString()));
        
        ((TagDataSQLite)tagData).close(); 
    }

    private static double evaluateHM3Sites() throws Exception{
        PositionList hm3List=PositionListIOUtils.readSNPConserveFile(HM3SitesFile);
        System.out.println("hm3List.numberOfSites() = " + hm3List.numberOfSites());
        PositionListBuilder plb=new PositionListBuilder();
        hm3List.stream()
        .forEach(p -> plb.add(getRandomOffsetSite(hm3List,p)));
        PositionList shiftList=plb.build();
        TagData tagData =new TagDataSQLite(GBSConstants.GBS_GBS2DB_FILE);
        long overlap=tagData.getSNPPositions().stream().filter(p -> hm3List.contains(p)).count();
        System.out.println("hm3List overlap = " + overlap);
        long overlapShift=tagData.getSNPPositions().stream().filter(p -> shiftList.contains(p)).count();
        System.out.println("hm3List overlapShift = " + overlapShift);
        double power=(double)overlap/(double)hm3List.numberOfSites();
        System.out.println("Power of finding HM3 sites:" + power);

        ((TagDataSQLite)tagData).close();
        return power;
    }

    private static Position getRandomOffsetSite(PositionList hm3List, Position initPosition) {
        int offset=256;
        while(true) {
            int totalOffset=offset+random.nextInt(offset)-(offset/2);
            Position shiftPosition=new GeneralPosition.Builder(initPosition.getChromosome(), initPosition.getPosition() + totalOffset).build();
            if(!hm3List.contains(shiftPosition)) return (shiftPosition);
        }
    }

}
