/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import static org.junit.Assert.*;

import java.io.IOException;
import java.lang.reflect.Method;
import java.lang.reflect.Type;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.Ignore;
import org.junit.Test;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Tuple;

/**
 * Aligns the {@link RepGenSimData} tags against its reference with both RepGen aligners. Every amplicon
 * should yield exactly one reference tag, and every tag should align best, on the forward strand, to
 * the reference tag of its own amplicon starting at the amplicon start.
 * <p>
 * The ignored tests below are the original developer runs on the ANP68 rAmpSeq data and the full AGPv4
 * genome. They are kept for reference until #43 decides whether the RepGen pipeline stays.
 * 
 * @author lcj34
 *
 */
public class RepGenAlignerPluginTest {

    private static final int REF_KMER_LEN = 300;
    private static final int PERFECT_SCORE = 2 * RepGenSimData.KMER_LENGTH;
    private static final int ONE_MISMATCH_SCORE = 2 * (RepGenSimData.KMER_LENGTH - 1) - 1;

    private static final String nextGenDir = GeneralConstants.DATA_DIR + "RepGen/";
    private static final String repGenInputDir = nextGenDir + "single_fastq/";
    private static final String referenceGenome = GeneralConstants.DATA_DIR + "CandidateTests/Zea_mays.AGPv4.dna.toplevel.fa.gz";
    private static final String repGenDB_FILE = GeneralConstants.TEMP_DIR + "RepGen/repGen_testOrigAligner.db";
    private static final String repGen_TESTING_KEY_FILE = nextGenDir + "repGen_keyFiles/anp68_locus9_key.txt";

    /**
     * Test of performFunction method, of class RepGenAlignerPlugin.
     */
    @Test
    public void testRepGenAlignerPlugin() throws Exception {
        RepGenSimData data = RepGenSimData.createUnder("aligner");
        String db = RepGenLoadSeqToDBPluginTest.loadDB(data, 5);

        new RepGenAlignerPlugin()
                .inputDB(db)
                .refGenome(data.referenceFile.toString())
                .seedLen(17)
                .kmerLen(RepGenSimData.KMER_LENGTH)
                .minTagCount(1)
                .minHitCount(5)
                .refKmerLen(REF_KMER_LEN)
                .performFunction(null);

        RepGenSQLite repGenData = new RepGenSQLite(db);
        try {
            List<RefTagData> refTags = new ArrayList<>(repGenData.getRefTags());
            assertEquals("One reference tag per amplicon", data.loci.size(), refTags.size());
            Map<RepGenSimData.Locus, RefTagData> refTagByLocus = new HashMap<>();
            for (RepGenSimData.Locus locus : data.loci) {
                List<RefTagData> covering = refTags.stream()
                        .filter(rtd -> rtd.chromosome().equals(locus.chromosome))
                        .filter(rtd -> rtd.tag().sequence().contains(locus.tag(false)))
                        .collect(Collectors.toList());
                assertEquals("Reference tags covering the amplicon at " + locus.chromosome + ":" + locus.start, 1, covering.size());
                RefTagData refTag = covering.get(0);
                String sequence = refTag.tag().sequence();
                assertEquals(REF_KMER_LEN, sequence.length());
                assertEquals("Reference tag is the reference sequence at its position",
                        data.referenceSequence(locus.chromosome, refTag.position(), refTag.position() + sequence.length() - 1), sequence);
                refTagByLocus.put(locus, refTag);
            }

            List<Tag> tags = new ArrayList<>(repGenData.getTags());
            int numTags = tags.size();
            assertEquals(data.expectedTags.size(), numTags);
            assertEquals("Every tag is aligned against every other tag",
                    numTags * (numTags - 1), repGenData.getAllNonRefTagAlignments(0).size());
            assertEquals("Every reference tag is aligned against every other reference tag",
                    refTags.size() * (refTags.size() - 1), repGenData.getAllRefTagAlignments(0).size());
            assertAlleleIsBestTagMatch(data, repGenData.getTagAlignmentsForTags(tags, 0));

            Multimap<Tag, AlignmentInfo> refAlignments = repGenData.getRefAlignmentsForTags(tags, 0);
            for (Tag tag : tags) {
                Collection<AlignmentInfo> alignments = refAlignments.get(tag);
                assertEquals("Each tag is aligned to both strands of every reference tag", 2 * refTags.size(), alignments.size());
                RepGenSimData.ExpectedTag expected = data.expectedTags.get(tag.sequence());
                RefTagData ownRefTag = refTagByLocus.get(expected.locus);
                AlignmentInfo best = alignments.stream().max(Comparator.comparingInt(AlignmentInfo::score)).get();
                assertEquals(ownRefTag.tag(), best.tag2());
                assertEquals(ownRefTag.position(), best.tag2pos());
                assertEquals("Forward strand", 1, best.ref_strand());
                assertEquals("Alignment starts at the amplicon start", expected.locus.start, best.alignmentPos());
                assertEquals(expected.alt ? ONE_MISMATCH_SCORE : PERFECT_SCORE, best.score());
            }
        } finally {
            repGenData.close();
        }
    }

    @Test
    public void testRepGenPhase2AlignerPlugin() throws Exception {
        RepGenSimData data = RepGenSimData.createUnder("phase2Aligner");
        String db = RepGenLoadSeqToDBPluginTest.loadDB(data, 5);
        Path primers = Files.write(data.dir.resolve("primers.txt"), List.of(
                "forward\treverse",
                "GATCGGAAGAGCACACGTCT\tAATGATACGGCGACCACCGA",
                RepGenSimData.FORWARD_PRIMER + "\t" + RepGenSimData.REVERSE_PRIMER));

        runPhase2(data, db, primers);

        RepGenSQLite repGenData = new RepGenSQLite(db);
        try {
            Set<RefTagData> refTags = repGenData.getRefTags();
            assertEquals("One reference tag per amplicon", data.loci.size(), refTags.size());
            for (RepGenSimData.Locus locus : data.loci) {
                List<RefTagData> atLocus = refTags.stream()
                        .filter(rtd -> rtd.chromosome().equals(locus.chromosome) && rtd.position() == locus.start)
                        .collect(Collectors.toList());
                assertEquals("Reference tag starting at the forward primer of " + locus.chromosome + ":" + locus.start, 1, atLocus.size());
                String sequence = atLocus.get(0).tag().sequence();
                assertTrue(sequence.startsWith(locus.refAmplicon));
                assertEquals(data.referenceSequence(locus.chromosome, locus.start, locus.start + REF_KMER_LEN), sequence);
            }
        } finally {
            repGenData.close();
        }
    }

    @Test
    public void testRepGenPhase2AlignerPluginPrimersNotFound() throws Exception {
        RepGenSimData data = RepGenSimData.createUnder("phase2AlignerNoPrimers");
        String db = RepGenLoadSeqToDBPluginTest.loadDB(data, 5);
        Path primers = Files.write(data.dir.resolve("primers.txt"), List.of(
                "forward\treverse",
                "GATCGGAAGAGCACACGTCT\tAATGATACGGCGACCACCGA"));

        runPhase2(data, db, primers);

        RepGenSQLite repGenData = new RepGenSQLite(db);
        try {
            assertTrue("Seed hits without a primer pair nearby make no reference tags", repGenData.getRefTags().isEmpty());
        } finally {
            repGenData.close();
        }
    }

    private static void runPhase2(RepGenSimData data, String db, Path primers) {
        new RepGenPhase2AlignerPlugin()
                .inputDB(db)
                .refGenome(data.referenceFile.toString())
                .seedLen(31)
                .kmerLen(RepGenSimData.KMER_LENGTH)
                .minTagCount(1)
                .refKmerLen(REF_KMER_LEN)
                .primers(primers.toString())
                .performFunction(null);
    }

    /** Each tag's highest-scoring tag alignment is the other allele of its amplicon, one mismatch apart. */
    static void assertAlleleIsBestTagMatch(RepGenSimData data, Multimap<Tag, AlignmentInfo> tagAlignments) {
        for (Tag tag : tagAlignments.keySet()) {
            RepGenSimData.ExpectedTag expected = data.expectedTags.get(tag.sequence());
            AlignmentInfo best = tagAlignments.get(tag).stream().max(Comparator.comparingInt(AlignmentInfo::score)).get();
            assertEquals(expected.locus.tag(!expected.alt), best.tag2().sequence());
            assertEquals(ONE_MISMATCH_SCORE, best.score());
        }
    }

    @Ignore("#43: duplicate of RepGenAlignerUnitTest.testCreateKmerSeedsFromDBTags")
    @Test
    public void testCreateKmerSeedsFromDBTags() {
        // Create list of tags.  Send this list into the DB.
        // Verify we get the tags we expect.
        String seq1 = "ACGTCGATCTAGGGGGTCTCGACGAAGGCAACCATCTTGGTATGCCGGAGGACGGTGATCCTCCTAGGCCTGCGCCTCGCGTTGACATCCTTCGGGAGCTAGCTGTGGTCCCAGTCCCTGCGAGGGGTCAGGACGCATAGCTCGAGCAAA";
        String seq2 = "ACGTCGATCTAGGGGGTCTCGACGAAGGAAACCATTTTGGTATGCCGGAGGACGGTGATCCCCCTAGGCCCGCGCCTCGCGTTGACATCCTTTGGGAGCTAGCTATGGTCCCAGTCCCTGTGGGGGGTCAGGACGACGCACAGCTCGAGC";
        String seq3 = "GCACAAGTTGTCCTGCTTCCTCGTCGAGCCTGGCGTGCATCTCGCGGATTTTCTCGAGCTGTGCGTCCTGACCCCCCGCAGGGACTGGGACCACAGCTAGCTCCCGAAGGATGTCAACGCGAGGCGCAGGCCTAGAGGGATCGCCGTCCT";
        String seq4 = "GCACAAGTTGTCCTGCTTCCTCGTCGAGCCTGGCCTGTACCTCGCGGATTTGCTCGAGCTGTGCGTCCTGACCTTCCGCAGGGACTGGGACCACAGCTAGCTCCCGAAGGATGTCAACGCGAGGGGCAGACCTAGGGGGATCACCGTCCT";
        String seq5 = "GCACAAGTTGTCCTGCTTCCTCGTCGAGCTTGGTCTGCATCTCACGGATTTGCTCGAGCTGTGTGTCCTGACCCCCCGCAGGGACTGGGACCACAGCTAGCTCCCGAAGGATGTCAACACGAGGTGCAGGCCTAGGGGGATCGTCGTCCT";
        //String seq6 = "GCACAGGTTGTCCTGCTTCCNNGTCGAGCTTGGTCTGCATATCACGGATTTGCTCGAGCTGTGGGTCCTGACCCCCCGCAGGGACTGGGACCACAGCTAGCTCCCGAAGGATGTCAACACGAGGTGCAGGCCTAGGGGGATCGTCGTCCT";

        // Sequences with null:  While the RepGenALignerPlugin method checks for N's in the sequence,
        // there won't be any be cause the sequenes are all of type object Tag.  The Tag builder returns
        // null if there are N's in the sequence.  So it will never have gotten onto the tag list.
        Set<String> sequences = new HashSet<String>();
        sequences.add(seq1);
        sequences.add(seq2);
        sequences.add(seq3);
        sequences.add(seq4);
        sequences.add(seq5);
       // sequences.add(seq6);

        // create list of kmers
        Multimap<String,String> kmerSequenceMap = HashMultimap.create();
        int window = 20;
        int seedlen = 17;
        for (String seq: sequences) {
            int idx = 0;
            for (idx = 0; idx < seq.length()-window;){
                // Add the shortened kmer and the reverse Complement of it
                String kmer = seq.substring(idx, idx+seedlen);
                kmerSequenceMap.put(kmer, seq);
                byte[] kmerRC = NucleotideAlignmentConstants.reverseComplementAlleleByteArray(kmer.getBytes());
                String kmerRCString = new String(kmerRC);
                if (kmerRCString.contains("N")) {
                    idx+= window;
                    continue;
                }
                kmerSequenceMap.put(kmerRCString, seq);
                idx += window;
            }
        }

        Set<Tag> tags = new HashSet<Tag>();
        tags.add(TagBuilder.instance(seq1).build());
        tags.add(TagBuilder.instance(seq2).build());
        tags.add(TagBuilder.instance(seq3).build());
        tags.add(TagBuilder.instance(seq4).build());
        tags.add(TagBuilder.instance(seq5).build());
        //tags.add(TagBuilder.instance(seq6).build());


        // Call RepGenAlignerPlugin method to create kmers
        Multimap<String,Tag> kmerTagMap = HashMultimap.create();     
        RepGenAlignerPlugin repGenAlignerPlugin = new RepGenAlignerPlugin();
        Class theClass = repGenAlignerPlugin.getClass();
        try {
            Method theMethod = theClass.getDeclaredMethod("createKmerSeedsFromDBTags", new Class[] {Set.class, Multimap.class, int.class} );
            theMethod.setAccessible(true);
            // Ideally, I'd set the plugin parameters to match what is in the test case.
            // here is an example of how to set private values:
            //  http://stackoverflow.com/questions/29901351/can-we-mock-private-fields-of-a-class-for-junit-testing
            // Unfortunately, This does not work well when your type is PluginParamter<>.  So the test
            // will ensure the seedLen matches the default for the plugin parameter.  NOte that
            // seedWindow is also a plugin parameter, but the createkmerSeedsFromDBTags method has
            // that passed as a parameter, so I can pass that here.
            theMethod.invoke(repGenAlignerPlugin,tags,kmerTagMap,20);

        } catch (Exception exc) {
            exc.printStackTrace();
        }
        System.out.println("testCreateKmerSeedsFromDBTags: test kmer size: " + kmerSequenceMap.keySet().size() +
                ", createKmerSeedsFromDBTags size: " + kmerTagMap.keySet().size());
        assertEquals(kmerSequenceMap.keySet().size(),kmerTagMap.keySet().size());

        // Verify the same kmers are in each string
        for (String kmer: kmerSequenceMap.keySet()) {
            if (!(kmerTagMap.keySet().contains(kmer))) {
                System.out.println("Missing from create method:" + kmer);
                fail("createKmerSeedsFromDBTags does not contain kmer from junit");
            }
        }
        for (String kmer: kmerTagMap.keySet()) {
            if (!(kmerSequenceMap.keySet().contains(kmer))) {
                System.out.println("In createKmerSeedsFromDBTags array, not in junit array:" + kmer);
                fail("createKmerSeedsFromDBTags contain kmer not in junit array");
            }
        }      
    }

    @Ignore("#43: needs the ANP68 reads and the full AGPv4 genome, which are not in tassel_test_data; "
            + "testRepGenAlignerPlugin stores and reads back tag-tag alignments on RepGenSimData")
    @Test
    public void testPutGetTagTagAlignment() {
        // This test calls RepGenALignerPlugin method to put tags to the
        // database, then calls getAlignmentsForTag to get them back out.a

        // Needs to be augmented to put/get nonref to reference tags in/out.
        LoggingUtils.setupDebugLogging();

        try{
            Files.deleteIfExists(Paths.get(repGenDB_FILE));
        } catch (IOException e) {
            e.printStackTrace();
        }
        long time=System.nanoTime();
        System.out.println();

        System.out.println("testRepGenAlignerPlugin, calling RepGenLoadSeqToDBPlugin to load the db ....");
        // First load the DB
        new RepGenLoadSeqToDBPlugin()
        .inputDirectory(repGenInputDir)
        .outputDatabaseFile(repGenDB_FILE)
        .keyFile(repGen_TESTING_KEY_FILE)
        .kmerLength(150)
        .minKmerCount(1) // should be higher when running for real, this is for running just 1 file
        //.deleteOldData(true)
        .performFunction(null);
        System.out.println("TotalTime for repGenLoadSeqToDBPlugin was " + (System.nanoTime() - time) / 1e9 + " seconds");

        RepGenDataWriter repGenData=new RepGenSQLite(repGenDB_FILE);
        Set<Tag> tagsToAlign = repGenData.getTags();
        List<Tag> tagList = new ArrayList<Tag>(tagsToAlign);
        List<Tag> shortTagList = new ArrayList<Tag>();

        int numTags = 0;
        System.out.println("LCJ - first 10 tags from db:");
        while(numTags < 10) { // want to test with a small number
            shortTagList.add(tagList.get(numTags));
            System.out.println(tagList.get(numTags).sequence());
            numTags++;
        }
        System.out.println("LCJ - created shortTagLIst, size is : " + shortTagList.size());
        Multimap<Tag,Tuple<Tag,Integer>> tagTagAlignMap = Multimaps.synchronizedMultimap(HashMultimap.<Tag,Tuple<Tag,Integer>>create());
        time = System.nanoTime();
        int totalCount = 0;

        repGenData.addReferenceGenome(referenceGenome); // add reference genome - needed to be present for test below
        RepGenAlignerPlugin repGenLoadSeqToDBPlugin = new RepGenAlignerPlugin();
        Class theClass = repGenLoadSeqToDBPlugin.getClass();
        try {
            Method theMethod = theClass.getDeclaredMethod("calculateTagTagAlignment", new Class[] {List.class, Multimap.class} );
            theMethod.setAccessible(true);
            Method[] allMethods = theClass.getDeclaredMethods();
            for (Method m : allMethods) {
                if (m.getName().equals("calculateTagTagAlignment")) {
                    Type[] pType = m.getParameterTypes();
                    for (Type pt : pType) {
                        System.out.println("  ptype: " + pt.toString());
                    }
                }
  //              System.out.println("Method name: " + m.getName() + ", parameters: " + m.getParameters() +
  //                       ", parameter types: " + m.getParameterTypes());
            }

            Multimap<Tag,AlignmentInfo> tagAlignInfoMap = Multimaps.synchronizedMultimap(HashMultimap.<Tag,AlignmentInfo>create());
            theMethod.invoke(repGenLoadSeqToDBPlugin,(Object)shortTagList,tagAlignInfoMap);
            // tag1_isref, tag2_isref, tag1chrom, tag1pos, refGenome
            repGenData.putTagTagAlignments(tagAlignInfoMap);

            // Now - verify we can get them back out!
            Multimap<Tag,AlignmentInfo> returnTagAlignInfoMap = repGenData.getTagAlignmentsForTags(shortTagList,0);
            System.out.println("Size of tagAlignInfoMap: " + tagAlignInfoMap.size() 
            + ", size of returnTagAlignInfoMap: " + returnTagAlignInfoMap.size());

            assertEquals(tagAlignInfoMap.keySet().size(),returnTagAlignInfoMap.keySet().size());

        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } 
        try {
            ((RepGenSQLite)repGenData).close();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    @Ignore("#43: needs the ANP68 reads and the full AGPv4 genome, which are not in tassel_test_data; "
            + "testRepGenAlignerPlugin checks the alignment counts on RepGenSimData")
    @Test
    public void testGetTagAndRefTagAlignments() {
        // This loads the db with tags and alignments,
        // then calls getALignmentsForTag() to retrieve
        // tag alignments for 10 non-ref tags from the db..  
        // There will be both tag-tag and tag-refTag alignments in the  returned batch

        //When run with the initial files loaded here, using just
        // chrom9 (IE - uncomment debug line 
        //    if (chrom.getChromosomeNumber() != 9) return; 
        // in RepGenAlignerPlugin:processData() to test with just chrom9)
        // there are 3051 tags in the db, and 551 reference tags
        // created.  Tags are not aligned against themselves, so for each
        // tag there should be 3050 + 1102 alignments, equals 4152 alignments
        // per tag, or a total of 41520 alignments returned for the 10 tags.

        LoggingUtils.setupDebugLogging();

        try{
            Files.deleteIfExists(Paths.get(repGenDB_FILE));
        } catch (IOException e) {
            e.printStackTrace();
        }
        long time=System.nanoTime();
        System.out.println();

        System.out.println("testRepGenAlignerPlugin, calling RepGenLoadSeqToDBPlugin to load the db ....");
        // First load the DB
        new RepGenLoadSeqToDBPlugin()
        .inputDirectory(repGenInputDir)
        .outputDatabaseFile(repGenDB_FILE)
        .keyFile(repGen_TESTING_KEY_FILE)
        .kmerLength(150)
        .minKmerCount(1) // should be higher when running for real, this is for running just 1 file
        //.deleteOldData(true)
        .performFunction(null);
        System.out.println("TotalTime for repGenLoadSeqToDBPlugin was " + (System.nanoTime() - time) / 1e9 + " seconds");

        System.out.println("\ntestRepGenAlignerPlugin ... calling RepGenAlignerPlugin");
        time=System.nanoTime();

        // Align against a reference genome
        new RepGenAlignerPlugin()
        .inputDB(repGenDB_FILE)
        .refGenome(referenceGenome)
        .seedLen(17)
        .kmerLen(150)
        .minTagCount(1)
        .minHitCount(10)
        .refKmerLen(300)
        .performFunction(null);
  System.out.println("TotalTime for RepGenLoadSeqToDBPlugin: " + (System.nanoTime() - time) / 1e9 + " seconds");

        RepGenDataWriter repGenData=new RepGenSQLite(repGenDB_FILE);
        Set<Tag> tagsToAlign = repGenData.getTags();
        List<Tag> tagList = new ArrayList<Tag>(tagsToAlign);
        List<Tag> shortTagList = new ArrayList<Tag>();

        Set<RefTagData> refTags = repGenData.getRefTags();
        int refTagSize = refTags.size();
        int tagSize = tagsToAlign.size();
        refTags = null;
        tagsToAlign = null;

        int numTags = 0;
        System.out.println("LCJ - first 10 tags from db:");
        while(numTags < 10) { // want to test with a small number
            shortTagList.add(tagList.get(numTags));
            System.out.println(tagList.get(numTags).sequence());
            numTags++;
        }
        System.out.println("LCJ - created shortTagLIst, size is : " + shortTagList.size());
        Multimap<Tag,Tuple<Tag,Integer>> tagTagAlignMap = Multimaps.synchronizedMultimap(HashMultimap.<Tag,Tuple<Tag,Integer>>create());
        time = System.nanoTime();

        repGenData.addReferenceGenome(referenceGenome); // add reference genome - needed to be present for test below
        RepGenAlignerPlugin repGenLoadSeqToDBPlugin = new RepGenAlignerPlugin();
        try {
            // Now - verify we can get them back out!
            System.out.println("\n Calling repGenData.getAlignmentsForTags with shortTagList .... ");
            Multimap<Tag,AlignmentInfo> returnTagAlignInfoMap = repGenData.getTagAlignmentsForTags(shortTagList,0);
            System.out.println("Size of returnTagAlignInfoMap: " + returnTagAlignInfoMap.size());
            System.out.println("TotalTime for getALignmentsForTags: " + (System.nanoTime() - time) / 1e9 + " seconds");

            // When run with the initial files loaded here, using just
            // chrom9 (uncomment debug line in RepGenAlignerPlugin to test with
            // just chrom9) there are 3051 tags in the db, and 551 reference tags
            // created.  Tags are not aligned against themselves, so for each
            // tag there should be 3050 + 1102 alignments, equals 4152 alignments
            // per tag, or a total of 41520 alignments returned for the 10 tags.

            // This should be computed programmatically - get the nubmer of tags and
            // the number of reference tags from the DB, then do the math and compare.
            // Will add that later.
            int numShortTags = shortTagList.size();
            int expectedAlignments = ((tagSize * numShortTags) - numShortTags) + (numShortTags * (refTagSize * 2));
            assertEquals(returnTagAlignInfoMap.size(),expectedAlignments);

        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } 

        try {
            time = System.nanoTime();
            // Test grabbing all the non-reference tag alignments.

            refTags = null; // clear some memory
            System.out.println("\n Calling repGenData.getAllNonRefTagAlignments - getting all non-ref tag alignments");
            Multimap<Tag,AlignmentInfo> returnTagAlignInfoMap = repGenData.getAllNonRefTagAlignments(0);          
            System.out.println("TotalTime for getAllNonRefTagAlignments: " + (System.nanoTime() - time) / 1e9 + " seconds");

            // Total alignments is number of non-ref tags times number of non-ref Tags minus the 
            // number of tags (because all tags are aligned against all others, but not against themselves)
            // plus twice the number of reftags times the number of tags.  Value is twice the number
            // of reftags because we aligned against both the refTag and the refTag reverse complement.
            int totalAlignments = (tagSize * tagSize) - tagSize + (refTagSize * 2 * tagSize);
            System.out.println("RefTagAlianmgnets: refTagSize: " + tagSize +", totalExpectedAlignments: " 
               + totalAlignments + ", size of returnTagAlignInfoMap: " + returnTagAlignInfoMap.size());
            assertEquals(returnTagAlignInfoMap.size(),totalAlignments);
            returnTagAlignInfoMap = null; // clearing memory - shouldn't need as this is within a block

        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } 

        try {
            time = System.nanoTime();

            // Test grabbing the reference tags.
            System.out.println("\n Calling repGenData.getAllRefTagAlignments- getting all refTag alignments");
            Multimap<RefTagData,AlignmentInfo> returnRefTagAlignInfoMap = repGenData.getAllRefTagAlignments(0);
            System.out.println("TotalTime for getAllRefTagAlignments: " + (System.nanoTime() - time) / 1e9 + " seconds");

            // The number of refTags should be (refTagSize * refTagSize) - refTagSize.  This is because
            // the ref tags were aligned against all other reftags, but not against themselves.
            int totalAlignments = (refTagSize * refTagSize) - refTagSize;
            System.out.println("RefTagAlianmgnets: refTagSize: " + refTagSize +", totalExpectedAlignments: " 
               + totalAlignments + ", size of returnRefTagAlignInfoMap: " + returnRefTagAlignInfoMap.size());
            assertEquals(returnRefTagAlignInfoMap.size(),totalAlignments);
            returnRefTagAlignInfoMap = null; // clearing memory

        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } 

        try {
            ((RepGenSQLite)repGenData).close();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    @Ignore("#43: needs a developer's Locus9 R1 RepGen database, which is not in tassel_test_data; "
            + "an unfinished sketch that asserts nothing")
    @Test
    public void testRefTagCreationAlgrithms() {
        String myDBFile = nextGenDir + "dbs/repGen_Locus9R1_rs.db";
        RepGenDataWriter repGenData=new RepGenSQLite(myDBFile);
        Set<Tag> tagsToAlign = repGenData.getTags();
        List<Tag> tagList = new ArrayList<Tag>(tagsToAlign);

        // get list of tags with references having minimum quality of 291 (3 bp mismatch/gap or less)
        // 291 assumes aligner was run with 2,-1.-1 and tag length=150;
        // if tag length < 150, will likely be tossed.  A perfect score would keep tags
        // of length 144 and up, but otherwise they are dropped.
        int minscore = 291;
        Multimap<Tag, AlignmentInfo> tagRefTagsMap = repGenData.getRefAlignmentsForTags(tagList,  minscore);

        // Above we got a list of tags with reftag alignment info
        // below we get list of refTags and the tags that align to them.
        // Need to write this up inSQlite -look for alignments were tag2=ref
        // and tag1 is non-ref.  then check score.
    }
}
