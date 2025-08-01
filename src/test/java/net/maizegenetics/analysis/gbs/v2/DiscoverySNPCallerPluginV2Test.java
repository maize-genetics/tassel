package net.maizegenetics.analysis.gbs.v2;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Multimap;
import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;

import junit.framework.Assert;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.GenomeSequence;
import net.maizegenetics.dna.map.GenomeSequenceBuilder;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TagData;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TaxaDistBuilder;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Tuple;

import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.Assert.*;

public class DiscoverySNPCallerPluginV2Test {
    private Map<Tag,TaxaDistribution> tagTaxaDistMap;
    private Map<Tag,TaxaDistribution> tagTaxaDistMapR;
    private Map<Tag, TaxaDistribution> tagTaxaDistMap2;
    private Map<Tag, TaxaDistribution> tagTaxaDistMap2R;
    List<Boolean> direction= ImmutableList.of(true, true, true, true,false);
    public static final String REFERENCE_GENOME_EXPECTED_DIR = "dataFiles/CandidateTests/";
    public static final String REFERENCE_GENOME_FILE = REFERENCE_GENOME_EXPECTED_DIR + "ZmB73_RefGen_v3_chr9_10_1st25MB.fasta.gz";

    private final static String outTagFasta=GBSConstants.GBS_TEMP_DIR + "tagsForAlign.fa.gz";
   // private final static String inTagSAM=GBSConstants.GBS_TEMP_DIR + "tagsForAlign910auto.sam"; // used if running bowtie during test
    private final static String inTagSAM=GBSConstants.GBS_EXPECTED_DIR + "bowtie2.1/tagsForAlign910auto.sam";
    @Before
    public void setUp() throws Exception {
        LoggingUtils.setupLogging();
        // WIth changes made to alignTags from TAS-1001, the sequences below are changes
        // to all be forward sequences.  alignTags previously looked at the strand direction
        // and converted reverse oriented strands to forward, then aligned.  
        // DiscoverySNPCallerPluginV2:alignTags() algorightm now takes all tags compiled in a
        // single direction and aligns them.  Forward and Reverse strands are aligned separately.
        
        // Test cases need to be added to include reverse strand examples
        // The initial reverse strand of TGTCGTCGT was translated to forward strand of ACGACGACA
        // via the Tag method tag.toReverseCompliment().
        List<String> sequences= ImmutableList.of("ACGACGACG","ACtACGACG","ACGACGtCG","ACtACGGtCG",
                "ACGACGACA"); // originally last tag was reverse and was "TGTCGTCGT"
        List<String> sequences2= ImmutableList.of("ACGACGACGTT","ACtACGACGTT","ACGACGtCGTT","AtACGGtCGAAGG",
                "ACGACGACA"); // originally last tag was reverse and was "TGTCGTCGT"
        //List<Boolean> direction= ImmutableList.of(true, true, true, true,false);
        List<TaxaDistribution> taxaDists=ImmutableList.of(
                TaxaDistBuilder.create(10,new int[]{0,1,2,3},new int[]{1,2,1,3}),
                TaxaDistBuilder.create(10,new int[]{4,5},new int[]{1,1}),
                TaxaDistBuilder.create(10,new int[]{3,6,7},new int[]{3,2,1}),
                TaxaDistBuilder.create(10,8),
                TaxaDistBuilder.create(10,1)
        );
        
        List<Tag> tags= sequences.stream().map(s -> TagBuilder.instance(s).build()).collect(Collectors.toList());
        tags.set(0,TagBuilder.instance(tags.get(0).sequence()).reference().build());  //make first tag the reference tag
        tagTaxaDistMap=new HashMap<>();
        for (int i = 0; i < tags.size(); i++) {  //This should all be doable with a Stream zip function, but it was removed from beta versions, and I don't know how to do it now.
            tagTaxaDistMap.put(tags.get(i),taxaDists.get(i));
        }
        
        // Create second tagTaxaDistMap for testing filtered alignments.  new sequences, first one is
        // still the reference, same taxDist for maps 1 and 2
        List<Tag> tags2= sequences2.stream().map(s -> TagBuilder.instance(s).build()).collect(Collectors.toList());
        tags2.set(0,TagBuilder.instance(tags2.get(0).sequence()).reference().build());  //make first tag the reference tag
        tagTaxaDistMap2=new HashMap<>();
        for (int i = 0; i < tags2.size(); i++) {  //This should all be doable with a Stream zip function, but it was removed from beta versions, and I don't know how to do it now.
            tagTaxaDistMap2.put(tags2.get(i),taxaDists.get(i));
        }
    }

    @Test
    public void testAlignTags() throws Exception {
        Map<Tag,String> alignedTags= usePrivAlignTags(tagTaxaDistMap);
        alignedTags.forEach((t,s) -> System.out.println(s+ ": from "+t.toString()+ tagTaxaDistMap.get(t).toString()));
        // 10 is the length of sequence after aligning, e.g. ACGACGACG becomes ACGAC-GACG
        alignedTags.forEach((t, s) -> assertEquals("Alignment Length incorrect", 10, s.length()));
        Position startPos=new GeneralPosition.Builder(new Chromosome("1"),40000).build();

        Table<Position, Byte, List<TagTaxaDistribution>> tAlign=usePrivConvertAlignmentToTagTable(alignedTags,
                tagTaxaDistMap,  startPos);

        Table<Position, Byte, List<TagTaxaDistribution>> alignT2= TreeBasedTable.create();
        tAlign.rowMap().forEach((p, mBTTD) -> {
            if(mBTTD.entrySet().size()>1) {
                mBTTD.forEach((a,ttdl) -> alignT2.put(p,a,ttdl));
            }
        } );
        alignT2.rowMap().forEach((p,ttd) -> System.out.println(p.getPosition()+" :"+ttd.entrySet().size()));
        System.out.println(tAlign.toString());
        
    }
    
	@Test
    public void testFilterAlignedTags() throws Exception {
		System.out.println("Running test TestFilterAlignedtags");
    	System.out.println("Begin sequence1 tests:\n");
        Map<Tag,String> alignedTagsUnfiltered= usePrivAlignTags(tagTaxaDistMap);
        alignedTagsUnfiltered.forEach((t,s) -> System.out.println(s+ ": from "+t.toString()+ tagTaxaDistMap.get(t).toString()));
        alignedTagsUnfiltered.forEach((t, s) -> assertEquals("Alignment Length incorrect", 10, s.length()));
        Position startPos=new GeneralPosition.Builder(new Chromosome("1"),40000).build();
        
        // The filtered tags and their alignment look as below.  note that all are below 0.5
        // ref tag:     ACGAC-GACG
        // aligned tag: ACGAC-GACG value: 0.0 for  IC=0,NC=9
        // aligned tag: ACTAC-GACG value: 0.0 for  IC=0,NC=9
        // aligned tag: ACTACGGTCG value: 0.1 for  IC=1,NC=9
        // aligned tag: ACGACG-TCG value: 0.2 for  IC=2,NC=8
        // aligned tag: ACGACG-ACA value: 0.2 for  IC=2,NC=8
        Map<Tag,String> alignedTags = usePrivFilterAlignTags(alignedTagsUnfiltered, startPos, 0.5);

        assertEquals(alignedTags.size(),5);
        System.out.println("\nsequence1: Tags successfully aligned within threshold range of 0.5");
        convertAlignedTagsToTable(alignedTags,tagTaxaDistMap, startPos);
        
        // Run again with filter  0.12
        // Same aligned tags as above, the last 2 cause the loci to be tossed
        alignedTags = usePrivFilterAlignTags(alignedTagsUnfiltered, startPos, 0.12); 
        assertEquals(alignedTags,null);
        System.out.println("\nsequence1: Null returned for threshold of 0.12 - tags at loci tossed\n");
        
    	// Repeat tests with second sequence of tags with 0.8 and 0.12 threshold 
        // The filtered tags and their alignments should look as below.
        // ref tag:     A------CGACGACGTT 
        // aligned tag: A------CGACGTCGTT value: 0.0    for IC=0, NC=11
        // aligned tag: A------CGACGACA-- value: 0.1818 for IC=2, NC=9
        // aligned tag: ATACGGTCGAAG--G-- value: 0.5882 for IC=10, NC=7
        // aligned tag: A------CGACGACGTT value: 0.0 for  IC=0,  NC=11
        // aligned tag: A------CTACGACGTT value: 0.0 for  IC=0, NC=11
        System.out.println("Begin sequence2 tests:\n");
        alignedTagsUnfiltered= usePrivAlignTags(tagTaxaDistMap2);
        alignedTagsUnfiltered.forEach((t,s) -> System.out.println(s+ ": from "+t.toString()+ tagTaxaDistMap2.get(t).toString()));
        alignedTagsUnfiltered.forEach((t, s) -> assertEquals("Alignment Length incorrect", 17, s.length()));
        
        //  Run with 0.8 as filter - all tags fall within range
        alignedTags = usePrivFilterAlignTags(alignedTagsUnfiltered, startPos, 0.8);
        assertEquals(alignedTags.size(),5);   // All tags pass threshold 
        System.out.println("\nsequence2: Tags successfully aligned within threshold range of 0.8\n");
        convertAlignedTagsToTable(alignedTags,tagTaxaDistMap2, startPos);
        
        //Run again with filter  0.12 - tags are tossed. the 0.1818 and 0.5882 above exceed the threshold
        alignedTags = usePrivFilterAlignTags(alignedTagsUnfiltered, startPos, 0.12);
        assertEquals(alignedTags,null); // 2 tags fall out, so null is returned.
        System.out.println("\nsequence2: Null returned for threshold of 0.12 - tags at loci tossed");
    }

    @Test
	public void testFullSNPCaller() throws Exception {
	    System.out.println("Running testFullSNPCaller");
	    Chromosome nine = new Chromosome("9");
	    Chromosome ten = new Chromosome("10");
	    new DiscoverySNPCallerPluginV2()
	            .inputDB(GBSConstants.GBS_GBS2DB_FILE)
	            .minMinorAlleleFreq(0.01)
	            .startChromosome(nine)
	            .endChromosome(ten)
	            .performFunction(null);
	}

    @Test
    public void testSinglePositionFiltering() throws Exception {
        System.out.println("DiscoverySNPCallerPluginV2Test.testSinglePositionFiltering");
        Position startPos=new GeneralPosition.Builder(new Chromosome("1"),10000).build();
        DiscoverySNPCallerPluginV2 caller=new DiscoverySNPCallerPluginV2();
        
        // LCJ - this must be changed to pass in the strand info (true=forward, false=reverse)
        // This probably means the data set must change to include only forward in one set, and
        // only backward in the other set.  Jeff has a direction table above with 4 values,
        // 3 trues, 1 false.  Re-work that to create 2 different maps, and then pass each
        // map in separately.
        Multimap<Tag,Allele> alignT=caller.findAlleleByAlignment(startPos,tagTaxaDistMap,new Chromosome("1"), true);
        alignT.asMap().forEach((t, alleleCollection) -> {
                    System.out.print(t.toString());
                    alleleCollection.forEach(a -> System.out.print(a.position().getPosition() + "=" + a.allele() + ","));
                    System.out.println();
                }
        );
    }

    @Test
    public void testChromPosition() throws Exception {
        System.out.println("LCJ - testChromPosition begin");
        Chromosome elevenC = new Chromosome("11c");
        Chromosome oneA = new Chromosome("1a");
        if (elevenC.compareTo(oneA) > 0) {
            String message = "The start chromosome " + elevenC.getName() 
                    + " is larger than the end chromosome " + oneA.getName();
            System.out.println(message);

        }  else {
            System.out.println("LCJ - 11c is less than 1a");
        }
        Chromosome tenA = new Chromosome("10a");
        if (tenA.compareTo(oneA) > 0) {
            String message = "The start chromosome " + tenA.getName() 
                    + " is larger than the end chromosome " + oneA.getName();
            System.out.println(message);

        }  else {
            System.out.println("LCJ - 10a is less than 1a");
        }
        Chromosome tenB = new Chromosome("10b");
        if (tenB.compareTo(oneA) > 0) {
            String message = "The start chromosome " + tenB.getName() 
                    + " is larger than the end chromosome " + oneA.getName();
            System.out.println(message);

        }  else {
            System.out.println("LCJ - 10b is less than 1a");
        }
    }
    // Test creating a reference
    @Test
    public void testCreateReferenceTag() throws Exception {
        //Create tag list
        System.out.println("Begin testCreateReferenceTag ...");
        Map<Tag,TaxaDistribution> tagTDMap;
        Map<Tag,TaxaDistribution> tagTDMap2;
        List<String> sequences= ImmutableList.of("ACGACGACG","ACtACGACG","ACGACGtCG","ACtACGGtCG",
                "TGTCGTCGT");
//        List<String> sequences2= ImmutableList.of("ACGACGACGTTCT","ACtACGACGTT","ACGACGtCGTT","AtACGGtCGAAGG",
//                "TGTCGTCGT");
        // These are reverse strands
        List<String> sequences2= ImmutableList.of("AGAACGTCGTCGT","AACGTCGTAGT","AACGACGTCGT","CCTTCGACCGTAT",
                "TGTCGTCGT");
        List<Boolean> direction= ImmutableList.of(true, true, true, true,false);
        List<TaxaDistribution> taxaDists=ImmutableList.of(
                TaxaDistBuilder.create(10,new int[]{0,1,2,3},new int[]{1,2,1,3}),
                TaxaDistBuilder.create(10,new int[]{4,5},new int[]{1,1}),
                TaxaDistBuilder.create(10,new int[]{3,6,7},new int[]{3,2,1}),
                TaxaDistBuilder.create(10,8),
                TaxaDistBuilder.create(10,1)
        );
        List<Tag> tags= sequences.stream().map(s -> TagBuilder.instance(s).build()).collect(Collectors.toList());
        //tags.set(0,TagBuilder.instance(tags.get(0).sequence()).reference().build());  //make first tag the reference tag
        tagTDMap=new HashMap<>();
        for (int i = 0; i < tags.size(); i++) {  //This should all be doable with a Stream zip function, but it was removed from beta versions, and I don't know how to do it now.
            tagTDMap.put(tags.get(i),taxaDists.get(i));
        }
        
        // Create second tagTaxaDistMap for testing filtered alignments. There exists a tag, not marked as
        // reference that will match the reference we create.
        List<Tag> tags2= sequences2.stream().map(s -> TagBuilder.instance(s).build()).collect(Collectors.toList());
        System.out.println("tags2 size: " + tags2.size());
        tagTDMap2=new HashMap<>();
        for (int i = 0; i < tags2.size(); i++) {  //This should all be doable with a Stream zip function, but it was removed from beta versions, and I don't know how to do it now.
            tagTDMap2.put(tags2.get(i),taxaDists.get(i));
            System.out.println("seq2 forloop: added tag, i=" + i);
        }
        System.out.println("tagTDMap2 size after creating: " + tagTDMap2.keySet().size());
        Chromosome chrom = new Chromosome("9");
        Position cutPos=new GeneralPosition
                .Builder(chrom,355420)
                .strand((byte)1) 
                .build();
        
        // Pretest - this verifies that position 355420 in the stored reference
        // genome file begins a sequence reasonable to be the cutposition site for
        // out test taxa tags.
        GenomeSequence myRefGenome = GenomeSequenceBuilder.instance(REFERENCE_GENOME_FILE);
        byte[] seqForTaxa = myRefGenome.chromosomeSequence(chrom,355420,355429);
        String seqForTaxaString = NucleotideAlignmentConstants.nucleotideBytetoString(seqForTaxa);
        System.out.println("seqForTaxaString is: " + seqForTaxaString + ", expecting ACGACGACGC");
        assertTrue(seqForTaxaString.equals("ACGACGACGC"));
        // Create a reference that contains the tags used above from tagTaxaDistMap
        // or has something close to these tags in the position below.  Verify what
        // comes back in tagTaxaDistMap is the tags with a reference added if we set
        // the reference, or the same tags if we did not.
        
        // Check current number of tags in the map.  
        int tagCount = tagTDMap.keySet().size();
        System.out.println("initial tagCount from map: " + tagCount);
        // Use the TASSEL-5-TEST stored reference file for partial chrom 9-10 as a reference.
        Map<Tag,TaxaDistribution> tagTDMapRet = usePrivCreateReferenceTag(cutPos,tagTDMap,chrom,192, REFERENCE_GENOME_FILE,true);
        
        //Verify tagTaxaDistMap has same number of entries, contains a reference tag, 
        // and ref tag sequence matches our expected reference sequence
        int newTagCount = tagTDMapRet.keySet().size();
        System.out.println("initial tag count:" + tagCount + " new tag count:" + newTagCount);
        assertEquals(tagCount+1,newTagCount);
        
        //List<Tag> tagTaxaDistList = tagTaxaDistMap.keySet().stream().collect(Collectors.toList());

        Tag refTag=tagTDMapRet.keySet().stream()
                .filter(Tag::isReference)
                .findFirst().orElseThrow(() -> new IllegalStateException("Reference not found"));
        System.out.println("LCJ - refTag sequence from tagTDMap: " + refTag.sequence() + " expecting:" + seqForTaxaString);
        assertTrue("Created reference tag doesnot match expected reference sequence",refTag.sequence().equals(seqForTaxaString));
        
        // Create new tagTaxaDistMap that has no sequence that is exact match of reference
        System.out.println("\nLCJ - testCreateReferenceTag - begin reverse strand test ...");
        Chromosome chrom10 = new Chromosome("10");
        cutPos=new GeneralPosition
                .Builder(chrom10,5291778)
                .strand((byte)1) 
                .build();
        byte[] seqForTaxa2 = myRefGenome.chromosomeSequence(chrom,5291778,5291790);
        String seqForTaxaString2 = NucleotideAlignmentConstants.nucleotideBytetoString(seqForTaxa2);
        Tag seq2Tag = TagBuilder.instance(seqForTaxaString2).build();
        String reverseSeq2Tag = seq2Tag.toReverseComplement();
        System.out.println("LCJ - reverseSeq2Tag is: " + reverseSeq2Tag + ", expecting AGAACGTCGTCGT");
        assertTrue(reverseSeq2Tag.equals("AGAACGTCGTCGT"));
        tagCount = tagTDMap2.keySet().size();
        
        // Test creating reverse strand tag
        tagTDMapRet = usePrivCreateReferenceTag(cutPos,tagTDMap2,chrom,192, REFERENCE_GENOME_FILE,false);
        // When I check tagTDMap2.keySet().size AFTER the usePrivCreateReferencTag call, it comes
        // back as 4 instead of 5.   However, the new map returned has the correct value of 5
        newTagCount = tagTDMapRet.keySet().size();
        System.out.println("old tagCount " + tagCount + " newTagCount " + newTagCount);
        // Tag count should be one greater as the reference tag was not already 
        // in the map
        assertEquals(tagCount+1,newTagCount);
                
        refTag=tagTDMapRet.keySet().stream()
                .filter(Tag::isReference)
                .findFirst().orElseThrow(() -> new IllegalStateException("Reference not found"));
        
        // This test needs fixing since we've added changes to finding the reference.
        
//        System.out.println("LCJ - refTag sequence from tagTDMap: " + refTag.sequence() + " expecting:" + seqForTaxaString2);
//        assertTrue("Created reference tag doesnot match expected reference sequence",refTag.sequence().equals(seqForTaxaString2));
//        System.out.println("refTag sequence from tagTDMap: " + refTag.sequence() + " expecting:" + reverseSeq2Tag);
//        assertTrue("Created reference tag doesnot match expected reference sequence",refTag.sequence().equals(reverseSeq2Tag));

        System.out.println("\nFinished testCreateReferenceTag");
    }   
    
    @Test
    public void testNumTagsWithReference() throws Exception {
        // run the pipeline from GBSSeqToTagPlugin through Discovery.
        // Verify at the end of Discovery the number of tags in the DB
        // has not increased, nor descreased from what GBSSeqToTagDBPlugin stored.
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

        // Get tag list
        TagData tagData =new TagDataSQLite(GBSConstants.GBS_GBS2DB_FILE);
        Set<Tag> tags = tagData.getTags();
        int numGBSTags = tags.size();
        System.out.println("JunitTest - number of tags in tag table after GBSSeq: " + tags.size());
        ((TagDataSQLite)tagData).close();
        
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
        .minMinorAlleleFreq(0.1)
        .startChromosome(new Chromosome("9"))
        .endChromosome(new Chromosome("10"))
        .referenceGenomeFile(REFERENCE_GENOME_FILE)
        .deleteOldData(true)
        //.gapAlignmentThreshold(0.12) //  default is 1.0 (keep all alignments)
        .performFunction(null);
        
        // Get tags again - verify the same number
        tagData =new TagDataSQLite(GBSConstants.GBS_GBS2DB_FILE);
        tags = tagData.getTags();
        int numDiscoveryTags = tags.size();
        System.out.println("JunitTest - number of tags in tag table after Discovery: " + tags.size());
        ((TagDataSQLite)tagData).close();
        
        assertEquals(numDiscoveryTags,numGBSTags);

    }
    private void convertAlignedTagsToTable (Map<Tag,String> alignedTags, Map<Tag,TaxaDistribution> distMap, Position startPos) throws Exception {
        Table<Position, Byte, List<TagTaxaDistribution>> tAlign=usePrivConvertAlignmentToTagTable(alignedTags,
                distMap,  startPos);
        Table<Position, Byte, List<TagTaxaDistribution>> alignT2= TreeBasedTable.create();
        tAlign.rowMap().forEach((p, mBTTD) -> {
            if(mBTTD.entrySet().size()>1) {
                mBTTD.forEach((a,ttdl) -> alignT2.put(p,a,ttdl));
            }
        } );
        alignT2.rowMap().forEach((p,ttd) -> System.out.println(p.getPosition()+" :"+ttd.entrySet().size()));
        System.out.println(tAlign.toString());
    }
    
    private Map<Tag,String> usePrivAlignTags(Map<Tag,TaxaDistribution> tags) throws Exception {
        DiscoverySNPCallerPluginV2 discoverySNPCallerPluginV2 = new DiscoverySNPCallerPluginV2();
        Class theClass = discoverySNPCallerPluginV2.getClass();
        Method theMethod = theClass.getDeclaredMethod("alignTags", new Class[] {Map.class, int.class,byte.class,boolean.class} );
        theMethod.setAccessible(true);
        // Hard coding true - need to fix this
        return (Map<Tag,String>)theMethod.invoke(discoverySNPCallerPluginV2,(Object)tags, 64,(byte)1,true);
    }
    
    private Map<Tag,String> usePrivFilterAlignTags(Map<Tag,String> tags, Position refStartPosition, double threshold) throws Exception {
        DiscoverySNPCallerPluginV2 discoverySNPCallerPluginV2 = new DiscoverySNPCallerPluginV2().gapAlignmentThreshold(threshold);
        Class theClass = discoverySNPCallerPluginV2.getClass();
        Method theMethod = theClass.getDeclaredMethod("filterAlignedTags", new Class[] {Map.class, Position.class, double.class} );
        theMethod.setAccessible(true);
        return (Map<Tag,String>)theMethod.invoke(discoverySNPCallerPluginV2,(Object)tags, refStartPosition, threshold);
    }

    private Table<Position, Byte, List<TagTaxaDistribution>>  usePrivConvertAlignmentToTagTable(Map<Tag,String> alignedTags,
                              Map<Tag,TaxaDistribution> tagTaxaDistMap, Position refStartPosition) throws Exception {
        DiscoverySNPCallerPluginV2 discoverySNPCallerPluginV2 = new DiscoverySNPCallerPluginV2();
        Class theClass = discoverySNPCallerPluginV2.getClass();
        Method theMethod = theClass.getDeclaredMethod("convertAlignmentToTagTable", new Class[] {Map.class, Map.class, Position.class} );
        theMethod.setAccessible(true);
        return (Table<Position, Byte, List<TagTaxaDistribution>>)theMethod.invoke(discoverySNPCallerPluginV2,
                (Object)alignedTags, (Object)tagTaxaDistMap,(Object)refStartPosition);
        //org.apache.commons.lang3.reflect.MethodUtils - this could be used to simplify this.
    }
    
    private Map<Tag,TaxaDistribution> usePrivCreateReferenceTag(Position cutPos, Map<Tag,
            TaxaDistribution> tags, Chromosome chrom, int numberOfTaxa, String refFile, boolean direction) throws Exception {
        DiscoverySNPCallerPluginV2 discoverySNPCallerPluginV2 = new DiscoverySNPCallerPluginV2();
        Class theClass = discoverySNPCallerPluginV2.getClass();
        // Set reference
        System.out.println("usePrivCreateReferenceTag, refFile is: " + refFile);

        discoverySNPCallerPluginV2.includeReference = true;
        discoverySNPCallerPluginV2.myRefSequence = GenomeSequenceBuilder.instance(refFile);

        Method theMethod = theClass.getDeclaredMethod("createReferenceTag", new Class[] {Position.class,Map.class,
                Chromosome.class,int.class,boolean.class});
        theMethod.setAccessible(true);
        return (Map<Tag,TaxaDistribution>)theMethod.invoke(discoverySNPCallerPluginV2,(Object)cutPos,(Object)tags, (Object)chrom,numberOfTaxa,direction);


    }

}