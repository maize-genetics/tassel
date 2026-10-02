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
	    // Build the whole pipeline (GBSSeqToTagDB -> TagExport -> generated SAM -> Discovery) from
	    // simulated data and assert on invariants: tag count does not change through Discovery and a
	    // nonzero (in fact the injected) number of SNPs is called.
	    GBSSimData sim = GBSSimData.createUnder("Discovery_full");
	    sim.buildDatabaseThroughSam();
	    int tagsBefore = sim.dbTagCount();

	    sim.runDiscovery(false); // no reference genome -> setCommonToReference path

	    int tagsAfter = sim.dbTagCount();
	    assertEquals("Discovery must not change the number of tags", tagsBefore, tagsAfter);

	    TagData tagData = new TagDataSQLite(sim.dbFile.toString());
	    int snpCount = tagData.getSNPPositions().size();
	    ((TagDataSQLite) tagData).close();
	    System.out.println("testFullSNPCaller: SNP positions called = " + snpCount);
	    assertTrue("Discovery should call at least one SNP", snpCount > 0);
	    assertEquals("Should call one SNP per injected SNP locus", sim.numSnpLoci, snpCount);
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
    // Test creating a reference tag from the simulated reference genome.
    @Test
    public void testCreateReferenceTag() throws Exception {
        System.out.println("Begin testCreateReferenceTag ...");
        GBSSimData sim = GBSSimData.createUnder("Discovery_createRef");

        // Choose a known reference-allele locus from the simulated data.
        GBSSimData.TagInfo refInfo = sim.tagInfos.stream()
                .filter(t -> !t.variant)
                .findFirst().orElseThrow(() -> new IllegalStateException("No reference-allele tag"));
        Chromosome chrom = new Chromosome(refInfo.chrom);
        Position cutPos = new GeneralPosition.Builder(chrom, refInfo.cutPosition).strand((byte) 1).build();

        // The expected reference tag is exactly the sequence stored at the cut position.
        GenomeSequence refGenome = GenomeSequenceBuilder.instance(sim.referenceFasta.toString());
        byte[] refBytes = refGenome.chromosomeSequence(chrom, refInfo.cutPosition,
                refInfo.cutPosition + GBSSimData.TAG_LENGTH - 1);
        String expectedRefSeq = NucleotideAlignmentConstants.nucleotideBytetoString(refBytes);
        assertEquals("Simulated reference should contain the planted reference-allele tag",
                refInfo.sequence, expectedRefSeq);

        // Build a tag map that does NOT already contain the reference-allele tag, using a full-length
        // tag from another locus so the longest-tag length is TAG_LENGTH.
        GBSSimData.TagInfo otherInfo = sim.tagInfos.stream()
                .filter(t -> !t.sequence.equals(refInfo.sequence))
                .findFirst().orElseThrow(() -> new IllegalStateException("No other tag"));
        Map<Tag, TaxaDistribution> tagTDMap = new HashMap<>();
        tagTDMap.put(TagBuilder.instance(otherInfo.sequence).build(),
                TaxaDistBuilder.create(sim.taxa.size(), new int[]{0, 1}, new int[]{2, 2}));
        int tagCount = tagTDMap.keySet().size();

        Map<Tag, TaxaDistribution> tagTDMapRet = usePrivCreateReferenceTag(cutPos, tagTDMap, chrom,
                sim.taxa.size(), sim.referenceFasta.toString(), true);

        // The reference tag was not already present, so the map grows by one.
        assertEquals(tagCount + 1, tagTDMapRet.keySet().size());

        Tag refTag = tagTDMapRet.keySet().stream()
                .filter(Tag::isReference)
                .findFirst().orElseThrow(() -> new IllegalStateException("Reference not found"));
        assertEquals("Created reference tag must match the simulated reference sequence",
                expectedRefSeq, refTag.sequence());

        System.out.println("\nFinished testCreateReferenceTag");
    }
    
    @Test
    public void testNumTagsWithReference() throws Exception {
        // Run the pipeline from GBSSeqToTagDBPlugin through Discovery WITH a reference genome.
        // Verify Discovery neither adds nor removes tags, and that SNPs are still called.
        GBSSimData sim = GBSSimData.createUnder("Discovery_numTagsRef");
        sim.buildDatabaseThroughSam();

        int numGBSTags = sim.dbTagCount();
        System.out.println("JunitTest - number of tags in tag table after GBSSeq+SAM: " + numGBSTags);

        sim.runDiscovery(true); // supply the simulated reference genome

        int numDiscoveryTags = sim.dbTagCount();
        System.out.println("JunitTest - number of tags in tag table after Discovery: " + numDiscoveryTags);
        assertEquals(numGBSTags, numDiscoveryTags);

        TagData tagData = new TagDataSQLite(sim.dbFile.toString());
        int snpCount = tagData.getSNPPositions().size();
        ((TagDataSQLite) tagData).close();
        assertTrue("Discovery with reference should still call SNPs", snpCount > 0);
        assertEquals("Should call one SNP per injected SNP locus", sim.numSnpLoci, snpCount);
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