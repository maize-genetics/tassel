/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import javax.swing.ImageIcon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;

import net.maizegenetics.analysis.gbs.neobio.BasicScoringScheme;
import net.maizegenetics.analysis.gbs.neobio.IncompatibleScoringSchemeException;
import net.maizegenetics.analysis.gbs.neobio.InvalidSequenceException;
import net.maizegenetics.analysis.gbs.neobio.PairwiseAlignment;
import net.maizegenetics.analysis.gbs.neobio.PairwiseAlignmentAlgorithm;
import net.maizegenetics.analysis.gbs.neobio.ScoringScheme;
import net.maizegenetics.analysis.gbs.neobio.SmithWaterman;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.GenomeSequence;
import net.maizegenetics.dna.map.GenomeSequenceBuilder;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.Utils;

/**
 * This plugin takes an existing repGen db, grabs the tags
 * whose depth meets that specified in the minCount parameter,
 * makes kmer seeds from these tags.  Window for kmer seeds is default 20.
 * 
 * The ref genome is walked with a sliding window of 1. Reference tags are created
 * based on peaks within clusters where kmer seeds align.  
 * 
 *  The kmerLen field should match the length of the kmers stored
 *  as tags during the RepGenLoadSeqToDBPlugin step.  The default is 150.
 *  
 *  The refKmerLen() should minimally be the length of the db kmer
 *  tags, but can be longer.  Our defaults are 150 for kmer tags, and
 *  twice this length (300) for the refKmerLen.  
 *  
 *  There are 2 count parameters:  minTagCount specifies the minimum depth
 *  of a tag for it to be used when creating seed kmers.
 *  
 *  minHitCount specifies the number of "hits" that must occur within a sliding 
 *  window (window size = refKmerLen) for the window to remain part of a cluser.
 *  When sliding, if the hit count drops below the minhitCount threshold, the
 *  cluster ends.  A new cluster does not begin until the hit count is back up
 *  to threshold.  See createRefTagsForAlignment() and storeRefTagPositions()
 *  for specifics.
 *  
 *  This plugin creates and stores the reference tags in the refTag table
 *  in the database.   Both the tagMapping and the physicalMapPosition table will
 *  we populated with the reference tag information.
 *  
 *  Once the tables have been populated with the reference information,
 *  Smith Waterman is run to align all the nonreference tags in the db 
 *  against each other; each non-reference tag against the reference tags;
 *  finally each refTag against all other refTags.
 *  
 * ALignment data is stored in the tagAlignments table.
 * 
 * Smith Waterman from SourceForge neobio project is used
 * to determine alignment score.  Settings for match rewards, mismatch penalty
 * and gap penalty may be changed by user via plugin parameters.
 * 
 * 
 * @author lcj34
 *
 */
public class RepGenAlignerPlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(RepGenAlignerPlugin.class);
    
    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with tags and taxa distribution").build();
    private PluginParameter<String> refGenome = new PluginParameter.Builder<String>("ref", null, String.class).guiName("Reference Genome File").required(true)
            .description("Referemce Genome File for aligning against ").build();
    private PluginParameter<Integer> minTagCount = new PluginParameter.Builder<Integer>("minTagCount", 1, Integer.class).guiName("Min Tag Count")
            .description("Minimum count of reads for a tag to be aligned").build();
    private PluginParameter<Integer> minHitCount = new PluginParameter.Builder<Integer>("minHitCount", 4, Integer.class).guiName("Minimum Hit Count")
            .description("Minimum number of hits in a cluster for a reference tag to be created from that location.").build();
    private PluginParameter<Integer> seedLen = new PluginParameter.Builder<Integer>("seedLen", 17, Integer.class).guiName("Seed Kmer Length")
            .description("Length of kmer seed created from DB tags and used as seeds for aligning against the reference genome.").build();
    private PluginParameter<Integer> seedWindow = new PluginParameter.Builder<Integer>("seedWindow", 20, Integer.class).guiName("Window Length for Seed Creation")
            .description("Length of window between positions when creating seed from DB tags.").build();
    private PluginParameter<Integer> kmerLen = new PluginParameter.Builder<Integer>("kmerLen", 150, Integer.class).guiName("Kmer Length")
            .description("Length of kmer from fastq reads stored as the tag sequence in the DB.").build();
    private PluginParameter<Integer> refKmerLen = new PluginParameter.Builder<Integer>("refKmerLen", 300, Integer.class).guiName("Reference Kmer Length")
            .description("Length of kmers created from reference genome to store in the DB. \nThis should be at as long or longer than the kmerLen parameter used for storing input sequence tags.").build();
    private PluginParameter<Integer> match_reward = new PluginParameter.Builder<Integer>("match_reward", 2, Integer.class).guiName("Match Reward Amount")
            .description("Parameter sent to Smith Waterman aligner for use in calculating reward when base pairs match.").build();
    private PluginParameter<Integer> mismatch_penalty = new PluginParameter.Builder<Integer>("mismatch_penalty", -1, Integer.class).guiName("Mismatch Penalty Amount")
            .description("Parameter sent to Smith Waterman aligner for use in calculating penalty when base pairs are mis-matched.").build();
    private PluginParameter<Integer> gap_penalty = new PluginParameter.Builder<Integer>("gap_penalty", -1, Integer.class).guiName("Gap Penalty Amount")
            .description("Parameter sent to Smith Waterman aligner for use in calculating penalty when when a gap is identified.").build();
    
    static GenomeSequence myRefSequence = null;

    public RepGenAlignerPlugin() {
        super(null, false);
    }

    public RepGenAlignerPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public RepGenAlignerPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    @Override
    public void postProcessParameters() {

        if (myDBFile.isEmpty() || !Files.exists(Paths.get(inputDB()))) {
            throw new IllegalArgumentException("RepGenAlignerPlugin: postProcessParameters: Input DB not set or found");
        }
        if (!refGenome.isEmpty()) {
            myRefSequence = GenomeSequenceBuilder.instance(refGenome());
        } else {
            throw new IllegalArgumentException("RepGenAlignerPlugin: postProcessParameters: reference genome not set or found");
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        long totalTime = System.nanoTime();
        long time=System.nanoTime();
 
        try {           
            System.out.println("RepGenAlignerPlugin:processData begin"); 
            RepGenDataWriter repGenData=new RepGenSQLite(inputDB());
            Map<Tag, Integer> tagsWithDepth = repGenData.getTagsWithDepth(minTagCount());
            if (tagsWithDepth.isEmpty()) {
                System.out.println("\nNo tags found with minimum depth " + minTagCount() + ". Halting Run.\n");
                ((RepGenSQLite)repGenData).close();
                return null;
            }
            
//            BufferedWriter tagbw = Utils.getBufferedWriter("/home/lcj34/repGen_outputFiles/tagsLoaded.txt");
//            
//            try {
//                List<String> tags = new ArrayList<String>();
//                for (Tag tag: tagsWithDepth.keySet()){
//                    tags.add(tag.sequence());
//                }
//                Collections.sort(tags);
//                for (String seq : tags) {
//                    tagbw.write(seq);
//                    tagbw.write("\n");
//                }
//                tagbw.close();
//            } catch (IOException ioe) {
//                ioe.printStackTrace();
//            }
                    
            Multimap<String,Tag> kmerTagMap = HashMultimap.create();
                          
            System.out.println("Calling createKmerSeedsFromDBTags with window size: " + seedWindow());
           // Create map of kmer seeds from db tags
            createKmerSeedsFromDBTags(tagsWithDepth.keySet(), kmerTagMap,  seedWindow());
            System.out.println("Num distinct kmerSeeds created: " + kmerTagMap.keySet().size() + 
                    ", kmerTagMap size is:" + kmerTagMap.size() + ",TotalTime for createKmerSeedsFromDBTags was " + (System.nanoTime() - time) / 1e9 + " seconds");
 
            System.out.println("Size of tagsWithDepth: " + tagsWithDepth.size());
            System.out.println("Size of kmerTagMap keyset: " + kmerTagMap.keySet().size());
            
            time = System.nanoTime();
            // myRefSequence populated in postProcessParameters 
            Set<Chromosome> chromsInRef = myRefSequence.chromosomes();
            
            // Create hash of bitmaps for the chromosomes:
            Map<String, OpenBitSet> chromBitMaps = new ConcurrentHashMap<String,OpenBitSet>();
            
            System.out.println("Start making array of chrom bitmaps ");
            // For each chrom in refSequence, walk the refSequence looking for kmers
            // matching those in the kmerTagMap.  Store hit positions in per-chromosome
            // bitmaps to be used for creating ref tags for aligning.
            
            chromsInRef.parallelStream().forEach(chrom -> { 
                  // Turn this on/off for debug purposes
                  //if (chrom.getChromosomeNumber() != 9) return; // just for initial testing !!! - remove
                                    
                  int kmersForChrom = 0;
                  int chromLength = myRefSequence.chromosomeSize(chrom);
                  System.out.println("\nChecking reference chrom " + chrom.getName() + ", size: " + chromLength + " for tag kmer matches.");
                  
                  // create bitmap for this chrom
                  OpenBitSet chromBits = new OpenBitSet(chromLength);
                  
                  for (int chromIdx = 0; chromIdx < chromLength;) {                                       
                      // chromosomeSequence:  start and end are inclusive and 1-based, adjust for this in call below
                      int end = Math.min(chromLength, chromIdx+seedLen());
                      byte[] chromSeedBytes = myRefSequence.chromosomeSequence(chrom,chromIdx+1, end);
                      // we're not processing from fastq, we're processing from tag table in db
                      int badValue = checkForN(chromSeedBytes);
                      if (badValue >=0) {
                          // found a non-ACGT character,  re-set index and
                          // grab new kmer from the chromosome
                          chromIdx += badValue+1; // returned value was at bad base - move beyond it
                          continue;
                      }
                                         
                      // Convert back to allele string
                      String chromKmerString = NucleotideAlignmentConstants.nucleotideBytetoString(chromSeedBytes);
                      // chromKmerBytes were good - check if this kmer exists among our seeds
                      if (kmerTagMap.containsKey(chromKmerString)){
                          kmersForChrom++;                         
                          chromBits.fastSet(chromIdx);  // set position in the bitmap
                      }                                                 
                      chromIdx++; // after processing, slide up by 1  
                  }

                  chromBitMaps.put(chrom.getName(), chromBits); 
                  System.out.println("Total tag seeds matching to kmers in chrom " + chrom.getName() + ": " 
                      + kmersForChrom + ", total fastBits set via cardinality: " + chromBits.cardinality());             
              });
                         
            // The bitmaps of positions matching a kmer seed start have been calculated for each chrom.
            // Take these maps and find where hits are clustered.  Find the maxima in each cluster,
            // use the maxima to determine a start position for the reference tag to be created.
            
            // Hashmap holds start positions, on each chromosome, where clusters
            // of hits occur, ie kmers map to these regions.  THe map is keyed by chromosome name,
            // map values are the position at the center of hit clusters.  SHould this map include the tag?  
            
            System.out.println("Entering loop to find clusters" );
            Multimap<String,Integer> chromMaximaMap = HashMultimap.create(); // holds chrom/positions for creating refTag
            Multimap<Tag,Position> refTagPositionMap =  HashMultimap.create(); // holds ref tags and positions

            for (String chrom : chromBitMaps.keySet()) {
                OpenBitSet chromHits = chromBitMaps.get(chrom);
                createRefTagsForAlignment(chromHits, chrom, chromMaximaMap, refTagPositionMap);
            }
            
            // LCJ - below is debug!  Writes to my cbsu dir
            String refTagOutFile = "/home/lcj34/repGen_outputFiles/allChroms_refTagsOrigAligner.txt";
            // LCJ - this one writes to laptop
            //String refTagOutFile = "/Users/lcj34/notes_files/repgen/junit_phase2Aligner/chrom" + chrom.getName() + "_refTagsFromPrimers.txt";
            BufferedWriter refwr = Utils.getBufferedWriter(refTagOutFile);                  
            try {
                for (Map.Entry<Tag, Position> entry : refTagPositionMap.entries()) {
                    int position = entry.getValue().getPosition();
                    String chromname = entry.getValue().getChromosome().getName();
                    Tag refTag = entry.getKey();
                    String refData = chromname + "\t" + position + "\t" + refTag.sequence() + "\n";
                    refwr.write(refData);
                }
                refwr.close();
            } catch (IOException e1) {
                // TODO Auto-generated catch block
                e1.printStackTrace();
            }
            System.out.println("\nNumber of distinct refTags to be loaded into db: " 
              + refTagPositionMap.keySet().size() + " total refTags: " + refTagPositionMap.size());
            // add ref and mapping approach to db, then reference tags
            repGenData.addMappingApproach("SmithWaterman");
            repGenData.addReferenceGenome(refGenome());
            
            // Add ref tag and tags separately. Add to tagAlignment map,
            repGenData.putRefTagMapping(refTagPositionMap, refGenome());
            
            // Get tags stored from RepGenLoadSeqToDB        
            Set<Tag> tagsToAlign = repGenData.getTags();
            List<Tag> tagList = new ArrayList<Tag>(tagsToAlign);
            
            // Create synchronized map for use in parallel streams 
            Multimap<Tag,AlignmentInfo> tagAlignInfoMap = Multimaps.synchronizedMultimap(HashMultimap.<Tag,AlignmentInfo>create());
            Multimap<RefTagData,AlignmentInfo> refTagAlignInfoMap = Multimaps.synchronizedMultimap(HashMultimap.<RefTagData,AlignmentInfo>create());
            
            // First align the tags against each other
            // Output of this is stored in db table tagAlignments
            System.out.println("Calling calculateTagTagAlignment for tags without refs");
            
            calculateTagTagAlignment( tagList, tagAlignInfoMap);
            System.out.println("Number of tag-tag alignments: " + tagAlignInfoMap.size() + ", store to db.");
            // neither tag is reference, null and -1 for tag1 chrom/pos
            //repGenData.putTagTagAlignments(tagAlignInfoMap, false,false, null,-1,refGenome());
            repGenData.putTagTagAlignments(tagAlignInfoMap);
 
            System.out.println("Calling calculateTagRefTagALignment for tags WITH ref");
            Set<RefTagData> refTags = repGenData.getRefTags();
            List<RefTagData> refTagList = new ArrayList<RefTagData>(refTags);
            
            tagAlignInfoMap.clear(); // remove old data
            calculateTagRefTagAlignment(tagList,refTagList,tagAlignInfoMap,refGenome());
            
            // Add the tag-refTag info to the tagAlignments table.           
            // tag1=nonref, tag2=ref, null and -1 for tag1 .  Alignment will be
            // done twice:  once for tag/refTag-fwd, once for tag/refTag-reverse
            System.out.println("Number of tag-refTag alignments, includes aligning to ref fwd and reverse strands: " + tagAlignInfoMap.size() + ", store to db.");
            //repGenData.putTagAlignments(tagAlignInfoMap,false,true,null,-1,refGenome());
            repGenData.putTagRefTagAlignments(tagAlignInfoMap,refGenome());
            
            System.out.println("Calling calculateRefRefAlignment");            
            calculateRefRefAlignment(refTagList,refTagPositionMap,refTagAlignInfoMap);
            System.out.println("Number of reftag-reftag alignments: " + refTagAlignInfoMap.size() + ", store to db.");
            repGenData.putRefRefAlignments(refTagAlignInfoMap, refGenome()); // CREATE putRefRefAlignments -  different parameters
            ((RepGenSQLite)repGenData).close();

            myLogger.info("Finished RepGenAlignerPlugin\n");
        } catch (Exception e) {
            myLogger.info("Catch in prcessing REpGEnALignerPlugin file e=" + e);
            e.printStackTrace();
        }
        System.out.println("Process took " + (System.nanoTime() - totalTime)/1e9 + " seconds.\n");
        return null;
    }
    
    private short calculateFirstRangeCount(int rangeSize, OpenBitSet obs) {
        short totalhits = 0;
        for (int idx = 0; idx < rangeSize; idx++) {
            if (obs.get(idx)) totalhits++;
        }
        return totalhits;
    }
        
    // THis method determines the peak maxima from the list of peak number/peak positions
    //  Input is a peak number, and all the positions that fall in this peak.
    //  From that list, we find the maxima,  create the reference tag and its position information
    // This is stored in a data structure that will be used to populate the DB.  We store a position
    // so the physical site, and any annotations can be used.
    private void storeRefTagPositions(int peaknum,List<Integer> positionsInPeak,Multimap<String,Integer> chromMaximaMap, 
            String chrom, long chromSize, Multimap<Tag, Position> refTagPositionMap) {
        // positionsInPeak is a sorted list.
        int size = positionsInPeak.size();
        
        if (size < minHitCount()) {
        //if (size < 400) { // LCJ - testing to get fewer refTags.  hit count still valid for movingSums
            System.out.println("storeRefTagPosition: dropping peak number " + peaknum + ", only " + size + " positions in cluster.");
            return;
        }
        int start = positionsInPeak.get(0);
        int end = positionsInPeak.get(size-1);
        int maxima = start + (end-start)/2 + seedLen()/2; // add half the kmer seed len, we want to store the midpoint of the kmer

        chromMaximaMap.put(chrom, maxima);
        // Create the reference tag, add to tag/position list for adding to db.
        int lowpos = maxima - (refKmerLen()/2);
        int startRefPos = lowpos >0 ? lowpos : 1;
        int endpos = startRefPos == 1 ? refKmerLen() : startRefPos + (refKmerLen()-1); // -1 so we aren't 1 over the specified length
        if (endpos > chromSize-1) {
            endpos = (int)(chromSize-1);
        }
        Chromosome chromObject = new Chromosome(chrom);
        byte[] refTagBytes = myRefSequence.chromosomeSequence(chromObject,startRefPos, endpos);

        // Convert to allele string
        String refString = NucleotideAlignmentConstants.nucleotideBytetoString(refTagBytes);
        if (refString.contains("N") || refString.contains("null")) {
            System.out.println("StoreRefTagPositions... refString contains N or null: " + refString); 
            // If the ref string has N's in it, cut it back to the same size as 
            // the kmer and re-grab it.  if it still has N's, discard
            lowpos = maxima - (kmerLen()/2);
            startRefPos = lowpos >0 ? lowpos : 1;
            endpos = startRefPos == 1 ? kmerLen() : startRefPos + (kmerLen()-1); // -1 so we aren't 1 over specified length
            refTagBytes = myRefSequence.chromosomeSequence(chromObject,startRefPos, endpos);
            refString = NucleotideAlignmentConstants.nucleotideBytetoString(refTagBytes);
            if (refString.contains("N") || refString.contains("null")) {
                System.out.println(" - after adjusting, refString still contains N, drop this one");
            } else {
                System.out.println(" - new refTag after adjusting, len=" + refString.length() + ", " + refString);
            }
        }
        // RefString could be null for some other reason - e.g. start/end pos is bad
        if (refString == null || refString.length() == 0) {
            System.out.println(" storeRefTagPositions - refString is NULL");
        } else {
            Tag refTag = TagBuilder.instance(refString).reference().build();
            if (refTag != null ) {
                Position refPos=new GeneralPosition
                        .Builder(chromObject,startRefPos)
                        .strand((byte)1)
                        .addAnno("mappingapproach", "SmithWaterman") // this wasn't done with SW !!  Is just a reference
                        .addAnno("forward", "true") // reference, so always forward strand.
                        .build();
                // This is a multimap because it is possible for a tag sequence to
                // show up in multiple places in the genome.  If multiple places with
                // this tag are at maxima sites, then the multimap has all positions
                // for this tag.
                refTagPositionMap.put(refTag, refPos);
            } else {
                System.out.println("- refTag is NULL for refString: " + refString);
            } 
        } 
        // DEBUG - COMMENT THIS OUT !!!
        //writePeakPositions(chrom,peaknum,positionsInPeak);
    }
    
    private void createKmerSeedsFromDBTags(Set<Tag> tagWithDepth,
            Multimap<String,Tag> kmerTagMap, int window) {
        for (Tag tag : tagWithDepth) {
            String tagSequence = tag.sequence();
            int maxIdx = tagSequence.length() - window;
            for (int seqIdx = 0; seqIdx < maxIdx;) {
                String kmer = tagSequence.substring(seqIdx, seqIdx + seedLen());  
                byte[] kmerBytesAsNum = NucleotideAlignmentConstants.convertHaplotypeStringToAlleleByteArray(kmer);
                
                int badValues = checkForN(kmerBytesAsNum);
                if (badValues >=0) {
                    // found a non-ACGT character,  re-set index and
                    // grab new kmer
                    //seqIdx += badValues+1;
                    // just increase by window size, don't care where the bad value is
                    seqIdx += window;
                    continue;
                }
                // good values - add unique kmer to map
                kmerTagMap.put(kmer, tag);
                // add reverse complement of kmer
                // THis isn't right ... this tag does not have this KMER.
                // And we need to know it is a reverse complement.  Ed's code
                // creates rc from the ref genome, and stores the position as negative
                // We are creating seeds from tags, not the ref genomes.  how to include RC of kmer?
                byte[] kmerRC = NucleotideAlignmentConstants.reverseComplementAlleleByteArray(kmer.getBytes());
                String kmerRCString = new String(kmerRC);
                kmerTagMap.put(kmerRCString, tag); 
                // loop end - increment seqIdx for next kmer
                seqIdx += window;
            }               
        }
    }

    private void createRefTagsForAlignment(OpenBitSet chromHits, String chrom,
            Multimap<String,Integer> chromMaximaMap, Multimap<Tag,Position> refTagPositionMap) {
        // This doesn't give same value as myRefSequence.chromosomeSize(chrom)
        // FOr chrom 9, myRefSequence.chromosomeSize(9) says 159769782
        // but this code below says size is 159769792, which is 10 more. WHY?
        // oh - it returns the capacity of the set, which might be more than
        // the size I set it to if we're storing bits in words.
        
       // String peakOutFile = "/Users/lcj34/notes_files/repgen/junit2_out/chrom" + chrom + "_peaklist400.txt";
        String peakOutFile = "/home/lcj34/repGen_outputFiles/chrom" + chrom + "_peaklist400.txt";
        BufferedWriter peakwr = Utils.getBufferedWriter(peakOutFile);
        String header = "peakStart,numHits\n"; // csv for R input
        try {
            peakwr.write(header);
        } catch (IOException e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        }
        long chromSize = chromHits.size(); 
        
        // This holds the moving sums for specific positions.  The key is
        // the sum amount.  Only sums that are equal or greater than the minCount()
        // parameter are stored. 
        Multimap<Short,Integer> movingSumMap = HashMultimap.create();

        System.out.println("Size of OpenBitSet for chrom " + chrom + ": " 
        + chromSize + ", cardiality is: " + chromHits.cardinality());

        // loop through the chromosome. 
        int rangeSize = refKmerLen();               
        // THis is the position where we store the first kmer
        int firstRange = rangeSize/2; // refKmerLen should be even number, but is ok if it isn't                
        short firstRangeCount = calculateFirstRangeCount(rangeSize,chromHits);
        System.out.println("firstRangeCount: " + firstRangeCount);
        
        short currentTotal = firstRangeCount;
        int lastRange = (int)chromSize - firstRange ;
        // calculate the rest of the values
        short peak_num = 0;
        boolean belowMin = true;
        for (int idx = firstRange+1; idx < lastRange; idx++) {
            // we move up one.  If the bitMap was set for the position
            // we drop, then drop one from the count.  If the bitMap
            // is set for the position we add, add 1 to the count.
            // if our range is 10, the first range is 0-9
            // We want to store the value in the middle, which is range/2 = 5
            // That was done above.  Here, we start writing to position 6.
            // We want the totals now in range 1-10.  We dropped position 0,
            // which is 6 - (10/2) -1 = 6-5-1 = 0.
            // We want to add position 10, which is 6 + (10/2) -1,
            // which is 6 + 5 - 1 = 10
            int posToDrop = idx - (rangeSize/2) - 1;
            int posToAdd = idx + (rangeSize/2) -1;
            if (chromHits.fastGet(posToDrop)) currentTotal--;
            if (chromHits.fastGet(posToAdd)) currentTotal++;

            // Only store positions with hit counts where the hit count meets our minimum
            // May make sense to create the tag here.  Would mean we don't have to
            // traverse this map later to find positions to create the tags                                       
            if (currentTotal >= minHitCount()) {
                if (belowMin) {
                    // We dropped out of minimum, but
                    // now the counts are back up. So start a new peak
                    peak_num++;
                    belowMin = false;
                }
                // This stores the middle position of a range
                // of positions where the total kmerhit count is >= minCount()
                movingSumMap.put(peak_num,idx);
            } else {
                belowMin = true;
            }
        }   
        
        // Algorithm to find maximas (per Peter):
        //  Only store the positions that have at least the minHitCount number of hits.
        //  Whenever the hit count drops below the min, we stop adding positions for this
        //  peak
        //  When the hit count get back up to min count, we start a new peak and
        //  add positions again until the hit count at a position drops below minCount
        //  Store the peak number and its positions in the movingSumMap.
        //  
        
        Set<Short> peaks = movingSumMap.keySet();
        System.out.println("Total number of peaks calculated for chrom: " + chrom + ": "+ peaks.size());
        
        List<Short> sumList = new ArrayList<Short>(peaks); // does this need to be sorted?
        
        // Create reference tag for each peak on this chromosome
        // NOTE: chromMaximaMap is used for nothing other than debug at the end of this for loop.
        // It could be removed.

        sumList.stream().forEach(peak-> {           
            Collection<Integer> values = movingSumMap.get(peak);
            List<Integer> positionsInPeak = new ArrayList<Integer>(values);
            Collections.sort(positionsInPeak);
            // LCJ - debug
//            if ( positionsInPeak.size() >= 400) {
//                String peakData = positionsInPeak.size() + "\tpeak " + peak + "\tstart idx: " + positionsInPeak.get(0) + "\tend idx: " 
//                        + positionsInPeak.get(positionsInPeak.size()-1) + "\n";
//                   try {
//                       peakwr.write(peakData);
//                   } catch (Exception e) {
//                       // TODO Auto-generated catch block
//                       e.printStackTrace();
//                   }
//            }
            String peakData =  positionsInPeak.get(0) + "," + positionsInPeak.size() + "\n";

            try {
                peakwr.write(peakData);
            } catch (Exception e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            // LCJ - end debug
 
            storeRefTagPositions( peak,positionsInPeak, chromMaximaMap, chrom,chromSize,refTagPositionMap);
        });
               
        try {
            peakwr.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        // write out to a file for debug.  The cluster maxima positions are printed
        // for each cluster.  These used to create ref tags.  File output shows  the 
        // number of reference tags created for each chromosome, as well as the spacing 
        // between them. 
        //writeToFile( chrom, chromMaximaMap, minHitCount()); // writes to Lynn's directory 
    }
    
    // this calculates tag against tag alignment from the tags in the DB
    // from the tag table (not the refTag table)
    private void calculateTagTagAlignment(List<Tag> tags, Multimap<Tag,AlignmentInfo> tagAlignInfoMap){
        long totalTime = System.nanoTime();
        // For each tag on the tags list, run SW against it and store in tagTagAlignMap 
        tags.parallelStream().forEach(tag1 -> {
            for (Tag tag2: tags) {
                if (tag1.equals(tag2)) continue; // don't aligne against yourself
                String seq1 = tag1.sequence();
                String seq2 = tag2.sequence();
                Reader reader1 = new StringReader(seq1);
                Reader reader2 = new StringReader(seq2);

                PairwiseAlignmentAlgorithm  algorithm;
                ScoringScheme  scoring;
                algorithm = new SmithWaterman();
                scoring = new BasicScoringScheme(match_reward(),mismatch_penalty(),gap_penalty());
                algorithm.setScoringScheme(scoring);
                int score = 0;
 
                try {
                    algorithm.loadSequences(reader1, reader2);
                    // for tag-tag alignment, we are only computing the score
                    score = algorithm.getScore();
                } catch (IOException ioe) {
                    ioe.printStackTrace();
                } catch (InvalidSequenceException ise) { 
                    ise.printStackTrace();
                } catch (IncompatibleScoringSchemeException isse) {
                    isse.printStackTrace();
                }
                // for tag/tag, we have no chrom or position or strand or alignment position.  Store "null" and -1
                AlignmentInfo tagAI = new AlignmentInfo(tag2, null, -1, -1, -1, refGenome(),score);
                tagAlignInfoMap.put(tag1,tagAI);
            }
        });
        System.out.println("Number of tags: " + tags.size() + ", TotalTime for calculateTagTagAlignment was " + (System.nanoTime() - totalTime) / 1e9 + " seconds");
    }
    
    // This calculates alignment of each tag against each reference tag.
    private void calculateTagRefTagAlignment(List<Tag> tags, List<RefTagData> refTagDataList,
            Multimap<Tag,AlignmentInfo> tagAlignInfoMap, String refGenome){
        long totalTime = System.nanoTime();
        // For each tag on the tags list, run SW against it and store in tagAlignInfoMap  
        tags.parallelStream().forEach(tag1 -> {          
            for (RefTagData rtd : refTagDataList) {
                // Create alignment against both refTag and reverse complement of refTag
                Tag tag2 = rtd.tag();
                
                String seq1 = tag1.sequence();
                String seq2 = tag2.sequence();
                Reader reader1 = new StringReader(seq1);
                Reader reader2 = new StringReader(seq2);

                PairwiseAlignmentAlgorithm  algorithm;
                ScoringScheme  scoring;
                PairwiseAlignment alignment;

                algorithm = new SmithWaterman();
                scoring = new BasicScoringScheme(match_reward(),mismatch_penalty(),gap_penalty());
                algorithm.setScoringScheme(scoring);
                int score;
                int refAlignStartPos = rtd.position();
                int tagAlignOffset = 0; // ajust incase SW sligns from somewhere in the middle of the tag
                try {
                    algorithm.loadSequences(reader1, reader2);
                    score = algorithm.getScore(); // get score
                    alignment = algorithm.getPairwiseAlignment(); // compute alignment
                  
                    // The first sequence given in loadSequences() is loaded into rows. THis is non-refTag
                    // The second is loaded into columns (matrix[seq1][seq2] - this is the refTag
                    tagAlignOffset = alignment.getRowStart();
                    refAlignStartPos += alignment.getColStart();
                    
                    if (tagAlignOffset > 0) {
                        // Tag was not aligned from the beginning,
                        // add back the bps that were skipped so alignment begins at start of the tag
                        refAlignStartPos -= tagAlignOffset;
                    }
                    // If clipping has dropped us below the start of the reference genome, skip it
                    if (refAlignStartPos >= 0) {
                        // The ref tag start position is needed in RepGenSQLite to create a
                        // RefTagData object.  This is stored in the BiMap and used along with chrom to distinguish
                        // one tag from another.  The actual alignment position is also needed (refAlignStartPos)
                        // for the tagAlignments table.
                        AlignmentInfo tagAI = new AlignmentInfo(tag2,rtd.chromosome(),rtd.position(),refAlignStartPos, 1,refGenome,score);
                        tagAlignInfoMap.put(tag1, tagAI); // data to be stored into tagAlignments table
                    }
                                       
                    // Now align against the reverse complement of the refTag
                    reader1 = new StringReader(seq1); // readers must be reset
                    seq2 = tag2.toReverseComplement();
                    reader2 = new StringReader(seq2);
                    algorithm.unloadSequences(); 
                    algorithm.loadSequences(reader1, reader2);
                    score = algorithm.getScore(); // get score
                    alignment = algorithm.getPairwiseAlignment(); // compute alignment                 
                    tagAlignOffset = alignment.getRowStart();
                    refAlignStartPos += alignment.getColStart();                    
                            
                    if (tagAlignOffset > 0) {
                        // Tag1 was not aligned from the beginning,
                        // add back the bps that were skipped so alignment begins at start of the tag
                        refAlignStartPos -= tagAlignOffset;
                    }
                    // If clipping has dropped us below the start of the reference genome, skip it
                    if (refAlignStartPos >= 0) {
                        AlignmentInfo tagAI = new AlignmentInfo(tag2,rtd.chromosome(),rtd.position(),refAlignStartPos, 0, refGenome(),score);
                        tagAlignInfoMap.put(tag1, tagAI); // data to be stored into tagAlignments table
                    }
                                       
                } catch (IOException e) {
                    e.printStackTrace();
                } catch (InvalidSequenceException e) {
                    e.printStackTrace();
                } catch (IncompatibleScoringSchemeException e) {
                    e.printStackTrace();
                }                                        

            } 
        });
        System.out.println("Num tags: " + tags.size() + ", Num refTags: " + refTagDataList.size() + ", TotalTime for calculateTagRefTagAlignment was " + (System.nanoTime() - totalTime) / 1e9 + " seconds");
    }
    
    // This calculates alignment of  reference tags against each other.
    private void calculateRefRefAlignment(List<RefTagData> refTags, Multimap<Tag,Position> refTagPosMap,
            Multimap<RefTagData,AlignmentInfo> refTagAlignInfoMap){
        long totalTime = System.nanoTime();
        // For each tag on the reftags list, run SW against all other tags in the list
       refTags.parallelStream().forEach(tag1 -> {
            for (RefTagData tag2: refTags) {
                if (tag1.equals(tag2)) continue; // don't align against yourself
                String seq1 = tag1.tag().sequence();
                String seq2 = tag2.tag().sequence();
                Reader reader1 = new StringReader(seq1);
                Reader reader2 = new StringReader(seq2);

                PairwiseAlignmentAlgorithm  algorithm;
                ScoringScheme  scoring;
                algorithm = new SmithWaterman();
                scoring = new BasicScoringScheme(match_reward(),mismatch_penalty(),gap_penalty());
                algorithm.setScoringScheme(scoring);
                int score = 0;

                try {
                    algorithm.loadSequences(reader1, reader2);
                    // for reftag-reftag alignment, we are only computing the score
                    score = algorithm.getScore();
                } catch (IOException ioe) {
                    ioe.printStackTrace();
                } catch (InvalidSequenceException ise) {
                    ise.printStackTrace();
                } catch (IncompatibleScoringSchemeException isse) {
                    isse.printStackTrace();
                }
                // for reftag/reftag, we have no alignment position .  Store -1.  
                // both alignment positions and reference strand (which is 1 for both) are ignored params
                // for ref-ref alignment.
                AlignmentInfo tagAI = new AlignmentInfo(tag2.tag(),tag2.chromosome(),tag2.position(),-1, 1,refGenome(),score);
                refTagAlignInfoMap.put(tag1,tagAI);
            }
        }); 
        System.out.println("Number of refTags: " + refTags.size() + ", TotalTime for calculateREfRefAlignment was " + (System.nanoTime() - totalTime) / 1e9 + " seconds");
    }
    // Checks for non ACGT character in kmer.  If found, returns the deviate
    // position in the array.  If not found, -1 is returned. -1 is good return
    private int checkForN(byte[] kmer) {
        boolean flag = false;
        int idx;
       // byte[] kmerInBytes = NucleotideAlignmentConstants.convertHaplotypeStringToAlleleByteArray(kmer);
        for (idx = 0; idx < kmer.length; idx++) {
            if (kmer[idx] >3) { // non ACGT value
                flag = true;
                break;
            }
        }
        if (flag) return idx; // toss sequences with non ACGT value
        else return -1;
    }

    
    // write to a file for debug
    public static void writeToFile(String chrom, Multimap<String,Integer> chromMaximaMap, int minCount) {
        String outFile = "/Users/lcj34/notes_files/repgen/junit2_out/chrom" + chrom + "_peakMaximaPositions_minCount"
                    + minCount + ".txt";
        BufferedWriter bw = Utils.getBufferedWriter(outFile);
        // Print the positions for this chrom
        Collection<Integer> maximaForChrom = chromMaximaMap.get(chrom);
        List<Integer> chromValues = new ArrayList<Integer>(maximaForChrom);
        Collections.sort(chromValues);
        System.out.println("Total number of maxima/peak positions for chrom " + chrom 
                + " is " + maximaForChrom.size() + ". Writing them to file " + outFile + " .... ");                              
        try {
            StringBuilder sb = new StringBuilder();
            int count = 0;
            for (int idx = 0; idx < chromValues.size(); idx++) {
                String hits = Integer.toString(chromValues.get(idx));
                sb.append(hits);
                sb.append("\n");
                count++;
                if (count > 100000) {
                    bw.write(sb.toString());
                    count = 0;
                    sb.setLength(0);
                }                       
            }
            if (count > 0) {
                // write last lines
                bw.write(sb.toString());
            }
            bw.close();
        } catch (IOException ioe) {
            System.out.println("LCJ - exception writing file " + outFile);
            ioe.printStackTrace();
        } 
    }
    
    public static void writePeakPositions(String chrom, int peak, List<Integer> peakPositions) {
        String outFile = "/Users/lcj34/notes_files/repgen/junit2_out/chrom" + chrom + "_peak" + peak + "_positions.txt";
                    
        BufferedWriter bw = Utils.getBufferedWriter(outFile);
        // Print the positions for this chrom
        
        System.out.println("Total number of  positions for peak " + peak
                + " is " + peakPositions.size() + ". Writing them to file " + outFile + " .... ");                              
        try {
            StringBuilder sb = new StringBuilder();
            int count = 0;
            for (int idx = 0; idx < peakPositions.size(); idx++) {
                String hits = Integer.toString(peakPositions.get(idx));
                sb.append(hits);
                sb.append("\n");
                count++;
                if (count > 100000) {
                    bw.write(sb.toString());
                    count = 0;
                    sb.setLength(0);
                }                       
            }
            if (count > 0) {
                // write last lines
                bw.write(sb.toString());
            }
            bw.close();
        } catch (IOException ioe) {
            System.out.println("LCJ - exception writing file " + outFile);
            ioe.printStackTrace();
        } 
    }

    public static void writeChrom9Bits(List<Integer> peakPositions) {
        String outFile = "/Users/lcj34/notes_files/repgen/junit_out/chrom9bit_positions.txt";
                    
        BufferedWriter bw = Utils.getBufferedWriter(outFile);
        // Print the positions for this chrom
        
        System.out.println("Total number of  positions for chrom9: " + peakPositions.size());                            
        try {
            StringBuilder sb = new StringBuilder();
            int count = 0;
            for (int idx = 0; idx < peakPositions.size(); idx++) {
                String hits = Integer.toString(peakPositions.get(idx));
                sb.append(hits);
                sb.append("\n");
                count++;
                if (count > 100000) {
                    bw.write(sb.toString());
                    count = 0;
                    sb.setLength(0);
                }                       
            }
            if (count > 0) {
                // write last lines
                bw.write(sb.toString());
            }
            bw.close();
        } catch (IOException ioe) {
            System.out.println("LCJ - exception writing file " + outFile);
            ioe.printStackTrace();
        } 
    }

    @Override
    public ImageIcon getIcon() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getButtonName() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getToolTipText() {
        // TODO Auto-generated method stub
        return null;
    }
 
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
     public static void main(String[] args) {
         GeneratePluginCode.generate(RepGenAlignerPlugin.class);
     }

    /**
     * Input database file with tags and taxa distribution
     *
     * @return Input DB
     */
    public String inputDB() {
        return myDBFile.value();
    }

    /**
     * Set Input DB. Input database file with tags and taxa
     * distribution
     *
     * @param value Input DB
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin inputDB(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }

    /**
     * Output fastq file to use as input for BWA or bowtie2
     *
     * @return Output File
     */
    public String refGenome() {
        return refGenome.value();
    }

    /**
     * Set Output File. Output fastq file to use as input
     * for BWA or bowtie2
     *
     * @param value Output File
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin refGenome(String value) {
        refGenome = new PluginParameter<>(refGenome, value);
        return this;
    }

    /**
     * Minimum count of reads for a tag to be output
     *
     * @return Min Count
     */
    public Integer minTagCount() {
        return minTagCount.value();
    }

    /**
     * Set Min Count. Minimum count of reads for a tag to
     * be output
     *
     * @param value Min Count
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin minTagCount(Integer value) {
        minTagCount = new PluginParameter<>(minTagCount, value);
        return this;
    }
    
    /**
     * Minimum count of hits in a cluster for a reference tag
     * to be created for this cluster
     *
     * @return Min Count
     */
    public Integer minHitCount() {
        return minHitCount.value();
    }

    /**
     * Set Min Count. Minimum count of hits in a cluster for a 
     * reference tag to be created fro this cluster
     *
     * @param value Min Count
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin minHitCount(Integer value) {
        minHitCount = new PluginParameter<>(minHitCount, value);
        return this;
    }
    /**
     * Length of seed kmers
     *
     * @return seed len
     */
    public Integer seedLen() {
        return seedLen.value();
    }

    /**
     * Set Seed kmer length. Length of seeds to 
     * use in aligning.
     *
     * @param value seed length
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin seedLen(Integer value) {
        seedLen = new PluginParameter<>(seedLen, value);
        return this;
    }
    
    /**
     * Length of window between positions 
     * when creating seed from DB tags
     *
     * @return seed window
     */
    public Integer seedWindow() {
        return seedWindow.value();
    }

    /**
     * Set Seed window length. Length of window 
     * between positions when creating seed from DB tags
     *
     * @param value seed length
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin seedWindow(Integer value) {
        seedWindow = new PluginParameter<>(seedWindow, value);
        return this;
    }
    /**
     * Length of  kmers as tag sequences in the db
     *
     * @return kmerLen
     */
    public Integer kmerLen() {
        return kmerLen.value();
    }

    /**
     * Set Kmer length. Length of kmers to be stored
     * as tag sequences in the db
     *
     * @param value kmer length
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin kmerLen(Integer value) {
        kmerLen = new PluginParameter<>(kmerLen, value);
        return this;
    }
    /**
     * Length of  kmers as tag sequences in the db
     *
     * @return kmerLen
     */
    public Integer refKmerLen() {
        return refKmerLen.value();
    }

    /**
     * Set Kmer length. Length of kmers to be stored
     * as tag sequences in the db
     *
     * @param value kmer length
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin refKmerLen(Integer value) {
        refKmerLen = new PluginParameter<>(refKmerLen, value);
        return this;
    }
    
    /**
     * Parameter sent to Smith Waterman aligner for use in
     * calculating reward when base pairs match.
     *
     * @return Match Reward Amount
     */
    public Integer match_reward() {
        return match_reward.value();
    }

    /**
     * Set Match Reward Amount. Parameter sent to Smith Waterman
     * aligner for use in calculating reward when base pairs
     * match.
     *
     * @param value Match Reward Amount
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin match_reward(Integer value) {
        match_reward = new PluginParameter<>(match_reward, value);
        return this;
    }

    /**
     * Parameter sent to Smith Waterman aligner for use in
     * calculating penalty when base pairs are mis-matched.
     *
     * @return Mismatch Penalty Amount
     */
    public Integer mismatch_penalty() {
        return mismatch_penalty.value();
    }

    /**
     * Set Mismatch Penalty Amount. Parameter sent to Smith
     * Waterman aligner for use in calculating penalty when
     * base pairs are mis-matched.
     *
     * @param value Mismatch Penalty Amount
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin mismatch_penalty(Integer value) {
        mismatch_penalty = new PluginParameter<>(mismatch_penalty, value);
        return this;
    }

    /**
     * Parameter sent to Smith Waterman aligner for use in
     * calculating penalty when when a gap is identified.
     *
     * @return Gap Penalty Amount
     */
    public Integer gap_penalty() {
        return gap_penalty.value();
    }

    /**
     * Set Gap Penalty Amount. Parameter sent to Smith Waterman
     * aligner for use in calculating penalty when when a
     * gap is identified.
     *
     * @param value Gap Penalty Amount
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin gap_penalty(Integer value) {
        gap_penalty = new PluginParameter<>(gap_penalty, value);
        return this;
    }
}
