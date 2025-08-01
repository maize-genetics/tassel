/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import java.awt.Frame;
import java.io.BufferedReader;
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
import java.util.HashMap;
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
import net.maizegenetics.dna.BaseEncoder;
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
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;

/**
 * This plugin takes an existing repGen db, grabs the tags
 * whose depth meets that specified in the minCount parameter,
 * makes kmer seeds from these tags.
 * 
 * Forward and reverse primer sequences are added as an input parameter.
 * When a kmer seed is found on a reference chromosome, a ref sequence
 * is created from 300bp before the hit, to 300 bp after the hit.  This
 * value is half the refKmerLen parameter passed by user.  Default
 * refKmerLen is 600.
 * 
 * From the ref sequence created, a search is made for the primer
 * pairs within this sequence.  IF either both forward primer and the
 * reverse complement of the reverse primer; or reverse primer and the
 * reverse complement of the forward primer are found, a reference tag
 * is created starting at the start of the first occurring primer from
 * the primer pair found in the sequence.  If both forward and reverse
 * pairs are found, the ref tag is created based on the best match,
 * defaulting to the forward primer if both are found.
 * 
 * Search for additional kmer matches on the chromosome begins at
 * the position on the ref chrom following the end of the second primer in the matched
 * pair.
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
public class RepGenPhase2AlignerPlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(RepGenPhase2AlignerPlugin.class);
    
    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with tags and taxa distribution").build();
    private PluginParameter<String> refGenome = new PluginParameter.Builder<String>("ref", null, String.class).guiName("Reference Genome File").required(true)
            .description("Referemce Genome File for aligning against ").build();
    private PluginParameter<Integer> minTagCount = new PluginParameter.Builder<Integer>("minTagCount", 1, Integer.class).guiName("Min Tag Count")
            .description("Minimum count of reads for a tag to be aligned").build();
    private PluginParameter<Integer> seedLen = new PluginParameter.Builder<Integer>("seedLen", 31, Integer.class).guiName("Seed Kmer Length")
            .description("Length of kmer seed created from DB tags and used as seeds for aligning against the reference genome.").build();
    private PluginParameter<Integer> seedWindow = new PluginParameter.Builder<Integer>("seedWindow", 1, Integer.class).guiName("Window Length for Seed Creation")
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
    private PluginParameter<String> primers = new PluginParameter.Builder<String>("primers", null, String.class).guiName("Primers").required(true).inFile()
            .description("Tab delimited file that contains a list of forward,reverse primer pairs.  \nThe values in each column are the forward primer sequence and the reverse primer sequence.").build();
 
    
    static GenomeSequence myRefSequence = null;
    // length of ref tag sequence from which to search for primer strings
    // This is half the length
    //static int refAlignLen = 600; // could be a plugin parameter - hard code for testing
    static int refAlignLen = 1000;
    static int bestScoreBad = 0;

    public enum BestScore {
        none(0), forward(1),reverse(2), forwardRC(3), reverseRC(4);
        private int primer;
        
        private BestScore(int value) {
            this.primer = value;
        }
    }
    public RepGenPhase2AlignerPlugin() {
        super(null, false);
    }

    public RepGenPhase2AlignerPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public RepGenPhase2AlignerPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    @Override
    public void postProcessParameters() {

        if (myDBFile.isEmpty() || !Files.exists(Paths.get(inputDB()))) {
            throw new IllegalArgumentException("RepGenPhase2AlignerPlugin: postProcessParameters: Input DB not set or found");
        }
        if (primers.isEmpty() || !Files.exists(Paths.get(inputDB()))) {
            throw new IllegalArgumentException("RepGenPhase2AlignerPlugin: postProcessParameters: Primer file not set or found");
        }
        if (!refGenome.isEmpty()) {
            myRefSequence = GenomeSequenceBuilder.instance(refGenome());
        } else {
            throw new IllegalArgumentException("RepGenPhase2AlignerPlugin: postProcessParameters: reference genome not set or found");
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        long totalTime = System.nanoTime();
        long time=System.nanoTime();
 
        try {           
            System.out.println("RepGenPhase2AlignerPlugin:processData begin"); 
            RepGenDataWriter repGenData=new RepGenSQLite(inputDB());
            Map<Tag, Integer> tagsWithDepth = repGenData.getTagsWithDepth(minTagCount());
            if (tagsWithDepth.isEmpty()) {
                System.out.println("\nNo tags found with minimum depth " + minTagCount() + ". Halting Run.\n");
                ((RepGenSQLite)repGenData).close();
                return null;
            }
            
            // DEBUG _ PRINTING TAGS BELOW
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
            // myRefSequence populated in post process parameters           
            Set<Chromosome> chromsInRef = myRefSequence.chromosomes();
            
            List<Tuple<String,String>> primerList = createPrimerList(); // get list of primers from input file
            if (primerList == null || primerList.isEmpty()) {
                ((RepGenSQLite)repGenData).close();
                System.out.println("Failed to process entries from primers file " + primers() 
                + ".\nPlease ensure file exists and a header line with 3 tab-delimited columns for chrom,forward,reverse\n");
                return null;
            }
                      
           //  Create synchronized map for use in parallel streams 
            Multimap<Tag,Position> refTagPositionMap = Multimaps.synchronizedMultimap(HashMultimap.<Tag,Position>create());

            chromsInRef.parallelStream().forEach(chrom -> { 
                  // Turn this on/off for debug purposes
                  //if (chrom.getChromosomeNumber() != 9) return; // just for initial testing !!! - remove
                                    
                  int kmersForChrom = 0;
                  int chromLength = myRefSequence.chromosomeSize(chrom);
                  System.out.println("\nChecking reference chrom " + chrom.getName() + ", size: " + chromLength + " for tag kmer matches.");
                  
                  int refTagsCreated = 0;
                  int noRefTagCreated = 0;
                  for (int chromIdx = 0; chromIdx < chromLength;) {                                       
                      // chromosomeSequence:  start and end are inclusive and 1-based, adjust for this in call below
                      // grab seq from ref the length of the kmer seed - check if refSeq matches one of the kmer seeds
                      int end = Math.min(chromLength, chromIdx+seedLen());
                      byte[] chromSeedBytes = myRefSequence.chromosomeSequence(chrom,chromIdx+1, end);
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
                          // Should chromIdx or chromIdx+1 be passed? 
                          // grab ref bytes to use when looking for primers
                          // This is a longer string than above - it is 300 above and 300 below
                          // where we found a kmer seed match (300 assuming the refAlignLen was 600)
                          //
                          // Update from Dan:  expect there to be 200 or less bps' between the
                          // primer ends.  Here we tested with 300 up and 300 below.  That should
                          // have covered it.  Let's try 1000 to see if we got better
                          // then the 600 we originally hard coded
                          int refHalfLen = refAlignLen/2;
                          int first= (chromIdx < (refHalfLen)) ? 0 : chromIdx;
                          int last = Math.min(chromLength, chromIdx + refHalfLen);
                          byte[] refBytes = myRefSequence.chromosomeSequence(chrom,first+1, last);
                          String refString = NucleotideAlignmentConstants.nucleotideBytetoString(refBytes);
                          if (refString == null) {
                              // bad ref - 
                              System.out.println("repGenPhase2ALigner:procesData - NULL returned for refString - continue");
                              continue; 
                          }
                          int nextChromIdx = createRefTagsForAlignment(chrom, refString, first+1, primerList, refTagPositionMap);
                          // if we found and created a tag, nextChromIdx tells us at what position within
                          // the reference string the "end" primer ends.  Add this value to the current
                          // chromIdx to begin the next kmer seed search.
                          // If we could not create a ref tag from this position, slide by 1 and restart
                          // the process of looking for seed kmer matches.  If we didn't find start/end primers
                          // in the 600bp window surrounding this match, will we find them by moving up just 1?
                          // Should the default when refTag can't be created be to increment by something greater than 1?
                          if (nextChromIdx > 0) {
                              chromIdx += nextChromIdx;
                              refTagsCreated++;
                          }
                          else {
                              chromIdx++;
                              noRefTagCreated++;
                          }
                      } else {
                          chromIdx++; // look for hit in next position
                      }                       
                  }
                  
                  System.out.println("Total tag seeds matching to kmers in chrom " + chrom.getName() + ": " 
                      + kmersForChrom + ",total refTagsCreated " + refTagsCreated + ", no ref tag created " + noRefTagCreated);  
              });
                   
            // LCJ - below is debug!  Writes to my cbsu dir
            String refTagOutFile = "/home/lcj34/repGen_outputFiles/allChroms_refTagsFromPrimers.txt";
//            // LCJ - this one writes to laptop
//            //String refTagOutFile = "/Users/lcj34/notes_files/repgen/junit_phase2Aligner/chrom" + chrom.getName() + "_refTagsFromPrimers.txt";
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
            // LCJ - end debug
            
            System.out.println("Finished with refs from primers, total time:" + (System.nanoTime() - time)/1e9 + " seconds.\n");
            time = System.nanoTime();
            
            System.out.println("\nNumber of distinct refTags to be loaded into db: " 
              + refTagPositionMap.keySet().size() + " total refTags: " + refTagPositionMap.size());
            // add ref and mapping approach to db, then reference tags
            repGenData.addMappingApproach("SmithWaterman");
            repGenData.addReferenceGenome(refGenome());
            
            // Add ref tag and tags separately. Add to tagAlignment map,
            repGenData.putRefTagMapping(refTagPositionMap, refGenome());
            
//            // Get tags stored from RepGenLoadSeqToDB        
//            Set<Tag> tagsToAlign = repGenData.getTags();
//            List<Tag> tagList = new ArrayList<Tag>(tagsToAlign);
//            
//            // Create synchronized map for use in parallel streams 
//            Multimap<Tag,AlignmentInfo> tagAlignInfoMap = Multimaps.synchronizedMultimap(HashMultimap.<Tag,AlignmentInfo>create());
//            Multimap<RefTagData,AlignmentInfo> refTagAlignInfoMap = Multimaps.synchronizedMultimap(HashMultimap.<RefTagData,AlignmentInfo>create());
//            
//            // First align the tags against each other
//            // Output of this is stored in db table tagAlignments
//            System.out.println("Calling calculateTagTagAlignment for tags without refs");
//            
//            calculateTagTagAlignment( tagList, tagAlignInfoMap);
//            System.out.println("Number of tag-tag alignments: " + tagAlignInfoMap.size() + ", store to db.");
//            // neither tag is reference, null and -1 for tag1 chrom/pos
//            repGenData.putTagAlignments(tagAlignInfoMap, false,false, null,-1,refGenome());
// 
//            System.out.println("Calling calculateTagRefTagALignment for tags WITH ref");
//            Set<RefTagData> refTags = repGenData.getRefTags();
//            List<RefTagData> refTagList = new ArrayList<RefTagData>(refTags);
//            
//            tagAlignInfoMap.clear(); // remove old data
//            calculateTagRefTagAlignment(tagList,refTagList,tagAlignInfoMap,refGenome());
//            
//            // Add the tag-refTag info to the tagAlignments table.           
//            // tag1=nonref, tag2=ref, null and -1 for tag1 .  Alignment will be
//            // done twice:  once for tag/refTag-fwd, once for tag/refTag-reverse
//            System.out.println("Number of tag-refTag alignments, includes aligning to ref fwd and reverse strands: " + tagAlignInfoMap.size() + ", store to db.");
//            repGenData.putTagAlignments(tagAlignInfoMap,false,true,null,-1,refGenome());
//            
//            System.out.println("Calling calculateRefRefAlignment");            
//            calculateRefRefAlignment(refTagList,refTagPositionMap,refTagAlignInfoMap);
//            System.out.println("Number of reftag-reftag alignments: " + refTagAlignInfoMap.size() + ", store to db.");
//            repGenData.putRefRefAlignments(refTagAlignInfoMap, refGenome()); // CREATE putRefRefAlignments -  different parameters
            ((RepGenSQLite)repGenData).close();

            myLogger.info("Finished RepGenPhase2AlignerPlugin\n");
        } catch (Exception e) {
            myLogger.info("Catch in prcessing RepGenPhase2AlignerPlugin file e=" + e);
            e.printStackTrace();
        }
        System.out.println("Process took " + (System.nanoTime() - totalTime)/1e9 + " seconds.\n");
        return null;
    }
    
    private List<Tuple<String,String>> createPrimerList() {
        // forward and reverse primer remain together,  store in hashmap
        // could be stored as String[][]
        List<Tuple<String, String>> pmap = new ArrayList<Tuple<String,String>>();
        BufferedReader br = Utils.getBufferedReader(primers());
        
        try {
            String line = br.readLine(); // ditch the header line
            while((line = br.readLine()) != null) {
                String[] primerTokens = line.split("\t");
                if (primerTokens.length != 2) {
                    System.out.println("Bad file format for primers.  Should be 2 tab delimited columns for forward, reverse");
                    return null;
                }
                Tuple<String, String> primers = new Tuple<String,String>(primerTokens[0],primerTokens[1]);
                pmap.add(primers);
            }
        } catch (IOException ioe) {
            System.out.println("createChromPrimerMap failed processing file: " + primers());
            ioe.printStackTrace();
            return null;
        }
        return pmap;
    }
    
    private void createKmerSeedsFromDBTags(Set<Tag> tagWithDepth,
            Multimap<String,Tag> kmerTagMap, int window) {
        for (Tag tag : tagWithDepth) {
            String tagSequence = tag.sequence();
            // maxIdx ensures we don't substring beyond the end of the tag.
            int maxIdx = window > seedLen() ? tagSequence.length() - window : tagSequence.length() - seedLen();
            for (int seqIdx = 0; seqIdx < maxIdx;) {
                String kmer = tagSequence.substring(seqIdx, seqIdx + seedLen());  
                byte[] kmerBytesAsNum = NucleotideAlignmentConstants.convertHaplotypeStringToAlleleByteArray(kmer);
                
                int badValues = checkForN(kmerBytesAsNum);
                if (badValues >=0) {
                    // found a non-ACGT character,  re-set index and
                    // grab new kmer
                    seqIdx += badValues+1; // skip past the bad value
                    // just increase by window size, don't care where the bad value is
                    //seqIdx += window;
                    continue;
                }
                // good values - add unique kmer to map
                kmerTagMap.put(kmer, tag);
                // add reverse complement of kmer
                byte[] kmerRC = NucleotideAlignmentConstants.reverseComplementAlleleByteArray(kmer.getBytes());
                String kmerRCString = new String(kmerRC);
                kmerTagMap.put(kmerRCString, tag); 
                // loop end - increment seqIdx for next kmer
                seqIdx += window;
            }               
        }
    }

    private int createRefTagsForAlignment( Chromosome chrom, String refString, int refOffset, List<Tuple<String,String>> primers,
             Multimap<Tag,Position> refTagPositionMap) {
        
        int nextChromIdx = -1; // what is returned
        
        // A kmer match was found on the reference chromosome.  The refSTring is 300 down
        //  to 300 up from the start of the seed kmer match.  (half the refKmerLen parameter)
        // Within this string, look for primer matches.
        // From each primer pair on the list, we need to match the
        // forward primer to the reverse complement of the reverse primer.  And we match
        // the reverse primer string to the reverse complement of the forward primer string.
        // If both matches of a pair are found, with =< 3bp mismatch, then create a reference.  
        // Use SW to find the matches
       
        for (Tuple<String,String> primer: primers) {
            String forward = primer.x;
            String reverse = primer.y;
            String forwardRC = BaseEncoder.getReverseComplement(forward);
            String reverseRC = BaseEncoder.getReverseComplement(reverse);
            
            // Assuming all primers have the same length!  For this method,
            // match value=2, both mismatch and gap get a -1. With these
            // values, the perfect scores and 1,2,3 mismatches are calculated below
            int perfectScore = forward.length() * 2;
            int oneMisMatch = ((forward.length() - 1) * 2) -1;
            int twoMisMatch = ((forward.length() - 2) * 2) -2;
            int threeMisMatch = ((forward.length() - 3) * 2) -3;

            Tuple<Integer,Integer> forwardPrimer = null;
            Tuple<Integer,Integer> reversePrimer = null;
            Tuple<Integer,Integer> forwardRCPrimer = null;
            Tuple<Integer,Integer> reverseRCPrimer = null;
            
            // Lynn - does this make sense, or is it just extra
            // computing?  If forward is best, then see if we have
            // reverseRC as not null.  If reverse is best, check if
            // forwardRC is not null, etc.  Maybe a case to handle them
            // But we still need to know if other side appears further
            // up the chrom.  if it doesn't, don't we skip it?
            // Since the reference is always on the forward strand, should
            // we look for values in that order?
            // Need to look at PCR videos again.  WHy is forward mapped to reverseRC?
            BestScore bestPrimer = BestScore.none;
            int bestScore = -1;
            Map<String,PairwiseAlignment> pAlign = new HashMap<String,PairwiseAlignment>();
            forwardPrimer = computePrimerSW(forward,  refString, 2,-1,-1,pAlign);
            if (forwardPrimer != null && forwardPrimer.x >= threeMisMatch) {
                bestPrimer = BestScore.forward;
                bestScore = forwardPrimer.x;
                reverseRCPrimer = computePrimerSW(reverseRC,  refString, 2,-1,-1,pAlign);
                if (reverseRCPrimer !=null && reverseRCPrimer.x >= threeMisMatch) {
                    if (reverseRCPrimer.x > bestScore) {
                        bestScore = reverseRCPrimer.x;
                        bestPrimer = BestScore.reverseRC;
                    }
                } else {
                    reverseRCPrimer = null;
                    bestScore = -1; // reset it - we don't have a primer pair
                }
            } else forwardPrimer = null; // not a good score!

            reversePrimer = computePrimerSW(reverse,  refString, 2,-1,-1,pAlign);
            if (reversePrimer != null && reversePrimer.x >= threeMisMatch) {
                int tempBest = bestScore;
                if (reversePrimer.x > tempBest) {
                    tempBest = reversePrimer.x;
                    bestPrimer = BestScore.reverse;
                }
                forwardRCPrimer = computePrimerSW(forwardRC,  refString, 2,-1,-1,pAlign);
                if (forwardRCPrimer !=null && forwardRCPrimer.x >= threeMisMatch) {
                    if (forwardRCPrimer.x > tempBest) {
                        tempBest = forwardRCPrimer.x;
                        bestPrimer = BestScore.forwardRC;
                    }
                    if (tempBest > bestScore) bestScore = tempBest;
                } else forwardRCPrimer = null;
            } else reversePrimer = null;
 
            if (bestScore < 0) {
                // no pairs of primers found with 3 or less mismatches from ref.
                continue; // check next primer pair
            }
            // At this point, whatever is recorded as bestPrimer can be used
            // if bestScore is > 0.  If we didn't have a pair, bestScore would be -1
            
            // Create a reference.  Remember - once we have the reference on the
            // refTag list, we set the refIndex to beyond the end primer and we
            // start looking from there to find next seed tag match, and once
            // that is found we call into here again.
            int refBegin = -1;
            // The beginning of the reference tag is the index in the ref where
            // the primer aligned.  The end of the reference tag is beginning plus
            // refTagLength input parameter from user.
            // Where we begin looking for the next seed (the return value) is at the end 
            // of the end primer
            switch (bestPrimer) {
            case forward:
            case reverseRC:
                // debug
//                System.out.println("BestPrimers: forwardPrimer.score: " + forwardPrimer.x +
//                        ", forwardPrimer.refBegin: " + forwardPrimer.y + ", reverseRC.score: "
//                        + reverseRCPrimer.x + ", reverseRC.refBegin: " + reverseRCPrimer.y);
//                if (forwardPrimer.y < 0 || reverseRCPrimer.y < 0) {
//                    System.out.println("ForwardPrimer: " + forward);
//                    System.out.println("ReverseRCPrim: " + reverseRC);
//                    System.out.println("Alignment offset, forward:" + forwardPrimer.y 
//                            + ", reverseRC:" + reverseRCPrimer.y + ", refOffset param:" + refOffset);
//                    if (forwardPrimer.y < 0) {
//                        PairwiseAlignment testAlign = pAlign.get(forward);
//                        System.out.println("Forward Alignment: \n" + testAlign);
//                    }
//                    if (reverseRCPrimer.y < 0) {
//                        PairwiseAlignment testAlign = pAlign.get(reverseRC);
//                        System.out.println("REverseRC Alignment: \n" + testAlign);
//                    }
//                }
                // end debug
                if (forwardPrimer.y < reverseRCPrimer.y) {
                    refBegin =  forwardPrimer.y;
                    nextChromIdx = reverseRCPrimer.y + reverseRC.length();
                } else {
                    refBegin =  reverseRCPrimer.y;
                    nextChromIdx = forwardPrimer.y + forward.length();
                }
                break;
            case reverse:
            case forwardRC:               
                // debug
//                System.out.println("BestPrimers: reversePrimer.score: " + reversePrimer.x +
//                        ", reversePrimer.refBegin: " + reversePrimer.y + ", forwardRC.score: "
//                        + forwardRCPrimer.x + ", forwardRC.refBegin: " + forwardRCPrimer.y);
//                if (forwardRCPrimer.y < 0 || reversePrimer.y < 0) {
//                    System.out.println("ForwardRCPrimer: " + forwardRC);
//                    System.out.println("ReversePrim: " + reverse);
//                    System.out.println("Alignment offset, forwardRC:" + forwardRCPrimer.y 
//                            + ", reverse:" + reversePrimer.y + ", refOffset param:" + refOffset);
//                    if (forwardRCPrimer.y < 0) {
//                        PairwiseAlignment testAlign = pAlign.get(forwardRC);
//                        System.out.println("ForwardRC Alignment: \n" + testAlign);
//                    }
//                    if (reversePrimer.y < 0) {
//                        PairwiseAlignment testAlign = pAlign.get(reverse);
//                        System.out.println("Reverse Alignment: \n" + testAlign);
//                    }
//                }
                // end debug
                if (reversePrimer.y < forwardRCPrimer.y) {
                    refBegin =  reversePrimer.y;
                    nextChromIdx = forwardRCPrimer.y + forwardRC.length();
                } else {
                    refBegin =  forwardRCPrimer.y;
                    nextChromIdx = reversePrimer.y + reverseRC.length();
                }
                break;
            default:
                System.out.println("ERROR: best Primer is not one of the 4 !! We should have ditched earlier!");
                 return -2;
            }
            // Create reference tag - store to map
            // refBeginIdx holds the starting position in the reference chromosome.  It is a 1-based number
            // from where we created the refString that was passed in.  Add to this
            // the offset from within the string where the primer was found.  Create the
            // reference tag from there.
            int refStartIdx = refOffset + refBegin;
            int refEndIdx = Math.min(refStartIdx + refKmerLen(), myRefSequence.chromosomeSize(chrom));
            byte[] refTagBytes = myRefSequence.chromosomeSequence(chrom,refStartIdx, refEndIdx);
           // System.out.println("LCJ - after finding best primer for chrom " + chrom.getName() + " calling myRefSequence to get refString with refBegin value " + refBegin);

            // Convert to allele string
            String refTagString = NucleotideAlignmentConstants.nucleotideBytetoString(refTagBytes);
            
            // RefString could be null for some other reason - e.g. start/end pos is bad
            if (refTagString == null || refTagString.length() == 0 || refTagString.contains("N") || refTagString.contains("null")) {
                System.out.println(" storeRefTagPositions - refString is NULL or contains N or null or is length 0");
            } else {
                Tag refTag = TagBuilder.instance(refTagString).reference().build();
                if (refTag != null ) {
                    Position refPos=new GeneralPosition
                            .Builder(chrom,refStartIdx)
                            .strand((byte)1)
                            .addAnno("mappingapproach", "SmithWaterman") // this wasn't done with SW !!  Is just a reference
                            .addAnno("forward", "true") // reference, so always forward strand.
                            .build();
                    // This is a multimap because it is possible for a tag sequence to
                    // show up in multiple places in the genome.  The multimap has all positions
                    // for this tag.
                    refTagPositionMap.put(refTag, refPos);
                } else {
                    System.out.println("- refTag is NULL for refString: " + refString);
                } 
            }           
        }
        return nextChromIdx;
    }
    
    // this calculates tag against tag alignment from the tags in the DB
    // tag table (not the refTag table)
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
                    alignment = algorithm.getPairwiseAlignment(); // compute alignment
                    score = algorithm.getScore(); // get score - this is done in getPairwiseAlignment
                  
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

    // This method computes the SW score of a primer or primer-rc sequence
    // against a 600 bp ref tag window.  It returns a Tuple<Integer,Integer>
    // were Tuple.x is the alignment score, and tuple.y is the starting position
    // on the reference tag.  These values will be used to determine whether
    // the primer sequence was found.
    private Tuple<Integer,Integer> computePrimerSW(String primerSeq, String refSeq, 
            int match, int mismatch, int gap, Map<String, PairwiseAlignment> pAlignment) {
        
        Reader reader1 = new StringReader(primerSeq);
        Reader reader2 = new StringReader(refSeq);

        PairwiseAlignmentAlgorithm  algorithm;
        ScoringScheme  scoring;
        PairwiseAlignment alignment;

        algorithm = new SmithWaterman();
        scoring = new BasicScoringScheme(match,mismatch,gap);
        algorithm.setScoringScheme(scoring);
        int score;

        int primerOffset = 0;
        int refAlignStartPos = 0;
        try {
            algorithm.loadSequences(reader1, reader2);
            alignment = algorithm.getPairwiseAlignment(); // compute alignment
            score = algorithm.getScore();

            // REMOVE THIS LINE and the map parameter
             pAlignment.put(primerSeq, alignment); // THIS IS DEBUG - remove it and the parameter
              
            // The first sequence given in loadSequences() is loaded into rows. This is primerSeq
            // The second is loaded into columns (matrix[seq1][seq2] - this is the refSeq
            primerOffset = alignment.getRowStart(); // if not 0, SW clipped the primer
            refAlignStartPos += alignment.getColStart(); // This probably isn't 0 as refString is longer than primer
            
            if (primerOffset > 0) {
                // Tag was not aligned from the beginning,
                // add back the bps that were skipped so alignment begins at start of the tag
                refAlignStartPos -= primerOffset;
            }
            return (new Tuple<Integer,Integer>(score,refAlignStartPos));
                                                  
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InvalidSequenceException e) {
            e.printStackTrace();
        } catch (IncompatibleScoringSchemeException e) {
            e.printStackTrace();
        }                                           
        return null; // error computing the values
    }
    
    // write to a file for debug
    public static void writeToFile(String chrom, Multimap<String,Integer> chromMaximaMap, int minCount) {
        String outFile = "/Users/lcj34/notes_files/repgen/junit_out/chrom" + chrom + "_peakMaximaPositions_minCount"
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
        String outFile = "/Users/lcj34/notes_files/repgen/junit_out/chrom" + chrom + "_peak" + peak + "_positions.txt";
                    
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
         GeneratePluginCode.generate(RepGenPhase2AlignerPlugin.class);
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
    public RepGenPhase2AlignerPlugin inputDB(String value) {
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
    public RepGenPhase2AlignerPlugin refGenome(String value) {
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
    public RepGenPhase2AlignerPlugin minTagCount(Integer value) {
        minTagCount = new PluginParameter<>(minTagCount, value);
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
    public RepGenPhase2AlignerPlugin seedLen(Integer value) {
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
    public RepGenPhase2AlignerPlugin seedWindow(Integer value) {
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
    public RepGenPhase2AlignerPlugin kmerLen(Integer value) {
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
    public RepGenPhase2AlignerPlugin refKmerLen(Integer value) {
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
    public RepGenPhase2AlignerPlugin match_reward(Integer value) {
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
    public RepGenPhase2AlignerPlugin mismatch_penalty(Integer value) {
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
    public RepGenPhase2AlignerPlugin gap_penalty(Integer value) {
        gap_penalty = new PluginParameter<>(gap_penalty, value);
        return this;
    }
    
    /**
     * Tab delimited file that contains the column headers
     * chrom,forward,reverse.  
     * The values in each column are the chromosone name,
     * the forward primer sequence and the reverse primer
     * sequence for the specified chromosome.
     *
     * @return Primers
     */
    public String primers() {
        return primers.value();
    }

    /**
     * Set Primers. Tab delimited file that contains the column
     * headers chrom,forward,reverse.  
     * The values in each column are the chromosone name,
     * the forward primer sequence and the reverse primer
     * sequence for the specified chromosome.
     *
     * @param value Primers
     *
     * @return this plugin
     */
    public RepGenPhase2AlignerPlugin primers(String value) {
        primers = new PluginParameter<>(primers, value);
        return this;
    }
}

