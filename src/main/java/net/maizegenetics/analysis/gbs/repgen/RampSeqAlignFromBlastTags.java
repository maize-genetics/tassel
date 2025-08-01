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
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;

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
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;

/**
 * This class takes a rAmpSeq database populated with tags,a reference Genome, a filtered
 * BLAST file output and primers, then:
 *   (1) creates  reference tags using the blast alignment output
 *   (2) runs tag/tag alignment for each tag in the DB,
 *   (3) runs tag/refTag alignment for each tag against each ref tag
 *   (4) runs refTag/refTag alignment for each refTag in the DB>
 *   (5) all alignments are stored in the db tagAlignments table.
 * 
 * Alignments are performed and stored in groups to prevent overwhelming
 * the DB with massive load commands.
 * 
 * Blast was run on CBSU using these parameters:
 * Make the reference files:
 *    makeblastdb -dbtype nucl -in <agpv4 referencene genome fasta> -parse_seqids -out maizeAGPV4.db
 *    
 * Run blast using maizeAGPV4.db from command above:
 *   blastn -num_threads 24 -db maizeAGPV4.db -query anp68R1Tags.fasta -evalue 1e-60 -max_target_seqs 5 -max_hsps 1 -outfmt 6 -out anp68TagsR1Result/blastANP68_R1.txt
 *   
 * The blast output file was filtered using the 3 commands below.  The first filters
 * identity down to 98 %,  the second gets alignment lengths that were at least 148,
 * the 3rd filters it down to just the chrom, start, end positions:
 * 
 * awk '$3 >= 98.000 {print $0}' blastANP68_R1.txt > blastANP68_R1_98per.txt
 * awk '$4 >= 148 {print $0}' blastANP68_R1_98per.txt > blastANP68_R1_98per_148align.txt
 * awk {'printf ("%s\t%s\t%s\n", $2, $9, $10)'} blastANP68_R1_98per_148align.txt > blastANP68_R1_98per_148align_3cols.txt
 * 
 * It is the last file from awk, blastANP68_R1_98per_148align_3cols.txt, that is given as a parameter here.
 * 
 * @author lcj34
 *
 */
public class RampSeqAlignFromBlastTags extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(RepGenPhase2AlignerPlugin.class);
    
    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with tags and taxa distribution").build();
    private PluginParameter<String> refGenome = new PluginParameter.Builder<String>("ref", null, String.class).guiName("Reference Genome File").required(true)
            .description("Referemce Genome File for aligning against ").build();
    private PluginParameter<Integer> minTagCount = new PluginParameter.Builder<Integer>("minTagCount", 1, Integer.class).guiName("Min Tag Count")
            .description("Minimum count of reads for a tag to be aligned").build();
    private PluginParameter<Integer> match_reward = new PluginParameter.Builder<Integer>("match_reward", 2, Integer.class).guiName("Match Reward Amount")
            .description("Parameter sent to Smith Waterman aligner for use in calculating reward when base pairs match.").build();
    private PluginParameter<Integer> mismatch_penalty = new PluginParameter.Builder<Integer>("mismatch_penalty", -1, Integer.class).guiName("Mismatch Penalty Amount")
            .description("Parameter sent to Smith Waterman aligner for use in calculating penalty when base pairs are mis-matched.").build();
    private PluginParameter<Integer> gap_penalty = new PluginParameter.Builder<Integer>("gap_penalty", -1, Integer.class).guiName("Gap Penalty Amount")
            .description("Parameter sent to Smith Waterman aligner for use in calculating penalty when when a gap is identified.").build();
    private PluginParameter<String> forwardp = new PluginParameter.Builder<String>("forwardp", null, String.class).guiName("Forward Primer").required(true)
            .description("String containing the forward primer sequence.").build();
    private PluginParameter<String> reversep = new PluginParameter.Builder<String>("reversep", null, String.class).guiName("Reverse Primer").required(true)
            .description("String containing the reverse primer sequence.").build();
    private PluginParameter<String> blastFile = new PluginParameter.Builder<String>("blastFile", null, String.class).guiName("Blast File").required(true).inFile()
            .description("Tab delimited Blast file output with NO header line, that contains only the data from columns for chrom, start postion, end position. \nThis data should be filtered to contain only entries whose identiy value was 98% or greater.").build();
    
    static GenomeSequence myRefSequence = null;
    // length of ref tag sequence from which to search for primer strings
    // This is half the length
    //static int refAlignLen = 600; // could be a plugin parameter - hard code for testing
    static int refAlignLen = 1000;

    public RampSeqAlignFromBlastTags() {
        super(null, false);
    }

    public RampSeqAlignFromBlastTags(Frame parentFrame) {
        super(parentFrame, false);
    }

    public RampSeqAlignFromBlastTags(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    @Override
    public void postProcessParameters() {

        if (myDBFile.isEmpty() || !Files.exists(Paths.get(dBFile()))) {
            throw new IllegalArgumentException("RepGenPhase2AlignerPlugin: postProcessParameters: Input DB not set or found");
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
            System.out.println("RampSeqAlignFromBlastTags:processData begin"); 
            RepGenDataWriter repGenData=new RepGenSQLite(dBFile());   
            time = System.nanoTime();
            
           //  Create synchronized map for use in parallel streams 
           //Multimap<Tag,Position> refTagPositionMap = Multimaps.synchronizedMultimap(HashMultimap.<Tag,Position>create());
            Multimap<Tag,Position> refTagPositionMap = HashMultimap.<Tag,Position>create();
            boolean success = createRefTagsFromBlast(blastFile(), forwardp(), reversep(), refTagPositionMap);
            if (!success) {
                ((RepGenSQLite)repGenData).close();
                return null;
            }
            
            // add to db
            // add ref and mapping approach to db, then reference tags
            repGenData.addMappingApproach("SmithWaterman");
            repGenData.addReferenceGenome(refGenome());
            
            // This adds to refTag table and creates/adds entries to physicalMapPosition table
            repGenData.putRefTagMapping(refTagPositionMap, refGenome());
            
            // Get tags stored from RepGenLoadSeqToDB        
            Set<Tag> tagsToAlign = repGenData.getTags();
            List<Tag> tagList = new ArrayList<Tag>(tagsToAlign);
            
            // Create synchronized map for use in parallel streams int the calculate alignment methods
            Multimap<Tag,AlignmentInfo> tagAlignInfoMap = Multimaps.synchronizedMultimap(HashMultimap.<Tag,AlignmentInfo>create());
            Multimap<RefTagData,AlignmentInfo> refTagAlignInfoMap = Multimaps.synchronizedMultimap(HashMultimap.<RefTagData,AlignmentInfo>create());
            
            // First align the tags against each other
            // Output of this is stored in db table tagAlignments
            System.out.println("Calling calculateTagTagAlignment for tags without refs");
            
            calculateTagTagAlignment( tagList, tagAlignInfoMap,repGenData);
            
            // neither tag is reference, null and -1 for tag1 chrom/pos
            //repGenData.putTagAlignments(tagAlignInfoMap, false,false, null,-1,refGenome());
 
            System.out.println("Calling calculateTagRefTagALignment for tags WITH ref");
            Set<RefTagData> refTags = repGenData.getRefTags();
            List<RefTagData> refTagList = new ArrayList<RefTagData>(refTags);
            
            tagAlignInfoMap.clear(); // remove old data
            calculateTagRefTagAlignment(tagList,refTagList,tagAlignInfoMap,refGenome(), repGenData);
            
            // Add the tag-refTag info to the tagAlignments table.           
            // tag1=nonref, tag2=ref, null and -1 for tag1 .  Alignment will be
            // done twice:  once for tag/refTag-fwd, once for tag/refTag-reverse complement
            System.out.println("Number of tag-refTag alignments, includes aligning to ref fwd and reverse strands: " + tagAlignInfoMap.size() + ", store to db.");
            //repGenData.putTagAlignments(tagAlignInfoMap,false,true,null,-1,refGenome());
            
            System.out.println("Calling calculateRefRefAlignment");            
            calculateRefRefAlignment(refTagList,refTagPositionMap,refTagAlignInfoMap, repGenData);
            
            //repGenData.putRefRefAlignments(refTagAlignInfoMap, refGenome()); // CREATE putRefRefAlignments -  different parameters
            ((RepGenSQLite)repGenData).close();

            myLogger.info("Finished RepGenPhase2AlignerPlugin\n");
        } catch (Exception exc) {
            exc.printStackTrace();
        }
        
        
        System.out.println("Process took " + (System.nanoTime() - totalTime)/1e9 + " seconds.\n");
        return null; 
    }

    public boolean createRefTagsFromBlast(String chromStartEndFile, String forwardP, String reverseP, Multimap<Tag,Position> refTagPositionMap) {
        System.out.println("Begin createRefTagsFromBlast");
        
        BufferedReader rdChroms = Utils.getBufferedReader(chromStartEndFile);      
        // HOlds 1-based start/end chromosome
        Multimap<String,Tuple<Integer,Integer>> chromStartEndMap = HashMultimap.create();
        
        // read the BLAST file data into the map
        try {
            String line;
            
            // FIle has been filtered via awk to contain just the
            // chrom, source-start and source-end columns.  There are often 
            // duplicates for each tag as they can align with 98% identity
            // on multiple chromosomes.
            // First column is chrom, seconds is startpos, third is endpos
            int lineCount = 0;
            while ( (line = rdChroms.readLine()) != null ) {
                lineCount++;
                String tokens[] = line.split("\t");
                String chrom = tokens[0];
                Tuple<Integer,Integer> startEnd = new Tuple<Integer,Integer>(Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
                chromStartEndMap.put(chrom, startEnd);
            }
            System.out.println("Read " + lineCount + " lines from BLAST input file.");
        } catch (Exception exc) {
            exc.printStackTrace();
            return false;
        }
        
        byte[] forwardRCBytes = NucleotideAlignmentConstants.reverseComplementAlleleByteArray(forwardP.getBytes());
        String forwardPRC = new String(forwardRCBytes);
        byte[] reverseRCBytes = NucleotideAlignmentConstants.reverseComplementAlleleByteArray(reverseP.getBytes());
        String reversePRC = new String(reverseRCBytes);
        // we have the map from which to create ref tags.
        
        int count = 0;
        System.out.println("Finished reading Blast file, num of entries added to chromStartEndMap: " + chromStartEndMap.size());

        // number of entries in map is quite a bit lower than number of entries in file.
        // I awk'd and sorted - there are many duplicate tags from the fastq files, so
        // this is fine.
        for (int idx = 1; idx < 11; idx++){
            String chrom = Integer.toString(idx);
            System.out.println("chromStartENdMap entries for " + idx + ": " + chromStartEndMap.get(chrom).size());
        }
        
        // These values are for printing metrics
        int perfectScore = 0;
        int oneMM = 0;
        int twoMM = 0;
        int notFound = 0;
        
        int fstart = 0;
        int frcstart = 0;
        int rstart = 0;
        int rrcstart = 0;
        int ostart = 0;
        int refLen150more = 0;
        int refLen150 = 0;
        int refLen150less = 0;
        // end metric values
        
        for (Map.Entry entry : chromStartEndMap.entries()) {
            count++;
            
            String chromString = (String)entry.getKey();
            Chromosome chrom = new Chromosome(chromString);
            Tuple<Integer,Integer> startEnd = (Tuple<Integer, Integer>) entry.getValue();
            int chromSize = myRefSequence.chromosomeSize(chrom);
            int origSize = Math.abs(startEnd.x - startEnd.y); // need length of original string for calcuating best primer end below
            
            int refStartIdx = startEnd.x;
            int refEndIdx = startEnd.y;
            if (refStartIdx > refEndIdx) { // reverse strand
                int temp = refStartIdx;
                refStartIdx = refEndIdx;
                refEndIdx = temp;
            }
            refEndIdx = (refEndIdx + 200) <= chromSize ? refEndIdx + 200 : chromSize;
 
            // Look for the reverse primer
            // Grab reference that is 200 greater than the end so can look for primer end
            byte[] refTagBytes = myRefSequence.chromosomeSequence(chrom,refStartIdx, refEndIdx);
            // Convert to allele string
            String refTagString = NucleotideAlignmentConstants.nucleotideBytetoString(refTagBytes);
             
            // lcj this is debug
            if (refTagString.startsWith(forwardp())) fstart++;
            else if (refTagString.startsWith(reversep())) rstart++;
            else if (refTagString.startsWith(forwardPRC)) frcstart++;
            else if (refTagString.startsWith(reversePRC)) rrcstart++;
            else ostart++;
            // lcj - end debug 
            
            // Check that ref string is good.
            if (refTagString == null || refTagString.length() == 0 || refTagString.contains("N") || refTagString.contains("null")) {
                System.out.println(" storeRefTagPositions - refString is NULL");
            } else {
                // look for end primer sequence in the string.
                // We don't know what the start primer is, so must check for all
                Tuple<Integer,Integer> forwardPrimer = null; // holds score and refAlignStartPos
                Tuple<Integer,Integer> reversePrimer = null;
                Tuple<Integer,Integer> forwardRCPrimer = null;
                Tuple<Integer,Integer> reverseRCPrimer = null;

                forwardPrimer = computePrimerSW(forwardP,  refTagString, 2,-1,-1); // Returns score and refAlignStartPos
                reversePrimer = computePrimerSW(reverseP,  refTagString, 2,-1,-1);
                forwardRCPrimer = computePrimerSW(forwardPRC,  refTagString, 2,-1,-1);
                reverseRCPrimer = computePrimerSW(reversePRC,  refTagString, 2,-1,-1);

                // Want to truncate the refTagSTring at the end of the ending primer.
                // Find which primer has the best score, and is at least 100 past the
                // beginning of the line (so we aren't matching on the forward primer).  
                int bestEndPos = origSize;
                
                int bestScore = 33; // 40 is best, 37 is 1 off, 34 is 2 off, will allow no more than 2 mismatch or gap
                // by being more stringent, we get longer refStrings for matching when primers aren't found
                if (forwardPrimer.y > 100 && forwardPrimer.x > bestScore) {
                    bestEndPos = forwardPrimer.y ;
                    bestScore = forwardPrimer.x;
                }
                if (reversePrimer.y > 100 && reversePrimer.x > bestScore) {
                    if (reversePrimer.x > bestScore) {
                        bestEndPos = reversePrimer.y ; 
                        bestScore = reversePrimer.x;
                    }
                    if (reversePrimer.x == bestScore) {
                        // if scores are equal, take the one furthest out
                        if (reversePrimer.y > bestEndPos) {
                            bestEndPos = reversePrimer.y ; 
                        }
                    }
                }
                if (reverseRCPrimer.y > 100 && reverseRCPrimer.x > bestScore) {
                    if (reverseRCPrimer.x > bestScore) {
                        bestEndPos = reverseRCPrimer.y ; 
                        bestScore = reverseRCPrimer.x;
                    }
                    if (reverseRCPrimer.x == bestScore) {
                        // if scores are equal, take the one with the furthest position out
                        if (reverseRCPrimer.y > bestEndPos) {
                            bestEndPos = reverseRCPrimer.y ; 
                        }
                    }
                }
                if (forwardRCPrimer.y > 100 && forwardRCPrimer.x > bestScore) {
                    if (forwardRCPrimer.x > bestScore) {
                        bestEndPos = forwardRCPrimer.y ; 
                        bestScore = forwardRCPrimer.x;
                    }
                    if (forwardRCPrimer.x == bestScore) {
                        // if scores are equal, take the one with the farthest position out
                        if (forwardRCPrimer.y > bestEndPos) {
                            bestEndPos = forwardRCPrimer.y ; 
                        }
                    }
                }
                if (bestScore == 40) perfectScore++;
                else if (bestScore > 36) oneMM++;
                else if (bestScore > 33) twoMM++;
                else notFound++;
                // Add forwardP.length to bestEnd POS !!!!!  then substr the string, and make the tag
                // Assumption is primers are all the same length, so doesn't matter which primer length is added
                int endPos = forwardP.length() + bestEndPos;
                String finalRef = refTagString.substring(0,endPos); // bestEndPos plus length of primer
                if (finalRef.length() > 150) refLen150more++;
                else if (finalRef.length() == 150) {
                    refLen150++;
                }
                else refLen150less++;
                Tag refTag = TagBuilder.instance(finalRef).reference().build();
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
                    System.out.println("- null refTag created - skipping !!");
                } 
            }           
        } 

        System.out.println("Number of tags processed: " + count );
        System.out.println("Finsished creating tags, refTagPositionMap size: " 
                + refTagPositionMap.size() + ", number refTags with length greater than 150: " + 
                refLen150more + ", equal 150: " + refLen150 + ", less than 150: " + refLen150less);
        System.out.println("Number of distinct ref tags: " + refTagPositionMap.keySet().size());
        System.out.println("Start seqences: forwardP " + fstart + ", reverseP " + rstart + ", forwardPRC " 
                + frcstart + ", reversePRC " + rrcstart + ", other " + ostart);
        System.out.println("End primers perfectscore: " + perfectScore + ", oneMM: " + oneMM 
                + ", twoMM: " + twoMM + ", end primer not found: " + notFound);
        
        return true;
    }

    private Tuple<Integer,Integer> computePrimerSW(String primerSeq, String refSeq, 
            int match, int mismatch, int gap) {
        
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
    
    // this calculates tag against tag alignment from the tags in the DB
    // tag table (not the refTag table)
    private void calculateTagTagAlignment(List<Tag> tags, Multimap<Tag,AlignmentInfo> tagAlignInfoMap, RepGenDataWriter repGenData){
        long time = System.nanoTime();
        long totalTime = System.nanoTime();
        // For each tag on the tags list, run SW against it and store in tagTagAlignMap 
        int totalProcessedTags = 0;
        int tagCount = 0;
        for (int tidx=0; tidx < tags.size(); tidx++) {
            Tag tag1 = tags.get(tidx);
        
            totalProcessedTags++;
            tagCount++;
            IntStream.range(tidx+1, tags.size()).parallel().forEach(item -> {
                //tagTagSW(tag1, tags.get(item), tagAlignInfoMap);
                Tag tag2 = tags.get(item);
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

            });

            if (tagCount > 999) {
                System.out.println("FInished aligning " + totalProcessedTags + " tags, this set took " + (System.nanoTime() - time)/1e9 + " seconds, now load to db ..." );
                time = System.nanoTime();
                repGenData.putTagTagAlignments(tagAlignInfoMap);
                System.out.println("Loading DB took " + (System.nanoTime() - time)/1e9 + " seconds.\n");
                tagCount = 0;
                tagAlignInfoMap.clear();
                time = System.nanoTime();
            }
        }
        if (tagCount > 0 ) {
            System.out.println("Finished processing last alignments, load to DB");
            time = System.nanoTime();
            repGenData.putTagTagAlignments(tagAlignInfoMap);
            System.out.println("Loading DB took " + (System.nanoTime() - time)/1e9 + " seconds.\n");
            tagCount = 0;
            tagAlignInfoMap.clear(); // start fresh with next 1000
        }

        System.out.println("Number of tags: " + tags.size() + ", TotalTime for calculateTagTagAlignment was " + (System.nanoTime() - totalTime) / 1e9 + " seconds");
    }

    // This calculates alignment of each tag against each reference tag.
    // For tag to refTag, we DO want every tag on tag list against every
    // tag on refTag list.
    private void calculateTagRefTagAlignment(List<Tag> tags, List<RefTagData> refTagDataList,
            Multimap<Tag,AlignmentInfo> tagAlignInfoMap, String refGenome, RepGenDataWriter repGenData){
        long totalTime = System.nanoTime();
        long time = System.nanoTime();
        // For each tag on the tags list, run SW against it and all refTags and store in tagAlignInfoMap 
        int totalProcessedTags = 0;
        int tagCount = 0;
        for (int tidx = 0; tidx < tags.size(); tidx++) {
            totalProcessedTags++;
            tagCount++;
            Tag tag1 = tags.get(tidx);

            refTagDataList.parallelStream().forEach (rtd -> {
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
            }); 
            if (tagCount > 999) {
                System.out.println("FInished aligning " + totalProcessedTags + " tags, this set took " + (System.nanoTime() - time)/1e9 + " seconds, now load to db ..." );
                time = System.nanoTime();
                // tag1 not ref tag (false), tag2 is ref tag (true)
                repGenData.putTagRefTagAlignments(tagAlignInfoMap,refGenome());               
                System.out.println("Loading DB took " + (System.nanoTime() - time)/1e9 + " seconds.\n");
                tagCount = 0;
                tagAlignInfoMap.clear();
                time = System.nanoTime();
            }
        }
        if (tagCount > 0) {
            System.out.println("FInished last alignments " + totalProcessedTags + " tags, this set took " + (System.nanoTime() - time)/1e9 + " seconds, now load to db ..." );
            time = System.nanoTime();
            // tag1 not ref tag (false), tag2 is ref tag (true)
            repGenData.putTagRefTagAlignments(tagAlignInfoMap,refGenome());               
            System.out.println("Loading DB took " + (System.nanoTime() - time)/1e9 + " seconds.\n");
            tagCount = 0;
            tagAlignInfoMap.clear();
            time = System.nanoTime();
        }
        System.out.println("Num tags: " + tags.size() + ", Num refTags: " + refTagDataList.size() + ", TotalTime for calculateTagRefTagAlignment was " + (System.nanoTime() - totalTime) / 1e9 + " seconds");
    }
    
    // This calculates alignment of  reference tags against each other.
    // Only need to caluclate refTag to other refTags on list below it - not
    // against a ref tag that has already calculated the alignment against the current.
    private void calculateRefRefAlignment(List<RefTagData> refTags, Multimap<Tag,Position> refTagPosMap,
            Multimap<RefTagData,AlignmentInfo> refTagAlignInfoMap, RepGenDataWriter repGenData){
        long totalTime = System.nanoTime();
        long time = System.nanoTime();
        int tagCount = 0;
        int totalProcessedTags = 0;
        // For each tag on the reftags list, run SW against all other tags in the list
        for (int tidx = 0; tidx < refTags.size(); tidx++) {
            tagCount++;
            totalProcessedTags++;
            RefTagData tag1 = refTags.get(tidx);
            final int tIdxFinal = tidx;
            IntStream.range(tidx+1, refTags.size()).parallel().forEach(item -> {
                RefTagData tag2 = refTags.get(tIdxFinal+1);
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
            });
            if (tagCount > 999) {
                System.out.println("FInished aligning " + totalProcessedTags + " tags, this set took " + (System.nanoTime() - time)/1e9 + " seconds, now load to db ..." );
                time = System.nanoTime();
                // tag1 not ref tag (false), tag2 is ref tag (true)
                repGenData.putRefRefAlignments(refTagAlignInfoMap, refGenome());
                               
                System.out.println("Loading DB took " + (System.nanoTime() - time)/1e9 + " seconds.\n");
                tagCount = 0;
                refTagAlignInfoMap.clear();
                time = System.nanoTime();
            }
        }
        if (tagCount > 0) {
            System.out.println("FInished last refRef alignments " + totalProcessedTags + " tags, this set took " + (System.nanoTime() - time)/1e9 + " seconds, now load to db ..." );
            time = System.nanoTime();
            // tag1 not ref tag (false), tag2 is ref tag (true)
            repGenData.putRefRefAlignments(refTagAlignInfoMap, refGenome());
                           
            System.out.println("Loading DB took " + (System.nanoTime() - time)/1e9 + " seconds.\n");
            tagCount = 0;
            refTagAlignInfoMap.clear();
            time = System.nanoTime();
        }
        System.out.println("Number of refTags: " + refTags.size() + ", TotalTime for calculateREfRefAlignment was " + (System.nanoTime() - totalTime) / 1e9 + " seconds");
    }
    
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
//     public static void main(String[] args) {
//         GeneratePluginCode.generate(RampSeqAlignFromBlastTags.class);
//     }
     
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

    /**
     * Input database file with tags and taxa distribution
     *
     * @return Input DB
     */
    public String dBFile() {
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
    public RampSeqAlignFromBlastTags dBFile(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }

    /**
     * Referemce Genome File for aligning against 
     *
     * @return Reference Genome File
     */
    public String refGenome() {
        return refGenome.value();
    }

    /**
     * Set Reference Genome File. Referemce Genome File for
     * aligning against 
     *
     * @param value Reference Genome File
     *
     * @return this plugin
     */
    public RampSeqAlignFromBlastTags refGenome(String value) {
        refGenome = new PluginParameter<>(refGenome, value);
        return this;
    }

    /**
     * Minimum count of reads for a tag to be aligned
     *
     * @return Min Tag Count
     */
    public Integer minTagCount() {
        return minTagCount.value();
    }

    /**
     * Set Min Tag Count. Minimum count of reads for a tag
     * to be aligned
     *
     * @param value Min Tag Count
     *
     * @return this plugin
     */
    public RampSeqAlignFromBlastTags minTagCount(Integer value) {
        minTagCount = new PluginParameter<>(minTagCount, value);
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
    public RampSeqAlignFromBlastTags match_reward(Integer value) {
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
    public RampSeqAlignFromBlastTags mismatch_penalty(Integer value) {
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
    public RampSeqAlignFromBlastTags gap_penalty(Integer value) {
        gap_penalty = new PluginParameter<>(gap_penalty, value);
        return this;
    }

    /**
     * String containing the forward primer sequence.
     *
     * @return Forward Primer
     */
    public String forwardp() {
        return forwardp.value();
    }

    /**
     * Set Forward Primer. String containing the forward primer
     * sequence.
     *
     * @param value Forward Primer
     *
     * @return this plugin
     */
    public RampSeqAlignFromBlastTags forwardp(String value) {
        forwardp = new PluginParameter<>(forwardp, value);
        return this;
    }

    /**
     * String containing the reverse primer sequence.
     *
     * @return Reverse Primer
     */
    public String reversep() {
        return reversep.value();
    }

    /**
     * Set Reverse Primer. String containing the reverse primer
     * sequence.
     *
     * @param value Reverse Primer
     *
     * @return this plugin
     */
    public RampSeqAlignFromBlastTags reversep(String value) {
        reversep = new PluginParameter<>(reversep, value);
        return this;
    }

    /**
     * Tab delimited Blast file output with NO header line,
     * that contains only the data from columns for chrom,
     * start postion, end position. 
     * This data should be filtered to contain only entries
     * whose identiy value was 98% or greater.
     *
     * @return Blast File
     */
    public String blastFile() {
        return blastFile.value();
    }

    /**
     * Set Blast File. Tab delimited Blast file output with
     * NO header line, that contains only the data from columns
     * for chrom, start position, end position. 
     * This data should be filtered to contain only entries
     * whose identify value was 98% or greater.
     *
     * @param value Blast File
     *
     * @return this plugin
     */
    public RampSeqAlignFromBlastTags blastFile(String value) {
        blastFile = new PluginParameter<>(blastFile, value);
        return this;
    }
}
