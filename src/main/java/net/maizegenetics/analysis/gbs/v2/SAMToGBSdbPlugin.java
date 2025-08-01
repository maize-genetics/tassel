/*
 * SAMConverterPlugin
 */
package net.maizegenetics.analysis.gbs.v2;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.tag.*;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;

import java.awt.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Optional;
import java.util.Set;

/**
 * Reads SAM file formats to determine the potential positions of Tags against the reference genome.
 *
 * @author Ed Buckler
 *
 */
public final class SAMToGBSdbPlugin extends AbstractPlugin {

    boolean cleanCutSites = true;
    private static final Logger myLogger = LogManager.getLogger(SAMToGBSdbPlugin.class);

    private PluginParameter<String> myInputFile = new PluginParameter.Builder<String>("i", null, String.class).guiName("SAM Input File").required(true).inFile()
            .description("Name of input file in SAM text format").build();
    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("GBS DB File").required(true).outFile()
            .description("Name of output file (e.g. GBSv2.db)").build();
    private PluginParameter<Double> alignProportion = new PluginParameter.Builder<Double>("aProp", 0.0, Double.class).guiName("SAM Min Align Proportion").required(false)
            .range(Range.closed(0.0, 1.0) ).description("Minimum proportion of sequence that must align to store the SAM entry").build();
    private PluginParameter<Integer> minAlignLength = new PluginParameter.Builder<Integer>("aLen", 0, Integer.class).guiName("SAM Min Align Length").required(false)
            .range(Range.closed(0, 1000) ).description("Minimum length of bps aligning to store the SAM entry").build();
    private PluginParameter<String> mappingApproach = new PluginParameter.Builder<String>("mapper", "BWA", String.class).guiName("Mapper").required(false)
            .description("Mapping approach (one of Bowtie2, BWA, or bwaMem)").build();
    private PluginParameter<Boolean> myDeleteOldData = new PluginParameter.Builder<Boolean>("deleteOldData",true,Boolean.class).guiName("Delete Old Data")
            .description("Delete existing SNP data from db tables").build();
    private PluginParameter<Integer> minMAPQ = new PluginParameter.Builder<Integer>("minMAPQ", 0, Integer.class).guiName("SAM Min MAPQ value").required(false)
            .description("Minimum value of MAPQ to store the SAM entry").build();
    
    private enum tagPresence {
        present,
        originalPresent,
        notPresent;       
    }

    private static int mapQBad = 0;
    private static int mapQGood = 0;
    public SAMToGBSdbPlugin() {
        super(null, false);
    }

    public SAMToGBSdbPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        int tagsNotFoundInDB=0, tagsNotMapped=0;       
        try {
            BufferedReader bw=Utils.getBufferedReader(sAMInputFile());
            TagDataWriter tagData=new TagDataSQLite(gBSDBFile());
            
            if (deleteOldData()) {
                myLogger.info("deleteOldData is TRUE: Clearing existing Alignment, Discovery and SNPQuality data");
                tagData.clearSNPQualityData();
                tagData.clearDiscoveryData();
                tagData.clearAlignmentData();
            }
            Set<Tag> knownTags=tagData.getTags();
            Multimap<Tag,Position> tagPositions= HashMultimap.create(knownTags.size(),2);
            String inputLine;
            while((inputLine=bw.readLine())!=null) {
                if(inputLine.startsWith("@")) {
                    // this is header - check if it is bowtie
                    if (inputLine.contains("bowtie2")) mappingApproach("Bowtie2");
                    continue;
                }
                Tuple<Tag,Optional<Position>> tagPositionTuple=parseRow(inputLine);
                if (tagPositionTuple == null) continue; 
                // Must use method "isKnownTag" to handle BWA-Mem issue
                // from Jesse Poland (TAS-722)   
                tagPresence tagP = isKnownTag(inputLine, tagPositionTuple.x, knownTags);
                if (tagP == tagPresence.notPresent) {
                    tagsNotFoundInDB++;
                }           
                if(tagPositionTuple.y.isPresent()) {
                    if (tagP == tagPresence.originalPresent) {                       
                        String[] stringTokens=inputLine.split("\\s");
                        String origSeq = stringTokens[0].split("=")[1]; //stringTokens[0] is tagSeg=<original sequence here>
                        Tag oTag = TagBuilder.instance(origSeq).build();
                        tagPositions.put(oTag, tagPositionTuple.y.get());
                    } else {
                        tagPositions.put(tagPositionTuple.x,tagPositionTuple.y.get());
                    }                             
                } else {
                    tagsNotMapped++;
                }
            }
            bw.close();

            if(tagsNotFoundInDB==0) {tagData.putTagAlignments(tagPositions);
                myLogger.info("Finished reading SAM file and adding tags to DB."
                    + "\nTotal number of tags mapped: " + tagPositions.keySet().size() + " (total mappings " + tagPositions.size() + ")"
                        + "\nTags not mapped: " + tagsNotMapped  
                        + "\nTags dropped due to minimum mapq value: " + mapQBad + "\n\n");}
            else {
                System.out.println("Unobserved tags were found in the SAM file count= " + tagsNotFoundInDB);
                myLogger.info("Finished reading SAM file.  No Tags added to DB as "+tagsNotFoundInDB+" unobserved tags were found.\n" +
                        "Please ensure all tags in the SAM file already exist in the DB.\n\n");
            }
            ((TagDataSQLite)tagData).close();  //todo autocloseable should do this but it is not working.


        } catch (Exception e) {
            myLogger.info("Catch in reading TagCount file e=" + e);
            e.printStackTrace();
        }
        return null;
    }


    //should this be converted to a stream?
    private Tuple<Tag,Optional<Position>> parseRow(String inputLine) {
        final int name = 0, flag = 1, chr = 2, pos = 3, cigar = 5, tagS = 9; // column indices in inputLine
        String[] s=inputLine.split("\\s+");
        // Use the tag we stored in the header field.  This handles forward vs
        // reverse strand (no longer need to reverse complement to our original)
        // and it contains seqments that were "clipped" by the aligner.
        
        String origSeq = s[0].split("=")[1]; //s[0] is tagSeg=<original sequence here
        Tag tag = TagBuilder.instance(origSeq).build();        
        int mapq = Integer.parseInt(s[4]);
        if (mapq < minMAPQ()) {
            // Dropping anything with quality less than user specified mapQ. 
            // Sites  discusing this:
            // http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
            // http://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42
            mapQBad++;
            return new Tuple<>(tag,Optional.<Position>empty());
        } else mapQGood++; 
        
        // A tag consisting of 32 T's become -1 in "getLongFromSequence", which results in a "null" tag
        // This was seen in the Zea_mays.AGPv3 chromosome files.  Shouldn't happen since we're now using
        // the sequence from the header line
        if (tag == null) return null;
        // The two lines need to be here to make sure the sequence can be found in the DB
        boolean forwardStrand=isForwardStrand(s[flag]);
        int samPos = Integer.parseInt(s[pos]); 
        
        int[] alignSpan = SAMUtils.adjustCoordinates(s[cigar], samPos);        
        int cutPos;
        if (forwardStrand == true) {
            cutPos = alignSpan[0];
        } else  {
            cutPos = alignSpan[1];
        }
        // adjustCpoordinates()  considers indels and clipping when re-calculating where this
        // sequence's first bp matches to the refence genome.  The result should not be negative.
        // It can become negative when clipping moves the alignment out of range of the
        // reference.  For example:  the .sam file says the leftmost mapping position of 
        // the first matched bp is 1 (ie "pos" = 1), but the cigar is 32S32M (32 soft clipped,
        // then 32 match).  Drop these alignments.
        
        // These alignments should be dropped whenever alignSpan[0] is < 0.
        // FOr a forward strand, that is the start and it should be dropped as agreed to
        // above.  For a reverse strand, when we hit DiscoverySNPCallerPluginV2, we will translate
        // the reverse strands into a forward strand.  The strand end becomes the strand start, and
        // we can't have a negative start position.
        if (alignSpan[0] < 0 || alignSpan[1] < 0){
            return new Tuple<>(tag,Optional.<Position>empty());
        }            

        if (!hasAlignment(s[flag])) return new Tuple<>(tag,Optional.<Position>empty());
        // Check for minimum alignment length and proportion
        if (!hasMinAlignLength(s)) return new Tuple<> (tag,Optional.<Position>empty());
        if (!hasMinAlignProportion(s)) return new Tuple<> (tag,Optional.<Position>empty());
        Chromosome chromosome = new Chromosome(s[chr]); // Chromosome class parses the chromosome
        String alignmentScore=getAlignmentScore(s);        
        // TASSEL defines forward/reverse as the following:
        // public interface Position extends Comparable<Position> {

        // public static final byte STRAND_PLUS = (byte) 1;
        // public static final byte STRAND_MINUS = (byte) 0;
        // public static final byte STRAND_UNKNOWN = Byte.MIN_VALUE;
        
        // The aligner determined if
        // the tag sequence from the read was forward or backwards. It marks it as
        // such in the strand bit of the flag.  If it was reverse, the aligner reverse
        // compliments the sequence.  We check that above and if the aligner reversed it,
        // we reverse it back so that it matches what is stored in our DB.  For this
        // position, store whether it was forward.
        byte strand = (byte)(forwardStrand ? 1 : 0);
        Position position=new GeneralPosition
                .Builder(chromosome,cutPos)
                .strand(strand)
                //.strand((byte)1)
                .addAnno("forward", forwardStrand?"true":"false")
                .addAnno("mappingapproach", mappingApproach())
                .addAnno("cigar", s[cigar])
                .addAnno("supportvalue", alignmentScore)  //todo include again
                .build();
        return new Tuple<>(tag,Optional.of(position));
    }

    /**Bit 3 of the flag (16) indicates whether forward strand*/
    private boolean hasAlignment(String samFlag){
        int flag=Integer.parseInt(samFlag);
        return ((flag & 4) == 0);
    }

    /**Bit 5 of the flag (16) indicates whether forward strand*/
    private boolean isForwardStrand(String samFlag){
        int flag=Integer.parseInt(samFlag);
        return ((flag & 16) == 0);
    }
    
    /**Check optional Tags for CIGAR and MD to verify min align length*/
    private boolean hasMinAlignLength(String[] samReadParsed){
    	// Optional user tags may or may not be present, and may be
    	// present in different order in different agligner outputs
    	// Tags start at position 11 (0 based), so start looking from here   	
    	if (minAlignLength() == 0) return true; // 0 is default - no minimum length
    	int matchLen = calculateNumberAligned(samReadParsed);
    	if (matchLen >= minAlignLength()) return true;
    	else return false;
    }
    
    /**Check optional Tags for CIGAR and MD to verify min align proportion*/
    private boolean hasMinAlignProportion(String[] samRead){  	
    	if (minAlignProportion() == 0) return true; // 0 is default - no minimum proportion 
    	float seqLength = samRead[9].length(); // sequence is in 9th position of 0 based array
    	int matchLen = calculateNumberAligned(samRead);
    	float matchProportion = (float)matchLen/seqLength;
    	if (matchProportion >= minAlignProportion()) return true;
    	else return false;
    }

    private int calculateNumberAligned(String[] samRead){   	
    	// Look for MD in the SAM output. Optional fields start at field 12 (11 in 0 based array)
    	String mdField = "";
    	for (int index=11; index < samRead.length; index++) {
    		String[] optField = samRead[index].split(":");
    		if (optField[0].equals("MD")) {
    			mdField = samRead[index];
    			break;
    		}
    	}
    	
		int matchLen = 0;
    	// Calculate minimum alignment length.  MD field is first choice, then CIGAR
    	// MD field specifies for sequence match 
    	if (!mdField.equals("")) {
    		String mdVal = mdField.split(":")[2]; // we want 3rd value in MD:Z:3T5^AG6
    		int curNum = 0;
    		for (int mdIdx = 0; mdIdx < mdVal.length(); mdIdx++) {
    			char currChar = mdVal.charAt(mdIdx);
    			if (Character.isDigit(currChar)) {
    				curNum = curNum *10 + Character.getNumericValue(currChar); // convert byte to unsigned int
    			} else {
    				matchLen += curNum;  // DO we need to consider num chars that don't match?
    				curNum = 0;
    			}
    		} 
    		matchLen += curNum; // takes care of case where MD is a single number string
    	} else { // use CIGAR value gives alignment match, not sequence match
        	String cigar = samRead[5];
    		int curNum = 0;
    		for (int cIdx = 0; cIdx < cigar.length(); cIdx++) {
    			char currChar = cigar.charAt(cIdx);
    			if (Character.isDigit(currChar)) {
    				curNum = curNum *10 + Character.getNumericValue(currChar); // convert byte to int, add to total
    			} else {
    				if (currChar == 'M' || currChar == 'm') {
    					// M indicates a match, should be proceeded by a number
    					// Don't count gaps indicated by 'S' or 'H' 
    	   				matchLen += curNum;
    				} 
    				curNum = 0;
    			}
    		} 
    		matchLen += curNum; // takes care of case where CIGAR is a single number string
    	}
    	return matchLen;
    }

    private String getAlignmentScore(String[] samRead) {
    	String asField = null;
    	for (int index=11; index < samRead.length; index++) {
    		String[] optField = samRead[index].split(":");
    		if (optField[0].equals("AS")) {
    			asField = samRead[index].split(":")[2];
    			break;
    		}
    	}
    	if (asField == null) {
    	    // Too many warning messages.  If the AS field is absent it is most probably absent
    	    // for all the entries in the file.
    		//myLogger.info("SAMToGBSDbPluginV2: warning: alignmentScore not present in Sam File, defaulting to 0");
    		asField = "0";
    	} 
    	return asField;
    }
    
    private tagPresence isKnownTag(String inputLine, Tag tag, Set knownTags){
        // 1.  Check if tag made from the aligner's sequence occurs in the db, if yes, return "present" 
        // 2.  Check if tag made from the original sequence occurs in the db, if yes, return "originalPresent"
        // 3.  If neither sequence can be found, return "notPresent"
        if (knownTags.contains(tag)) {
            return tagPresence.present; // good - no processing needed
        }

        String[] stringTokens=inputLine.split("\\s");
        String origSeq = stringTokens[0].split("=")[1];

        Tag oTag = TagBuilder.instance(origSeq).build();
        if (knownTags.contains(oTag)) {
            // The tag created from the aligner's sequence does not appear in the database,
            // However the original tag sequence DOES appear in the db, so store the position
            // against this tag rather than the tag created from the aligner sequence.  When the aligner
            // performs "hard-clipping" the "clipped" portions of the tag are removed. BWA-MEM does this.
            //System.out.println("LCJ - SAMToGBSDb:isKnownTag - CHANGE THE TAG !!");
            return tagPresence.originalPresent;
        }
        return tagPresence.notPresent;
    }

    /**
     * Reads SAM files output from BWA or Bowtie2
     */
//    private void readSAMFile(String inputFileName, int tagLengthInLong) {
//        System.out.println("Reading SAM format tag alignment from: " + inputFileName);
//        this.tagLengthInLong = tagLengthInLong;
//        String inputStr = "Nothing has been read from the file yet";
//        int nHeaderLines = countTagsInSAMfile(inputFileName); // detects if the file is Bowtie2, initializes topm matrices
//        int tagIndex = Integer.MIN_VALUE;
//        try {
//            BufferedReader br;
//            if (inputFileName.endsWith(".gz")) {
//                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(inputFileName)))));
//            } else {
//                br = new BufferedReader(new FileReader(new File(inputFileName)), 65536);
//            }
//            for (int i = 0; i < nHeaderLines; i++) {
//                br.readLine();
//            } // Skip over the header
//            for (tagIndex = 0; tagIndex < myNumTags; tagIndex++) {
//                inputStr = br.readLine();
//                parseSAMAlignment(inputStr, tagIndex);
//                if (tagIndex % 1000000 == 0) {
//                    System.out.println("Read " + tagIndex + " tags.");
//                }
//            }
//            br.close();
//        } catch (Exception e) {
//            System.out.println("\n\nCatch in reading SAM alignment file at tag " + tagIndex + ":\n\t" + inputStr + "\nError: " + e + "\n\n");
//            e.printStackTrace();
//            System.exit(1);
//        }
//    }
//
//    private int countTagsInSAMfile(String inputFileName) {
//        mySAMFormat = SAMFormat.BWA;  // format is BWA by default
//        myNumTags = 0;
//        int nHeaderLines = 0;
//        String currLine = null;
//        try {
//            String[] inputLine;
//            ArrayList<String> chrNames = new ArrayList<String>();
//            BufferedReader br;
//            if (inputFileName.endsWith(".gz")) {
//                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(inputFileName)))));
//            } else {
//                br = new BufferedReader(new FileReader(new File(inputFileName)), 65536);
//            }
//            while ((currLine = br.readLine()) != null) {
//                inputLine = currLine.split("\\s");
//                if (inputLine[0].contains("@")) {
//                    //SAM files produced by Bowtie2 contain the string "@PG     ID:bowtie2      PN:bowtie2 "
//                    if (inputLine[1].contains("bowtie2")) {
//                        mySAMFormat = SAMFormat.BOWTIE2;
//                    }
//                    nHeaderLines++;
//                } else {
//                    String chr = inputLine[2];
//                    if (!chrNames.contains(chr)) {
//                        chrNames.add(chr);
//                    }
//                    myNumTags++;
//                    if (myNumTags % 1000000 == 0) {
//                        System.out.println("Counted " + myNumTags + " tags.");
//                    }
//                }
//            }
//            br.close();
//            System.out.println("Found " + myNumTags + " tags in SAM file.  Assuming " + mySAMFormat + " file format.");
//        } catch (Exception e) {
//            System.out.println("Catch in counting lines of alignment file at line " + currLine + ": " + e);
//            e.printStackTrace();
//            System.exit(1);
//        }
//        initMatrices(myNumTags);
//        return nHeaderLines;
//    }
//
//
//    private void parseSAMAlignment(String inputStr, int tagIndex) {
//        String[] inputLine = inputStr.split("\t");
//        int name = 0, flag = 1, chr = 2, pos = 3, cigar = 5, tagS = 9; // column indices in inputLine
//        String nullS = this.getNullTag();
//        byte currStrand = ((Integer.parseInt(inputLine[flag]) & 16) == 16) ? (byte) -1 : (byte) 1; // bit 0x10 (= 2^4 = 16) is set: REVERSE COMPLEMENTED
//        if ((Integer.parseInt(inputLine[flag]) & 4) == 4) {  // bit 0x4 (= 2^2 = 4) is set: NO ALIGNMENT
//            recordLackOfSAMAlign(tagIndex, inputLine[tagS], inputLine[name], nullS, currStrand);
//        } else {  // aligns to one or more positions
//            HashMap<String, Integer> SAMFields = parseOptionalFieldsFromSAMAlignment(inputLine);
//            byte bestHits = (byte) Math.min(SAMFields.get("nBestHits"), Byte.MAX_VALUE);
//            byte editDist = (byte) Math.min(SAMFields.get("editDist"), Byte.MAX_VALUE);
//            recordSAMAlign(
//                    tagIndex,
//                    inputLine[tagS],
//                    inputLine[name],
//                    nullS,
//                    bestHits,
//                    inputLine[chr],
//                    currStrand,
//                    Integer.parseInt(inputLine[pos]),
//                    inputLine[cigar],
//                    editDist);
//        }
//    }
//
//    private HashMap<String, Integer> parseOptionalFieldsFromSAMAlignment(String[] inputLine) {
//        HashMap<String, Integer> SAMFields = new HashMap<String, Integer>();
//        if (mySAMFormat == SAMFormat.BWA) {
//            for (int field = 11; field < inputLine.length; field++) { // Loop through all the optional field of the SAM alignment
//                if (inputLine[field].regionMatches(0, "X0", 0, 2)) {        // X0 = SAM format for # of "high-quality" alignments of this query.  Specific to BWA.
//                    SAMFields.put("nBestHits", Integer.parseInt(inputLine[field].split(":")[2]));
//                } else if (inputLine[field].regionMatches(0, "NM", 0, 2)) { // NM = SAM format for edit distance to the reference.  Common to BWA and Bowtie2.
//                    SAMFields.put("editDist", Integer.parseInt(inputLine[field].split(":")[2]));
//                }
//            }
//        } else {  // bowtie2 -M format
//            for (int field = 11; field < inputLine.length; field++) { // Loop through all the optional field of the SAM alignment
//                if (inputLine[field].regionMatches(0, "AS", 0, 2)) {        // AS = SAM format for alignment score of the best alignment.  Specific to bowtie2.
//                    SAMFields.put("bestScore", Integer.parseInt(inputLine[field].split(":")[2]));
//                } else if (inputLine[field].regionMatches(0, "XS", 0, 2)) { // XS = SAM format for alignment score of 2nd best alignment.  Specific to bowtie2.
//                    SAMFields.put("nextScore", Integer.parseInt(inputLine[field].split(":")[2]));
//                } else if (inputLine[field].regionMatches(0, "NM", 0, 2)) { // NM = SAM format for edit distance to the reference.  Common to BWA and Bowtie2.
//                    SAMFields.put("editDist", Integer.parseInt(inputLine[field].split(":")[2]));
//                }
//            }
//            if (SAMFields.containsKey("bestScore")) {
//                if (SAMFields.containsKey("nextScore")) {
//                    if (SAMFields.get("bestScore") > SAMFields.get("nextScore")) {
//                        SAMFields.put("nBestHits", 1);
//                    } else {
//                        SAMFields.put("nBestHits", 99);  // 99 will stand for an unknown # of multiple hits
//                    }
//                } else {
//                    SAMFields.put("nBestHits", 1);
//                }
//            }
//        }
//        return SAMFields;
//    }
//
//    private void writeLogFile(TagsOnPhysicalMap topm) {
//        try {
//            DataOutputStream report = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile() + ".log"), 65536));
//            int[] aligned = topm.mappedTags();
//            int unique = 0, multi = 1;  // the indices of aligned
//            int unaligned = topm.getTagCount() - aligned[unique] - aligned[multi];
//            report.writeBytes(
//                    "Input file: " + inputFile() + "\n"
//                    + "Output file: " + outputFile() + "\n"
//                    + "Total " + topm.getTagCount() + " tags\n\t"
//                    + aligned[unique] + " were aligned to unique postions\n\t"
//                    + aligned[multi] + " were aligned to multiple postions\n\t"
//                    + unaligned + " could not be aligned.\n\n");
//            int[] dist = topm.mappingDistribution();
//            report.writeBytes("nPositions  nTags\n");
//            for (int i = 0; i < dist.length; i++) {
//                if (dist[i] > 0) {
//                    if (i < 10) {
//                        report.writeBytes(i + "           " + dist[i] + "\n");
//                    } else if (i < 100) {
//                        report.writeBytes(i + "          " + dist[i] + "\n");
//                    } else if (i < 1000) {
//                        report.writeBytes(i + "         " + dist[i] + "\n");
//                    }
//                }
//            }
//            report.close();
//        } catch (Exception e) {
//            myLogger.warn("Caught exception while writing log file: " + e);
//        }
//    }

    /**
     * Set Mapper (one of 'Bowtie2', 'BWA', or 'bwaMem').
     *
     * @param value Mapper type
     *
     * @return this plugin
     */
    public SAMToGBSdbPlugin mappingApproach(String value) {
        mappingApproach = new PluginParameter<>(mappingApproach, value);
        return this;
    }

    /**
     * Get Mapper (one of 'Bowtie2', 'BWA', or 'bwaMem').
     *
     * @param value Mapper type
     *
     * @return String for mapper type
     */
    public String mappingApproach() {
        return mappingApproach.value();
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(SAMToGBSdbPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public TagData runPlugin(DataSet input) {
        return (TagData) performFunction(input).getData(0).getData();
    }

    /**
     * Name of input file in SAM text format
     *
     * @return SAM Input File
     */
    public String sAMInputFile() {
        return myInputFile.value();
    }

    /**
     * Set SAM Input File. Name of input file in SAM text
     * format
     *
     * @param value SAM Input File
     *
     * @return this plugin
     */
    public SAMToGBSdbPlugin sAMInputFile(String value) {
        myInputFile = new PluginParameter<>(myInputFile, value);
        return this;
    }

    /**
     * Name of output file (e.g. GBSv2.db)
     *
     * @return GBS DB File
     */
    public String gBSDBFile() {
        return myOutputFile.value();
    }

    /**
     * Set GBS DB File. Name of output file (e.g. GBSv2.db)
     *
     * @param value GBS DB File
     *
     * @return this plugin
     */
    public SAMToGBSdbPlugin gBSDBFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
        return this;
    }

    /**
     * Minimum proportion of sequence that must align to store
     * the SAM entry
     *
     * @return SAM Min Align Proportion
     */
    public Double minAlignProportion() {
        return alignProportion.value();
    }

    /**
     * Set SAM Min Align Proportion. Minimum proportion of
     * sequence that must align to store the SAM entry
     *
     * @param value SAM Min Align Proportion
     *
     * @return this plugin
     */
    public SAMToGBSdbPlugin minAlignProportion(Double value) {
        alignProportion = new PluginParameter<>(alignProportion, value);
        return this;
    }

    /**
     * Minimum length of bps aligning to store the SAM entry
     *
     * @return SAM Min Align Length
     */
    public Integer minAlignLength() {
        return minAlignLength.value();
    }

    /**
     * Set SAM Min Align Length. Minimum length of bps aligning
     * to store the SAM entry
     *
     * @param value SAM Min Align Length
     *
     * @return this plugin
     */
    public SAMToGBSdbPlugin minAlignLength(Integer value) {
        minAlignLength = new PluginParameter<>(minAlignLength, value);
        return this;
    }
    /**
     * Minimum value of MAPQ to store the SAM entry
     *
     * @return SAM minimum MAPQ value
     */
    public Integer minMAPQ() {
        return minMAPQ.value();
    }

    /**
     * Set SAM minimum MAPQ value. Minimum value of MAPQ
     * to store the SAM entry
     *
     * @param value SAM minimum MAPQ value
     *
     * @return this plugin
     */
    public SAMToGBSdbPlugin minMAPQ(Integer value) {
        minMAPQ = new PluginParameter<>(minMAPQ, value);
        return this;
    }

    /**
     * Delete exisiting Alignment data from DB
     *
     * @return deleteOldData
     */
    public Boolean deleteOldData() {
        return myDeleteOldData.value();
    }

    /**
     * Set Delete old data flag.  True indicates we want the
     * db tables cleared
     *
     * @param value true/false - whether to delete data
     *
     * @return this plugin
     */
    public SAMToGBSdbPlugin deleteOldData(Boolean value) {
        myDeleteOldData = new PluginParameter<>(myDeleteOldData, value);
        return this;
    }
    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "SAM to TOPM Converter";
    }

    @Override
    public String getToolTipText() {
        return "SAM to TOPM Converter";
    }
}
