/*
 * DiscoverySNPCallerPluginV2
 */
package net.maizegenetics.analysis.gbs.v2;

import static net.maizegenetics.dna.snp.NucleotideAlignmentConstants.GAP_ALLELE;
import static net.maizegenetics.dna.snp.NucleotideAlignmentConstants.GAP_ALLELE_CHAR;
import static net.maizegenetics.dna.snp.NucleotideAlignmentConstants.getNucleotideAlleleByte;

import java.awt.Frame;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.GenomeSequence;
import net.maizegenetics.dna.map.GenomeSequenceBuilder;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.SimpleAllele;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TagDataWriter;
import net.maizegenetics.dna.tag.TaxaDistBuilder;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Tuple;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.util.ConcurrencyTools;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;

/**
 * This class aligns tags at the same physical location against one another,
 * calls SNPs, and then outputs the SNPs to a HapMap file.
 *
 * It is multi-threaded, as there are substantial speed increases with it.
 *
 * @author Ed Buckler
 * @author Jeff Glaubitz
 */
public class DiscoverySNPCallerPluginV2 extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(DiscoverySNPCallerPluginV2.class);

    private PluginParameter<String> myInputDB = new PluginParameter.Builder<>("db", null, String.class).guiName("Input GBS Database").required(true).inFile()
            .description("Input Database file if using SQLite").build();
    private PluginParameter<Double> myMinMinorAlleleFreq = new PluginParameter.Builder<>("mnMAF", 0.01, Double.class).guiName("Min Minor Allele Freq")
            .description("Minimum minor allele frequency").build();
    private PluginParameter<Double> myMinLocusCoverage = new PluginParameter.Builder<>("mnLCov", 0.1, Double.class).guiName("Min Locus Coverage")
            .description("Minimum locus coverage (proportion of Taxa with a genotype)").build();
    private PluginParameter<String> myRefGenome = new PluginParameter.Builder<>("ref", null, String.class).guiName("Reference Genome File").inFile()
            .description("Path to reference genome in fasta format. Ensures that a tag from the reference genome is always included "
                    + "when the tags at a locus are aligned against each other to call SNPs. The reference allele for each site "
                    + "is then provided in the output HapMap files, under the taxon name \"REFERENCE_GENOME\" (first taxon). "
                    + "DEFAULT: Don't use reference genome.").build();
    private PluginParameter<Chromosome> myStartChr = new PluginParameter.Builder<>("sC", null, Chromosome.class).guiName("Start Chromosome").required(false)
            .description("Start Chromosome").build();
    private PluginParameter<Chromosome> myEndChr = new PluginParameter.Builder<>("eC", null, Chromosome.class).guiName("End Chromosome").required(false)
            .description("End Chromosome").build();
//    private PluginParameter<Boolean> myIncludeRareAlleles = new PluginParameter.Builder<>("inclRare", false, Boolean.class).guiName("Include Rare Alleles")
//            .description("Include the rare alleles at site (3 or 4th states)").build();
    private PluginParameter<Boolean> myIncludeGaps = new PluginParameter.Builder<>("inclGaps", false, Boolean.class).guiName("Include Gaps")
            .description("Include sites where major or minor allele is a GAP").build();
//    private PluginParameter<Boolean> myCallBiSNPsWGap = new PluginParameter.Builder<>("callBiSNPsWGap", false, Boolean.class).guiName("Call Biallelic SNPs with Gap")
//            .description("Include sites where the third allele is a GAP (mutually exclusive with inclGaps)").build();
    private PluginParameter<Double> myGapAlignmentThreshold = new PluginParameter.Builder<>("gapAlignRatio", 1.0, Double.class).guiName("Gap Alignment Threshold")
            .description("Gap alignment threshold ratio of indel contrasts to non indel contrasts: IC/(IC + NC)."
            		+ " Any loci with a tag alignment value above this threshold will be excluded from the pool.").build();
    private PluginParameter<Integer> maxTagsPerCutSite = new PluginParameter.Builder<Integer>("maxTagsCutSite", 64, Integer.class).guiName("Max Number of Cut Sites").required(false)
            .description("Maximum number of tags per cut site").build();
    private PluginParameter<Boolean> myDeleteOldData = new PluginParameter.Builder<>("deleteOldData", true, Boolean.class).guiName("Delete Previous Discovery Data")
            .description("Delete existing SNP data from tables").build();
    
    private TagDataWriter tagDataWriter = null;
    boolean includeReference = false;
    static GenomeSequence myRefSequence = null; 
    private boolean customSNPLogging = true;  // a custom SNP log that collects useful info for filtering SNPs through machine learning criteria
//    private CustomSNPLog myCustomSNPLog = null;
    private boolean customFiltering = false;
    private int refTagsNotAdded = 0;

    public DiscoverySNPCallerPluginV2() {
        super(null, false);
    }

    public DiscoverySNPCallerPluginV2(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        myLogger.info("Finding SNPs in " + inputDB() + ".");
        myLogger.info(String.format("StartChr:%s EndChr:%s %n", startChromosome(), endChromosome()));
        tagDataWriter =new TagDataSQLite(inputDB());
 
        if (deleteOldData()) {
            myLogger.info("deleteOldData is TRUE: Clearing previous Discovery and SNPQuality data");
            tagDataWriter.clearSNPQualityData();
            tagDataWriter.clearDiscoveryData();
        }
        // Get list of stored chromosomes, we'll process a subset of this list
        List<Chromosome> myChroms = tagDataWriter.getChromosomesFromCutPositions();
        if (myChroms == null || myChroms.size() == 0) {
            myLogger.error("No Chromosomes found in cutPosition tables");
            try{
                ((TagDataSQLite)tagDataWriter).close();
            } catch (Exception ex) {ex.printStackTrace();}
            return null;
        }
        Collections.sort(myChroms); // put in order
       
        Chromosome startChrom = myStartChr.isEmpty() ? myChroms.get(0) : startChromosome();
        Chromosome endChrom = myEndChr.isEmpty() ? myChroms.get(myChroms.size()-1) : endChromosome();
        if (startChrom.compareTo(endChrom) > 0) {
            String message = "The start chromosome " + startChrom.getName() 
                    + " is larger than the end chromosome " + endChrom.getName();
            myLogger.error(message);
            try{
                ((TagDataSQLite)tagDataWriter).close();
            } catch (Exception ex) {ex.printStackTrace();}
            return null;
        } 
        List<Chromosome> chromsToProcess=myChroms.stream()
                .filter(chrom -> { 
                    // Only keep chromosomes within our range.  The chromosome 
                    // constructor stored the string as a name, but also attempted
                    // to create an integer from the name.  
                    //
                    // The "compareTo" method of Chromosome compares the chrom number if
                    // available, otherwise does a string compare on the name. 
                    return (chrom.compareTo(startChrom) >=0 && chrom.compareTo(endChrom) <=0);
                })
                .collect(Collectors.toList());

        chromsToProcess.stream()
            .forEach(chr -> {
            	        myLogger.info("Start processing chromosome " + chr + "\n");
                        Multimap<Tag, Allele> chromosomeAllelemap = HashMultimap.create();
//                        tagDataWriter.getCutPositionTagTaxaMap(chr, -1, -1).entrySet().stream()
//                                .forEach(emp -> {
//                                    Multimap<Tag, Allele> tm = findAlleleByAlignment(emp.getKey(), emp.getValue(), chr);
//                                    if (tm != null) chromosomeAllelemap.putAll(tm);
//                                }); 
                        

                        // Get Tags on forward strands that map to cut positions 
                        myLogger.info("Calling getCutPosForStrand FORWARD strands...");
                        tagDataWriter.getCutPosForStrandTagTaxaMap(chr, -1, -1,true).entrySet().stream()
                        .forEach(emp -> {
                            Multimap<Tag, Allele> tm = findAlleleByAlignment(emp.getKey(), emp.getValue(), chr, true);
                            if (tm != null) chromosomeAllelemap.putAll(tm);
                        });
                        // Get Tags on reverse strands that map to cut positions
                        myLogger.info("\nCalling getCutPosForStrand REVERSE strands...");
                        tagDataWriter.getCutPosForStrandTagTaxaMap(chr, -1, -1,false).entrySet().stream()
                        .forEach(emp -> {
                            Multimap<Tag, Allele> tm = findAlleleByAlignment(emp.getKey(), emp.getValue(), chr, false);
                            if (tm != null) chromosomeAllelemap.putAll(tm);
                        });
                        
                        myLogger.info("Finished processing chromosome " + chr + "\n\n");
                        tagDataWriter.putTagAlleles(chromosomeAllelemap);                       
                    }
            );
        ConcurrencyTools.shutdown();
        //System.out.println("DIscoverySNPCaller - number of referencd tags created but not added: " + refTagsNotAdded);
        try{
            ((TagDataSQLite)tagDataWriter).close();
        } catch (Exception ex) {ex.printStackTrace();}
        return null;
    }

    @Override
    public void postProcessParameters() {

        if (myInputDB.isEmpty() || !Files.exists(Paths.get(inputDB()))) {
            throw new IllegalArgumentException("DiscoverySNPCallerPlugin: postProcessParameters: Input DB not set or found");
        }
        if (!myRefGenome.isEmpty()) {
            includeReference = true;
            myRefSequence = GenomeSequenceBuilder.instance(referenceGenomeFile());
        }
//        if (callBiallelicSNPsWithGap() && includeGaps()) {
//            throw new IllegalArgumentException("The callBiSNPsWGap option is mutually exclusive with the inclGaps option.");
//        }
        if (!myStartChr.isEmpty() && !myEndChr.isEmpty()) {
            if (startChromosome().compareTo(endChromosome()) > 0) {
                throw new IllegalArgumentException("The start chromosome is larger than the end chromosome.");
            } 
        } 
        myLogger.info(String.format("MinMAF:%g %n", minMinorAlleleFreq()));
//        myLogger.info(String.format("includeRare:%s includeGaps:%s %n", includeRareAlleles(), includeGaps()));
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Discovery SNP Caller";
    }

    @Override
    public String getToolTipText() {
        return "Discovery SNP Caller";
    }


    /**
     * Method takes all tags and their taxa distributions at a single cut Position and then identifies the segregating
     * SNPs relative to the reference genome.  Each tag can result in multiple allele being called (Allele contain SNP
     * Position as one of their attributes).  It does not record, which taxa have which SNPs.
     * The steps are:
     * 1. Ensure more than 1 tag (otherwise would be monomorphic).
     * 2. Get the reference sequence
     * 3. Ensure meets minimal coverage
     * 4. Align using BioJava 
     * 5. Filters the aligned tags against a user defined gap alignment ratio threshold (default is include all tags)
     * 6. convert remaining tags to Guava Table 
     * 7. Evaluate whether SNPs meet MAF and coverage levels.
     * 8. Return SNPs in Tag -> Allele map
     * genome.
     * @param cutPosition the cut position that all tags start with
     * @param tagTaxaMap  map of Tag -> Tuple (Boolean if reference, TaxaDistribution)
     * @return multimap of tag -> allele
     */
    Multimap<Tag,Allele> findAlleleByAlignment(Position cutPosition, Map<Tag,TaxaDistribution> tagTaxaMap, Chromosome chromosome, boolean forwardStrand) {
        if(tagTaxaMap.isEmpty()) return null;  //todo why would this be empty?
        final int numberOfTaxa=tagTaxaMap.values().stream().findFirst().get().maxTaxa();
        if((minMinorAlleleFreq()>0) && (tagTaxaMap.size()<2)) {//homozygous
            return null;  //consider reporting homozygous
        }
 
        if(!tagTaxaMap.keySet().stream().anyMatch(Tag::isReference)) {
        	if (myRefSequence != null) {
               	   tagTaxaMap = createReferenceTag(cutPosition, tagTaxaMap, chromosome, numberOfTaxa, forwardStrand);
               	   if (tagTaxaMap == null) return null;
        	} else {
        		tagTaxaMap=setCommonToReference(tagTaxaMap);
        	}
        }
        boolean printDebug=(cutPosition.getPosition()>179_000 && cutPosition.getPosition()<1_500_000);  //todo remove after debugging
        printDebug = false;
        if(printDebug) System.out.println("\nTagLocus: "+cutPosition.toString());  //todo remove after debugging
        final double taxaCoverage=tagTaxaMap.values().stream().mapToInt(t -> t.numberOfTaxaWithTag()).sum()/(double)numberOfTaxa;  //todo this could be changed to taxa with tag
        if(printDebug) System.out.println("taxaCoverage = " + taxaCoverage+" myMinLocusCoverage:"+myMinLocusCoverage.value());
        if(taxaCoverage < myMinLocusCoverage.value()) {
            return null;  //consider reporting low coverage
        }
        
        // This aligns the tags against each other - it doesn't call SNPs
        Map<Tag,String> alignedTagsUnfiltered=alignTags(tagTaxaMap,maxTagsPerCutSite(),cutPosition.getStrand(),printDebug);
        if (alignedTagsUnfiltered == null || alignedTagsUnfiltered.size() == 0) {
        	// Errors related to CompoundNotFound were logged in alignTags. 
        	return null;
        }        
        // Filter the aligned tags.  Throw out all tags from a loci
        // that has any tag with a gap ratio that exceeds the threshold
        Map<Tag,String> alignedTags = filterAlignedTags(alignedTagsUnfiltered, cutPosition, myGapAlignmentThreshold.value());
        if (alignedTags == null || alignedTags.size() == 0) {
            return null;
        }
        // Convert the Position, allele, and tagtaxadist values to a table.  "position" is
        // the key (first item in the row).  
        Table<Position, Byte, List<TagTaxaDistribution>> tAlign=convertAlignmentToTagTable(alignedTags, tagTaxaMap,  cutPosition);
        // Positions returned in table now have a reference allele stored
        List<Position> positionToKeep=tAlign.rowMap().entrySet().stream()
                .filter(entry -> (double) numberTaxaAtSiteIgnoreGaps(entry.getValue()) / (double) numberOfTaxa > minLocusCoverage())
                .filter(entry -> {
                    List<Tuple<Byte,Integer>> aC=alleleTaxaCounts(entry.getValue());
                    if(!includeGaps()) {
                        if(aC.get(0).x==GAP_ALLELE) return false;
                        if(aC.size()>1 && (aC.get(1).x==GAP_ALLELE)) return false;
                    }
                    if(minMinorAlleleFreq()<=0) return true;  //permits export of monomorphic SNPs
                    //if(aC.size()>1) System.out.printf("%s %d %g %g %g%n",entry.getKey().toString(),aC.get(1).y, minMinorAlleleFreq(),taxaCoverage,(minMinorAlleleFreq()*taxaCoverage*(double)numberOfTaxa));
                    return (aC.size()>1 && aC.get(1).y>(minMinorAlleleFreq()*taxaCoverage*(double)numberOfTaxa));
                })
                .map(Map.Entry::getKey) //get Position
//                .peek(p -> {if(printDebug) System.out.println("SNP:"+p.toString());})  //todo remove after debugging
                .collect(Collectors.toList());
        //todo convert to stream
        Multimap<Tag,Allele> tagAllelemap= HashMultimap.create();
        
        for (Position position: positionToKeep) {
            for (Map.Entry<Byte, List<TagTaxaDistribution>> entry : tAlign.row(position).entrySet()) {
                Allele allele=new SimpleAllele(entry.getKey(),position);
                for (TagTaxaDistribution tagTaxaDistribution : entry.getValue()) {
                	Tag currentTag = tagTaxaDistribution.tag();
                	TaxaDistribution td = tagTaxaDistribution.taxaDist();
                	// Address the problem of adding more tags than should be in the db.
                	// don't want to add the reference tag if it wasn't already there.
                	// Don't want to skip it if it was in the DB.
                	// Total depth should distinguish those that were previously present
                	// from those tags that were not.  The problem with adding it is there
                	// is no taxa distribution for it, meaning we created and empty TD, with
                	// total depth of null.  Storing this in the db can cause problems
                	// down the line if a method attempts to grab the td values.
                	// Jira tasks for this were TAS-761 (added), TAS-1006 (removed),
                	// TAS-1138 (changed to totalDepth check)
                	if (td.totalDepth() == 0) {                	    
                	    refTagsNotAdded++;
                	} 
                	else
                	    tagAllelemap.put(currentTag, allele);   
                	//The code below has been in and out. It is replaced by the
                	// code above, which is a better way of determining which
                	// tags were created just for aligning with a reference.

//                	if (currentTag.isReference() && myRefSequence != null) {
//                	    refTagsNotAdded++;
//                	} //don't add ref tag from ref genome to map
//                	else
//                        tagAllelemap.put(currentTag,allele);
                }
            }
        }
        //System.out.printf("%s SNPNum:%d \n",cutPosition.toString(),tagAllelemap.size());
        return tagAllelemap;
    }

    private static Map<Tag,TaxaDistribution> setCommonToReference(Map<Tag,TaxaDistribution> tagTaxaMap) {
        Tag commonTag=tagTaxaMap.entrySet().stream()
                .max(Comparator.comparingInt(e -> e.getValue().numberOfTaxaWithTag()))
                .map(e -> e.getKey())
                .get();
        TaxaDistribution commonTD=tagTaxaMap.get(commonTag);
        Tag refTag=TagBuilder.instance(commonTag.seq2Bit(),commonTag.seqLength()).reference().build();
        ImmutableMap.Builder<Tag,TaxaDistribution> tagTaxaMapBuilder=new ImmutableMap.Builder<>();
        for (Map.Entry<Tag, TaxaDistribution> entry : tagTaxaMap.entrySet()) {
            if(entry.getKey()!=commonTag) {tagTaxaMapBuilder.put(entry);}
            else {tagTaxaMapBuilder.put(refTag,commonTD);}
        }
        return tagTaxaMapBuilder.build();
    }
    private static Map<Tag,TaxaDistribution> createReferenceTag(Position cutPos, 
            Map<Tag,TaxaDistribution> tagTaxaMap, Chromosome myChrom, int numberOfTaxa, boolean forwardStrand) {

        Tag refTag = null;
        int cutPosition = cutPos.getPosition();
         
        // Find the longest sequence length of all tags
//        short longestTag = tagTaxaMap.keySet().stream()
//                .max(Comparator.comparingInt(key -> key.seqLength()))
//                .map(key -> key.seqLength())
//                .get();
        Tag longestTag = tagTaxaMap.keySet().stream()
                .max(Comparator.comparingInt(key -> key.seqLength()))
                .get();
        short longestTagLen = longestTag.seqLength();
        
        // Aligners create the SAM files with the "pos" variable and the "seq" value relative
        // to the forward strand.

        int seqStart =  cutPosition ;
        int seqEnd = cutPosition + longestTagLen-1 ;
        if (!forwardStrand) {
            seqStart =  cutPosition - (longestTagLen-1);
            seqEnd = cutPosition;
        }
 
        byte[] seqInBytes;
        try {
            // chromosomeSequence() now throws an exception on failure.
            // Catch error and record message.
            seqInBytes = myRefSequence.chromosomeSequence(myChrom, seqStart, seqEnd );
        } catch (IllegalArgumentException iae) {
            String msg = "Error creating reference tag at position " + cutPosition + " with length " + longestTagLen
                    + ", seqStart " + seqStart + " seqEnd " + seqEnd + " isForwardStrand " + forwardStrand
                    + ". Position not found in reference file.  " 
                    + ". Please verify the reference file used for the plugin matches the reference file used for the aligner.";
            myLogger.error(msg);
            return null;
        }
 
        String seqInBytesString = NucleotideAlignmentConstants.nucleotideBytetoString(seqInBytes); 
        if (seqInBytesString == null) {
            System.out.println(" createReferenceTag: seqInBytesString is null, seqInBytes value: " + seqInBytes);
        }
        // "null" may be returned from NucleotideBytetoString() for anything that isn't
        // AGCTN.  The Wheat fasta file has R,Y,K,M values, which get translated to 
        // UNDEFINED_ALLELE, which is "6".  These are translated to "null" by nucleotideBytetoString() above
        if (seqInBytesString.contains("N") || seqInBytesString.contains("null")) {
            // This gets lots of hits in wheat - commenting out, but leave for debugging
            //System.out.println("createReferenceTag: reftag contains non-ACGT character, returning Null for cutPosition " + cutPosition + " forwardStrand:" + forwardStrand);
            return null;
        }
        refTag=TagBuilder.instance(seqInBytesString).reference().build();
        if (refTag== null) {
            // We caught bad data above, so this probably means we hit a long string of T's.
            // Two examples are from barley:
            //   CAGCAGGAGAAAGTATGATACTTTGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
            //   AATGTGTATGTCTATGCCAACAAACGGTGCTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
            // These strings have 32 T's, and BaseEncoder.getLongFromSeq() returns -1 for this, as a result of continually
            // shifting and adding 3's as per this code:
            // long v=0; - then in a loop that keeps hitting T's:
            // v = (v << 2) + (byte) 3;
            // after 32 T's, we have -1
            System.out.println("\ncreateReferenceTag: forward strand refTag is null, seqInBytesString: " + seqInBytesString + "\n");
            return null;
        }

        if(!forwardStrand) {           
            refTag=TagBuilder.reverseComplement(refTag).reference().build();
            if (refTag== null) {
                // See problem with 32 T's as above.
                System.out.println("\ncreateReferenceTag: reverse complemented refTag is null for seqInBytesString: " + seqInBytesString + "\n");
                return null;
            }
        }
                
        ImmutableMap.Builder<Tag,TaxaDistribution> tagTaxaMapBuilder=new ImmutableMap.Builder<>();
        //TaxaDistribution refTD=TaxaDistBuilder.create(numberOfTaxa); // is numberOfTaxa an appropriate value to use? 
        TaxaDistribution refTD;
        // If reference tag is not represented, create it and put it on
        // the map.  If it is represented, set this tag to be reference, grab the
        // taxa distribution from the existing tag, and put them on the map
        int mapCount = 0;
        if (!tagTaxaMap.containsKey(refTag)) {
            refTD = TaxaDistBuilder.create(numberOfTaxa);
            tagTaxaMapBuilder.put(refTag,refTD);
            mapCount++;
        }  else {
            // Remove the tag, re-add it with reference set
            refTD=tagTaxaMap.get(refTag);
            tagTaxaMap.remove(refTag);
            tagTaxaMapBuilder.put(refTag,refTD);
            mapCount++;
        }

        for (Map.Entry<Tag, TaxaDistribution> entry : tagTaxaMap.entrySet()) {
            if(!(entry.getKey().equals(refTag))) {
                tagTaxaMapBuilder.put(entry);
                mapCount++;
            }
            // else - do nothing - refTag is already on the map
        }
        return tagTaxaMapBuilder.build();
    }

    private static List<Tuple<Byte,Integer>> alleleTaxaCounts(Map<Byte, List<TagTaxaDistribution>> alleleDistMap) {
        List<Tuple<Byte,Integer>> alleleCnt = new ArrayList<>();
        for (Map.Entry<Byte, List<TagTaxaDistribution>> entry : alleleDistMap.entrySet()) {
            alleleCnt.add(new Tuple<>(entry.getKey(),numberOfTaxa(entry.getValue())));
        }
        Collections.sort(alleleCnt, (Tuple<Byte, Integer> a, Tuple<Byte, Integer> b) -> b.y.compareTo(a.y));  //a,b reversed so common allele  first
        return alleleCnt;
    }

    /**
     * Aligns a set of tags anchored to the same reference position.
     * Tags have been pre-sorted to align
     * @return map with tag(values) mapping to String with alignment
     */
    private static Map<Tag,String> alignTags(Map<Tag, TaxaDistribution> tags, int maxTagsPerCutSite, byte strand, boolean printDebug) {
        List<DNASequence> lst=new ArrayList<>();
        
        for (Map.Entry<Tag,  TaxaDistribution> entry : tags.entrySet())
        {
            Tag tag = entry.getKey();
            if(printDebug) System.out.println(tag.toString());
            String sequence = (strand == 1) ? tag.sequence() : tag.toReverseComplement();
            try {
                DNASequence ds = new DNASequence(sequence);
                ds.setUserCollection(ImmutableList.of(tag));
                lst.add(ds);
            } catch (CompoundNotFoundException ex) {
                // This shouldn't occur and indicates a coding error in previous processing
                myLogger.error("DSNPCaller:alignTags, compoundNotFound exception from DNASequence call for: " + sequence);
                myLogger.debug(ex.getMessage(),ex); 
                return null;
            }
        }        
        ImmutableMap.Builder<Tag,String> result=new ImmutableMap.Builder<>();
        if(lst.size()==1) {
            Tag tag=(Tag)((ImmutableList)lst.get(0).getUserCollection()).get(0);
            result.put(tag,tag.sequence());
            return result.build();
        }
        
        if (lst.size() > maxTagsPerCutSite) {
            // biojava getMultipleSequenceAligment() can handle aligning only so many tags.
            return null;
        }
        // Alignments.getmultipleSequenceAlignment aligns the tags against each other using
        // the ClustalW algorithm
        Profile<DNASequence, NucleotideCompound> profile = Alignments.getMultipleSequenceAlignment(lst);
        if(printDebug) System.out.printf("Clustalw:%n%s%n", profile);
        for (AlignedSequence<DNASequence, NucleotideCompound> compounds : profile) {
            ImmutableList tagList=(ImmutableList)compounds.getOriginalSequence().getUserCollection();
            result.put((Tag)tagList.get(0),compounds.getSequenceAsString());
        }
        return result.build();
    }
    
    /**
     * Takes the aligned tags and compares to the reference genome. If there are too many indels
     * for the tag to be a good match,throw out the alignment. Threshold is user defined via the  
     * myGapAlignmentThreshold plugin parameter.  Default value is 1.0 (keep all tags)
     *   Algorithm:  gap alignment ratio = IC/(IC + NC)
     */
    private static Map<Tag,String> filterAlignedTags(Map<Tag,String> alignedTags, Position refStartPosition, double threshold) {
    	Map<Tag,String> filteredTags = new HashMap<Tag, String>();
        final List<Optional<Position>> referencePositions=referencePositions(refStartPosition,alignedTags);
        // Java 8 "alignedTags.forEach(tag,value) -> {...} " is not used as it does not support "break".
        for (Map.Entry<Tag, String> entry : alignedTags.entrySet())
        {
            Tag tag = entry.getKey();
            String value = entry.getValue();
            int IC = 0;
            int NC = 0;
            for (int index = 0; index < value.length(); index++) {
                byte allele=getNucleotideAlleleByte(value.charAt(index));
                Optional<Position> p=referencePositions.get(index);
                if(!p.isPresent()){
                    // If both are missing, continue.  OTherwise increment indel-contrast
                    if (allele == GAP_ALLELE) continue;
                    IC++;
                } else {
                    if (allele == GAP_ALLELE) IC++;
                    else NC++;
                }                
            }
            
            double alignValue = (double)(IC / (double)(IC + NC));
            if (alignValue <= threshold) {
            	filteredTags.put(tag, value);
            } else {
            	// Toss the entire loci if any of the tags fail the alignment threshold.
            	filteredTags.clear();
            	return null;
            }
        }
        return filteredTags;
    }
    /**
     * Converts at TagAlignment to a Guava Table with
     * @param alignedTags
     * @param tagTaxaDistMap
     * @param refStartPosition
     * @return
     */
    private static Table<Position, Byte, List<TagTaxaDistribution>> convertAlignmentToTagTable(Map<Tag,String> alignedTags,
                        Map<Tag,TaxaDistribution> tagTaxaDistMap, Position cutPosition) {
        Table<Position, Byte, List<TagTaxaDistribution>> alignT= TreeBasedTable.create(); //These could be sorted by depth
        final List<Optional<Position>> referencePositions=referencePositions(cutPosition,alignedTags);
        alignedTags.forEach((t,s) -> {
            TagTaxaDistribution td=new TagTaxaDistribution(t,tagTaxaDistMap.get(t));
            for (int i = 0; i < s.length(); i++) {
                byte allele=getNucleotideAlleleByte(s.charAt(i));
                Optional<Position> p=referencePositions.get(i);
                if(!p.isPresent()) continue;  //skip alignment sites not in reference
                List tdList=alignT.get(p.get(),allele);
                if(tdList!=null) {tdList.add(td);} // Adding more taxaDist for this allele at this position.
                else {
                    List<TagTaxaDistribution> ttdL=new ArrayList<>();
                    ttdL.add(td);
                    alignT.put(p.get(),allele, ttdL);
                }
            }
        });
        return alignT;
    }

    private static String toString(Table<Position, Byte, List<TagTaxaDistribution>> snpTagTable) {
        StringBuilder sb=new StringBuilder();
        sb.append(String.format("Rows %d Columns %d\n", snpTagTable.rowKeySet().size(), snpTagTable.columnKeySet().size()));
        sb.append("Columns bases:");
        snpTagTable.columnKeySet().stream().forEach(b -> sb.append(b+","));
        sb.append("\n");
        snpTagTable.rowMap().forEach((p, ttd) -> {
            sb.append(String.format("Position: %d", p.getPosition()));
            for (Map.Entry<Byte, List<TagTaxaDistribution>> byteListEntry : ttd.entrySet()) {
                sb.append(String.format("\tAllele: %d ->", byteListEntry.getKey()));
                byteListEntry.getValue().forEach(tt -> sb.append(String.join(",",tt.tag().sequence()+"#"+tt.taxaDist().numberOfTaxaWithTag()+" ")));//not doing this correctly
            }
            sb.append("\n");

        });
        return sb.toString();
    }

    private static int numberOfTaxa(List<TagTaxaDistribution> tagTaxaDistributions) {
        return tagTaxaDistributions.stream().mapToInt(t -> t.taxaDist().numberOfTaxaWithTag()).sum();
    }

    private static int numberTaxaAtSite(Map<Byte,List<TagTaxaDistribution>> tagTaxaDistributions) {
        return tagTaxaDistributions.values().stream()
                .flatMap(tList -> tList.stream())
                .mapToInt(t -> t.taxaDist().numberOfTaxaWithTag())
                .sum();
    }

    private static int numberTaxaAtSiteIgnoreGaps(Map<Byte,List<TagTaxaDistribution>> tagTaxaDistributions) {
        return tagTaxaDistributions.entrySet().stream()
                .filter(e -> e.getKey()!=GAP_ALLELE)
                .flatMap(tList -> tList.getValue().stream())
                .mapToInt(tag -> tag.taxaDist().numberOfTaxaWithTag())
                .sum();
    }

    private static List<Optional<Position>> referencePositions(Position cutPosition, Map<Tag,String> alignedTags){
        // Find the tag marked as the reference.
        Tag refTag=alignedTags.keySet().stream()
                .filter(Tag::isReference)
                .findFirst().orElseThrow(() -> new IllegalStateException("Reference not found"));
        AtomicInteger start=new AtomicInteger(cutPosition.getPosition());
        // This takes all the characters in the reference tag.  As long as the character
        // is not equal to a GAP character, it creates a Position object.  The initial
        // int position in the object is the cut position.  As we move through the sequence,
        // the position is incremented by 1 and that is the position stored in the Position object

 
        // Create position list that contains the reference allele
        // The cut position is relative to the forward strand
        String refChars = alignedTags.get(refTag);
        byte strand = cutPosition.getStrand();
        List <Optional<Position>> positionList = new ArrayList<Optional<Position>>();
        if (strand == 1) {
            for (int cidx = 0; cidx < refChars.length(); cidx++) {
                char currChar = refChars.charAt(cidx);
                if (currChar == GAP_ALLELE_CHAR) {
                    positionList.add(Optional.<Position>empty());
                } else {
                    byte refAllele = NucleotideAlignmentConstants.getNucleotideAlleleByte(currChar);
                    // These positions are always relative to the forward strand
                    // if (strand == -1) refAllele = NucleotideAlignmentConstants.getNucleotideComplement(refAllele);
                    positionList.add(Optional.<Position>of(new GeneralPosition.Builder(cutPosition.getChromosome(), 
                            start.getAndIncrement()).allele(WHICH_ALLELE.Reference, refAllele).build()));
                }
            } 
        } else {
            List <Optional<Position>> revPositionList = new ArrayList<Optional<Position>>();
            for (int revCidx = refChars.length()-1; revCidx >= 0; revCidx--) {
                char currChar = refChars.charAt(revCidx);
                if (currChar == GAP_ALLELE_CHAR) {
                    revPositionList.add(Optional.<Position>empty());
                } else {
                    byte refAllele = NucleotideAlignmentConstants.getNucleotideAlleleByte(currChar);
                    revPositionList.add(Optional.<Position>of(new GeneralPosition.Builder(cutPosition.getChromosome(), 
                            start.getAndDecrement()).allele(WHICH_ALLELE.Reference, refAllele).build()));
                }
            }
            for (int revCidx = revPositionList.size()-1; revCidx >= 0; revCidx--) {
                positionList.add(revPositionList.get(revCidx));
            }
        }
        return positionList;
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(DiscoverySNPCallerPluginV2.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
//    // TODO: Replace <Type> with specific type.
//    public <TagDataWriter> runPlugin(DataSet input) {
//        return (TagDataWriter) performFunction(input).getData(0).getData();
//    }

    /**
     * Input TagsByTaxa file (if hdf5 format, use .hdf or
     * .h5 extension)
     *
     * @return Input Tags by Taxa File
     */
    public String inputDB() {
        return myInputDB.value();
    }

    /**
     * Set Input Tags by Taxa File. Input TagsByTaxa file
     * (if hdf5 format, use .hdf or .h5 extension)
     *
     * @param value Input Tags by Taxa File
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 inputDB(String value) {
        myInputDB = new PluginParameter<>(myInputDB, value);
        return this;
    }

    /**
     * Minimum minor allele frequency
     *
     * @return Min Minor Allele Freq
     */
    public Double minMinorAlleleFreq() {
        return myMinMinorAlleleFreq.value();
    }

    /**
     * Set Min Minor Allele Freq. Minimum minor allele frequency
     *
     * @param value Min Minor Allele Freq
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 minMinorAlleleFreq(Double value) {
        myMinMinorAlleleFreq = new PluginParameter<>(myMinMinorAlleleFreq, value);
        return this;
    }

    /**
     * Minimum locus coverage (proportion of Taxa with a genotype)
     *
     * @return Min Locus Coverage
     */
    public Double minLocusCoverage() {
        return myMinLocusCoverage.value();
    }

    /**
     * Set Min Locus Coverage. Minimum locus coverage (proportion
     * of Taxa with a genotype)
     *
     * @param value Min Locus Coverage
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 minLocusCoverage(Double value) {
        myMinLocusCoverage = new PluginParameter<>(myMinLocusCoverage, value);
        return this;
    }

    /**
     * Path to reference genome in fasta format. Ensures that
     * a tag from the reference genome is always included
     * when the tags at a locus are aligned against each other
     * to call SNPs. The reference allele for each site is
     * then provided in the output HapMap files, under the
     * taxon name "REFERENCE_GENOME" (first taxon). DEFAULT:
     * Don't use reference genome.
     *
     * @return Reference Genome File
     */
    public String referenceGenomeFile() {
        return myRefGenome.value();
    }

    /**
     * Set Reference Genome File. Path to reference genome
     * in fasta format. Ensures that a tag from the reference
     * genome is always included when the tags at a locus
     * are aligned against each other to call SNPs. The reference
     * allele for each site is then provided in the output
     * HapMap files, under the taxon name "REFERENCE_GENOME"
     * (first taxon). DEFAULT: Don't use reference genome.
     *
     * @param value Reference Genome File
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 referenceGenomeFile(String value) {
        myRefGenome = new PluginParameter<>(myRefGenome, value);
        return this;
    }

    /**
     * Start Chromosome
     *
     * @return Start Chromosome
     */
    public Chromosome startChromosome() {
        return myStartChr.value();
    }

    /**
     * Set Start Chromosome. Start Chromosome
     *
     * @param value Start Chromosome
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 startChromosome(Chromosome value) {
        myStartChr = new PluginParameter<>(myStartChr, value);
        return this;
    }

    /**
     * End Chromosome
     *
     * @return End Chromosome
     */
    public Chromosome endChromosome() {
        return myEndChr.value();
    }

    /**
     * Set End Chromosome. End Chromosome
     *
     * @param value End Chromosome
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 endChromosome(Chromosome value) {
        myEndChr = new PluginParameter<>(myEndChr, value);
        return this;
    }

//    /**
//     * Include the rare alleles at site (3 or 4th states)
//     *
//     * @return Include Rare Alleles
//     */
//    public Boolean includeRareAlleles() {
//        return myIncludeRareAlleles.value();
//    }
//
//    /**
//     * Set Include Rare Alleles. Include the rare alleles
//     * at site (3 or 4th states)
//     *
//     * @param value Include Rare Alleles
//     *
//     * @return this plugin
//     */
//    public DiscoverySNPCallerPluginV2 includeRareAlleles(Boolean value) {
//        myIncludeRareAlleles = new PluginParameter<>(myIncludeRareAlleles, value);
//        return this;
//    }

    /**
     * Include sites where major or minor allele is a GAP
     *
     * @return Include Gaps
     */
    public Boolean includeGaps() {
        return myIncludeGaps.value();
    }

    /**
     * Set Include Gaps. Include sites where major or minor
     * allele is a GAP
     *
     * @param value Include Gaps
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 includeGaps(Boolean value) {
        myIncludeGaps = new PluginParameter<>(myIncludeGaps, value);
        return this;
    }

//    /**
//     * Include sites where the third allele is a GAP (mutually
//     * exclusive with inclGaps)
//     *
//     * @return Call Biallelic SNPs with Gap
//     */
//    public Boolean callBiallelicSNPsWithGap() {
//        return myCallBiSNPsWGap.value();
//    }
//
//    /**
//     * Set Call Biallelic SNPs with Gap. Include sites where
//     * the third allele is a GAP (mutually exclusive with
//     * inclGaps)
//     *
//     * @param value Call Biallelic SNPs with Gap
//     *
//     * @return this plugin
//     */
//    public DiscoverySNPCallerPluginV2 callBiallelicSNPsWithGap(Boolean value) {
//        myCallBiSNPsWGap = new PluginParameter<>(myCallBiSNPsWGap, value);
//        return this;
//    }

    /**
     * Maximum gap alignment allowed from the equation:
     *   IndelContrast / (IndelContrast/Non-IndelContrast)
     *   
     *   IC=Indel contrasts=Sum the number ACGT vs - 
     *   NC=non-indel constrasts = Sum the number of ACGT vs ACGT
     *   ignore = - vs -
     * Gapped Alignment ratio = IC/(IC+NC)
     *
     * @return Maxmimum Gap alignment ratio
     */
    public Double gapAlignmentThreshold() {
        return myGapAlignmentThreshold.value();
    }

    /**
     * Maximum gap alignment allowed from the equation:
     *   IndelContrast / (IndelContrast/Non-IndelContrast)
     *   
     *   IC=Indel contrasts=Sum the number ACGT vs - 
     *   NC=non-indel constrasts = Sum the number of ACGT vs ACGT
     *   ignore = - vs -
     * Gapped Alignment ratio = IC/(IC+NC)
     *
     * @param value Max gap alignment ratio
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 gapAlignmentThreshold(Double value) {
        myGapAlignmentThreshold = new PluginParameter<>(myGapAlignmentThreshold, value);
        return this;
    }
    
    /**
     * Maximum number of tags per cut site (for alignment)
     *
     * @return MaxTagsPerCutSite
     */
    public Integer maxTagsPerCutSite() {
        return maxTagsPerCutSite.value();
    }

    /**
     * Set maxTagsPerCutSite. This is the maximum number of tags
     * allowed per cute site when performaing an alignment.  Too
     * many tags and biojava 4 getMultipleSequenceAlignment grinds to a halt.
     *
     * @param value Max number of tags per cut site
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 maxTagsPerCutSite(Integer value) {
    	maxTagsPerCutSite = new PluginParameter<>(maxTagsPerCutSite, value);
        return this;
    }

    private static Map<String,Integer> keyFileStringToInt = null;

    // For junit testing.  Used for ReferenceGenomeSequence:readReferenceGenomeChr() tests
    // that want to read from a key file.
    public static void setKeyFileStringToInt (Map<String,Integer> inputKeyFile) {
        keyFileStringToInt = inputKeyFile;
    }

    /**
     * Delete exisiting Discovery data from DB
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
    public DiscoverySNPCallerPluginV2 deleteOldData(Boolean value) {
        myDeleteOldData = new PluginParameter<>(myDeleteOldData, value);
        return this;
    }
    /**
     * GIven a chromosome string value, search for it's corresponding
     * number from the keyFileStringToInt map
     *
     * @return integer version of Chromosome
     */
    public static int keyFileReturnChromInt(String chromosome) {
        if (keyFileStringToInt != null) {
                if (keyFileStringToInt.containsKey(chromosome))
                        return(keyFileStringToInt.get(chromosome));
        } else {
                System.out.println("LCJ - DiscoverySNP .. keyFileStringToInt is NULL");
        }
        return -1; // failure case
    }
}

class TagTaxaDistribution {
    private final Tag myTag;
    private final TaxaDistribution myTaxaDist;

    public TagTaxaDistribution(Tag myTag, TaxaDistribution myTaxaDist) {
        this.myTag = myTag;
        this.myTaxaDist = myTaxaDist;
    }

    public Tag tag() {
        return myTag;
    }

    public TaxaDistribution taxaDist() {
        return myTaxaDist;
    }
}

