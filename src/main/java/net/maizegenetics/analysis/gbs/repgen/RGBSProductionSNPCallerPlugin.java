/*
 * ProductionSNPCallerPlugin
 */
package net.maizegenetics.analysis.gbs.repgen;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import net.maizegenetics.analysis.gbs.Barcode;
import net.maizegenetics.analysis.gbs.v2.GBSUtils;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.dna.snp.score.AlleleDepthUtil;
import net.maizegenetics.dna.snp.genotypecall.BasicGenotypeMergeRule;
import net.maizegenetics.dna.snp.genotypecall.GenotypeMergeRule;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TagData;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;
import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.Collectors;

/**
 * This plugin converts all of the fastq (and/or qseq) files in the input folder
 * and keyfile to genotypes and adds these to a genotype file in HDF5 format.
 *
 * We refer to this step as the "Production Pipeline".
 *
 * The output format is either HDF5 or VCF genotypes with allelic depth stored. 
 * Output file type is determined by presence of the ".h5" suffix.  SNP calling is
 * quantitative with the option of using either the Glaubitz/Buckler binomial
 * method (pHet/pErr > 1 = het) (=default), or the Stacks method.
 *
 * Merging of samples with the same LibraryPrepID is handled by
 * GenotypeTableBuilder.addTaxon(), with the genotypes re-called based upon the
 * new depths. Therefore, if you want to keep adding genotypes to the same
 * target HDF5 file in subsequent runs, use the -ko (keep open) option so that
 * the output GenotypeTableBuilder will be mutable, using closeUnfinished()
 * rather than build().
 *
 * If the target output is HDF5, and that GenotypeTable file doesn't exist, it will be
 * created.  
 *
 * Each taxon in the output file is named "ShortName:LibraryPrepID" and is
 * annotated with "Flowcell_Lanes" (=source seq data for current genotype).
 *
 * Requires a database with variants added from a previous "Discovery Pipeline" run.
 * 
 * References to "tag" are being replaced by references to "kmer" as the pipeline
 * is really a kmer alignment process.
 *
 * TODO add the Stacks likelihood method to BasicGenotypeMergeRule
 *
 * @author Ed Buckler
 * @author Jeff Glaubitz
 */
public class RGBSProductionSNPCallerPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(RGBSProductionSNPCallerPlugin.class);

    private PluginParameter<String> myInputDir = new PluginParameter.Builder<>("i", null, String.class).guiName("Input Directory").required(true).inDir()
            .description("Input directory containing fastq AND/OR qseq files.").build();
    private PluginParameter<String> myKeyFile = new PluginParameter.Builder<>("k", null, String.class).guiName("Key File").required(true).inFile()
            .description("Key file listing barcodes distinguishing the samples").build();
    private PluginParameter<String> myInputDB = new PluginParameter.Builder<>("db", null, String.class).guiName("Input GBS Database").required(true).inFile()
            .description("Input Database file if using SQLite").build();
    private PluginParameter<String> myOutputGenotypes = new PluginParameter.Builder<>("o", null, String.class).guiName("Output Genotypes File").required(true).outFile()
            .description("Output (target) genotypes file to produce.  Default output file type is VCF.  If file suffix is .h5, an hdf5 file will be created instead.")
            .build();
    private PluginParameter<Double> myAveSeqErrorRate = new PluginParameter.Builder<>("eR", 0.01, Double.class).guiName("Ave Seq Error Rate")
            .description("Average sequencing error rate per base (used to decide between heterozygous and homozygous calls)").build();
    private PluginParameter<Integer> myMaxDivergence = new PluginParameter.Builder<>("d", 0, Integer.class).guiName("Max Divergence")
            .description("Maximum divergence (edit distance) between new read and previously mapped read (Default: 0 = perfect matches only)").build();
    private PluginParameter<Boolean> myDepthOutput = new PluginParameter.Builder<>("do", true, Boolean.class).guiName("Write Depths to Output")
            .description("Depth output: True means write depths to the output hdf5 genotypes file, false means do NOT write depths to the hdf5 file").build();
    private PluginParameter<Integer> myKmerLength = new PluginParameter.Builder<>("kmerLength", 64, Integer.class).guiName("Maximum Kmer Length")
            .description("Length of kmers to process").build();
    private PluginParameter<Double> posQualityScore = new PluginParameter.Builder<>("minPosQS", 0.0, Double.class).guiName("Minimun snp quality score")
            .description("Minimum quality score for snp position to be included").build();
    private PluginParameter<Integer> myBatchSize = new PluginParameter.Builder<>("batchSize", 8, Integer.class).guiName("Batch size of fastq files").required(false)
            .description("Number of flow cells being processed simultaneously").build();
    private PluginParameter<Integer> myMinQualScore = new PluginParameter.Builder<>("mnQS", 0, Integer.class).guiName("Minimum quality score").required(false)
            .description("Minimum quality score within the barcode and read length to be accepted").build();
    //private PluginParameter<Boolean> myStacksLikelihood = new PluginParameter.Builder<>("sL", false, Boolean.class).guiName("Use Stacks Likelihood")
    //        .description("Use STACKS likelihood method to call heterozygotes (default: use tasselGBS likelihood ratio method)").build();

    private String myOutputDir = null;
    private static boolean isHDF5 = false; // default is VCF
    private TagData tagDataReader = null;
    Multimap<Taxon,Tag> tagCntMap=Multimaps.synchronizedMultimap(ArrayListMultimap.create(384, 500_000));
    private Set<String> seqFilesInKeyAndDir = new TreeSet<>(); // fastq (or qseq) file names present in input directory that have a "Flowcell_Lane" in the key file

    protected static int readEndCutSiteRemnantLength;
    private Trie ahoCorasickTrie; // import from ahocorasick-0.2.1.jar
    private String[] likelyReadEndStrings;

    //Documentation of read depth per sample (one recorded per replicate)
    // Treemap is synchronized as multiple threads may increment values.
    private Map<String, Integer> rawReadCountsMap = new TreeMap<>();
    private Map<String, Integer> rawReadCountsForFullSampleName = Collections.synchronizedMap(rawReadCountsMap);
    private Map<String, Integer> matchedReadCountsMap = new TreeMap<>();
    private Map<String, Integer> matchedReadCountsForFullSampleName = Collections.synchronizedMap(matchedReadCountsMap);

    private GenotypeMergeRule genoMergeRule = null;
    private boolean taglenException;

    public RGBSProductionSNPCallerPlugin() {
        super(null, false);
    }

    public RGBSProductionSNPCallerPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public void postProcessParameters() {
        try {
            myOutputDir = (new File(outputGenotypesFile())).getCanonicalFile().getParent();
        } catch (IOException e) {
            throw new IllegalStateException("Problem resolving output directory:" + e);
        }
        genoMergeRule = new BasicGenotypeMergeRule(aveSeqErrorRate());
        
        if (!myOutputGenotypes.isEmpty()) {
            if (outputGenotypesFile().endsWith(".h5")) {
                isHDF5 = true;
            }
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        int batchSize = batchSize();
        List<Path> directoryFiles= DirectoryCrawler.listPaths(GBSUtils.inputFileGlob, Paths.get(myInputDir.value()).toAbsolutePath());
        if(directoryFiles.isEmpty()) {
            myLogger.warn("No files matching:"+GBSUtils.inputFileGlob);
            return null;
        }


        tagDataReader =new TagDataSQLite(myInputDB.value());
        TaxaList masterTaxaList= TaxaListIOUtils.readTaxaAnnotationFile(keyFile(), GBSUtils.sampleNameField, new HashMap<>(), true);
        Map<String,Taxon> fileTaxaMap=TaxaListIOUtils.getUniqueMapOfTaxonByAnnotation(masterTaxaList,GBSUtils.fileNameField)
                .orElseThrow(() -> new IllegalArgumentException("Error: Same file points more than one taxon in the KeyFile"));
        List<Path> inputSeqFiles =  directoryFiles.stream()
                .peek(path -> System.out.println(path.getFileName().toString()))
                .filter(path -> fileTaxaMap.containsKey(path.getFileName().toString()))
                .collect(Collectors.toList());
        if (inputSeqFiles.size() == 0) return null; // no files to process
        writeInitialTaxaReadCounts(masterTaxaList); // initialize synchronized maps
        //todo perhaps subset the masterTaxaList based on the files in there, but it seems like it will all be figure out.
        Map<Tag,Tag> canonicalTag=new HashMap<>();  //canonicalize them OR eventually we will use a Trie
        tagDataReader.getTags().stream().forEach(t -> canonicalTag.put(t,t));
        int batchNum = inputSeqFiles.size()/batchSize;
       
        if (inputSeqFiles.size() % batchSize !=0) batchNum++;
        System.out.println("ProductionSNPCallerPluginV2: Total batches to process: " + batchNum);

        final PositionList positionList=tagDataReader.getSNPPositions(positionQualityScore());
        if (positionList == null || positionList.size() == 0) {
        	String errMsg = "\nNo snp positons found with quality score of " + positionQualityScore() + ".\n"
        			+ "Please run UpdateSNPPositionQualityPlugin to add quality scores for your positions,\n"
        			+ " then select snp positions within a quality range you have specified.\n";
        	myLogger.error(errMsg);
        	return null;
        }
               
        GenotypeTableBuilder gtb=setUpGenotypeTableBuilder(outputGenotypesFile(), positionList, genoMergeRule);
        final Multimap<Tag,AlleleWithPosIndex> tagsToIndex=ArrayListMultimap.create();
        tagDataReader.getAlleleMap().entries().stream()
                .forEach(e -> {
                	// indexOf returns -1 if the list doesn't contain the element, which it won't
                	// if there are snpposition entries with a quality score less than minimumQualityScore 
                    int posIndex=positionList.indexOf(e.getValue().position());
                    if (posIndex >= 0) {
                    	tagsToIndex.put(e.getKey(),new AlleleWithPosIndex(e.getValue(),posIndex));
                    }                   
                });
        
        taglenException = false;
        for (int idx = 0; idx < inputSeqFiles.size(); idx+=batchSize) {
        	tagCntMap.clear(); // start fresh with each new batch
            int end = idx+batchSize;
            if (end > inputSeqFiles.size()) end = inputSeqFiles.size();
            ArrayList<Path> sub = new ArrayList<Path>();
            for (int jdx = idx; jdx < end; jdx++) sub.add(inputSeqFiles.get(jdx));
            System.out.println("\nStart processing batch " + String.valueOf(idx/batchSize+1));
            sub.parallelStream()
            .forEach(inputSeqFile -> {
                try {
                    int taxaIndex=masterTaxaList.indexOf(fileTaxaMap.get(inputSeqFile.getFileName().toString()));
                    //processFastQFile(masterTaxaList,keyPath, inputSeqFile, enzyme(),canonicalTag,kmerLength(), minimumQualityScore());
                    processFastQ(inputSeqFile,taxaIndex,masterTaxaList,canonicalTag,kmerLength(),minimumQualityScore());
                } catch (StringIndexOutOfBoundsException oobe) {
                    oobe.printStackTrace();
                    myLogger.error(oobe.getMessage());
                    setTagLenException();
                    return;
                }              
            });
            if (taglenException == true) return null; // Tag length failure from processFastQ - halt processing
         
            tagCntMap.asMap().entrySet().stream()
            .forEach(e -> {
                callGenotypes(e.getKey(), e.getValue(), tagsToIndex, positionList, genoMergeRule,gtb,depthToOutput());
                //System.out.println(e.x.getName()+ Arrays.toString(Arrays.copyOfRange(e.y,0,10)))); 
            });
            System.out.println("\nFinished processing batch " + String.valueOf(idx/batchSize+1));
        }
        GenotypeTable myGt = gtb.build();
        ExportUtils.writeToVCF(myGt, outputGenotypesFile(), depthToOutput());
        return null;
    }

    private static void callGenotypes(Taxon taxon, Collection<Tag> tags, Multimap<Tag,AlleleWithPosIndex> tagsToIndex,
                   PositionList positionList, GenotypeMergeRule genoMergeRule, GenotypeTableBuilder gtb, boolean outputDepths) {
        int[][] alleleDepths = new int[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][positionList.numberOfSites()];
        tags.stream().map(t -> tagsToIndex.get(t)).flatMap(c -> c.stream())
                .forEach(a -> alleleDepths[a.allele()][a.positionIndex()]++);
        if (outputDepths) {
            byte[][] byteDepths = AlleleDepthUtil.depthIntToByte(alleleDepths);
            gtb.addTaxon(taxon, resolveGenosForTaxon(alleleDepths, genoMergeRule),byteDepths);
        } else {
        	gtb.addTaxon(taxon, resolveGenosForTaxon(alleleDepths, genoMergeRule));
        }
    }

    private class AlleleWithPosIndex extends SimpleAllele {
        private int positionIndex;

        private AlleleWithPosIndex(Allele myAllele, int positionIndex) {
            super(myAllele.allele(), myAllele.position());
            this.positionIndex=positionIndex;
        }

        public int positionIndex() {
            return positionIndex;
        }
    }

    private class CountOfReadQuality {
        LongAdder allReads=new LongAdder();
        LongAdder goodBarcodedReads=new LongAdder();
        LongAdder goodMatched=new LongAdder();
        LongAdder perfectMatches=new LongAdder();
        LongAdder imperfectMatches=new LongAdder();
        LongAdder singleImperfectMatches=new LongAdder();
    }

//    private void processFastQ(Path fastqFile, int taxaIndex, TaxaList masterTaxaList,
//                              TagDistributionMap masterTagTaxaMap, int preferredTagLength, int minQual) throws StringIndexOutOfBoundsException{

    private void processFastQ(Path fastqFile, int taxaIndex, TaxaList masterTaxaList, Map<Tag,Tag> canonicalTags,
            int preferredTagLength, int minQual) throws StringIndexOutOfBoundsException {
        int allReads=0, goodBarcodedReads = 0, lowQualityReads = 0;
        try {
            int qualityScoreBase=GBSUtils.determineQualityScoreBase(fastqFile);
            BufferedReader br = Utils.getBufferedReader(fastqFile.toString(), 1 << 22);
            long time=System.nanoTime();
            String[] seqAndQual;
            Taxon taxon=masterTaxaList.get(taxaIndex);
            while ((seqAndQual=GBSUtils.readFastQBlock(br, allReads)) != null) {
                allReads++;
                // Decode barcode using the current sequence & quality  score

                if(minQual>0) {
                    //todo move getFirstLowQualityPos into this class?
                    if(BaseEncoder.getFirstLowQualityPos(seqAndQual[1],minQual, qualityScoreBase)<(preferredTagLength)){
                        lowQualityReads++;
                        continue;
                    }
                }
//                rawReadCountsForFullSampleName.put(barcode.getTaxaName(), rawReadCountsForFullSampleName.get(barcode.getTaxaName()) + 1);
                if (seqAndQual[0].length() < preferredTagLength) {
                    String errMsg = "\n\nERROR processing " + fastqFile.toString() + "\n" +
                            "Reading entry number " + allReads + " fails the length test.\n" +
                            "Sequence length " + seqAndQual[0].length() +
                            " is less then maxKmerLength " + preferredTagLength + ".\n" +
                            "Re-run your files with either a shorter mxKmerL value or a higher minimum quality score.\n";
                    throw new StringIndexOutOfBoundsException(errMsg);
                }

                Tag tag = TagBuilder.instance(seqAndQual[0].substring(0,preferredTagLength)).build();
                //Tag tag = removeSecondCutSiteAhoC(seqAndQual[0].substring(barcodeLen),preferredTagLength);
                //Tag tag= TagBuilder.instance(seqAndQual[0].substring(barcode.getBarLength(), barcode.getBarLength() + preferredTagLength)).build();
                if(tag==null) continue;   //null occurs when any base was not A, C, G, T
                goodBarcodedReads++;
                Tag canonicalTag=canonicalTags.get(tag);
                if(canonicalTag!=null) {
                    tagCntMap.put(taxon,canonicalTag);
                    //matchedReadCountsForFullSampleName.put(barcode.getTaxaName(), matchedReadCountsForFullSampleName.get(barcode.getTaxaName()) + 1);
                }
                if (allReads % 1000000 == 0) {
                    myLogger.info("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads
                            + " rate:" + (System.nanoTime()-time)/allReads +" ns/read");
                }
            }
            myLogger.info("Total number of reads in lane=" + allReads);
            myLogger.info("Total number of good barcoded reads=" + goodBarcodedReads);
            myLogger.info("Total number of low quality reads=" + lowQualityReads);
            myLogger.info("Timing process (sorting, collapsing, and writing TagCount to file).");
            myLogger.info("Process took " + (System.nanoTime() - time)/1e6 + " milliseconds for file " + fastqFile.toString());
            br.close();
        } catch (Exception e) {
            myLogger.error("Good Barcodes Read: " + goodBarcodedReads);
            e.printStackTrace();
        }
    }



    private void reportProgress(int[] counters, long readSeqReadTime, long ifRRNotNullTime) {
        myLogger.info(
                "totalReads:" + counters[0]
                + "  goodBarcodedReads:" + counters[1]
                + "  goodMatchedToTOPM:" + counters[2]
                //            + "  perfectMatches:" + counters[3]
                //            + "  nearMatches:" + counters[4]
                //            + "  uniqueNearMatches:" + counters[5]
                + "  cumulReadSequenceTime: " + ((double) (readSeqReadTime) / 1_000_000_000.0) + " sec"
                + "  cumulProcessSequenceTime: " + ((double) (ifRRNotNullTime) / 1_000_000_000.0) + " sec"
        );
    }

    private void reportTotals(Path fileName, int[] counters, int nFilesProcessed) {
        myLogger.info("Total number of reads in lane=" + counters[0]);
        myLogger.info("Total number of good, barcoded reads=" + counters[1]);
        myLogger.info("Total number of good, barcoded reads matched to the TOPM=" + counters[2]);
        myLogger.info("Finished reading " + nFilesProcessed + " of " + seqFilesInKeyAndDir.size() + " sequence files: " + fileName + "\n");
    }


    private static GenotypeTableBuilder setUpGenotypeTableBuilder(String anOutputFile, PositionList positionList, GenotypeMergeRule mergeRule) {
        if (isHDF5) {
            File hdf5File = new File(anOutputFile);
            if (hdf5File.exists()) {
                myLogger.info("\nGenotypes will be added to existing HDF5 file:\n  " + anOutputFile + "\n");
                return GenotypeTableBuilder.mergeTaxaIncremental(anOutputFile, mergeRule);
            } else {
                myLogger.info("\nThe target HDF5 file:\n  " + anOutputFile
                        + "\ndoes not exist. A new HDF5 file of that name will be created \nto hold the genotypes from this run.");
                return GenotypeTableBuilder.getTaxaIncrementalWithMerging(anOutputFile, positionList, mergeRule);
            }
        } else { // create genotype table for VCF
            GenotypeTableBuilder gtb = GenotypeTableBuilder.getTaxaIncremental(positionList, mergeRule);
            myLogger.info("\nOutput VCF file: \n" + anOutputFile +
                    " \ncreated for genotypes from this run.");
            return gtb;
        }
    }

    private static byte[] resolveGenosForTaxon(int[][] depthsForTaxon, GenotypeMergeRule genoMergeRule) {
        int nAlleles = depthsForTaxon.length;
        int[] depthsAtSite = new int[nAlleles];
        int nSites = depthsForTaxon[0].length;
        byte[] genos = new byte[nSites];
        for (int site = 0; site < nSites; site++) {
            for (int allele = 0; allele < nAlleles; allele++) {
                depthsAtSite[allele] = depthsForTaxon[allele][site];
            }
            genos[site] = genoMergeRule.callBasedOnDepth(depthsAtSite);
        }
        return genos;
    }

    
    private void writeInitialTaxaReadCounts(TaxaList tl) {
    	tl.stream() // Add initial taxa names with count of 0 to synchronized maps
    	.forEach(taxon -> {
    		 rawReadCountsForFullSampleName.put(taxon.getName(), 0); 
    	     matchedReadCountsForFullSampleName.put(taxon.getName(), 0);
    	});
    }

    public void setTagLenException() {
        taglenException = true;
    }
    
    private void printFileNameConventions(String actualFileName) {
        String message
                = "\n\n"
                + "Error in parsing file name:"
                + "\n   The raw sequence filename does not contain either 3, 4, or 5 underscore-delimited values."
                + "\n   Acceptable file naming conventions include the following (where FLOWCELL indicates the flowcell name and LANE is an integer):"
                + "\n       FLOWCELL_LANE_fastq.gz"
                + "\n       FLOWCELL_s_LANE_fastq.gz"
                + "\n       code_FLOWCELL_s_LANE_fastq.gz"
                + "\n       FLOWCELL_LANE_fastq.txt.gz"
                + "\n       FLOWCELL_s_LANE_fastq.txt.gz"
                + "\n       code_FLOWCELL_s_LANE_fastq.txt.gz"
                + "\n       FLOWCELL_LANE_qseq.txt.gz"
                + "\n       FLOWCELL_s_LANE_qseq.txt.gz"
                + "\n       code_FLOWCELL_s_LANE_qseq.txt.gz"
                + "\n"
                + "\n   Actual Filename: " + actualFileName
                + "\n\n";

        myLogger.error(message);
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Production SNP Caller";
    }

    @Override
    public String getToolTipText() {
        return "Production SNP Caller";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(ProductionSNPCallerPluginV2.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public TagData runPlugin(DataSet input) {
        return (TagData) performFunction(input).getData(0).getData();
    }

    /**
     * Input directory containing fastq AND/OR qseq files.
     *
     * @return Input Directory
     */
    public String inputDirectory() {
        return myInputDir.value();
    }

    /**
     * Set Input Directory. Input directory containing fastq
     * AND/OR qseq files.
     *
     * @param value Input Directory
     *
     * @return this plugin
     */
    public RGBSProductionSNPCallerPlugin inputDirectory(String value) {
        myInputDir = new PluginParameter<>(myInputDir, value);
        return this;
    }

    /**
     * Key file listing barcodes distinguishing the samples
     *
     * @return Key File
     */
    public String keyFile() {
        return myKeyFile.value();
    }

    /**
     * Set Key File. Key file listing barcodes distinguishing
     * the samples
     *
     * @param value Key File
     *
     * @return this plugin
     */
    public RGBSProductionSNPCallerPlugin keyFile(String value) {
        myKeyFile = new PluginParameter<>(myKeyFile, value);
        return this;
    }


    /**
     * Input Database file if using SQLite
     *
     * @return Input GBS Database
     */
    public String inputGBSDatabase() {
        return myInputDB.value();
    }

    /**
     * Set Input GBS Database. Input Database file if using
     * SQLite
     *
     * @param value Input GBS Database
     *
     * @return this plugin
     */
    public RGBSProductionSNPCallerPlugin inputGBSDatabase(String value) {
        myInputDB = new PluginParameter<>(myInputDB, value);
        return this;
    }

    /**
     * Output (target) genotypes file to add new genotypes
     * to (new file created if it doesn't exist)
     *
     * @return Output file Genotypes File
     */
    public String outputGenotypesFile() {
        return myOutputGenotypes.value();
    }

    /**
     * Set Output  Genotypes File. Output (target)
     * genotypes file to add new genotypes to (new file created
     * if it doesn't exist)
     *
     * @param value Output Genotypes File
     *
     * @return this plugin
     */
    public RGBSProductionSNPCallerPlugin outputGenotypesFile(String value) {
        myOutputGenotypes = new PluginParameter<>(myOutputGenotypes, value);
        return this;
    }

    /**
     * Average sequencing error rate per base (used to decide
     * between heterozygous and homozygous calls)
     *
     * @return Ave Seq Error Rate
     */
    public Double aveSeqErrorRate() {
        return myAveSeqErrorRate.value();
    }

    /**
     * Set Ave Seq Error Rate. Average sequencing error rate
     * per base (used to decide between heterozygous and homozygous
     * calls)
     *
     * @param value Ave Seq Error Rate
     *
     * @return this plugin
     */
    public RGBSProductionSNPCallerPlugin aveSeqErrorRate(Double value) {
        myAveSeqErrorRate = new PluginParameter<>(myAveSeqErrorRate, value);
        return this;
    }

    /**
     * Maximum divergence (edit distance) between new read
     * and previously mapped read (Default: 0 = perfect matches
     * only)
     *
     * @return Max Divergence
     */
    public Integer maxDivergence() {
        return myMaxDivergence.value();
    }

    /**
     * Set Max Divergence. Maximum divergence (edit distance)
     * between new read and previously mapped read (Default:
     * 0 = perfect matches only)
     *
     * @param value Max Divergence
     *
     * @return this plugin
     */
    public RGBSProductionSNPCallerPlugin maxDivergence(Integer value) {
        myMaxDivergence = new PluginParameter<>(myMaxDivergence, value);
        return this;
    }

    /**
     * Output depth: write depths to the output
     * hdf5 genotypes file
     *
     * @return Depth to Output - true or false
     */
    public Boolean depthToOutput() {
        return myDepthOutput.value();
    }

    /**
     * User sets true or false, indicating if they do
     * or do not want depth information written to the
     * HDF5 file.
     *
     * @param value Write depth to output file
     *
     * @return this plugin
     */
    public RGBSProductionSNPCallerPlugin depthToOutput(Boolean value) {
        myDepthOutput = new PluginParameter<>(myDepthOutput, value);
        return this;
    }
    /**
     * Maximum Tag Length
     *
     * @return Maximum Tag Length
     */
    public Integer kmerLength() {
        return myKmerLength.value();
    }

    /**
     * Set Maximum Tag Length:  User should set this value
     * equivalent to what was used in GBSSeqToTagDBPlugin
     * for maximum tag length when creating the database.
     * If the two values are not equal inconsistent results
     * may occur.
     *
     * @param value Maximum Tag Length
     *
     * @return this plugin
     */
    public RGBSProductionSNPCallerPlugin kmerLength(Integer value) {
        myKmerLength = new PluginParameter<>(myKmerLength, value);
        return this;
    }
    /**
     *  Minimum Position Quality Score
     *
     * @return Minimum position quality score
     */
    public Double positionQualityScore() {
        return posQualityScore.value();
    }

    /**
     * Set Minimum quality score for position:  This value is used to pull
     * SNPs out of the snpposition table.  Only snps with quality
     * scores meeting or exceeding the specified value will be 
     * processed.
     *
     * @param value Minimum position quality score
     *
     * @return this plugin
     */
    public RGBSProductionSNPCallerPlugin positionQualityScore(Double value) {
        posQualityScore = new PluginParameter<>(posQualityScore, value);
        return this;
    }
    
    /**
     *  Batch size for processing fastq files
     *
     * @return batchSize
     */
    public Integer batchSize() {
        return myBatchSize.value();
    }
    /**
     * Set number of Fastq files processed simultaneously
     * @param value
     * @return
     */
    public RGBSProductionSNPCallerPlugin batchSize(Integer value) {
        myBatchSize = new PluginParameter<>(myBatchSize, value);
        return this;
    }
    /**
     * Minimum quality score within the barcode and read length
     * to be accepted
     *
     * @return Minimum quality score
     */
    public Integer minimumQualityScore() {
        return myMinQualScore.value();
    }

    /**
     * Set Minimum quality score. Minimum quality score within
     * the barcode and read length to be accepted
     *
     * @param value Minimum quality score
     *
     * @return this plugin
     */
    public RGBSProductionSNPCallerPlugin minimumQualityScore(Integer value) {
        myMinQualScore = new PluginParameter<>(myMinQualScore, value);
        return this;
    }
    /**
     * Use STACKS likelihood method to call heterozygotes (default: use
     * tasselGBS likelihood ratio method)
     *
     * @return Use Stacks Likelihood
     */
    //public Boolean useStacksLikelihood() {
    //    return myStacksLikelihood.value();
    //}
    /**
     * Set Use Stacks Likelihood. Use STACKS likelihood method to call
     * heterozygotes (default: use tasselGBS likelihood ratio method)
     *
     * @param value Use Stacks Likelihood
     *
     * @return this plugin
     */
    //public ProductionSNPCallerPlugin useStacksLikelihood(Boolean value) {
    //    myStacksLikelihood = new PluginParameter<>(myStacksLikelihood, value);
    //    return this;
    //}
}
