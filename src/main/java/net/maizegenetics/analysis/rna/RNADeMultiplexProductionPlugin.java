package net.maizegenetics.analysis.rna;

import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import net.maizegenetics.analysis.gbs.v2.GBSUtils;
import net.maizegenetics.analysis.gbs.v2.ProductionSNPCallerPluginV2;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.tag.*;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import net.maizegenetics.analysis.gbs.v2.EnzymeList;


public class RNADeMultiplexProductionPlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(RNADeMultiplexProductionPlugin.class);

    private PluginParameter<String> myInputDir = new PluginParameter.Builder<>("i", null, String.class).guiName("Input Directory").required(true).inDir()
            .description("Input directory containing fastq AND/OR qseq files.").build();
    private PluginParameter<String> myKeyFile = new PluginParameter.Builder<>("k", null, String.class).guiName("Key File").required(true).inFile()
            .description("Key file listing barcodes distinguishing the samples").build();
    private PluginParameter<String> myInputDB = new PluginParameter.Builder<>("db", null, String.class).guiName("Input GBS Database").required(true).inFile()
            .description("Input Database file if using SQLite").build();
    private PluginParameter<FindMatchByWordHash.MatchType> myMatchingType = new PluginParameter.Builder<>("matchType", null, FindMatchByWordHash.MatchType.class)
            .guiName("Matching approach").required(true)
            .range(Range.encloseAll(Arrays.asList(FindMatchByWordHash.MatchType.values())))
            .description("Approach used for matching the reads to the contigs").build();
    private PluginParameter<Integer> myWordSizeMatching = new PluginParameter.Builder<>("word", 16, Integer.class).guiName("Word size used by match")
            .description("Word size used to find matches with reference sequences").build();
    private PluginParameter<Integer> myMaxWordRepeats = new PluginParameter.Builder<>("wordRep", 10, Integer.class).guiName("Maximum repeats of word")
            .description("Maximum repetitiveness of the word").build();
    private PluginParameter<Boolean> mySearchReverseComplement = new PluginParameter.Builder<>("searchRevComp",true,Boolean.class).guiName("Search reverse complement")
            .description("Search of the reverse complements of the reference sequences").build();
    private PluginParameter<Integer> myBatchSize = new PluginParameter.Builder<>("batchSize", 8, Integer.class).guiName("Batch size of fastq files").required(false)
            .description("Number of flow cells being processed simultaneously").build();
    private PluginParameter<Integer> myMinQualScore = new PluginParameter.Builder<>("mnQS", 0, Integer.class).guiName("Minimum quality score").required(false)
            .description("Minimum quality score within the barcode and read length to be accepted").build();



    private TagDataWriter tdw = null;
    private Map<String,int[]> taxatissueCntMap;
    private FindMatchByWordHash findMatchByWordHash;
 
    protected static int readEndCutSiteRemnantLength;
    private static String myEnzyme="ignore";
    private static int maxTaxa;
    private static int maxTissue;

    private boolean taglenException;
    
    public RNADeMultiplexProductionPlugin() {
        super(null, false);
    }

    public RNADeMultiplexProductionPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public void postProcessParameters() {
        EnzymeList.Enzyme enzyme = EnzymeList.defaultCache.getEnzyme(myEnzyme);
        readEndCutSiteRemnantLength = enzyme.readEndCutSiteRemnantLength;
    }

    @Override
    public DataSet processData(DataSet input) {
        int batchSize = batchSize();
        List<Path> directoryFiles= DirectoryCrawler.listPaths(GBSUtils.inputFileGlob, Paths.get(myInputDir.value()).toAbsolutePath());
        if(directoryFiles.isEmpty()) {
            myLogger.warn("No files matching:"+GBSUtils.inputFileGlob);
            return null;
        }
        List<Path> inputSeqFiles =directoryFiles;
//        if (inputSeqFiles.size() == 0) return null; // no files to process

        tdw =new TagDataSQLite(myInputDB.value());



        TaxaList masterTaxaList= TaxaListIOUtils.readTaxaAnnotationFile(keyFile(), GBSUtils.sampleNameField, new HashMap<>(), true);
        Multimap<String,Taxon> temp=TaxaListIOUtils.getMapOfTaxonByAnnotation(masterTaxaList, "FileName");
        Map<String,Taxon> fileToTaxonMap=new HashMap<>();
        temp.entries().stream().forEach(entries -> {
            Taxon t=fileToTaxonMap.put(entries.getKey(), entries.getValue());
            if(t!=null) System.err.println("Duplicate file");
        });  //wrote over duplicated

        ArrayList<String> masterTissueList = (ArrayList)TaxaListIOUtils.readTissueAnnotationFile(keyFile(), GBSUtils.tissueNameField);
        // Write tissue list to the db
        tdw.putTaxaList(masterTaxaList);
        tdw.putAllTissue(masterTissueList);

        maxTaxa = masterTaxaList.size();
        maxTissue = masterTissueList.size();

        //get all the possible tags (=sequences = contigs)
        findMatchByWordHash = FindMatchByWordHash.getBuilder(tdw.getTags())
                .matchType(matchingType())
                .wordLength(wordSizeMatching())
                .maxWordCopies(maxWordRepeats())
                .searchBiDirectional(searchReverseComplement())
                .build();
        taxatissueCntMap=new ConcurrentHashMap<>(maxTissue*maxTaxa*2);

        int batchNum = inputSeqFiles.size()/batchSize;
       
        if (inputSeqFiles.size() % batchSize !=0) batchNum++;
        System.out.println("RNADeMultiplexProductionPlugin: Total batches to process: " + batchNum);


        taglenException = false;
        for (int idx = 0; idx < inputSeqFiles.size(); idx+=batchSize) {
            int end = idx+batchSize;
            if (end > inputSeqFiles.size()) end = inputSeqFiles.size();
            ArrayList<Path> sub = new ArrayList<Path>();
            for (int jdx = idx; jdx < end; jdx++) sub.add(inputSeqFiles.get(jdx));
            System.out.println("\nStart processing batch " + String.valueOf(idx/batchSize+1));
            sub.parallelStream()
            .forEach(inputSeqFile -> {
                try {
//                    System.err.println(inputSeqFile);
//                    System.err.println(fileToTaxonMap.get(inputSeqFile.getFileName().toString()));
                    if(fileToTaxonMap.get(inputSeqFile.getFileName().toString())==null) return;
                    processFastQ(inputSeqFile,fileToTaxonMap.get(inputSeqFile.getFileName().toString()),minQualScore());

                } catch (StringIndexOutOfBoundsException oobe) {
                    oobe.printStackTrace();
                    myLogger.error(oobe.getMessage());
                    setTagLenException();
                    return;
                }              
            });
            if (taglenException == true) return null; // Tag length failure from processFastQ - halt processing
         
            System.out.println("\nFinished processing batch " + String.valueOf(idx/batchSize+1));
        }
        
        // put counts to database.
        taxatissueCntMap.entrySet().stream().forEach(entry -> {
            String[] taxaTissueNames = entry.getKey().split(",");
            tdw.putTaxaTissueDistribution(taxaTissueNames[0],taxaTissueNames[1], findMatchByWordHash.getTags(),entry.getValue());
        });

//        writeReadsPerSampleReports(tagTaxaTissueMap.size());
        return new DataSet(new Datum("DB Name",inputDB(),""),null);
    }


    private void processFastQ(Path fastqFile, Taxon taxon, int minQual) throws StringIndexOutOfBoundsException {
        String tissuetaxaKey=taxon.getName()+","+taxon.getAnnotation().getTextAnnotation("Tissue")[0];
        taxatissueCntMap.putIfAbsent(tissuetaxaKey, new int[findMatchByWordHash.totalNumberOfTags()]);
        int allReads=0, readsMatchedToContigs = 0, lowQualityReads = 0;
        int shortReads = 0;
        try {
            int qualityScoreBase=GBSUtils.determineQualityScoreBase(fastqFile);
            BufferedReader br = Utils.getBufferedReader(fastqFile.toString(), 1 << 22);
            long time=System.nanoTime();
            String[] seqAndQual;
            String likelyReadEnd = "AGATCGGA";
            while ((seqAndQual=GBSUtils.readDeMultiPlexFastQBlock(br, allReads)) != null) {
                allReads++;
                int adapterStart = 0;
                // Decode barcode using the current sequence & quality  score

                // Find barcode for smallRNA (taxa comes from barcode)
                // Strip off the barcode - original sequence did NOT have it,
                // this was added when we read the fastQ block.
                // We to revert the sequence to what is in the fastQ file
                // so the quality scores match up correctly.
                String sequence = seqAndQual[0];


                adapterStart = sequence.indexOf(likelyReadEnd);
                if (adapterStart > 0) {
                    sequence = sequence.substring(0, adapterStart-1);
                }

                if(minQual>0) {
                    //todo move getFirstLowQualityPos into this class?
                    if(BaseEncoder.getFirstLowQualityPos(seqAndQual[1],minQual, qualityScoreBase) < sequence.length()){
                        lowQualityReads++;
                        continue;
                    }
                }
                FindMatchByWordHash.Match hitIndex= findMatchByWordHash.match(sequence);
                if(!hitIndex.isEmpty()) {
                    readsMatchedToContigs++;
                    taxatissueCntMap.get(tissuetaxaKey)[hitIndex.tagIndex()]++;
                } else {

                   // if(allReads%10000==0) System.out.println(">"+allReads+"\n"+seqAndQual[0]+"\n"+sequence);
                }

            }
            myLogger.info("Total number of reads in FASTQ file=" + allReads  + " rate:" + (System.nanoTime()-time)/allReads +" ns/read");
            myLogger.info("Total number of reads matched to <= " + maxWordRepeats() + " contigs/genes=" + readsMatchedToContigs + " (adjust upper limit using maxWordRepeats)");
            myLogger.info("Total number of rejected low quality reads and reads with non-AGCT bases=" + lowQualityReads);
            myLogger.info("Timing process (sorting, collapsing, and writing TagCount to file).");
            myLogger.info("Process took " + (System.nanoTime() - time)/1e9 + " seconds for file " + fastqFile.toString());
            br.close();
        } catch (Exception e) {
            myLogger.error("Matched reads: " + readsMatchedToContigs);
            e.printStackTrace();
        }
    }

    public void setTagLenException() {
        taglenException = true;
    }
    

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "RNA Production Counter";
    }

    @Override
    public String getToolTipText() {
        return "RNA Production Counter";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
     public static void main(String[] args) {
         GeneratePluginCode.generate(RNADeMultiplexProductionPlugin.class);
     }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public String runPlugin(DataSet input) {
        return (String) performFunction(input).getData(0).getData();
    }

    /**
     * Input directory containing fastq AND/OR qseq files.
     *
     * @return Input Directory
     */
    public String inputDir() {
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
    public RNADeMultiplexProductionPlugin inputDir(String value) {
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
    public RNADeMultiplexProductionPlugin keyFile(String value) {
        myKeyFile = new PluginParameter<>(myKeyFile, value);
        return this;
    }

    /**
     * Input Database file if using SQLite
     *
     * @return Input GBS Database
     */
    public String inputDB() {
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
    public RNADeMultiplexProductionPlugin inputDB(String value) {
        myInputDB = new PluginParameter<>(myInputDB, value);
        return this;
    }

    /**
     * Approach used for matching the reads to the contigs
     *
     * @return Matching approach
     */
    public FindMatchByWordHash.MatchType matchingType() {
        return myMatchingType.value();
    }

    /**
     * Set Matching approach. Approach used for matching the
     * reads to the contigs
     *
     * @param value Matching approach
     *
     * @return this plugin
     */
    public RNADeMultiplexProductionPlugin matchingType(FindMatchByWordHash.MatchType value) {
        myMatchingType = new PluginParameter<>(myMatchingType, value);
        return this;
    }

    /**
     * Word size used to find matches with reference sequences
     *
     * @return Word size used by match
     */
    public Integer wordSizeMatching() {
        return myWordSizeMatching.value();
    }

    /**
     * Set Word size used by match. Word size used to find
     * matches with reference sequences
     *
     * @param value Word size used by match
     *
     * @return this plugin
     */
    public RNADeMultiplexProductionPlugin wordSizeMatching(Integer value) {
        myWordSizeMatching = new PluginParameter<>(myWordSizeMatching, value);
        return this;
    }

    /**
     * Maximum repetitiveness of the word
     *
     * @return Maximum repeats of word
     */
    public Integer maxWordRepeats() {
        return myMaxWordRepeats.value();
    }

    /**
     * Set Maximum repeats of word. Maximum repetitiveness
     * of the word
     *
     * @param value Maximum repeats of word
     *
     * @return this plugin
     */
    public RNADeMultiplexProductionPlugin maxWordRepeats(Integer value) {
        myMaxWordRepeats = new PluginParameter<>(myMaxWordRepeats, value);
        return this;
    }

    /**
     * Search of the reverse complements of the reference
     * sequences
     *
     * @return Search reverse complement
     */
    public Boolean searchReverseComplement() {
        return mySearchReverseComplement.value();
    }

    /**
     * Set Search reverse complement. Search of the reverse
     * complements of the reference sequences
     *
     * @param value Search reverse complement
     *
     * @return this plugin
     */
    public RNADeMultiplexProductionPlugin searchReverseComplement(Boolean value) {
        mySearchReverseComplement = new PluginParameter<>(mySearchReverseComplement, value);
        return this;
    }

    /**
     * Number of flow cells being processed simultaneously
     *
     * @return Batch size of fastq files
     */
    public Integer batchSize() {
        return myBatchSize.value();
    }

    /**
     * Set Batch size of fastq files. Number of flow cells
     * being processed simultaneously
     *
     * @param value Batch size of fastq files
     *
     * @return this plugin
     */
    public RNADeMultiplexProductionPlugin batchSize(Integer value) {
        myBatchSize = new PluginParameter<>(myBatchSize, value);
        return this;
    }

    /**
     * Minimum quality score within the barcode and read length
     * to be accepted
     *
     * @return Minimum quality score
     */
    public Integer minQualScore() {
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
    public RNADeMultiplexProductionPlugin minQualScore(Integer value) {
        myMinQualScore = new PluginParameter<>(myMinQualScore, value);
        return this;
    }

}
