package net.maizegenetics.analysis.gbs.repgen;

import java.awt.Frame;
import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javax.swing.ImageIcon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.analysis.gbs.v2.GBSUtils;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TaxaDistBuilder;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;

/**
 * Develops a discovery rGBS database based on a folder of sequencing files
 *
 * Keeps only good reads having no N's in the
 * useful part of the sequence. Trims off the barcodes and truncates sequences
 * that (1) have a second cut site, or (2) read into the common adapter.
 * 
 * Originally the reference throughout was to "tag". This is being changed
 * to "kmer" as the pipeline is a kmer alignment process.  
 *
 * @author Ed Buckler
 */
public class RepGenLoadSeqToDBPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(RepGenLoadSeqToDBPlugin.class);

    private PluginParameter<String> myInputDir = new PluginParameter.Builder<>("i", null, String.class).guiName("Input Directory").required(true).inDir()
            .description("Input directory containing FASTQ files in text or gzipped text.\n"
                    + "     NOTE: Directory will be searched recursively and should\n"
                    + "     be written WITHOUT a slash after its name.").build();
    private PluginParameter<String> myKeyFile = new PluginParameter.Builder<>("k", null, String.class).guiName("Key File").required(true).inFile()
            .description("Key file listing barcodes distinguishing the samples").build();
    private PluginParameter<Integer> myKmerLength = new PluginParameter.Builder<>("kmerLength", 150, Integer.class).guiName("Maximum Kmer Length")
            .description("Specified length for each kmer to process").build();
    private PluginParameter<Integer> myMinKmerLength = new PluginParameter.Builder<>("minKmerL", 120, Integer.class).guiName("Minimum Kmer Length")
            .description("Minimum kmer Length after second cut site is removed").build();
    private PluginParameter<Integer> myMinKmerCount = new PluginParameter.Builder<>("c", 10, Integer.class).guiName("Min Kmer Count")
            .description("Minimum kmer count").build();
    private PluginParameter<Integer> minTaxa = new PluginParameter.Builder<>("minTaxa",2, Integer.class).guiName("Min Taxa Represented")
            .description("Minimum numer of taxa where kmer is found").build();
    private PluginParameter<String> myOutputDB = new PluginParameter.Builder<>("db", null, String.class).guiName("Output Database File").required(true).outFile()
            .description("Output Database File").build();
    private PluginParameter<Integer> myMinQualScore = new PluginParameter.Builder<>("mnQS", 0, Integer.class).guiName("Minimum quality score").required(false)
            .description("Minimum quality score within the barcode and read length to be accepted").build();
    private PluginParameter<Integer> myMaxKmerNumber = new PluginParameter.Builder<>("mxKmerNum", 50000000, Integer.class).guiName("Maximum Kmer Number").required(false)
            .description("Maximum number of kmers").build();
    private PluginParameter<Integer> myBatchSize = new PluginParameter.Builder<>("batchSize", 100, Integer.class).guiName("Batch size of fastq files").required(false)
            .description("Number of flow cells being processed simultaneously").build();
//    private PluginParameter<Boolean> myDeleteOldData = new PluginParameter.Builder<Boolean>("deleteOldData",true,Boolean.class).guiName("Delete Old Data")
//            .description("Delete existing SNP quality data from db tables").build();
    LongAdder roughTagCnt = new LongAdder();

    private TagDistributionMap tagCntMap;
    private TagCountQualityScoreMap tagCntQSMap;
    private boolean taglenException;
    protected static int readEndCutSiteRemnantLength;
    protected static int qualityScoreBase;
    String[] likelyReadEndStrings;
    
    public RepGenLoadSeqToDBPlugin() {
        super(null, false);
    }

    public RepGenLoadSeqToDBPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }


    private long[] calcTagMapStats(TagDistributionMap tagCntMap) {
        long totalDepth=0, memory=0;
        int cnt=0;
        for (Map.Entry<Tag, TaxaDistribution> entry : tagCntMap.entrySet()) {
            memory+=entry.getValue().memorySize();
            memory+=25;
            totalDepth+=entry.getValue().totalDepth();
            cnt++;
        }
        int currentSize = tagCntMap.size();
        memory+=tagCntMap.size()*2*16;  //estimate for the map size
        long[] stats={currentSize,memory, totalDepth,totalDepth/cnt};
        System.out.printf("Map Tags:%,d  Memory:%,d  TotalDepth:%,d  AvgDepthPerTag:%d%n",stats[0],stats[1],stats[2],stats[3]);
        return stats;
    }
    @Override
    public void postProcessParameters() {

    }
    @Override
    public DataSet processData(DataSet input) {
        int batchSize = myBatchSize.value();
        float loadFactor = 0.95f; //TODO: Seems high to Ed
        tagCntMap = new TagDistributionMap (myMaxKmerNumber.value(),loadFactor, 128, this.minKmerCount());
        tagCntQSMap = new TagCountQualityScoreMap(myMaxKmerNumber.value(),loadFactor, 128, this.minKmerCount()); // must add Collections.synchonizedList<String> as the list! 
        try {
            TaxaList masterTaxaList= TaxaListIOUtils.readTaxaAnnotationFile(keyFile(), GBSUtils.sampleNameField, new HashMap<>(), true);
            Map<String,Taxon> fileTaxaMap=TaxaListIOUtils.getUniqueMapOfTaxonByAnnotation(masterTaxaList,GBSUtils.fileNameField)
                    .orElseThrow(() -> new IllegalArgumentException("Error: Same file points more than one taxon in the KeyFile"));
            //Get the list of fastq files
            List<Path> directoryFiles= DirectoryCrawler.listPaths(GBSUtils.inputFileGlob, Paths.get(myInputDir.value()).toAbsolutePath());
            if(directoryFiles.isEmpty()) {
                myLogger.warn("No files matching:"+ GBSUtils.inputFileGlob);
                return null;
            }
            // Cull files that are not represented in the given key file
            List<Path> inputSeqFiles =  directoryFiles.stream()
                    .peek(path -> System.out.println(path.getFileName().toString()))
                    .filter(path -> fileTaxaMap.containsKey(path.getFileName().toString()))
                    .collect(Collectors.toList());
            if (inputSeqFiles.size() == 0) {
                System.out.println("RepGenLoadSeqToDB:processData - found NO files represented in key file.");
                System.out.println("Please verify your file names are formatted correctly and that your key file contains the required headers.");
                return null; // no files in this directory to process
            }
            System.out.println("Found " + inputSeqFiles.size() + " files to process");
            //Files in a batch have roughly the same size
            //inputSeqFiles = this.sortFastqBySize(inputSeqFiles);
            int batchNum = inputSeqFiles.size()/batchSize;
            if (inputSeqFiles.size()%batchSize != 0) batchNum++;

            // ALways deleting old DB - need to work on code to re-create
            // array of tag quality scores if we add to old DB.
            RepGenDataWriter tdw = null;
            if (Files.exists(Paths.get(myOutputDB.value()))) {
                try {
                    Files.delete(Paths.get(myOutputDB.value()));
                } catch (Exception exc){
                    System.out.println("Error when trying to delete database file: " + myOutputDB.value());
                    System.out.println("File delete error: " + exc.getMessage());
                    return null;
                }
            }

//                if (deleteOldData()) {
// 
//                }
//                } else {
//                 // We'll append data to new DB if all new taxa are contained in db
//                    tdw=new RepGenSQLite(myOutputDB.value());
//                    TaxaList oldTaxaList = tdw.getTaxaList();
//                    boolean sameMasterTaxaList = true;
//                    Taxon badTaxon = null;
//                    for (Taxon taxa : masterTaxaList) {
//                        if (!oldTaxaList.contains(taxa)) {
//                            sameMasterTaxaList = false;
//                            badTaxon = taxa;
//                            break;
//                        }
//                    }
//                    if (!sameMasterTaxaList) {
//                        myLogger.error("Taxa list in keyfile contains taxa not in db:  " + badTaxon +
//                            ".\nEither delete existing DB, use a new DB output file, or use keyfile with entries matching existing taxa list.\n");
//                        ((TagDataSQLite)tdw).close();
//                        return null;
//                    }
//                    // Grab existing data from db, append to empty tagCntMap
//                    Map<Tag, TaxaDistribution> existingTDM = tdw.getAllTagsTaxaMap();
//                    tagCntMap.putAll(existingTDM);
//                    tdw.clearTagTaxaDistributionData(); // clear old data - it will be re-added at the end.
//                }
//            }

            if (tdw == null) tdw=new RepGenSQLite(myOutputDB.value());
            taglenException = false;
            for (int i = 0; i < inputSeqFiles.size(); i+=batchSize) {
                int end = i+batchSize;
                if (end > inputSeqFiles.size()) end = inputSeqFiles.size();
                ArrayList<Path> sub = new ArrayList<Path>();
                for (int j = i; j < end; j++) sub.add(inputSeqFiles.get(j));
                System.out.println("\nStart processing batch " + String.valueOf(i/batchSize+1));
                sub//.parallelStream()
                    .forEach(inputSeqFile -> {
                        try {
                            int taxaIndex=masterTaxaList.indexOf(fileTaxaMap.get(inputSeqFile.getFileName().toString()));
                            processFastQ(inputSeqFile, taxaIndex, masterTaxaList, tagCntMap, kmerLength(), minimumQualityScore());
                        } catch (StringIndexOutOfBoundsException oobe) {
                            oobe.printStackTrace();
                            myLogger.error(oobe.getMessage());
                            setTagLenException();
                            return;
                        }
                });
                if (taglenException == true) {
                    ((TagDataSQLite)tdw).close();
                    return null; // Tag length failure from processFastQ - halt processing
                }

                System.out.println("\nKmers are added from batch "+String.valueOf(i/batchSize+1) + ". Total batch number: " + batchNum);
                int currentSize = tagCntMap.size();
                System.out.println("Current number: " + String.valueOf(currentSize) + ". Max kmer number: " + String.valueOf(myMaxKmerNumber.value()));
                System.out.println(String.valueOf((float)currentSize/(float)myMaxKmerNumber.value()) + " of max tag number");

                if (currentSize > 0) { // calcTagMapStats() gets "divide by 0" error when size == 0
                    this.calcTagMapStats(tagCntMap);
                    System.out.println();
                    //make sure don't lose rare ones, need to set maxTagNumber large enough
                    System.out.println("BEFORE removeTagsWihtoutReplication, tagCntMap.size= " + tagCntMap.keySet().size()
                      + ", tagCntQSMap.size= " + tagCntQSMap.keySet().size());
                    removeTagsWithoutReplication(tagCntMap,tagCntQSMap,2,4);  //TODO lows values that keep the map under control
                    if (tagCntMap.size() == 0) {
                        System.out.println("WARNING:  After removing tags without replication, there are NO  tags left in the database");
                    } else {
                        this.calcTagMapStats(tagCntMap);
                        System.out.println();
                        System.out.println("After removeTagsWithoutReplication: Kmer number is reduced to " + tagCntMap.size()+"\n");
                    }
                    this.roughTagCnt.reset();
                    this.roughTagCnt.add(tagCntMap.size());
                } else {
                    System.out.println("WARNING: Current tagcntmap size is 0 after processing batch " + String.valueOf(i/batchSize+1) );
                }

                System.out.println("Total memory: "+ String.valueOf((double)(Runtime.getRuntime().totalMemory()/1024/1024/1024))+" Gb");
                System.out.println("Free memory: "+ String.valueOf((double)(Runtime.getRuntime().freeMemory()/1024/1024/1024))+" Gb");
                System.out.println("Max memory: "+ String.valueOf((double)(Runtime.getRuntime().maxMemory()/1024/1024/1024))+" Gb");
                System.out.println("\n");
            }
            System.out.println("\nAll the batch are processed");
            removeTagsWithoutReplication(tagCntMap,tagCntQSMap,minTaxa(),minKmerCount());
//            tagCntMap.removeTagByCount(myMinKmerCount.value());
//            tagCntQSMap.removeTagByCount(myMinKmerCount.value()); // should remove the same tags as in call above
           // System.out.println("By removing kmers with minCount of " + myMinKmerCount.value() + "Kmer number is reduced to " + tagCntMap.size()+"\n");
            System.out.printf("Filter by kmerCnt:%d taxaCount:%d results in MapSize:%d\n",minKmerCount(), minTaxa(), tagCntMap.size());

            // Now create map of count/qualityString
            // THis must have the same number of keys as does tagCntMap.
            
            if (tagCntMap.keySet().size() != tagCntQSMap.keySet().size()) {
                System.out.println("Mismatch between size of tagCntMap " + tagCntMap.keySet().size() 
                        + ", and tagCntQSMap " + tagCntQSMap.keySet().size());
                System.out.println("quitting - nothing added to DB!");
                ((RepGenSQLite)tdw).close();
                return null;    
            } else {
                System.out.println("\ntagCntMap and tagCntQSMap have same size: " + tagCntMap.keySet().size());
            }
            
            Map<Tag,Tuple<Integer,String>> tagInstanceAverageQS = calculateTagAveQS(tagCntQSMap,qualityScoreBase);
            System.out.println("Before add to DB: sizeof tagInstanceAverageQS.keySet = " + tagInstanceAverageQS.keySet().size());
            tdw.putTaxaList(masterTaxaList);
            // 2 fields to putALlTag as some callers of this method don't hvae the qs/num-instances data
            //tdw.putAllTag(tagCntMap.keySet(), null);
            tdw.putAllTag(tagCntMap.keySet(),tagInstanceAverageQS); // includes tag, the number of instances of this tag, and the average quality string
            tdw.putTaxaDistribution(tagCntMap);
            ((RepGenSQLite)tdw).close();  //todo autocloseable should do this but it is not working.
        } catch(Exception e) {
            e.printStackTrace();
        }
        return new DataSet(new Datum("TagMap",tagCntMap,""),this);
    }

    private void processFastQ(Path fastqFile, int taxaIndex, TaxaList masterTaxaList,
                              TagDistributionMap masterTagTaxaMap, int preferredTagLength, int minQual) throws StringIndexOutOfBoundsException{
        int allReads=0, goodBarcodedReads = 0, lowQualityReads = 0, tooShortReads = 0;
        int maxTaxaNumber=masterTaxaList.size();
        int checkSize = 10000000;
        myLogger.info("processing file " + fastqFile.toString());
        try {
            //int qualityScoreBase=GBSUtils.determineQualityScoreBase(fastqFile);
            qualityScoreBase=GBSUtils.determineQualityScoreBase(fastqFile);
            BufferedReader br = Utils.getBufferedReader(fastqFile.toString(), 1 << 22);
            long time=System.nanoTime();
            String[] seqAndQual;

            while ((seqAndQual=GBSUtils.readFastQBlock(br,allReads)) != null) {
                allReads++;
                int tagEnd = seqAndQual[0].length() < preferredTagLength ? seqAndQual[0].length() : preferredTagLength;
                // Check for and toss low quality reads
                if(minQual>0) {
                    //todo move getFirstLowQualityPos into this class?
                    if(BaseEncoder.getFirstLowQualityPos(seqAndQual[1],minQual, qualityScoreBase)<(tagEnd)){
                    	lowQualityReads++;
                    	continue;
                    }
                }

                // preferredTagLength is actually maximum tag length.  We need to check if 
                // tag is less than minimumKmerLength, toss that.  And truncate if is greater
                // than kmerLength.
                if (seqAndQual[0].length() < minimumKmerLength()) {
                //if (seqAndQual[0].length() < preferredTagLength) {
                    //System.out.println("Tossing kmer, length too short: " + seqAndQual[0].length());
//                	String errMsg = "\n\nERROR processing " + fastqFile.toString() + "\n" +
//                			"Reading entry number " + allReads + " fails the length test.\n" +
//                			"Sequence length " + seqAndQual[0].length() +
//                			" is less then maxKmerLength " + minimumKmerLength() + ".\n" +
//                			"Re-run your files with either a shorter mxKmerL value or a higher minimum quality score.\n";
                	//throw new StringIndexOutOfBoundsException(errMsg);
                    tooShortReads++;
                	continue; // for combined reads, just continue, don't error out.
                }
                //Tag tag = TagBuilder.instance(seqAndQual[0].substring(0,preferredTagLength)).build();

                Tag tag = TagBuilder.instance(seqAndQual[0].substring(0,tagEnd)).build();
                

                if(tag==null) continue;   //null occurs when any base was not A, C, G, T
                goodBarcodedReads++;
                TaxaDistribution taxaDistribution=masterTagTaxaMap.get(tag);
                if(taxaDistribution==null) {
                    masterTagTaxaMap.put(tag,TaxaDistBuilder.create(maxTaxaNumber,taxaIndex));
                    this.roughTagCnt.increment();
                } else {
                    taxaDistribution.increment(taxaIndex);
                }
                // add to qualityScoreMap
               // String tagQS = seqAndQual[1].substring(0,preferredTagLength);
                String tagQS = seqAndQual[1].substring(0,tagEnd);
                List<String> tagScores = tagCntQSMap.get(tag);
                if (tagScores==null) {
                    // create new synchronizedList
                    tagScores =  Collections.synchronizedList(new ArrayList<String>());
                    tagScores.add(tagQS);
                    tagCntQSMap.put(tag, tagScores);
                } else {
                    tagScores.add(tagQS);
                }
                if (allReads % checkSize == 0) {
                    myLogger.info("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads
                            + " rate:" + (System.nanoTime()-time)/allReads +" ns/read. Current tag count:" + this.roughTagCnt);
                }
            }
            myLogger.info("Summary for "+fastqFile.toString()+"\n"+
                    "Total number of reads in lane=" + allReads +"\n"+
                    "Total number of good barcoded reads=" + goodBarcodedReads+"\n"+
                    "Total number of low quality reads=" + lowQualityReads+"\n"+
                    "Total number of too short reads=" + tooShortReads + "\n"+
                    "Timing process (sorting, collapsing, and writing TagCount to file)."+"\n"+
                    "Process took " + (System.nanoTime() - time)/1e6 + " milliseconds.");
            System.out.println("tagCntMap size: "+masterTagTaxaMap.size());
            br.close();
        } catch (StringIndexOutOfBoundsException oobe) {
        	throw oobe; // pass it up to print error and stop processing
        } catch (Exception e) {
            myLogger.error("Good Barcodes Read: " + goodBarcodedReads);

            e.printStackTrace();
        }
    }


    /**
     * This method removes all tags are are never repeated in a single sample (taxa).  The concept is that
     * all biologically real tag should show up twice somewhere.  This could be called at the end of every
     * flowcell to test all the novel tags.
     * 
     * In addition, it removes the corresponding entry from the tag-qualityscore map.  They must remain in synch
     */
    private static void removeTagsWithoutReplication (TagDistributionMap masterTagTaxaMap, TagCountQualityScoreMap tagCntQSMap,int minTaxa, int minDepth) {
        int currentSize = masterTagTaxaMap.size();
       // int minTaxa=2;
        System.out.println("Starting removeTagsWithoutReplication. Current tag number: " + currentSize);
        LongAdder tagsRemoved=new LongAdder();
        masterTagTaxaMap.entrySet().parallelStream().forEach(t -> {
            TaxaDistribution td = t.getValue();
            if(td.totalDepth()<minDepth) {  //depth
                masterTagTaxaMap.remove(t.getKey());
                tagCntQSMap.remove(t.getKey());
                tagsRemoved.increment();
            }
            else if(IntStream.of(td.depths()).filter(depth -> depth>1).count()<minTaxa) {  //taxa test
                masterTagTaxaMap.remove(t.getKey());
                tagCntQSMap.remove(t.getKey());
                tagsRemoved.increment();
            }
        });
        System.out.println("Finished removeTagsWithoutReplication.  tagsRemoved = " + tagsRemoved + ". Current tag number: " + String.valueOf(currentSize-tagsRemoved.intValue()));
    }

    public void setTagLenException() {
    	taglenException = true;
    }
    
    public static Map<Tag,Tuple<Integer,String>> calculateTagAveQS(TagCountQualityScoreMap tagCntQSMap, int qualityScoreBase) {
        
        Map<Tag,Tuple<Integer,String>> tagInstanceAverageQS = new ConcurrentHashMap<Tag,Tuple<Integer,String>>();
        
        // the number of instances for each tag equals the number of values on the arraylist
        // iterate through the list, add values to the tagInstanceAverageQS map
        tagCntQSMap.entrySet().parallelStream().forEach(entry -> {
            Tag tag = entry.getKey();
            List<String> scores = entry.getValue();
            // add all the scores, divide by number of scores, convert back to string
            int numInstances = scores.size();
            int scoreLen = scores.get(0).length();
            // the scores should all be the same length as the tag
            int[] scoreArray = new int[scoreLen];
            
            for (int idx = 0; idx < numInstances; idx++) {
                String currentScore = scores.get(idx);
                for (int jdx = 0; jdx < scoreLen; jdx++) {
                    scoreArray[jdx] += (currentScore.charAt(jdx) - qualityScoreBase);
                }
            }
            // positions for all scores are added.  now divide by number of scores to get
            // the average, translate value back to string
            StringBuilder sb = new StringBuilder();
            for (int jdx = 0; jdx < scoreLen; jdx++) {
                scoreArray[jdx] /= numInstances;
                char ch = (char)(scoreArray[jdx] + qualityScoreBase);
                sb.append(ch);
            }
            // Store to map
            Tuple<Integer,String> instanceScore = new Tuple<Integer,String>(numInstances,sb.toString());
            tagInstanceAverageQS.put(tag, instanceScore);           
        });
        
        return tagInstanceAverageQS;
    }

// The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
//     public static void main(String[] args) {
//         GeneratePluginCode.generate(RepGenLoadSeqToDBPlugin.class);
//     }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
//    public <Type> runPlugin(DataSet input) {
//        return (<Type>) performFunction(input).getData(0).getData();
//    }

    /**
     * Input directory containing FASTQ files in text or gzipped
     * text.
     *      NOTE: Directory will be searched recursively and
     * should
     *      be written WITHOUT a slash after its name.
     *
     * @return Input Directory
     */
    public String inputDirectory() {
        return myInputDir.value();
    }

    /**
     * Set Input Directory. Input directory containing FASTQ
     * files in text or gzipped text.
     *      NOTE: Directory will be searched recursively and
     * should
     *      be written WITHOUT a slash after its name.
     *
     * @param value Input Directory
     *
     * @return this plugin
     */
    public RepGenLoadSeqToDBPlugin inputDirectory(String value) {
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
    public RepGenLoadSeqToDBPlugin keyFile(String value) {
        myKeyFile = new PluginParameter<>(myKeyFile, value);
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
     * Set Maximum Tag Length. Maximum Tag Length
     *
     * @param value Maximum Tag Length
     *
     * @return this plugin
     */
    public RepGenLoadSeqToDBPlugin kmerLength(Integer value) {
        myKmerLength = new PluginParameter<>(myKmerLength, value);
        return this;
    }

    /**
     * Minimum Tag Length
     *
     * @return Minimum Tag Length
     */
    public Integer minimumKmerLength() {
        return myMinKmerLength.value();
    }

    /**
     * Set Minimum Tag Length. Minimum Tag Length
     *
     * @param value Minimum Tag Length
     *
     * @return this plugin
     */
    public RepGenLoadSeqToDBPlugin minimumKmerLength(Integer value) {
        myMinKmerLength = new PluginParameter<>(myMinKmerLength, value);
        return this;
    }

    /**
     * Minimum tag count
     *
     * @return Min Tag Count
     */
    public Integer minKmerCount() {
        return myMinKmerCount.value();
    }

    /**
     * Set Min Tag Count. Minimum tag count
     *
     * @param value Min Tag Count
     *
     * @return this plugin
     */
    public RepGenLoadSeqToDBPlugin minKmerCount(Integer value) {
        myMinKmerCount = new PluginParameter<>(myMinKmerCount, value);
        return this;
    }

    /**
     * Minimum taxa containing tag
     *
     * @return Min Taxa Count
     */
    public Integer minTaxa() {
        return minTaxa.value();
    }

    /**
     * Set Min Taxa Count. Minimum Taxa count
     *
     * @param value Min Taxa containing a tag
     *
     * @return this plugin
     */
    public RepGenLoadSeqToDBPlugin minTaxa(Integer value) {
        minTaxa = new PluginParameter<>(minTaxa, value);
        return this;
    }
    /**
     * Output Database File
     *
     * @return Output Database File
     */
    public String outputDatabaseFile() {
        return myOutputDB.value();
    }

    /**
     * Set Output Database File. Output Database File
     *
     * @param value Output Database File
     *
     * @return this plugin
     */
    public RepGenLoadSeqToDBPlugin outputDatabaseFile(String value) {
        myOutputDB = new PluginParameter<>(myOutputDB, value);
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
    public RepGenLoadSeqToDBPlugin minimumQualityScore(Integer value) {
        myMinQualScore = new PluginParameter<>(myMinQualScore, value);
        return this;
    }

    /**
     * Set maximum number of kmers
     * @param value
     * @return
     */
    public RepGenLoadSeqToDBPlugin maximumKmerNumber(Integer value) {
        myMaxKmerNumber = new PluginParameter<>(myMaxKmerNumber, value);
        return this;
    }

    /**
     * Set number of Fastq files processed simultaneously
     * @param value
     * @return
     */
    public RepGenLoadSeqToDBPlugin batchSize(Integer value) {
        myBatchSize = new PluginParameter<>(myBatchSize, value);
        return this;
    }
//    /**
//     * Delete exisiting  DB
//     *
//     * @return deleteOldData
//     */
//    public Boolean deleteOldData() {
//        return myDeleteOldData.value();
//    }
//
//    /**
//     * Set Delete old data flag.  True indicates we want the
//     * db tables cleared
//     *
//     * @param value true/false - whether to delete data
//     *
//     * @return this plugin
//     */
//    public RepGenLoadSeqToDBPlugin deleteOldData(Boolean value) {
//        myDeleteOldData = new PluginParameter<>(myDeleteOldData, value);
//        return this;
//    }
    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Discovery Tags By Taxa";
    }

    @Override
    public String getToolTipText() {
        return "Discovery Tags By Taxa";
    }

    /**
     * This ConcurrentHashMap constrain the size of the map, and purges low distribution count tags when the size needs
     * to be reduced.
     */
    static class TagDistributionMap extends ConcurrentHashMap<Tag,TaxaDistribution> {
        private final int maxTagNum;
        private int minDepthToRetainInMap=2;
        private final int minCount;

        TagDistributionMap (int maxTagNumber, float loadFactor, int concurrencyLevel, int minCount) {
            super((maxTagNumber*2), loadFactor, concurrencyLevel);
            maxTagNum = maxTagNumber;
            this.minCount = minCount;
        }

        @Override
        public TaxaDistribution put(Tag key, TaxaDistribution value) {
            return super.put(key, value);
        }

        public synchronized void removeTagByCount(int minCnt) {
            entrySet().parallelStream()
                    .filter(e -> e.getValue().totalDepth()<minCnt)
                    .forEach(e -> remove(e.getKey()));
        }


        public long estimateMapMemorySize() {
            long size=0;
            int cnt=0;
            for (Entry<Tag, TaxaDistribution> entry : entrySet()) {
                size+=8+16+1; //Tag size
                size+=16; //Map references
                size+=entry.getValue().memorySize();
                cnt++;
                if(cnt>10000) break;
        }
            long estSize=(size()/cnt)*size;
            return estSize;
        }

        public long[] depthDistribution() {
            long[] base2bins=new long[34];
            int cnt=0;
            for (Entry<Tag, TaxaDistribution> entry : entrySet()) {
                base2bins[31-Integer.numberOfLeadingZeros(entry.getValue().totalDepth())]++;
                cnt++;
               // if(cnt>100000) break;
            }
            return base2bins;
        }
    }
    
    // LCJ - THe easiest and fastest method is to keep track of all the
    // quality strings seen for each iteration of the tag.  SOme of these
    // tags will be dropped when we removeTagWIthoutReplication() and other
    // methods.  Those need to be synced up with this list.
    
    //If we create a hashmap of tag/list-of-quality-strings, then after
    // culling the list at the end, we can create the average score.  The
    // number of strings on the list give us the number of tag instances.
    
    // THe map itself is a concurrent hashmap, but the list also must
    // be concurrently.  It can be defined as "List<Sting>" but the value
    // added must be a synchonizedList
    
    // the list that is added needs to be a synchonized list.  
    // http://stackoverflow.com/questions/5923405/concurrenthashmap-with-arraylist-as-value
    // Cannot make this a multimap as multimap will not store duplicate quality scores,
    public static class TagCountQualityScoreMap extends ConcurrentHashMap<Tag,List<String>> {
        private  int tagCount;

        TagCountQualityScoreMap (int maxTagNumber, float loadFactor, int concurrencyLevel, int minCount) {
            super((maxTagNumber*2), loadFactor, concurrencyLevel);

        }

        public synchronized void removeTagByCount(int minCnt) {
            // remove tags which are not represented a minimum number of times
            // THis is determimed by how many quality strings occur for the tag
            // should be consistent with removeTagByCount for TagDistributionMap
            entrySet().parallelStream()
                    .filter(e -> e.getValue().size() <minCnt)
                    .forEach(e -> remove(e.getKey()));
        }
    }
    
}


