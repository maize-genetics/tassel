package net.maizegenetics.analysis.rna;

import java.awt.Frame;
import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.LongAdder;

import javax.swing.ImageIcon;

import net.maizegenetics.analysis.gbs.Barcode;
import net.maizegenetics.analysis.gbs.v2.BarcodeTrie;
import net.maizegenetics.analysis.gbs.v2.EnzymeList;
import net.maizegenetics.analysis.gbs.v2.GBSUtils;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TagDataWriter;
import net.maizegenetics.dna.tag.TaxaDistBuilder;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class RNADeMultiPlexSeqToDBPlugin extends AbstractPlugin{
    private static final Logger myLogger = LogManager.getLogger(RNADeMultiPlexSeqToDBPlugin.class);
    static LongAdder roughTagCnt = new LongAdder();
    private PluginParameter<String> myInputDir = new PluginParameter.Builder<>("i", null, String.class).guiName("Input Directory").required(true).inDir()
            .description("Input directory containing FASTQ files in text or gzipped text.\n"
                    + "     NOTE: Directory will be searched recursively and should\n"
                    + "     be written WITHOUT a slash after its name.").build();
    private PluginParameter<String> myKeyFile = new PluginParameter.Builder<>("k", null, String.class).guiName("Key File").required(true).inFile()
            .description("Key file listing barcodes distinguishing the samples").build();
    private PluginParameter<Integer> myMinKmerLength = new PluginParameter.Builder<>("minKmerL", 20, Integer.class).guiName("Minimum Kmer Length")
            .description("Minimum kmer Length after second cut site is removed").build();
    private PluginParameter<Integer> myMinKmerCount = new PluginParameter.Builder<>("c", 10, Integer.class).guiName("Min Kmer Count")
            .description("Minimum kmer count").build();
    private PluginParameter<String> myOutputDB = new PluginParameter.Builder<>("db", null, String.class).guiName("Output Database File").required(true).outFile()
            .description("Output Database File").build();
    private PluginParameter<Integer> myMinQualScore = new PluginParameter.Builder<>("mnQS", 0, Integer.class).guiName("Minimum quality score").required(false)
            .description("Minimum quality score within the barcode and read length to be accepted").build();
    
  
    private static String myEnzyme = "ignore";
    private static Integer myMaxKmerNumber = 50000000;
    private static Integer myBatchSize = 8;
    static final String inputFileGlob="glob:*{.fq,fq.gz,fastq,fastq.txt,fastq.gz,fastq.txt.gz,_sequence.txt,_sequence.txt.gz}";
    static final String sampleNameField="FullSampleName";
    static final String flowcellField="Flowcell";
    static final String laneField="Lane";
    static final String barcodeField="Barcode";

    private static TagDistributionMap tagCntMap;
    private static boolean taglenException;

    @Override
    public DataSet processData(DataSet input) {
        float loadFactor = 0.95f;
        tagCntMap = new TagDistributionMap (myMaxKmerNumber,loadFactor, 128, minKmerCount());  
        try {
            //Get the list of fastq files
            Path keyPath= Paths.get(keyFile()).toAbsolutePath();
            List<Path> directoryFiles= DirectoryCrawler.listPaths(inputFileGlob, Paths.get(inputDirectory()).toAbsolutePath());
            if(directoryFiles.isEmpty()) { 
                myLogger.warn("No files matching:"+ inputFileGlob);
                System.out.println("RNADeMultiPlex - no files matching " + inputFileGlob);
                return null;
            } 

            // Cull files that are not represented in the given key file 
            // This is useful when a user keeps all fastq files in a single
            // directory, but uses the key file to determine which files to run.
            // Commented out here as culledFiles expects to find the flowcell and
            // lane from the file name, required 3,4 or 5 underscores with the
            // the number of underscores determining the flowcell/lane position/value.
            // NOTE: When we "cull" the files, we don't read them, we merely look at the names
            // to determine whether to include them or not.
           // List<Path> inputSeqFiles = GBSUtils.culledFiles(directoryFiles,keyPath);
            List<Path> inputSeqFiles = directoryFiles;
            if (inputSeqFiles.size() == 0) return null; // no files in this directory to process
            int batchNum = inputSeqFiles.size()/myBatchSize;
            if (inputSeqFiles.size()%myBatchSize != 0) batchNum++;
            TaxaList masterTaxaList= TaxaListIOUtils.readTaxaAnnotationFile(keyFile(), sampleNameField, new HashMap<>(), true);
            
            //  clear existing db.
            if (Files.exists(Paths.get(myOutputDB.value()))) {
                try {
                    Files.delete(Paths.get(outputDatabaseFile()));
                 } catch (Exception exc){
                    System.out.println("Error when trying to delete database file: " + outputDatabaseFile());
                    System.out.println("File delete error: " + exc.getMessage());
                    return null;
                 }
            }
  
            TagDataWriter tdw=new TagDataSQLite(outputDatabaseFile());
            taglenException = false;
            for (int i = 0; i < inputSeqFiles.size(); i+=myBatchSize) {
                int end = i+myBatchSize;
                if (end > inputSeqFiles.size()) end = inputSeqFiles.size();
                ArrayList<Path> sub = new ArrayList<Path>();
                for (int j = i; j < end; j++) sub.add(inputSeqFiles.get(j));
                System.out.println("\nStart processing batch " + String.valueOf(i/myBatchSize+1));
                sub.parallelStream()
                .forEach(inputSeqFile -> {
                    try {
                        processFastQFile(masterTaxaList,keyPath, inputSeqFile, myEnzyme,
                                minimumQualityScore(), minimumKmerLength(), tagCntMap);
                    } catch (StringIndexOutOfBoundsException oobe) {
                        oobe.printStackTrace();
                        myLogger.error(oobe.getMessage());
                        setTagLenException();
                        return;
                    }
                });
                if (taglenException == true) return null; // Tag length failure from processFastQ - halt processing
                System.out.println("\nKmers are added from batch "+String.valueOf(i/myBatchSize+1) + ". Total batch number: " + batchNum);
                int currentSize = tagCntMap.size();
                System.out.println("Current number: " + String.valueOf(currentSize) + ". Max kmer number: " + String.valueOf(myMaxKmerNumber));
                System.out.println(String.valueOf((float)currentSize/(float)myMaxKmerNumber) + " of max tag number");

                if (currentSize > 0) { // calcTagMapStats() gets "divide by 0" error when size == 0
                    calcTagMapStats(tagCntMap);
                    System.out.println();
                    roughTagCnt.reset();
                    roughTagCnt.add(tagCntMap.size());
                } else {
                    System.out.println("WARNING: Current tagcntmap size is 0 after processing batch " + String.valueOf(i/myBatchSize+1) );
                }

                System.out.println("Total memory: "+ String.valueOf((double)(Runtime.getRuntime().totalMemory()/1024/1024/1024))+" Gb");
                System.out.println("Free memory: "+ String.valueOf((double)(Runtime.getRuntime().freeMemory()/1024/1024/1024))+" Gb");
                System.out.println("Max memory: "+ String.valueOf((double)(Runtime.getRuntime().maxMemory()/1024/1024/1024))+" Gb");
                System.out.println("\n");
            }
            System.out.println("\nAll the batch are processed");
            tagCntMap.removeTagByCount(minKmerCount());
            System.out.println("By removing kmers with minCount of " + myMinKmerCount + " Kmer number is reduced to " + tagCntMap.size()+"\n");
            
            tdw.putTaxaList(masterTaxaList);
            tdw.putAllTag(tagCntMap.keySet());
         // For RNA, taxaTissueDist will be stored when Production is run
            //tdw.putTaxaDistribution(tagCntMap);
            ((TagDataSQLite)tdw).close();
        } catch(Exception e) {
            e.printStackTrace();
        }
        return null;
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
    private static void processFastQFile(TaxaList masterTaxaList, Path keyPath, Path fastQPath, String enzymeName,
            int minQuality, int minKmerLen, TagDistributionMap masterTagTaxaMap) throws StringIndexOutOfBoundsException {
       // ArrayList<Taxon> tl=GBSUtils.getLaneAnnotatedTaxaList(keyPath, fastQPath);
        // Call readTaxaANnotationFileAL directly to avoid needing to name files
        // with a specific number of underscores.
        ArrayList<Taxon> tl=TaxaListIOUtils.readTaxaAnnotationFileAL(keyPath.toAbsolutePath().toString(), GBSUtils.sampleNameField, new HashMap<>());
        if (tl.size() == 0) return; 
        BarcodeTrie barcodeTrie=GBSUtils.initializeBarcodeTrie(tl, masterTaxaList, EnzymeList.defaultCache.getEnzyme(enzymeName));
        try {
                processFastQ(fastQPath,barcodeTrie,masterTaxaList,masterTagTaxaMap,
                        minQuality, minKmerLen);
        } catch (StringIndexOutOfBoundsException oobe) {
                throw oobe; // Let processData() handle it - we want to stop processing on this error
        }        
    }

    private static void processFastQ(Path fastqFile, BarcodeTrie barcodeTrie, TaxaList masterTaxaList,
            TagDistributionMap masterTagTaxaMap, int minQual, int minKmerLen) throws StringIndexOutOfBoundsException{
        int allReads=0, goodBarcodedReads = 0, lowQualityReads = 0, badNoBarcode = 0, nullTags = 0;
        int shortReads = 0;
        int maxTaxaNumber=masterTaxaList.size();
        int checkSize = 10000000;
        myLogger.info("processing file " + fastqFile.toString());
        try {
            int qualityScoreBase=GBSUtils.determineQualityScoreBase(fastqFile);
            BufferedReader br = Utils.getBufferedReader(fastqFile.toString(), 1 << 22);
            long time=System.nanoTime();
            String[] seqAndQual;
            String likelyReadEnd = "AGATCGGA";
 
            while ((seqAndQual=GBSUtils.readDeMultiPlexFastQBlock(br,allReads)) != null) {
                allReads++;
                //After quality score is read, decode barcode using the current sequence & quality  score
                int adapterStart = 0;
                String tmpSeq = seqAndQual[2] + seqAndQual[0];
                Barcode barcode =barcodeTrie.longestPrefix(tmpSeq); 
                if(barcode==null) {
                    System.out.println("BC not found: " + seqAndQual[0]);
                    badNoBarcode++;
                    continue;
                }

                String sequence = seqAndQual[0];
                // Substring to remove the common adaptor
                // and everything beyond it               
                adapterStart = sequence.indexOf(likelyReadEnd);
                if (adapterStart > 0) {
                    sequence = sequence.substring(0, adapterStart-1);
                }
 
                if (sequence.length() < minKmerLen){
                     shortReads++;
                    continue;
                }

                if(minQual>0) {
                    //todo move getFirstLowQualityPos into this class?
                    if(BaseEncoder.getFirstLowQualityPos(seqAndQual[1],minQual, qualityScoreBase) < sequence.length()){
                        lowQualityReads++;
                        continue;
                    }
                } 
 
                Tag tag = null;
                tag=TagBuilder.instance(sequence).build();
                               
                if(tag==null) {
                    nullTags++;
                    continue;   //null occurs when any base was not A, C, G, T
                }
                goodBarcodedReads++;
                TaxaDistribution taxaDistribution=masterTagTaxaMap.get(tag);
                if(taxaDistribution==null) {
                    masterTagTaxaMap.put(tag,TaxaDistBuilder.create(maxTaxaNumber,barcode.getTaxaIndex()));
                    roughTagCnt.increment();
                } else {
                    taxaDistribution.increment(barcode.getTaxaIndex());
                }
                if (allReads % checkSize == 0) {
                    myLogger.info("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads
                            + " rate:" + (System.nanoTime()-time)/allReads +" ns/read. Current tag count:" + roughTagCnt);
                }
                if (allReads % checkSize == 0) {
                    myLogger.info("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads
                            + " rate:" + (System.nanoTime()-time)/allReads +" ns/read. Current tag count:" + roughTagCnt);
                }
            }
            myLogger.info("Summary for "+fastqFile.toString()+"\n"+
                    "Total number of reads in lane=" + allReads +"\n"+
                    "Total number of good barcoded reads=" + goodBarcodedReads+"\n"+
                    "Total number of low quality reads=" + lowQualityReads+"\n"+
                    "Total number of short reads=" + shortReads+"\n"+
                    "Total number of bad or no barcode found=" + badNoBarcode+"\n" +
                    "Total number of null tags created=" + nullTags+"\n" +
                    "Timing process (sorting, collapsing, and writing TagCount to file)."+"\n"+
                    "Process took " + (System.nanoTime() - time)/1e6 + " milliseconds.");
            br.close();
        } catch (StringIndexOutOfBoundsException oobe) {
                throw oobe; // pass it up to print error and stop processing
        } catch (Exception e) {
            myLogger.error("Good Barcodes Read: " + goodBarcodedReads);            
            e.printStackTrace();
        }
    }   

    public static void setTagLenException() {
        taglenException = true;
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
            for (Map.Entry<Tag, TaxaDistribution> entry : entrySet()) {
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
            for (Map.Entry<Tag, TaxaDistribution> entry : entrySet()) {
                base2bins[31-Integer.numberOfLeadingZeros(entry.getValue().totalDepth())]++;
                cnt++;
               // if(cnt>100000) break;
            }
            return base2bins;
        }
    }


    public RNADeMultiPlexSeqToDBPlugin() {
        super(null, false);
    }

    public RNADeMultiPlexSeqToDBPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public RNADeMultiPlexSeqToDBPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
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
    public RNADeMultiPlexSeqToDBPlugin inputDirectory(String value) {
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
    public RNADeMultiPlexSeqToDBPlugin keyFile(String value) {
        myKeyFile = new PluginParameter<>(myKeyFile, value);
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
    public RNADeMultiPlexSeqToDBPlugin minimumKmerLength(Integer value) {
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
    public RNADeMultiPlexSeqToDBPlugin minKmerCount(Integer value) {
        myMinKmerCount = new PluginParameter<>(myMinKmerCount, value);
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
    public RNADeMultiPlexSeqToDBPlugin outputDatabaseFile(String value) {
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
    public RNADeMultiPlexSeqToDBPlugin minimumQualityScore(Integer value) {
        myMinQualScore = new PluginParameter<>(myMinQualScore, value);
        return this;
    }
}
