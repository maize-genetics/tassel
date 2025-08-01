/*
 *  AdjustPhasingPlugin
 * 
 *  Created on Feb 17, 2017
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class AdjustPhasingPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(AdjustPhasingPlugin.class);

    private static final int NUM_VCF_HEADER_COLUMNS = 9;

    private PluginParameter<String> myInputVCFFile = new PluginParameter.Builder<>("inputVCFFile", null, String.class)
            .description("Input VCF file")
            .inFile()
            .guiName("Input VCF File")
            .required(true)
            .build();

    private PluginParameter<String> myHapcutDir = new PluginParameter.Builder<>("hapcutDir", null, String.class)
            .description("Directory containing Hapcut output files.")
            .inDir()
            .required(true)
            .build();

    private PluginParameter<String> myOutputVCFFile = new PluginParameter.Builder<>("outputVCFFile", null, String.class)
            .description("Output VCF file")
            .outFile()
            .guiName("Output VCF File")
            .required(true)
            .build();

    private int myNumTaxa = 0;
    private String[] myTaxaNames = null;
    private int myNumSites = 0;

    public AdjustPhasingPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        myNumSites = Utils.getNumberLinesNotHashOrBlank(inputVCFFile());

        try {

            ExecutorService pool = Executors.newWorkStealingPool();

            List<Future<ProcessHapcut>> futures = new ArrayList<>();

            File temp = new File(hapcutDir());
            for (File current : temp.listFiles()) {
                String name = current.getCanonicalPath();
                if (name.endsWith("_haplotypes")) {
                    ProcessHapcut process = new ProcessHapcut(name, myNumSites);
                    futures.add(pool.submit(process));
                }
            }

            String[] positionsFromHapcut = new String[myNumSites];
            Map<String, ProcessHapcut> processedHapcutFiles = new HashMap<>();
            for (Future<ProcessHapcut> future : futures) {
                ProcessHapcut processed = future.get();
                for (int s = 0; s < myNumSites; s++) {
                    if (positionsFromHapcut[s] == null) {
                        positionsFromHapcut[s] = processed.myPositions[s];
                    } else if (processed.myPositions[s] != null && !positionsFromHapcut[s].equals(processed.myPositions[s])) {
                        throw new IllegalStateException("AdjustPhasingPlugin: position at site: " + s + " in: " + processed.myFilename + ": " + processed.myPositions[s] + " doesn't match other Hapcut files: " + positionsFromHapcut[s]);
                    }
                }
                processedHapcutFiles.put(processed.taxaName(), processed);
                myLogger.info("finished processing: " + processed.filename() + " taxa: " + processed.taxaName());
            }

            try (BufferedReader reader = Utils.getBufferedReader(inputVCFFile());
                    BufferedWriter writer = Utils.getBufferedWriter(outputVCFFile())) {

                String line = reader.readLine();
                while (line != null && line.startsWith("##")) {
                    writer.write(line);
                    line = reader.readLine();
                }

                if (line == null || !line.startsWith("#CHROM")) {
                    throw new IllegalArgumentException("AdjustPhasingPlugin: processData: First line after ## lines should be header #CHROM...");
                }

                writer.write(line);
                writer.write('\n');

                // Taxa Names
                String[] tokens = line.split("\t");
                myNumTaxa = tokens.length - NUM_VCF_HEADER_COLUMNS;

                myLogger.info("Number of Taxa: " + myNumTaxa);
                myLogger.info("Number of Sites: " + myNumSites);

                myTaxaNames = new String[myNumTaxa];
                Switch[] switches = new Switch[myNumTaxa];
                int numIntersectingTaxa = 0;
                int numTaxaWithOutHapcut = 0;
                for (int i = NUM_VCF_HEADER_COLUMNS; i < tokens.length; i++) {
                    myTaxaNames[i - NUM_VCF_HEADER_COLUMNS] = tokens[i];
                    ProcessHapcut current = processedHapcutFiles.remove(tokens[i]);
                    if (current == null) {
                        switches[i - NUM_VCF_HEADER_COLUMNS] = FALSE_SWITCH;
                        numTaxaWithOutHapcut++;
                    } else {
                        switches[i - NUM_VCF_HEADER_COLUMNS] = current.result();
                        numIntersectingTaxa++;
                    }
                }
                int numHapcutNotInVCF = processedHapcutFiles.size();

                myLogger.info("Number intersecting taxa between Hapcut files and VCF: " + numIntersectingTaxa);
                myLogger.info("Number Taxa in VCF without Hapcut files: " + numTaxaWithOutHapcut);
                myLogger.info("Number Hapcut files without taxa in VCF: " + numHapcutNotInVCF);

                BlockingQueue<Future<ProcessLines>> queue = new LinkedBlockingQueue<>();

                Future<?> readFuture = pool.submit(new ReadLines(reader, queue, switches, pool));

                Future<?> writeFuture = pool.submit(new WriteLines(writer, queue));

                readFuture.get();

                writeFuture.get();

            } catch (Exception ex) {
                myLogger.debug(ex.getMessage(), ex);
                throw new IllegalStateException("AdjustPhasingPlugin: processData: problem converting file: " + inputVCFFile() + " to file: " + outputVCFFile() + " error: " + ex.getMessage());
            }

            pool.shutdown();

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
        }

        return null;

    }

    /**
     * Input VCF file
     *
     * @return Input VCF File
     */
    public String inputVCFFile() {
        return myInputVCFFile.value();
    }

    /**
     * Set Input VCF File. Input file
     *
     * @param value Input VCF File
     *
     * @return this plugin
     */
    public AdjustPhasingPlugin inputVCFFile(String value) {
        myInputVCFFile = new PluginParameter<>(myInputVCFFile, value);
        return this;
    }

    /**
     * Directory containing Hapcut output files.
     *
     * @return Hapcut Dir
     */
    public String hapcutDir() {
        return myHapcutDir.value();
    }

    /**
     * Set Hapcut Dir. Directory containing Hapcut output files.
     *
     * @param value Hapcut Dir
     *
     * @return this plugin
     */
    public AdjustPhasingPlugin hapcutDir(String value) {
        myHapcutDir = new PluginParameter<>(myHapcutDir, value);
        return this;
    }

    /**
     * Output VCF file
     *
     * @return Output VCF File
     */
    public String outputVCFFile() {
        return myOutputVCFFile.value();
    }

    /**
     * Set Output VCF File. Output VCF file
     *
     * @param value Output VCF File
     *
     * @return this plugin
     */
    public AdjustPhasingPlugin outputVCFFile(String value) {
        myOutputVCFFile = new PluginParameter<>(myOutputVCFFile, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Adjust Phasing";
    }

    @Override
    public String getToolTipText() {
        return "Adjust Phasing";
    }

    private class ProcessHapcut implements Callable<ProcessHapcut> {

        private final String myFilename;
        private final String myLogfile;
        private final String myTaxaName;
        private final BitSet myResult;
        private final String[] myPositions;
        private final List<String> myChromosomes = new ArrayList<>();
        private final List<Integer> myOffsets = new ArrayList<>();

        public ProcessHapcut(String filename, int numSites) {
            myFilename = filename;
            myLogfile = filename + ".log";
            myTaxaName = Utils.getFilename(filename).split("_")[0];
            myResult = new OpenBitSet(numSites);
            myPositions = new String[numSites];
        }

        @Override
        public ProcessHapcut call() throws Exception {

            try (BufferedReader reader = Utils.getBufferedReader(myFilename);
                    BufferedWriter log = Utils.getBufferedWriter(myLogfile)) {

                String line = reader.readLine();
                if (line == null) {
                    myLogger.warn("Hapcut file: " + myFilename + " is empty.");
                }

                int blockNum = 0;
                while (line != null) {

                    if (line.startsWith("BLOCK")) {
                        blockNum++;
                        processBlock(reader, line, blockNum, log);
                    } else {
                        throw new IllegalStateException("AdjustPhasingPlugin: expected BLOCK statement: " + line);
                    }

                    line = reader.readLine();

                }

            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("AdjustPhasingPlugin: ProcessHapcut: problem reading file: " + myFilename + ": " + e.getMessage());
            }

            return this;

        }

        private void processBlock(BufferedReader reader, String blockLine, int blockNum, BufferedWriter log) {

            String currentChr = null;
            int highestChrIndex = -1;
            List<Integer>[] phased = new List[2];
            phased[0] = new ArrayList<>();
            phased[1] = new ArrayList<>();
            Map<Integer, String> blockLines = new HashMap<>();

            try {

                String line = reader.readLine();
                while (line != null && !line.startsWith("********")) {

                    String[] tokens = line.split("\t");

                    if (currentChr == null) {
                        // tokens[3] is chromosome
                        currentChr = tokens[3];
                    } else if (!currentChr.equals(tokens[3])) {
                        throw new IllegalStateException("AdjustPhasingPlugin: ProcessHapcut: Different chr: " + tokens[3] + " within block that started with chr: " + currentChr);
                    }

                    // tokens[0] site number plus 1
                    int tempIndex = Integer.parseInt(tokens[0]) - 1;
                    blockLines.put(tempIndex, line);
                    if (highestChrIndex >= tempIndex) {
                        throw new IllegalStateException("AdjustPhasingPlugin: ProcessHapcut: index out of order: " + tempIndex);
                    } else {
                        highestChrIndex = tempIndex;
                    }

                    // tokens[4] is physical position
                    myPositions[tempIndex] = tokens[4];

                    // tokens[1] haploid A
                    // tokens[2] haploid B
                    // tokens[7] VCF genotype field
                    String[] fromVCF = tokens[7].split(":")[0].split("\\|");
                    if (tokens[1].equals(fromVCF[0]) || tokens[2].equals(fromVCF[1])) {
                        phased[0].add(tempIndex);
                    } else {
                        phased[1].add(tempIndex);
                    }

                    line = reader.readLine();

                }

            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("AdjustPhasingPlugin: ProcessHapcut: problem reading file: " + myFilename + ": " + e.getMessage());
            }

            int chrIndex = myChromosomes.indexOf(currentChr);
            if (chrIndex == -1) {
                myChromosomes.add(currentChr);
                myOffsets.add(highestChrIndex);
            } else if (myOffsets.get(chrIndex) < highestChrIndex) {
                myOffsets.add(chrIndex, highestChrIndex);
            }

            if (!phased[0].isEmpty() || !phased[1].isEmpty()) {
                List<Integer> smallest = null;
                if (phased[0].size() < phased[1].size()) {
                    smallest = phased[0];
                } else {
                    smallest = phased[1];
                }
                for (Integer index : smallest) {
                    myResult.fastSet(index);
                    try {
                        log.write(String.valueOf(blockNum));
                        log.write("\t");
                        log.write(blockLines.get(index));
                        log.write("\n");
                    } catch (Exception e) {
                        myLogger.debug(e.getMessage(), e);
                    }
                }
            }

        }

        public String filename() {
            return myFilename;
        }

        public String taxaName() {
            return myTaxaName;
        }

        public Switch result() {
            if (myResult.cardinality() == 0) {
                return FALSE_SWITCH;
            } else {
                return new SwitchBitSet(myResult);
            }
        }

        public String[] positions() {
            return myPositions;
        }

        public List<String> chromosomes() {
            return myChromosomes;
        }

        public List<Integer> offsets() {
            return myOffsets;
        }

    }

    private final Switch FALSE_SWITCH = new Switch();

    private class Switch {

        public boolean alleles(int site) {
            return false;
        }
    }

    private class SwitchBitSet extends Switch {

        private final BitSet myBitSet;

        public SwitchBitSet(BitSet bitSet) {
            myBitSet = bitSet;
        }

        public boolean alleles(int site) {
            return myBitSet.fastGet(site);
        }
    }

    private class ProcessLines implements Callable<ProcessLines> {

        private final int myStartSite;
        private final String[] myLines;
        private final Switch[] mySwitches;
        private final StringBuilder myBuilder = new StringBuilder();
        private final boolean myIsFinalProcessLines;
        private final int[] myNumSwitchesTaxa;
        private final int[] myNumSwitchesSites;

        public ProcessLines(int startSite, String[] lines, Switch[] switches) {
            myStartSite = startSite;
            myLines = lines;
            mySwitches = switches;
            myIsFinalProcessLines = false;
            myNumSwitchesTaxa = new int[myNumTaxa];
            myNumSwitchesSites = new int[myLines.length];
        }

        public ProcessLines() {
            myStartSite = 0;
            myLines = new String[0];
            mySwitches = null;
            myIsFinalProcessLines = true;
            myNumSwitchesTaxa = new int[0];
            myNumSwitchesSites = new int[0];
        }

        @Override
        public ProcessLines call() throws Exception {

            int currentSite = myStartSite;
            for (String current : myLines) {
                char[] temp = current.toCharArray();
                int numChars = temp.length;
                int copyToNthTab = NUM_VCF_HEADER_COLUMNS;
                int numTabsFound = 0;
                int startIndex = 0;
                int charIndex = 0;
                for (int t = 0; t < myNumTaxa; t++) {
                    if (mySwitches[t].alleles(currentSite)) {
                        myNumSwitchesTaxa[t]++;
                        myNumSwitchesSites[currentSite - myStartSite]++;
                        while (true) {
                            if (temp[charIndex] == '\t') {
                                numTabsFound++;
                                if (numTabsFound == copyToNthTab) {
                                    myBuilder.append(temp, startIndex, charIndex - startIndex + 1);
                                    charIndex++;
                                    char first = temp[charIndex];
                                    charIndex++;
                                    if (temp[charIndex] != '|') {
                                        throw new IllegalStateException("AdjustPhasingPlugin: ProcessLines: | char expected.");
                                    }
                                    charIndex++;
                                    myBuilder.append(temp[charIndex]);
                                    myBuilder.append('|');
                                    myBuilder.append(first);
                                    charIndex++;
                                    while (charIndex < numChars && temp[charIndex] != '\t') {
                                        myBuilder.append(temp[charIndex]);
                                        charIndex++;
                                    }
                                    if (charIndex < numChars) {
                                        myBuilder.append('\t');
                                    } else {
                                        myBuilder.append('\n');
                                    }
                                    startIndex = charIndex + 1;
                                    break;
                                }
                            }
                            charIndex++;
                        }
                        copyToNthTab++;
                    } else {
                        copyToNthTab++;
                    }
                }

                if (startIndex < numChars) {
                    myBuilder.append(temp, startIndex, numChars - startIndex);
                    myBuilder.append('\n');
                }

                currentSite++;
            }

            return this;

        }

    }

    private static final int NUM_LINES_PER_BLOCK = 100;

    private class ReadLines implements Runnable {

        private final BufferedReader myReader;
        private final BlockingQueue<Future<ProcessLines>> myQueue;
        private final Switch[] mySwitches;
        private final ExecutorService myPool;

        public ReadLines(BufferedReader reader, BlockingQueue<Future<ProcessLines>> queue, Switch[] switches, ExecutorService pool) {
            myReader = reader;
            myQueue = queue;
            mySwitches = switches;
            myPool = pool;
        }

        @Override
        public void run() {

            try {

                List<String> temp = new ArrayList<>();
                String line = myReader.readLine();
                int count = 0;
                int startSite = 0;
                while (line != null) {
                    temp.add(line);
                    count++;
                    if (count == NUM_LINES_PER_BLOCK) {
                        String[] lines = new String[NUM_LINES_PER_BLOCK];
                        myQueue.add(myPool.submit(new ProcessLines(startSite, temp.toArray(lines), mySwitches)));
                        temp.clear();
                        count = 0;
                        startSite += NUM_LINES_PER_BLOCK;
                    }
                    line = myReader.readLine();
                }

                if (!temp.isEmpty()) {
                    String[] lines = new String[temp.size()];
                    myQueue.add(myPool.submit(new ProcessLines(startSite, temp.toArray(lines), mySwitches)));
                }

                myQueue.add(myPool.submit(new ProcessLines()));

            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("AdjustPhasingPlugin: ReadLines: problem reading file: " + inputVCFFile());
            }
        }

    }

    private class WriteLines implements Runnable {

        private final BufferedWriter myWriter;
        private final BlockingQueue<Future<ProcessLines>> myQueue;

        public WriteLines(BufferedWriter writer, BlockingQueue<Future<ProcessLines>> queue) {
            myWriter = writer;
            myQueue = queue;
        }

        @Override
        public void run() {

            int[] numSwitchesTaxa = new int[myNumTaxa];
            int[] numSwitchesSites = new int[myNumSites];

            try {
                Future<ProcessLines> future = myQueue.take();
                ProcessLines processed = future.get();
                while (!processed.myIsFinalProcessLines) {
                    myWriter.write(processed.myBuilder.toString());

                    int percent = Math.min(Math.abs((int) ((double) processed.myStartSite / (double) myNumSites * 100.0)), 100);
                    progress(percent, this);

                    for (int t = 0; t < myNumTaxa; t++) {
                        numSwitchesTaxa[t] += processed.myNumSwitchesTaxa[t];
                    }
                    System.arraycopy(processed.myNumSwitchesSites, 0, numSwitchesSites, processed.myStartSite, processed.myLines.length);
                    future = myQueue.take();
                    processed = future.get();
                }
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("AdjustPhasingPlugin: WriteLines: problem writing file: " + outputVCFFile());
            }

            for (int t = 0; t < myNumTaxa; t++) {
                System.out.println(t + ": " + myTaxaNames[t] + ": alleles switched: " + numSwitchesTaxa[t]);
            }

        }

    }

}
