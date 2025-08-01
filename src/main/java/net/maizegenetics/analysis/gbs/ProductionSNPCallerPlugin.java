/*
 * ProductionSNPCallerPlugin
 */
package net.maizegenetics.analysis.gbs;

import cern.colt.list.IntArrayList;
import com.google.common.collect.ImmutableMap;

import com.google.common.collect.ImmutableTable;
import com.google.common.collect.Multimap;
import com.google.common.collect.Table;
import com.google.common.collect.TreeMultimap;

import net.maizegenetics.dna.map.*;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.MultiMemberGZIPInputStream;
import net.maizegenetics.dna.snp.score.AlleleDepthUtil;
import net.maizegenetics.dna.snp.genotypecall.BasicGenotypeMergeRule;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.Taxon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.*;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.util.GeneralAnnotation;

/**
 * This plugin converts all of the fastq (and/or qseq) files in the input folder
 * and keyfile to genotypes and adds these to a genotype file in HDF5 format.
 *
 * We refer to this step as the "Production Pipeline".
 *
 * The output format is HDF5 genotypes with allelic depth stored. SNP calling is
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
 * If the target output HDF5 GenotypeTable file doesn't exist, it will be
 * created.
 *
 * Each taxon in the HDF5 file is named "ShortName:LibraryPrepID" and is
 * annotated with "Flowcell_Lanes" (=source seq data for current genotype).
 *
 * Requires a TOPM with variants added from a previous "Discovery Pipeline" run.
 * In binary topm or HDF5 format (TOPMInterface).
 *
 * TODO add the Stacks likelihood method to BasicGenotypeMergeRule
 *
 * @author Jeff Glaubitz
 * @author Ed Buckler
 */
public class ProductionSNPCallerPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(ProductionSNPCallerPlugin.class);

    private PluginParameter<String> myInputDirectory = new PluginParameter.Builder<>("i", null, String.class).guiName("Input Directory").required(true).inDir()
            .description("Input directory containing fastq AND/OR qseq files.").build();
    private PluginParameter<String> myKeyFile = new PluginParameter.Builder<>("k", null, String.class).guiName("Key File").required(true).inFile()
            .description("Key file listing barcodes distinguishing the samples").build();
    private PluginParameter<String> myEnzyme = new PluginParameter.Builder<>("e", null, String.class).guiName("Enzyme").required(true)
            .description("Enzyme used to create the GBS library").build();
    private PluginParameter<String> myProductionTOPM = new PluginParameter.Builder<>("m", null, String.class).guiName("Input TOPM File").required(true).inFile()
            .description("Physical map file containing tags and corresponding variants (production TOPM)").build();
    private PluginParameter<String> myOutputGenotypes = new PluginParameter.Builder<>("o", null, String.class).guiName("Output HDF5 Genotypes File").required(true).outFile()
            .description("Output (target) HDF5 genotypes file to add new genotypes to (new file created if it doesn't exist)").build();
    private PluginParameter<Double> myAveSeqErrorRate = new PluginParameter.Builder<>("eR", 0.01, Double.class).guiName("Ave Seq Error Rate")
            .description("Average sequencing error rate per base (used to decide between heterozygous and homozygous calls)").build();
//    private PluginParameter<Integer> myMaxDivergence = new PluginParameter.Builder<>("d", 0, Integer.class).guiName("Max Divergence")
//            .description("Maximum divergence (edit distance) between new read and previously mapped read (Default: 0 = perfect matches only)").build();
    private PluginParameter<Boolean> myKeepGenotypesOpen = new PluginParameter.Builder<>("ko", false, Boolean.class).guiName("Keep Genotypes Open")
            .description("Keep hdf5 genotypes open for future runs that add more taxa or more depth").build();
    private PluginParameter<Boolean> myNoDepthOutput = new PluginParameter.Builder<>("ndo", false, Boolean.class).guiName("No Depth to Output")
            .description("No depth output: do not write depths to the output hdf5 genotypes file").build();
    //private PluginParameter<Boolean> myStacksLikelihood = new PluginParameter.Builder<>("sL", false, Boolean.class).guiName("Use Stacks Likelihood")
    //        .description("Use STACKS likelihood method to call heterozygotes (default: use tasselGBS likelihood ratio method)").build();

    private String[] myRawSeqFileNames = null;
    private String myOutputDir = null;
    private TOPMInterface topm = null;
    private int[] chromosomes = null;
    private boolean fastq = true;
    private Map<String, Integer> keyFileColumns = new HashMap<>();
    private Set<String> flowcellLanesInKey = new TreeSet<>();  // flowcellLanes present in key file (stored here as "Flowcell_Lane")
    private Map<String, String> seqFileNameToFlowcellLane = new HashMap<>(); // map from name of fastq or qseq file to "Flowcell_Lane"
    private Set<String> seqFilesInKeyAndDir = new TreeSet<>(); // fastq (or qseq) file names present in input directory that have a "Flowcell_Lane" in the key file
    private Map<String, String> fullNameToHDF5Name = new TreeMap<>();
    private Multimap<String, String> flowCellLaneToLibPrepIDs = TreeMultimap.create();
    private Map<String, String> libraryPrepIDToSampleName = new TreeMap<String, String>();
    
    private GenotypeTableBuilder genos = null; //output genotype table
    private TaxaList taxaList = null;
    private PositionList myPositionList = null;
    private IntArrayList[] obsTagsForEachTaxon = null;
    private Table<Integer, Integer, Integer> positionToSite = null;  // indices = chrIndices.  For a given position (key), each HashMap provides the site in the MutableNucleotideDepthAlignment (value)

    //Documentation of read depth per sample (one recored per replicate)
    private Map<String, Integer> rawReadCountsForFullSampleName = new TreeMap<>();
    private Map<String, Integer> matchedReadCountsForFullSampleName = new TreeMap<>();

    private BasicGenotypeMergeRule genoMergeRule = null;

    public ProductionSNPCallerPlugin() {
        super(null, false);
    }

    public ProductionSNPCallerPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public void postProcessParameters() {

        File rawSeqDirectory = new File(inputDirectory());
        myRawSeqFileNames = DirectoryCrawler.listFileNames(rawSeqFileNameRegex, rawSeqDirectory.getAbsolutePath());
        if (myRawSeqFileNames.length == 0 || myRawSeqFileNames == null) {
            throw new IllegalArgumentException(noMatchingRawSeqFileNamesMessage + inputDirectory());
        } else {
            Arrays.sort(myRawSeqFileNames);
            String message = "\nProductionSNPCallerPlugin:\n\nThe following GBS raw sequence data files were found in the input folder (and sub-folders):";
            for (String filename : myRawSeqFileNames) {
                message += "\n   " + filename;
            }
            myLogger.info(message + "\n");
        }

        topm = TOPMUtils.readTOPM(inputTOPMFile());
        if (topm.getSize() == 0) {
            throw new IllegalStateException("TagsOnPhysicalMap file not available or is empty");
        }

        try {
            myOutputDir = (new File(outputHDF5GenotypesFile())).getCanonicalFile().getParent();
        } catch (IOException e) {
            throw new IllegalStateException("Problem resolving output directory:" + e);
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        long previous, current;
        readKeyFile();  // TODO: read/write full set of metadata
        matchKeyFileToAvailableRawSeqFiles();
        myPositionList = getUniquePositions();
        generateFastSiteLookup(myPositionList);
        setUpGenotypeTableBuilder();
        int nFilesProcessed = 0;
        for (int fileNum = 0; fileNum < myRawSeqFileNames.length; fileNum++) {
            if (!seqFilesInKeyAndDir.contains(myRawSeqFileNames[fileNum])) {
                continue;  // skip fastq/qseq files that are not in the key file
            }
            int[] counters = {0, 0, 0, 0, 0, 0}; // 0:allReads 1:goodBarcodedReads 2:goodMatched 3:perfectMatches 4:imperfectMatches 5:singleImperfectMatches
            myLogger.info("\nLooking for known SNPs in sequence reads from file:");
            myLogger.info("  " + myRawSeqFileNames[fileNum] + "\n");
            previous = System.nanoTime();
            readRawSequencesAndRecordDepth(fileNum, counters);  // TODO: read the machine name from the fastq/qseq file
            current = System.nanoTime();
            System.out.println("ProductionSNPCallerPlugin: performFunction: readRawSequencesAndRecordDepth: " + myRawSeqFileNames[fileNum] + ": " + ((double) (current - previous) / 1_000_000_000.0) + " sec");
            previous = System.nanoTime();
            callGenotypes();
            current = System.nanoTime();
            System.out.println("ProductionSNPCallerPlugin: performFunction: callGenotypes: " + myRawSeqFileNames[fileNum] + ": " + ((double) (current - previous) / 1_000_000_000.0) + " sec");
            ++nFilesProcessed;
            reportTotals(fileNum, counters, nFilesProcessed);
        }
        if (keepGenotypesOpen()) {
            previous = System.nanoTime();
            genos.closeUnfinished();
            current = System.nanoTime();
            System.out.println("ProductionSNPCallerPlugin: performFunction: genos.closeUnfinished(): " + ((double) (current - previous) / 1_000_000_000.0) + " sec");
        } else {
            previous = System.nanoTime();
            genos.build();
            current = System.nanoTime();
            System.out.println("ProductionSNPCallerPlugin: performFunction: genos.build(): " + ((double) (current - previous) / 1_000_000_000.0) + " sec");
        }
        writeReadsPerSampleReports();
        return null;
    }

    private void readRawSequencesAndRecordDepth(int fileNum, int[] counters) {
        long previous, current;
        ParseBarcodeRead thePBR = setUpBarcodes(fileNum);
        if (thePBR == null || thePBR.getBarCodeCount() == 0) {
            myLogger.info("No barcodes found. Skipping this raw sequence file.");
            return;
        }
        taxaList = new TaxaListBuilder().addAll(getHDF5Taxa(fileNum)).sortTaxaAlphabetically().build();
        obsTagsForEachTaxon = new IntArrayList[taxaList.numberOfTaxa()];
        for (int t = 0; t < obsTagsForEachTaxon.length; t++) {
            obsTagsForEachTaxon[t] = new IntArrayList(750_000); // initial capacity
        }
        String temp = "Nothing has been read from the raw sequence file yet";
        BufferedReader br = getBufferedReaderForRawSeqFile(fileNum);
        try {
            long readSeqReadTime = 0;
            long ifRRNotNullTime = 0;
            while ((temp = br.readLine()) != null) {
                if (counters[0] % 1000000 == 0) {
                    reportProgress(counters, readSeqReadTime, ifRRNotNullTime);
                }
                previous = System.nanoTime();
                ReadBarcodeResult rr = readSequenceRead(br, temp, thePBR, counters);
                current = System.nanoTime();
                readSeqReadTime += (current - previous);
                previous = System.nanoTime();
                if (rr != null) {
                    counters[1]++;  // goodBarcodedReads
                    rawReadCountsForFullSampleName.put(rr.getTaxonName(), rawReadCountsForFullSampleName.get(rr.getTaxonName()) + 1);
                    int tagIndex = topm.getTagIndex(rr.getRead());
                    if (tagIndex >= 0) {
                        counters[3]++;  // perfectMatches
                    }
//                    if (tagIndex < 0 && maxDivergence() > 0) {
//                        tagIndex = findBestImperfectMatch(rr.getRead(), counters);
//                    }
                    if (tagIndex < 0) {
                        continue;
                    }
                    counters[2]++;  // goodMatched++;
                    matchedReadCountsForFullSampleName.put(rr.getTaxonName(), matchedReadCountsForFullSampleName.get(rr.getTaxonName()) + 1);
                    int taxonIndex = taxaList.indexOf(fullNameToHDF5Name.get(rr.getTaxonName()));
                    if (taxonIndex == -1) {
                        System.out.println("We're here!");
                    }
                    obsTagsForEachTaxon[taxonIndex].add(tagIndex);
                }
                current = System.nanoTime();
                ifRRNotNullTime += (current - previous);
            }
            System.out.println("ProductionSNPCallerPlugin: readRawSequencesAndRecordDepth: readSequenceTime: " + ((double) (readSeqReadTime) / 1_000_000_000.0) + " sec");
            System.out.println("ProductionSNPCallerPlugin: readRawSequencesAndRecordDepth: processSequenceTime: " + ((double) (ifRRNotNullTime) / 1_000_000_000.0) + " sec");
            br.close();
        } catch (Exception e) {
            myLogger.error("Catch in readRawSequencesAndRecordDepth() at nReads=" + counters[0] + " e=" + e);
            myLogger.error("Last line read: " + temp);
            e.printStackTrace();
            System.exit(1);
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

    private void reportTotals(int fileNum, int[] counters, int nFilesProcessed) {
        myLogger.info("Total number of reads in lane=" + counters[0]);
        myLogger.info("Total number of good, barcoded reads=" + counters[1]);
        myLogger.info("Total number of good, barcoded reads matched to the TOPM=" + counters[2]);
        myLogger.info("Finished reading " + nFilesProcessed + " of " + seqFilesInKeyAndDir.size() + " sequence files: " + myRawSeqFileNames[fileNum] + "\n");
    }

    private void readKeyFile() {
        flowcellLanesInKey.clear();
        seqFileNameToFlowcellLane.clear();
        seqFilesInKeyAndDir.clear();
        fullNameToHDF5Name.clear();
        flowCellLaneToLibPrepIDs.clear();
        libraryPrepIDToSampleName.clear();
        String inputLine = "Nothing has been read from the keyfile yet";
        try {
            BufferedReader br = new BufferedReader(new FileReader(keyFile()), 65536);
            int currLine = 0;
            while ((inputLine = br.readLine()) != null) {
                if (currLine == 0) {
                    parseKeyFileHeader(inputLine);
                } else {
                    populateKeyFileFields(inputLine);
                }
                currLine++;
            }
        } catch (Exception e) {
            myLogger.error("Couldn't read key file: " + e);
            myLogger.error("Last line read from key file: " + inputLine);
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void parseKeyFileHeader(String headerLine) {
        headerLine.trim();
        String[] header = headerLine.split("\\t");
        keyFileColumns.clear();
        for (int col = 0; col < header.length; col++) {
            if (header[col].equalsIgnoreCase("Flowcell")) {
                keyFileColumns.put("Flowcell", col);
            } else if (header[col].equalsIgnoreCase("Lane")) {
                keyFileColumns.put("Lane", col);
            } else if (header[col].equalsIgnoreCase("Barcode")) {
                keyFileColumns.put("Barcode", col);
            } else if (header[col].equalsIgnoreCase("DNASample") || header[col].equalsIgnoreCase("Sample")) {
                keyFileColumns.put("Sample", col);
            } else if (header[col].equalsIgnoreCase("LibraryPrepID")) {
                keyFileColumns.put("LibPrepID", col);
            }
        }
        if (!confirmKeyFileHeader()) {
            throwBadKeyFileError();
        }
    }

    private boolean confirmKeyFileHeader() {
        // check the headers
        if (!keyFileColumns.containsKey("Flowcell")) {
            return false;
        }
        if (!keyFileColumns.containsKey("Lane")) {
            return false;
        }
        if (!keyFileColumns.containsKey("Barcode")) {
            return false;
        }
        if (!keyFileColumns.containsKey("Sample")) {
            return false;
        }
        if (!keyFileColumns.containsKey("LibPrepID")) {
            return false;
        }

        // check the column order (needed for ParseBarcodeReads)
        if (!keyFileColumns.get("Flowcell").equals(0)) {
            return false;
        }
        if (!keyFileColumns.get("Lane").equals(1)) {
            return false;
        }
        if (!keyFileColumns.get("Barcode").equals(2)) {
            return false;
        }
        if (!keyFileColumns.get("Sample").equals(3)) {
            return false;
        }
        if (!keyFileColumns.get("LibPrepID").equals(7)) {
            return false;
        }
        return true;
    }

    private void throwBadKeyFileError() {
        String badKeyFileMessage
                = "\n\nThe keyfile does not conform to expections.\n"
                + "It must contain columns with the following (exact) headers, and\n"
                + "in the indicated columns:\n"
                + "   \"Flowcell\"                 (column A)\n"
                + "   \"Lane\"                     (column B)\n"
                + "   \"Barcode\"                  (column C)\n"
                + "   \"DNASample\" or \"Sample\"    (column D)\n"
                + "   \"LibraryPrepID\"            (column H)\n"
                + "\n\n";
        try {
            throw new IllegalStateException(badKeyFileMessage);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void populateKeyFileFields(String keyFileLine) {
        keyFileLine.trim();
        String[] cells = keyFileLine.split("\\t");

        String sample = cells[keyFileColumns.get("Sample")];
        String libPrepID = cells[keyFileColumns.get("LibPrepID")];
        String fullName = sample + ":" + cells[keyFileColumns.get("Flowcell")] + ":" + cells[keyFileColumns.get("Lane")] + ":" + libPrepID;
        rawReadCountsForFullSampleName.put(fullName, 0);
        matchedReadCountsForFullSampleName.put(fullName, 0);
        fullNameToHDF5Name.put(fullName, sample + ":" + libPrepID);

        String flowCell_Lane = cells[keyFileColumns.get("Flowcell")] + "_" + cells[keyFileColumns.get("Lane")];
        flowcellLanesInKey.add(flowCell_Lane);
        flowCellLaneToLibPrepIDs.put(flowCell_Lane, libPrepID);

        String prevSample = libraryPrepIDToSampleName.get(libPrepID);
        if (prevSample == null) {
            libraryPrepIDToSampleName.put(libPrepID, sample);
        } else if (!prevSample.contentEquals(sample)) {
            try {
                throw new IllegalStateException("\nThe key file contains different Sample names (\"" + prevSample + "\" and \"" + sample + "\") for the sample LibraryPrepID (" + libPrepID + ")\n\n");
            } catch (Exception e) {
                myLogger.error("Error in key file: " + e);
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    private void matchKeyFileToAvailableRawSeqFiles() {
        String message
                = "\nThe following raw sequence files in the input directory conform to one of our file naming conventions and have corresponding samples in the barcode key file:";
        for (int fileNum = 0; fileNum < myRawSeqFileNames.length; fileNum++) {
            String[] flowcellLane = parseRawSeqFileName(myRawSeqFileNames[fileNum]);
            if (flowcellLane != null && flowcellLanesInKey.contains(flowcellLane[0] + "_" + flowcellLane[1])) {
                seqFileNameToFlowcellLane.put(myRawSeqFileNames[fileNum], flowcellLane[0] + "_" + flowcellLane[1]);
                seqFilesInKeyAndDir.add(myRawSeqFileNames[fileNum]);
                message += "\n  " + myRawSeqFileNames[fileNum];
            }
        }
        message += "\n";
        myLogger.info(message);
    }

    private ParseBarcodeRead setUpBarcodes(int fileNum) {
        System.gc();
        String message = "\nWorking on GBS raw sequence file: " + myRawSeqFileNames[fileNum];
        fastq = true;
        if (new File(myRawSeqFileNames[fileNum]).getName().contains("qseq")) {
            fastq = false;
        }
        if (fastq) {
            message += "\n\tThis file is assumed to be in fastq format";
        } else {
            message += "\n\tThis file contains 'qseq' in its name so is assumed to be in qseq format";
        }
        String[] flowcellLane = parseRawSeqFileName(myRawSeqFileNames[fileNum]);
        if (flowcellLane == null) {
            myLogger.info(message);
            return null;
        } else {
            ParseBarcodeRead thePBR = new ParseBarcodeRead(keyFile(), enzyme(), flowcellLane[0], flowcellLane[1]);
            message += "\nTotal barcodes found in key file for this lane:" + thePBR.getBarCodeCount() + "\n";
            myLogger.info(message);
            return thePBR;
        }
    }

    /**
     * Parses out the flowcell and lane from the raw GBS sequence filename
     * (fastq or qseq file)
     *
     * @param rawSeqFileName
     * @return String[2] where element[0]=flowcell and element[1]=lane
     */
    private String[] parseRawSeqFileName(String rawSeqFileName) {
        File rawSeqFile = new File(rawSeqFileName);
        String[] FileNameParts = rawSeqFile.getName().split("_");
        if (FileNameParts.length == 3) {
            return new String[]{FileNameParts[0], FileNameParts[1]};
        } else if (FileNameParts.length == 4) {
            return new String[]{FileNameParts[0], FileNameParts[2]};
        } else if (FileNameParts.length == 5) {
            return new String[]{FileNameParts[1], FileNameParts[3]};
        } else {
            printFileNameConventions(rawSeqFileName);
            return null;
        }
    }

    private void setUpGenotypeTableBuilder() {
        genoMergeRule = new BasicGenotypeMergeRule(aveSeqErrorRate());
        File hdf5File = new File(outputHDF5GenotypesFile());
        if (hdf5File.exists()) {
            myLogger.info("\nGenotypes will be added to existing HDF5 file:\n  " + outputHDF5GenotypesFile() + "\n");
            genos = GenotypeTableBuilder.mergeTaxaIncremental(outputHDF5GenotypesFile(), genoMergeRule);
        } else {
            myLogger.info("\nThe target HDF5 file:\n  " + outputHDF5GenotypesFile()
                    + "\ndoes not exist. A new HDF5 file of that name will be created \nto hold the genotypes from this run.");
            genos = GenotypeTableBuilder.getTaxaIncrementalWithMerging(outputHDF5GenotypesFile(), myPositionList, genoMergeRule);
        }
    }

    /**
     * Gets an ArrayList of taxa, each named "Sample:LibraryPrepID" and annotated
     * with Flowcell_Lane as well as all other annotations in the key file, for 
     * the corresponding fastq file.
     *
     * @return ArrayList<Taxon>
     */
    private ArrayList<Taxon> getHDF5Taxa(int fileNum) {
        String currFlowcellLane = seqFileNameToFlowcellLane.get(myRawSeqFileNames[fileNum]);
        String[] flowcellLane = currFlowcellLane.split("_");
        TaxaList annoTL = TaxaListIOUtils.readTaxaAnnotationFile(keyFile(), "LibraryPrepID",
                ImmutableMap.of("Flowcell", flowcellLane[0], "Lane", flowcellLane[1]), false);
        ArrayList<Taxon> taxaAL = new ArrayList();
        for (Taxon tax : annoTL) {
            GeneralAnnotation annotation = tax.getAnnotation();
            String newName = annotation.getTextAnnotation("Sample").length == 0 ? annotation.getTextAnnotation("DNASample")[0] : annotation.getTextAnnotation("Sample")[0];
            String libPrepID = tax.getName();
            newName += ":" + libPrepID;
            Taxon gbsTaxon = new Taxon.Builder(tax).name(newName).addAnno("Flowcell_Lane", currFlowcellLane)
                    .addAnno("LibraryPrepID", libPrepID).addAnno("Status", "private").build();
            taxaAL.add(gbsTaxon);
        }
        return taxaAL;
    }

    private PositionList getUniquePositions() {
        //todo Move this code to TOPM
        myLogger.info("\nCounting sites in TOPM file");
        PositionListBuilder plb = new PositionListBuilder();
        chromosomes = topm.getChromosomes();
        for (Integer chrNum : chromosomes) {
            Chromosome chr = new Chromosome(chrNum.toString());
            for (int pos : topm.getUniquePositions(topm.getChromosomeIndex(chrNum))) {
                Position p = new GeneralPosition.Builder(chr, pos).build();
                plb.add(p);
            }
        }
        PositionList pl = plb.sortPositions().build();
        myLogger.info("In total, the TOPM contains " + chromosomes.length + " chromosomes and " + pl.numberOfSites() + " sites.");
        return pl;
    }

    private void generateFastSiteLookup(PositionList pl) {
        ImmutableTable.Builder ptsB = new ImmutableTable.Builder<Chromosome, Integer, Integer>();
        for (int i = 0; i < pl.numberOfSites(); i++) {
            Position position = pl.get(i);
            ptsB.put(position.getChromosome().getChromosomeNumber(), position.getPosition(), i);
        }
        positionToSite = ptsB.build();
    }

    private BufferedReader getBufferedReaderForRawSeqFile(int fileNum) {
        BufferedReader br = null;
        try {
            if (myRawSeqFileNames[fileNum].endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(myRawSeqFileNames[fileNum]))));
            } else {
                br = new BufferedReader(new FileReader(myRawSeqFileNames[fileNum]), 65536);
            }
        } catch (Exception e) {
            myLogger.error("Catch in getBufferedReader(): e=" + e);
            e.printStackTrace();
        }
        return br;
    }

    private ReadBarcodeResult readSequenceRead(BufferedReader br, String temp, ParseBarcodeRead thePBR, int[] counters) {
        ReadBarcodeResult rr = null;
        String sl = "";
        try {
            if (fastq) {
                sl = br.readLine();    // read the 2nd line in the set of 4 lines = sequence
                temp = br.readLine();  // skip the 3rd line
                temp = br.readLine();  // skip the 4th in the set of 4 lines = quality score (note that the QS is scaled differently in Cassava 1.8 - we don't use it so it is not corrected here)
                rr = thePBR.parseReadIntoTagAndTaxa(sl, null, true, 0);
            } else {  // qseq
                String[] jj = temp.split("\\s");
                sl = jj[8];
                rr = thePBR.parseReadIntoTagAndTaxa(sl, null, false, 0);
            }
        } catch (Exception e) {
            myLogger.error("Catch in readSequenceRead() at nReads=" + counters[0] + " e=" + e);
            myLogger.error("Last line read:\n" + temp + "\n");
            e.printStackTrace();
        }
        counters[0]++;  // allReads
        return rr;
    }

    private int findBestImperfectMatch(long[] read, int[] counters) {
        // this method is not ready for prime time -- to resolve a tie, it currently chooses a random tag out of the tied tags
        int tagIndex = -1;
        TagMatchFinder tmf = new TagMatchFinder(topm);
//        TreeMap<Integer, Integer> bestHitsAndDiv = tmf.findMatchesWithIntLengthWords(read, maxDivergence(), true);
        TreeMap<Integer, Integer> bestHitsAndDiv = tmf.findMatchesWithIntLengthWords(read, 0, true);
        if (bestHitsAndDiv.size() > 0) {
            counters[4]++; // imperfectMatches
            if (bestHitsAndDiv.size() == 1) {
                counters[5]++; // singleImperfectMatches
            }
            tagIndex = bestHitsAndDiv.firstKey();  // a random tag (firstKey) chosen to resolve the tie = suboptimal behavior
        }
        return tagIndex;
    }

    private void incrementDepthForTagVariants(int tagIndex, int[][] alleleDepths, int increment) {
        int chromosome = topm.getChromosome(tagIndex);
        if (chromosome == TOPMInterface.INT_MISSING) {
            return;
        }
        int startPos = topm.getStartPosition(tagIndex);
        for (int variant = 0; variant < topm.getMaxNumVariants(); variant++) {
            byte newBase = topm.getVariantDef(tagIndex, variant);
            if ((newBase == TOPMInterface.BYTE_MISSING) || (newBase == GenotypeTable.UNKNOWN_ALLELE)) {
                continue;
            }
            int offset = topm.getVariantPosOff(tagIndex, variant);
            int pos = startPos + offset;
//            int currSite = genos.getSiteOfPhysicalPosition(pos, locus);
            int currSite = positionToSite.get(chromosome, pos);
            if (currSite < 0) {
                continue;
            }
            alleleDepths[newBase][currSite] += increment;
        }
    }

    private void callGenotypes() {
        myLogger.info("\nCalling genotypes...");
        for (int currTaxonIndex = 0; currTaxonIndex < obsTagsForEachTaxon.length; currTaxonIndex++) {
            IntArrayList currTagList = obsTagsForEachTaxon[currTaxonIndex];
            currTagList.sort();
            int[][] alleleDepths = new int[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][myPositionList.numberOfSites()];
            int prevTag = currTagList.getQuick(0);
            int currInc = 0;
            for (int t = 0; t < currTagList.size(); t++) {
                int tag = currTagList.getQuick(t);
                if (tag == prevTag) {
                    currInc++;
                } else {
                    incrementDepthForTagVariants(prevTag, alleleDepths, currInc);
                    prevTag = tag;
                    currInc = 1;
                }
            }
            incrementDepthForTagVariants(prevTag, alleleDepths, currInc);
            byte[][] byteDepths = AlleleDepthUtil.depthIntToByte(alleleDepths);
            byte[] taxonGenos = resolveGenosForTaxon(byteDepths);
            if (noDepthToOutput()) {
                genos.addTaxon(taxaList.get(currTaxonIndex), taxonGenos, null);
            } else {
                genos.addTaxon(taxaList.get(currTaxonIndex), taxonGenos, byteDepths);
            }
            myLogger.info("  finished calling genotypes for " + taxaList.get(currTaxonIndex).getName());
        }
        myLogger.info("Finished calling genotypes for " + obsTagsForEachTaxon.length + " taxa\n");
    }

    private byte[] resolveGenosForTaxon(byte[][] depthsForTaxon) {
        int nAlleles = depthsForTaxon.length;
        byte[] depthsAtSite = new byte[nAlleles];
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

    private void writeReadsPerSampleReports() {
        myLogger.info("\nWriting ReadsPerSample log file...");
        String outFileS = myOutputDir + File.separator + (new File(keyFile())).getName();
        outFileS = outFileS.replaceAll(".txt", "_ReadsPerSample.log");
        outFileS = outFileS.replaceAll("_key", "");
        try {
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outFileS))), 65536);
            bw.write("FullSampleName\tgoodBarcodedReads\tgoodReadsMatchedToTOPM\n");
            for (String fullSampleName : rawReadCountsForFullSampleName.keySet()) {
                bw.write(fullSampleName + "\t" + rawReadCountsForFullSampleName.get(fullSampleName) + "\t" + matchedReadCountsForFullSampleName.get(fullSampleName) + "\n");
            }
            bw.close();
        } catch (Exception e) {
            myLogger.error("Couldn't write to ReadsPerSample log file: " + e);
            e.printStackTrace();
            System.exit(1);
        }
        myLogger.info("   ...done\n");
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

    public static final String rawSeqFileNameRegex
            = "(?i)" + // case insensitve
            ".*\\.fq" + "$|"
            + ".*\\.fq\\.gz" + "$|"
            + ".*\\.fastq" + "$|"
            + ".*_fastq\\.txt" + "$|"
            + ".*_fastq\\.gz" + "$|"
            + ".*_fastq\\.txt\\.gz" + "$|"
            + ".*_sequence\\.txt" + "$|"
            + ".*_sequence\\.txt\\.gz" + "$|"
            + ".*_qseq\\.txt" + "$|"
            + ".*_qseq\\.txt\\.gz" + "$";
    //                \\. denotes escape . so it doesn't mean 'any char'

    private String noMatchingRawSeqFileNamesMessage
            = "Couldn't find any files that end with "
            + "\".fq\", "
            + "\".fq.gz\", "
            + "\".fastq\", "
            + "\"_fastq.txt\", "
            + "\"_fastq.gz\", "
            + "\"_fastq.txt.gz\", "
            + "\"_sequence.txt\", "
            + "\"_sequence.txt.gz\", "
            + "\"_qseq.txt\", or "
            + "\"_qseq.txt.gz\" "
            + "in the supplied directory: ";

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
    //     GeneratePluginCode.generate(ProductionSNPCallerPlugin.class);
    // }
    /**
     * Input directory containing fastq AND/OR qseq files.
     *
     * @return Input Directory
     */
    public String inputDirectory() {
        return myInputDirectory.value();
    }

    /**
     * Set Input Directory. Input directory containing fastq AND/OR qseq files.
     *
     * @param value Input Directory
     *
     * @return this plugin
     */
    public ProductionSNPCallerPlugin inputDirectory(String value) {
        myInputDirectory = new PluginParameter<>(myInputDirectory, value);
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
     * Set Key File. Key file listing barcodes distinguishing the samples
     *
     * @param value Key File
     *
     * @return this plugin
     */
    public ProductionSNPCallerPlugin keyFile(String value) {
        myKeyFile = new PluginParameter<>(myKeyFile, value);
        return this;
    }

    /**
     * Enzyme used to create the GBS library
     *
     * @return Enzyme
     */
    public String enzyme() {
        return myEnzyme.value();
    }

    /**
     * Set Enzyme. Enzyme used to create the GBS library
     *
     * @param value Enzyme
     *
     * @return this plugin
     */
    public ProductionSNPCallerPlugin enzyme(String value) {
        myEnzyme = new PluginParameter<>(myEnzyme, value);
        return this;
    }

    /**
     * Physical map file containing tags and corresponding variants (production
     * TOPM)
     *
     * @return Input TOPM File
     */
    public String inputTOPMFile() {
        return myProductionTOPM.value();
    }

    /**
     * Set Input TOPM File. Physical map file containing tags and corresponding
     * variants (production TOPM)
     *
     * @param value Input TOPM File
     *
     * @return this plugin
     */
    public ProductionSNPCallerPlugin inputTOPMFile(String value) {
        myProductionTOPM = new PluginParameter<>(myProductionTOPM, value);
        return this;
    }

    /**
     * Output (target) HDF5 genotypes file to add new genotypes to (new file
     * created if it doesn't exist)
     *
     * @return Output HDF5 Genotypes File
     */
    public String outputHDF5GenotypesFile() {
        return myOutputGenotypes.value();
    }

    /**
     * Set Output HDF5 Genotypes File. Output (target) HDF5 genotypes file to
     * add new genotypes to (new file created if it doesn't exist)
     *
     * @param value Output HDF5 Genotypes File
     *
     * @return this plugin
     */
    public ProductionSNPCallerPlugin outputHDF5GenotypesFile(String value) {
        myOutputGenotypes = new PluginParameter<>(myOutputGenotypes, value);
        return this;
    }

    /**
     * Average sequencing error rate per base (used to decide between
     * heterozygous and homozygous calls)
     *
     * @return Ave Seq Error Rate
     */
    public Double aveSeqErrorRate() {
        return myAveSeqErrorRate.value();
    }

    /**
     * Set Ave Seq Error Rate. Average sequencing error rate per base (used to
     * decide between heterozygous and homozygous calls)
     *
     * @param value Ave Seq Error Rate
     *
     * @return this plugin
     */
    public ProductionSNPCallerPlugin aveSeqErrorRate(Double value) {
        myAveSeqErrorRate = new PluginParameter<>(myAveSeqErrorRate, value);
        return this;
    }

//    /**
//     * Maximum divergence (edit distance) between new read and previously mapped
//     * read (Default: 0 = perfect matches only)
//     *
//     * @return Max Divergence
//     */
//    public Integer maxDivergence() {
//        return myMaxDivergence.value();
//    }
//
//    /**
//     * Set Max Divergence. Maximum divergence (edit distance) between new read
//     * and previously mapped read (Default: 0 = perfect matches only)
//     *
//     * @param value Max Divergence
//     *
//     * @return this plugin
//     */
//    public ProductionSNPCallerPlugin maxDivergence(Integer value) {
//        myMaxDivergence = new PluginParameter<>(myMaxDivergence, value);
//        return this;
//    }

    /**
     * Keep hdf5 genotypes open for future runs that add more taxa or more depth
     *
     * @return Keep Genotypes Open
     */
    public Boolean keepGenotypesOpen() {
        return myKeepGenotypesOpen.value();
    }

    /**
     * Set Keep Genotypes Open. Keep hdf5 genotypes open for future runs that
     * add more taxa or more depth
     *
     * @param value Keep Genotypes Open
     *
     * @return this plugin
     */
    public ProductionSNPCallerPlugin keepGenotypesOpen(Boolean value) {
        myKeepGenotypesOpen = new PluginParameter<>(myKeepGenotypesOpen, value);
        return this;
    }

    /**
     * No depth output: do not write depths to the output hdf5 genotypes file
     *
     * @return No Depth to Output
     */
    public Boolean noDepthToOutput() {
        return myNoDepthOutput.value();
    }

    /**
     * Set No Depth to Output. No depth output: do not write depths to the
     * output hdf5 genotypes file
     *
     * @param value No Depth to Output
     *
     * @return this plugin
     */
    public ProductionSNPCallerPlugin noDepthToOutput(Boolean value) {
        myNoDepthOutput = new PluginParameter<>(myNoDepthOutput, value);
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
