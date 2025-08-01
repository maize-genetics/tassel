/*
 * DiscoverySNPCallerPlugin
 */
package net.maizegenetics.analysis.gbs;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import net.maizegenetics.dna.map.TOPMInterface;
import net.maizegenetics.dna.map.TagLocus;
import net.maizegenetics.dna.map.TagsOnPhysicalMap;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxaByteFileMap;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.nbio.core.util.ConcurrencyTools;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.Arrays;
import java.util.HashMap;

import static net.maizegenetics.dna.snp.NucleotideAlignmentConstants.*;
import net.maizegenetics.plugindef.PluginParameter;

/**
 * This class aligns tags at the same physical location against one another,
 * calls SNPs, and then outputs the SNPs to a HapMap file.
 *
 * It is multi-threaded, as there are substantial speed increases with it.
 *
 * @author Jeff Glaubitz
 * @author Ed Buckler
 */
public class DiscoverySNPCallerPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(DiscoverySNPCallerPlugin.class);

    private PluginParameter<String> myInputTagsByTaxa = new PluginParameter.Builder<>("i", null, String.class).guiName("Input Tags by Taxa File").required(true).inFile()
            .description("Input TagsByTaxa file (if hdf5 format, use .hdf or .h5 extension)").build();
    private PluginParameter<Boolean> myUseByteFormat = new PluginParameter.Builder<>("y", false, Boolean.class).guiName("Use Byte Format")
            .description("Use byte format TagsByTaxa file (*.tbt.byte)").build();
    private PluginParameter<String> myInputTOPM = new PluginParameter.Builder<>("m", null, String.class).guiName("Input TOPM File").required(true).inFile()
            .description("TagsOnPhysicalMap (TOPM) file containing genomic positions of tags").build();
    private PluginParameter<String> myOutputTOPM = new PluginParameter.Builder<>("o", null, String.class).guiName("Output TOPM File").required(true).outFile()
            .description("Output TagsOnPhysicalMap (TOPM) file with allele calls (variants) added (for use with ProductionSNPCallerPlugin").build();
    private PluginParameter<String> myLogFile = new PluginParameter.Builder<>("log", null, String.class).guiName("Log File").outFile()
            .description("TagLocus log file name. (Default: TagLocusLog.txt)").build();
    private PluginParameter<Double> myMinF = new PluginParameter.Builder<>("mnF", -2.0, Double.class).guiName("Minimum F")
            .description("Minimum F (inbreeding coefficient)").build();
    private PluginParameter<String> myPedigreeFile = new PluginParameter.Builder<>("p", null, String.class).guiName("Pedigree File").inFile()
            .description("Pedigree file containing full sample names (or expected names after merging) & expected inbreeding "
                    + "coefficient (F) for each.  Only taxa with expected F >= mnF used to calculate F = 1-Ho/He. "
                    + "(default: use ALL taxa to calculate F").build();
    private PluginParameter<Double> myMinMinorAlleleFreq = new PluginParameter.Builder<>("mnMAF", 0.01, Double.class).guiName("Min Minor Allele Freq")
            .description("Minimum minor allele frequency").build();
    private PluginParameter<Integer> myMinMinorAlleleCount = new PluginParameter.Builder<>("mnMAC", 10, Integer.class).guiName("Min Minor Allele Count")
            .description("Minimum minor allele count").build();
    private PluginParameter<Double> myMinLocusCoverage = new PluginParameter.Builder<>("mnLCov", 0.1, Double.class).guiName("Min Locus Coverage")
            .description("Minimum locus coverage (proportion of Taxa with a genotype)").build();
    private PluginParameter<Double> myAveSeqErrorRate = new PluginParameter.Builder<>("eR", 0.01, Double.class).guiName("Ave Seq Error Rate")
            .description("Average sequencing error rate per base (used to decide between heterozygous and homozygous calls)").build();
    private PluginParameter<String> myRefGenome = new PluginParameter.Builder<>("ref", null, String.class).guiName("Reference Genome File").inFile()
            .description("Path to reference genome in fasta format. Ensures that a tag from the reference genome is always included "
                    + "when the tags at a locus are aligned against each other to call SNPs. DEFAULT: Don't use reference genome.").build();
    private PluginParameter<Integer> myStartChr = new PluginParameter.Builder<>("sC", null, Integer.class).guiName("Start Chromosome").required(true)
            .description("Start Chromosome").build();
    private PluginParameter<Integer> myEndChr = new PluginParameter.Builder<>("eC", null, Integer.class).guiName("End Chromosome").required(true)
            .description("End Chromosome").build();
    private PluginParameter<Boolean> myIncludeRareAlleles = new PluginParameter.Builder<>("inclRare", false, Boolean.class).guiName("Include Rare Alleles")
            .description("Include the rare alleles at site (3 or 4th states)").build();
    private PluginParameter<Boolean> myIncludeGaps = new PluginParameter.Builder<>("inclGaps", false, Boolean.class).guiName("Include Gaps")
            .description("Include sites where major or minor allele is a GAP").build();
    private PluginParameter<Boolean> myCallBiSNPsWGap = new PluginParameter.Builder<>("callBiSNPsWGap", false, Boolean.class).guiName("Call Biallelic SNPs with Gap")
            .description("Include sites where the third allele is a GAP (mutually exclusive with inclGaps)").build();

    private TagsOnPhysicalMap theTOPM = null;
    private TagsByTaxa theTBT = null;
    private boolean usePedigree = false;
    private HashMap<String, Double> taxaFs = null;
    private boolean[] useTaxaForMinF = null;
    private int nInbredTaxa = Integer.MIN_VALUE;
    private int minTaxaWithLocus;
    private boolean includeReference = false;
    private long[] refGenomeChr = null;
    private boolean fuzzyStartPositions = false;
    private int locusBorder = 0;
    private final static int CHR = 0, STRAND = 1, START_POS = 2;  // indices of these position attributes in array returned by theTOPM.getPositionArray(a)
    private boolean customSNPLogging = true;  // a custom SNP log that collects useful info for filtering SNPs through machine learning criteria
    private CustomSNPLog myCustomSNPLog = null;
    private boolean customFiltering = false;

    public DiscoverySNPCallerPlugin() {
        super(null, false);
    }

    public DiscoverySNPCallerPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        myLogger.info("Finding SNPs in " + inputTagsByTaxaFile() + ".");
        myLogger.info(String.format("StartChr:%d EndChr:%d %n", startChromosome(), endChromosome()));
        theTOPM.sortTable(true);
        myLogger.info("\nAs a check, here are the first 5 tags in the TOPM (sorted by position):");
        theTOPM.printRows(5, true, true);
        DataOutputStream locusLogDOS = openLocusLog(logFile());
        if (customSNPLogging) {
            myCustomSNPLog = new CustomSNPLog(logFile());
        }
        for (int chr = startChromosome(); chr <= endChromosome(); chr++) {
            myLogger.info("\n\nProcessing chromosome " + chr + "...");
            if (includeReference) {
                refGenomeChr = readReferenceGenomeChr(referenceGenomeFile(), chr);
                if (refGenomeChr == null) {
                    myLogger.info("  WARNING: chromosome " + chr + " not found in the reference genome file. Skipping this chromosome.");
                    continue;
                }
            }
            discoverSNPsOnChromosome(chr, locusLogDOS);
            myLogger.info("Finished processing chromosome " + chr + "\n\n");
        }

        if (outputTOPMFile().endsWith(".txt")) {
            theTOPM.writeTextFile(new File(outputTOPMFile()));
        } else {
            theTOPM.writeBinaryFile(new File(outputTOPMFile()));
        }

        try {
            locusLogDOS.close();
        } catch (Exception e) {
            catchLocusLogException(e);
        }
        if (customSNPLogging) {
            myCustomSNPLog.close();
        }
        ConcurrencyTools.shutdown();
        return null;
    }

    @Override
    public void postProcessParameters() {

        if (myInputTagsByTaxa.isEmpty()) {
            throw new IllegalArgumentException("DiscoverySNPCallerPlugin: postProcessParameters: Input Tags by Taxa File not Set.");
        } else {
            if (inputTagsByTaxaFile().endsWith(".hdf") || inputTagsByTaxaFile().endsWith(".h5")) {
                theTBT = new TagsByTaxaByteHDF5TagGroups(inputTagsByTaxaFile());
            } else if (useByteFormat()) {
                theTBT = new TagsByTaxaByteFileMap(inputTagsByTaxaFile());
            }
        }
        if (theTBT == null) {
            throw new IllegalArgumentException("DiscoverySNPCallerPlugin: postProcessParameters: Problem reading Tags by Taxa File: " + inputTagsByTaxaFile());
        }

        boolean loadBinary = (inputTOPMFile().endsWith(".txt")) ? false : true;
        theTOPM = new TagsOnPhysicalMap(inputTOPMFile(), loadBinary);

        if (myLogFile.isEmpty()) {
            try {
                File outDir = (new File(outputTOPMFile())).getCanonicalFile().getParentFile();
                logFile(outDir.getCanonicalPath() + File.separator + "TagLocusLog.txt");
            } catch (IOException e) {
                throw new IllegalArgumentException("Problem creating the tagLocusLog file. Program aborted: "+e);
            }
        }

        if (!myPedigreeFile.isEmpty()) {
            taxaFs = readTaxaFsFromFile(new File(pedigreeFile()));
            if (taxaFs == null) {
                throw new IllegalArgumentException("Problem reading the pedigree file. Progam aborted.");
            }
            if (!maskNonInbredTaxa()) {
                throw new IllegalArgumentException("Mismatch between taxa names in the pedigree file and TBT. Progam aborted.");
            }
            usePedigree = true;
        }

        minTaxaWithLocus = (int) Math.round(theTBT.getTaxaCount() * minLocusCoverage());

        if (!myRefGenome.isEmpty()) {
            includeReference = true;
        }

        if (callBiallelicSNPsWithGap() && includeGaps()) {
            throw new IllegalArgumentException("The callBiSNPsWGap option is mutually exclusive with the inclGaps option.");
        }

        if (endChromosome() - startChromosome() < 0) {
            throw new IllegalArgumentException("The start chromosome is larger than the end chromosome.");
        }

        myLogger.info(String.format("minTaxaWithLocus:%d MinF:%g MinMAF:%g MinMAC:%d %n", minTaxaWithLocus, minimumF(), minMinorAlleleFreq(), minMinorAlleleCount()));
        myLogger.info(String.format("includeRare:%s includeGaps:%s %n", includeRareAlleles(), includeGaps()));
    }

    //
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(DiscoverySNPCallerPlugin.class);
    // }
    //
    /**
     * Input TagsByTaxa file (if hdf5 format, use .hdf or .h5 extension)
     *
     * @return Input Tags by Taxa File
     */
    public String inputTagsByTaxaFile() {
        return myInputTagsByTaxa.value();
    }

    /**
     * Set Input Tags by Taxa File. Input TagsByTaxa file (if hdf5 format, use
     * .hdf or .h5 extension)
     *
     * @param value Input Tags by Taxa File
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin inputTagsByTaxaFile(String value) {
        myInputTagsByTaxa = new PluginParameter<>(myInputTagsByTaxa, value);
        return this;
    }

    /**
     * Use byte format TagsByTaxa file (*.tbt.byte)
     *
     * @return Use Byte Format
     */
    public Boolean useByteFormat() {
        return myUseByteFormat.value();
    }

    /**
     * Set Use Byte Format. Use byte format TagsByTaxa file (*.tbt.byte)
     *
     * @param value Use Byte Format
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin useByteFormat(Boolean value) {
        myUseByteFormat = new PluginParameter<>(myUseByteFormat, value);
        return this;
    }

    /**
     * TagsOnPhysicalMap (TOPM) file containing genomic positions of tags
     *
     * @return Input TOPM File
     */
    public String inputTOPMFile() {
        return myInputTOPM.value();
    }

    /**
     * Set Input TOPM File. TagsOnPhysicalMap (TOPM) file containing genomic
     * positions of tags
     *
     * @param value Input TOPM File
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin inputTOPMFile(String value) {
        myInputTOPM = new PluginParameter<>(myInputTOPM, value);
        return this;
    }

    /**
     * Output TagsOnPhysicalMap (TOPM) file with allele calls (variants) added
     * (for use with ProductionSNPCallerPlugin
     *
     * @return Output TOPM File
     */
    public String outputTOPMFile() {
        return myOutputTOPM.value();
    }

    /**
     * Set Output TOPM File. Output TagsOnPhysicalMap (TOPM) file with allele
     * calls (variants) added (for use with ProductionSNPCallerPlugin
     *
     * @param value Output TOPM File
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin outputTOPMFile(String value) {
        myOutputTOPM = new PluginParameter<>(myOutputTOPM, value);
        return this;
    }

    /**
     * TagLocus log file name
     *
     * @return Log File
     */
    public String logFile() {
        return myLogFile.value();
    }

    /**
     * Set Log File. TagLocus log file name
     *
     * @param value Log File
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin logFile(String value) {
        myLogFile = new PluginParameter<>(myLogFile, value);
        return this;
    }

    /**
     * Minimum F (inbreeding coefficient)
     *
     * @return Minimum F
     */
    public Double minimumF() {
        return myMinF.value();
    }

    /**
     * Set Minimum F. Minimum F (inbreeding coefficient)
     *
     * @param value Minimum F
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin minimumF(Double value) {
        myMinF = new PluginParameter<>(myMinF, value);
        return this;
    }

    /**
     * Pedigree file containing full sample names (or expected names after
     * merging) & expected inbreeding coefficient (F) for each. Only taxa with
     * expected F >= mnF used to calculate F = 1-Ho/He. (default: use ALL taxa
     * to calculate F
     *
     * @return Pedigree File
     */
    public String pedigreeFile() {
        return myPedigreeFile.value();
    }

    /**
     * Set Pedigree File. Pedigree file containing full sample names (or
     * expected names after merging) & expected inbreeding coefficient (F) for
     * each. Only taxa with expected F >= mnF used to calculate F = 1-Ho/He.
     * (default: use ALL taxa to calculate F
     *
     * @param value Pedigree File
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin pedigreeFile(String value) {
        myPedigreeFile = new PluginParameter<>(myPedigreeFile, value);
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
    public DiscoverySNPCallerPlugin minMinorAlleleFreq(Double value) {
        myMinMinorAlleleFreq = new PluginParameter<>(myMinMinorAlleleFreq, value);
        return this;
    }

    /**
     * Minimum minor allele count
     *
     * @return Min Minor Allele Count
     */
    public Integer minMinorAlleleCount() {
        return myMinMinorAlleleCount.value();
    }

    /**
     * Set Min Minor Allele Count. Minimum minor allele count
     *
     * @param value Min Minor Allele Count
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin minMinorAlleleCount(Integer value) {
        myMinMinorAlleleCount = new PluginParameter<>(myMinMinorAlleleCount, value);
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
     * Set Min Locus Coverage. Minimum locus coverage (proportion of Taxa with a
     * genotype)
     *
     * @param value Min Locus Coverage
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin minLocusCoverage(Double value) {
        myMinLocusCoverage = new PluginParameter<>(myMinLocusCoverage, value);
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
    public DiscoverySNPCallerPlugin aveSeqErrorRate(Double value) {
        myAveSeqErrorRate = new PluginParameter<>(myAveSeqErrorRate, value);
        return this;
    }

    /**
     * Path to reference genome in fasta format. Ensures that a tag from the
     * reference genome is always included when the tags at a locus are aligned
     * against each other to call SNPs. The reference allele for each site is
     * then provided in the output HapMap files, under the taxon name
     * "REFERENCE_GENOME" (first taxon). DEFAULT: Don't use reference genome.
     *
     * @return Reference Genome File
     */
    public String referenceGenomeFile() {
        return myRefGenome.value();
    }

    /**
     * Set Reference Genome File. Path to reference genome in fasta format.
     * Ensures that a tag from the reference genome is always included when the
     * tags at a locus are aligned against each other to call SNPs. The
     * reference allele for each site is then provided in the output HapMap
     * files, under the taxon name "REFERENCE_GENOME" (first taxon). DEFAULT:
     * Don't use reference genome.
     *
     * @param value Reference Genome File
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin referenceGenomeFile(String value) {
        myRefGenome = new PluginParameter<>(myRefGenome, value);
        return this;
    }

    /**
     * Include the rare alleles at site (3 or 4th states)
     *
     * @return Include Rare Alleles
     */
    public Boolean includeRareAlleles() {
        return myIncludeRareAlleles.value();
    }

    /**
     * Set Include Rare Alleles. Include the rare alleles at site (3 or 4th
     * states)
     *
     * @param value Include Rare Alleles
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin includeRareAlleles(Boolean value) {
        myIncludeRareAlleles = new PluginParameter<>(myIncludeRareAlleles, value);
        return this;
    }

    /**
     * Include sites where major or minor allele is a GAP
     *
     * @return Include Gaps
     */
    public Boolean includeGaps() {
        return myIncludeGaps.value();
    }

    /**
     * Set Include Gaps. Include sites where major or minor allele is a GAP
     *
     * @param value Include Gaps
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin includeGaps(Boolean value) {
        myIncludeGaps = new PluginParameter<>(myIncludeGaps, value);
        return this;
    }

    /**
     * Include sites where the third allele is a GAP (mutually exclusive with
     * inclGaps)
     *
     * @return Call Biallelic SNPs with Gap
     */
    public Boolean callBiallelicSNPsWithGap() {
        return myCallBiSNPsWGap.value();
    }

    /**
     * Set Call Biallelic SNPs with Gap. Include sites where the third allele is
     * a GAP (mutually exclusive with inclGaps)
     *
     * @param value Call Biallelic SNPs with Gap
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin callBiallelicSNPsWithGap(Boolean value) {
        myCallBiSNPsWGap = new PluginParameter<>(myCallBiSNPsWGap, value);
        return this;
    }

    /**
     * Start Chromosome
     *
     * @return Start Chromosome
     */
    public Integer startChromosome() {
        return myStartChr.value();
    }

    /**
     * Set Start Chromosome. Start Chromosome
     *
     * @param value Start Chromosome
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin startChromosome(Integer value) {
        myStartChr = new PluginParameter<>(myStartChr, value);
        return this;
    }

    /**
     * End Chromosome
     *
     * @return End Chromosome
     */
    public Integer endChromosome() {
        return myEndChr.value();
    }

    /**
     * Set End Chromosome. End Chromosome
     *
     * @param value End Chromosome
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPlugin endChromosome(Integer value) {
        myEndChr = new PluginParameter<>(myEndChr, value);
        return this;
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

    public void discoverSNPsOnChromosome(int targetChromo, DataOutputStream locusLogDOS) {
        long time = System.currentTimeMillis();
        TagLocus currTAL = new TagLocus(Integer.MIN_VALUE, Byte.MIN_VALUE, Integer.MIN_VALUE, Integer.MIN_VALUE, includeReference, fuzzyStartPositions, aveSeqErrorRate());
        int[] currPos = null;
        int countLoci = 0;
        int siteCnt = 0;
        for (int i = 0; i < theTOPM.getSize(); i++) {
            int ri = theTOPM.getReadIndexForPositionIndex(i);  // process tags in order of physical position
            int[] newPos = theTOPM.getPositionArray(ri);
            if (newPos[CHR] != targetChromo) {
                continue;    //Skip tags from other chromosomes
            }
            if ((fuzzyStartPositions && nearbyTag(newPos, currPos)) || Arrays.equals(newPos, currPos)) {
                currTAL.addTag(ri, theTOPM, theTBT, includeReference, fuzzyStartPositions);
            } else {
                int nTaxaCovered = currTAL.getNumberTaxaCovered();
                if (currTAL.getSize() > 1 && nTaxaCovered >= minTaxaWithLocus) {  // finish the current TAL
                    siteCnt += findSNPsInTagLocus(currTAL, locusLogDOS);  // note that with fuzzyStartPositions there may be no overlapping tags!!
                    countLoci++;
                    if (siteCnt % 100 == 0) {
                        double rate = (double) siteCnt / (double) (System.currentTimeMillis() - time);
                        myLogger.info(String.format(
                                "Chr:%d Pos:%d Loci=%d SNPs=%d rate=%g SNP/millisec", currPos[CHR], currPos[START_POS], countLoci, siteCnt, rate));
                    }
                } else if (currPos != null) {
                    logRejectedTagLocus(currTAL, locusLogDOS);
                }
                currPos = newPos; // start a new TAL with the current tag
                if ((currPos[STRAND] != TagsOnPhysicalMap.BYTE_MISSING) && (currPos[START_POS] != TOPMInterface.INT_MISSING)) {  // we already know that currPos[CHR]==targetChromo
                    currTAL = new TagLocus(currPos[CHR], (byte) currPos[STRAND], currPos[START_POS], theTOPM.getTagLength(ri), includeReference, fuzzyStartPositions, aveSeqErrorRate());
                    currTAL.addTag(ri, theTOPM, theTBT, includeReference, fuzzyStartPositions);
                } else {
                    currPos = null;  // invalid position
                }
            }
        }
        if ((currTAL.getSize() > 1) && (currTAL.getNumberTaxaCovered() >= minTaxaWithLocus)) { // then finish the final TAL for the targetChromo
            findSNPsInTagLocus(currTAL, locusLogDOS);
        } else if (currPos != null) {
            logRejectedTagLocus(currTAL, locusLogDOS);
        }
        myLogger.info("Number of marker sites recorded for chr" + targetChromo + ": " + siteCnt);
    }

    boolean nearbyTag(int[] newTagPos, int[] currTagPos) {
        if (newTagPos == null || currTagPos == null) {
            return false;
        }
        // because we move through the TOPM in positional order, the newTag startPosition is guaranteed to be >= that of the current tag
        if (newTagPos[CHR] == currTagPos[CHR] && newTagPos[START_POS] - currTagPos[START_POS] < locusBorder) {  // &&newTagPos[STRAND]==currTagPos[STRAND]
            // grab all of the tags that align to a local region (until a gap > tolerance is reached)
            currTagPos[START_POS] = newTagPos[START_POS];
            return true;
        }
        return false;
    }

    private synchronized int findSNPsInTagLocus(TagLocus theTAL, DataOutputStream locusLogDOS) {
        boolean refTagUsed = false;
        if (theTAL.getSize() < 2) {
            logRejectedTagLocus(theTAL, locusLogDOS);
            return 0;  // need at least two (overlapping!) sequences to make an alignment
        }
        byte[][] callsBySite;
        if (includeReference) {
            if (fuzzyStartPositions) {
                String refSeqInRegion = getRefSeqInRegion(theTAL);
                callsBySite = theTAL.getSNPCallsQuant(refSeqInRegion, callBiallelicSNPsWithGap());
            } else {
                addRefTag(theTAL);
                refTagUsed = true;
                callsBySite = theTAL.getSNPCallsQuant(callBiallelicSNPsWithGap(), includeReference);
            }
        } else {
            callsBySite = theTAL.getSNPCallsQuant(callBiallelicSNPsWithGap(), includeReference);
        }
        if (callsBySite == null) {
            logAcceptedTagLocus(theTAL.getLocusReport(minTaxaWithLocus, null), locusLogDOS);
            return 0;
        }
        int[] positionsInLocus = theTAL.getPositionsOfVariableSites();
        int strand = theTAL.getStrand();
        boolean[] varSiteKept = new boolean[callsBySite.length];  // initializes to false
        TagLocusSiteQualityScores SiteQualityScores = new TagLocusSiteQualityScores(callsBySite.length);
        for (int s = 0; s < callsBySite.length; s++) {
            byte[] alleles = null;
            if ((alleles = isSiteGood(callsBySite[s])) == null) {
                continue;
            }
            if (includeReference && !fuzzyStartPositions && theTAL.getRefGeno(s) == NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) {
                continue;
            }
            int position = (strand == -1) ? theTAL.getMinStartPosition() - positionsInLocus[s] : theTAL.getMinStartPosition() + positionsInLocus[s];
            CustomSNPLogRecord mySNPLogRecord;
            if (customSNPLogging) {
                mySNPLogRecord = new CustomSNPLogRecord(s, theTAL, position, useTaxaForMinF, refTagUsed);
                myCustomSNPLog.writeEntry(mySNPLogRecord.toString());
                SiteQualityScores.addSite(s, mySNPLogRecord.getInbredCoverage(), mySNPLogRecord.getInbredHetScore(), alleles, position);
                if (customFiltering && !mySNPLogRecord.isGoodSNP()) {
                    continue;
                }
            }
            varSiteKept[s] = true;
            if (!customFiltering) {
                updateTOPM(theTAL, s, position, strand, alleles);
            }
        }
        logAcceptedTagLocus(theTAL.getLocusReport(minTaxaWithLocus, varSiteKept), locusLogDOS);
        if (customFiltering) {
            updateTOPM(theTAL, varSiteKept, SiteQualityScores);
        }
        return callsBySite.length;
    }

    private void updateTOPM(TagLocus myTAL, boolean[] varSiteKept, TagLocusSiteQualityScores SiteQualityScores) {
        SiteQualityScores.sortByQuality();
        byte strand = myTAL.getStrand();
        for (int s = 0; s < SiteQualityScores.getSize(); s++) {
            int siteInTAL = SiteQualityScores.getSiteInTAL(s);
            if (varSiteKept[siteInTAL]) {
                updateTOPM(myTAL, siteInTAL, SiteQualityScores.getPosition(s), strand, SiteQualityScores.getAlleles(s));
            }
        }
    }

    private void updateTOPM(TagLocus myTAL, int variableSite, int position, int strand, byte[] alleles) {
        for (int tg = 0; tg < myTAL.getSize(); tg++) {
            int topmTagIndex = myTAL.getTOPMIndexOfTag(tg);
            if (topmTagIndex == Integer.MIN_VALUE) {
                continue; // skip the reference genome tag (which may not be in the TOPM)
            }
            byte baseToAdd = myTAL.getCallAtVariableSiteForTag(variableSite, tg);
            boolean matched = false;
            for (byte cb : alleles) {
                if (baseToAdd == cb) {
                    matched = true;
                    break;
                }
            }
            // so that all tags in the tagAlignment have the same corresponding variants in the TOPM, add a variant no matter what (set to missing if needed)
            byte offset = (byte) (position - myTAL.getMinStartPosition());
            if (!matched) {
                baseToAdd = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
            }
            if (strand == -1) {
                baseToAdd = NucleotideAlignmentConstants.getNucleotideComplement(baseToAdd);  // record everything relative to the plus strand
            }
            // convert from allele from 0-15 style to IUPAC ASCII character value (e.g., (byte) 'A') (maintains compatibility with Tassel3 TOPM)
            baseToAdd = getIUPACAllele(baseToAdd);
            theTOPM.addVariant(topmTagIndex, offset, baseToAdd);
        }
    }

    /**
     *
     * @param calls
     * @return
     */
    private byte[] isSiteGood(byte[] calls) {
        int[][] allelesByFreq = GenotypeTableUtils.getAllelesSortedByFrequency(calls);
        if (allelesByFreq[0].length < 2) {
            return null; // quantitative SNP calling rendered the site invariant
        }
        int aCnt = allelesByFreq[1][0] + allelesByFreq[1][1];
        double theMAF = (double) allelesByFreq[1][1] / (double) aCnt;
        if ((theMAF < minMinorAlleleFreq()) && (allelesByFreq[1][1] < minMinorAlleleCount())) {
            return null;  // note that a site only needs to pass one of the criteria, minMAF &/or minMAC
        }
        byte majAllele = (byte) allelesByFreq[0][0];
        byte minAllele = (byte) allelesByFreq[0][1];
        if (!includeGaps() && ((majAllele == NucleotideAlignmentConstants.GAP_ALLELE) || (minAllele == NucleotideAlignmentConstants.GAP_ALLELE))) {
            return null;
        }
        byte homMaj = (byte) ((majAllele << 4) | majAllele);
        byte homMin = (byte) ((minAllele << 4) | minAllele);
        byte hetG1 = GenotypeTableUtils.getDiploidValue(majAllele, minAllele);
        byte hetG2 = GenotypeTableUtils.getDiploidValue(minAllele, majAllele);
        if (minimumF() > -1.0) { // only test for minF if the parameter has been set above the theoretical minimum
            double obsF = calculateF(calls, allelesByFreq, hetG1, hetG2, theMAF);
            if (obsF < minimumF()) {
                return null;
            }
        }
        return getGoodAlleles(calls, allelesByFreq, homMaj, homMin, hetG1, hetG2, majAllele, minAllele);
    }

    private double calculateF(byte[] calls, int[][] alleles, byte hetG1, byte hetG2, double theMAF) {
        boolean report = false;
        double obsF;
        int hetGCnt = 0;
        if (usePedigree) {
            byte[] callsToUse = filterCallsForInbreds(calls);
            //int[][] allelesToUse = getSortedAlleleCounts(callsToUse);
            int[][] allelesToUse = GenotypeTableUtils.getAllelesSortedByFrequency(callsToUse);
            if (allelesToUse[0].length < 2) {
                return 1.0;  // lack of variation in the known inbreds will NOT reject a SNP
            }
            int aCnt = allelesToUse[1][0] + allelesToUse[1][1];
            double newMAF = (double) allelesToUse[1][1] / (double) aCnt;
            if (newMAF <= 0.0) {
                return 1.0;  // lack of variation in the known inbreds will NOT reject a SNP
            }
            byte majAllele = (byte) allelesToUse[0][0];
            byte minAllele = (byte) allelesToUse[0][1];
            //byte newHetG = IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(majAllele, minAllele);
            byte newHetG1 = GenotypeTableUtils.getDiploidValue(majAllele, minAllele);
            byte newHetG2 = GenotypeTableUtils.getDiploidValue(minAllele, majAllele);
            for (byte i : callsToUse) {
                if (i == newHetG1 || i == newHetG2) {
                    hetGCnt++;
                }
            }
            int majGCnt = (allelesToUse[1][0] - hetGCnt) / 2; // number of homozygous major allele genotypes
            int minGCnt = (allelesToUse[1][1] - hetGCnt) / 2; // number of homozygous minor allele genotypes
            double propHets = (double) hetGCnt / (double) (hetGCnt + majGCnt + minGCnt);
            double expHets = 2.0 * newMAF * (1 - newMAF);
            obsF = 1.0 - (propHets / expHets);
            if (report) {
                System.out.printf("%d %d %d propHets:%g expHets:%g obsF:%g %n", majGCnt, minGCnt, hetGCnt, propHets, expHets, obsF);
            }
            return obsF;
        } else {
            for (byte i : calls) {
                if (i == hetG1 || i == hetG2) {
                    hetGCnt++;
                }
            }
            int majGCnt = (alleles[1][0] - hetGCnt) / 2; // number of homozygous major allele genotypes
            int minGCnt = (alleles[1][1] - hetGCnt) / 2; // number of homozygous minor allele genotypes
            double propHets = (double) hetGCnt / (double) (hetGCnt + majGCnt + minGCnt);
            double expHets = 2.0 * theMAF * (1 - theMAF);
            obsF = 1.0 - (propHets / expHets);
            if (report) {
                System.out.printf("%d %d %d propHets:%g expHets:%g obsF:%g %n", majGCnt, minGCnt, hetGCnt, propHets, expHets, obsF);
            }
            return obsF;
        }
    }

    private byte[] getGoodAlleles(byte[] calls, int[][] allelesByFreq, byte homMaj, byte homMin, byte hetG1, byte hetG2, byte majAllele, byte minAllele) {
        if (includeRareAlleles()) {
            byte[] byteAlleles = new byte[allelesByFreq[0].length];
            for (int a = 0; a < allelesByFreq[0].length; a++) {
                byteAlleles[a] = (byte) allelesByFreq[0][a];
            }
            return byteAlleles;
        } else {
            setBadCallsToMissing(calls, homMaj, homMin, hetG1, hetG2, majAllele, minAllele);
            allelesByFreq = GenotypeTableUtils.getAllelesSortedByFrequency(calls); // the allele frequency & number of alleles may have been altered by setBadCallsToMissing()
            if (allelesByFreq[0].length < 2) {
                return null; // setBadCallsToMissing() rendered the site invariant
            } else if (allelesByFreq[0].length == 2) {
                return getMajMinAllelesOnly(allelesByFreq);
            } else {
                if (callBiallelicSNPsWithGap()) {
                    boolean hasGap = false;
                    for (int a = 2; a < allelesByFreq[0].length; a++) {  // NOTE: the maj & min alleles are not checked (we know that they are NOT gap (inclGaps mutually exclusive with callBiallelicSNPsWithGap)
                        if (((byte) allelesByFreq[0][a]) == NucleotideAlignmentConstants.GAP_ALLELE) {
                            hasGap = true;
                            break;
                        }
                    }
                    if (hasGap) {
                        byte[] byteAlleles = new byte[3];
                        byteAlleles[0] = (byte) allelesByFreq[0][0];
                        byteAlleles[1] = (byte) allelesByFreq[0][1];
                        byteAlleles[2] = NucleotideAlignmentConstants.GAP_ALLELE;
                        return byteAlleles;
                    } else {
                        return getMajMinAllelesOnly(allelesByFreq);
                    }
                } else {
                    return getMajMinAllelesOnly(allelesByFreq);
                }
            }
        }
    }

    private byte[] getMajMinAllelesOnly(int[][] alleles) {
        byte[] byteAlleles = new byte[2];
        byteAlleles[0] = (byte) alleles[0][0];
        byteAlleles[1] = (byte) alleles[0][1];
        return byteAlleles;
    }

    private void setBadCallsToMissing(byte[] calls, byte homMaj, byte homMin, byte hetG1, byte hetG2, byte majAllele, byte minAllele) {
        if (callBiallelicSNPsWithGap()) {
            for (int i = 0; i < calls.length; i++) {
                if (isGoodBiallelicWithGapCall(calls[i], homMaj, homMin, hetG1, hetG2, majAllele, minAllele)) {
                    continue;
                } else {
                    calls[i] = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
                }
            }
        } else {
            for (int i = 0; i < calls.length; i++) {
                if ((calls[i] == homMaj) || (calls[i] == homMin) || (calls[i] == hetG1) || (calls[i] == hetG2)) {
                    continue;
                } else {
                    calls[i] = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
                }
            }
        }
    }

    private boolean isGoodBiallelicWithGapCall(byte call, byte homMaj, byte homMin, byte hetG1, byte hetG2, byte majAllele, byte minAllele) {
        if (call == homMaj) {
            return true;
        }
        if (call == homMin) {
            return true;
        }
        if (call == hetG1) {
            return true;
        }
        if (call == hetG2) {
            return true;
        }
        if (call == GenotypeTableUtils.getDiploidValue(majAllele, NucleotideAlignmentConstants.GAP_ALLELE)) {
            return true;
        }
        if (call == GenotypeTableUtils.getDiploidValue(NucleotideAlignmentConstants.GAP_ALLELE, majAllele)) {
            return true;
        }
        if (call == GenotypeTableUtils.getDiploidValue(minAllele, NucleotideAlignmentConstants.GAP_ALLELE)) {
            return true;
        }
        if (call == GenotypeTableUtils.getDiploidValue(NucleotideAlignmentConstants.GAP_ALLELE, minAllele)) {
            return true;
        }
        if (call == NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) {
            return true;
        }
        return false;
    }

    private DataOutputStream openLocusLog(String logFileName) {
        try {
            DataOutputStream locusLogDOS
                    = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(new File(logFileName)), 65536));
            locusLogDOS.writeBytes(
                    "chr\tstart\tend\tstrand\ttotalbp\tnTags\tnReads\tnTaxaCovered\tminTaxaCovered\tstatus\tnVariableSites\tposVariableSites\tnVarSitesKept\tposVarSitesKept\trefTag?\tmaxTagLen\tminTagLen\n");
            return locusLogDOS;
        } catch (Exception e) {
            catchLocusLogException(e);
        }
        return null;
    }

    private void logRejectedTagLocus(TagLocus currTAL, DataOutputStream locusLogDOS) {
        int start, end;
        if (currTAL.getStrand() == -1) {
            end = currTAL.getMinStartPosition();
            start = currTAL.getMinStartPosition() - currTAL.getMaxTagLength() + 1;
        } else {
            start = currTAL.getMinStartPosition();
            end = currTAL.getMinStartPosition() + currTAL.getMaxTagLength() - 1;
        }
        int totalbp = end - start + 1;
        String status, refTag;
        if (currTAL.getSize() == 1) {
            status = "invariant\t0";
            refTag = currTAL.getDivergenceOfTag(0) == 0 ? "1" : "0";
        } else {
            status = "tooFewTaxa\tNA";
            boolean refTagFound = false;
            int t = -1;
            while (!refTagFound && t < currTAL.getSize() - 1) {
                t++;
                if (currTAL.getDivergenceOfTag(t) == 0) {
                    refTagFound = true;
                }
            }
            refTag = refTagFound ? "1" : "0";
        }
        try {
            locusLogDOS.writeBytes(
                    currTAL.getChromosome() + "\t"
                    + start + "\t"
                    + end + "\t"
                    + currTAL.getStrand() + "\t"
                    + totalbp + "\t"
                    + currTAL.getSize() + "\t"
                    + currTAL.getTotalNReads() + "\t"
                    + currTAL.getNumberTaxaCovered() + "\t"
                    + minTaxaWithLocus + "\t"
                    + status + "\t"
                    + "NA" + "\t"
                    + "0" + "\t"
                    + "NA" + "\t"
                    + refTag + "\t"
                    + currTAL.getMaxTagLength() + "\t"
                    + currTAL.getMinTagLength() + "\n"
            );
        } catch (Exception e) {
            catchLocusLogException(e);
        }
    }

    private void logAcceptedTagLocus(String locusLogRecord, DataOutputStream locusLogDOS) {
        try {
            locusLogDOS.writeBytes(locusLogRecord);
        } catch (Exception e) {
            catchLocusLogException(e);
        }
    }

    private void catchLocusLogException(Exception e) {
        System.out.println("ERROR: Unable to write to locus log file: " + e);
        e.printStackTrace();
        System.exit(1);
    }

    private byte[] filterCallsForInbreds(byte[] calls) {
        byte[] callsForInbredsOnly = new byte[nInbredTaxa];
        int inbred = 0;
        for (int taxon = 0; taxon < calls.length; taxon++) {
            if (useTaxaForMinF[taxon]) {
                callsForInbredsOnly[inbred] = calls[taxon];
                inbred++;
            }
        }
        return callsForInbredsOnly;
    }

    public static HashMap<String, Double> readTaxaFsFromFile(File pedigreeFile) {
        HashMap<String, Double> taxaFs = new HashMap<String, Double>();
        String inputLine = "Nothing has been read from the pedigree input file yet";
        int nameCol = -1, fCol = -1, nTaxa = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(pedigreeFile), 65536);
            inputLine = br.readLine();  // header line
            String[] cells = inputLine.split("\t");  // headers
            for (int col = 0; col < cells.length; col++) {
                if (cells[col].equalsIgnoreCase("Name")) {
                    nameCol = col;
                }
                if (cells[col].equalsIgnoreCase("F")) {
                    fCol = col;
                }
            }
            if (nameCol > -1 && fCol > -1) {
                while ((inputLine = br.readLine()) != null) {
                    cells = inputLine.split("\t");
                    if (cells[fCol].equals("NA")) {
                        taxaFs.put(cells[nameCol], -2.0);
                    } else {
                        taxaFs.put(cells[nameCol], Double.parseDouble(cells[fCol]));
                    }
                    ++nTaxa;
                }
            } else {
                throw new Exception("Name and/or F column not found in header");
            }
        } catch (Exception e) {
            myLogger.error("Catch in reading pedigree file e=" + e);
            e.printStackTrace();
            System.out.println(inputLine);
            return null;
        }
        myLogger.info(nTaxa + " taxa read from the pedigree file");
        return taxaFs;
    }

    private boolean maskNonInbredTaxa() {
        useTaxaForMinF = new boolean[theTBT.getTaxaCount()];  // initialized to false
        nInbredTaxa = 0;
        try {
            for (int taxon = 0; taxon < theTBT.getTaxaCount(); taxon++) {
                if (taxaFs.containsKey(theTBT.getTaxaName(taxon))) {
                    if (taxaFs.get(theTBT.getTaxaName(taxon)) >= minimumF()) {
                        useTaxaForMinF[taxon] = true;
                        nInbredTaxa++;
                    }
                } else {
                    throw new Exception("Taxon " + theTBT.getTaxaName(taxon) + " not found in the pedigree file");
                }
            }
            myLogger.info(nInbredTaxa + " taxa with an Expected F >= the mnF of " + minimumF() + " were found in the input TBT");
            return true;
        } catch (Exception e) {
            myLogger.error("Mismatch between TBT and pedigree file e=" + e);
            e.printStackTrace();
            return false;
        }
    }

    private long[] readReferenceGenomeChr(String inFileStr, int targetChr) {
        int nBases = getLengthOfReferenceGenomeChr(inFileStr, targetChr);
        if (nBases == 0) {
            return null;
        }
        int basesPerLong = BaseEncoder.chunkSize;
        int nLongs = (nBases % basesPerLong == 0) ? nBases / basesPerLong : (nBases / basesPerLong) + 1;
        long[] refGenomeChrAsLongs = new long[nLongs];
        myLogger.info("\n\nReading in the target chromosome " + targetChr + " from the reference genome fasta file: " + inFileStr);
        String temp = "Nothing has been read yet from the reference genome fasta file";
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(inFileStr)));
            StringBuilder currStrB = new StringBuilder();
            int currChr = Integer.MIN_VALUE, chunk = 0;
            while (br.ready()) {
                temp = br.readLine().trim();
                if (temp.startsWith(">")) {
                    if (chunk > 0) {
                        break;  // finished reading the targetChr (no need to read the rest of the file)
                    }
                    String chrS = temp.replace(">", "");
                    chrS = chrS.replace("chr", "");
                    currChr = Integer.parseInt(chrS);  // don't need to catch exception because getLengthOfReferenceGenomeChr() would have caught it already
                    myLogger.info("Currently reading chromosome " + currChr + " (target chromosome = " + targetChr + ")");
                } else if (currChr == targetChr) {
                    currStrB.append(temp.replace("N", "A")); // BaseEncoder encodes sequences with N's as (long) -1
                    while (currStrB.length() >= basesPerLong) {
                        refGenomeChrAsLongs[chunk] = BaseEncoder.getLongFromSeq(currStrB.substring(0, basesPerLong));
                        currStrB = (currStrB.length() > basesPerLong) ? new StringBuilder(currStrB.substring(basesPerLong)) : new StringBuilder();
                        chunk++;
                        if (chunk % 1000000 == 0) {
                            myLogger.info(chunk + " chunks of " + basesPerLong + " bases read from the reference genome fasta file for chromosome " + targetChr);
                        }
                    }
                }
            }
            if (currStrB.length() > 0) {
                refGenomeChrAsLongs[chunk] = BaseEncoder.getLongFromSeq(currStrB.toString());
                chunk++;
            }
            myLogger.info("\n\nFinished reading target chromosome " + targetChr + " into a total of " + chunk + " " + basesPerLong + "bp chunks\n\n");
            if (chunk != nLongs) {
                throw new Exception("The number of 32 base chunks read (" + chunk + ") was not equal to the expected number (" + nLongs + ")");
            }
            br.close();
        } catch (Exception e) {
            myLogger.error("Exception caught while reading the reference genome fasta file at line.  Error=" + e);
            e.printStackTrace();
            System.exit(1);
        }
        return refGenomeChrAsLongs;
    }

    private int getLengthOfReferenceGenomeChr(String inFileStr, int targetChr) {
        myLogger.info("\n\nDetermining the length (in bases) of target chromosome " + targetChr + " in the reference genome fasta file: " + inFileStr);
        String temp = "Nothing has been read yet from the reference genome fasta file";
        int line = 0, nBases = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(inFileStr)));
            int currChr = Integer.MIN_VALUE;
            while (br.ready()) {
                temp = br.readLine().trim();
                line++;
                if (line % 1000000 == 0) {
                    myLogger.info(line + " lines read from the reference genome fasta file");
                }
                if (temp.startsWith(">")) {
                    if (nBases > 0) {
                        break;  // finished reading the targetChr (no need to read the rest of the file)
                    }
                    String chrS = temp.replace(">", "");
                    chrS = chrS.replace("chr", "");
                    try {
                        currChr = Integer.parseInt(chrS);
                    } catch (NumberFormatException e) {
                        myLogger.error("\n\nTagsToSNPByAlignment detected a non-numeric chromosome name in the reference genome sequence fasta file: " + chrS
                                + "\n\nPlease change the FASTA headers in your reference genome sequence to integers "
                                + "(>1, >2, >3, etc.) OR to 'chr' followed by an integer (>chr1, >chr2, >chr3, etc.)\n\n");
                        System.exit(1);
                    }
                    myLogger.info("Currently reading chromosome " + currChr + " (target chromosome = " + targetChr + ")");
                } else if (currChr == targetChr) {
                    nBases += temp.length();
                }
            }
            if (nBases == 0) {
                throw new Exception("Target chromosome (" + targetChr + ") not found");
            }
            myLogger.info("The target chromosome " + targetChr + " is " + nBases + " bases long");
            br.close();
        } catch (Exception e) {
            if (nBases == 0) {
                myLogger.warn("Exception caught while reading the reference genome fasta file at line " + line + "\n   e=" + e + "\n   Skipping this chromosome...");
            } else {
                myLogger.error("Exception caught while reading the reference genome fasta file at line " + line + "\n   e=" + e);
                e.printStackTrace();
                System.exit(1);
            }
        }
        return nBases;
    }

    private String getRefSeqInRegion(TagLocus theTAL) {
        int basesPerLong = BaseEncoder.chunkSize;
        int refSeqStartPosition = theTAL.getMinStartPosition() - 128;
        int startIndex = Math.max((refSeqStartPosition / basesPerLong) - 1, 0);
        int refSeqEndPosition = theTAL.getMaxStartPosition() + 128;
        int endIndex = Math.min((refSeqEndPosition / basesPerLong) + 1, refGenomeChr.length - 1);
        StringBuilder sb = new StringBuilder();
        for (int i = startIndex; i <= endIndex; ++i) {
            sb.append(BaseEncoder.getSequenceFromLong(refGenomeChr[i]));
        }
        theTAL.setMinStartPosition(startIndex * basesPerLong + 1);
        return sb.toString();
    }

    private void addRefTag(TagLocus theTAL) {
        String refTag;
        int basesPerLong = BaseEncoder.chunkSize;
        int refSeqStartPos, refSeqEndPos;
        if (theTAL.getStrand() == -1) {
            refSeqEndPos = theTAL.getMinStartPosition();
            refSeqStartPos = refSeqEndPos - theTAL.getMaxTagLength() + 1;
        } else {
            refSeqStartPos = theTAL.getMinStartPosition();
            refSeqEndPos = refSeqStartPos + theTAL.getMaxTagLength() - 1;
        }
        int startIndex = Math.max((refSeqStartPos / basesPerLong) - 1, 0);
        int endIndex = Math.min((refSeqEndPos / basesPerLong), refGenomeChr.length - 1);
        StringBuilder sb = new StringBuilder();
        for (int i = startIndex; i <= endIndex; ++i) {
            sb.append(BaseEncoder.getSequenceFromLong(refGenomeChr[i]));
        }
        refTag = sb.substring(Math.max(refSeqStartPos - startIndex * basesPerLong - 1, 0),
                Math.min(refSeqStartPos - startIndex * basesPerLong - 1 + theTAL.getMaxTagLength(), sb.length()));
        if (theTAL.getStrand() == -1) {
            refTag = revComplement(refTag);
        }
        theTAL.addRefTag(refTag, theTOPM.getTagSizeInLong(), theTOPM.getNullTag());
    }

    //TODO remove this.  It should be in some other class, like BaseEncoder
    public static byte getIUPACAllele(byte allele) {
        byte iupacAllele = (byte) 'N';
        switch (allele) {
            case 0x00:
                iupacAllele = (byte) 'A';
                break;
            case 0x01:
                iupacAllele = (byte) 'C';
                break;
            case 0x02:
                iupacAllele = (byte) 'G';
                break;
            case 0x03:
                iupacAllele = (byte) 'T';
                break;
            case 0x05:
                iupacAllele = (byte) '-';
                break;
            default:
                iupacAllele = (byte) 'N';
                break;
        }
        return iupacAllele;
    }

    public static String revComplement(String seq) {
        StringBuilder sb = new StringBuilder();
        for (int i = seq.length() - 1; i >= 0; i--) {
            sb.append(NucleotideAlignmentConstants.getNucleotideDiploidIUPACComplement(seq.charAt(i)));
        }
        return sb.toString();
    }

}

class CustomSNPLog {

    private BufferedWriter myWriter;
    private final String HEADER
            = "Chr" + "\t"
            + "TagLocusStartPos" + "\t"
            + "TagLocusStrand" + "\t"
            + "SNPPosition" + "\t"
            + "Alleles" + "\t"
            + "nTagsAtLocus" + "\t"
            + "nReads" + "\t"
            + "nTaxa" + "\t"
            + "nTaxaCovered" + "\t"
            + "nInbreds" + "\t"
            + "nInbredsCovered" + "\t"
            + "nInbreds1Read" + "\t"
            + "nInbreds1ReadMaj" + "\t"
            + "nInbreds1ReadMin" + "\t"
            + "nInbredsGT1Read" + "\t"
            + "nInbredsGT1ReadHomoMaj" + "\t"
            + "nInbredsGT1ReadHomoMin" + "\t"
            + "nInbredHets" + "\t"
            + "inbredCoverage" + "\t"
            + "inbredHetScore" + "\t"
            + "nOutbreds" + "\t"
            + "nOutbredsCovered" + "\t"
            + "nOutbreds1Read" + "\t"
            + "nOutbreds1ReadMaj" + "\t"
            + "nOutbreds1ReadMin" + "\t"
            + "nOutbredsGT1Read" + "\t"
            + "nOutbredsGT1ReadHomoMaj" + "\t"
            + "nOutbredsGT1ReadHomoMin" + "\t"
            + "nOutbredHets" + "\t"
            + "passed?" + "\n";

    public CustomSNPLog(String locusLogFileName) {
        String locusLogFileDir;
        try {
            locusLogFileDir = (new File(locusLogFileName)).getCanonicalFile().getParent();
            String snpLogFileName = locusLogFileDir + File.separator + "CustomSNPLog.txt";
            myWriter = Utils.getBufferedWriter(snpLogFileName, false);
            myWriter.append(HEADER);
        } catch (IOException e) {
            System.out.println("\n\nERROR: problem creating custom SNP log file: " + e + "\n\n");
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void writeEntry(String entry) {
        try {
            myWriter.append(entry);
        } catch (IOException e) {
            System.out.println("\n\nERROR: problem writing to custom SNP log file: " + e + "\n\n");
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void close() {
        try {
            myWriter.close();
        } catch (Exception e) {
            // do nothing;
        }
    }
}

class CustomSNPLogRecord {

    private int chr;
    private int tagLocusStartPos;
    private byte tagLocusStrand;
    private int snpPosition;
    private byte majAllele;
    private byte minAllele;
    private String alleles;
    private int nTagsAtLocus;
    private int nReads;
    private int nTaxa;
    private int nTaxaCovered;
    private int nInbreds;
    private int nInbredsCovered;
    private int nInbreds1Read;
    private int nInbreds1ReadMaj;
    private int nInbreds1ReadMin;
    private int nInbredsGT1Read;
    private int nInbredsGT1ReadHomoMaj;
    private int nInbredsGT1ReadHomoMin;
    private int nInbredHets;
    private int nOutbreds;
    private int nOutbredsCovered;
    private int nOutbreds1Read;
    private int nOutbreds1ReadMaj;
    private int nOutbreds1ReadMin;
    private int nOutbredsGT1Read;
    private int nOutbredsGT1ReadMajHomo;
    private int nOutbredsGT1ReadMinHomo;
    private int nOutbredHets;
    private double inbredCoverage;
    private double inbredHetScore;
    private boolean pass;
    private static final String DELIM = "\t";

    public CustomSNPLogRecord(int site, TagLocus myTAL, int position, boolean[] isInbred, boolean includeReference) {
        chr = myTAL.getChromosome();
        tagLocusStartPos = myTAL.getMinStartPosition();
        tagLocusStrand = myTAL.getStrand();
        snpPosition = position;
        byte[][] byteAlleles = myTAL.getCommonAlleles();
        majAllele = tagLocusStrand == -1 ? getNucleotideComplement(byteAlleles[0][site]) : byteAlleles[0][site];
        minAllele = tagLocusStrand == -1 ? getNucleotideComplement(byteAlleles[1][site]) : byteAlleles[1][site];
        alleles = NUCLEOTIDE_ALLELES[0][majAllele] + "/"
                + NUCLEOTIDE_ALLELES[0][minAllele];
        nTagsAtLocus = (includeReference) ? myTAL.getSize() - 1 : myTAL.getSize();
        nReads = myTAL.getTotalNReads();
        nTaxaCovered = myTAL.getNumberTaxaCovered();
        getCounts(site, myTAL.getAlleleDepthsInTaxa(), isInbred);
    }

    private void getCounts(int site, byte[][][] alleleDepthsInTaxa, boolean[] isInbred) {
        nTaxa = alleleDepthsInTaxa[0][site].length;
        int genoDepth, nAlleles;
        boolean majPresent;
        for (int tx = 0; tx < nTaxa; tx++) {
            genoDepth = 0;
            nAlleles = 0;
            majPresent = false;
            if (isInbred == null || isInbred[tx]) {  // if no pedigree file was used, assume that all taxa are inbred
                ++nInbreds;
                for (int a = 0; a < 2; a++) {
                    int alleleDepth = alleleDepthsInTaxa[a][site][tx];
                    if (alleleDepth > 0) {
                        genoDepth += alleleDepth;
                        ++nAlleles;
                        if (a == 0) {
                            majPresent = true;
                        }
                    }
                }
                if (nAlleles > 0) {
                    ++nInbredsCovered;
                    if (genoDepth > 1) {
                        ++nInbredsGT1Read;
                        if (nAlleles > 1) {
                            ++nInbredHets;
                        } else if (majPresent) {
                            ++nInbredsGT1ReadHomoMaj;
                        } else {
                            ++nInbredsGT1ReadHomoMin;
                        }
                    } else {
                        ++nInbreds1Read;
                        if (majPresent) {
                            ++nInbreds1ReadMaj;
                        } else {
                            ++nInbreds1ReadMin;
                        }
                    }
                }
            } else {
                ++nOutbreds;
                for (int a = 0; a < 2; a++) {
                    int alleleDepth = alleleDepthsInTaxa[a][site][tx];
                    if (alleleDepth > 0) {
                        genoDepth += alleleDepth;
                        ++nAlleles;
                        if (a == 0) {
                            majPresent = true;
                        }
                    }
                }
                if (nAlleles > 0) {
                    ++nOutbredsCovered;
                    if (genoDepth > 1) {
                        ++nOutbredsGT1Read;
                        if (nAlleles > 1) {
                            ++nOutbredHets;
                        } else if (majPresent) {
                            ++nOutbredsGT1ReadMajHomo;
                        } else {
                            ++nOutbredsGT1ReadMinHomo;
                        }
                    } else {
                        ++nOutbreds1Read;
                        if (majPresent) {
                            ++nOutbreds1ReadMaj;
                        } else {
                            ++nOutbreds1ReadMin;
                        }
                    }
                }
            }
        }
        inbredCoverage = (double) nInbredsCovered / nInbreds;
        inbredHetScore = (double) nInbredHets / (nInbredsGT1ReadHomoMin + nInbredHets + 0.5);
        if (inbredCoverage > 0.15 && inbredHetScore < 0.21) {
            pass = true;  // machine learning cutoffs set by Ed
        }
    }

    public boolean isGoodSNP() {
        return pass;
    }

    public double getInbredCoverage() {
        return inbredCoverage;
    }

    public double getInbredHetScore() {
        return inbredHetScore;
    }

    @Override
    public String toString() {
        StringBuilder sBuilder = new StringBuilder();
        sBuilder.append(String.valueOf(chr));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(tagLocusStartPos));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(tagLocusStrand));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(snpPosition));
        sBuilder.append(DELIM);
        sBuilder.append(alleles);
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nTagsAtLocus)).append(DELIM);
        sBuilder.append(String.valueOf(nReads));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nTaxa));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nTaxaCovered));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nInbreds));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nInbredsCovered));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nInbreds1Read));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nInbreds1ReadMaj));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nInbreds1ReadMin));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nInbredsGT1Read));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nInbredsGT1ReadHomoMaj));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nInbredsGT1ReadHomoMin));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nInbredHets));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(inbredCoverage));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(inbredHetScore));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nOutbreds));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nOutbredsCovered));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nOutbreds1Read));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nOutbreds1ReadMaj));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nOutbreds1ReadMin));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nOutbredsGT1Read));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nOutbredsGT1ReadMajHomo));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nOutbredsGT1ReadMinHomo));
        sBuilder.append(DELIM);
        sBuilder.append(String.valueOf(nOutbredHets));
        sBuilder.append(DELIM);
        if (pass) {
            sBuilder.append(String.valueOf(1));
        } else {
            sBuilder.append(String.valueOf(0));
        }
        sBuilder.append("\n");
        return sBuilder.toString();
    }
}

class TagLocusSiteQualityScores {

    private int[] siteIndicesInTAL;
    private double[] inbredCoverage;
    private double[] inbredHetScore;
    private byte[][] alleles;  // [nSites][nAlleles]
    private int[] position;
    private int currSize;

    public TagLocusSiteQualityScores(int nSites) {
        siteIndicesInTAL = new int[nSites];
        inbredCoverage = new double[nSites];
        inbredHetScore = new double[nSites];
        alleles = new byte[nSites][];
        position = new int[nSites];
        currSize = 0;
    }

    public void addSite(int siteIndex, double inbredCov, double inbredHetS, byte[] alleles, int position) {
        siteIndicesInTAL[currSize] = siteIndex;
        inbredCoverage[currSize] = inbredCov;
        inbredHetScore[currSize] = inbredHetS;
        this.alleles[currSize] = alleles;
        this.position[currSize] = position;
        ++currSize;
    }

    public int getSize() {
        return currSize;
    }

    public int getSiteInTAL(int site) {
        return siteIndicesInTAL[site];
    }

    public byte[] getAlleles(int site) {
        return alleles[site];
    }

    public int getPosition(int site) {
        return position[site];
    }

    public void sortByQuality() {

        Swapper swapperQual = new Swapper() {
            public void swap(int a, int b) {
                int tempInt;
                tempInt = siteIndicesInTAL[a];
                siteIndicesInTAL[a] = siteIndicesInTAL[b];
                siteIndicesInTAL[b] = tempInt;

                double score = inbredCoverage[a];
                inbredCoverage[a] = inbredCoverage[b];
                inbredCoverage[b] = score;

                score = inbredHetScore[a];
                inbredHetScore[a] = inbredHetScore[b];
                inbredHetScore[b] = score;

                byte[] tempAlleles = alleles[a];
                alleles[a] = alleles[b];
                alleles[b] = tempAlleles;

                tempInt = position[a];
                position[a] = position[b];
                position[b] = tempInt;
            }
        };

        IntComparator compQual = new IntComparator() {
            public int compare(int a, int b) {

                // reverse sort (high inbredCoverage is good)
                if (inbredCoverage[a] > inbredCoverage[b]) {
                    return -1;
                }
                if (inbredCoverage[a] < inbredCoverage[b]) {
                    return 1;
                }

                // normal sort (low inbredHetScore is good)
                if (inbredHetScore[a] < inbredHetScore[b]) {
                    return -1;
                }
                if (inbredHetScore[a] > inbredHetScore[b]) {
                    return 1;
                }

                // normal sort (low site indices are better because closer to start of read)
                if (siteIndicesInTAL[a] < siteIndicesInTAL[b]) {
                    return -1;
                }
                if (siteIndicesInTAL[a] > siteIndicesInTAL[b]) {
                    return 1;
                }
                return 0;
            }
        };

        GenericSorting.quickSort(0, currSize, compQual, swapperQual);
    }

}
