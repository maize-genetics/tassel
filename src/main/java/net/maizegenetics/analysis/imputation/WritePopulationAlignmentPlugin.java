package net.maizegenetics.analysis.imputation;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import javax.swing.ImageIcon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.analysis.imputation.NucleotideImputationUtils;
import net.maizegenetics.analysis.imputation.PopulationData;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.CombineGenotypeTable;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

public class WritePopulationAlignmentPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(WritePopulationAlignmentPlugin.class);
    public static final String brkptComment1 = "#Donor Haplotypes\n";
    public static final String brkptComment2 = "#Taxa Breakpoints\n";
    public static final String brkptComment3 = "#Block are defined chr:startPos:endPos:donor1:donor2 (-1 means no hypothesis)\n";
    private static final byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
    private static final byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
    private static final byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
    private static final byte AC = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
    private static final byte CA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
    boolean mergeAlignments = false;
    boolean writeParentCalls = true;
    boolean writeNucleotides = true;
    boolean outputDiploid = false;
    boolean writeBreakpoints = false;
    double minSnpCoverage = Double.NaN;//0.1;
    double maxMafForMono = Double.NaN;//0.01;
    boolean outputAlternateNucleotides = true;
    boolean writeToFile = false;
    String baseFile = "";
    String breakpointBase = "";
    String algorithm = "";
    Boolean breakpointHetsToMissing = false;
    
    public WritePopulationAlignmentPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        List<Datum> theData = input.getDataOfType(PopulationData.class);
        List<Datum> theResult = new ArrayList<Datum>();
        if (theData.size() > 0) {
            if (writeParentCalls) {
                theResult.addAll(writeOutput(theData, false));
            }
            if (writeNucleotides) {
                theResult.addAll(writeOutput(theData, true));
            }
            if (writeBreakpoints && breakpointBase.length() > 1) {
                for (Datum familyDatum : theData) {
                    PopulationData popdata = (PopulationData) familyDatum.getData();
                    writeBreakpoints(popdata);
                }
            }
            DataSet resultDataSet = new DataSet(theResult, this);
            fireDataSetReturned(new PluginEvent(resultDataSet, this));
            return resultDataSet;
        } else {
            return null;
        }
    }

    private List<Datum> writeOutput(List<Datum> theData, boolean asNucleotides) {
        List<Datum> theResult = new ArrayList<Datum>();
        if (mergeAlignments) {
            GenotypeTable[] allOfTheAlignments = new GenotypeTable[theData.size()];
            int count = 0;
            for (Datum datum : theData) {
                PopulationData family = (PopulationData) datum.getData();
                if (asNucleotides) {
                    allOfTheAlignments[count++] =
                            createOutputAlignmentImputingAllNucleotides(family);
                } else {
                    allOfTheAlignments[count++] = createOutputAlignment(family, asNucleotides);
                }
            }
            GenotypeTable myImputedGenotypes =
                    CombineGenotypeTable.getInstance(allOfTheAlignments, true);

            String myDatumName;
            String myDatumComment;
            String filepath;

            if (asNucleotides) {
                myDatumName = "imputed_genotypes";
                myDatumComment = "imputed genotypes, merged";
                filepath = baseFile + ".nuc.hmp.txt.gz";
            } else {
                myDatumName = "imputed_parents";
                myDatumComment =
                        "Imputed parents were coded as A and C.\nA and C were assigned at random for each chromosome independently.";
                filepath = baseFile + ".parents.hmp.txt.gz";
            }

            Datum myDatum = new Datum(myDatumName, myImputedGenotypes, myDatumComment);
            theResult.add(myDatum);
            if (writeToFile)
                ExportUtils.writeToHapmap(myImputedGenotypes, filepath);
        } else {
            for (Datum datum : theData) {
                PopulationData family = (PopulationData) datum.getData();
                String familyName = family.name.replace('/', '.');
                String chrName = family.original.chromosomeName(0);
                GenotypeTable myImputedGenotypes;
                StringBuilder myDatumName;
                StringBuilder myDatumComment;
                StringBuilder filepath = new StringBuilder(baseFile);
                if (asNucleotides) {
                    myImputedGenotypes = createOutputAlignmentImputingAllNucleotides(family);
                    myDatumName = new StringBuilder("imputed_genotypes_Chr");
                    myDatumName.append(chrName).append("_").append(familyName);
                    myDatumComment = new StringBuilder("imputed genotypes");
                    myDatumComment.append("\nchromosome ").append(chrName);
                    myDatumComment.append("\nfamily = ").append(familyName);
                    filepath.append(".chr").append(chrName).append(".").append(familyName).append(".nuc.hmp.txt.gz");
                } else {
                    myImputedGenotypes = createOutputAlignment(family, asNucleotides);
                    myDatumName = new StringBuilder("imputed_parents_Chr");
                    myDatumName.append(chrName).append("_").append(familyName);
                    myDatumComment = new StringBuilder("imputed parents");
                    myDatumComment.append("\nchromosome ").append(chrName);
                    myDatumComment.append("\nfamily = ").append(familyName);
                    myDatumComment.append("\nImputed parents have been coded as A and C.");
                    myDatumComment.append("\nA and C were assigned to parents at random for each chromosome independently.");
                    filepath.append(".chr").append(chrName).append(".").append(familyName).append(".parents.hmp.txt.gz");
                }

                Datum myDatum =
                        new Datum(myDatumName.toString(), myImputedGenotypes, myDatumComment.toString());
                theResult.add(myDatum);
                if (writeToFile)
                    ExportUtils.writeToHapmap(myImputedGenotypes, filepath.toString());
            }
        }
        return theResult;
    }

    private GenotypeTable createOutputAlignment(PopulationData popdata, boolean asNucleotides) {
        GenotypeTable out = null;

        if (!asNucleotides) {
            out = popdata.imputed;
        } else {
            //change the parent calls to original nucleotides
            GenotypeTable outPoly =
                    NucleotideImputationUtils.convertParentCallsToNucleotides(popdata);

            if (!Double.isNaN(minSnpCoverage) && !Double.isNaN(maxMafForMono)) {
                int nsnps = popdata.original.numberOfSites();
                double ngametes = 2 * popdata.original.numberOfTaxa();

                int[] monomorphicSnps = new int[nsnps];
                int snpCount = 0;
                for (int s = 0; s < nsnps; s++) {
                    double coverage = popdata.original.totalGametesNonMissingForSite(s) / ngametes;
                    if (!popdata.snpIndex.fastGet(s)
                            && popdata.original.minorAlleleFrequency(s) <= maxMafForMono
                            && coverage >= minSnpCoverage) {
                        monomorphicSnps[snpCount++] = s;
                    }
                }
                monomorphicSnps = Arrays.copyOf(monomorphicSnps, snpCount);
                GenotypeTable fa =
                        FilterGenotypeTable.getInstance(popdata.original, monomorphicSnps);
                if (fa.numberOfSites() == 0) {	//If there are no monomorphic sites (e.g, have been pre-filtered), just return polymorphic ones
                    out = outPoly;
                } else { //Return both monomorphic and polymorphic sites
                    GenotypeTableBuilder builder =
                            GenotypeTableBuilder.getSiteIncremental(fa.taxa());

                    // fill in all values with the major allele
                    nsnps = fa.numberOfSites();
                    int ntaxa = fa.numberOfTaxa();
                    for (int s = 0; s < nsnps; s++) {
                        byte majorAllele = fa.majorAllele(s);
                        byte major = (byte) ((majorAllele << 4) | majorAllele);
                        byte[] snpgeno = new byte[ntaxa];
                        Arrays.fill(snpgeno, major);
                        builder.addSite(fa.positions().get(s), snpgeno);
                    }
                    out = builder.build();
                }
            } else
                out = outPoly;
        }
        return out;
    }

    private GenotypeTable createOutputAlignmentImputingAllNucleotides(PopulationData family) {
        GenotypeTable filledImputedGenotypes = NucleotideImputationUtils.fillGapsInImputedAlignment(family);
        GenotypeTable filteredOriginalGenotypes = FilterGenotypeTable.getInstance(family.original, filledImputedGenotypes.taxa());
        int nsites = filteredOriginalGenotypes.numberOfSites();
        int nImputedSites = filledImputedGenotypes.numberOfSites();
        int[] imputedPos = filledImputedGenotypes.physicalPositions();
        int[] origPos = filteredOriginalGenotypes.physicalPositions();

        //first fill in gaps flanked by the same parent in the imputed alignment
        int ntaxa = filteredOriginalGenotypes.numberOfTaxa();
        GenotypeTableBuilder genoBuilder =
                GenotypeTableBuilder.getSiteIncremental(filteredOriginalGenotypes.taxa());
        for (int s = 0; s < nsites; s++) {

            byte[] nuc = filteredOriginalGenotypes.alleles(s);
            int nalleles = nuc.length;
            if (nalleles == 0) {
                //do nothing
                genoBuilder.addSite(filteredOriginalGenotypes.positions().get(s), filteredOriginalGenotypes.genotypeAllTaxa(s));
            } else if (nalleles > 0) {
                //find flanking markers in imputed
                int ndx = Arrays.binarySearch(imputedPos, origPos[s]);

                //do not impute if ndx is before the first or after the last marker
                if (ndx != -1 && ndx > -nImputedSites - 1) {
                    //deal with original sites that fall before the first imputed site or after the last one
                    //assume no recombination and assign them to the haplotype of the first or last site
                    if (ndx == -1)
                        ndx = 0;
                    if (-ndx - 1 >= nImputedSites)
                        ndx = nImputedSites - 1;

                    OpenBitSet mjImputed;
                    OpenBitSet mnImputed;

                    byte[] genotype = new byte[ntaxa];
                    if (ndx >= 0) {
                        mjImputed = new OpenBitSet(filledImputedGenotypes.allelePresenceForAllTaxa(ndx, WHICH_ALLELE.Major));
                        mnImputed = new OpenBitSet(filledImputedGenotypes.allelePresenceForAllTaxa(ndx, WHICH_ALLELE.Minor));
                    } else {
                        ndx = -ndx - 1; //the original site falls between ndx-1 and ndx
                        mjImputed = new OpenBitSet(filledImputedGenotypes.allelePresenceForAllTaxa(ndx, WHICH_ALLELE.Major));
                        mnImputed = new OpenBitSet(filledImputedGenotypes.allelePresenceForAllTaxa(ndx, WHICH_ALLELE.Minor));
                        OpenBitSet flankingSame = new OpenBitSet(filledImputedGenotypes.allelePresenceForAllTaxa(ndx, WHICH_ALLELE.Major));
                        flankingSame.notXor(filledImputedGenotypes.allelePresenceForAllTaxa(ndx - 1, WHICH_ALLELE.Major));
                        OpenBitSet mnImputedSame = new OpenBitSet(filledImputedGenotypes.allelePresenceForAllTaxa(ndx, WHICH_ALLELE.Minor));
                        mnImputedSame.notXor(filledImputedGenotypes.allelePresenceForAllTaxa(ndx - 1, WHICH_ALLELE.Minor));
                        flankingSame.and(mnImputedSame); //are the flanking imputed markers the same?

                        mjImputed.and(flankingSame);
                        mnImputed.and(flankingSame);
                    }

                    BitSet mjOrig = filteredOriginalGenotypes.allelePresenceForAllTaxa(s, WHICH_ALLELE.Major);
                    BitSet mnOrig = filteredOriginalGenotypes.allelePresenceForAllTaxa(s, WHICH_ALLELE.Minor);
                    OpenBitSet imj = new OpenBitSet(mjImputed);
                    imj.andNot(mnImputed); //homozygous major allele
                    OpenBitSet imn = new OpenBitSet(mnImputed);
                    imn.andNot(mjImputed); //homozygous minor allele
                    OpenBitSet omj = new OpenBitSet(mjOrig);
                    omj.andNot(mnOrig); //homozygous major allele
                    OpenBitSet omn = new OpenBitSet(mnOrig);
                    omn.andNot(mjOrig); //homozygous minor allele
                    int[][] counts = new int[2][2];
                    counts[0][0] = (int) OpenBitSet.intersectionCount(imj, omj);
                    counts[0][1] = (int) OpenBitSet.intersectionCount(imj, omn);
                    counts[1][0] = (int) OpenBitSet.intersectionCount(imn, omj);
                    counts[1][1] = (int) OpenBitSet.intersectionCount(imn, omn);
                    byte[] alleles = getMajorAndMinorAllelesAtSite(counts, nuc);
                    byte[] genotypes = filteredOriginalGenotypes.genotypeAllTaxa(s);
                    if (alleles != null) {
                        byte major = GenotypeTableUtils.getUnphasedDiploidValue(alleles[0], alleles[0]);
                        byte minor = GenotypeTableUtils.getUnphasedDiploidValue(alleles[1], alleles[1]);
                        byte het;
                        if (alleles[0] == GenotypeTable.UNKNOWN_ALLELE || alleles[1] == GenotypeTable.UNKNOWN_ALLELE)
                            het = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
                        else
                            het = GenotypeTableUtils.getUnphasedDiploidValue(alleles[0], alleles[1]);
                        for (int t = 0; t < ntaxa; t++) {

                            if (mjImputed.fastGet(t)) {
                                if (mnImputed.fastGet(t))
                                    genotypes[t] = het;
                                else
                                    genotypes[t] = major;
                            } else if (mnImputed.fastGet(t)) {
                                genotypes[t] = minor;
                            }

                        }
                    }
                    genoBuilder.addSite(filteredOriginalGenotypes.positions().get(s), genotypes);
                } else {
                    genoBuilder.addSite(filteredOriginalGenotypes.positions().get(s), filteredOriginalGenotypes.genotypeAllTaxa(s));
                }

            }

        }

        return genoBuilder.build();
    }

    private static byte[] getMajorAndMinorAllelesAtSite(int[][] counts, byte[] OriginalNucleotides) {
        int minratio = 4;
        int minPresent = 7;
        int totalPresent = counts[0][0] + counts[0][1] + counts[1][0] + counts[1][1];

        if ((counts[0][1] == 0 || counts[0][0] / counts[0][1] >= minratio)
                && (counts[1][0] == 0 || counts[1][1] / counts[1][0] >= minratio)) {
            if (OriginalNucleotides.length == 1) {
                if (totalPresent < minPresent)
                    return null; //not enough data to call
                return new byte[] { OriginalNucleotides[0], GenotypeTable.UNKNOWN_ALLELE }; //biologically missing
            }
            return OriginalNucleotides;
        }

        if ((counts[0][0] == 0 || counts[0][1] / counts[0][0] >= minratio)
                && (counts[1][1] == 0 || counts[1][0] / counts[1][1] >= minratio)) {
            if (OriginalNucleotides.length == 1) {
                if (totalPresent < minPresent)
                    return null;  //not enough data to call
                return new byte[] { GenotypeTable.UNKNOWN_ALLELE, OriginalNucleotides[0] }; //biologically missing
            }
            return new byte[] { OriginalNucleotides[1], OriginalNucleotides[0] };
        }

        //the site may be monomorphic
        //if so, one of the parents may be biologically missing
        //if a column sum is 1 or 0, then the site is monomorphic
        //if a row sum is 1 or 0, then that parent is biologically missing
        minratio = 10;
        int col0 = counts[0][0] + counts[1][0];
        int col1 = counts[0][1] + counts[1][1];
        int row0 = counts[0][0] + counts[0][1];
        int row1 = counts[1][0] + counts[1][1];

        if (row0 > 1 && row1 > 1) { 		//not biologically missing case
            if (col1 <= 1 || col0 / col1 >= minratio)
                return new byte[] { OriginalNucleotides[0], OriginalNucleotides[0] };
            if (col0 <= 1 || col1 / col0 >= minratio)
                return new byte[] { OriginalNucleotides[1], OriginalNucleotides[1] };
        } else {	//biologically missing case
            if (row0 <= 1 && totalPresent >= minPresent)
                return new byte[] { GenotypeTable.UNKNOWN_ALLELE, OriginalNucleotides[0] };
            if (row1 <= 1 && totalPresent >= minPresent)
                return new byte[] { OriginalNucleotides[0], GenotypeTable.UNKNOWN_ALLELE };
        }
        return null;
    }

    private void writeBreakpoints(PopulationData family) {
//        String filepath = String.format("%s_%s_%s_%s.pa.txt.gz", breakpointBase, family.chromosome, family.name, algorithm);
        String filepath = String.format("%s_%s_%s_%s.pa.txt", breakpointBase, family.chromosome, family.name, algorithm);
//        try (BufferedWriter bw =
//                new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(filepath.toString()))))) {
        
        //count the number of taxa to be written. does not include parents
        int nsites = family.imputed.numberOfSites();
        int ntaxa = family.imputed.numberOfTaxa();
        int nOutputTaxa = 0;
        for (int t = 0; t < ntaxa; t++) {
            String tname = family.imputed.taxaName(t);
            if (!tname.equals(family.parent1) && !tname.equals(family.parent2)) nOutputTaxa++;;

        }
        
        try (BufferedWriter bw =
                new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filepath.toString())))) {
            bw.write(String.format("%d\t%d\n", 2, nOutputTaxa));
            bw.write(brkptComment1);
            bw.write(String.format("0\t%s\n", family.parent1));
            bw.write(String.format("1\t%s\n", family.parent2));

            bw.write(brkptComment2);
            bw.write(brkptComment3);
            for (int t = 0; t < ntaxa; t++) {
                String tname = family.imputed.taxaName(t);
                if (tname.equals(family.parent1) || tname.equals(family.parent2)) continue;
                bw.write(tname);
                byte[] taxaGeno = family.imputed.genotypeAllSites(t);
                int segStartPos = family.imputed.chromosomalPosition(0);
                int segStart = 0;
                int segGeno = taxaGeno[0];
                for (int s = 0; s < nsites; s++) {
                    byte tGeno = taxaGeno[s];
                    if (tGeno == CA)
                        tGeno = AC;
                    if (tGeno != segGeno) {
                        int segEndPos = family.imputed.chromosomalPosition(s - 1);
                        if (segGeno == AA)
                            bw.write(String.format("\t%s:%d:%d:%d:%d", family.chromosome, segStartPos, segEndPos, 0, 0));
                        if (segGeno == CC)
                            bw.write(String.format("\t%s:%d:%d:%d:%d", family.chromosome, segStartPos, segEndPos, 1, 1));
                        if (segGeno == AC && !breakpointHetsToMissing)
                            bw.write(String.format("\t%s:%d:%d:%d:%d", family.chromosome, segStartPos, segEndPos, 0, 1));
                        segStart = s;
                        segGeno = tGeno;
                        segStartPos = family.imputed.chromosomalPosition(s);
                    }
                }
                int segEndPos = family.imputed.chromosomalPosition(nsites - 1);
                if (segGeno == AA)
                    bw.write(String.format("\t%s:%d:%d:%d:%d", family.chromosome, segStartPos, segEndPos, 0, 0));
                if (segGeno == CC)
                    bw.write(String.format("\t%s:%d:%d:%d:%d", family.chromosome, segStartPos, segEndPos, 1, 1));
                if (segGeno == AC && !breakpointHetsToMissing)
                    bw.write(String.format("\t%s:%d:%d:%d:%d", family.chromosome, segStartPos, segEndPos, 0, 1));

                bw.write("\n");
            }

        } catch (IOException e) {
            throw new RuntimeException("Error writing to file " + filepath, e);
        }
    }

    @Override
    public void setParameters(String[] args) {
        if (args == null || args.length == 0) {
            myLogger.info(getUsage());
            return;
        }

        int narg = args.length;
        for (int i = 0; i < narg; i++) {
            if (args[i].equals("-m") || args[i].equalsIgnoreCase("-merge")) {
                String val = args[++i];
                if (val.toUpperCase().startsWith("T")) {
                    mergeAlignments = true;
                } else {
                    mergeAlignments = false;
                }
            } else if (args[i].equals("-f") || args[i].equalsIgnoreCase("-file")) {
                setBaseFile(args[++i]);
            } else if (args[i].equals("-p") || args[i].equalsIgnoreCase("-parentCalls")) {
                String val = args[++i];
                if (val.toUpperCase().startsWith("T")) {
                    writeParentCalls = true;
                    writeNucleotides = false;
                } else {
                    writeParentCalls = false;
                    writeNucleotides = true;
                }
            } else if (args[i].equals("-o") || args[i].equalsIgnoreCase("-outputType")) {
                String val = args[++i];
                if (val.toUpperCase().startsWith("P")) {
                    writeParentCalls = true;
                    writeNucleotides = false;
                } else if (val.toUpperCase().startsWith("N")) {
                    writeParentCalls = false;
                    writeNucleotides = true;
                } else if (val.toUpperCase().startsWith("B")) {
                    writeParentCalls = true;
                    writeNucleotides = true;
                } else {
                    writeParentCalls = true;
                    writeNucleotides = false;
                }
            } else if (args[i].equals("-d") || args[i].equalsIgnoreCase("-diploid")) {
                String val = args[++i];
                if (val.toUpperCase().startsWith("T")) {
                    outputDiploid = true;
                } else {
                    outputDiploid = false;
                }
            } else if (args[i].equals("-c") || args[i].equalsIgnoreCase("-minCoverage")) {
                minSnpCoverage = Double.parseDouble(args[++i]);
            } else if (args[i].equals("-x") || args[i].equalsIgnoreCase("-maxMono")) {
                maxMafForMono = Double.parseDouble(args[++i]);
            } else if (args[i].equals("-pa") || args[i].equalsIgnoreCase("-breakpointFile")) {
                writeBreakpoints = true;
                breakpointBase = args[++i];
            } else if (args[i].equals("-al") || args[i].equalsIgnoreCase("-algorithm")) {
                algorithm = args[++i];
            } else if (args[i].equalsIgnoreCase("-bpHetsMissing")) {
                breakpointHetsToMissing = true;
            } else if (args[i].equals("?")) {
                myLogger.info(getUsage());
            }
        }
    }

    public void setMergeAlignments(boolean mergeAlignments) {
        this.mergeAlignments = mergeAlignments;
    }

    public void setWriteParentCalls(boolean writeParentCalls) {
        this.writeParentCalls = writeParentCalls;
    }

    public void setOutputDiploid(boolean outputDiploid) {
        this.outputDiploid = outputDiploid;
    }

    //    public void setBaseFileName(String baseFileName) {
    //        this.baseFileName = baseFileName;
    //    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Write Populations";
    }

    @Override
    public String getToolTipText() {
        return null;
    }

    public String getUsage() {
        StringBuilder usage =
                new StringBuilder("The WritePopulationAlignmentPlugin can take the following parameters:\n");
        usage.append("-f or -file : the base file name and path to which output will be written");
        usage.append("-m or -merge : if true families are merged into a data set, if false each family is output to a separate data set (default = false)\n");
        usage.append("-o or -outputType : parents = output parent calls, nucleotides = output nucleotides, both = output both (the default)\n");
        usage.append("-c or -minCoverage : the minimum coverage for a monomorphic snp to be included in the nucleotide output\n");
        usage.append("-x or -maxMono : the maximum minor allele frequency used to call monomorphic snps\n");
        usage.append("if -c or -x equals NaN and merge is true, then missing values at monomorphic sites (within a family) will be left missing\n");
        usage.append("-pa or -breakpointFile : the stem for the basepoint file names. Chromosome, family, algorithm and .pa.txt.gz will be appended.");
        usage.append("-al or -algorithm : the algorithm used for imputation. Only used to construct the full breakpoint file names.");
        usage.append("? : print the parameter list.\n");

        return usage.toString();
    }

    public void setWriteNucleotides(boolean writeNucleotides) {
        this.writeNucleotides = writeNucleotides;
    }

    public void setMinSnpCoverage(double minSnpCoverage) {
        this.minSnpCoverage = minSnpCoverage;
    }

    public void setMaxMafForMono(double maxMafForMono) {
        this.maxMafForMono = maxMafForMono;
    }

    public void setBaseFile(String baseFile) {
        this.baseFile = baseFile;
        writeToFile = true;
    }

    public void setBreakpointBaseName(String bkptBase) {
        breakpointBase = bkptBase;
        writeBreakpoints = true;
    }

    public void setAlgorithmName(String algorithm) {
        this.algorithm = algorithm;
    }
    
    public void setDeleteHetsInBreakpoints(boolean hetsToMissing) {
        breakpointHetsToMissing = hetsToMissing;
    }
}
