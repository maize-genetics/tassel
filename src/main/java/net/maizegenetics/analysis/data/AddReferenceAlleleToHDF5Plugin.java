/*
 * AddReferenceAlleleToHDF5Plugin
 */
package net.maizegenetics.analysis.data;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;

import java.util.List;

import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

/**
 *
 * @author Jeff Glaubitz (jcg233@cornell.edu)
 */
public class AddReferenceAlleleToHDF5Plugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(AddReferenceAlleleToHDF5Plugin.class);

    private PluginParameter<String> myInputGenotypes = new PluginParameter.Builder<>("i", null, String.class)
        .guiName("Input HDF5 Genotype File")
        .required(true)
        .inFile()
        .description("Input HDF5 genotype (*.h5) file to be annotated with the reference allele")
        .build();
    private PluginParameter<String> myRefGenome = new PluginParameter.Builder<>("ref", null, String.class)
        .guiName("Reference Genome File")
        .required(true)
        .inFile()
        .description("Reference genome file in fasta format")
        .build();
    private PluginParameter<String> myRefGenomeVersion = new PluginParameter.Builder<>("ver", null, String.class)
        .guiName("Reference Genome Version")
        .required(true)
        .description("Version of the reference genome")
        .build();
    private PluginParameter<String> myOutputGenotypes = new PluginParameter.Builder<>("o", null, String.class)
        .guiName("Output HDF5 Genotype File")
        .required(false)
        .outFile()
        .description("Output HDF5 genotype file annotated with the reference allele (Default: write to same folder as input, with '*.h5' replaced '*_withRef.h5')")
        .build();
    
    private BufferedReader refReader = null;
    private PositionListBuilder newPosListBuilder = null;
    private Chromosome currChr = null;
    private int currPos = Integer.MIN_VALUE;
    private String contextSeq;
    private final boolean writePositions = false;

    public AddReferenceAlleleToHDF5Plugin() {
        super(null, false);
    }

    public AddReferenceAlleleToHDF5Plugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public void postProcessParameters() {
        if (myOutputGenotypes.isEmpty()) {
            String outFile = (new File( inputHDF5GenotypeFile() )).getAbsolutePath().replaceFirst("\\.h5$", "_withRef.h5");
            outputHDF5GenotypeFile(outFile);
        }
        refReader = Utils.getBufferedReader(referenceGenomeFile());
    }


    @Override
    public DataSet processData(DataSet input) {
        String message = addRefAlleleToHDF5GenoTable();
        if(message != null) {
            myLogger.error(message);
            try {Thread.sleep(500);} catch(Exception e) {}
            throw new IllegalStateException(message);
        }
        try {refReader.close();} catch (IOException e) {}
        fireProgress(100);
        return null;
    }

    private String addRefAlleleToHDF5GenoTable() {
        String message = populatePositionsWithRefAllele();
        if (message != null) return message;
        PositionList newPos = newPosListBuilder.build();
        String genomeVer = newPos.hasReference() ? newPos.genomeVersion() : "unknown";
        myLogger.info("\nGenome version: "+genomeVer+"\n");
        GenotypeTableBuilder newGenos = GenotypeTableBuilder.getTaxaIncremental(newPos, outputHDF5GenotypeFile());
        IHDF5Reader h5Reader = HDF5Factory.open(inputHDF5GenotypeFile());
        List<String> taxaNames = HDF5Utils.getAllTaxaNames(h5Reader);
        int nTaxaWritten = 0;
        for (String taxonName : taxaNames) {
            Taxon taxon = HDF5Utils.getTaxon(h5Reader, taxonName);
            byte[] genos = HDF5Utils.getHDF5GenotypesCalls(h5Reader, taxonName);
            byte[][] depth = HDF5Utils.getHDF5GenotypesDepth(h5Reader, taxonName);
            newGenos.addTaxon(taxon, genos, depth);
            ++nTaxaWritten;
            if(nTaxaWritten % 100 == 0) {
                myLogger.info("...finished writing genotypes and depth for "+nTaxaWritten+" taxa ");
            }
        }
        newGenos.build();
        myLogger.info("\n\nFinished adding reference alleles to file:");
        myLogger.info("  "+outputHDF5GenotypeFile()+"\n\n");
        return null;
    }
    
    private String populatePositionsWithRefAllele() {
        IHDF5Writer h5Writer = HDF5Factory.open(inputHDF5GenotypeFile());
        PositionList oldPosList = PositionListBuilder.getInstance(h5Writer);
//        h5Writer.close();
        newPosListBuilder = new PositionListBuilder();
        if (writePositions) {
            String genomeVer = oldPosList.hasReference() ? oldPosList.genomeVersion() : "unkown";
            myLogger.info("\ncurrent genome version: "+genomeVer+"\n");
            myLogger.info("SNPID\tchr\tpos\tstr\tmaj\tmin\tref\tmaf\tcov\tcontext");
        }
        for (Position oldPos : oldPosList) {
            Chromosome chr = oldPos.getChromosome();
            int pos = oldPos.getPosition();
            byte strand = oldPos.getStrand();
            byte refAllele = retrieveRefAllele(chr, pos, strand);
            if (refAllele == Byte.MIN_VALUE) {
                return "\nCould not find position "+pos+" on chromosome "+chr+" in the reference genome fasta file.\n\n\n";
            }
            Position newPos = new GeneralPosition.Builder(oldPos) 
                .allele(WHICH_ALLELE.Reference, refAllele)
                .build();
            if (writePositions) writePosition(newPos, contextSeq);
            newPosListBuilder.add(newPos);
        }
        newPosListBuilder.genomeVersion(referenceGenomeVersion());
        myLogger.info("Finished populating positions with RefAllele");
        return null;
    }
    
    private byte retrieveRefAllele(Chromosome chr, int pos, int strand) {
        findChrInRefGenomeFile(chr);
        char currChar = findPositionInRefGenomeFile(pos);
        if (currPos == pos) {
            byte refAllele = NucleotideAlignmentConstants.getNucleotideAlleleByte(currChar);
            if (strand == -1) refAllele = NucleotideAlignmentConstants.getNucleotideComplement(refAllele);
            return refAllele;
        } else {
            myLogger.warn("currPos:"+currPos);
            return Byte.MIN_VALUE;
        }
    }
    
    private void findChrInRefGenomeFile(Chromosome chr) {
        String temp = "Nothing has been read from the reference genome fasta file yet";
        try {
            while (refReader.ready() && (currChr == null || currChr.compareTo(chr) < 0)){
                temp = refReader.readLine().trim();
                if (temp.startsWith(">")) {
                    String chrS = temp.replace(">", "");
                    currChr = new Chromosome(chrS); // Chromosome class removes leading "chr" or "chromosome"
                    if (!writePositions) myLogger.info("\nCurrently reading chromosome "+currChr.getName()+" from reference genome fasta file\n\n");
                }
                currPos = 0;
            }
        } catch (IOException e) {
            myLogger.error("Exception caught while reading the reference genome fasta file:\n  "+e+"\nLast line read:\n  "+temp);
            try {Thread.sleep(500);} catch(InterruptedException iE) {}
            e.printStackTrace();
            throw new IllegalStateException("Problem reading reference genome file");
        }
        if (!currChr.equals(chr)) {
            myLogger.error("\nCould not find chromosome "+chr+" in the reference genome fasta file.\nMake sure that the chromosomes are in numerical order in that file\n\n\n");
            try {Thread.sleep(500);} catch(InterruptedException iE) {}
            throw new IllegalStateException("Problem reading reference genome file: Make sure that the chromosomes are in numerical order");
        }
    }
    
    private char findPositionInRefGenomeFile(int pos) {
        char currChar = Character.MAX_VALUE;
        contextSeq = "";
        try {
            while (currPos < pos) {
                int intChar = refReader.read();
                if (intChar == -1) {
                    currPos = Integer.MAX_VALUE;
                    return Character.MAX_VALUE;
                }
                currChar = (char) intChar;
                if (!Character.isWhitespace(currChar)) {
                    currPos++;
                    if (pos - currPos < 60) {
                        contextSeq += currChar;
                    }
                }
            }
        } catch (IOException e) {
            myLogger.error("\n\nError reading reference genome file:\n  "+e+"\n\n");
            throw new IllegalStateException("Problem reading reference genome file");
        }
        return currChar;
    }
    
    private void writePosition(Position pos, String contextSeq) {
        myLogger.info(
            pos.getSNPID()+
            "\t"+pos.getChromosome().getChromosomeNumber()+
            "\t"+pos.getPosition()+
            "\t"+pos.getStrand()+
            "\t"+pos.getAllele(WHICH_ALLELE.GlobalMajor)+
            "\t"+pos.getAllele(WHICH_ALLELE.GlobalMinor)+
            "\t"+pos.getAllele(WHICH_ALLELE.Reference)+
            "\t"+pos.getGlobalMAF()+
            "\t"+pos.getGlobalSiteCoverage()+
            "\t"+contextSeq
        );
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Add reference allele";
    }

    @Override
    public String getToolTipText() {
        return "Add reference allele to HDF5 genotypes";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(AddReferenceAlleleToHDF5Plugin.class);
    // }

    /**
     * Input HDF5 genotype (*.h5) file to be annotated with
     * the reference allele
     *
     * @return Input HDF5 Genotype File
     */
    public String inputHDF5GenotypeFile() {
        return myInputGenotypes.value();
    }

    /**
     * Set Input HDF5 Genotype File. Input HDF5 genotype (*.h5)
     * file to be annotated with the reference allele
     *
     * @param value Input HDF5 Genotype File
     *
     * @return this plugin
     */
    public AddReferenceAlleleToHDF5Plugin inputHDF5GenotypeFile(String value) {
        myInputGenotypes = new PluginParameter<>(myInputGenotypes, value);
        return this;
    }

    /**
     * Reference genome file in fasta format
     *
     * @return Reference Genome File
     */
    public String referenceGenomeFile() {
        return myRefGenome.value();
    }

    /**
     * Set Reference Genome File. Reference genome file in
     * fasta format
     *
     * @param value Reference Genome File
     *
     * @return this plugin
     */
    public AddReferenceAlleleToHDF5Plugin referenceGenomeFile(String value) {
        myRefGenome = new PluginParameter<>(myRefGenome, value);
        return this;
    }

    /**
     * Version of the reference genome
     *
     * @return Reference Genome Version
     */
    public String referenceGenomeVersion() {
        return myRefGenomeVersion.value();
    }

    /**
     * Set Reference Genome Version. Version of the reference
     * genome
     *
     * @param value Reference Genome Version
     *
     * @return this plugin
     */
    public AddReferenceAlleleToHDF5Plugin referenceGenomeVersion(String value) {
        myRefGenomeVersion = new PluginParameter<>(myRefGenomeVersion, value);
        return this;
    }

    /**
     * Output HDF5 genotype file annotated with the reference
     * allele (Default: write to same folder as input, with
     * '*.h5' replaced '*_withRef.h5')
     *
     * @return Output HDF5 Genotype File
     */
    public String outputHDF5GenotypeFile() {
        return myOutputGenotypes.value();
    }

    /**
     * Set Output HDF5 Genotype File. Output HDF5 genotype
     * file annotated with the reference allele (Default:
     * write to same folder as input, with '*.h5' replaced
     * '*_withRef.h5')
     *
     * @param value Output HDF5 Genotype File
     *
     * @return this plugin
     */
    public AddReferenceAlleleToHDF5Plugin outputHDF5GenotypeFile(String value) {
        myOutputGenotypes = new PluginParameter<>(myOutputGenotypes, value);
        return this;
    }
}
