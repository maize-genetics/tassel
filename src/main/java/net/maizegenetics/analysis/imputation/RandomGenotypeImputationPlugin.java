package net.maizegenetics.analysis.imputation;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;

public class RandomGenotypeImputationPlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(RandomGenotypeImputationPlugin.class);
    private static final byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
    public static final String tab = "\t";
    private final Random randomizer = new Random();
    
    private PluginParameter<Boolean> writeFile =
            new PluginParameter.Builder<>("tofile", false, Boolean.class)
                    .description("Should the imputed data be written directly to a hapmap file (in IUPAC format) rather than stored in memory. (Default = false). "
                            + "Data should be filtered on proportion missing prior to imputation.")
                    .guiName("Write Hapmap File")
                    .build();

    private PluginParameter<String> outfile =
            new PluginParameter.Builder<>("filename", null, String.class)
                    .description("The name of the Hapmap file to be written.")
                    .guiName("Hapmap File Name")
                    .dependentOnParameter(writeFile)
                    .outFile()
                    .build();

    public RandomGenotypeImputationPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        List<Datum> genotypeList = input.getDataOfType(GenotypeTable.class);
        if (genotypeList.size() != 1) {
            String errmsg =
                    "Error in random imputation: exactly one input genotype data set must be selected.";
            if (isInteractive())
                JOptionPane.showMessageDialog(getParentFrame(), errmsg, "Error", JOptionPane.ERROR_MESSAGE);
            else
                myLogger.error(errmsg);
            return null;
        }
        GenotypeTable myGeno = (GenotypeTable) genotypeList.get(0).getData();
        String dataname = genotypeList.get(0).getName();
        
        if (writeFile.value()) {
            if (outfile.value() == null) throw new IllegalArgumentException("Must supply an output filename when tofile is true.");
            return writeImputationToHapmapFile(myGeno, outfile.value());
        } else {
            return storeImputationInMemory(myGeno, dataname);
        }
        
    }

    private DataSet storeImputationInMemory(GenotypeTable myGenotype, String dataname) {

        GenotypeTableBuilder myBuilder = GenotypeTableBuilder.getSiteIncremental(myGenotype.taxa());
        int nsites = myGenotype.numberOfSites();
        int progressAt = nsites / 50;
        progressAt = Math.max(1,  progressAt);
        
        for (int s = 0; s < nsites; s++) {
            byte[] sitegeno = myGenotype.genotypeAllTaxa(s);
            Optional<byte[]> imputedgeno = imputeRandomGenotypes(sitegeno);
            if (imputedgeno.isPresent()) {
                myBuilder.addSite(myGenotype.positions().get(s), sitegeno);
            }
            if (s % progressAt == 0) {
                fireProgress(s * 100 / nsites);
            }
        }

        String comment =
                "Missing genotypes imputed randomly\nImputed genotypes selected from genotype distribution for each site.";
        Datum outDatum = new Datum("Imputed_" + dataname, myBuilder.build(), comment);

        return new DataSet(outDatum, this);
    }

    private DataSet writeImputationToHapmapFile(GenotypeTable myGenotype, String filename) {
        
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(filename))) {
            bw.write(hapmapHeader(myGenotype.taxa()));
            int nsites = myGenotype.numberOfSites();
            int ntaxa = myGenotype.numberOfTaxa();
            int progressAt = nsites / 50;
            progressAt = Math.max(1,  progressAt);
            for (int s = 0; s < nsites; s++) {
                byte[] sitegeno = myGenotype.genotypeAllTaxa(s);
                Optional<byte[]> imputedgeno = imputeRandomGenotypes(sitegeno);
                if (imputedgeno.isPresent()) {
                    Position pos = myGenotype.positions().get(s);
                    StringBuilder sb = new StringBuilder(pos.getSNPID());
                    String alleles = String.format("%s/%s",
                            NucleotideAlignmentConstants.getHaplotypeNucleotide(pos.getAllele(WHICH_ALLELE.Major)),
                            NucleotideAlignmentConstants.getHaplotypeNucleotide(pos.getAllele(WHICH_ALLELE.Minor)));
                    sb.append(tab).append(alleles);
                    sb.append(tab).append(pos.getChromosome().getName());
                    sb.append(tab).append(pos.getPosition());
                    if (pos.getStrand() == Position.STRAND_PLUS) {
                        sb.append(tab).append("+");
                    } else if (pos.getStrand() == Position.STRAND_MINUS) {
                        sb.append(tab).append("-");
                    } else {
                        sb.append(tab).append("NA");
                    }
                    sb.append(tab).append("NA");
                    sb.append(tab).append("NA");
                    sb.append(tab).append("NA");
                    sb.append(tab).append("NA");
                    sb.append(tab).append("NA");
                    sb.append(tab).append("NA");
                    
                    for (int t = 0; t < ntaxa; t++) sb.append(tab).append(NucleotideAlignmentConstants.getNucleotideIUPAC(sitegeno[t]));
                    sb.append("\n");
                    bw.write(sb.toString());
                }
                if (s % progressAt == 0) {
                    fireProgress(s * 100 / nsites);
                }
            }

        } catch (IOException ioe) {
            throw new RuntimeException("Error writing to " + filename, ioe);
        }
        return null;
    }

    private String hapmapHeader(TaxaList taxa) {
        StringBuilder sb = new StringBuilder("rs#");
        sb.append(tab).append("alleles");
        sb.append(tab).append("chrom");
        sb.append(tab).append("pos");
        sb.append(tab).append("strand");
        sb.append(tab).append("assembly#");
        sb.append(tab).append("center");
        sb.append(tab).append("protLSID");
        sb.append(tab).append("assayLSID");
        sb.append(tab).append("panel");
        sb.append(tab).append("QCcode");
        for (Taxon taxon : taxa) {
            sb.append(tab).append(taxon.getName());
        }
        sb.append("\n");
        return sb.toString();
    }
    
    private Optional<byte[]> imputeRandomGenotypes(byte[] sitegeno) {
        Object[] genotypeCounts = byteCounts(sitegeno);
        byte[] siteGenotypes = (byte[]) genotypeCounts[0];
        int[] genoCounts = (int[]) genotypeCounts[1];
        int nCounts = genoCounts.length;
        if (nCounts > 0) {
            int ntaxa = sitegeno.length;
            int maxCount = genoCounts[nCounts - 1];
            for (int t = 0; t < ntaxa; t++) {
                if (sitegeno[t] == NN) {
                    int ranval = randomizer.nextInt(maxCount);
                    int genoIndex = 0;
                    while (ranval > genoCounts[genoIndex])
                        genoIndex++;
                    sitegeno[t] = siteGenotypes[genoIndex];
                }
            }
            return Optional.of(sitegeno);
        } else return Optional.empty();
    }
    
    public static Object[] byteCounts(byte[] genotypes) {
        //Object[0] is a byte[] array of genotypes
        //Object[1] is an int[] array of the cumulative counts of the genotypes
        byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        int n = genotypes.length;
        Map<Byte, Long> byteCounts = IntStream.range(0, n).filter(i -> genotypes[i] != NN).boxed()
                .collect(Collectors.groupingBy(i -> new Byte(genotypes[i]), Collectors.counting()));
        int mapSize = byteCounts.size();
        byte[] geno = new byte[mapSize];
        int[] genocount = new int[mapSize];
        int ndx = 0;
        int sum = 0;
        for (Map.Entry<Byte, Long> me : byteCounts.entrySet()) {
            geno[ndx] = me.getKey().byteValue();
            sum += me.getValue().intValue();
            genocount[ndx] = sum;
            ndx++;
        }
        return new Object[] { geno, genocount };
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Random Imputation";
    }

    @Override
    public String getToolTipText() {
        return "Replace missing genotypes with a value drawn from the site genotype distribution.";
    }

}
