/**
 *
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Plugin to convert genotypes to parental comparisons and store them in an
 * output file. The user can get outputs as A/H/B, 0/1/2, or 0/0.5/1 for a
 * genotype matching parent A's genotype, heterozygous, or parent B's genotype
 * respectively. If the genotype is neither A or B or a het combination of A/B
 * (B/A) then it is coded as "NA".
 *
 * @author Stefan Reuscher
 * @author jeff Glaubitz
 * @author lynn Johnson
 * @author Josh Guy
 *
 */
public class GenosToABHPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(GenosToABHPlugin.class);
    private ArrayList<Integer> parentAIndices = null;
    private ArrayList<Integer> parentBIndices = null;

    public static enum OUTPUT_CHECK {
        c, i, r
    }

    private PluginParameter<String> outfile = new PluginParameter.Builder<>("o", null, String.class)
            .required(true).outFile().guiName("Output file").description("Output genotype file").build();
    private PluginParameter<String> parentA = new PluginParameter.Builder<>("parentA", null, String.class)
            .required(true).guiName("Parent A").inFile()
            .description("The full name of file containing list of taxa names for parent A").build();
    private PluginParameter<String> parentB = new PluginParameter.Builder<>("parentB", null, String.class)
            .required(true).guiName("Parent B").inFile()
            .description("The full name of file containing list of taxa names for parent B").build();
    private PluginParameter<OUTPUT_CHECK> outputFormat = new PluginParameter.Builder<>("outputFormat", OUTPUT_CHECK.c, OUTPUT_CHECK.class)
            .required(false).guiName("Output Format")
            .description("Enter \"c\" for A,H,B to code for parent A, het, and parent B; Enter \"i\" for 0, 1, 2 to code for parent A, het, and parent B; Enter \"r\" for 0, 0.5, 1 to code for parent A, het, and parent B")
            .range(OUTPUT_CHECK.values()).build();

    private GenotypeTable myInput = null;

    public GenosToABHPlugin() {
        super(null, false);
    }

    public GenosToABHPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        if (input == null) {
            throw new IllegalArgumentException("GenosToABHPlugin: preProcessParameters: Please select one Genotype Table.");
        }
        List<Datum> genotypeTables = input.getDataOfType(GenotypeTable.class);
        if (genotypeTables.size() == 1) {
            myInput = (GenotypeTable) genotypeTables.get(0).getData();
        } else {
            throw new IllegalArgumentException("GenosToABHPlugin: preProcessParameters: Please select one Genotype Table.");
        }
    }

    /**
     * The main method. Plugin to convert genotypes to parental comparisons and
     * store them in an output file. The user can get outputs as A/H/B, 0/1/2,
     * or 0/0.5/1 for a genotype matching parent A's genotype, heterozygous, or
     * parent B's genotype respectively. If the genotype is neither A or B or a
     * het combination of A/B (B/A) then it is coded as "NA".
     */
    @Override
    public DataSet processData(DataSet input) {

        parentAIndices = getParentIndex(myInput, parentA());
        parentBIndices = getParentIndex(myInput, parentB());
        if (parentAIndices == null || parentAIndices.size() == 0 || parentBIndices == null
                || parentBIndices.size() == 0) {
            return null;
        }

        byte[] parentAGenos = new byte[myInput.numberOfSites()];
        byte[] parentBGenos = new byte[myInput.numberOfSites()];
        byte[][] hets = new byte[2][myInput.numberOfSites()];

        createParentalByteGenos(myInput, parentAGenos, parentBGenos, hets);

        myLogger.info(String.format("GenosToABHPlugin: number Of sites:%d  number of taxa:%d %n", myInput.numberOfSites(), myInput.numberOfTaxa()));

        // Write converted genotypes to a text file
        writeConvertedGenos(myInput, parentAGenos, parentBGenos, hets);

        return null;
    }

    // Assumes a file format that contains NO header line.  Each line
    // contains one parental taxon name.
    private ArrayList<Integer> getParentIndex(GenotypeTable genos, String parentFile) {
        try (BufferedReader reader = Utils.getBufferedReader(parentFile)) {
            ArrayList<String> parentNames = new ArrayList<String>();
            String line;
            // Read each line of the parentFile, add the parent taxon
            // to the ArrayList of parent names
            while ((line = reader.readLine()) != null) {
                parentNames.add(line);
            }
            // if there are parentNames, find their index in the genotype table
            if (parentNames.size() > 0) {
                ArrayList<Integer> parentIndices = new ArrayList<Integer>();
                for (String parent : parentNames) {
                    int index = genos.taxa().indexOf(parent);
                    if (index >= 0) {
                        parentIndices.add(index);
                    } else {
                        myLogger.error("Parent " + parent + " missing from the genotype file");
                        return null;
                    }
                }
                System.out.println("File " + parentFile + " has " + parentIndices.size() + " elements");
                return parentIndices;
            } else {
                myLogger.error("No parent names found in file:  " + parentFile);
            }

        } catch (IOException exc) {
            exc.printStackTrace();
        }
        return null;
    }

    // Create byte arrays to hold consensus parental genotypes
    // Rejected sites will have parentAGenos[site] set to N (don't output these to converted ABH genotype file)
    // Also create a 2-D byte array holding the corresponding heterozygous call for each site (in both possible phases)
    public void createParentalByteGenos(GenotypeTable genos, byte[] parentAGenos, byte[] parentBGenos, byte[][] hets) {
        int numAccepted = 0;
        int numRejected = 0;
        // Populate arrays with genotypes
        for (int site = 0; site < genos.numberOfSites(); site++) {

            // StringBuilder for output of parental genotypes of accepted and rejected sites to log
            StringBuilder strB = new StringBuilder(genos.siteName(site));

            // Resolve consensus genotypes for parents A and B
            byte parentAGeno = getParentalGenotype(genos, parentAIndices, site);
            byte parentBGeno = getParentalGenotype(genos, parentBIndices, site);

            // reject sites that are not polymorphic or unknown or heterozygous in either parent
            if (parentAGeno == parentBGeno
                    || parentAGeno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE
                    || parentBGeno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE
                    || GenotypeTableUtils.isHeterozygous(parentAGeno)
                    || GenotypeTableUtils.isHeterozygous(parentBGeno)) {
                // reject
                parentAGenos[site] = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
                numRejected++;
                strB.append("\treject");
            } else {
                // accept
                parentAGenos[site] = parentAGeno;
                parentBGenos[site] = parentBGeno;

                // Takes the consensus parental genotypes and gives us the 2 possible heterozygotes           
                hets[0][site] = GenotypeTableUtils.getDiploidValuePhased(parentAGeno, parentBGeno);
                hets[1][site] = GenotypeTableUtils.getDiploidValuePhased(parentBGeno, parentAGeno);
                strB.append("\taccept");
                numAccepted++;
            }
            strB.append("\t" + NucleotideAlignmentConstants.getNucleotideIUPAC(parentAGeno) + ":");
            for (Integer parentAIdx : parentAIndices) {
                strB.append(genos.genotypeAsString(parentAIdx, site));
            }
            strB.append("\t" + NucleotideAlignmentConstants.getNucleotideIUPAC(parentBGeno) + ":");
            for (Integer parentBIdx : parentBIndices) {
                strB.append(genos.genotypeAsString(parentBIdx, site));
            }
            myLogger.info(strB.toString());
        }
        myLogger.info("Number of accepted sites: " + numAccepted + ", number of Rejected sites: " + numRejected);
    }

    public byte getParentalGenotype(GenotypeTable genos, ArrayList<Integer> parentIndices, int site) {
        boolean finalIsUnknown = false;

        byte finalParent = genos.genotype(parentIndices.get(0), site);
        if (finalParent == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
            finalIsUnknown = true;
        }
        for (int idx = 1; idx < parentIndices.size(); idx++) {
            byte nextGeno = genos.genotype(parentIndices.get(idx), site);
            if (nextGeno != finalParent) {
                if (finalIsUnknown) {
                    finalParent = nextGeno;
                    finalIsUnknown = false;
                } else if (nextGeno != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    // Mismatch - what do we set it to?
                    finalParent = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
                    return finalParent;
                }
            }
        }
        return finalParent;
    }

    private void writeConvertedGenos(GenotypeTable genos, byte[] parentAGenos, byte[] parentBGenos, byte[][] hets) {
        try (BufferedWriter bw = Utils.getBufferedWriter(outfile())) {

            // create the first line of the file - this is the header line,
            // it contains the site names.  Skip all sites where parentAGenos of 
            // that site == UNKNOWN_DIPLOID_ALLELE
            StringBuilder strB = new StringBuilder("id");
            for (int site = 0; site < genos.numberOfSites(); site++) {
                if (parentAGenos[site] == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    continue;
                }
                strB.append("," + genos.siteName(site));
            }
            writeLine(strB, bw); // writeline appends the newline

            // Second line contains the chromosome names
            strB = new StringBuilder("NA");
            for (int site = 0; site < genos.numberOfSites(); site++) {
                if (parentAGenos[site] == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    continue;
                }
                strB.append("," + genos.chromosomeName(site));
            }
            writeLine(strB, bw);

            // write the converted genotype for each taxon/site
            for (int taxon = 0; taxon < genos.numberOfTaxa(); taxon++) {
                if (parentAIndices.contains(taxon) || parentBIndices.contains(taxon)) {
                    continue; // skip the parents
                }
                strB = new StringBuilder(genos.taxaName(taxon));
                for (int site = 0; site < genos.numberOfSites(); site++) {
                    if (parentAGenos[site] == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                        continue;
                    }
                    byte geno = genos.genotype(taxon, site);

                    //Check desired output and append the corresponding value
                    if (outputFormat.value() == OUTPUT_CHECK.c || outputFormat == null) {
                        if (geno == parentAGenos[site]) {
                            strB.append(",A");
                        } else if (geno == parentBGenos[site]) {
                            strB.append(",B");
                        } else if (geno == hets[0][site] || geno == hets[1][site]) {
                            strB.append(",H");
                        } else {
                            strB.append(",NA");
                        }
                    } else if (outputFormat.value() == OUTPUT_CHECK.i) {
                        if (geno == parentAGenos[site]) {
                            strB.append(",0");
                        } else if (geno == parentBGenos[site]) {
                            strB.append(",2");
                        } else if (geno == hets[0][site] || geno == hets[1][site]) {
                            strB.append(",1");
                        } else {
                            strB.append(",NA");
                        }
                    } else if (outputFormat.value() == OUTPUT_CHECK.r) {
                        if (geno == parentAGenos[site]) {
                            strB.append(",0");
                        } else if (geno == parentBGenos[site]) {
                            strB.append(",1");
                        } else if (geno == hets[0][site] || geno == hets[1][site]) {
                            strB.append(",0.5");
                        } else {
                            strB.append(",NA");
                        }
                    }
                }
                writeLine(strB, bw);
            }

        } catch (IOException e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("GenosToABHPlugin: writeConvertedGenos: Problem writing file: " + outfile() + ": " + e.getMessage());
        }
    }

    private void writeLine(StringBuilder strB, BufferedWriter writer) {
        try {
            writer.write(strB.toString() + "\n");
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    @Override
    public String getToolTipText() {
        return "Convert Genotypes to A/B/H";
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = GenosToABHPlugin.class.getResource("/net/maizegenetics/analysis/images/homozygous.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "ABH Genotype";
    }

    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(GenosToABHPlugin.class);
    // }
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(GenosToABHPlugin.class);
    // }
    /**
     * Convenience method to run plugin with one return object.
     */
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Output genotype file with ABH encoding
     *
     * @return Output file
     */
    public String outfile() {
        return outfile.value();
    }

    /**
     * Set Output file. Output genotype file with ABH encoding
     *
     * @param value Output file
     *
     * @return this plugin
     */
    public GenosToABHPlugin outfile(String value) {
        outfile = new PluginParameter<>(outfile, value);
        return this;
    }

    /**
     * The full name of parent to be encoded as A
     *
     * @return Parent A
     */
    public String parentA() {
        return parentA.value();
    }

    /**
     * Set Parent A. The full name of parent to be encoded as A
     *
     * @param value Parent A
     *
     * @return this plugin
     */
    public GenosToABHPlugin parentA(String value) {
        parentA = new PluginParameter<>(parentA, value);
        return this;
    }

    /**
     * The full name of parent to be encoded as B
     *
     * @return Parent B
     */
    public String parentB() {
        return parentB.value();
    }

    /**
     * Set Parent B. The full name of parent to be encoded as B
     *
     * @param value Parent B
     *
     * @return this plugin
     */
    public GenosToABHPlugin parentB(String value) {
        parentB = new PluginParameter<>(parentB, value);
        return this;
    }

    /**
     * The output format enum
     *
     * @return Output Format
     */
    public OUTPUT_CHECK outputFormat() {
        return outputFormat.value();
    }

    /**
     * Set Output Format. The output format enum
     *
     * @param value Output Format
     *
     * @return this plugin
     */
    public GenosToABHPlugin outputFormat(OUTPUT_CHECK value) {
        outputFormat = new PluginParameter<>(outputFormat, value);
        return this;
    }

    @Override
    public String getCitation() {
        return "Stefan Reuscher, Jeff Glaubitz, Lynn Johnson (2015) First Annual TASSEL Hackathon";
    }

    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/GenosToABH/GenosToABHPlugin";
    }

}
