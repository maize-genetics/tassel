/**
 * 
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Plugin to thin out sites based on their physical position on the chromosome.
 * Expects the user to enter a minimum distance between sites in base pair.
 * Accepts and writes .hmp.txt, .hmp.txt.gz, .vcf, .vcf.gz and .h5 files.
 * 
 * @author StefanReuscher
 * @author JeffGlaubitz
 * @author lcj34
 *
 */
public class ThinSitesByPositionPlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(GenosToABHPlugin.class);

    private PluginParameter<String> outfile= new PluginParameter.Builder<>("o", null, String.class)
            .required(true)
            .outFile()
            .guiName("Output file")
            .description("Output genotype file").build();
    private PluginParameter<Integer> minDist= new PluginParameter.Builder<>("minDist", null, Integer.class)
            .required(true)
            .guiName("Minimum distance")
            .description("Minimum distance in bp between adjacent sites")
            .build();
    
    private GenotypeTable myInput = null;
    
    public ThinSitesByPositionPlugin() {
        super(null, false);
    }

    public ThinSitesByPositionPlugin(Frame parentFrame, boolean isInteractive) {
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

    public DataSet processData(DataSet input) {
        String[] sitesToKeep = findSiteNamesToKeep(myInput); 
        GenotypeTable filteredGenos = null;
        if ((sitesToKeep != null) && (sitesToKeep.length != 0)) {
            filteredGenos = FilterGenotypeTable.getInstance(myInput, sitesToKeep);
            writeOutputGenos(filteredGenos);
        } else {
            myLogger.warn("WARNING - no sites kept, no output file written !!");
        }
        return null;
    }
    
    private String[] findSiteNamesToKeep(GenotypeTable genos){
        // This method traverses the Positions in each chromosome,
        // adding to a list those positions within the chromosome
        // that are the user specified distance apart.
        ArrayList<String> sitesToKeep = new ArrayList<String>();
        PositionList positions = genos.positions();
        Chromosome chrom = new Chromosome("-1");
        int prevPos = -1;
        for (Position currPos : positions){
            if (!currPos.getChromosome().equals(chrom)) {
                // Always keep the first position in a chromosome
                sitesToKeep.add(currPos.getSNPID());
                chrom = currPos.getChromosome();
                prevPos = currPos.getPosition();
            } else {
                if (currPos.getPosition() - prevPos >= minDist()) {;
                    sitesToKeep.add(currPos.getSNPID());
                    prevPos = currPos.getPosition();
                } // do nothing if distance between positions < minDistß
            }           
        }  
        myLogger.info("ThinSitesByPosition: original number of sites: " + positions.size() 
                + " number of sites kept: " + sitesToKeep.size());
        String[] stringArray = sitesToKeep.toArray(new String[sitesToKeep.size()]);
        return stringArray;
    }
    
    private void writeOutputGenos (GenotypeTable genos){
        // Output format is based on name given for output file
        // If the extension is not a TASSEL supported extension, default to hapmap
        String fileName = outFile();
        if (fileName.endsWith(".h5")) {
            ExportUtils.writeGenotypeHDF5(genos,fileName,true);            
        } else if (fileName.endsWith("hmp.txt.gz") || fileName.endsWith("hmp.txt")) {
            ExportUtils.writeToHapmap(genos, fileName);
        } else if (fileName.endsWith(".vcf") || fileName.endsWith(".vcf.gz")) {
            ExportUtils.writeToVCF(genos, fileName, true);
        } else {
            fileName = fileName + ".hmp.txt";
            myLogger.warn("File extension not recognized, writing to hapmap file:\n    " + fileName);
            ExportUtils.writeToHapmap(genos, fileName);
        }
    }
    
    @Override
    public ImageIcon getIcon() {
        URL imageURL = HetsToUnknownPlugin.class.getResource("/net/maizegenetics/analysis/images/homozygous.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Thin Sites by Position";
    }

    @Override
    public String getToolTipText() {
        return "This sites based on their physical position on the chromosome";
    }

    /**
     * Output genotype file
     *
     * @return Output file
     */
    public String outFile() {
        return outfile.value();
    }

    /**
     * Set Output file. Output genotype file
     *
     * @param value Output file
     *
     * @return this plugin
     */
    public ThinSitesByPositionPlugin outfile(String value) {
        outfile = new PluginParameter<>(outfile, value);
        return this;
    }

    /**
     * Minimum distance in bp between adjacent sites
     *
     * @return Minimum distance
     */
    public Integer minDist() {
        return minDist.value();
    }

    /**
     * Set Minimum distance. Minimum distance in bp between
     * adjacent sites
     *
     * @param value Minimum distance
     *
     * @return this plugin
     */
    public ThinSitesByPositionPlugin minDist(Integer value) {
        minDist = new PluginParameter<>(minDist, value);
        return this;
    }

    @Override
    public String getCitation() {
        return "Stefan Reuscher, Jeff Glaubitz, Lynn Johnson (2015) First Annual TASSEL Hackathon";
    }

//     The following getters and setters were auto-generated.
//     Please use this method to re-generate.
    
//     public static void main(String[] args) {
//         GeneratePluginCode.generate(ThinSitesByPositionPlugin.class);
//     }

    /**
     * Convenience method to run plugin with one return object.
     */
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }
}
