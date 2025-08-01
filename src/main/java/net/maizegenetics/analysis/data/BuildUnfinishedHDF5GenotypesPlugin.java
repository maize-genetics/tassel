package net.maizegenetics.analysis.data;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.SimpleDateFormat;
import java.util.Date;
import javax.swing.*;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
//import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.util.HDF5Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Opens an "unfinished" HDF5 genotypes file, on which .closeUnfinished() was called rather than .build()
 * and finalizes it by calling .build().
 * 
 * If an output file name is provided, then the input file is first copied to that, and then build() is called on that copy.
 * 
 * If provided, root level annotations DataSetName and DataSetDescription are also added. "__DATE__" in the provided
 * dataSet name is replaced a current date stamp ("_yyyyMMdd"). "__SNPS__" and "__TAXA__" in the dataSetDescription are
 * replaced with the number of sites and taxa, respectively.
 *
 * @author Jeff Glaubitz (jcg233@cornell.edu)
 */
public class BuildUnfinishedHDF5GenotypesPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(BuildUnfinishedHDF5GenotypesPlugin.class);
    
    String dataSetDescrip, date;

    private PluginParameter<String> inputGenotypes = new PluginParameter.Builder<>("i", null, String.class)
        .guiName("Input file")
        .required(true)
        .inFile()
        .description("Input, unfinished HDF5 genotype (*.h5) file to be fininalized")
        .build();
    private PluginParameter<String> outputGenotypes = new PluginParameter.Builder<>("o", null, String.class)
        .guiName("Output file")
        .required(false)
        .outFile()
        .description("Output, finished HDF5 genotype (*.h5) file which can be opened with the TASSEL5 GUI. __DATE__ is replaced with a _yyyyMMdd date stamp.")
        .build();
    private PluginParameter<String> dataSetName = new PluginParameter.Builder<>("name", null, String.class)
        .guiName("Data set name")
        .required(false)
        .description("(Optional) Short data set name to be added as an root level annotation under \"dataSetName\"")
        .build();
    private PluginParameter<String> dataSetDescription = new PluginParameter.Builder<>("desc", null, String.class)
        .guiName("Data set description")
        .required(false)
        .description("(Optional) Short data set description to be added as an root level annotation under \"dataSetDescription\"")
        .build();

    public BuildUnfinishedHDF5GenotypesPlugin() {
        super(null, false);
    }

    public BuildUnfinishedHDF5GenotypesPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, false);
    }
    
    @Override
    protected void preProcessParameters(DataSet input) {
        date = new SimpleDateFormat("yyyyMMdd").format(new Date());
        String outfile = outputFile();
        outputFile(outfile.replace("__DATE__", "_"+date));
    }

    @Override
    public DataSet processData(DataSet input) {
        String message = buildUnfinishedHDF5Genotypes();
        if(message != null) {
            myLogger.error(message);
            try {Thread.sleep(500);} catch(Exception e) {}
            throw new IllegalStateException(message);
        }
        fireProgress(100);
        return null;
    }

    private String buildUnfinishedHDF5Genotypes() {
        if (outputFile() == null) {
            myLogger.info("\n\nBuildUnfinishedHDF5GenotypesPlugin:\nFinalizing the following HDF5 genotype file:\n   "
                    +inputFile()+"\n\n");
        } else {
            myLogger.info("\n\nBuildUnfinishedHDF5GenotypesPlugin: Copying the HDF5 genotypes from the file:\n   "
                    +inputFile()
                    +"\n"+"and finalizing them in this output file:\n   "
                    +outputFile()+"\n\n");
        }
        if (dataSetDescription() != null) {
            dataSetDescrip = parseDataSetDescription(dataSetDescription()); // get the number of SNPs and taxa
        }
        GenotypeTableBuilder genoTable;
        if (outputFile() == null) {
            genoTable = GenotypeTableBuilder.getBuilder(inputFile());
        } else {
            String message = copyInputFile();
            if (message != null) return message;
            genoTable = GenotypeTableBuilder.getBuilder(outputFile());
        }
        
        if (dataSetName() != null) {
            genoTable.dataSetName(parseDataSetName(dataSetName()));
        }
        if (dataSetDescription() != null) {
            genoTable.dataSetDescription(dataSetDescrip);
        }
        
        genoTable.build();
        myLogger.info("\n\nFinished finalizing HDF5 genotypes file\n\n");
        return null;
    }
    
    private String parseDataSetName(String dataSetName) {
        return dataSetName.replace("__DATE__", "_"+date);
    }
    
    private String parseDataSetDescription(String dataSetDescrip) {
        IHDF5Reader h5Reader = HDF5Factory.openForReading(inputFile());
        int nSNPs = HDF5Utils.getHDF5PositionNumber(h5Reader);
        TaxaList tL = new TaxaListBuilder().buildFromHDF5Genotypes(h5Reader);
        int nTaxa = tL.numberOfTaxa();
        h5Reader.close();
        return dataSetDescrip.replace("__SNPS__", ""+nSNPs).replace("__TAXA__", ""+nTaxa).replace("__DATE__", date);
    }
    
    private String copyInputFile() {
        try {
            Files.copy(new File(inputFile()).toPath(), new File(outputFile()).toPath());
        } catch (IOException e) {
            myLogger.error("Could not copy file: " + e.getMessage());
            return "\n\nERROR: Could not copy input file to output file: " + e.getMessage() +"\n\n";
        }
        return null;
    }

    @Override
    public ImageIcon getIcon() {
//        URL imageURL = SeparatePlugin.class.getResource("/net/maizegenetics/analysis/images/Merge.gif");
//        if (imageURL == null) {
            return null;
//        } else {
//            return new ImageIcon(imageURL);
//        }
    }

    @Override
    public String getButtonName() {
        return "Split chromosomes from HDF5 genotype file";
    }

    @Override
    public String getToolTipText() {
        return "Split chromosomes from HDF5 genotype file";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(BuildUnfinishedHDF5GenotypesPlugin.class);
    // }

    /**
     * Input, unfinished HDF5 genotype (*.h5) file to be fininalized
     *
     * @return Input file
     */
    public String inputFile() {
        return inputGenotypes.value();
    }

    /**
     * Set Input file. Input, unfinished HDF5 genotype (*.h5)
     * file to be fininalized
     *
     * @param value Input file
     *
     * @return this plugin
     */
    public BuildUnfinishedHDF5GenotypesPlugin inputFile(String value) {
        inputGenotypes = new PluginParameter<>(inputGenotypes, value);
        return this;
    }

    /**
     * Output, finished HDF5 genotype (*.h5) file which can
     * be opened with the TASSEL5 GUI
     *
     * @return Output file
     */
    public String outputFile() {
        return outputGenotypes.value();
    }

    /**
     * Set Output file. Output, finished HDF5 genotype (*.h5)
     * file which can be opened with the TASSEL5 GUI
     *
     * @param value Output file
     *
     * @return this plugin
     */
    public BuildUnfinishedHDF5GenotypesPlugin outputFile(String value) {
        outputGenotypes = new PluginParameter<>(outputGenotypes, value);
        return this;
    }

    /**
     * (Optional) Short data set name to be added as an root
     * level annotation under "/DataSetName"
     *
     * @return Data set name
     */
    public String dataSetName() {
        return dataSetName.value();
    }

    /**
     * Set Data set name. (Optional) Short data set name to
     * be added as an root level annotation under "dataSetName"
     *
     * @param value Data set name
     *
     * @return this plugin
     */
    public BuildUnfinishedHDF5GenotypesPlugin dataSetName(String value) {
        dataSetName = new PluginParameter<>(dataSetName, value);
        return this;
    }

    /**
     * (Optional) Short data set description to be added as
     * an root level annotation under "dataSetDescription"
     *
     * @return Data set description
     */
    public String dataSetDescription() {
        return dataSetDescription.value();
    }

    /**
     * Set Data set description. (Optional) Short data set
     * description to be added as an root level annotation
     * under "/DataSetDescription"
     *
     * @param value Data set description
     *
     * @return this plugin
     */
    public BuildUnfinishedHDF5GenotypesPlugin dataSetDescription(String value) {
        dataSetDescription = new PluginParameter<>(dataSetDescription, value);
        return this;
    }

}
