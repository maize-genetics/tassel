/*
 * SplitHDF5ByChromosomePlugin
 */
package net.maizegenetics.analysis.data;

// standard imports for plugins
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import net.maizegenetics.util.ArgsEngine;
import javax.swing.*;
import java.awt.*;

// specifically needed for this plugin
import java.io.File;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.plugindef.PluginParameter;

/**
 *
 * @author Jeff Glaubitz (jcg233@cornell.edu)
 */
public class SplitHDF5ByChromosomePlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(SplitHDF5ByChromosomePlugin.class);

    private PluginParameter<String> myInputGenotypes = new PluginParameter.Builder<>("i", null, String.class)
        .guiName("Input HDF5 Genotype File")
        .required(true)
        .inFile()
        .description("Input HDF5 genotype (*.h5) file to be split by chromosome. The output files will be named *_chr#.h5 (where # = chromosome number)")
        .build();
    private PluginParameter<Boolean> myIgnoreDepth = new PluginParameter.Builder<>("iD", false, Boolean.class)
        .guiName("Ignore Depth")
        .required(false)
        .description("If genotypic depth information is present in the input file, ignore it (i.e., do not write depth information to the output files). Default: keep depth.")
        .build();
    
    private ArgsEngine myArgsEngine = null;
    private String inHDF5FileName = null;
    private boolean keepDepth = true;

    public SplitHDF5ByChromosomePlugin() {
        super(null, false);
    }

    public SplitHDF5ByChromosomePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, false);
    }
//
//    private void printUsage() {
//        myLogger.info(
//            "\n\n"
//            +"The options for the TASSEL5 SplitHDF5ByChromosomePlugin are as follows:\n"
//            +"  -i   Input HDF5 genotype (*.h5) file to be split by chromosome. The\n"
//            +"       output files will be named *_chr#.h5 (where # = chromosome number).\n"
//                     +"\n"
//            +"  -iD  Ignore depth. If genotypic depth information is present in the input\n"
//            +"       input files, ignore it (i.e., do not write depth information to the \n"
//            +"       output files).  Default: keep depth.\n"
//            +"\n\n"
//        );
//    }
//
//    @Override
//    public void setParameters(String[] args) {
//        if (args == null || args.length == 0) {
//            printUsage();
//            try {Thread.sleep(500);} catch(InterruptedException e) {}
//            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
//        }
//        if (myArgsEngine == null) {
//            myArgsEngine = new ArgsEngine();
//            myArgsEngine.add("-i", "--input-file", true);
//            myArgsEngine.add("-iD", "--ignore-depth", false);
//        }
//        myArgsEngine.parse(args);
//        inHDF5FileName = myArgsEngine.getString("-i");
//        File inFile = null;
//        if (inHDF5FileName != null) {
//            inFile = new File(inHDF5FileName);
//            if (!inFile.exists()) {
//                printUsage();
//                try {Thread.sleep(500);} catch(InterruptedException e) {}
//                String errMessage = "\nThe input HDF5 genotype (*.h5) file name you supplied does not exist:\n"+inHDF5FileName+"\n\n";
//                myLogger.error(errMessage);
//                try {Thread.sleep(500);} catch(InterruptedException e) {}
//                throw new IllegalArgumentException(errMessage);
//            }
//        } else {
//            printUsage();
//            try {Thread.sleep(500);} catch(InterruptedException e) {}
//            String errMessage = "\nPlease specify an input HDF5 genotype (*.h5) file to be split by chromosome (-i option)\n\n";
//            myLogger.error(errMessage);
//            try {Thread.sleep(500);} catch(InterruptedException e) {}
//            throw new IllegalArgumentException(errMessage);
//        }
//        if (myArgsEngine.getBoolean("-iD")) {
//            keepDepth = false;
//        }
//    }

    @Override
    public DataSet processData(DataSet input) {
        String message = splitHDF5GenoTableByChr();
        if(message != null) {
            myLogger.error(message);
            try {Thread.sleep(500);} catch(Exception e) {}
            throw new IllegalStateException(message);
        }
        fireProgress(100);
        return null;
    }

    private String splitHDF5GenoTableByChr() {
        System.out.println("\n\nSplitHDF5ByChromosomePlugin:\nSplitting the following input file by chromosome:");
        System.out.println("  "+inputHDF5GenotypeFile()+"\n\n");
        GenotypeTable inGenos = ImportUtils.readGuessFormat(inputHDF5GenotypeFile());
        Chromosome[] chrs = inGenos.chromosomes();
        for (Chromosome chr : chrs) {
            GenotypeTable genosForChr = FilterGenotypeTable.getInstance(inGenos, chr);
            String outFileName;
            if (ignoreDepth()) { 
                outFileName = inputHDF5GenotypeFile().replaceFirst("\\.h5$", "_noDepth_chr"+chr.getName()+".h5");
            } else {
                outFileName = inputHDF5GenotypeFile().replaceFirst("\\.h5$", "_chr"+chr.getName()+".h5");
            }
            ExportUtils.writeGenotypeHDF5(genosForChr, outFileName, !ignoreDepth());
            System.out.println("    Genotypes from chromosome "+chr.getName()+" written to file");
            System.out.println("        "+outFileName);
        }
        System.out.println("\n\nSplitHDF5ByChromosomePlugin: Finished splitting chromosomes.\n\n");
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
    //     GeneratePluginCode.generate(SplitHDF5ByChromosomePlugin.class);
    // }

    /**
     * Input HDF5 genotype (*.h5) file to be split by chromosome.
     * The output files will be named *_chr#.h5 (where # =
     * chromosome number)
     *
     * @return Input HDF5 Genotype File
     */
    public String inputHDF5GenotypeFile() {
        return myInputGenotypes.value();
    }

    /**
     * Set Input HDF5 Genotype File. Input HDF5 genotype (*.h5)
     * file to be split by chromosome. The output files will
     * be named *_chr#.h5 (where # = chromosome number)
     *
     * @param value Input HDF5 Genotype File
     *
     * @return this plugin
     */
    public SplitHDF5ByChromosomePlugin inputHDF5GenotypeFile(String value) {
        myInputGenotypes = new PluginParameter<>(myInputGenotypes, value);
        return this;
    }

    /**
     * If genotypic depth information is present in the input
     * file, ignore it (i.e., do not write depth information
     * to the output files). Default: keep depth.
     *
     * @return Ignore Depth
     */
    public Boolean ignoreDepth() {
        return myIgnoreDepth.value();
    }

    /**
     * Set Ignore Depth. If genotypic depth information is
     * present in the input file, ignore it (i.e., do not
     * write depth information to the output files). Default:
     * keep depth.
     *
     * @param value Ignore Depth
     *
     * @return this plugin
     */
    public SplitHDF5ByChromosomePlugin ignoreDepth(Boolean value) {
        myIgnoreDepth = new PluginParameter<>(myIgnoreDepth, value);
        return this;
    }
}
