/*
 * MergeHDF5GenotypesSameSitesPlugin
 */
package net.maizegenetics.analysis.data;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.DirectoryCrawler;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.net.URL;
import java.util.Arrays;
import java.util.List;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.BasicGenotypeMergeRule;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;

/**
 *
 * @author Jeff Glaubitz (jcg233@cornell.edu)
 */
public class MergeHDF5GenotypesSameSitesPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(MergeHDF5GenotypesSameSitesPlugin.class);
    private ArgsEngine myArgsEngine = null;
    private String[] inHD5F5GenoFileNames = null;
    private String mergedHDF5GenoFileName;
    private boolean ignoreDepth = false;

    public MergeHDF5GenotypesSameSitesPlugin() {
        super(null, false);
    }

    public MergeHDF5GenotypesSameSitesPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        myLogger.info(
            "\n\n"
            +"The options for the TASSEL5 MergeHDF5GenotypesSameSitesPlugin are as follows:\n"
            +"  -i   Input folder containing HDF5 genotype (*.h5) files to be merged. These\n"
            +"       must all contain exactly the same sites.\n"
                     +"\n"
            +"  -o   Output, merged HDF5 genotype file name (*.h5). Must in a different\n"
            +"       folder than the input folder, that lies outside of the input folder.\n"
                     +"\n"
            +"  -iD  Ignore depth. If genotypic depth information is present in any of the\n"
            +"       input files, ignore it. If none of the files contain depth, then this\n"
            +"       option is not strictly necessary, but it will likely improve speed.\n"
            +"       Default: use depth (under the assumption that it is present in all\n"
            +"       input files)."
            +"\n\n"
        );
    }

    @Override
    public void setParameters(String[] args) {
        if (args == null || args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--input-folder", true);
            myArgsEngine.add("-o", "--output-file", true);
            myArgsEngine.add("-iD", "--ignore-depth", false);
        }
        myArgsEngine.parse(args);
        String inDirName = myArgsEngine.getString("-i");
        File inDirectory = null;
        if (inDirName != null) {
            inDirectory = new File(inDirName);
            if (!inDirectory.isDirectory()) {
                printUsage();
                String errMessage = "\nThe input folder name you supplied is not a directory:\n"+inDirName+"\n\n";
                myLogger.error(errMessage);
                throw new IllegalArgumentException(errMessage);
            }
            inHD5F5GenoFileNames = DirectoryCrawler.listFileNames(".*\\.h5$", inDirectory.getAbsolutePath());
            if (inHD5F5GenoFileNames.length == 0 || inHD5F5GenoFileNames == null) {
                printUsage();
                String errMessage = "\nNo files with the extension .h5 (HDF5 genotype files) were found in the input folder:\n  "
                        +inDirName+"\n"+"or in any of its sub-folders\n\n";
                myLogger.error(errMessage);
                throw new IllegalArgumentException(errMessage);
            } else if (inHD5F5GenoFileNames.length < 2) {
                printUsage();
                String errMessage = "\nThere must be at least 2 files with the extension .h5 (HDF5 genotype files) in the input folder:\n  "
                        +inDirName+"\n"+"or its sub-folders\n\n";
                myLogger.error(errMessage);
                throw new IllegalArgumentException(errMessage);
            } else {
                Arrays.sort(inHD5F5GenoFileNames);
                myLogger.info("\n\nMergeHDF5GenotypesSameSitesPlugin:\n\nThe following HDF5 genotype files were found in the input folder (and sub-folders):");
                for (String filename : inHD5F5GenoFileNames) {
                    System.out.println("   "+filename);
                }
                System.out.println("\n");
            }
        } else {
            printUsage();
            String errMessage = "\nPlease specify an input directory containing HDF5 genotype files to be merged (option -i)\n\n";
            myLogger.error(errMessage);
            throw new IllegalArgumentException(errMessage);
        }

        if (myArgsEngine.getBoolean("-o")) {
            mergedHDF5GenoFileName = myArgsEngine.getString("-o");
            File outFile = new File(mergedHDF5GenoFileName);
            if (outFile.exists()) {
                printUsage();
                String errMessage = "\nThe output file (option -o) must not already exist.\n\n";
                myLogger.error(errMessage);
                throw new IllegalArgumentException(errMessage);
            }
            if (outFile.getAbsolutePath().startsWith(inDirectory.getAbsolutePath())) {
                printUsage();
                String errMessage = "\nThe output file (option -o) must not lie within the input folder.\n\n";
                myLogger.error(errMessage);
                throw new IllegalArgumentException(errMessage);
            }
            if (!mergedHDF5GenoFileName.endsWith(".h5")) {
                printUsage();
                String errMessage = "\nThe output file name (option -o) must end with the extension \".h5\"\n\n";
                myLogger.error(errMessage);
                throw new IllegalArgumentException(errMessage);
            }
        } else {
            printUsage();
            String errMessage = "\nPlease specify the output file name (option -o), using the extension \".h5\"\n\n";
            myLogger.error(errMessage);
            throw new IllegalArgumentException(errMessage);
        }
        
        if (myArgsEngine.getBoolean("-iD")) {
            ignoreDepth = true;
        }
    }

    @Override
    public DataSet performFunction(DataSet input) {
        String message = mergeHDF5GenoTables();
        if(message != null) {
            printUsage();
            try {Thread.sleep(500);} catch(Exception e) {}
            myLogger.error(message);
            try {Thread.sleep(500);} catch(Exception e) {}
            throw new UnsupportedOperationException(message);
        }
        fireProgress(100);
        return null;
    }

    // todo TAS-54 covers some of these topics, as we want general merging rule sets
    public String mergeHDF5GenoTables() {
        IHDF5Reader reader = HDF5Factory.openForReading(inHD5F5GenoFileNames[0]);
        PositionList firstPosList = PositionListBuilder.getInstance(reader);
        GenotypeTableBuilder mrgGenos 
            = GenotypeTableBuilder.getTaxaIncrementalWithMerging(mergedHDF5GenoFileName,firstPosList,new BasicGenotypeMergeRule(0.01));
        for (int i = 0; i < inHD5F5GenoFileNames.length; i++) {
            if (i > 0) {
                reader = HDF5Factory.openForReading(inHD5F5GenoFileNames[i]);
                if (!samePositions(firstPosList, PositionListBuilder.getInstance(reader))) {
                    return "Error: "+inHD5F5GenoFileNames[i]+" does not contain the same positions as "+inHD5F5GenoFileNames[0];
                }
            }
            List<String> taxaNames = HDF5Utils.getAllTaxaNames(reader);
            for (String taxonName : taxaNames) {
                Taxon taxon = HDF5Utils.getTaxon(reader, taxonName);
                byte[] genos = HDF5Utils.getHDF5GenotypesCalls(reader,taxonName);
                byte[][] depth = null;
                if (!ignoreDepth) {
                    depth = HDF5Utils.getHDF5GenotypesDepth(reader,taxonName);
                }
                try {
                    mrgGenos.addTaxon(taxon, genos, depth);
                } catch (Exception e) {
                    return "Error: "+ e;
                }
            }
            System.out.println("...finished adding "+inHD5F5GenoFileNames[i]+" to the output file");
        }
        mrgGenos.build();
        System.out.println("\n\nFinished merging genotypes to output file:");
        System.out.println("  "+mergedHDF5GenoFileName+"\n\n");
        return null;
    }
    
    // if PositionList had an equals() method this wouldn't be necessary
    private boolean samePositions(PositionList PosListA, PositionList PosListB) {
        if (PosListA.size() != PosListB.size()) {
            return false;
        }
        for (int j = 0; j < PosListA.size(); j++) {
            // only GeneralPosition seems to have an appropriate equals() method
            Position genPos1 = new GeneralPosition.Builder(PosListA.get(j)).build();
            Position genPos2 = new GeneralPosition.Builder(PosListB.get(j)).build();
            if(!genPos1.equals(genPos2)) return false;
        }
        return true;
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = SeparatePlugin.class.getResource("/net/maizegenetics/analysis/images/Merge.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Merge GenotypeTables";
    }

    @Override
    public String getToolTipText() {
        return "Merge GenotypeTables";
    }
}
