/*
 * ProductionPipeline
 */
package net.maizegenetics.analysis.gbs;

import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.util.LoggingUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FilenameFilter;

import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;

/**
 *
 * @author Terry Casstevens
 */
public class ProductionPipeline extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(ProductionPipeline.class);
    private static final SimpleDateFormat LOGGING_DATE_FORMAT = new SimpleDateFormat("yyyyMMdd HH:mm:ss");
    private static final String READY_FILE_NAME = "ready.txt";
    private static final String LOCK_FILE_NAME = "lock.txt";
    private static final String SUMMARY_LOG_FILE = "ProductionPipeline.log";

    private PluginParameter<String> myInputDirectory = new PluginParameter.Builder<>("inputDirectory", null, String.class).required(true).inDir()
            .description("Input directory containing subdirectories with fastq AND/OR qseq files").build();
    private PluginParameter<String> myEnzyme = new PluginParameter.Builder<>("enzyme", null, String.class).required(true)
            .description("Enzyme used to create the GBS library").build();
    private PluginParameter<String> myProductionTOPM = new PluginParameter.Builder<>("productionTOPM", null, String.class).required(true).inFile()
            .description("Physical map file containing tags and corresponding variants (production TOPM)").build();
    private PluginParameter<String> myOutputGenotypeFile = new PluginParameter.Builder<>("outputGenotypeFile", null, String.class).required(true).outFile()
            .description("Output (target) HDF5 genotypes file to add new genotypes to (new file created if it doesn't exist)").build();
    private PluginParameter<String> myArchiveDirectory = new PluginParameter.Builder<>("archiveDirectory", null, String.class).required(true).outDir()
            .description("Archive directory where to move processed files").build();

    private String myOutputDirectory;
    private BufferedWriter mySummaryLogFile;
    private String date;

    public ProductionPipeline(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public void postProcessParameters() {
        if (myOutputGenotypeFile.value() != null) {
            myOutputDirectory = Utils.getDirectory(myOutputGenotypeFile.value());
            setupLogfile();
        }
        myLogger.info(getTimeStamp());
        date = new SimpleDateFormat("yyyyMMdd").format(new Date());
    }

    @Override
    public DataSet processData(DataSet input) {

        String lockFilename = inputDirectory() + File.separator + LOCK_FILE_NAME;

        try {

            if (!new File(lockFilename).createNewFile()) {
                myLogger.warn("Production Pipeline already running.  File exists: " + lockFilename + "  Aborting...");
                return null;
            }

            String logFilename = inputDirectory() + File.separator + SUMMARY_LOG_FILE;
            mySummaryLogFile = Utils.getBufferedWriter(logFilename, true);

            File inputDirectory = new File(inputDirectory());
            String[] directories = inputDirectory.list(new FilenameFilter() {
                @Override
                public boolean accept(File current, String name) {
                    return new File(current, name).isDirectory();
                }
            });

            if(directories.length==0){
                myLogger.warn("No subdirectories found within " + inputDirectory() + "; no SNPs will be called");
                return null;
            }

            for (String current : directories) {
                String fullDirName = inputDirectory() + File.separator + current;
                try {
                    if (processSubDirectory(fullDirName)) {
                        File currentLocation = new File(fullDirName);
                        File newLocation = new File(archiveDirectory()+current+"_"+date);  // there may be more lanes from the same flowcell released at a later date
                        currentLocation.renameTo(newLocation);
                        myLogger.info("Moved : " + currentLocation.getAbsolutePath() + " to: " + newLocation.getAbsolutePath());
                        writeToSummaryLogFile("Moved : " + currentLocation.getAbsolutePath() + " to: " + newLocation.getAbsolutePath());
                    }
                } catch (Exception e) {
                    writeToSummaryLogFile("Production Pipeline Failed: " + fullDirName);
                    myLogger.error(e.getMessage(), e);
                }

            }

            return null;

        } catch (Exception ex) {
            writeToSummaryLogFile("Problem Running Production Pipeline: " + ex.getMessage());
            myLogger.error(ex.getMessage(), ex);
            return null;
        } finally {
            LoggingUtils.closeLogfile();
            try {
                mySummaryLogFile.close();
            } catch (Exception e) {
                // do nothing
            }
            new File(lockFilename).delete();
        }

    }

    private boolean processSubDirectory(String subDirectory) {

        writeToSummaryLogFile("----------- Production Pipeline Started: " + subDirectory);

        String readyFilename = subDirectory + File.separator + READY_FILE_NAME;
        File readyFile = new File(readyFilename);
        if (readyFile.exists()) {
            myLogger.info("Processing directory: " + subDirectory);
            readyFile.delete();
        } else {
            myLogger.warn("This directory is not ready yet: " + subDirectory);
            writeToSummaryLogFile("Ready File not Found; " + readyFilename);
            return false;
        }

        String keyFile = subDirectory + File.separator + Utils.getFilename(subDirectory) + "_key.txt";
        if (!new File(keyFile).exists()) {
            myLogger.error("Keyfile doesn't exist: " + keyFile);
            writeToSummaryLogFile("Keyfile doesn't exist: " + keyFile);
            return false;
        }

        String[] rawSeqFileNames = DirectoryCrawler.listFileNames(ProductionSNPCallerPlugin.rawSeqFileNameRegex, subDirectory);
        if ((rawSeqFileNames == null) || (rawSeqFileNames.length == 0)) {
            myLogger.warn("No sequence files in directory: " + subDirectory);
            writeToSummaryLogFile("No sequence files in directory");
            return false;
        }

        String[] args = getPluginArgs(subDirectory, keyFile);

        myLogger.info("Raw Sequence Files: " + Arrays.deepToString(rawSeqFileNames));
        myLogger.info("Parameters Passed to ProductionSNPCallerPlugin: " + Arrays.deepToString(args));
        writeToSummaryLogFile("Raw Sequence Files: " + Arrays.deepToString(rawSeqFileNames));
        writeToSummaryLogFile("Parameters Passed to ProductionSNPCallerPlugin: " + Arrays.deepToString(args));

        ProductionSNPCallerPlugin plugin = new ProductionSNPCallerPlugin();

        plugin.setParameters(args);

        printParameterValues();
        plugin.performFunction(null);

        writeToSummaryLogFile("Production Pipeline Finished: " + subDirectory);
        return true;
    }

    private String[] getPluginArgs(String inputDir, String keyFile) {
        String[] args = {
            "-i", inputDir,
            "-k", keyFile,
            "-e", enzyme(),
            "-o", outputGenotypeFile(),
            "-m", productionTOPM(),
            "-ko"
        };
        return args;
    }

    private void setupLogfile() {

        String todayDate = new SimpleDateFormat("yyyyMMdd").format(new Date());
        String logFileName = todayDate + "_" + "ProductionPipeline" + ".log";
        logFileName = myOutputDirectory + "/" + logFileName;

        myLogger.info("Log File: " + logFileName);
        try {
            LoggingUtils.setupLogfile(logFileName);
        } catch (Exception e) {
            throw new IllegalArgumentException("ProductionPipeline: setupLogfile: " + logFileName + " doesn't exist.");
        }

    }

    private void writeToSummaryLogFile(String str) {
        try {
            mySummaryLogFile.write(getTimeStamp());
            mySummaryLogFile.write(str);
            mySummaryLogFile.write("\n");
        } catch (Exception e) {
            myLogger.error("writeToSummaryLogFile: Problem writing to Summary Log File.");
        }
    }

    /**
     * Convenience method to provide uniformly labeled timestamps
     */
    private static String getTimeStamp() {
        return "Timestamp: " + LOGGING_DATE_FORMAT.format(new Date()) + ": ";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Production Pipeline";
    }

    @Override
    public String getToolTipText() {
        return "Production Pipeline";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(ProductionPipeline.class);
    // }
    /**
     * Input directory containing subdirectories with fastq AND/OR qseq files
     *
     * @return Input Directory
     */
    public String inputDirectory() {
        return myInputDirectory.value();
    }

    /**
     * Set Input Directory. Input directory containing subdirectories with fastq
     * AND/OR qseq files
     *
     * @param value Input Directory
     *
     * @return this plugin
     */
    public ProductionPipeline inputDirectory(String value) {
        myInputDirectory = new PluginParameter<>(myInputDirectory, value);
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
    public ProductionPipeline enzyme(String value) {
        myEnzyme = new PluginParameter<>(myEnzyme, value);
        return this;
    }

    /**
     * Physical map file containing tags and corresponding variants (production
     * TOPM)
     *
     * @return Production T O P M
     */
    public String productionTOPM() {
        return myProductionTOPM.value();
    }

    /**
     * Set Production T O P M. Physical map file containing tags and
     * corresponding variants (production TOPM)
     *
     * @param value Production T O P M
     *
     * @return this plugin
     */
    public ProductionPipeline productionTOPM(String value) {
        myProductionTOPM = new PluginParameter<>(myProductionTOPM, value);
        return this;
    }

    /**
     * Output (target) HDF5 genotypes file to add new genotypes to (new file
     * created if it doesn't exist)
     *
     * @return Output Genotype File
     */
    public String outputGenotypeFile() {
        return myOutputGenotypeFile.value();
    }

    /**
     * Set Output Genotype File. Output (target) HDF5 genotypes file to add new
     * genotypes to (new file created if it doesn't exist)
     *
     * @param value Output Genotype File
     *
     * @return this plugin
     */
    public ProductionPipeline outputGenotypeFile(String value) {
        myOutputGenotypeFile = new PluginParameter<>(myOutputGenotypeFile, value);
        return this;
    }

    /**
     * Archive directory where to move processed files
     *
     * @return Archive Directory
     */
    public String archiveDirectory() {
        return myArchiveDirectory.value();
    }

    /**
     * Set Archive Directory. Archive directory where to move processed files
     *
     * @param value Archive Directory
     *
     * @return this plugin
     */
    public ProductionPipeline archiveDirectory(String value) {
        myArchiveDirectory = new PluginParameter<>(myArchiveDirectory, value);
        return this;
    }
}
