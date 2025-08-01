/*
 *  VCAPScanPlugin
 * 
 *  Created on Apr 4, 2016
 */
package net.maizegenetics.analysis.distance;

import com.google.common.collect.Range;
import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.TimeUnit;
import javax.swing.ImageIcon;
import net.maizegenetics.analysis.data.ExportPlugin;
import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.FilterSite;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.io.JSONUtils;
import net.maizegenetics.dna.snp.io.ReadBedfile;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.DefaultPluginListener;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.plugindef.ThreadedPluginListener;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class VCAPScanPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(VCAPScanPlugin.class);

    public static enum SCAN_METHOD {

        Chromosome,
        Site_Blocks,
        Directory_Of_Files,
        Bed_File
    };

    private PluginParameter<SCAN_METHOD> myMethod = new PluginParameter.Builder<>("method", SCAN_METHOD.Chromosome, SCAN_METHOD.class)
            .guiName("Scan method")
            .range(SCAN_METHOD.values())
            .description("")
            .build();

    private PluginParameter<Integer> myBlockingWindowSize = new PluginParameter.Builder<>("blockingWindowSize", 0, Integer.class)
            .description("Blocking window size. Number of sites (Site_Blocks) or physical positions (Bed_File) on each side of regions.")
            .range(Range.atLeast(0))
            .build();

    private PluginParameter<Integer> myNumSitesPerBlock = new PluginParameter.Builder<>("numSitesPerBlock", 10000, Integer.class)
            .description("For Site_Blocks method, this sets number of sites per block. Blocks do not span chromosomes.")
            .range(Range.atLeast(1))
            .build();

    private PluginParameter<Integer> myStepSize = new PluginParameter.Builder<>("stepSize", 0, Integer.class)
            .description("Step Size. If step size set to default 0, the step size will equal number sites per block.")
            .range(Range.atLeast(0))
            .build();

    private PluginParameter<String> myDirOfFiles = new PluginParameter.Builder<>("dirOfFiles", null, String.class)
            .description("Directory contains files for sub-matrices. Can be text (.txt) containing chromosome / positions, Bed (.bed), or Position List (.json or .json.gz)")
            .inDir()
            .build();

    private PluginParameter<String> myBedFile = new PluginParameter.Builder<>("bedFile", null, String.class)
            .description("For Bed_File method, this specifies bed file to use.  Each line / range in the bed file will be a block in the scan.")
            .inFile()
            .build();

    private PluginParameter<String> myBlockingBedFile = new PluginParameter.Builder<>("blockingBedFile", null, String.class)
            .description("For Bed_File method, this specifies optional bed file to define blocking window for each block in the main bed file.")
            .inFile()
            .build();

    private PluginParameter<String> myPhenotypeFile = new PluginParameter.Builder<>("phenotypeFile", null, String.class)
            .description("The phenotype file to use with LDAK.  Must have only one phenotype.")
            .inFile()
            .required(true)
            .build();

    private PluginParameter<String> myLDAKCommand = new PluginParameter.Builder<>("ldakCommand", "ldak.4.9.fast", String.class)
            .description("Command to call LDAK.  If not on PATH, full pathname needed.")
            .build();

    private PluginParameter<String> myOutputDir = new PluginParameter.Builder<>("outputDir", ".", String.class)
            .description("Directory to output kinship matrices and LDAK reml results.")
            .outDir()
            .build();

    private PluginParameter<String> myWholeMatrix = new PluginParameter.Builder<>("wholeMatrix", null, String.class)
            .description("Kinship matrix of whole genotype dataset.  This will be created if not specified.  Must have been exported from Tassel with export type SqrMatrixBin.")
            .inFile()
            .build();

    public VCAPScanPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
        if ((alignInList == null) || (alignInList.size() != 1)) {
            throw new IllegalArgumentException("VCAPScanPlugin: Must input exactly one genotype table.");
        }
    }

    @Override
    protected void postProcessParameters() {
        if (outputDir().charAt(outputDir().length() - 1) != '/') {
            outputDir(outputDir() + "/");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        GenotypeTable genotype = (GenotypeTable) input.getDataOfType(GenotypeTable.class).get(0).getData();
        List<String> matrixFiles = new ArrayList<>();

        if (method() == SCAN_METHOD.Chromosome) {

            ForkJoinPool threadPool = new ForkJoinPool();

            KinshipPlugin kinshipPlugin = new KinshipPlugin(null, false)
                    .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS);
            kinshipPlugin.addListener(DefaultPluginListener.getInstance());

            SubtractDistanceMatrixPlugin subtractPlugin = new SubtractDistanceMatrixPlugin(null, false)
                    .wholeMatrix(createWholeMatrixIfNeeded(input));

            for (Chromosome currentChr : genotype.chromosomes()) {

                String chrStr = currentChr.getName();
                FilterSiteBuilderPlugin filterChr = new FilterSiteBuilderPlugin(null, false)
                        .startChr(currentChr)
                        .endChr(currentChr);
                DataSet genotypeDataSet = filterChr.performFunction(input);
                GenotypeTable genotypeChr = (GenotypeTable) genotypeDataSet.getData(0).getData();
                int numSites = genotypeChr.numberOfSites();

                DataSet part = kinshipPlugin.performFunction(genotypeDataSet);

                DataSet rest = subtractPlugin.performFunction(part);

                ExportPlugin exportPlugin = new ExportPlugin(null, false);
                exportPlugin.fileType(FileLoadPlugin.TasselFileType.SqrMatrixBin);
                String startPosStr = String.format("%012d", genotypeChr.chromosomalPosition(0));
                String endPosStr = String.format("%012d", genotypeChr.chromosomalPosition(numSites - 1));
                String saveFilename = outputDir() + "Kinship_" + chrStr + "_" + startPosStr + "_" + endPosStr;
                matrixFiles.add(saveFilename);
                exportPlugin.saveFile(saveFilename);
                threadPool.submit(new ThreadedPluginListener(exportPlugin, new PluginEvent(new DataSet(part.getData(0), part.getCreator()))));

                ExportPlugin exportPlugin1 = new ExportPlugin(null, false);
                exportPlugin1.fileType(FileLoadPlugin.TasselFileType.SqrMatrixBin);
                exportPlugin1.saveFile(saveFilename + "Rest");
                threadPool.submit(new ThreadedPluginListener(exportPlugin1, new PluginEvent(new DataSet(rest.getData(0), rest.getCreator()))));

            }

            threadPool.shutdown();
            try {
                threadPool.awaitTermination(20, TimeUnit.MINUTES);
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("VCAPScanPlugin: processData: problem: " + e.getMessage());
            }

        } else if (method() == SCAN_METHOD.Site_Blocks) {

            ForkJoinPool threadPool = new ForkJoinPool();

            int numSitesPerBlock = numSitesPerBlock();

            KinshipPlugin kinshipPlugin = new KinshipPlugin(null, false)
                    .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS);
            kinshipPlugin.addListener(DefaultPluginListener.getInstance());

            SubtractDistanceMatrixPlugin subtractPlugin = new SubtractDistanceMatrixPlugin(null, false)
                    .wholeMatrix(createWholeMatrixIfNeeded(input));

            Chromosome[] chromosomes = genotype.chromosomes();

            for (Chromosome c : chromosomes) {

                String chrStr = c.getName();
                FilterSiteBuilderPlugin filterChr = new FilterSiteBuilderPlugin(null, false)
                        .startChr(c)
                        .endChr(c);
                DataSet genotypeDataSet = filterChr.performFunction(input);
                GenotypeTable genotypeChr = (GenotypeTable) genotypeDataSet.getData(0).getData();
                int numSites = genotypeChr.numberOfSites();

                int winBufferSize = blockingWindowSize();
                int stepSize = stepSize();
                if (stepSize < 1) {
                    stepSize = numSitesPerBlock;
                }

                for (int startSite = 0; startSite < numSites; startSite += stepSize) {

                    int endSite = Math.min(startSite + numSitesPerBlock - 1, numSites - 1);

                    int startSiteBlocking = Math.max(startSite - winBufferSize, 0);
                    int endSiteBlocking = Math.min(endSite + winBufferSize, numSites - 1);

                    String startPosStr = String.format("%012d", genotypeChr.chromosomalPosition(startSite));
                    String endPosStr = String.format("%012d", genotypeChr.chromosomalPosition(endSite));
                    String saveFilename = outputDir() + "Kinship_" + chrStr + "_" + startPosStr + "_" + endPosStr;
                    matrixFiles.add(saveFilename);
                    if (new File(saveFilename + "Rest.grm.bin").isFile() && new File(saveFilename + ".grm.bin").isFile()) {
                        myLogger.info(saveFilename + " already exists");
                        continue;
                    }

                    FilterSiteBuilderPlugin filter = new FilterSiteBuilderPlugin(null, false)
                            .startSite(startSite)
                            .endSite(endSite);
                    DataSet filteredGenotype = filter.performFunction(genotypeDataSet);

                    DataSet part = kinshipPlugin.performFunction(filteredGenotype);

                    DataSet blocking;
                    if (startSite != startSiteBlocking || endSite != endSiteBlocking) {
                        FilterSiteBuilderPlugin filterBlocking = new FilterSiteBuilderPlugin(null, false)
                                .startSite(startSiteBlocking)
                                .endSite(endSiteBlocking);
                        DataSet filteredBlocking = filterBlocking.performFunction(genotypeDataSet);

                        blocking = kinshipPlugin.performFunction(filteredBlocking);
                    } else {
                        blocking = part;
                    }
                    DataSet rest = subtractPlugin.performFunction(blocking);

                    ExportPlugin exportPlugin = new ExportPlugin(null, false);
                    exportPlugin.fileType(FileLoadPlugin.TasselFileType.SqrMatrixBin);
                    exportPlugin.saveFile(saveFilename);
                    threadPool.submit(new ThreadedPluginListener(exportPlugin, new PluginEvent(new DataSet(part.getData(0), part.getCreator()))));

                    ExportPlugin exportPlugin1 = new ExportPlugin(null, false);
                    exportPlugin1.fileType(FileLoadPlugin.TasselFileType.SqrMatrixBin);
                    exportPlugin1.saveFile(saveFilename + "Rest");
                    threadPool.submit(new ThreadedPluginListener(exportPlugin1, new PluginEvent(new DataSet(rest.getData(0), rest.getCreator()))));

                }

            }

            threadPool.shutdown();
            try {
                threadPool.awaitTermination(20, TimeUnit.MINUTES);
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("VCAPScanPlugin: processData: problem: " + e.getMessage());
            }

        } else if (method() == SCAN_METHOD.Directory_Of_Files) {

            ForkJoinPool threadPool = new ForkJoinPool();

            KinshipPlugin kinshipPlugin = new KinshipPlugin(null, false)
                    .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS);
            kinshipPlugin.addListener(DefaultPluginListener.getInstance());

            SubtractDistanceMatrixPlugin subtractPlugin = new SubtractDistanceMatrixPlugin(null, false)
                    .wholeMatrix(createWholeMatrixIfNeeded(input));

            String[] files = new File(dirOfFiles()).list();
            for (String current : files) {

                String file = dirOfFiles() + "/" + current;

                DataSet genotypeDataSet = null;
                if (file.endsWith(".txt")) {
                    FilterSiteBuilderPlugin filterChr = new FilterSiteBuilderPlugin(null, false)
                            .chrPosFile(file);
                    genotypeDataSet = filterChr.performFunction(input);
                } else if (file.endsWith(".bed")) {
                    FilterSiteBuilderPlugin filterChr = new FilterSiteBuilderPlugin(null, false)
                            .bedFile(file);
                    genotypeDataSet = filterChr.performFunction(input);
                } else if (file.endsWith(".json") || file.endsWith(".json.gz")) {
                    FilterSiteBuilderPlugin filterChr = new FilterSiteBuilderPlugin(null, false)
                            .positionList(JSONUtils.importPositionListFromJSON(file));
                    genotypeDataSet = filterChr.performFunction(input);
                } else {
                    continue;
                }

                DataSet part = kinshipPlugin.performFunction(genotypeDataSet);
                DataSet rest = subtractPlugin.performFunction(part);
                GenotypeTable genotypeChr = (GenotypeTable) genotypeDataSet.getData(0).getData();
                int numSites = genotypeChr.numberOfSites();
                ExportPlugin exportPlugin = new ExportPlugin(null, false);
                exportPlugin.fileType(FileLoadPlugin.TasselFileType.SqrMatrixBin);
                String startPosStr = String.format("%012d", genotypeChr.chromosomalPosition(0));
                String endPosStr = String.format("%012d", genotypeChr.chromosomalPosition(numSites - 1));
                String saveFilename = outputDir() + "Kinship_0_" + startPosStr + "_" + endPosStr + "_" + Utils.getFilename(file);
                matrixFiles.add(saveFilename);
                exportPlugin.saveFile(saveFilename);
                threadPool.submit(new ThreadedPluginListener(exportPlugin, new PluginEvent(new DataSet(part.getData(0), part.getCreator()))));

                ExportPlugin exportPlugin1 = new ExportPlugin(null, false);
                exportPlugin1.fileType(FileLoadPlugin.TasselFileType.SqrMatrixBin);
                exportPlugin1.saveFile(saveFilename + "Rest");
                threadPool.submit(new ThreadedPluginListener(exportPlugin1, new PluginEvent(new DataSet(rest.getData(0), rest.getCreator()))));

            }

            threadPool.shutdown();
            try {
                threadPool.awaitTermination(20, TimeUnit.MINUTES);
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("VCAPScanPlugin: processData: problem: " + e.getMessage());
            }

        } else if (method() == SCAN_METHOD.Bed_File) {

            if (bedFile() == null || bedFile().isEmpty()) {
                throw new IllegalArgumentException("VCAPScanPlugin: processData: bed file must be specified for method Bed_File.");
            } else if (!new File(bedFile()).isFile()) {
                throw new IllegalArgumentException("VCAPScanPlugin: processData: bed file doesn't exist: " + bedFile());
            }

            ForkJoinPool threadPool = new ForkJoinPool();

            KinshipPlugin kinshipPlugin = new KinshipPlugin(null, false)
                    .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS);
            kinshipPlugin.addListener(DefaultPluginListener.getInstance());

            SubtractDistanceMatrixPlugin subtractPlugin = new SubtractDistanceMatrixPlugin(null, false)
                    .wholeMatrix(createWholeMatrixIfNeeded(input));

            List<ReadBedfile.BedFileRange> ranges = ReadBedfile.getRanges(bedFile());
            int numRanges = ranges.size();

            List<ReadBedfile.BedFileRange> blockingRanges = null;
            if (blockingBedFile() == null || blockingBedFile().isEmpty()) {
                myLogger.info("processData: no blocking windows used with bed file: " + bedFile());
            } else if (!new File(blockingBedFile()).isFile()) {
                throw new IllegalArgumentException("VCAPScanPlugin: processData: blocking bed file doesn't exist: " + blockingBedFile());
            } else {
                blockingRanges = ReadBedfile.getRanges(blockingBedFile());
                if (ranges.size() != blockingRanges.size()) {
                    throw new IllegalArgumentException("VCAPScanPlugin: processData: must be same number of ranges in bed file and blocking bed file.");
                }
                for (int i = 0; i < numRanges; i++) {
                    ReadBedfile.BedFileRange range = ranges.get(i);
                    ReadBedfile.BedFileRange blockingRange = blockingRanges.get(i);
                    if (!blockingRange.chr().equals(range.chr())) {
                        throw new IllegalArgumentException("VCAPScanPlugin: processData: block range chr: " + blockingRange.chr() + " should equal range chr: " + range.chr());
                    }
                    if (blockingRange.start() > range.start()) {
                        throw new IllegalArgumentException("VCAPScanPlugin: processData: blocking range start: " + blockingRange.start() + " should be less than or equal to range start: " + range.start());
                    }
                    if (blockingRange.end() < range.end()) {
                        throw new IllegalArgumentException("VCAPScanPlugin: processData: blocking range end: " + blockingRange.end() + " should be greater than or equal to range end: " + range.end());
                    }
                }
            }

            for (int i = 0; i < numRanges; i++) {

                ReadBedfile.BedFileRange range = ranges.get(i);

                int startSite = genotype.siteOfPhysicalPosition(range.start(), new Chromosome(range.chr()));
                if (startSite < 0) {
                    startSite = -startSite - 1;
                }

                int endSite = genotype.siteOfPhysicalPosition(range.end(), new Chromosome(range.chr()));
                if (endSite < 0) { // end position doesn't exist, so already excluded
                    endSite = -endSite - 2;
                } else { // end position is exclusive
                    endSite--;
                }

                if (endSite <= startSite) {
                    myLogger.warn("No sites in region chr: " + range.chr() + " start: " + range.start() + " end: " + range.end());
                    continue;
                }

                int startSiteBlocking = startSite;
                int endSiteBlocking = endSite;
                if (blockingRanges != null) {
                    ReadBedfile.BedFileRange blockingRange = blockingRanges.get(i);
                    startSiteBlocking = genotype.siteOfPhysicalPosition(blockingRange.start(), new Chromosome(blockingRange.chr()));
                    if (startSiteBlocking < 0) {
                        startSiteBlocking = -startSiteBlocking - 1;
                    }

                    endSiteBlocking = genotype.siteOfPhysicalPosition(blockingRange.end(), new Chromosome(blockingRange.chr()));
                    if (endSiteBlocking < 0) { // end position doesn't exist, so already excluded
                        endSiteBlocking = -endSiteBlocking - 2;
                    } else { // end position is exclusive
                        endSiteBlocking--;
                    }
                } else if (blockingWindowSize() != 0) {
                    startSiteBlocking = genotype.siteOfPhysicalPosition(Math.max(0, range.start() - blockingWindowSize()), new Chromosome(range.chr()));
                    if (startSiteBlocking < 0) {
                        startSiteBlocking = -startSiteBlocking - 1;
                    }

                    endSiteBlocking = genotype.siteOfPhysicalPosition(range.end() + blockingWindowSize(), new Chromosome(range.chr()));
                    if (endSiteBlocking < 0) { // end position doesn't exist, so already excluded
                        endSiteBlocking = -endSiteBlocking - 2;
                    } else { // end position is exclusive
                        endSiteBlocking--;
                    }
                }

                String startPosStr = String.format("%012d", genotype.chromosomalPosition(startSite));
                String endPosStr = String.format("%012d", genotype.chromosomalPosition(endSite));
                String saveFilename = outputDir() + "Kinship_" + range.chr() + "_" + startPosStr + "_" + endPosStr;
                if (range.name() != null) {
                    saveFilename += "_" + range.name();
                }
                matrixFiles.add(saveFilename);
                if (new File(saveFilename + "Rest.grm.bin").isFile() && new File(saveFilename + ".grm.bin").isFile()) {
                    myLogger.info(saveFilename + " already exists");
                    continue;
                }

                FilterSiteBuilderPlugin filter = new FilterSiteBuilderPlugin(null, false)
                        .siteFilter(FilterSite.SITE_RANGE_FILTER_TYPES.SITES)
                        .startSite(startSite)
                        .endSite(endSite);
                DataSet filteredGenotype = filter.performFunction(input);

                DataSet part = kinshipPlugin.performFunction(filteredGenotype);

                DataSet blocking;
                if (startSite != startSiteBlocking || endSite != endSiteBlocking) {
                    FilterSiteBuilderPlugin filterBlocking = new FilterSiteBuilderPlugin(null, false)
                            .startSite(startSiteBlocking)
                            .endSite(endSiteBlocking);
                    DataSet filteredBlocking = filterBlocking.performFunction(input);

                    blocking = kinshipPlugin.performFunction(filteredBlocking);
                } else {
                    blocking = part;
                }
                DataSet rest = subtractPlugin.performFunction(blocking);

                ExportPlugin exportPlugin = new ExportPlugin(null, false);
                exportPlugin.fileType(FileLoadPlugin.TasselFileType.SqrMatrixBin);
                exportPlugin.saveFile(saveFilename);
                threadPool.submit(new ThreadedPluginListener(exportPlugin, new PluginEvent(new DataSet(part.getData(0), part.getCreator()))));

                ExportPlugin exportPlugin1 = new ExportPlugin(null, false);
                exportPlugin1.fileType(FileLoadPlugin.TasselFileType.SqrMatrixBin);
                exportPlugin1.saveFile(saveFilename + "Rest");
                threadPool.submit(new ThreadedPluginListener(exportPlugin1, new PluginEvent(new DataSet(rest.getData(0), rest.getCreator()))));

            }

            threadPool.shutdown();
            try {
                threadPool.awaitTermination(20, TimeUnit.MINUTES);
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("VCAPScanPlugin: processData: problem: " + e.getMessage());
            }

        }

        runLDAK(matrixFiles);
        getResults();

        return null;

    }

    private String createWholeMatrixIfNeeded(DataSet input) {
        if ((wholeMatrix() == null) || (wholeMatrix().isEmpty())) {
            String saveFilename = outputDir() + "kinship_whole.txt";
            if (!isFileEmpty(saveFilename)) {
                myLogger.info("Whole Kinship already exists: " + saveFilename);
                return saveFilename;
            }

            KinshipPlugin kinshipPlugin = new KinshipPlugin(null, false)
                    .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS);
            kinshipPlugin.addListener(DefaultPluginListener.getInstance());
            DataSet whole = kinshipPlugin.performFunction(input);

            ExportPlugin exportPlugin = new ExportPlugin(null, false);
            exportPlugin.fileType(FileLoadPlugin.TasselFileType.SqrMatrix);
            exportPlugin.saveFile(saveFilename);
            exportPlugin.performFunction(whole);
            return saveFilename;
        } else {
            return wholeMatrix();
        }
    }

    private void runLDAK(List<String> matrixFiles) {

        // ldak.4.9.fast --reml results_chr1 --mgrm kinship_list.txt --pheno NAM_ap_dta_multiblup.txt --kinship-details NO
        int numProcesses = 0;
        Process[] runningProcesses = new Process[3];
        BufferedReader[] readers = new BufferedReader[3];
        String[] commands = new String[3];
        for (String filename : matrixFiles) {

            String resultFilename = outputDir() + "Results" + Utils.getFilename(filename);
            if (!isFileEmpty(resultFilename + ".reml")) {
                continue;
            }

            String listFilename = "kinship_list" + numProcesses + ".txt";
            try (BufferedWriter writer = Utils.getBufferedWriter(listFilename)) {
                writer.write(filename);
                writer.write("\n");
                writer.write(filename + "Rest");
                writer.write("\n");
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
            }

            commands[numProcesses] = ldakCommand()
                    + " --reml " + resultFilename
                    + " --mgrm " + listFilename
                    + " --pheno " + phenotypeFile()
                    + " --kinship-details NO";

            try {
                runningProcesses[numProcesses] = Runtime.getRuntime().exec(commands[numProcesses]);
                readers[numProcesses] = new BufferedReader(
                        new InputStreamReader(runningProcesses[numProcesses].getInputStream()));
            } catch (Exception e) {
                myLogger.error(e.getMessage(), e);
            }

            numProcesses++;
            if (numProcesses == 3) {
                for (int i = 0; i < numProcesses; i++) {
                    try {
                        runningProcesses[i].waitFor();
                        myLogger.info("command: " + commands[i]);
                        String line = null;
                        while ((line = readers[i].readLine()) != null) {
                            System.out.println(line);
                        }
                        readers[i].close();
                    } catch (Exception e) {
                        myLogger.error(e.getMessage(), e);
                    }
                }
                numProcesses = 0;
            }

        }

        for (int i = 0; i < numProcesses; i++) {
            try {
                runningProcesses[i].waitFor();
                myLogger.info("command: " + commands[i]);
                String line = null;
                while ((line = readers[i].readLine()) != null) {
                    System.out.println(line);
                }
                readers[i].close();
            } catch (Exception e) {
                myLogger.error(e.getMessage(), e);
            }
        }

    }

    private void getResults() {

        TreeMap<ChrPos, String[]> rows = new TreeMap<>();
        String output = outputDir() + Utils.getFilename(phenotypeFile()) + "Results.txt";

        String[] files = new File(outputDir()).list();
        for (String file : files) {
            if (file.endsWith(".reml")) {
                String[] result = new String[6];
                String[] tokens = file.substring(0, file.indexOf(".reml")).split("_");
                result[0] = tokens[1]; // chr
                result[1] = tokens[2]; // start pos
                result[2] = tokens[3]; // end pos
                if (tokens.length > 4) {
                    result[3] = tokens[4]; // comment
                } else {
                    result[3] = "";
                }
                String file3 = outputDir() + file;
                try (BufferedReader reader = Utils.getBufferedReader(file3)) {
                    String line = reader.readLine();
                    while (line != null) {
                        if (line.startsWith("Her_K1")) {
                            String[] temp = line.split(" ");
                            result[4] = temp[1];
                        } else if (line.startsWith("Her_K2")) {
                            String[] temp = line.split(" ");
                            result[5] = temp[1];
                        }
                        line = reader.readLine();
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
                rows.put(new ChrPos(Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2])), result);

            }
        }

        myLogger.info("Writing Results file: " + output);
        try (BufferedWriter writer = Utils.getBufferedWriter(output)) {

            writer.write("Chromosome\tStart Position\tEnd Position\tComment\tHeritably Subset\tHeritably Rest\n");

            for (Map.Entry<ChrPos, String[]> current : rows.entrySet()) {
                boolean first = true;
                for (int i = 0; i < current.getValue().length; i++) {
                    if (!first) {
                        writer.write("\t");
                    }
                    first = false;
                    writer.write(current.getValue()[i]);
                }
                writer.write("\n");
            }

        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    private boolean isFileEmpty(String filename) {
        if (new File(filename).isFile()) {
            try (BufferedReader reader = Utils.getBufferedReader(filename)) {
                return reader.readLine() == null;
            } catch (Exception e) {
                return true;
            }
        } else {
            return true;
        }
    }

    private class ChrPos implements Comparable<ChrPos> {

        private final int myChr;
        private final int myPos;

        public ChrPos(int chr, int pos) {
            myChr = chr;
            myPos = pos;
        }

        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof ChrPos)) {
                return false;
            }
            ChrPos other = (ChrPos) obj;
            if ((myChr == other.myChr) && (myPos == other.myPos)) {
                return true;
            } else {
                return false;
            }
        }

        @Override
        public int hashCode() {
            int hash = 5;
            hash = 37 * hash + myChr;
            hash = 37 * hash + myPos;
            return hash;
        }

        @Override
        public int compareTo(ChrPos o) {
            if (myChr < o.myChr) {
                return -1;
            } else if (myChr > o.myChr) {
                return 1;
            }
            if (myPos < o.myPos) {
                return -1;
            } else if (myPos > o.myPos) {
                return 1;
            } else {
                return 0;
            }
        }

    }

    /**
     * Scan method
     *
     * @return Scan method
     */
    public SCAN_METHOD method() {
        return myMethod.value();
    }

    /**
     * Set Scan method. Scan method
     *
     * @param value Scan method
     *
     * @return this plugin
     */
    public VCAPScanPlugin method(SCAN_METHOD value) {
        myMethod = new PluginParameter<>(myMethod, value);
        return this;
    }

    /**
     * Blocking window size
     *
     * @return Blocking Window Size
     */
    public Integer blockingWindowSize() {
        return myBlockingWindowSize.value();
    }

    /**
     * Set Blocking Window Size. Blocking window size
     *
     * @param value Blocking Window Size
     *
     * @return this plugin
     */
    public VCAPScanPlugin blockingWindowSize(Integer value) {
        myBlockingWindowSize = new PluginParameter<>(myBlockingWindowSize, value);
        return this;
    }

    /**
     * Num Sites Per Block
     *
     * @return Num Sites Per Block
     */
    public Integer numSitesPerBlock() {
        return myNumSitesPerBlock.value();
    }

    /**
     * Set Num Sites Per Block. Num Sites Per Block
     *
     * @param value Num Sites Per Block
     *
     * @return this plugin
     */
    public VCAPScanPlugin numSitesPerBlock(Integer value) {
        myNumSitesPerBlock = new PluginParameter<>(myNumSitesPerBlock, value);
        return this;
    }

    /**
     * Step Size
     *
     * @return Step Size
     */
    public Integer stepSize() {
        return myStepSize.value();
    }

    /**
     * Set Step Size. Step Size
     *
     * @param value Step Size
     *
     * @return this plugin
     */
    public VCAPScanPlugin stepSize(Integer value) {
        myStepSize = new PluginParameter<>(myStepSize, value);
        return this;
    }

    /**
     * Directory contains files for sub-matrices. Can be text (.txt) containing
     * chromosome / positions, Bed (.bed), or Position List (.json or .json.gz)
     *
     * @return Dir Of Files
     */
    public String dirOfFiles() {
        return myDirOfFiles.value();
    }

    /**
     * Set Dir Of Files. Directory contains files for sub-matrices. Can be text
     * (.txt) containing chromosome / positions, Bed (.bed), or Position List
     * (.json or .json.gz)
     *
     * @param value Dir Of Files
     *
     * @return this plugin
     */
    public VCAPScanPlugin dirOfFiles(String value) {
        myDirOfFiles = new PluginParameter<>(myDirOfFiles, value);
        return this;
    }

    /**
     * Bed File
     *
     * @return Bed File
     */
    public String bedFile() {
        return myBedFile.value();
    }

    /**
     * Set Bed File. Bed File
     *
     * @param value Bed File
     *
     * @return this plugin
     */
    public VCAPScanPlugin bedFile(String value) {
        myBedFile = new PluginParameter<>(myBedFile, value);
        return this;
    }

    /**
     * For Bed_File method, this specifies optional bed file to define blocking
     * window for each block in the main bed file.
     *
     * @return Blocking Bed File
     */
    public String blockingBedFile() {
        return myBlockingBedFile.value();
    }

    /**
     * Set Blocking Bed File. For Bed_File method, this specifies optional bed
     * file to define blocking window for each block in the main bed file.
     *
     * @param value Blocking Bed File
     *
     * @return this plugin
     */
    public VCAPScanPlugin blockingBedFile(String value) {
        myBlockingBedFile = new PluginParameter<>(myBlockingBedFile, value);
        return this;
    }

    /**
     * Phenotype File
     *
     * @return Phenotype File
     */
    public String phenotypeFile() {
        return myPhenotypeFile.value();
    }

    /**
     * Set Phenotype File. Phenotype File
     *
     * @param value Phenotype File
     *
     * @return this plugin
     */
    public VCAPScanPlugin phenotypeFile(String value) {
        myPhenotypeFile = new PluginParameter<>(myPhenotypeFile, value);
        return this;
    }

    /**
     * Ldak Command
     *
     * @return Ldak Command
     */
    public String ldakCommand() {
        return myLDAKCommand.value();
    }

    /**
     * Set Ldak Command. Ldak Command
     *
     * @param value Ldak Command
     *
     * @return this plugin
     */
    public VCAPScanPlugin ldakCommand(String value) {
        myLDAKCommand = new PluginParameter<>(myLDAKCommand, value);
        return this;
    }

    /**
     * Output Dir
     *
     * @return Output Dir
     */
    public String outputDir() {
        return myOutputDir.value();
    }

    /**
     * Set Output Dir. Output Dir
     *
     * @param value Output Dir
     *
     * @return this plugin
     */
    public VCAPScanPlugin outputDir(String value) {
        myOutputDir = new PluginParameter<>(myOutputDir, value);
        return this;
    }

    /**
     * Whole Matrix
     *
     * @return Whole Matrix
     */
    public String wholeMatrix() {
        return myWholeMatrix.value();
    }

    /**
     * Set Whole Matrix. Whole Matrix
     *
     * @param value Whole Matrix
     *
     * @return this plugin
     */
    public VCAPScanPlugin wholeMatrix(String value) {
        myWholeMatrix = new PluginParameter<>(myWholeMatrix, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = VCAPScanPlugin.class.getResource("/net/maizegenetics/analysis/images/VCAP.png");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "VCAP Scan";
    }

    @Override
    public String getToolTipText() {
        return "Variance Component Annotation Pipeline Scan";
    }

    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/VCAPScan/VCAPScan";
    }

}
