package net.maizegenetics.analysis.imputation;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.swing.ImageIcon;
import javax.swing.JFileChooser;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginListener;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Utils;

public class MergeBreakpointFilesPlugin extends AbstractPlugin {
    private static Logger myLogger = LogManager.getLogger(MergeBreakpointFilesPlugin.class);
    private static final int minSite = 0;
    private static final int maxSite = 350000000;
    static {
        LoggingUtils.setupDebugLogging();
    }
    
    private PluginParameter<Boolean> selectFiles =
            new PluginParameter.Builder<>("selectFiles", false, Boolean.class)
                    .description("Select the files to merge. Alternatively, merge all files in a directory. (Default = false)")
                    .guiName("Select Files to Merge")
                    .build();
    private PluginParameter<String> selectedFilesList =
            new PluginParameter.Builder<>("fileList", null, String.class)
                    .description("The name and path of a file containing the names of the breakpoint files to be merged.")
                    .guiName("List of files to merge")
                    .inFile()
                    .build();
    private PluginParameter<String> outputFile =
            new PluginParameter.Builder<>("outputFile", null, String.class)
                    .description("The output file. If the filename ends in gz it will be zipped. Required.")
                    .required(true)
                    .guiName("Output File")
                    .required(true)
                    .outFile()
                    .build();
    private PluginParameter<Boolean> fillToEnd = new PluginParameter.Builder<>("fillends", false, Boolean.class)
            .description("Should the first and last break point intervals in each chromosome be exended to the ends?")
            .guiName("Fill to Ends")
            .build();
            
    public MergeBreakpointFilesPlugin() {
        this(null, false);
    }

    public MergeBreakpointFilesPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        if (!selectFiles.value() && selectedFilesList.value() == null) {
            throw new IllegalArgumentException("Either a list of files must be provided or individual files selected.");
        } else if (selectFiles.value()) {
            if (!isInteractive())
                throw new IllegalArgumentException("Must supply a file list in command line mode.");
            JFileChooser myInputFileChooser = new JFileChooser();
            myInputFileChooser.setName("George");
            int result = myInputFileChooser.showOpenDialog(getParentFrame());
            File[] mySelectedFiles = new File[0];
            if (result == JFileChooser.APPROVE_OPTION) {
                mySelectedFiles = myInputFileChooser.getSelectedFiles();
            }
            if (mySelectedFiles.length > 0 || outputFile != null) {
                mergeBreakpoints(mySelectedFiles);
            } else if (mySelectedFiles.length == 0) {
                throw new IllegalArgumentException("No input files were selected.");
            } else if (outputFile == null) {
                throw new IllegalArgumentException("An output file name was not supplied.");
            }

        } else {
            try (BufferedReader br = new BufferedReader(new FileReader(selectedFilesList.value()))) {
                ArrayList<File> files = new ArrayList<>();
                String inputLine;
                while ((inputLine = br.readLine()) != null) {
                    if (inputLine.length() > 1) {
                        File nextFile = new File(inputLine);
                        if (nextFile.exists()) {
                            files.add(nextFile);
                        } else {
                            String msg = nextFile.getPath() + " does not exist.";
                            myLogger.error(msg);
                        }
                    }
                }
                
                if (files.size() > 0) {
                    File[] inputFiles = files.stream().toArray(File[]::new);
                    mergeBreakpoints(inputFiles);
                }
            } catch(IOException e) {
                String msg = "Unable to read from " + selectedFilesList.value() + ". Make certain that file exists.";
                throw new RuntimeException(msg , e);
            }
            
        }
        return null;
    }

    private void mergeBreakpoints(File[] inputFiles) {
        final Pattern white_space = Pattern.compile("\\s+");
        final Pattern colon = Pattern.compile(":");
        Map<String, Integer> allParentMap = readAndIndexParents(inputFiles);

        //Read each breakpoint file, convert the parent indexes using the allParentMap, and store in a Map keyed on Taxon name
        Map<String, List<Segment>> breakpointListMap = new HashMap<>();
        for (File input : inputFiles) {
            try (BufferedReader br = new BufferedReader(new FileReader(input))) {
                String inputLine = br.readLine();
                String[] parsedLine = white_space.split(inputLine);
                int nparents = Integer.parseInt(parsedLine[0]);
                int ntaxa = Integer.parseInt(parsedLine[1]);
                br.readLine();
                Map<Integer, Integer> parentIndexer = new HashMap<>();
                for (int i = 0; i < nparents; i++) {
                    parsedLine = white_space.split(br.readLine());
                    parentIndexer.put(Integer.decode(parsedLine[0]), allParentMap.get(parsedLine[1]));
                }
                br.readLine(); //skip two lines
                br.readLine();

                //for each taxon, add each segment to that taxon list in breakpointListMap
                for (int i = 0; i < ntaxa; i++) {
                    
                    try{
                        parsedLine = white_space.split(br.readLine());
                    } catch(Exception e) {
                        myLogger.debug(String.format("null pointer at %d in %s.", i, input.getName()));
                        System.out.println();
                    }
                    String taxon = parsedLine[0];
                    List<Segment> segmentList = Arrays.stream(parsedLine).skip(1)
                            .map(seg -> Segment.getInstance(colon.split(seg)))
                            .filter(opt -> opt.isPresent())
                            .map(opt -> {
                                Segment seg = opt.get();
                                seg.parent1 = parentIndexer.get(seg.parent1);
                                seg.parent2 = parentIndexer.get(seg.parent2);
                                return seg;
                            })
                            .collect(Collectors.toList());

                    List<Segment> bplist = breakpointListMap.get(taxon);
                    if (bplist == null) {
                        bplist = segmentList;
                        breakpointListMap.put(taxon, bplist);
                    } else {
                        bplist.addAll(segmentList);
                    }
                }

            } catch (IOException e) {
                throw new RuntimeException("Failed to read " + input.getName(), e);
            }
        }

        List<String> taxa = new ArrayList<>(breakpointListMap.keySet());
        Collections.sort(taxa);

        try (BufferedWriter bw = Utils.getBufferedWriter(outputFile.value())) {
            //write: #parents #taxa
            int nparents = allParentMap.size();
            bw.write(String.format("%d\t%d\n", nparents, taxa.size()));
            
            //write first comment
            bw.write(WritePopulationAlignmentPlugin.brkptComment1);
            
            //write parent indices
            String[] parents = new String[nparents];
            for (Map.Entry<String, Integer> parentIndex : allParentMap.entrySet()) {
                parents[parentIndex.getValue()] = parentIndex.getKey();
            }
            for (int i = 0; i < nparents; i++) {
                bw.write(String.format("%d\t%s\n", i, parents[i]));
            }
            
            //write second and third comments
            bw.write(WritePopulationAlignmentPlugin.brkptComment2);
            bw.write(WritePopulationAlignmentPlugin.brkptComment3);
            
            //sort the segments in each breakpoint list
            //output the sorted segments to the bp file
            for (String taxon : taxa) {
                List<Segment> segList = breakpointListMap.get(taxon);
                Collections.sort(segList);
                if (fillToEnd.value()) fillToChromosomeEnds(segList);
                String outline =
                        segList.stream().map(seg -> seg.toString())
                        .collect(Collectors.joining("\t", taxon + "\t", "\n"));
                bw.write(outline);
            }

        } catch (IOException e) {
            throw new RuntimeException("Unable to write to " + outputFile.value(), e);
        }

    }

    private Map<String, Integer> readAndIndexParents(File[] inputFiles) {
        //go through the files and read the parents
        final Pattern white_space = Pattern.compile("\\s+");
        TreeSet<String> parentSet = new TreeSet<>();
        for (File input : inputFiles) {
            try (BufferedReader br = new BufferedReader(new FileReader(input))) {
                //skip two lines
                br.readLine();
                br.readLine();
                String inputLine = br.readLine();
                while (!inputLine.startsWith("#")) {
                    if (inputLine.length() > 0) {
                        String[] parsedLine = white_space.split(inputLine);
                        if (parsedLine.length == 2)
                            parentSet.add(parsedLine[1]);
                    }
                    inputLine = br.readLine();
                }
            } catch (IOException e) {
                throw new RuntimeException("Failed to read " + input.getName(), e);
            }
        }

        Map<String, Integer> parentMap = new HashMap<>();
        int parentCount = 0;
        for (String parent : parentSet)
            parentMap.put(parent, parentCount++);
        return parentMap;
    }

    private void fillToChromosomeEnds(List<Segment> segList) {
        //assume the list has been sorted
        int nseg = segList.size();
        Segment first = segList.get(0);
        first.start = minSite;
        segList.get(nseg - 1).end = maxSite;
        String prevChr = "none";
        Segment prevSeg = first;
        for (int s = 1; s < nseg; s++) {
            Segment currSeg = segList.get(s);
            if (!currSeg.chr.equals(prevSeg.chr)) {
                currSeg.start = minSite;
                prevSeg.end = maxSite;
            }
            prevSeg = currSeg;
        }
    }
    
    @Override
    public ImageIcon getIcon() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getButtonName() {
        return "Merge breakpoint files";
    }

    @Override
    public String getToolTipText() {
        return "Merge selected breakpoint files.";
    }

    static class Segment implements Comparable<Segment> {
        static final Pattern colon = Pattern.compile(":");
        String chr;
        int chrnumber;
        int start;
        int end;
        int parent1;
        int parent2;

        Segment() {
        }

        Segment(String[] input) {
            chr = input[0];
            try {
                chrnumber = Integer.parseInt(chr);
            } catch (Exception e) {
                chrnumber = -1;
            }
            start = Integer.parseInt(input[1]);
            end = Integer.parseInt(input[2]);
            parent1 = Integer.parseInt(input[3]);
            parent2 = Integer.parseInt(input[4]);
        }

        public static Optional<Segment> getInstance(String[] input) {
            if (input.length != 5)
                return Optional.empty();
            return Optional.of(new Segment(input));
        }

        public String toString() {
            return String.format("%s:%d:%d:%d:%d", chr, start, end, parent1, parent2);
        }

        @Override
        public int compareTo(Segment seg) {
            if (chr != seg.chr) {
                if (chrnumber == seg.chrnumber)
                    return chr.compareTo(seg.chr);
                return chrnumber - chrnumber;
            }
            return start - seg.start;
        }
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(MergeBreakpointFilesPlugin.class);
    // }

    /**
     * Select the files to merge. Alternatively, merge all
     * files in a directory. (Default = false)
     *
     * @return Select Files to Merge
     */
    public Boolean selectFiles() {
        return selectFiles.value();
    }

    /**
     * Set Select Files to Merge. Select the files to merge.
     * Alternatively, merge all files in a directory. (Default
     * = false)
     *
     * @param value Select Files to Merge
     *
     * @return this plugin
     */
    public MergeBreakpointFilesPlugin selectFiles(Boolean value) {
        selectFiles = new PluginParameter<>(selectFiles, value);
        return this;
    }

    /**
     * The name and path of a file containing the names of
     * the breakpoint files to be merged.
     *
     * @return List of files to merge
     */
    public String selectedFilesList() {
        return selectedFilesList.value();
    }

    /**
     * Set List of files to merge. The name and path of a
     * file containing the names of the breakpoint files to
     * be merged.
     *
     * @param value List of files to merge
     *
     * @return this plugin
     */
    public MergeBreakpointFilesPlugin selectedFilesList(String value) {
        selectedFilesList = new PluginParameter<>(selectedFilesList, value);
        return this;
    }

    /**
     * The output file. If the filename ends in gz it will
     * be zipped. Required.
     *
     * @return Output File
     */
    public String outputFile() {
        return outputFile.value();
    }

    /**
     * Set Output File. The output file. If the filename ends
     * in gz it will be zipped. Required.
     *
     * @param value Output File
     *
     * @return this plugin
     */
    public MergeBreakpointFilesPlugin outputFile(String value) {
        outputFile = new PluginParameter<>(outputFile, value);
        return this;
    }

}
