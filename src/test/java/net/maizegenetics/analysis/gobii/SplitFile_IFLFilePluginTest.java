/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.util.LoggingUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.Assert.*;

/**
 * Splits a GOBII marker intermediate file built from the markers in a Tutorial HapMap file.
 *
 * @author lcj34
 *
 */
public class SplitFile_IFLFilePluginTest {

    private static final Path TEMP_DIR = Paths.get(GeneralConstants.TEMP_DIR, "SplitFile_IFLFilePluginTest");
    private static final String HEADER = "name\tchr\tposition";

    @Test
    public void testSplitFile_IFLFilePlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        Path outDir = cleanDir("split");
        List<String> rows = markerRows();
        Path infile = writeMarkerFile("DS_4.marker", rows);
        int maximumSize = 100;

        new SplitFile_IFLFilePlugin()
                .inFile(infile.toString())
                .outputDir(outDir.toString() + "/")
                .maxSize(maximumSize)
                .performFunction(null);

        int expectedFiles = (rows.size() + maximumSize - 1) / maximumSize;
        assertEquals("Test data should need more than one file", 3, expectedFiles);
        try (Stream<Path> files = Files.list(outDir)) {
            assertEquals(expectedFiles, files.count());
        }
        for (int fileNum = 1; fileNum <= expectedFiles; fileNum++) {
            Path split = outDir.resolve("DS_4_" + fileNum + ".marker");
            assertTrue("Missing " + split, Files.exists(split));
            List<String> lines = Files.readAllLines(split);
            assertEquals("Every file starts with the header", HEADER, lines.get(0));
            int from = (fileNum - 1) * maximumSize;
            int to = Math.min(rows.size(), from + maximumSize);
            assertEquals("Rows of file " + fileNum, rows.subList(from, to), lines.subList(1, lines.size()));
        }
    }

    @Test
    public void testFileNameNeedsOnePeriod() throws Exception {
        Path outDir = cleanDir("badName");
        Path infile = writeMarkerFile("DS_4.marker.txt", markerRows());

        new SplitFile_IFLFilePlugin()
                .inFile(infile.toString())
                .outputDir(outDir.toString() + "/")
                .maxSize(100)
                .performFunction(null);

        try (Stream<Path> files = Files.list(outDir)) {
            assertEquals("A name without exactly one period can't be split into dataset and table", 0, files.count());
        }
    }

    /** Marker name, chromosome and position for each site in the chromosome 7 Tutorial HapMap file. */
    private static List<String> markerRows() throws IOException {
        try (Stream<String> lines = Files.lines(Paths.get(TutorialConstants.TUTORIAL_DIR + "mdp_genotype_chrom7.hmp.txt"))) {
            return lines.skip(1)
                    .map(line -> line.split("\t"))
                    .map(tokens -> tokens[0] + "\t" + tokens[2] + "\t" + tokens[3])
                    .collect(Collectors.toList());
        }
    }

    private static Path writeMarkerFile(String name, List<String> rows) throws IOException {
        Path inDir = TEMP_DIR.resolve("input");
        Files.createDirectories(inDir);
        List<String> lines = new ArrayList<>();
        lines.add(HEADER);
        lines.addAll(rows);
        return Files.write(inDir.resolve(name), lines);
    }

    private static Path cleanDir(String name) throws IOException {
        Path dir = TEMP_DIR.resolve(name);
        if (Files.exists(dir)) {
            try (Stream<Path> files = Files.list(dir)) {
                for (Path file : files.collect(Collectors.toList())) {
                    Files.delete(file);
                }
            }
        }
        return Files.createDirectories(dir);
    }

}
