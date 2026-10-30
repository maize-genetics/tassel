/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.util.LoggingUtils;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.Assert.*;

/**
 * Builds breakpoint (.pa.txt) files and GID mapping files for taxa in the Tutorial HapMap
 * panel, then checks the hap_breakpoint and breakpoint_set intermediate files.
 *
 * @author lcj34
 *
 */
public class HapBreakpoints_IFLFilePluginTest {

    private static final Path TEMP_DIR = Paths.get(GeneralConstants.TEMP_DIR, "HapBreakpoints_IFLFilePluginTest");
    private static final String SET_NAME = "ames_unimp";
    private static final String MAPSET = "AGPV3";
    private static final String METHOD = "FILLIN";
    private static final String SRC_DATASET = "ZeaGBSv27_20150803_AGPv3_Ames";
    private static final String UNKNOWN_TAXON = "NotInTaxaMap";
    private static final String MISSING_GID = "49999999";

    @Test
    public void testHapBreakpointsIFLFilePlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        Path baseDir = cleanDir(TEMP_DIR);
        Path breakDir = Files.createDirectories(baseDir.resolve("breakpoints"));
        Path outDir = Files.createDirectories(baseDir.resolve("out"));

        List<String> panel = tutorialTaxa();
        List<String> donors = panel.subList(0, 4);
        List<String> taxa = panel.subList(4, 7);
        Map<String, String> donorGids = gids(donors, 42000000);
        Map<String, String> taxaGids = gids(taxa, 43000000);

        // The donor map appends library IDs to the names, which the plugin strips.
        Path donorMapFile = writeMapFile(baseDir.resolve("donorMap.txt"), "TaxaColumn", donorGids, ":LIB1");
        Path taxaMapFile = writeMapFile(baseDir.resolve("taxaMap.txt"), "Name", taxaGids, "");

        // Chromosome 1: two taxa from the map, one taxon missing from it, and a block with no second donor.
        writeBreakFile(breakDir.resolve("chr1.pa.txt"), donors, Arrays.asList(
                taxa.get(0) + "\t1:100:200:0:1\t1:200:300:2:-1",
                taxa.get(1) + "\t1:100:300:3:3",
                UNKNOWN_TAXON + "\t1:100:300:1:2"));
        writeBreakFile(breakDir.resolve("chr2.pa.txt"), donors, Collections.singletonList(
                taxa.get(2) + "\t2:500:900:2:0"));
        // Not a .pa.txt file, so it must be ignored.
        Files.write(breakDir.resolve("notes.txt"), Collections.singletonList("ignore me"));

        new HapBreakpoints_IFLFilePlugin()
                .outputDir(outDir.toString() + "/")
                .breakFile(breakDir.toString())
                .setName(SET_NAME)
                .mapset(MAPSET)
                .src_dataset(SRC_DATASET)
                .method(METHOD)
                .donorMapFile(donorMapFile.toString())
                .taxaMapFile(taxaMapFile.toString())
                .performFunction(null);

        List<String> expectedHap = Arrays.asList(
                "taxa\tchr\tposition_range\tdonor1\tdonor2\tstatus\tbreakpoint_set_name",
                hapLine(taxaGids.get(taxa.get(0)), "1", "[100,200)", donorGids.get(donors.get(0)), donorGids.get(donors.get(1))),
                hapLine(taxaGids.get(taxa.get(0)), "1", "[200,300)", donorGids.get(donors.get(2)), "-1"),
                hapLine(taxaGids.get(taxa.get(1)), "1", "[100,300)", donorGids.get(donors.get(3)), donorGids.get(donors.get(3))),
                hapLine(MISSING_GID, "1", "[100,300)", donorGids.get(donors.get(1)), donorGids.get(donors.get(2))),
                hapLine(taxaGids.get(taxa.get(2)), "2", "[500,900)", donorGids.get(donors.get(2)), donorGids.get(donors.get(0))));
        assertEquals(expectedHap, Files.readAllLines(outDir.resolve(SET_NAME + ".hapbreakpoints")));

        assertEquals(Arrays.asList(
                "name\tmethod\tmapset_name\tsource_dataset_name\tstatus",
                SET_NAME + "\t" + METHOD + "\t" + MAPSET + "\t" + SRC_DATASET + "\t1"),
                Files.readAllLines(outDir.resolve(SET_NAME + ".breakpoint_set")));

        assertEquals(Collections.singletonList(UNKNOWN_TAXON),
                Files.readAllLines(outDir.resolve(SET_NAME + ".taxaNotFound.txt")));
        assertEquals("Every taxon line is listed, including the one given a placeholder GID",
                Arrays.asList(taxa.get(0), taxa.get(1), UNKNOWN_TAXON, taxa.get(2)),
                Files.readAllLines(outDir.resolve(SET_NAME + ".taxaFound.txt")));
    }

    @Test
    public void testCreateNameGidMap() throws Exception {
        Path baseDir = cleanDir(TEMP_DIR.resolve("nameGidMap"));
        Map<String, String> gids = gids(tutorialTaxa().subList(0, 3), 42000000);
        Path mapFile = writeMapFile(baseDir.resolve("map.txt"), "TaxaColumn", gids, ":LIB1");

        HapBreakpoints_IFLFilePlugin plugin = new HapBreakpoints_IFLFilePlugin();
        assertEquals("Library IDs after a colon are stripped", gids, plugin.createNameGidMap(mapFile.toString(), "TAXACOLUMN"));
        assertNull("A missing name column returns null", plugin.createNameGidMap(mapFile.toString(), "NAME"));
    }

    private static String hapLine(String taxonGid, String chr, String range, String donor1, String donor2) {
        return taxonGid + "\t" + chr + "\t" + range + "\t" + donor1 + "\t" + donor2 + "\t1\t" + SET_NAME;
    }

    private static List<String> tutorialTaxa() throws IOException {
        try (BufferedReader br = Files.newBufferedReader(Paths.get(TutorialConstants.TUTORIAL_DIR + "mdp_genotype_chrom7.hmp.txt"))) {
            String[] header = br.readLine().split("\t");
            return Arrays.asList(header).subList(11, header.length);
        }
    }

    private static Map<String, String> gids(List<String> names, int firstGid) {
        Map<String, String> gids = new LinkedHashMap<>();
        for (int i = 0; i < names.size(); i++) {
            gids.put(names.get(i), Integer.toString(firstGid + i));
        }
        return gids;
    }

    private static Path writeMapFile(Path file, String nameHeader, Map<String, String> gids, String nameSuffix) throws IOException {
        List<String> lines = new ArrayList<>();
        lines.add(nameHeader + "\tGID\tSource");
        gids.forEach((name, gid) -> lines.add(name + nameSuffix + "\t" + gid + "\tTutorial"));
        return Files.write(file, lines);
    }

    private static void writeBreakFile(Path file, List<String> donors, List<String> taxaLines) throws IOException {
        List<String> lines = new ArrayList<>();
        lines.add(donors.size() + "\t" + taxaLines.size());
        lines.add("#Donor Haplotypes");
        for (int i = 0; i < donors.size(); i++) {
            lines.add(i + "\t" + donors.get(i));
        }
        lines.add("#Block are defined chr:startPos:endPos:donor1:donor2 (-1 means no hypothesis)");
        lines.addAll(taxaLines);
        Files.write(file, lines);
    }

    private static Path cleanDir(Path dir) throws IOException {
        if (Files.exists(dir)) {
            try (Stream<Path> paths = Files.walk(dir)) {
                for (Path path : paths.sorted(Comparator.reverseOrder()).collect(Collectors.toList())) {
                    Files.delete(path);
                }
            }
        }
        return Files.createDirectories(dir);
    }

}
