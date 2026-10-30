package net.maizegenetics.analysis.gobii;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.Assert.*;

/**
 * GWAS_IFLPlugin, MarkerDNARun_IFLFilePlugin, MonetDB_IFLFilePlugin and PreProcessGOBIIMappingFilePlugin
 * need a live GOBII PostgreSQL database, so these tests cover what happens before connecting: reading the
 * config file, and stopping without writing intermediate files when no connection can be made.
 */
public class GOBIIDbUtilsTest {

    private static final Path TEMP_DIR = Paths.get(GeneralConstants.TEMP_DIR, "GOBIIDbUtilsTest");
    private static final String NO_PASSWORD = "user=gobii\nDB=gobii_maize\n";

    @Test
    public void testConnectToDBRejectsIncompleteConfig() throws IOException {
        assertNull("missing user", GOBIIDbUtils.connectToDB(writeConfig("noUser.txt", "password=secret\nDB=gobii_maize\n")));
        assertNull("missing password", GOBIIDbUtils.connectToDB(writeConfig("noPassword.txt", NO_PASSWORD)));
        assertNull("missing DB", GOBIIDbUtils.connectToDB(writeConfig("noDB.txt", "user=gobii\npassword=secret\n")));
        assertNull("missing file", GOBIIDbUtils.connectToDB(TEMP_DIR.resolve("doesNotExist.txt").toString()));
    }

    @Test
    public void testGWAS_IFLPluginStopsWithoutDB() throws IOException {
        Path outDir = cleanDir("gwas");
        Path gwasFile = Files.write(TEMP_DIR.resolve("gwas.txt"), Arrays.asList(
                "Trait\tMarker\tChr\tPos\tp\tr2\tdf",
                "EarHT\tPZA01426.1\t7\t729478\t0.01\t0.2\t1"));

        new GWAS_IFLPlugin()
                .b4rConfigFile(writeConfig("gwasConfig.txt", NO_PASSWORD))
                .outputDir(outDir.toString() + "/DNDiff_")
                .inputFile(gwasFile.toString())
                .methodIds("1,2,3")
                .expID("3")
                .statNames("p,r2,df")
                .performFunction(null);

        assertEmpty(outDir);
    }

    @Test
    public void testMarkerDNARun_IFLFilePluginStopsWithoutDB() throws IOException {
        Path outDir = cleanDir("markerDNARun");

        new MarkerDNARun_IFLFilePlugin()
                .dbConfigFile(writeConfig("markerConfig.txt", NO_PASSWORD))
                .inputFile(TutorialConstants.TUTORIAL_DIR + "mdp_genotype_chrom7.hmp.txt")
                .outputFileDir(outDir.toString() + "/")
                .refFile(TEMP_DIR.resolve("unusedReference.fa").toString())
                .mappingFile(TEMP_DIR.resolve("unusedMapping.txt").toString())
                .mapsetName("AGPV3")
                .expName("282Maize_wgs321")
                .platformName("Hapmap_321")
                .refName("Maize_AGPv3")
                .datasetName("ZeaWGSraw_hmp321_20160628_AGPv3_282")
                .performFunction(null);

        assertEmpty(outDir);
    }

    @Test
    public void testPreProcessGOBIIMappingFilePluginStopsWithoutDB() throws IOException {
        Path outDir = cleanDir("preProcess");

        new PreProcessGOBIIMappingFilePlugin()
                .dbConfigFile(writeConfig("preProcessConfig.txt", NO_PASSWORD))
                .outputDir(outDir.toString() + "/")
                .datasetName("ZeaWGS_hmp321_imputed_AGPv3")
                .mappingFile(TEMP_DIR.resolve("unusedMapping.txt").toString())
                .performFunction(null);

        assertEmpty(outDir);
    }

    @Test
    public void testMonetDB_IFLFilePluginWritesNoIdsWithoutDB() throws IOException {
        Path outDir = cleanDir("monetDB");

        new MonetDB_IFLFilePlugin()
                .dbConfigFile(writeConfig("monetConfig.txt", NO_PASSWORD))
                .outputDir(outDir.toString() + "/DS_1")
                .datasetName("ZeaGBSv27_20160209_AGPv2_282")
                .performFunction(null);

        // The plugin opens both output files before connecting.
        for (String suffix : new String[]{".dnarun_id", ".marker_id"}) {
            Path idFile = outDir.resolve("DS_1" + suffix);
            assertTrue(Files.notExists(idFile) || Files.size(idFile) == 0);
        }
    }

    private static void assertEmpty(Path dir) throws IOException {
        try (Stream<Path> files = Files.list(dir)) {
            assertEquals("No intermediate files should be written without a database",
                    Arrays.asList(), files.collect(Collectors.toList()));
        }
    }

    private static String writeConfig(String name, String contents) throws IOException {
        Files.createDirectories(TEMP_DIR);
        return Files.write(TEMP_DIR.resolve(name), contents.getBytes()).toString();
    }

    private static Path cleanDir(String name) throws IOException {
        Path dir = TEMP_DIR.resolve(name);
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
