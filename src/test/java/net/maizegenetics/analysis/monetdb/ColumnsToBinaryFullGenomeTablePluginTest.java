/**
 * 
 */
package net.maizegenetics.analysis.monetdb;

import net.maizegenetics.constants.GeneralConstants;
import org.junit.Test;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.Assert.*;

/**
 * Writes MonetDB binary columns for range data on chromosomes 4 and 5 of the short test reference
 * ({@code dataFiles/ReferenceGenomeSequenceTests/testChrom3456_short.fa}, 239 and 1,560 bp).
 * The plugin writes one value per reference position, so every column must hold 1,799 values.
 *
 * @author lcj34
 *
 */
public class ColumnsToBinaryFullGenomeTablePluginTest {

    private static final Path TEMP_DIR = Paths.get(GeneralConstants.TEMP_DIR, "ColumnsToBinaryFullGenomeTablePluginTest");
    private static final String REFERENCE = GeneralConstants.DATA_DIR + "ReferenceGenomeSequenceTests/testChrom3456_short.fa";
    private static final int CHR4_SIZE = 239;
    private static final int CHR5_SIZE = 1560;
    private static final String HEADER = "chr\tstart\tend\tmean\tcorrected\tcount";
    // Ranges are 1-based, start inclusive and end exclusive.
    private static final List<String> CHR4_LINES = Arrays.asList(
            "4\t10\t20\t1.5\t0.25\t3",
            "4\t50\t51\t2.0\tNA\t7",
            "4\t200\t240\t-0.5\t1.0\t1");
    private static final List<String> CHR5_LINES = Arrays.asList(
            "5\t1\t5\t3.0\t3.5\t2",
            "5\t1000\t1010\t4.0\t4.5\tNA");

    @Test
    public void testColumnsToBinaryFullGenomeTest() throws Exception {
        Path baseDir = cleanDir(TEMP_DIR.resolve("singleFile"));
        Path inputFile = baseDir.resolve("xoDensity.txt");
        Files.write(inputFile, concat(Arrays.asList(HEADER), CHR4_LINES, CHR5_LINES));
        String outBase = baseDir.resolve("out").resolve("xoDensity").toString();

        runPlugin(inputFile.toString(), outBase);

        float missingFloat = -Float.MAX_VALUE;
        float[] mean = filled(missingFloat);
        float[] corrected = filled(missingFloat);
        int[] count = new int[CHR4_SIZE + CHR5_SIZE];
        Arrays.fill(count, Integer.MIN_VALUE);
        setRange(mean, 4, 10, 20, 1.5f);     setRange(corrected, 4, 10, 20, 0.25f);  setRange(count, 4, 10, 20, 3);
        setRange(mean, 4, 50, 51, 2.0f);                                             setRange(count, 4, 50, 51, 7);
        setRange(mean, 4, 200, 240, -0.5f);  setRange(corrected, 4, 200, 240, 1.0f); setRange(count, 4, 200, 240, 1);
        setRange(mean, 5, 1, 5, 3.0f);       setRange(corrected, 5, 1, 5, 3.5f);     setRange(count, 5, 1, 5, 2);
        setRange(mean, 5, 1000, 1010, 4.0f); setRange(corrected, 5, 1000, 1010, 4.5f);

        assertArrayEquals(mean, readFloats(outBase + "mean_real.bin"), 0f);
        assertArrayEquals(corrected, readFloats(outBase + "corrected_real.bin"), 0f);
        assertArrayEquals(count, readInts(outBase + "count_int.bin"));

        String copy = new String(Files.readAllBytes(Paths.get(outBase + "copy.sql")));
        for (String column : new String[]{"xoDensitymean_real.bin", "xoDensitycorrected_real.bin", "xoDensitycount_int.bin"}) {
            assertTrue("copy.sql should load " + column + ": " + copy, copy.contains("'" + column + "'"));
        }
        String create = new String(Files.readAllBytes(Paths.get(outBase + "create.sql")));
        assertTrue(create, create.startsWith("CREATE TABLE ... (chr int,pos int,"));
        assertTrue(create, create.contains("xoDensitycount int"));
    }

    @Test
    public void testInputDirectorySplitByChromosome() throws Exception {
        Path singleDir = cleanDir(TEMP_DIR.resolve("singleFileForCompare"));
        Path singleFile = singleDir.resolve("xoDensity.txt");
        Files.write(singleFile, concat(Arrays.asList(HEADER), CHR4_LINES, CHR5_LINES));
        String singleOut = singleDir.resolve("out").resolve("xoDensity").toString();
        runPlugin(singleFile.toString(), singleOut);

        Path splitDir = cleanDir(TEMP_DIR.resolve("byChromosome"));
        Path inputDir = Files.createDirectories(splitDir.resolve("input"));
        Files.write(inputDir.resolve("chr04.txt"), concat(Arrays.asList(HEADER), CHR4_LINES));
        Files.write(inputDir.resolve("chr05.txt"), concat(Arrays.asList(HEADER), CHR5_LINES));
        String splitOut = splitDir.resolve("out").resolve("xoDensity").toString();
        runPlugin(inputDir.toString(), splitOut);

        for (String column : new String[]{"mean_real.bin", "corrected_real.bin", "count_int.bin"}) {
            assertArrayEquals("Splitting the input by chromosome shouldn't change " + column,
                    Files.readAllBytes(Paths.get(singleOut + column)), Files.readAllBytes(Paths.get(splitOut + column)));
        }
    }

    private static void runPlugin(String inputFile, String outBase) throws IOException {
        Files.createDirectories(Paths.get(outBase).getParent());
        new ColumnsToBinaryFullGenomeTablePlugin()
                .inputFile(inputFile)
                .refFile(REFERENCE)
                .outBase(outBase)
                .colsFloat("mean,corrected")
                .colsInt("count")
                .range(true)
                .performFunction(null);
    }

    private static float[] filled(float value) {
        float[] values = new float[CHR4_SIZE + CHR5_SIZE];
        Arrays.fill(values, value);
        return values;
    }

    private static int offset(int chr) {
        return chr == 4 ? 0 : CHR4_SIZE;
    }

    private static void setRange(float[] values, int chr, int start, int end, float value) {
        Arrays.fill(values, offset(chr) + start - 1, offset(chr) + end - 1, value);
    }

    private static void setRange(int[] values, int chr, int start, int end, int value) {
        Arrays.fill(values, offset(chr) + start - 1, offset(chr) + end - 1, value);
    }

    private static float[] readFloats(String file) throws IOException {
        ByteBuffer buffer = ByteBuffer.wrap(Files.readAllBytes(Paths.get(file))).order(ByteOrder.LITTLE_ENDIAN);
        float[] values = new float[buffer.remaining() / Float.BYTES];
        buffer.asFloatBuffer().get(values);
        return values;
    }

    private static int[] readInts(String file) throws IOException {
        ByteBuffer buffer = ByteBuffer.wrap(Files.readAllBytes(Paths.get(file))).order(ByteOrder.LITTLE_ENDIAN);
        int[] values = new int[buffer.remaining() / Integer.BYTES];
        buffer.asIntBuffer().get(values);
        return values;
    }

    @SafeVarargs
    private static List<String> concat(List<String>... parts) {
        return Stream.of(parts).flatMap(List::stream).collect(Collectors.toList());
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
