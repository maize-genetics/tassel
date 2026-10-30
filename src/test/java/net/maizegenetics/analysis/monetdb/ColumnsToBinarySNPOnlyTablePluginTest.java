/**
 * 
 */
package net.maizegenetics.analysis.monetdb;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import org.junit.Test;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.Assert.*;

/**
 * Uses the sites in the chromosome 9 and 10 Tutorial HapMap file as the SNP position list
 * ({@code SNPPos_chrom<N>.txt}), then writes per-SNP columns for every third site. Sites
 * missing from the input must be written as MonetDB nulls so the columns line up with the list.
 *
 * @author lcj34
 *
 */
public class ColumnsToBinarySNPOnlyTablePluginTest {

    private static final Path TEMP_DIR = Paths.get(GeneralConstants.TEMP_DIR, "ColumnsToBinarySNPOnlyTablePluginTest");
    private static final int[] CHROMOSOMES = {9, 10};
    private static final String INT_COLUMN = "Tutorial_NNonMiss";
    private static final String ALLELE_COLUMN = "Tutorial_MajAll";

    @Test
    public void testColumnsToBinarySNPOnlyTest() throws Exception {
        Path baseDir = cleanDir(TEMP_DIR);
        Path refDir = Files.createDirectories(baseDir.resolve("refGenomeFiles"));
        Path outDir = Files.createDirectories(baseDir.resolve("output"));
        Map<Integer, List<String[]>> sitesByChrom = readTutorialSites();

        List<String> input = new ArrayList<>();
        input.add("chr\tpos\t" + INT_COLUMN + "\t" + ALLELE_COLUMN);
        List<Integer> expectedInts = new ArrayList<>();
        List<Byte> expectedAlleles = new ArrayList<>();
        for (int chrom : CHROMOSOMES) {
            List<String[]> sites = sitesByChrom.get(chrom);
            Files.write(refDir.resolve("SNPPos_chrom" + chrom + ".txt"),
                    sites.stream().map(site -> site[3]).collect(Collectors.toList()));
            for (int i = 0; i < sites.size(); i++) {
                String[] site = sites.get(i);
                if (i % 3 == 0) {
                    int nonMissing = nonMissingCalls(site);
                    char majorAllele = site[1].charAt(0);
                    input.add(chrom + "\t" + site[3] + "\t" + nonMissing + "\t" + majorAllele);
                    expectedInts.add(nonMissing);
                    expectedAlleles.add(alleleCode(majorAllele));
                } else {
                    expectedInts.add(Integer.MIN_VALUE);
                    expectedAlleles.add(Byte.MIN_VALUE);
                }
            }
        }
        Path inputFile = Files.write(baseDir.resolve("MAFAndCov.popInfo.txt"), input);

        String outBase = outDir.toString() + "/";
        new ColumnsToBinarySNPOnlyTablePlugin()
                .inputFile(inputFile.toString())
                .refDir(refDir.toString())
                .outBase(outBase)
                .colsInt(INT_COLUMN)
                .colsAllele(ALLELE_COLUMN)
                .performFunction(null);

        int numSites = sitesByChrom.get(9).size() + sitesByChrom.get(10).size();
        int[] ints = readInts(outBase + INT_COLUMN + "_int.bin");
        assertEquals("One value per SNP in the position lists", numSites, ints.length);
        assertArrayEquals(expectedInts.stream().mapToInt(Integer::intValue).toArray(), ints);

        byte[] alleles = Files.readAllBytes(Paths.get(outBase + ALLELE_COLUMN + "_allelebyte.bin"));
        assertEquals(numSites, alleles.length);
        for (int i = 0; i < numSites; i++) {
            assertEquals("Allele code at site " + i, expectedAlleles.get(i).byteValue(), alleles[i]);
        }

        String copy = new String(Files.readAllBytes(Paths.get(outBase + "copy.sql")));
        assertTrue(copy, copy.contains("'" + INT_COLUMN + "_int.bin'"));
        assertTrue(copy, copy.contains("'" + ALLELE_COLUMN + "_allelebyte.bin'"));
    }

    /** HapMap rows for chromosomes 9 and 10, in file order (sorted by position). */
    private static Map<Integer, List<String[]>> readTutorialSites() throws IOException {
        try (Stream<String> lines = Files.lines(Paths.get(TutorialConstants.HAPMAP_CHR_9_10_FILENAME))) {
            Map<Integer, List<String[]>> sites = lines.skip(1)
                    .map(line -> line.split("\t"))
                    .collect(Collectors.groupingBy(tokens -> Integer.parseInt(tokens[2])));
            for (int chrom : CHROMOSOMES) {
                assertTrue("Tutorial file should have chromosome " + chrom, sites.get(chrom).size() > 100);
            }
            return sites;
        }
    }

    private static int nonMissingCalls(String[] site) {
        int count = 0;
        for (int i = 11; i < site.length; i++) {
            if (!site[i].equals("N")) count++;
        }
        return count;
    }

    private static byte alleleCode(char allele) {
        return (byte) "ACGT".indexOf(allele);
    }

    private static int[] readInts(String file) throws IOException {
        ByteBuffer buffer = ByteBuffer.wrap(Files.readAllBytes(Paths.get(file))).order(ByteOrder.LITTLE_ENDIAN);
        int[] values = new int[buffer.remaining() / Integer.BYTES];
        buffer.asIntBuffer().get(values);
        return values;
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
