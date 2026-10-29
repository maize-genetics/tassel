/**
 * 
 */
package net.maizegenetics.util;

import net.maizegenetics.analysis.gbs.v2.GBSSimData;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

import static org.junit.Assert.*;

/**
 * This test calls ConvertOldFastqToModernFormatPlugin to
 * transform fastq files from an old illumina format which
 * contained multiple taxa with barcode added to each read sequence
 * into the newer format consisting of 1 fastq file per taxon
 * and no barcode in the read sequences.
 *
 * <p>The multiplexed input fastq and key file come from {@link GBSSimData}, so the
 * expected reads for each taxon are known exactly.</p>
 * 
 * @author lcj34
 *
 */
public class ConvertOldFastqToModernFormatTest {

    private static final String PROJECT_NAME = "chrom9-10junit";

    @Test
    public void testConvertFiles() throws IOException {
        GBSSimData sim = GBSSimData.createUnder("ConvertOldFastqToModernFormat");
        Path outputDir = sim.baseDir.resolve("converted");
        Files.createDirectories(outputDir);

        new ConvertOldFastqToModernFormatPlugin()
                .inputFile(sim.fastqDir.toString())
                .keyFile(sim.keyFile.toString())
                .outputDir(outputDir.toString())
                .projectName(PROJECT_NAME)
                .enzyme(GBSSimData.ENZYME)
                .performFunction(null);

        try (Stream<Path> files = Files.list(outputDir)) {
            assertEquals("One fastq file per taxon", sim.taxa.size(), files.count());
        }

        for (int taxonIdx = 0; taxonIdx < sim.taxa.size(); taxonIdx++) {
            Set<String> expectedTags = new HashSet<>();
            for (GBSSimData.TagInfo info : sim.tagInfos) {
                for (int carrier : info.taxaIndices) {
                    if (carrier == taxonIdx) expectedTags.add(info.sequence);
                }
            }

            String taxon = sim.taxa.get(taxonIdx);
            Path outFile = outputDir.resolve(PROJECT_NAME + "_" + taxon + "_" + GBSSimData.FLOWCELL + "_"
                    + GBSSimData.LANE + "_" + sim.barcodes.get(taxonIdx) + ".fastq");
            assertTrue("Missing output for " + taxon + ": " + outFile, Files.exists(outFile));

            List<String[]> records = readFastq(outFile);
            assertEquals("Read count for " + taxon,
                    expectedTags.size() * GBSSimData.DEPTH_PER_TAXON, records.size());

            Set<String> seenTags = new HashSet<>();
            for (String[] record : records) {
                assertTrue("Header should be copied unchanged: " + record[0],
                        record[0].startsWith("@" + GBSSimData.FLOWCELL + ":"));
                assertTrue("Barcode should be removed and the read should be one of " + taxon + "'s tags: " + record[1],
                        expectedTags.contains(record[1]));
                assertEquals("+", record[2]);
                assertEquals("Quality must be trimmed the same as the sequence",
                        record[1].length(), record[3].length());
                seenTags.add(record[1]);
            }
            assertEquals("Every tag carried by " + taxon + " should appear", expectedTags, seenTags);
        }
    }

    private static List<String[]> readFastq(Path file) throws IOException {
        List<String[]> records = new ArrayList<>();
        try (BufferedReader br = Utils.getBufferedReader(file.toString())) {
            String header;
            while ((header = br.readLine()) != null) {
                records.add(new String[]{header, br.readLine(), br.readLine(), br.readLine()});
            }
        }
        return records;
    }

}
