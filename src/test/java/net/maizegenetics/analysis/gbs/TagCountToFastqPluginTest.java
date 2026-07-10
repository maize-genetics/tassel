/*
 * TagCountToFastqPluginTest
 */
package net.maizegenetics.analysis.gbs;

import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Set;

import net.maizegenetics.analysis.gbs.v2.GBSSimData;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.util.Utils;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Property-based rehabilitation of the legacy GBSv1 {@code TagCountToFastqPlugin} test.
 *
 * <p>The original test MD5-compared the produced {@code .fq.gz} against a downloaded golden fixture.
 * Instead this test self-generates a deterministic master tag list with {@link GBSv1SimData}, runs the
 * plugin, and asserts on structural properties of the exported FASTQ: one 4-line record per tag, each
 * carrying the {@code @length=..count=..} header, sequence/quality lengths matching the header length,
 * and the exported sequence set covering exactly the known simulated tags.</p>
 *
 * @author terry (original), rehabilitated for self-generated data
 */
public class TagCountToFastqPluginTest {

    @Test
    public void testPerformFunction() throws Exception {
        GBSv1SimData sim = GBSv1SimData.createUnder("TagCountToFastq");
        sim.buildMasterTagCounts();

        TagCounts masterTags = new TagCounts(sim.masterTagCounts.toString(), FilePacking.Byte);
        int expectedRecords = masterTags.getTagCount();
        assertTrue("Master tag list should be non-empty", expectedRecords > 0);

        Path fastqOut = sim.v1Dir.resolve("tagsForAlign.fq.gz");
        TagCountToFastqPlugin plugin = new TagCountToFastqPlugin();
        plugin.setParameters(new String[]{
                "-i", sim.masterTagCounts.toString(),
                "-o", fastqOut.toString(),
                "-c", "1"
        });
        plugin.performFunction(null);

        int records = 0;
        Set<String> exportedSeqs = new HashSet<>();
        try (BufferedReader br = Utils.getBufferedReader(fastqOut.toString())) {
            String header;
            while ((header = br.readLine()) != null) {
                String seq = br.readLine();
                br.readLine();           // "+" separator
                String qual = br.readLine();

                assertTrue("FASTQ header should carry the length= prefix: " + header,
                        header.startsWith("@length="));
                int countIdx = header.indexOf("count=");
                assertTrue("FASTQ header should carry the count= field: " + header, countIdx > 0);

                int declaredLength = Integer.parseInt(header.substring("@length=".length(), countIdx));
                assertEquals("Sequence length should match the header length", declaredLength, seq.length());
                assertEquals("Quality length should match the header length", declaredLength, qual.length());
                assertTrue("Exported tag should retain the ApeKI cut remnant: " + seq,
                        seq.startsWith(GBSSimData.CUT_REMNANT));

                exportedSeqs.add(seq);
                records++;
            }
        }

        assertEquals("One FASTQ record should be exported per master tag", expectedRecords, records);
        for (GBSSimData.TagInfo info : sim.sim.tagInfos) {
            assertTrue("Exported FASTQ should contain simulated tag " + info.sequence,
                    exportedSeqs.contains(info.sequence));
        }
    }
}
