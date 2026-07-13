/**
 *
 */
package net.maizegenetics.analysis.gbs.v2;

import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.Assert.*;

/**
 * Tests {@link GetTagSequenceFromDBPlugin} against a self-generated {@link GBSSimData} database.
 * Rather than depending on a pre-built {@code GBSv2.db}, each test builds its own and asserts on
 * known simulated tags and output line counts.
 *
 * @author lcj34 (original), rewritten for simulated data
 */
public class GetTagSequenceFromDBPluginTest {

    @Test
    public void testSingleTagInput() throws Exception {
        GBSSimData sim = GBSSimData.createUnder("GetTagSeq_single");
        sim.buildTagDB();

        Path out = sim.baseDir.resolve("dbSingleTag.txt");
        String knownTag = sim.aKnownTagSequence();
        new GetTagSequenceFromDBPlugin()
                .inputDB(sim.dbFile.toString())
                .outputFile(out.toString())
                .tagSequence(knownTag)
                .performFunction(null);

        List<String> lines = Files.readAllLines(out);
        assertTrue("A known simulated tag should be reported as found",
                lines.stream().anyMatch(l -> l.equals(knownTag)));
        assertFalse("A known tag should not be reported as NOT found",
                lines.stream().anyMatch(l -> l.contains("NOT found")));
    }

    @Test
    public void testSingleTagNOTFound() throws Exception {
        GBSSimData sim = GBSSimData.createUnder("GetTagSeq_notFound");
        sim.buildTagDB();

        Path out = sim.baseDir.resolve("dbSingleTag.txt");
        new GetTagSequenceFromDBPlugin()
                .inputDB(sim.dbFile.toString())
                .outputFile(out.toString())
                .tagSequence(sim.anUnknownTagSequence())
                .performFunction(null);

        List<String> lines = Files.readAllLines(out);
        assertTrue("A bogus tag should be reported as NOT found",
                lines.stream().anyMatch(l -> l.contains("NOT found")));
    }

    @Test
    public void testMultipleTagInput() throws Exception {
        GBSSimData sim = GBSSimData.createUnder("GetTagSeq_all");
        sim.buildTagDB();

        Path out = sim.baseDir.resolve("dbMultipleTag.txt");
        new GetTagSequenceFromDBPlugin()
                .inputDB(sim.dbFile.toString())
                .outputFile(out.toString())
                .performFunction(null);

        // One header line ("Tags") plus one line per distinct tag.
        long tagLines = Files.readAllLines(out).stream()
                .filter(l -> !l.isEmpty() && !l.equals("Tags"))
                .count();
        assertEquals(sim.expectedDistinctTags, tagLines);
    }
}
