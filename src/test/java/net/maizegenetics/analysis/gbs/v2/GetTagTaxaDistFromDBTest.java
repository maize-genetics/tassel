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
 * Tests {@link GetTagTaxaDistFromDBPlugin} against a self-generated {@link GBSSimData} database.
 *
 * @author lcj34 (original), rewritten for simulated data
 */
public class GetTagTaxaDistFromDBTest {

    @Test
    public void testTagTaxaDistFromDBPlugin() throws Exception {
        GBSSimData sim = GBSSimData.createUnder("GetTagTaxaDist");
        sim.buildTagDB();

        Path out = sim.baseDir.resolve("tagsTaxaDepthOutput.txt");
        new GetTagTaxaDistFromDBPlugin()
                .inputDB(sim.dbFile.toString())
                .outputFile(out.toString())
                .performFunction(null);

        List<String> lines = Files.readAllLines(out);
        assertFalse("Output should not be empty", lines.isEmpty());
        assertTrue("First line should be the Tag/taxa header", lines.get(0).startsWith("Tag"));

        // One data row per distinct tag, each with a depth column per taxon.
        long dataRows = lines.stream().skip(1).filter(l -> !l.isEmpty()).count();
        assertEquals(sim.expectedDistinctTags, dataRows);

        int expectedColumns = sim.taxa.size() + 1; // tag sequence + one column per taxon
        assertEquals(expectedColumns, lines.get(1).split("\t").length);
    }
}
