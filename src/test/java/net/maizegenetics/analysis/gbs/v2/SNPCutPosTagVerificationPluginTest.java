/**
 *
 */
package net.maizegenetics.analysis.gbs.v2;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.tag.TagData;
import net.maizegenetics.dna.tag.TagDataSQLite;

import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.Assert.*;

/**
 * Tests {@link SNPCutPosTagVerificationPlugin} against a self-generated {@link GBSSimData} database.
 * The simulated pipeline injects a known number of cut positions (one per locus) and SNP positions
 * (one per SNP locus); this test verifies those counts and that the debugging plugin can dump the
 * tags at a cut and a SNP position.
 *
 * @author lcj34 (original), rewritten for simulated data
 */
public class SNPCutPosTagVerificationPluginTest {

    @Test
    public void SNPCutPosTagVerificationTest() throws Exception {
        System.out.println("Running SNPCutPosTagVerificationTest");
        GBSSimData sim = GBSSimData.createUnder("SNPCutPos");
        sim.buildDatabaseThroughSam();
        sim.runDiscovery(true);

        TagData tagData = new TagDataSQLite(sim.dbFile.toString());
        PositionList cutPositions = tagData.getTagCutPositions(true);
        PositionList snpPositions = tagData.getSNPPositions();
        ((TagDataSQLite) tagData).close();

        assertEquals("One cut position per simulated locus", sim.numLoci, cutPositions.size());
        assertEquals("One SNP position per simulated SNP locus", sim.numSnpLoci, snpPositions.size());

        // Dump the tags at the first cut position.
        Position cut = cutPositions.get(0);
        Path cutOut = sim.baseDir.resolve("cutPositionData.txt");
        new SNPCutPosTagVerificationPlugin()
                .inputDB(sim.dbFile.toString())
                .cutOrSnpPosition(cut.getPosition())
                .chrom(cut.getChromosome().getName())
                .positionType("cut")
                .strand(cut.getStrand())
                .outputFile(cutOut.toString())
                .performFunction(null);
        assertTrue("Cut position dump file should exist", Files.exists(cutOut));

        // Dump the tags at the first SNP position.
        Position snp = snpPositions.get(0);
        Path snpOut = sim.baseDir.resolve("snpPositionData.txt");
        new SNPCutPosTagVerificationPlugin()
                .inputDB(sim.dbFile.toString())
                .cutOrSnpPosition(snp.getPosition())
                .chrom(snp.getChromosome().getName())
                .positionType("snp")
                .strand(snp.getStrand())
                .outputFile(snpOut.toString())
                .performFunction(null);
        assertTrue("SNP position dump file should exist", Files.exists(snpOut));

        System.out.println("SNPCutPosTagVerificationTest finished successfully !!!");
    }
}
