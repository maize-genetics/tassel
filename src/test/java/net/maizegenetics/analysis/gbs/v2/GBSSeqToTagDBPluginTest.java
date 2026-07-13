/*
 * GBSSeqToTagDBPluginTest
 */
package net.maizegenetics.analysis.gbs.v2;

import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagData;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.util.LoggingUtils;

import org.junit.Test;

import static org.junit.Assert.*;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Set;

/**
 * Property-based tests for {@link GBSSeqToTagDBPlugin} that build their own tiny GBSv2 database
 * from {@link GBSSimData} (deterministic simulated FASTQ/key/reference) rather than relying on
 * large downloaded fixtures. Assertions check invariants (distinct tag counts, append/keep-old-data
 * stability) instead of byte-exact golden files.
 *
 * @author Ed Buckler (original), rewritten for simulated data
 */
public class GBSSeqToTagDBPluginTest {

    public GBSSeqToTagDBPluginTest() {
    }

    @Test
    public void testGBSSeqToTagDBPlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        GBSSimData sim = GBSSimData.createUnder("GBSSeqToTagDB_basic");
        sim.buildTagDB();

        assertEquals("DB should contain exactly the simulated distinct tags",
                sim.expectedDistinctTags, sim.dbTagCount());
    }

    @Test
    public void testKeepOldData() throws Exception {
        LoggingUtils.setupDebugLogging();
        GBSSimData sim = GBSSimData.createUnder("GBSSeqToTagDB_keepOld");

        // First run: clean build.
        sim.buildTagDB();
        int firstCount = sim.dbTagCount();
        assertEquals(sim.expectedDistinctTags, firstCount);

        // Second run over the same files WITHOUT deleting old data. This exercises the append
        // path (TAS-1120): depths are maintained so removeTagsWithoutReplication does not drop tags.
        new GBSSeqToTagDBPlugin()
                .enzyme(GBSSimData.ENZYME)
                .inputDirectory(sim.fastqDir.toString())
                .outputDatabaseFile(sim.dbFile.toString())
                .keyFile(sim.keyFile.toString())
                .kmerLength(GBSSimData.TAG_LENGTH)
                .minKmerCount(GBSSimData.MIN_KMER_COUNT)
                .minimumQualityScore(20)
                .deleteOldData(false)
                .performFunction(null);

        assertEquals("Tag count must be stable after keep-old-data append", firstCount, sim.dbTagCount());
    }

    @Test
    public void testTagExportPlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        GBSSimData sim = GBSSimData.createUnder("GBSSeqToTagDB_export");
        sim.buildTagDB();
        sim.exportAndAlign();

        // Every distinct tag should be exported (as an @tagSeq= FASTQ record).
        long exported = java.nio.file.Files.lines(java.nio.file.Paths.get(sim.samFile.toString()))
                .filter(l -> l.startsWith("tagSeq=")).count();
        assertEquals(sim.expectedDistinctTags, exported);
    }

    @Test
    public void testSAMImportPlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        GBSSimData sim = GBSSimData.createUnder("GBSSeqToTagDB_sam");
        sim.buildDatabaseThroughSam();

        // Importing the (aligner-free) SAM must not change the number of tags, and cut positions
        // should now be present for every simulated locus.
        assertEquals(sim.expectedDistinctTags, sim.dbTagCount());
        TagData tagData = new TagDataSQLite(sim.dbFile.toString());
        int cutPositions = tagData.getTagCutPositions(true).size();
        ((TagDataSQLite) tagData).close();
        assertEquals("Each simulated locus should yield one cut position", sim.numLoci, cutPositions);
    }

    @Test
    public void GBSSeqToTagDBPluginAppendTest() throws Exception {
        GBSSimData sim = GBSSimData.createUnder("GBSSeqToTagDB_append");

        sim.buildTagDB();
        TagData tdw = new TagDataSQLite(sim.dbFile.toString());
        Set<Tag> firstRunTags = tdw.getTags();
        int firstDepth = tdw.getAllTagsTaxaMap().values().stream()
                .mapToInt(td -> td.totalDepth()).sum();
        ((TagDataSQLite) tdw).close();

        // Append the same reads again (no delete). Same tags, but greater total depth.
        new GBSSeqToTagDBPlugin()
                .enzyme(GBSSimData.ENZYME)
                .inputDirectory(sim.fastqDir.toString())
                .outputDatabaseFile(sim.dbFile.toString())
                .keyFile(sim.keyFile.toString())
                .kmerLength(GBSSimData.TAG_LENGTH)
                .minKmerCount(GBSSimData.MIN_KMER_COUNT)
                .minimumQualityScore(20)
                .deleteOldData(false)
                .performFunction(null);

        tdw = new TagDataSQLite(sim.dbFile.toString());
        Set<Tag> secondRunTags = tdw.getTags();
        int secondDepth = tdw.getAllTagsTaxaMap().values().stream()
                .mapToInt(td -> td.totalDepth()).sum();
        ((TagDataSQLite) tdw).close();

        assertEquals(firstRunTags.size(), secondRunTags.size());
        assertTrue("Appending identical reads should increase total depth", secondDepth > firstDepth);
    }

    @Test
    public void testRemoveSecondCutSiteIndexOf() {

        GBSSeqToTagDBPlugin GBSSeqToTagdb = new GBSSeqToTagDBPlugin();
        Class GBSSeqClass = GBSSeqToTagdb.getClass();
        Method removeSecondCutSiteIndexOf;
        try {
            Class[] args = new Class[]{ String.class,Integer.TYPE};
            removeSecondCutSiteIndexOf = GBSSeqClass.getDeclaredMethod("removeSecondCutSiteIndexOf", args);
            removeSecondCutSiteIndexOf.setAccessible(true);
            GBSSeqToTagdb.enzyme("ApeKI");
            GBSSeqToTagdb.likelyReadEndStrings = new String[]{"ATGCAT", "ATGCAAGAT"};
            GBSSeqToTagdb.readEndCutSiteRemnantLength = 5;
            String seq1 = "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATGCAAGATAGGAACAGCGCTAGGGGAATG";
            String seq2 = "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATGCATGATAGGAACAGCGCTAGGGGAATG";
            String expectedTag = "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATGCA";
            // Test string with second cut site = likelyReadEndString[1]
            // ATGCAAGAT at position 41 in string
            Tag tag = (Tag) removeSecondCutSiteIndexOf.invoke(GBSSeqToTagdb, (Object)seq1,55);
            assertTrue(expectedTag.equals(tag.sequence()));

            // Test string with 2nd cut site = likelyReadEndString[0]
            // ATGCAT at position 41 in string
            tag = (Tag) removeSecondCutSiteIndexOf.invoke(GBSSeqToTagdb, (Object)seq2,55);
            assertTrue(expectedTag.equals(tag.sequence()));

            // Test sequence with NO second cut site appearing
            expectedTag = "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATTCATGATAGGAAC";
            String seq3 = "TAGGAACAGCGCTAGGGGAATGCTAAATTGCTAGCGCCATATTCATGATAGGAACAGCGCTAGGGGAATG";
            tag = (Tag) removeSecondCutSiteIndexOf.invoke(GBSSeqToTagdb, (Object)seq3,55);
            assertTrue(expectedTag.equals(tag.sequence()));
        } catch (NoSuchMethodException e) {
            e.printStackTrace();
        } catch (SecurityException e) {
            e.printStackTrace();
        } catch (InvocationTargetException ite) {
            ite.getCause();
            ite.printStackTrace();
        } catch (Exception exc)  {
            exc.printStackTrace();
        }
    }
}
