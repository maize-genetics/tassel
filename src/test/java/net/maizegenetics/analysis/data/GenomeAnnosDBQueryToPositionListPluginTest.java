/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.constants.GeneralConstants;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import static org.junit.Assert.*;

/**
 * Running a query needs a live PostgreSQL genome annotations DB, so these tests
 * cover query checking and config-file validation, which happen before connecting.
 *
 * @author jcg233
 */
public class GenomeAnnosDBQueryToPositionListPluginTest {

    private static final Path TEMP_DIR = Paths.get(GeneralConstants.TEMP_DIR, "GenomeAnnosDBQueryToPositionListPluginTest");

    /**
     * Test of checkQuery method, of class GenomeAnnosDBQueryToPositionListPlugin.
     */
    @Test
    public void testCheckQuery() {
        String checkedQuery = GenomeAnnosDBQueryToPositionListPlugin.checkQuery("SELECT chr, position alter where");
        assertNull(checkedQuery); 
        checkedQuery = GenomeAnnosDBQueryToPositionListPlugin.checkQuery("SELECT chr, position copy where");
        assertNull(checkedQuery); 
        checkedQuery = GenomeAnnosDBQueryToPositionListPlugin.checkQuery("SELECT chr, position create where");
        assertNull(checkedQuery); 
        checkedQuery = GenomeAnnosDBQueryToPositionListPlugin.checkQuery("SELECT chr, position delete where");
        assertNull(checkedQuery); 
        checkedQuery = GenomeAnnosDBQueryToPositionListPlugin.checkQuery("SELECT chr, position drop where");
        assertNull(checkedQuery); 
        checkedQuery = GenomeAnnosDBQueryToPositionListPlugin.checkQuery("SELECT chr, position insert where");
        assertNull(checkedQuery); 
        checkedQuery = GenomeAnnosDBQueryToPositionListPlugin.checkQuery("SELECT chr, position truncate where");
        assertNull(checkedQuery); 
        checkedQuery = GenomeAnnosDBQueryToPositionListPlugin.checkQuery("SELECT chr, position update where");
        assertNull(checkedQuery); 

        checkedQuery = GenomeAnnosDBQueryToPositionListPlugin.checkQuery("SELECT chr, position where");
        assertNotNull(checkedQuery); 
        checkedQuery = GenomeAnnosDBQueryToPositionListPlugin.checkQuery("SELECT chr, position, tear_drop where");
        assertNotNull(checkedQuery); 
        checkedQuery = GenomeAnnosDBQueryToPositionListPlugin.checkQuery("SELECT chr, position, kmer_copynumber where");
        assertNotNull(checkedQuery); 
    }

    @Test
    public void testCheckQueryRequiresSelectChrAndPosition() {
        assertNull(GenomeAnnosDBQueryToPositionListPlugin.checkQuery("chr, position from annos"));
        assertNull(GenomeAnnosDBQueryToPositionListPlugin.checkQuery("SELECT position from annos"));
        assertNull(GenomeAnnosDBQueryToPositionListPlugin.checkQuery("SELECT chr from annos"));
        assertEquals("select chr, position from annos",
                GenomeAnnosDBQueryToPositionListPlugin.checkQuery("  select chr, position from annos \n"));
    }

    @Test
    public void testConnectToDBRejectsIncompleteConfig() throws IOException {
        Files.createDirectories(TEMP_DIR);
        assertNull("missing user", GenomeAnnosDBQueryToPositionListPlugin.connectToDB(
                writeConfig("noUser.txt", "password=secret\nDB=maizeGenomeAnnos\n")));
        assertNull("missing password", GenomeAnnosDBQueryToPositionListPlugin.connectToDB(
                writeConfig("noPassword.txt", "user=postgres\nDB=maizeGenomeAnnos\n")));
        assertNull("missing DB", GenomeAnnosDBQueryToPositionListPlugin.connectToDB(
                writeConfig("noDB.txt", "user=postgres\npassword=secret\n")));
        assertNull("missing file", GenomeAnnosDBQueryToPositionListPlugin.connectToDB(
                TEMP_DIR.resolve("doesNotExist.txt").toString()));
    }

    private static String writeConfig(String name, String contents) throws IOException {
        Path config = TEMP_DIR.resolve(name);
        Files.write(config, contents.getBytes());
        return config.toString();
    }

}
