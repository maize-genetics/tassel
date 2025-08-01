/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.analysis.data;

import java.util.Map;
import net.maizegenetics.dna.map.PositionList;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author jcg233
 */
public class GenomeAnnosDBQueryToPositionListPluginTest {
    
    public GenomeAnnosDBQueryToPositionListPluginTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

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
    public void testPlugin() {
        PositionList positList = 
            new GenomeAnnosDBQueryToPositionListPlugin()
                .connConfigFile("/Users/jcg233/Documents/GenomeAnnotations/genomeAnnosConnectionConfigCBSU.txt")
                .queryFile("/Users/jcg233/Documents/GenomeAnnotations/high_impact_high_maf_query.sql")
                .runPlugin(null);
        System.out.println("\nPostionList:\n");
        System.out.println("chr\tposition\tannotations");
        positList.stream().forEach( 
            posit -> {
                System.out.print(posit.getChromosome()+"\t");
                System.out.print(posit.getPosition()+"\t");
                int entryCount = 0;
                for (Map.Entry<String,String> entry : posit.getAnnotation().getAllAnnotationEntries()) {
                    entryCount++;
                    if (entryCount > 1) {
                        System.out.print(";");
                    }
                    System.out.print(entry.toString());
                }
               System.out.print("\n");
            } 
        );
        assertEquals("Wrong number of positions", 820, positList.size());
    }
    
//    @Test
//    public void testPlugin() {
//        PositionList positList = 
//            new GenomeAnnosDBQueryToPositionListPlugin()
//                .connConfigFile("/Users/jcg233/Documents/GenomeAnnotations/genomeAnnosConnectionConfig.txt")
//                .queryFile("/Users/jcg233/Documents/GenomeAnnotations/height_above_ear_rmip_query.sql")
//                .runPlugin(null);
//        System.out.println("\nPostionList:\n");
//        System.out.println("chr\tposition\tannotations");
//        positList.stream().forEach( 
//            posit -> {
//                System.out.print(posit.getChromosome()+"\t");
//                System.out.print(posit.getPosition()+"\t");
//                int entryCount = 0;
//                for (Map.Entry<String,String> entry : posit.getAnnotation().getAllAnnotationEntries()) {
//                    entryCount++;
//                    if (entryCount > 1) {
//                        System.out.print(";");
//                    }
//                    System.out.print(entry.toString());
//                }
//               System.out.print("\n");
//            } 
//        );
//        assertEquals("Wrong number of positions", 820, positList.size());
//    }

//    /**
//     * Test of processData method, of class GenomeAnnosDBQueryToPositionListPlugin.
//     */
//    @Test
//    public void testProcessData() {
//    }
//
//    /**
//     * Test of getToolTipText method, of class GenomeAnnosDBQueryToPositionListPlugin.
//     */
//    @Test
//    public void testGetToolTipText() {
//    }
//
//    /**
//     * Test of getIcon method, of class GenomeAnnosDBQueryToPositionListPlugin.
//     */
//    @Test
//    public void testGetIcon() {
//    }
//
//    /**
//     * Test of getButtonName method, of class GenomeAnnosDBQueryToPositionListPlugin.
//     */
//    @Test
//    public void testGetButtonName() {
//    }
//
//    /**
//     * Test of runPlugin method, of class GenomeAnnosDBQueryToPositionListPlugin.
//     */
//    @Test
//    public void testRunPlugin() {
//    }
//
//    /**
//     * Test of connConfigFile method, of class GenomeAnnosDBQueryToPositionListPlugin.
//     */
//    @Test
//    public void testConnConfigFile_0args() {
//    }
//
//    /**
//     * Test of connConfigFile method, of class GenomeAnnosDBQueryToPositionListPlugin.
//     */
//    @Test
//    public void testConnConfigFile_String() {
//    }
//
//    /**
//     * Test of queryFile method, of class GenomeAnnosDBQueryToPositionListPlugin.
//     */
//    @Test
//    public void testQueryFile_0args() {
//    }
//
//    /**
//     * Test of queryFile method, of class GenomeAnnosDBQueryToPositionListPlugin.
//     */
//    @Test
//    public void testQueryFile_String() {
//    }
    
}
