package net.maizegenetics.dna.map;

import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import static org.junit.Assert.*;
import static org.junit.Assert.assertEquals;

public class GenomeFeatureBuilderTest {

    GenomeFeatureBuilder myBuilder = null;
    GenomeFeature baseFeature = null;

    @Before
    public void setUp() throws Exception {
        baseFeature = new GenomeFeatureBuilder().chromosome(1).id("noID").parentId("noParent").start(1).stop(1)
                .type("noType").build();
    }

    @After
    public void tearDown() throws Exception {

    }

    @Test
    public void testParseGffLine() throws Exception {

    }

    @Test
    public void testAssignmentMethods() throws Exception{
        //Test build from existing feature
        new GenomeFeatureBuilder(baseFeature).build();

        //Test ID assignment
        String myID="testID";
        assertEquals("Wrong self ID assignment", myID, new GenomeFeatureBuilder(baseFeature).id(myID).build().id());

        //Test parent ID assignment
        String parID="parentID";
        assertEquals("Wrong parent ID assignment", parID, new GenomeFeatureBuilder(baseFeature).parentId(parID).build().parentId());

        //Test type assignment
        String myType="testType";
        assertEquals("Wrong feature type assignment", myType, new GenomeFeatureBuilder(baseFeature).type(myType).build().type());

        //Test Chromosome assignment
        Chromosome testChrom = new Chromosome("chr1_testbuild");
        assertEquals("Wrong chromosome assignment (integer)", 1, Integer.parseInt(new GenomeFeatureBuilder(baseFeature).chromosome(1).build().chromosome()));
        assertEquals("Wrong chromosome assignment (string)", "chr1", new GenomeFeatureBuilder(baseFeature).chromosome("chr1").build().chromosome());
        assertEquals("Wrong chromosome assignment (Chromosome class)", testChrom.getName(), new GenomeFeatureBuilder(baseFeature).chromosome(testChrom).build().chromosome());

        //Test start and stop assigned correctly
        int pos=123;
        assertEquals("Wrong start assignment", pos, new GenomeFeatureBuilder(baseFeature).start(pos).stop(pos+1).build().start());
        assertEquals("Wrong stop assignment", pos+1, new GenomeFeatureBuilder(baseFeature).start(pos).stop(pos+1).build().stop());

        //Test that assigning a position also assigns the start and stop
        assertEquals("'Position' not translated to start", pos, new GenomeFeatureBuilder(baseFeature).position(pos).build().start());
        assertEquals("'Position' not translated to stop", pos, new GenomeFeatureBuilder(baseFeature).position(pos).build().start());

        //Test that assigning miscellaneous annotations works, and that trying to get a nonexistant one returns NA
        String testkey="testkey", testvalue="testvalue";
        assertEquals("Misc. annotation not added correctly", testvalue, new GenomeFeatureBuilder(baseFeature).addAnnotation(testkey, testvalue).build().getAnnotation(testkey));
        assertEquals("Unknown annotation should return NA", "NA", new GenomeFeatureBuilder(baseFeature).build().getAnnotation(testkey));
    }

    /**
     * Test if making a copy of a feature protects the original so changes to one don't affect the other
     */
    @Test
    public void testCopyProtection(){
        GenomeFeature newFeature = new GenomeFeatureBuilder(baseFeature).build();
        assertNotSame("Copy protection error on annotations", baseFeature.annotations(), newFeature.annotations());

        newFeature = new GenomeFeatureBuilder(baseFeature).id("newID").type("newType").parentId("newParent")
                .chromosome(99).start(99).stop(100).build();
        assertNotSame("Copy protection error on ID", baseFeature.id(), newFeature.id());
        assertNotSame("Copy protection error on parentID", baseFeature.parentId(), newFeature.parentId());
        assertNotSame("Copy protection error on chromosome", baseFeature.chromosome(), newFeature.chromosome());
        assertNotSame("Copy protection error on start position", baseFeature.start(), newFeature.start());
        assertNotSame("Copy protection error on stop position", baseFeature.stop(), newFeature.stop());
        assertNotSame("Copy protection error on feature type", baseFeature.type(), newFeature.type());

    }

    //Test that giving various synonyms of common annotations in gets the correct standard ones out
 
    @Test
    public void testSynonymizeKeys(){
        myBuilder = new GenomeFeatureBuilder();
        //id
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "id", myBuilder.synonymizeKeys("id"));
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "id", myBuilder.synonymizeKeys("name"));

        //chromosome
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "chromosome", myBuilder.synonymizeKeys("chr"));
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "chromosome", myBuilder.synonymizeKeys("chrom"));
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "chromosome", myBuilder.synonymizeKeys("chromosome"));

        //stop
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "stop", myBuilder.synonymizeKeys("stop"));
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "stop", myBuilder.synonymizeKeys("end"));

        //parent
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "parent_id", myBuilder.synonymizeKeys("parent"));
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "parent_id", myBuilder.synonymizeKeys("parentID"));
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "parent_id", myBuilder.synonymizeKeys("parent_id"));

        //position
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "position", myBuilder.synonymizeKeys("pos"));
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "position", myBuilder.synonymizeKeys("position"));

        //Test that undoing capitalization works
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "id", myBuilder.synonymizeKeys("ID"));

        //Test that values not in the list are returned unchanged/just lowercased
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "random123", myBuilder.synonymizeKeys("random123"));
        assertEquals("GenomeFeatureBuilder.synonymizeKeys error", "random123", myBuilder.synonymizeKeys("RANDOM123"));
    }

    /**
     * Tests if various exceptions are thrown for not having the required fields. While it's possible to construct these with
     * some form of ExpectedException, these use try-catch to (1) keep everything in one method easier, and (2) allow
     * more control over what's tested.
     */
    @Test
    public void testMissingFields(){
        //Test if no valid id entered
        try{
            new GenomeFeatureBuilder().type("notype").start(1).stop(1).chromosome(1).build();
            fail("Uncaught feature ID not assigned");
        }catch(Exception e){
            assertEquals("Wrong exception for unassigned feature ID", UnsupportedOperationException.class, e.getClass());
            assertTrue("Wrong exception message for unassigned feature ID", e.getMessage()
                    .contains("Cannot build a feature without a personal identifier"));
        }

        //Test if no valid type entered
        try{
            new GenomeFeatureBuilder().id("testID").start(1).stop(1).chromosome(1).build();
            fail("Uncaught feature type not assigned");
        }catch(Exception e){
            assertEquals("Wrong exception for unassigned feature type", UnsupportedOperationException.class, e.getClass());
            assertTrue("Wrong exception message for unassigned feature type", e.getMessage()
                    .contains("Cannot build a feature without a feature type"));
        }

        //Test if no valid chromosome entered
        try{
            new GenomeFeatureBuilder().id("testID").type("notype").start(1).stop(1).build();
            fail("Uncaught feature chromosome not assigned");
        }catch(Exception e){
            assertEquals("Wrong exception for unassigned chromosome", UnsupportedOperationException.class, e.getClass());
            assertTrue("Wrong exception message for unassigned chromosome", e.getMessage()
                    .contains("Cannot build a feature without a chromosome"));
        }

        //Test if no valid start position entered
        try{
            new GenomeFeatureBuilder().id("testID").type("notype").stop(1).chromosome(1).build();
            fail("Uncaught feature start position not assigned");
        }catch(Exception e){
            assertEquals("Wrong exception for unassigned start position", UnsupportedOperationException.class, e.getClass());
            assertTrue("Wrong exception message for unassigned start position", e.getMessage()
                    .contains("Cannot build a feature without a start position"));
        }

        //Test if no valid stop position entered
        try{
            new GenomeFeatureBuilder().id("testID").type("notype").start(1).chromosome(1).build();
            fail("Uncaught feature stop position not assigned");
        }catch(Exception e){
            assertEquals("Wrong exception for unassigned stop position", UnsupportedOperationException.class, e.getClass());
            assertTrue("Wrong exception message for unassigned stop position", e.getMessage()
                    .contains("Cannot build a feature without a stop position"));
        }

    }

    /**
     * Tests if various exceptions are thrown when they should. While it's possible to construct these with
     * some form of ExpectedException, these use try-catch to (1) keep everything in one method easier, and (2) allow
     * more control over what's tested.
     */
    @Test
    public void testExceptionsOfInvalidPositions(){


        //Test if start before stop
        try{
            new GenomeFeatureBuilder().id("testID").type("notype").start(2).stop(1).chromosome(1).build();
            fail("Uncaught stop position before start");
        }catch(Exception e){
            assertEquals("Wrong exception for stop before start positions", UnsupportedOperationException.class, e.getClass());
            assertTrue("Wrong exception message for stop before start positions", e.getMessage()
                    .contains("Start coordinate is greater than stop coordinate"));
        }

        //Test if start < 0
        try{
            new GenomeFeatureBuilder().id("testID").type("notype").start(-1).stop(1).chromosome(1).build();
            fail("Uncaught negative start position on build");
        }catch(Exception e){
            assertEquals("Wrong/missed exception for negative start value on build", UnsupportedOperationException.class, e.getClass());
            assertTrue("Wrong exception message for negative start value on build", e.getMessage()
                    .contains("coordinate is negative"));
        }

        //Test if stop < 0
        try{
            new GenomeFeatureBuilder().id("testID").type("notype").start(1).stop(-1).chromosome(1).build();
            fail("Uncaught negative stop position on build");
        }catch(Exception e){
            assertEquals("Wrong/missed exception for negative stop value on build", UnsupportedOperationException.class, e.getClass());
            assertTrue("Wrong exception message for negative stop value on build", e.getMessage()
                    .contains("coordinate is negative"));
        }
    }

    /*@Test
    public void testGetParentFromGffAttributes() throws Exception {
        myBuilder = new GenomeFeatureBuilder();
        String parent="TestGene123";
        String decoy = "WrongGene999";
        String noise = "exon_number \"19\"; seqedit \"false\"; protein_id \"GRMZM2G023858_P03\";";

        //Test that empty string returned if none of the appropriate flags present
        assertEquals("Empty list:", "", myBuilder.getParentFromGffAttributes(""));
        assertEquals("Parent without flag:", "", myBuilder.getParentFromGffAttributes(parent));
        assertEquals("Only different flags:", "", myBuilder.getParentFromGffAttributes(noise));

        //Test that finds 'Parent=GRMZM2G005232'
        assertEquals("Test1 for parse Parent=", parent, myBuilder.getParentFromGffAttributes("Parent=" + parent));
        assertEquals("Test2 for parse Parent=", parent, myBuilder.getParentFromGffAttributes("Parent=" + parent + "; " + noise));
        assertEquals("Test3 for parse Parent=", parent, myBuilder.getParentFromGffAttributes(noise + " Parent=" + parent));

        //Test that finds 'Parent=gene:GRMZM2G005232' properly ('gene' can be any identifier)
        assertEquals("Test1 for parse Parent=gene:", parent, myBuilder.getParentFromGffAttributes("Parent=gene:" + parent));
        assertEquals("Test2 for parse Parent=gene:",parent, myBuilder.getParentFromGffAttributes("Parent=transcript:" + parent));
        assertEquals("Test3 for parse Parent=gene:",parent, myBuilder.getParentFromGffAttributes("Parent=cds:" + parent));
        assertEquals("Test4 for parse Parent=gene:",parent, myBuilder.getParentFromGffAttributes("Parent=gobbeleygook:" + parent));
        assertEquals("Test5 for parse Parent=gene:",parent, myBuilder.getParentFromGffAttributes("Parent=gene:" + parent + "; " +  noise));
        assertEquals("Test6 for parse Parent=gene:",parent, myBuilder.getParentFromGffAttributes(noise + " Parent=gene:" + parent));

        //Test that finds 'parent_id "GRMZM2G005232"' properly
        assertEquals("Test1 for parse parent_id", parent, myBuilder.getParentFromGffAttributes("parent_id \"" + parent + "\""));
        assertEquals("Test2 for parse parent_id", parent, myBuilder.getParentFromGffAttributes("parent_id \"" + parent + "\"; " + noise));
        assertEquals("Test3 for parse parent_id", parent, myBuilder.getParentFromGffAttributes(noise + " parent_id \"" + parent + "\"; "));

        //Test that finds 'transcript_id "GRMZM2G005232_T01"' properly
        assertEquals("Test1 for parse transcript_id", parent, myBuilder.getParentFromGffAttributes("transcript_id \"" + parent + "\""));
        assertEquals("Test2 for parse transcript_id", parent, myBuilder.getParentFromGffAttributes("transcript_id \"" + parent + "\"; " + noise));
        assertEquals("Test3 for parse transcript_id", parent, myBuilder.getParentFromGffAttributes(noise + " transcript_id \"" + parent + "\"; "));

        //Test that finds 'gene_id "GRMZM2G005232"' properly
        assertEquals("Test1 for parse gene_id", parent, myBuilder.getParentFromGffAttributes("gene_id \"" + parent + "\""));
        assertEquals("Test2 for parse gene_id", parent, myBuilder.getParentFromGffAttributes("gene_id \"" + parent + "\"; " + noise));
        assertEquals("Test3 for parse gene_id", parent, myBuilder.getParentFromGffAttributes(noise + " gene_id \"" + parent + "\"; "));

        //Test that hierarchy works so that 'parent' > 'transcript' > 'gene'
        String realgene = "gene_id \"" + parent + "\"; ";
        String realtrans = "transcript_id \"" + parent + "\"; ";
        String realparent1 = "parent_id \"" + parent + "\"; ";
        String realparent2 = "Parent=" + parent + "; ";
        String decoygene = "gene_id \"" + decoy + "\"; ";
        String decoytrans = "transcript_id \"" + decoy + "\"; ";
        assertEquals("Hierarchy test 1", parent, myBuilder.getParentFromGffAttributes(realgene + noise));
        assertEquals("Hierarchy test 2", parent, myBuilder.getParentFromGffAttributes(realtrans + decoygene));
        assertEquals("Hierarchy test 3", parent, myBuilder.getParentFromGffAttributes(realparent1 + decoytrans));
        assertEquals("Hierarchy test 4", parent, myBuilder.getParentFromGffAttributes(realparent1 + decoygene));
        assertEquals("Hierarchy test 5", parent, myBuilder.getParentFromGffAttributes(realparent2 + decoytrans));
        assertEquals("Hierarchy test 6", parent, myBuilder.getParentFromGffAttributes(realparent2 + decoygene));
    }*/

    /*@Test
    public void testGetFeatureIdFromGffAttributes() throws Exception {
        myBuilder = new GenomeFeatureBuilder();
        String id="TestGene123";
        String noise = "exon_number \"19\"; seqedit \"false\"; protein_id \"GRMZM2G023858_P03\"; ;Parent=transcript:GRMZM5G856777_T01;rank=1;";

        //Test that returns null when appropriate flags not present
        assertEquals("Null test 1", null, myBuilder.getFeatureIdFromGffAttributes(""));
        assertEquals("Null test 2", null, myBuilder.getFeatureIdFromGffAttributes(id));
        assertEquals("Null test 3", null, myBuilder.getFeatureIdFromGffAttributes("Parent=" + id));

        //Test that finds ID in pattern ID=
        assertEquals("ID= test 1", id, myBuilder.getFeatureIdFromGffAttributes("ID=" + id));
        assertEquals("ID= test 2", id, myBuilder.getFeatureIdFromGffAttributes("ID=gene:" + id));
        assertEquals("ID= test 3", id, myBuilder.getFeatureIdFromGffAttributes("ID=gene:" + id + ";" + noise));
        assertEquals("ID= test 4", id, myBuilder.getFeatureIdFromGffAttributes(noise + "ID=gene:" + id));

        //Test that finds ID in pattern Name=
        assertEquals("Name= test 1", id, myBuilder.getFeatureIdFromGffAttributes("Name=" + id));
        assertEquals("Name= test 2", id, myBuilder.getFeatureIdFromGffAttributes("Name=gene:" + id));
        assertEquals("Name= test 3", id, myBuilder.getFeatureIdFromGffAttributes("Name=gene:" + id + ";" + noise));
        assertEquals("Name= test 4", id, myBuilder.getFeatureIdFromGffAttributes(noise + "Name=gene:" + id));
    }*/
}