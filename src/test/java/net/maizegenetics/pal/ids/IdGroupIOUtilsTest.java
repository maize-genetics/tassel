/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.ids;

import net.maizegenetics.constants.GeneralConstants;
import org.junit.*;

/**
 *
 * @author edbuckler
 */
public class IdGroupIOUtilsTest {
    String infile=GeneralConstants.DATA_DIR+"CandidateTests/AnnotationForGoodman282.txt";
    
    public IdGroupIOUtilsTest() {
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
     * Test of readPedigree method, of class TaxaListIOUtils.
     */
    @Test
    public void testReadPedigree() {
//        System.out.println("readPedigree");
//        TaxaList result = TaxaListIOUtils.readPedigree(infile);
//     //   System.out.println(result.toString());
//        for (int i = 0; i < result.getIdCount(); i++) {
//           AnnotatedIdentifier ai=(AnnotatedIdentifier)result.getIdentifier(i);
//           if(ai.getTextAnnotation("Slice").length>0) System.out.println(ai.toString());
//        }
//        TreeMultimap<String,String> memberOfSets=TaxaListIOUtils.getMapOfTextAnnotatedIds(result, "Slice");
//        for (String annoValue: memberOfSets.keySet()) {
//            System.out.println(annoValue+":"+memberOfSets.get(annoValue).size());
//        }
//        //System.out.println(memberOfSets);
//
////        assertEquals(expResult, result);
////        // TODO review the generated test code and remove the default call to fail.
////        fail("The test case is a prototype.");
    }
}