/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.dna.map;

import java.io.File;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 * Test for TagsOnPhysMapHDF5
 * <p>
 * TODO:  Do a comparison of the speed to load, compare against so default read speed
 * TODO:  Develop a test HDF5 the only includes the first 20Mb of chr 9 and 10.
 * TODO:  Test creation of this class from Tag
 * 
 * 
 * @author edbuckler
 */
public class TagsOnPhysMapHDF5Test {
    static String root="/Volumes/LaCie/build20120701/04_TOPM/";
//    static String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/04_TOPM/v26/";
//    static String inputFile="/Volumes/LaCie/build20120701/04_TOPM/AllZeaGBS_v2.6_MergedUnfiltProdTOPM_20130515.topm.h5";
    static String realFile=root+"AllZeaGBS_v2.6_MergedUnfiltProdTOPM_20130425.topm";
    static String realOutFile=root+"AllZeaGBSProdTOPM_20130603.topm.h5";
    static String inputFile=root+"SmallTOPM_20130515.topm.h5";
    static String textFile=root+"head100ktest.topm.txt";
//    static String inputFile=root+"FullTOPM_20130515.topm.h5";
//    static String outputFile="/Volumes/LaCie/build20120701/04_TOPM/TestTOPM_20130515.topm.h5";
    static String outputFile=root+"TestTOPM_20130515.topm.h5";
    
    static TagsOnPhysMapHDF5 instance;
    
    public TagsOnPhysMapHDF5Test() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        long time=System.currentTimeMillis();
        instance=(TagsOnPhysMapHDF5)TOPMUtils.readTOPM(inputFile);
        
//        AbstractTagsOnPhysicalMap instance2=(AbstractTagsOnPhysicalMap)TOPMUtils.readTOPM(realFile);
//        TagsOnPhysMapHDF5.createFile(instance2, realOutFile, 4, 16);
        
//        TagsOnPhysMapHDF5 newBig=new TagsOnPhysMapHDF5(realOutFile,false);
//        for (int i = 0; i < 10000; i++) {
//            System.out.println(newBig.printRow(i));
//            
//        }
//        System.exit(0);
        
        System.out.println(System.currentTimeMillis()-time);
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

//    /**
//     * Test of createFile method, of class TagsOnPhysMapHDF5.
//     */
//    @Test
//    public void testCreateFile() {
//        System.out.println("createFile");
//        TagsOnPhysicalMap inTags = null;
//        String newHDF5file = "";
//        int maxMapping = 0;
//        int maxVariants = 0;
//        TagsOnPhysMapHDF5.createFile(inTags, newHDF5file, maxMapping, maxVariants);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getFileReadyForClosing method, of class TagsOnPhysMapHDF5.
//     */
//    @Test
//    public void testGetFileReadyForClosing() {
//        System.out.println("getFileReadyForClosing");
//        TagsOnPhysMapHDF5 instance = null;
//        instance.getFileReadyForClosing();
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getAlternateTagMappingInfo method, of class TagsOnPhysMapHDF5.
//     */
//    @Test
//    public void testGetAlternateTagMappingInfo() {
//        System.out.println("getAlternateTagMappingInfo");
//        int index = 0;
//        int mapIndex = 0;
//        TagsOnPhysMapHDF5 instance = null;
//        TagMappingInfo expResult = null;
//        TagMappingInfo result = instance.getAlternateTagMappingInfo(index, mapIndex);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setAlternateTagMappingInfo method, of class TagsOnPhysMapHDF5.
//     */
//    @Test
//    public void testSetAlternateTagMappingInfo() {
//        System.out.println("setAlternateTagMappingInfo");
//        int index = 0;
//        int mapIndex = 0;
//        TagMappingInfo theTMI = null;
//        TagsOnPhysMapHDF5 instance = null;
//        instance.setAlternateTagMappingInfo(index, mapIndex, theTMI);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of swapTagMappingInfo method, of class TagsOnPhysMapHDF5.
//     */
//    @Test
//    public void testSwapTagMappingInfo() {
//        System.out.println("swapTagMappingInfo");
//        int index = 0;
//        int mapIndex = 0;
//        int mapIndex2 = 0;
//        TagsOnPhysMapHDF5 instance = null;
//        instance.swapTagMappingInfo(index, mapIndex, mapIndex2);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of addVariant method, of class TagsOnPhysMapHDF5.
//     */
//    @Test
//    public void testAddVariant() {
//        System.out.println("addVariant");
//        int tagIndex = 0;
//        byte offset = 0;
//        byte base = 0;
//        TagsOnPhysMapHDF5 instance = null;
//        int expResult = 0;
//        int result = instance.addVariant(tagIndex, offset, base);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }

//    /**
//     * Test of getChromosome method, of class TagsOnPhysMapHDF5.
//     */
//    @Test
//    public void testGetChromosome() {
//        System.out.println("getChromosome");
//        int index = 0;
//        TagsOnPhysMapHDF5 instance = null;
//        int expResult = 0;
//        int result = instance.getChromosome(index);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }

    /**
     * Test of getChromosomes method, of class TagsOnPhysMapHDF5.
     */
    @Ignore
    @Test
    public void testConversionChromosomes() {
        long time=System.currentTimeMillis();
        TOPMInterface basicTOPM=TOPMUtils.readTOPM(inputFile);
        System.out.println("Writing test TOPM");
//        TOPMUtils.writeTOPM(basicTOPM, outputFile);
        TOPMUtils.writeTOPM(basicTOPM, outputFile.replace(".h5", ".txt"));
    }
    
        /**
     * Test of many of the get method of class TagsOnPhysMapHDF5.
     */
    @Test
    public void testCompareTOPM() {
//        long time=System.currentTimeMillis();
//        TagsOnPhysicalMap theOld1=new TagsOnPhysicalMap(realFile,true);  //creates a text test file
//        theOld1.writeTextFile(new File(root+"head100ktest.topm.txt"));
//        if(theOld1!=null) System.exit(0);
        System.out.println("Comparing a range of TOPM against one another");
        TagsOnPhysicalMap theText=new TagsOnPhysicalMap(textFile,false);
        theText.writeBinaryFile(new File(textFile.replace(".txt", ".bin")));
        TagsOnPhysicalMap theTextToBin=new TagsOnPhysicalMap(textFile.replace(".txt", ".bin"),true);
        compareTOPMFields(theText, theTextToBin);
        assertEquals(theText.getTagCount(), theTextToBin.getTagCount());
        String hdf5File=textFile.replace(".txt", ".h5");
        TagsOnPhysMapHDF5.createFile(theText, hdf5File, 4, 16);
        TagsOnPhysMapHDF5 theNewHDF5=new TagsOnPhysMapHDF5(hdf5File);
        compareTOPMFields(theText, theNewHDF5);
        theNewHDF5.writeTextFile(new File(textFile.replace(".txt", ".hdf2.txt")));
//        String line;
//        while ((line = expected.readLine()) != null) {
//          assertEquals(line, actual.readLine());
//        }
    }
    
    
    
    
    public static void compareTOPMFields(AbstractTagsOnPhysicalMap theOld, AbstractTagsOnPhysicalMap theNew) {
        System.out.println("Comparing two TOPM");
        assertEquals(theOld.getTagCount(), theNew.getTagCount());
        for (int i = 0; i < theOld.getTagCount(); i+=1) {
//            System.out.println(i);
//            System.out.println(theOld.printRow(i));
//            System.out.println(theNew.printRow(i));
////            System.out.println("Tag:"+i);
            assertEquals(theOld.getMultiMaps(i), theNew.getMultiMaps(i));
            assertEquals(theOld.getChromosome(i), theNew.getChromosome(i));
            assertEquals(theOld.getStartPosition(i), theNew.getStartPosition(i));
            assertEquals(theOld.getStrand(i), theNew.getStrand(i));
            assertEquals(theOld.getTagLength(i), theNew.getTagLength(i));
            assertEquals(theOld.getEndPosition(i), theNew.getEndPosition(i));
            assertEquals(theOld.getReadIndexForPositionIndex(i), theNew.getReadIndexForPositionIndex(i));
        }
    }
    
         /**
     * Test of many of the get method of class TagsOnPhysMapHDF5.
     */
    @Ignore@Test
    public void testCompareBigTOPM() {
        System.out.println("Comparing a giant TOPM against one another");
 //       TagsOnPhysMapHDF5 oldBig=new TagsOnPhysMapHDF5("/Volumes/LaCie/build20120701/04_TOPM/AllZeaGBS_v2.6_MergedUnfiltProdTOPM_20130515.topm.h5",false);
        TagsOnPhysicalMap oldBig=new TagsOnPhysicalMap(realFile,true);
        TagsOnPhysMapHDF5 newBig=new TagsOnPhysMapHDF5("/Volumes/LaCie/build20120701/04_TOPM/AllZeaGBSProdTOPM_20130603.topm.h5",false);
        compareTOPMFields(oldBig, newBig);


    }
    
    /**
     * Test of getChromosomes method, of class TagsOnPhysMapHDF5.
     */
 //   @Ignore
    @Test
    public void testGetChromosomes() {
        System.out.println("getChromosomes");
        int[] expResult = {0,1,2,3,4,5,6,7,8,9,10,11,12};
        int[] result = instance.getChromosomes();
        assertArrayEquals(expResult, result);
        int[] chr1pos=instance.getUniquePositions(1);
//        System.out.printf("Chr 1 positions: %d ... %d %n",chr1pos[0],chr1pos[chr1pos.length-1]);       
    }
    
    
    
    @Ignore@Test
    public void testLoadSpeed() {
        long time=System.currentTimeMillis();
        System.out.println("File Loaded");
//        TOPMUtils.writeTOPM(instance, outputFile);
//        System.out.println(System.currentTimeMillis()-time);
//        time=System.currentTimeMillis();
//        TagsOnPhysMapHDF5 inst2=(TagsOnPhysMapHDF5)TOPMUtils.readTOPM(outputFile);
//        System.out.println(System.currentTimeMillis()-time);
//        System.out.println("Done writing TOPM");
//        System.out.println("getChromosomes");
//        int[] expResult = {0,1,2,3,4,5,6,7,8,9,10,11,12};
//        System.out.println("Try1");
//        int[] result = instance.getChromosomes();
//        assertArrayEquals(expResult, result);
//        System.out.println("Try2");
//        result = instance.getChromosomes();
//        assertArrayEquals(expResult, result);
    }

}