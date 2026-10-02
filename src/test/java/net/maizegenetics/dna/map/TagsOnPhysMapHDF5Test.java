/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.dna.map;

import java.io.File;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.constants.GeneralConstants;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Test for {@link TagsOnPhysMapHDF5}.
 *
 * <p>Originally this test read multi-GB TOPM fixtures from hard-coded {@code /Volumes/LaCie/...}
 * paths that are not part of the test-data release, and asserted chromosome sets tied to the old
 * AllZeaGBS build. It now derives its fixtures at runtime from the present binary TOPM
 * ({@link GBSConstants#GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE}): it writes a text
 * TOPM and an HDF5 TOPM under the git-ignored temp dir and asserts the round-trips agree. The
 * HDF5 write/read path exercises the native jhdf5 library.</p>
 *
 * @author edbuckler
 */
public class TagsOnPhysMapHDF5Test {

    private static final String TEST_DIR = GeneralConstants.TEMP_DIR + "TagsOnPhysMapHDF5Test/";
    private static final String sourceBinaryTOPM = GBSConstants.GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE;
    private static final String textFile = TEST_DIR + "test.topm.txt";
    private static final String hdf5File = TEST_DIR + "test.topm.h5";

    private static TagsOnPhysMapHDF5 instance;
    private static int[] expectedChromosomes;

    public TagsOnPhysMapHDF5Test() {
    }

    @BeforeClass
    public static void setUpClass() {
        new File(TEST_DIR).mkdirs();
        long time = System.currentTimeMillis();

        // Read the present binary TOPM fixture and derive text + HDF5 TOPMs from it.
        TagsOnPhysicalMap sourceBin = new TagsOnPhysicalMap(sourceBinaryTOPM, true);
        expectedChromosomes = sourceBin.getChromosomes();
        sourceBin.writeTextFile(new File(textFile));
        TagsOnPhysMapHDF5.createFile(sourceBin, hdf5File, 4, 16);

        instance = new TagsOnPhysMapHDF5(hdf5File);
        System.out.println("TagsOnPhysMapHDF5Test setUpClass took " + (System.currentTimeMillis() - time) + " ms");
    }

    /**
     * Round-trips the derived text TOPM through binary and HDF5 and confirms the fields agree.
     */
    @Test
    public void testCompareTOPM() {
        System.out.println("Comparing a range of TOPM against one another");
        TagsOnPhysicalMap theText = new TagsOnPhysicalMap(textFile, false);
        theText.writeBinaryFile(new File(textFile.replace(".txt", ".bin")));
        TagsOnPhysicalMap theTextToBin = new TagsOnPhysicalMap(textFile.replace(".txt", ".bin"), true);
        compareTOPMFields(theText, theTextToBin);
        assertEquals(theText.getTagCount(), theTextToBin.getTagCount());

        // Distinct from the HDF5 file setUpClass keeps open as {@code instance}, so this
        // createFile does not try to truncate a file that is already open.
        String hdf5FromText = TEST_DIR + "test_fromText.topm.h5";
        TagsOnPhysMapHDF5.createFile(theText, hdf5FromText, 4, 16);
        TagsOnPhysMapHDF5 theNewHDF5 = new TagsOnPhysMapHDF5(hdf5FromText);
        compareTOPMFields(theText, theNewHDF5);
        theNewHDF5.writeTextFile(new File(textFile.replace(".txt", ".hdf2.txt")));
    }

    public static void compareTOPMFields(AbstractTagsOnPhysicalMap theOld, AbstractTagsOnPhysicalMap theNew) {
        System.out.println("Comparing two TOPM");
        assertEquals(theOld.getTagCount(), theNew.getTagCount());
        for (int i = 0; i < theOld.getTagCount(); i+=1) {
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
     * The chromosomes reported by the HDF5 TOPM must match those in the source binary TOPM.
     */
    @Test
    public void testGetChromosomes() {
        System.out.println("getChromosomes");
        int[] result = instance.getChromosomes();
        assertArrayEquals(expectedChromosomes, result);
    }

}
