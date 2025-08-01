/*
 *  GenotypeTableBuilderTest
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author Terry Casstevens
 */
public class GenotypeTableBuilderTest {

    public GenotypeTableBuilderTest() {
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
     * Test of getHomozygousInstance method, of class GenotypeTableBuilder.
     */
    @Test
    public void testGetHomozygousInstance() {

        System.out.println("Testing GenotypeTableBuilder.getHomozygousInstance()...");

        GenotypeTable inputAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
        //GenotypeTable inputAlign = ImportUtils.readFromHapmap("/Users/terry/data/B73RefGenV2_1/maizeHapMapV2_B73RefGenV2_201203028_chr10.hmp.txt.gz", null);
        System.out.println("   Input Alignment: Taxa Count: " + inputAlign.numberOfTaxa() + "  Site Count: " + inputAlign.numberOfSites());

        long previous = System.nanoTime();
        GenotypeTable result = GenotypeTableBuilder.getHomozygousInstance(inputAlign);
        long current = System.nanoTime();
        System.out.println("   getHomozygousInstance time: " + ((double) (current - previous) / 1_000_000_000.0));

        String expectedResultsFilename = GeneralConstants.EXPECTED_RESULTS_DIR + "HomozygousGenotypeTable.hmp.txt";
        GenotypeTable expectedAlign = ImportUtils.readFromHapmap(expectedResultsFilename, null);

        AlignmentTestingUtils.alignmentsEqual(expectedAlign, result);

    }

}
