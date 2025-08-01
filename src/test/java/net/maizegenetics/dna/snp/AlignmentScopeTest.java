/*
 * AlignmentScopeTest
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import org.junit.*;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.fail;

/**
 *
 * @author terry
 */
public class AlignmentScopeTest {

    public AlignmentScopeTest() {
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

    @Test
    public void testAlignmentScopeMethods() throws IOException {
        System.out.println("Testing Alignment Scope Methods...");

        String expectedResultsFilename = GeneralConstants.EXPECTED_RESULTS_DIR + "FilterAlignmentTest.hmp.txt";
        GenotypeTable expectedAlign = ImportUtils.readFromHapmap(expectedResultsFilename, null);

        File input = new File(TutorialConstants.HAPMAP_FILENAME);
        System.out.println("   Input File: " + input.getCanonicalPath());
        if (!input.exists()) {
            fail("Input File: " + TutorialConstants.HAPMAP_FILENAME + " doesn't exist.");
        }
        GenotypeTable inputAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
        System.out.println("   Input Alignment: Taxa Count: " + inputAlign.numberOfTaxa() + "  Site Count: " + inputAlign.numberOfSites());

        TaxaList origTaxa = inputAlign.taxa();
        Taxon[] newIds = new Taxon[origTaxa.numberOfTaxa() / 3 + 1];
        int count = 0;
        for (int i = 0; i < origTaxa.numberOfTaxa(); i = i + 3) {
            newIds[count++] = origTaxa.get(i);
        }
        TaxaList taxa = new TaxaListBuilder().addAll(newIds).build();
        GenotypeTable filterTaxa = FilterGenotypeTable.getInstance(inputAlign, taxa);
        System.out.println("   Filtered Taxa Alignment: Taxa Count: " + filterTaxa.numberOfTaxa() + "  Site Count: " + filterTaxa.numberOfSites());

        double minFreq = 0.1;
        double maxFreq = 0.8;
        int minCount = 0;
        GenotypeTable filterTaxaSites = GenotypeTableUtils.removeSitesBasedOnFreqIgnoreMissing(filterTaxa, minFreq, maxFreq, minCount);
        System.out.println("   Filtered Taxa and Sites Alignment: Taxa Count: " + filterTaxaSites.numberOfTaxa() + "  Site Count: " + filterTaxaSites.numberOfSites());

        AlignmentTestingUtils.alignmentsScopeEqual(GenotypeTable.ALLELE_SORT_TYPE.Frequency, expectedAlign, filterTaxaSites);
    }
}
