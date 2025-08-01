/*
 * FilterGenotypeTableTest
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.analysis.data.UnionAlignmentPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.fail;

/**
 * @author Terry Casstevens
 */
public class FilterGenotypeTableTest {

    public FilterGenotypeTableTest() {
    }

    @Test
    public void testFilterTaxaAndSites() throws IOException {
        System.out.println("Testing Filtering Taxa and Sites...");

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

        double minFreq = 0.1;
        double maxFreq = 0.8;
        int minCount = 0;
        GenotypeTable filterTaxaSites = GenotypeTableUtils.removeSitesBasedOnFreqIgnoreMissing(filterTaxa, minFreq, maxFreq, minCount);
        System.out.println("   Filtered Alignment: Taxa Count: " + filterTaxaSites.numberOfTaxa() + "  Site Count: " + filterTaxaSites.numberOfSites());

        GenotypeTable copyOfFilter = GenotypeTableBuilder.getGenotypeCopyInstance(filterTaxaSites);

        AlignmentTestingUtils.alignmentsEqual(expectedAlign, filterTaxaSites);
        AlignmentTestingUtils.alignmentsEqual(expectedAlign, copyOfFilter);
    }

    @Test
    public void testFilterVCF() {

        System.out.println("Testing Filtering from VCF File");

        String input = TutorialConstants.TUTORIAL_DIR + "mdp_genotype_w_depth.vcf";
        DataSet orig = ImportUtils.readDataSet(input);

        FilterSiteBuilderPlugin first100 = new FilterSiteBuilderPlugin(null, false)
                .startSite(0)
                .endSite(99);
        GenotypeTable filtered = (GenotypeTable) first100.performFunction(orig).getData(0).getData();

        String expectedResultsFilename = GeneralConstants.EXPECTED_RESULTS_DIR + "mdp_genotype_w_depth_first_100_sites.vcf";
        GenotypeTable expected = ImportUtils.read(expectedResultsFilename);

        AlignmentTestingUtils.alignmentsEqual(expected, filtered);

    }

    @Test
    public void testFilterIndelsVCF() {

        System.out.println("Testing Filtering Indels from VCF File");

        String input = TutorialConstants.TUTORIAL_DIR + "mdp_genotype_w_depth_and_indels.vcf";
        DataSet orig = ImportUtils.readDataSet(input);

        FilterSiteBuilderPlugin filter = new FilterSiteBuilderPlugin(null, false)
                .removeSitesWithIndels(true);
        GenotypeTable filtered = filter.runPlugin(orig);

        String expectedResultsFilename = GeneralConstants.EXPECTED_RESULTS_DIR + "mdp_genotype_w_depth_no_indels.vcf";
        GenotypeTable expected = ImportUtils.read(expectedResultsFilename);

        AlignmentTestingUtils.alignmentsEqual(expected, filtered);

    }
}
