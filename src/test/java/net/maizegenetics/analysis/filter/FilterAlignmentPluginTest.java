/*
 *  FilterAlignmentPluginTest
 * 
 *  Created on Feb 19, 2015
 */
package net.maizegenetics.analysis.filter;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;

import net.maizegenetics.dna.snp.AlignmentTestingUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import org.junit.Test;

/**
 *
 * @author Terry Casstevens
 */
public class FilterAlignmentPluginTest {

    public FilterAlignmentPluginTest() {
    }

    /**
     * Test of setFilterMinorSNPs method, of class FilterAlignmentPlugin.
     */
    @Test
    public void testSetFilterMinorSNPs() {
        System.out.println("FilterAlignmentPluginTest: setFilterMinorSNPs");

        String inputFile = TutorialConstants.TUTORIAL_DIR + "mdp_genotype_w_rare_alleles.hmp.txt";
        GenotypeTable input = ImportUtils.readFromHapmap(inputFile);
        FilterAlignmentPlugin filter = new FilterAlignmentPlugin(null, false);
        filter.setMinFreq(0.0);
        filter.setMaxFreq(1.0);
        filter.setMinCount(0);
        filter.setFilterMinorSNPs(true);
        DataSet result = filter.performFunction(new DataSet(new Datum("input", input, null), null));
        GenotypeTable noRareAlleles = (GenotypeTable) result.getData(0).getData();

        String expectedPed = GeneralConstants.EXPECTED_RESULTS_DIR + "mdp_genotype_w_rare_alleles_removed.plk.ped";
        String expectedMap = GeneralConstants.EXPECTED_RESULTS_DIR + "mdp_genotype_w_rare_alleles_removed.plk.map";
        GenotypeTable expected = ImportUtils.readFromPLink(expectedPed, expectedMap, null);

        AlignmentTestingUtils.alignmentsEqual(expected, noRareAlleles);
    }

}
