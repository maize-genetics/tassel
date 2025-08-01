/*
 *  FilterSiteBuilderPluginTest
 * 
 *  Created on Oct 8, 2015
 */
package net.maizegenetics.analysis.filter;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.AlignmentTestingUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.io.JSONUtils;
import net.maizegenetics.plugindef.DataSet;
import static org.junit.Assert.*;
import org.junit.Test;

/**
 *
 * @author Terry Casstevens
 */
public class FilterSiteBuilderPluginTest {

    private static final String TEST_DIR = GeneralConstants.DATA_DIR + "FilterSiteBuilderPlugin/";

    public FilterSiteBuilderPluginTest() {
    }

    /**
     * Test of positionList method, of class FilterSiteBuilderPlugin.
     */
    @Test
    public void testFilterPositionList() {

        System.out.println("Testing Filtering Genotype Table by Position List...");

        String expectedResultsFilename = TEST_DIR + "GenotypeFilteredByPositionList.hmp.txt.gz";
        GenotypeTable expectedAlign = ImportUtils.readFromHapmap(expectedResultsFilename, null);

        DataSet input = ImportUtils.readDataSet(TutorialConstants.HAPMAP_FILENAME);

        PositionList positions = JSONUtils.importPositionListFromJSON(TEST_DIR + "position_list.json.gz");

        DataSet output = new FilterSiteBuilderPlugin(null, false)
                .positionList(positions)
                .performFunction(input);

        GenotypeTable result = (GenotypeTable) output.getDataOfType(GenotypeTable.class).get(0).getData();

        AlignmentTestingUtils.alignmentsEqual(expectedAlign, result);

    }

    @Test
    public void testFilterBedFile() {

        int[] expectedSites = new int[]{
            1947984,
            2914066,
            2914171,
            2915078,
            2915242,
            2973508,
            3706018,
            4175293,
            4429897,
            4429927};

        DataSet input = ImportUtils.readDataSet(TutorialConstants.HAPMAP_FILENAME);

        DataSet output = new FilterSiteBuilderPlugin(null, false)
                .bedFile(TEST_DIR + "test_filter.bed")
                .performFunction(input);

        GenotypeTable result = (GenotypeTable) output.getDataOfType(GenotypeTable.class).get(0).getData();

        assertEquals("Number of sites must be equal: ", expectedSites.length, result.numberOfSites());

        for (int i = 0; i < expectedSites.length; i++) {
            assertEquals("position should be the same: ", expectedSites[i], result.positions().get(i).getPosition());
        }

    }
    
    @Test
    public void testRemoveMinorSNPStates() {

        System.out.println("Testing Removing Minor SNP States...");

        String expectedResultsFilename = TEST_DIR + "GenotypeRemoveMinorSNPStates.hmp.txt.gz";
        GenotypeTable expectedAlign = ImportUtils.readFromHapmap(expectedResultsFilename, null);

        DataSet input = ImportUtils.readDataSet(TutorialConstants.HAPMAP_FILENAME);

        DataSet output = new FilterSiteBuilderPlugin(null, false)
                .removeMinorSNPStates(true)
                .performFunction(input);

        GenotypeTable result = (GenotypeTable) output.getDataOfType(GenotypeTable.class).get(0).getData();

        AlignmentTestingUtils.alignmentsEqual(expectedAlign, result);

    }

}
