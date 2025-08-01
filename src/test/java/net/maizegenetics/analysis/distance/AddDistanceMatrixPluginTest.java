/*
 *  AddDistanceMatrixPluginTest
 * 
 *  Created on Feb 22, 2016
 */
package net.maizegenetics.analysis.distance;

import java.util.ArrayList;
import java.util.List;
import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.plugindef.DataSet;
import org.junit.Test;

/**
 *
 * @author Terry Casstevens
 */
public class AddDistanceMatrixPluginTest {

    @Test
    public void testAddCenteredIBS() {

        System.out.println("Testing Add Centered IBS Matrix ...");

        DataSet input = ImportUtils.readDataSet(TutorialConstants.HAPMAP_FILENAME);
        int numSites = ((GenotypeTable) input.getData(0).getData()).numberOfSites();

        List<DataSet> submatrices = new ArrayList<>();

        FilterSiteBuilderPlugin filter1 = new FilterSiteBuilderPlugin(null, false);
        DataSet input1 = filter1
                .startSite(0)
                .endSite(500)
                .performFunction(input);

        KinshipPlugin kinship = new KinshipPlugin(null, false);
        submatrices.add(kinship
                .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS)
                .performFunction(input1));

        FilterSiteBuilderPlugin filter2 = new FilterSiteBuilderPlugin(null, false);
        DataSet input2 = filter2
                .startSite(501)
                .endSite(numSites - 1)
                .performFunction(input);

        submatrices.add(kinship
                .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS)
                .performFunction(input2));

        AddDistanceMatrixPlugin add = new AddDistanceMatrixPlugin(null, false);
        DataSet total = add.performFunction(DataSet.getDataSet(submatrices, null));

        FileLoadPlugin load = new FileLoadPlugin(null, false);
        load.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        load.setOpenFiles(new String[]{GeneralConstants.EXPECTED_RESULTS_DIR + "Centered_IBS_Matrix.txt"});
        DataSet expected = load.performFunction(null);

        DistanceMatrixTestingUtils.compare(expected, total, 0.0000008);

    }

    @Test
    public void testAddNormalizedIBS() {

        System.out.println("Testing Add Normalized IBS Matrix ...");

        DataSet input = ImportUtils.readDataSet(TutorialConstants.HAPMAP_FILENAME);
        int numSites = ((GenotypeTable) input.getData(0).getData()).numberOfSites();

        List<DataSet> submatrices = new ArrayList<>();

        FilterSiteBuilderPlugin filter1 = new FilterSiteBuilderPlugin(null, false);
        DataSet input1 = filter1
                .startSite(0)
                .endSite(500)
                .performFunction(input);

        KinshipPlugin kinship = new KinshipPlugin(null, false);
        submatrices.add(kinship
                .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Normalized_IBS)
                .performFunction(input1));

        FilterSiteBuilderPlugin filter2 = new FilterSiteBuilderPlugin(null, false);
        DataSet input2 = filter2
                .startSite(501)
                .endSite(numSites - 1)
                .performFunction(input);

        submatrices.add(kinship
                .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Normalized_IBS)
                .performFunction(input2));

        AddDistanceMatrixPlugin add = new AddDistanceMatrixPlugin(null, false);
        DataSet total = add.performFunction(DataSet.getDataSet(submatrices, null));

        FileLoadPlugin load = new FileLoadPlugin(null, false);
        load.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        load.setOpenFiles(new String[]{GeneralConstants.EXPECTED_RESULTS_DIR + "Normalized_IBS_Matrix.txt"});
        DataSet expected = load.performFunction(null);

        DistanceMatrixTestingUtils.compare(expected, total, 0.0000008);

    }

}
