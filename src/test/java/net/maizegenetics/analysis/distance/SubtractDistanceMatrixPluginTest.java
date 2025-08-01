/*
 *  SubtractDistanceMatrixPluginTest
 * 
 *  Created on Dec 18, 2015
 */
package net.maizegenetics.analysis.distance;

import java.util.ArrayList;
import java.util.List;
import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.plugindef.DataSet;
import org.junit.Test;

/**
 *
 * @author Terry Casstevens
 */
public class SubtractDistanceMatrixPluginTest {

    @Test
    public void testSubtractCenteredIBS() {

        System.out.println("Testing Subtract Centered IBS Matrix ...");

        DataSet input = ImportUtils.readDataSet(TutorialConstants.HAPMAP_FILENAME);

        FilterSiteBuilderPlugin filter = new FilterSiteBuilderPlugin(null, false);
        input = filter
                .startSite(100)
                .endSite(500)
                .performFunction(input);

        KinshipPlugin kinship = new KinshipPlugin(null, false);
        DataSet output = kinship
                .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS)
                .performFunction(input);

        SubtractDistanceMatrixPlugin subtract = new SubtractDistanceMatrixPlugin(null, false);
        DataSet rest = subtract
                .wholeMatrix(GeneralConstants.EXPECTED_RESULTS_DIR + "Centered_IBS_Matrix.txt")
                .performFunction(output);

        FileLoadPlugin load = new FileLoadPlugin(null, false);
        load.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        load.setOpenFiles(new String[]{GeneralConstants.EXPECTED_RESULTS_DIR + "Centered_IBS_Matrix_Minus_100-500.txt"});
        DataSet expected = load.performFunction(null);

        DistanceMatrixTestingUtils.compare(expected, rest, 0.000001);

    }

    @Test
    public void testSubtractMultipleCenteredIBS() {

        System.out.println("Testing Subtract Multiple Centered IBS Matrix ...");

        DataSet input = ImportUtils.readDataSet(TutorialConstants.HAPMAP_FILENAME);

        List<DataSet> subMatrices = new ArrayList<>();

        FilterSiteBuilderPlugin filter = new FilterSiteBuilderPlugin(null, false);
        DataSet input0 = filter
                .startSite(100)
                .endSite(300)
                .performFunction(input);

        KinshipPlugin kinship = new KinshipPlugin(null, false);
        DataSet output = kinship
                .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS)
                .performFunction(input0);
        subMatrices.add(output);

        FilterSiteBuilderPlugin filter1 = new FilterSiteBuilderPlugin(null, false);
        DataSet input1 = filter1
                .startSite(301)
                .endSite(400)
                .performFunction(input);

        KinshipPlugin kinship1 = new KinshipPlugin(null, false);
        DataSet output1 = kinship1
                .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS)
                .performFunction(input1);
        subMatrices.add(output1);
        
        FilterSiteBuilderPlugin filter2 = new FilterSiteBuilderPlugin(null, false);
        DataSet input2 = filter2
                .startSite(401)
                .endSite(500)
                .performFunction(input);

        KinshipPlugin kinship2 = new KinshipPlugin(null, false);
        DataSet output2 = kinship2
                .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS)
                .performFunction(input2);
        subMatrices.add(output2);

        SubtractDistanceMatrixPlugin subtract = new SubtractDistanceMatrixPlugin(null, false);
        DataSet rest = subtract
                .wholeMatrix(GeneralConstants.EXPECTED_RESULTS_DIR + "Centered_IBS_Matrix.txt")
                .performFunction(DataSet.getDataSet(subMatrices, null));

        FileLoadPlugin load = new FileLoadPlugin(null, false);
        load.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        load.setOpenFiles(new String[]{GeneralConstants.EXPECTED_RESULTS_DIR + "Centered_IBS_Matrix_Minus_100-500.txt"});
        DataSet expected = load.performFunction(null);

        DistanceMatrixTestingUtils.compare(expected, rest, 0.000008);

    }

}
