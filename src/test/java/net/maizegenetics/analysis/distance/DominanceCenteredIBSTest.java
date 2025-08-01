/*
 *  DominanceCenteredIBSTest
 * 
 *  Created on Nov 25, 2015
 */
package net.maizegenetics.analysis.distance;

import org.junit.Test;

import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.analysis.data.ExportPlugin;
import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.constants.GeneralConstants;

/**
 *
 * @author Terry Casstevens
 */
public class DominanceCenteredIBSTest {

    @Test
    public void testDominanceCenteredIBSTest() {

        System.out.println("Testing Dominance Centered IBS Matrix ...");

        DataSet input = ImportUtils.readDataSet(TutorialConstants.HAPMAP_FILENAME);
        KinshipPlugin plugin = new KinshipPlugin(null, false);
        DataSet output = plugin
                .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Dominance_Centered_IBS)
                .performFunction(input);

        String tempFile = GeneralConstants.TEMP_DIR + "Dominance_Centered_IBS_Matrix.txt";
        ExportPlugin export = new ExportPlugin(null, false);
        export.setSaveFile(tempFile);
        export.setAlignmentFileType(FileLoadPlugin.TasselFileType.SqrMatrix);
        export.performFunction(output);

        FileLoadPlugin load = new FileLoadPlugin(null, false);
        load.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        load.setOpenFiles(new String[]{GeneralConstants.EXPECTED_RESULTS_DIR + "Dominance_Centered_IBS_Matrix.txt"});
        DataSet expected = load.performFunction(null);

        DistanceMatrixTestingUtils.compare(expected, output, 0.000008);

        load.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        load.setOpenFiles(new String[]{GeneralConstants.TEMP_DIR + "Dominance_Centered_IBS_Matrix.txt"});
        DataSet temp = load.performFunction(null);

        DistanceMatrixTestingUtils.compare(expected, temp, 0.000008);

    }

    @Test
    public void testDominanceCenteredIBSProportionHeterozygousTest() {

        System.out.println("Testing Dominance Centered IBS Matrix w Proportion Heterozygous ...");

        DataSet input = ImportUtils.readDataSet(TutorialConstants.HAPMAP_FILENAME);
        KinshipPlugin plugin = new KinshipPlugin(null, false);
        DataSet output = plugin
                .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Dominance_Centered_IBS)
                .algorithmVariation(KinshipPlugin.ALGORITHM_VARIATION.Proportion_Heterozygous)
                .performFunction(input);

        FileLoadPlugin load = new FileLoadPlugin(null, false);
        load.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        load.setOpenFiles(new String[]{GeneralConstants.EXPECTED_RESULTS_DIR + "Dominance_Centered_IBS_Proportion_Heterozygous_Matrix.txt"});
        DataSet expected = load.performFunction(null);

        DistanceMatrixTestingUtils.compare(expected, output, 0.000008);

    }
}
