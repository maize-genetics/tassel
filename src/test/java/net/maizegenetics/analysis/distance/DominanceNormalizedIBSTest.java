/*
 *  DominanceNormalizedIBSTest
 * 
 *  Created on Mar 10, 2016
 */
package net.maizegenetics.analysis.distance;

import org.junit.Ignore;
import org.junit.Test;

import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.constants.GeneralConstants;

/**
 *
 * @author Terry Casstevens
 */
public class DominanceNormalizedIBSTest {
    @Ignore
    @Test
    public void testDominanceNormalizedIBSTest() {

        System.out.println("Testing Dominance Normalized IBS Matrix ...");

        DataSet input = ImportUtils.readDataSet(TutorialConstants.HAPMAP_FILENAME);
        KinshipPlugin plugin = new KinshipPlugin(null, false);
        DataSet output = plugin
                .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Dominance_Normalized_IBS)
                .performFunction(input);

        FileLoadPlugin load = new FileLoadPlugin(null, false);
        load.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        load.setOpenFiles(new String[]{GeneralConstants.EXPECTED_RESULTS_DIR + "Dominance_Normalized_IBS_Matrix.txt"});
        DataSet expected = load.performFunction(null);

        DistanceMatrixTestingUtils.compare(expected, output, 0.0000000008);

    }

}
