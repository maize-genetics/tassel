/*
 *  CenteredIBSTest
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
public class CenteredIBSTest {

    @Test
    public void testCenteredIBSTest() {

        System.out.println("Testing Centered IBS Matrix ...");

        DataSet input = ImportUtils.readDataSet(TutorialConstants.HAPMAP_FILENAME);
        KinshipPlugin plugin = new KinshipPlugin(null, false);
        DataSet output = plugin
                .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS)
                .performFunction(input);

        String tempFile = GeneralConstants.TEMP_DIR + "Centered_IBS_Matrix.txt";
        ExportPlugin export = new ExportPlugin(null, false);
        export.setSaveFile(tempFile);
        export.setAlignmentFileType(FileLoadPlugin.TasselFileType.SqrMatrix);
        export.performFunction(output);

        FileLoadPlugin load = new FileLoadPlugin(null, false);
        load.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        load.setOpenFiles(new String[]{GeneralConstants.EXPECTED_RESULTS_DIR + "Centered_IBS_Matrix.txt"});
        DataSet expected = load.performFunction(null);
        
        DistanceMatrixTestingUtils.compare(expected, output, 0.0000000008);
        
    }
}
