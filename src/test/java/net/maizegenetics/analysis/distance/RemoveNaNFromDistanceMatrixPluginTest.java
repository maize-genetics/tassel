/*
 *  RemoveNaNFromDistanceMatrixPluginTest
 * 
 *  Created on Mar 2, 2016
 */
package net.maizegenetics.analysis.distance;

import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.plugindef.DataSet;
import org.junit.Test;

/**
 *
 * @author Terry Casstevens
 */
public class RemoveNaNFromDistanceMatrixPluginTest {

    @Test
    public void testRemoveNaNCenteredIBS() {

        System.out.println("Testing Remove NaN Centered IBS Matrix ...");

        FileLoadPlugin load = new FileLoadPlugin(null, false);
        load.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        load.setOpenFiles(new String[]{GeneralConstants.DATA_DIR + "RemoveNaNFromDistanceMatrixPlugin/Centered_IBS_Matrix.txt"});
        DataSet original = load.performFunction(null);

        DataSet actual = RemoveNaNFromDistanceMatrixPlugin.runPlugin(original);

        load = new FileLoadPlugin(null, false);
        load.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        load.setOpenFiles(new String[]{GeneralConstants.DATA_DIR + "RemoveNaNFromDistanceMatrixPlugin/Centered_IBS_No_NaN.txt"});
        DataSet expected = load.performFunction(null);

        DistanceMatrixTestingUtils.compare(actual, expected, 0.0);

    }

}
