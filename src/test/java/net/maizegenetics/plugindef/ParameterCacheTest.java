package net.maizegenetics.plugindef;

import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.distance.DistanceMatrixTestingUtils;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.pipeline.TasselPipeline;
import net.maizegenetics.util.LoggingUtils;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @author Terry Casstevens
 * Created September 13, 2018
 */
public class ParameterCacheTest {

    @Before
    public void setUp() {
        LoggingUtils.setupDebugLogging();

        String configFile = "dataFiles/ParameterCache/config.txt";
        ParameterCache.load(configFile);
    }

    @Test
    public void testParameterCache() {

        assertEquals("", "VCF", ParameterCache.value("format").get());
        assertEquals("", "Dominance_Centered_IBS", ParameterCache.value("KinshipPlugin.method").get());

    }

    @Test
    public void testKinshipWithParameterCache() {

        String tempFile = GeneralConstants.TEMP_DIR + "DominanceCenteredIBSKinshipActual.txt";
        String[] args = new String[]{"-importGuess", TutorialConstants.HAPMAP_FILENAME, "-KinshipPlugin", "-endPlugin", "-export", tempFile};
        new TasselPipeline(args, null);

        FileLoadPlugin load = new FileLoadPlugin(null, false);
        load.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        load.setOpenFiles(new String[]{"dataFiles/ParameterCache/DominanceCenteredIBSKinshipExpected.txt"});
        DataSet expected = load.performFunction(null);

        load.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        load.setOpenFiles(new String[]{tempFile});
        DataSet actual = load.performFunction(null);

        DistanceMatrixTestingUtils.compare(expected, actual, 0.0001);

    }

}