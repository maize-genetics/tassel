package net.maizegenetics.analysis.gbs.v2;

import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.util.LoggingUtils;

import org.junit.Test;

public class SNPQualityProfilerPluginTest {

	private final static String statFileName = GBSConstants.GBS_TEMP_DIR + "snpQualityProfilerOutput.txt";
    @Test
    public void testSNPQualityProfilerPlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        System.out.println("Running SNPQualityProfilerPlugin");
//        new SNPQualityProfilerPlugin()
//                .dBFile(GBSConstants.GBS_GBS2DB_FILE)
//                .taxaListName("All")
//                .statFile(statFileName)
//                .performFunction(null);
        new SNPQualityProfilerPlugin()
                .dBFile(GBSConstants.GBS_GBS2DB_FILE)
                .taxaFile(GeneralConstants.DATA_DIR + "CandidateTests/IBMGBSTaxaList.txt")
                .taxaListName("IBM")
                .statFile(statFileName)
                .performFunction(null);
    }

    @Test
    public void testBasicOperationWithPipeline() throws Exception {
      //  (new GBSSeqToTagDBPluginTest()).pipeline();
        testSNPQualityProfilerPlugin();

    }
}