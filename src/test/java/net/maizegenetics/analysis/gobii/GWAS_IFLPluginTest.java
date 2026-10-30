/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import org.junit.Ignore;
import org.junit.Test;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.util.LoggingUtils;

/**
 * @author lcj34
 *
 */
public class GWAS_IFLPluginTest {
    @Ignore("#41: needs a live GOBII Postgres database with GWAS experiments and its B4R config file, plus the "
            + "DNDiff GWAS results; none are in tassel_test_data. GOBIIDbUtilsTest covers the failed-connection path")
    @Test
    public void testGWAS_IFLFilePlugin() throws Exception {

        String outfiledir = GeneralConstants.TEMP_DIR + "gobii_gwas/DNDiff_"; // adding prefix to output files
        String configFile = GeneralConstants.DATA_DIR + "GOBII/b4rconfig.txt";
        String inputFileDir = GeneralConstants.DATA_DIR + "GOBII/terry_DNDiff/"; 
        String expId = "3";
        String methodIds = "1,2,3"; // the experiment order needs to match the statname order, ie 1= experiment for p,, 2 for r2, 3=df
        String statnames = "p,r2,df"; // code will change "p" to "pvalue"
        System.out.println("testGWAS_IFLFIlePluginTest: begin GWAS_IFLPluginTest");
        // Call the plugin - verify files can be created
        LoggingUtils.setupDebugLogging();
        
        new GWAS_IFLPlugin()
        .b4rConfigFile(configFile) // ifl scripts do the id mapping for us
        .outputDir(outfiledir)
        .inputFile(inputFileDir)
        .methodIds(methodIds)
        .expID(expId)
        .statNames(statnames)
        .performFunction(null);
        
        System.out.println("\ntestGWAS_IFLFilePluginTest: finished GWAS_IFLPluginTest");
        
    }
}
