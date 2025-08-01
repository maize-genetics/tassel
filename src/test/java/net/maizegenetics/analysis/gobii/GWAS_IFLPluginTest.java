/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import org.junit.Test;

import net.maizegenetics.util.LoggingUtils;

/**
 * @author lcj34
 *
 */
public class GWAS_IFLPluginTest {
    @Test
    public void testGWAS_IFLFilePlugin() throws Exception {

        String outfiledir = "/Volumes/Samsung_T1/gobii_files/ifl_binaries/gobii_gwas/DNDiff_"; // adding prefix to output files
        String configFile = "/Users/lcj34/notes_files/gobiiANDBms/gwas/b4rconfig.txt";
        String inputFileDir = "/Volumes/Samsung_T1/gobii_files/gwas/terry_DNDiff/"; 
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
