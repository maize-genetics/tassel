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
public class SplitFile_IFLFilePluginTest {
    @Test
    public void testSplitFile_IFLFilePlugin() throws Exception {
        System.out.println("testSplitFile_IFLFilePlugin: begin testSplitFile_IFLFilePlugin");
        // Call the plugin - verify files can be created
        LoggingUtils.setupDebugLogging();
        
        //String infile = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/gobii_ifl_files/ds4_files/DS_4.marker_linkage_group";
        String infile = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/gobii_ifl_files/gobii_maize3/DS_4.marker";
        String outdir = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/gobii_ifl_files/gobii_maize3/";
        int maximumSize = 10000000;
        
        new SplitFile_IFLFilePlugin()
        .inFile(infile)
        .outputDir(outdir)
        .maxSize(maximumSize)
        .performFunction(null);
        
        System.out.println("\testSplitFile_IFLFilePlugin: finished testSplitFile_IFLFilePlugin");
    }

}
