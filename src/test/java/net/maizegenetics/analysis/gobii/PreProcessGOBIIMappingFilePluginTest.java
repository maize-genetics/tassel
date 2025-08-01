/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import org.junit.Test;

import net.maizegenetics.util.LoggingUtils;

/**
 * Calls PreprocessGOBIIMappingFilePlugin to create intermediate
 * dnasample and germplasm files for loading.
 * @author lcj34
 *
 */
public class PreProcessGOBIIMappingFilePluginTest {
    @Test
    public void testPreprocessGobiiMappingFile() throws Exception {
        // NOTE: the outfiledir should end with /
        
        System.out.println("LCJ - begin testPreprocessGobiiMappingFile");
      String outfiledir = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/gobii_ifl_files/gobii_lynntestoldschema/";
        //String outfiledir = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/gobii_ifl_files/gobii_hapbreakpoints/";
      // The dbConfigFile_gobiiTest.txt maps to a db that doesn't have germplasm administered, so all
      // germplasm should show up un the intermediate file
      String configFile = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/dbConfigFile_gobiilynntestoldschema.txt";
      //String configFile = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/dbConfigFile_hapbreakpoints.txt";
     // String datasetName = "TESTING_DS_282Panel_WGS";
      String datasetName = "ZeaWGS_hmp321_imputed_AGPv3";
      //String mapFile = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/mapping_files/germplasmInformation_hapmapv3_20160802.txt";
      //String mapFile = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/mapping_files/germplasmInformation_Ames_20160810.txt";
      String mapFile = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/mapping_files/germplasmInformation_hapmapv3_20160802_cintaFixesSep1.txt";
        
        System.out.println("PreProcessGOBIIMappingFilePluginTest: begin processing file: " + mapFile);
        // Call the plugin - verify files can be created
        LoggingUtils.setupDebugLogging();
        
        new PreProcessGOBIIMappingFilePlugin()
        .dbConfigFile(configFile)
        .outputDir(outfiledir)
        .datasetName(datasetName) // can be anything - is prefix to the *.germplasm and *.dnasample file. GOBII doesn't use it
        .mappingFile(mapFile)
        .performFunction(null);
        
        System.out.println("\nPreprocessGobiiMappingFilePluginTest: finished testPreprocessGobiiMappingFile");
        
    }
}
