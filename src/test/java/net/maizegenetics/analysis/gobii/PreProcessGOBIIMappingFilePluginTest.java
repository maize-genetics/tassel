/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import org.junit.Ignore;
import org.junit.Test;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.util.LoggingUtils;

/**
 * Calls PreprocessGOBIIMappingFilePlugin to create intermediate
 * dnasample and germplasm files for loading.
 * @author lcj34
 *
 */
public class PreProcessGOBIIMappingFilePluginTest {
    @Ignore("#41: needs a live GOBII Postgres database (dbConfigFile_gobiilynntestoldschema.txt) and the HapMap v3 "
            + "germplasm mapping file; neither is in tassel_test_data. GOBIIDbUtilsTest covers the failed-connection path")
    @Test
    public void testPreprocessGobiiMappingFile() throws Exception {
        // NOTE: the outfiledir should end with /
        
        System.out.println("LCJ - begin testPreprocessGobiiMappingFile");
      String outfiledir = GeneralConstants.TEMP_DIR + "gobii_ifl_files/gobii_lynntestoldschema/";
      // The dbConfigFile_gobiiTest.txt maps to a db that doesn't have germplasm administered, so all
      // germplasm should show up un the intermediate file
      String configFile = GeneralConstants.DATA_DIR + "GOBII/dbConfigFile_gobiilynntestoldschema.txt";
     // String datasetName = "TESTING_DS_282Panel_WGS";
      String datasetName = "ZeaWGS_hmp321_imputed_AGPv3";
      String mapFile = GeneralConstants.DATA_DIR + "GOBII/mapping_files/germplasmInformation_hapmapv3_20160802_cintaFixesSep1.txt";
        
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
