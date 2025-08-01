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
public class HapBreakpoints_IFLFilePluginTest {
    @Test
    public void testHapBreakpointsIFLFilePlugin() throws Exception {

        String outfiledir = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/gobii_ifl_files/gobii_hapbreakpoints/";
        String configFile = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/dbConfigFile_hapbreakpoints.txt";
        //String breakFile = "/Volumes/Samsung_T1/gobii_files/ames_hapmap_2and3/kelly_ames_unimputed_pa/Ames_HM32_Unimp_chr10.pa.txt";
        String breakFile = "/Volumes/Samsung_T1/gobii_files/ames_hapmap_2and3/kelly_ames_unimputed_pa/"; // all files in directory
        String setname = "ames_unimp";
        String mapset = "AGPV3";
        String method = "FILLIN";
        String src_dataset = "ZeaGBSv27_20150803_AGPv3_Ames";
        // Donor map file is the WGS dataset file
        String donorMapFile = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/mapping_files/germplasmInformation_hapmapv3_20160802_cintaFixesSep1.txt";
        // taxaMapFile is the map file matching the dataset for this breakpoint block.  Initially it is Ames unimputed
        String taxaMapFile = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/mapping_files/germplasmInformation_Ames_20160810.txt";
        System.out.println("HapBreakpoitns_IFLFIlePluginTest: begin testHapBreapointsIFLFilePlugin");
        // Call the plugin - verify files can be created
        LoggingUtils.setupDebugLogging();
        
        new HapBreakpoints_IFLFilePlugin()
        //.dbConfigFile(configFile)
        .outputDir(outfiledir)
        .breakFile(breakFile)
        .setName(setname)
        .mapset(mapset)
        .src_dataset(src_dataset)
        .method(method)
        .donorMapFile(donorMapFile)
        .taxaMapFile(taxaMapFile)
        .performFunction(null);
        
        System.out.println("\nHapBreakopints_IFLFilePluginTest: finished testHapBreakpointsIFLFilePlugin");
        
    }
}
