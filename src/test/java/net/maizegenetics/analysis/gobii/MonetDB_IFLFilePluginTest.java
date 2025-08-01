/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.FileOutputStream;

import org.junit.Test;

import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Utils;

/**
 * @author lcj34
 *
 */
public class MonetDB_IFLFilePluginTest {
    @Test
    public void testMonetdbIFLFileCreation() throws Exception {
        // NOTE: the outfiledir includes the prefix for the intermeidate files, e.g. DS_1 or DS_3
        
        // FOr gobii_maize3:
        // NOTE:  dataset 3 is Illumina, which is NOT yet populated !! (This used to be 282 Maize AGPv3,
        //       which is no longer populated via Ed instructions (only populating orginal mapset where
        //       SNPs called.
        // Dataset4 is the WGS for the 282 MAIZE, from a vcf file.
        // Dataset5 is the WGS for unimp.
      String outfiledir = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/gobii_ifl_files/gobii_maize3newschema/DS_1";
        //String outfiledir = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/gobii_ifl_files/gobii_hapbreakpoints/DS_3";
        String configFile = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/dbConfigFile_maize3newschema.txt";
        //String configFile = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/dbConfigFile_hapbreakpoints.txt";
     // String datasetName = "ZeaGBSv27impV5_20160209_AGPv2_282";
        //String datasetName = "ZeaWGSraw_hmp321_20160628_AGPv3_282";
        String datasetName = "ZeaGBSv27_20160209_AGPv2_282";
        
        System.out.println("MonetDB_IFLFilePluginTest: begin testMonetdbIFLFileCreation");
        // Call the plugin - verify files can be created
        LoggingUtils.setupDebugLogging();
        
        new MonetDB_IFLFilePlugin()
        .dbConfigFile(configFile)
        .outputDir(outfiledir)
        .datasetName(datasetName)
        .performFunction(null);
        
        System.out.println("\nMonetDB_IFLFilePluginTest: finished testMonetdbIFLFileCreation");
        System.out.println("If your dnarun_id and marker_id files are empty, check that your dataset name is correct");
        
    }
    
    @Test
    public void testMonetdbNumTaxa() throws Exception {
        DataOutputStream writerDNA = null;
        BufferedReader dnasampleBR = null;
        String headerFile = "/Volumes/Samsung_T1/gobii_files/wgs_282_vcf/fastqTest/vcfHeader.txt";
        String headerOut = "/Volumes/Samsung_T1/gobii_files/wgs_282_vcf/fastqTest/vcfHeader_perline.txt";
        try {
            dnasampleBR = Utils.getBufferedReader(headerFile);
            String mline = dnasampleBR.readLine(); 
            String[] tokens = mline.split("\\t");
            writerDNA = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(headerOut)));
            for (int idx=0; idx < tokens.length; idx++) {
                writerDNA.writeBytes(tokens[idx]);
                writerDNA.writeBytes("\n");
            }
            writerDNA.close();
        } catch (Exception exc) {
            exc.printStackTrace();
        }
    }
}
