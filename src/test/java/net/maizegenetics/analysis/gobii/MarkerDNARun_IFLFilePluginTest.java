/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.junit.Ignore;
import org.junit.Test;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Utils;

/**
 * @author lcj34
 *
 */
public class MarkerDNARun_IFLFilePluginTest {
    private static final String GOBII_DIR = GeneralConstants.DATA_DIR + "GOBII/";

    @Ignore("#41: needs a live GOBII Postgres database (dbConfigFile_maize2.txt), the 282 WGS germplasm mapping file, "
            + "the 282 WGS VCF and the full AGPv3 reference; none are in tassel_test_data. "
            + "GOBIIDbUtilsTest covers the failed-connection path")
    @Test
    public void testMarkerDNARunFileCreation() throws Exception {
      String configFile = GOBII_DIR + "dbConfigFile_maize2.txt";

      // LCJ - change this outfileDir when you change databases.
      String outfiledir = GeneralConstants.TEMP_DIR + "gobii_ifl_files/";
      
      String mappingFile = GOBII_DIR + "cinta_282Maize_files/Sample_germplasm_information_282_v3_forWGS.txt";
      String inputFile = GOBII_DIR + "wgs_282_vcf/fastq";
      
      String datasetName = "ZeaWGSraw_hmp321_20160628_AGPv3_282"; // pick correct one from db !!
      //String expName = "282Maize_gbsv2";
      String expName = "282Maize_wgs321";
     // String platformName = "GBS27.v2";
      String platformName = "Hapmap_321";
      String refName = "Maize_AGPv3"; // change reference as needed
      String projectId = "1";

      String agpv2Ref = GOBII_DIR + "reference_files/maize_agpv2.fa";
      String agpv3Ref = GOBII_DIR + "reference_files/maize_agpv3.fa";
        
        System.out.println("MarkerDNARunFromHMP_IFLFilePluginTest: begin testMarkerDNARunFileCreation");
        // Call the plugin - verify files can be created
        LoggingUtils.setupDebugLogging();
        
        new MarkerDNARun_IFLFilePlugin()
        .dbConfigFile(configFile)
        .inputFile(inputFile)
        .outputFileDir(outfiledir)
        .refFile(agpv3Ref)
        .mappingFile(mappingFile)
        .mapsetName("AGPV3") // vs AGPV2
        .expName(expName)
        .platformName(platformName)
        .refName(refName)
        .datasetName(datasetName)
        //.projectId(projectId) // no longer and input param - get from mapping file
        .performFunction(null);
        
        System.out.println("\nMarkerDNARunFromHMP_IFLFilePluginTest: finished testMarkerDNARunFileCreation");
        
    }

    @Ignore("#41: needs a live GOBII Postgres database (dbConfigFile_maize2.txt), its marker/DNA-run mapping file "
            + "and the Maize282 ZeaGBSv27 HapMap; none are in tassel_test_data. Prints taxa only and asserts nothing")
    @Test
    public void checkTaxaList() throws Exception {
        // Compare names with taxa list
        String dbConfigFile = GOBII_DIR + "dbConfigFile_maize2.txt";
        String mappingFile = GOBII_DIR + "gobii_cinta_filesToLoad/markerDNARUn_mapping_282raw_v2.txt";
        String inputFileComments = GOBII_DIR + "gobii_curator_training/Maize282_ZeaGBSv27raw_MAF02_AGPv2.hmp.txt";
        
        MarkerDNARun_IFLFilePlugin caller = new MarkerDNARun_IFLFilePlugin();
        
        // Get db connection.  needed to query for table ids based on name
        Connection dbConnection = GOBIIDbUtils.connectToDB(dbConfigFile);
        if (dbConnection == null) {
            throw new IllegalStateException("MOnetdb_IFLFilePLugin:processData: Problem connecting to database.");
        }
        HashMap<String,MarkerDNARun_IFLFilePlugin.HmpTaxaData> taxaDataMap = caller.createTaxaMap(dbConnection,mappingFile);
        //caller.mappingFile(mappingFile);
        
        boolean foundHeader = false;
 
        List<String> taxaListFound = new ArrayList<String>();
        List<String> taxaListNotFound = new ArrayList<String>();
        try {
            BufferedReader markerbr = Utils.getBufferedReader(inputFileComments, 1 << 22);
            String mline;
            int totalLines= 0;
            while ( (mline=markerbr.readLine() )!= null && !foundHeader) {
                totalLines++;
                
                if (mline.startsWith("##")) continue; // toss comments, assumes all comments are at top of file, followed by header
                if (!foundHeader) {
 
                    String[] hdrTokens = mline.split("\t");
                    for (int idx=11; idx < hdrTokens.length; idx++) {
                        MarkerDNARun_IFLFilePlugin.HmpTaxaData taxaData = taxaDataMap.get(hdrTokens[idx]);
                        if (taxaData != null) {
                            taxaListFound.add(hdrTokens[idx]);
                        } else {
                            taxaListNotFound.add(hdrTokens[idx]);
                        }
                    }
                    foundHeader = true;
                    System.out.println("LCJ - found header at line " + totalLines);
                    break;
                }
            }
        } catch (IOException ioe) {
            System.out.println("LCJ - error reading files for checkTaxaList !!");
            ioe.printStackTrace();
        }
        System.out.println("LCJ - taxas found ... ");
        for (int idx = 0;idx < taxaListFound.size(); idx++) {
            System.out.println("    " + taxaListFound.get(idx));
        }
        System.out.println("\nLCJ - taxas NOT found ... ");
        for (int idx = 0;idx < taxaListNotFound.size(); idx++) {
            System.out.println("    " + taxaListNotFound.get(idx));
        }
    }
    
    @Ignore("#41: queries a live GOBII Postgres database (dbConfigFile_maize2.txt) for platform 4 marker names; "
            + "a database export utility rather than a test, and it asserts nothing")
    @Test
    public void testMarkerNames() throws Exception {
      String configFile = GOBII_DIR + "dbConfigFile_maize2.txt";
      String outfile = GeneralConstants.TEMP_DIR + "hmpv3_markers.txt";
      
      // Grab marker names for v3 prior to adding vcf 321 for 282
      // Get db connection.  needed to query for table ids based on name
      Connection dbConnection = GOBIIDbUtils.connectToDB(configFile);
      if (dbConnection == null) {
          throw new IllegalStateException("MOnetdb_IFLFilePLugin:processData: Problem connecting to database.");
      }
      
      // get the dataset id;
      StringBuilder sb = new StringBuilder();
      sb.append("select name from marker where platform_id = 4;");
 
      String query = sb.toString();
      DataOutputStream writerMarker = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile)));
      try {
          ResultSet rs = dbConnection.createStatement().executeQuery(query);
          int count = 0;
          while (rs.next()) {
              String name = rs.getString("name");
              writerMarker.writeBytes(name);
              writerMarker.writeBytes("\n");
              count++;
          }
          writerMarker.close();
          System.out.println("LCJ - finished testMarkerNames, count " + count);
      } catch (SQLException sqle) {
          System.out.println("getTableId barfed on query: " + query);
          sqle.printStackTrace();
      }
      
    }
}
