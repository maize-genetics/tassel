package net.maizegenetics.analysis.gbs;


import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.util.CheckSum;
import org.junit.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import static org.junit.Assert.*;

/**

 *  @author Dallas Kroon

 */
public class ProductionPipelineMainTest {

    private static String tempDir = GeneralConstants.TEMP_DIR + "GBS/";
    private static String applicationConfiguration = "production_pipeline.properties";

    private static String exampleAppConfigFile =
            "runFileSuffix=.run\n" +
             "emailHost=appsmtp.mail.cornell.edu\n" +
             "emailAddress=dek29@cornell.edu;dek29@cornell.edu\n" +
             "runDirectory=" + tempDir + "\n" +
             "archiveDirectory=" + GeneralConstants.TEMP_DIR + "\n" +
             "haplosDirectory=/SSD/haplos/\n";

    private static String exampleRunFile =
            "inputFolder=" + GBSConstants.GBS_INPUT_DIR + "\n" +
            "enzyme=ApeKI\n" +
            "topmFile="+ GBSConstants.GBS_INPUT_TOPM + "\n" +
            "outputFolder=" + tempDir + "\n" +
            "keyFile=" + GBSConstants.GBS_DATA_DIR + "Pipeline_Testing_key.txt";


    @Test
    public void testProductionPipelineMain(){

            // create or use a temp dir
            File aTempDir = new File(tempDir);
            if(aTempDir.exists()) aTempDir.delete();
            if(!aTempDir.exists())  aTempDir.mkdir();


            // create a .run file in the temp dir
            // create an application configuration file pointing at the .run file
            BufferedWriter writer = null;
            BufferedWriter writer2 = null;
            try{
                writer = new BufferedWriter(new FileWriter(tempDir + "/temp.run"));
                writer.write(exampleRunFile);
                writer.flush();
                writer.close();

                writer2 = new BufferedWriter(new FileWriter("production_pipeline.properties"));
                writer2.write(exampleAppConfigFile);
                writer2.flush();
                writer2.close();
            }catch (IOException ioe) { ioe.printStackTrace(); }


        String expectedCheckSum = "10e75e612ade7979f210958933da4de9";

        boolean runCheckSum = false, runImputation = false;
        new ProductionPipelineMain(null, runCheckSum, runImputation, null);

        // get checksum for hmp.txt file
        String observedCheckSum = CheckSum.getChecksum(tempDir + "Pipeline_Testing.hmp.txt.gz", "MD5");

        assertEquals(expectedCheckSum, observedCheckSum);

    }
}
