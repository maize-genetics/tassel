package net.maizegenetics.analysis.gbs.v2;

import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.util.LoggingUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import static org.junit.Assert.*;

public class ProductionSNPCallerPluginV2Test {

    @Test
    public void testProcessData() throws Exception {
    	//  This test assumes the database has already been populated.
        // The test will fail if db does not already exist.  Run the 
        // DiscoveryPipeline test to create the database 
    	// (EvaluateSNPCallQualityOfPipelineTest.java)
        LoggingUtils.setupDebugLogging();
        System.out.println(Paths.get(GBSConstants.GBS_GBS2DB_FILE).toAbsolutePath().toString());
        
        // Delete output file if it exists. 
        try{
        	Files.deleteIfExists(Paths.get("junk1.vcf"));
        	Files.deleteIfExists(Paths.get("junk1.h5"));
        } catch (IOException e) {
        	e.printStackTrace();
        }
        
        long time=System.nanoTime();
        System.out.println();

        new ProductionSNPCallerPluginV2()
                .enzyme("ApeKI")
                .inputDirectory(GBSConstants.GBS_INPUT_DIR)
                .inputGBSDatabase(GBSConstants.GBS_GBS2DB_FILE)
                .keyFile(GBSConstants.GBS_TESTING_KEY_FILE)
                .outputGenotypesFile("junk1.vcf")
                .kmerLength(64)
                .minimumQualityScore(0) // tag read quality score
                //.positionQualityScore(5.0) // snp position quality score
                //.maximumMapMemoryInMb(5500)
                .performFunction(null);
        System.out.printf("TotalTime for vcf: %g sec%n", (double) (System.nanoTime() - time) / 1e9);
        
        time=System.nanoTime();
        System.out.println();

        new ProductionSNPCallerPluginV2()
                .enzyme("ApeKI")
                .inputDirectory(GBSConstants.GBS_INPUT_DIR)
                .inputGBSDatabase(GBSConstants.GBS_GBS2DB_FILE)
                .keyFile(GBSConstants.GBS_TESTING_KEY_FILE)
                .outputGenotypesFile("junk1.h5")
                .kmerLength(64)
                .minimumQualityScore(0) // tag read quality score
                //.positionQualityScore(5.0) // snp position quality score
                //.maximumMapMemoryInMb(5500)
                .performFunction(null);
        System.out.printf("TotalTime for h5: %g sec%n", (double) (System.nanoTime() - time) / 1e9);

    }
}