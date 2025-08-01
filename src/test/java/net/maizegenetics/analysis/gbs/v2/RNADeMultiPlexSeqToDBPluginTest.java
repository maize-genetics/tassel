package net.maizegenetics.analysis.gbs.v2;

import static org.junit.Assert.*;

import java.nio.file.Paths;
import java.util.List;
import java.util.Set;

import org.junit.Test;

import net.maizegenetics.analysis.rna.RNADeMultiPlexSeqToDBPlugin;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TagDataWriter;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.util.LoggingUtils;

public class RNADeMultiPlexSeqToDBPluginTest {
    private static String myInputDir = "/Users/lcj34/notes_files/hackathon_Oct2015/fastQ_sy";
    private static String myKeyFile = "/Users/lcj34/notes_files/hackathon_Oct2015/Leaf_key.txt";
    private static String myOutputDB = "/Users/lcj34/notes_files/hackathon_Oct2015/RNASeq.db";
    

    
    @Test
    public void testProcessData() throws Exception {
        //  This test assumes the database has already been populated.
        // The test will fail if db does not already exist.  Run the 
        // DiscoveryPipeline test to create the database 
        // (EvaluateSNPCallQualityOfPipelineTest.java)
        LoggingUtils.setupDebugLogging();
        
        new RNADeMultiPlexSeqToDBPlugin()
        .inputDirectory(myInputDir)
        .keyFile(myKeyFile)
        .outputDatabaseFile(myOutputDB)
        .minimumKmerLength(18)
        .minimumQualityScore(0) // tag read quality score
        .minKmerCount(1)
        .performFunction(null);
        
        // Verify values in DB
        TagDataWriter tdw=new TagDataSQLite(myOutputDB);
        
        Set<Tag> dbTags = tdw.getTags();
        assertEquals(8,dbTags.size());
        
    }
    
    @Test
    public void testTissueSet() throws Exception {
        // This test is to verify creating the tissue set works.
        LoggingUtils.setupDebugLogging();
        
        String keyfile2 = "/Users/lcj34/notes_files/hackathon_Oct2015/ProDBplugin_testFile/Leaf_Root_key_ProDBtest.txt";
        String keyfile3 = "/Users/lcj34/notes_files/hackathon_Oct2015/Leaf_key_noTissueValues.txt";
        
        List<String> tissueSet = TaxaListIOUtils.readTissueAnnotationFile(myKeyFile, GBSUtils.tissueNameField);
        
        if (tissueSet != null && tissueSet.size() > 0) {
            tissueSet.forEach(System.out::println);
        } else {
            System.out.println("LCJ - tissueSet ONE is NULL");
        }
        
        tissueSet = TaxaListIOUtils.readTissueAnnotationFile(keyfile2, GBSUtils.tissueNameField);       
        System.out.println("LCJ - now test with keyfile that contains tissue column, size of file: "
                + tissueSet.size());
        if (tissueSet != null && tissueSet.size() > 0) {
            tissueSet.forEach(System.out::println);
        } else {
            System.out.println("LCJ - tissueSet ONE is NULL");
        }
        
        tissueSet = TaxaListIOUtils.readTissueAnnotationFile(keyfile3, GBSUtils.tissueNameField);
        System.out.println("LCJ - now test with keyfile that contains tissue column, NO DATA, size of file: "
                + tissueSet.size()); 
        if (tissueSet != null && tissueSet.size() > 0) {
            tissueSet.forEach(System.out::println);
        } else {
            System.out.println("LCJ - tissueSet ONE is NULL");
        }
    }
}
