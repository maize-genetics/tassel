package net.maizegenetics.analysis.rna;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Set;

import net.maizegenetics.analysis.gbs.v2.GBSUtils;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TagDataWriter;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.util.LoggingUtils;

import org.junit.Test;

public class RNADeMultiplexSeqToDBTest {

    private static String myInputDir = GeneralConstants.DATA_DIR + "RNA";
    private static String myKeyFile = GeneralConstants.DATA_DIR + "RNA/Key_file_smallRNA_test_plugin.txt";
    private static String myOutputDB = GeneralConstants.TEMP_DIR + "RNASeq.db";
        
    @Test
    public void testProcessData() throws Exception {
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
        System.out.println("testProcessData:  Total number of tags in db: " + dbTags.size());
       // assertEquals(8,dbTags.size());
        assertEquals(1102,dbTags.size());
        
        ((TagDataSQLite)tdw).close();        
    }
    
    @Test
    public void testTissueSet() throws Exception {
        // This test is to verify creating the tissue set works.
        LoggingUtils.setupDebugLogging();
        
        String keyfile2 = GeneralConstants.DATA_DIR + "RNA/Leaf_Root_key_ProDBtest.txt";
        String keyfile3 = GeneralConstants.DATA_DIR + "RNA/Leaf_key_noTissueValues.txt";
        
        List<String> tissueSet = TaxaListIOUtils.readTissueAnnotationFile(myKeyFile, GBSUtils.tissueNameField);
        
        System.out.println("Tissues in key file:");
        assertEquals(2,tissueSet.size());
        if (tissueSet != null && tissueSet.size() > 0) {
            tissueSet.forEach(System.out::println);
        } else {
            System.out.println("TissueSet ONE is NULL");
        }
        
        tissueSet = TaxaListIOUtils.readTissueAnnotationFile(keyfile2, GBSUtils.tissueNameField);       
        System.out.println("Test with keyfile that contains tissue column, size of file: "
                + tissueSet.size());
        assertEquals(2,tissueSet.size());
        if (tissueSet != null && tissueSet.size() > 0) {
            tissueSet.forEach(System.out::println);
        } else {
            System.out.println("TissueSet TWO is NULL");
        }
        
        tissueSet = TaxaListIOUtils.readTissueAnnotationFile(keyfile3, GBSUtils.tissueNameField);
        assertEquals(0,tissueSet.size());
        System.out.println("Test with keyfile that contains tissue column, NO DATA, size of file: "
                + tissueSet.size()); 
        if (tissueSet != null && tissueSet.size() > 0) {
            tissueSet.forEach(System.out::println);
        } else {
            System.out.println("TissueSet THREE is NULL");
        }
    }
}
