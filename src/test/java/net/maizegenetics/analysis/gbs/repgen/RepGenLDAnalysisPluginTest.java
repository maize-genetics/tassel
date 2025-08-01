/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import static org.junit.Assert.*;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.junit.Test;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;

import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Tuple;

/**
 * @author lcj34
 *
 */
public class RepGenLDAnalysisPluginTest {
    public static final String nextGenDir="/Users/lcj34/notes_files/repGen/";
   // public static final String repGenInputDir="/Volumes/Samsung_T1/repgen/anp68_Locus9_fastq/";
    public static final String repGenInputDir="/Volumes/Samsung_T1/repgen/single_fastq/";
    public static final String repGen_TESTING_KEY_FILE=nextGenDir +"repGen_keyFiles/anp68_locus9_key.txt";

    public static final String repGenDB_FILE=nextGenDir + "repGen.db";
    @Test
    public void testRepGenLDAnalysisPlugin() {
        LoggingUtils.setupDebugLogging();
        System.out.println(Paths.get(repGenDB_FILE).toAbsolutePath().toString());
        try{
            Files.deleteIfExists(Paths.get(repGenDB_FILE));
        } catch (IOException e) {
            e.printStackTrace();
        }
        long time=System.nanoTime();
        System.out.println();
        
        new RepGenLoadSeqToDBPlugin()
        .inputDirectory(repGenInputDir)
        .outputDatabaseFile(repGenDB_FILE)
        .keyFile(repGen_TESTING_KEY_FILE)
        .kmerLength(150)
        .minKmerCount(1) // should be higher when running for real, this is for running just 1 file
        .performFunction(null);
        System.out.printf("Loading RepGen DB:  TotalTime %g sec%n", (double) (System.nanoTime() - time) / 1e9);
        
        System.out.println("\nBegin RepGenLDAnalysisPlugin");
        // DB is now loaded - run the correlation tests
        new RepGenLDAnalysisPlugin()
              .inputDB(repGenDB_FILE)
              .minTaxa(20) // 20 is also the default
              .performFunction(null);
        
        // Get correlations back out 
        System.out.println("\nTags have been correlated - grab a few values");
        RepGenDataWriter repGenData=new RepGenSQLite(repGenDB_FILE);
        Set<Tag> tagsToAlign = repGenData.getTags();
        List<Tag> tagList = new ArrayList<Tag>(tagsToAlign);
        
        List<Tag> shortTagList = new ArrayList<Tag>();

        int numTags = 0;
        System.out.println("LCJ - first 10 tags from db:");
        while(numTags < 10) { // want to test with a small number
            shortTagList.add(tagList.get(numTags));
            System.out.println(tagList.get(numTags).sequence());
            numTags++;
        }
        System.out.println("LCJ - created shortTagLIst, size is : " + shortTagList.size());
        time = System.nanoTime();

     // Now - verify we can get them back out!
        Multimap<Tag,TagCorrelationInfo> tagCorrelationsShortList = repGenData.getCorrelationsForTags(shortTagList);
        
        // Need to verify data within as well.
        // size of tags comes from tagList.  
        int expectedCorrelations = shortTagList.size() * tagList.size() - shortTagList.size();
        assertEquals(tagCorrelationsShortList.keySet().size(),shortTagList.size());
        System.out.println("Expected number of correlations: " + expectedCorrelations + ", actual returned: " + tagCorrelationsShortList.size());
        assertEquals(tagCorrelationsShortList.size(), expectedCorrelations);
        System.out.println("Finished test testRepGenLDAnalysisPlugin, took.");
    }
    
}
