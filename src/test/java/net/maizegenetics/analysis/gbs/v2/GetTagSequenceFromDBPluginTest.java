/**
 * 
 */
package net.maizegenetics.analysis.gbs.v2;

import java.util.Iterator;
import java.util.Set;
import java.util.concurrent.atomic.LongAdder;

import org.junit.Test;

import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagData;
import net.maizegenetics.dna.tag.TagDataSQLite;

/**
 * This class prints tag sequences from a specified database.
 * Each test case in this class assumes the GBSv2 pipeline has
 * been run and there exists a database file found in  
 *   {user directory path}/tempDir/GBS/Chr9_10-20000000/GBSv2.db
 * 
 * @author lcj34
 *
 */
public class GetTagSequenceFromDBPluginTest {

    private final static String dbSingleTagOutFile=GBSConstants.GBS_TEMP_DIR + "dbSingleTag.txt";
    private final static String dbAllTagsOutFile=GBSConstants.GBS_TEMP_DIR + "dbMultipleTag.txt";

    @Test
    public void testSingleTagInput() throws Exception {
        // Get tags from db, grab and entry, send to the plugin.
        // Plugin code should find this value in the db.
        // Output file should show one entry that is this tag.
        TagData tagData=new TagDataSQLite(GBSConstants.GBS_GBS2DB_FILE);
        Set<Tag> dbTags = tagData.getTags();
        ((TagDataSQLite)tagData).close();
        
        Tag myTag = null;
        for (Iterator<Tag> it = dbTags.iterator(); it.hasNext();) {
             myTag = it.next();
        }
            
        new GetTagSequenceFromDBPlugin()
        .inputDB(GBSConstants.GBS_GBS2DB_FILE)
        .outputFile(dbSingleTagOutFile)
        .tagSequence(myTag.sequence())
        .performFunction(null);
    }
    
    @Test
    public void testSingleTagNOTFound() throws Exception {           
        // Create dummy sequence that should not exist as tag in db
        // Output file should show one entry:  "NOT found: AAAAACCCCGGGGNNNN"
        new GetTagSequenceFromDBPlugin()
        .inputDB(GBSConstants.GBS_GBS2DB_FILE)
        .outputFile(dbSingleTagOutFile)
        .tagSequence("AAAAACCCCGGGGNNNN")
        .performFunction(null);
    }
        
    @Test
    public void testMultipleTagInput() throws Exception {       
        // Grab all tags - no need to access DB from the test -
        // source plugin handles it all.
        // Output file should show a list of tags from the specified DB
        new GetTagSequenceFromDBPlugin()
        .inputDB(GBSConstants.GBS_GBS2DB_FILE)
        .outputFile(dbAllTagsOutFile)
        .performFunction(null);
    }
}
