/**
 * 
 */
package net.maizegenetics.analysis.gbs.v2;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Set;
import java.util.concurrent.atomic.LongAdder;

import javax.swing.ImageIcon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.dna.map.GenomeSequenceBuilder;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TagData;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

/**
 * This plugin queries a GBSv2 sql database for existing tags.
 * 
 * If the user specifies a tag, and the tag is found, a tab-delimited file is written 
 * containing that tag sequence.  If the tag is NOT found, a tab-delimited file is 
 * written that shows the tag and "not found".
 * 
 * If no tag sequence is specified, the method will print out all tag sequences stored 
 * in the specified db.
 * 
 * If input db not found, plugin throws an exception.
 * If output directory doesn't exist, plugin throws an exception
 *
 * @author lcj34
 *
 */
public class GetTagSequenceFromDBPlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(GetTagSequenceFromDBPlugin.class);

    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with tags").build();
    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<String>("o", null, String.class).guiName("Output File").required(true).outFile()
            .description("Output txt file that can be imported to Excel").build();
    private PluginParameter<String> myTagSequence = new PluginParameter.Builder<String>("tagSequence", null, String.class).guiName("Tag Sequence")
            .description("Enter specific tag sequence to verify existence in database.  If no sequence is provided, all tags from the DB will be printed").build();

    public GetTagSequenceFromDBPlugin() {
        super(null, false);
    }

    public GetTagSequenceFromDBPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public GetTagSequenceFromDBPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public void postProcessParameters() {
        if (myDBFile.isEmpty() || !Files.exists(Paths.get(inputDB()))) {
            throw new IllegalArgumentException("GetTagSequenceFromDBPlugin: postProcessParameters: Input DB not set or found");
        }

        try {
            String myOutputDir = (new File(outputFile())).getCanonicalFile().getParent();
        } catch (IOException ioe) {
            throw new IllegalStateException("Problem resolving output directory:" + ioe);
        }
    }
    @Override
    public DataSet processData(DataSet input) {
        try {
            BufferedWriter fileWriter=Utils.getBufferedWriter(outputFile());
            TagData tagData=new TagDataSQLite(inputDB());
            LongAdder count=new LongAdder();
            Set<Tag> dbTags = tagData.getTags();
            try { // Write column header
                String columnHeader = "Tags";
                fileWriter.write(columnHeader); 
            } catch (Exception ioe) {
                myLogger.info("Catch in writing Tag file header error=" + ioe);
                ioe.printStackTrace();
            }
 
            if (tagSequence() != null) {
                Tag requestedTag = TagBuilder.instance(tagSequence()).build();
                if (dbTags.contains(requestedTag)) {
                    String tagString = "\n" + tagSequence();
                    writeTagSequence(fileWriter, tagString);
                } else {
                    String tagString = "\nNOT found: " + tagSequence();
                    writeTagSequence(fileWriter, tagString);
                }
            } else { // no tags specified, print all tag sequences from DB
                tagData.getTags()
                .forEach((tag)->{
                    String tagString = "\n" + tag.sequence();
                    writeTagSequence(fileWriter, tagString);
                    count.increment();
                });
            }

            fileWriter.close();
            ((TagDataSQLite)tagData).close(); 
            myLogger.info("Finished printing tag sequences from database."
                    + "\nTotal number of tags written: " + count.longValue() 
                    + "\nOuput tab-delimited file: " + outputFile() + "\n\n");
        } catch (Exception e) {
            myLogger.info("Catch in reading database error=" + e);
            e.printStackTrace();
        }
        return null;
    }

    private void writeTagSequence(BufferedWriter fileWriter, String tagSequence) {
        try {  
            //String tagString=tagSequence + "\n";
            fileWriter.write(tagSequence);
        } catch (IOException ioe) {
            myLogger.info("Catch in writing TagSequence file error=" + ioe);
            ioe.printStackTrace();
        }
    }

    @Override
    public String getToolTipText() {
        return "Verify single tag sequence exists in DB, or get list of all tag sequences in data base";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Get Tag Sequence from DB";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(GetTagSequenceFromDBPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public String runPlugin(DataSet input) {
        return (String) performFunction(input).getData(0).getData();
    }

    /**
     * Input database file with tags and taxa distribution
     *
     * @return Input DB
     */
    public String inputDB() {
        return myDBFile.value();
    }

    /**
     * Set Input DB. Input database file with tags and taxa
     * distribution
     *
     * @param value Input DB
     *
     * @return this plugin
     */
    public GetTagSequenceFromDBPlugin inputDB(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }

    /**
     * Output tab-delimited file showing tag sequences found in DB.
     *
     * @return Output File
     */
    public String outputFile() {
        return myOutputFile.value();
    }

    /**
     * Set Output File. This is tab delimited file with a list
     * of tag sequences found in the database (either a single tag
     * if the user is checking for a particular tag's presence, or
     * a list of all tags found in the db (if no user tag specified).
     *
     * @param value Output File
     *
     * @return this plugin
     */
    public GetTagSequenceFromDBPlugin outputFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
        return this;
    }

    /**
     * Tag sequence user would like to verify exists in DB.
     *
     * @return tag sequence 
     */
    public String tagSequence() {
        return myTagSequence.value();
    }

    /**
     * Set User tag sequence. This is tag user wants to verify exists in DB.
     * @param value value sequence
     *
     * @return this plugin
     */
    public GetTagSequenceFromDBPlugin tagSequence(String value) {
        myTagSequence = new PluginParameter<>(myTagSequence, value);
        return this;
    }
}
