/**
 * 
 */
package net.maizegenetics.analysis.gbs.v2;

import java.awt.Frame;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.tag.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;

/**
 * This plugin takes a GBSv2 database as input, queries for the tags
 * and their taxa distribution, and creates a tab-delimited file of tag
 * and taxa-distribution.  IT can be used for verifying data in the db.
 *
 * The user may supply a file with tab-delimited columns specifying on which tags to present data.
 * The input "tag" file must have at least 1 column named "TAGS" (any case - Tags/TAGS/tags).  Other
 * data in the file will be ignored.
 *
 * Output:  A tab-delimited file where the first column is the tag, and subsequent columns
 * indicate taxa found in the db.  The values in the "taxa" columns are the number of occurranes of
 * that tag for that taxa.
 * 
 * @author lcj34
 *
 */
public class GetTagTaxaDistFromDBPlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(GetTagSequenceFromDBPlugin.class);

    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with tags").build();
    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<String>("o", null, String.class).guiName("Output File").required(true).outFile()
            .description("Output txt file that can be imported to Excel").build();
    private PluginParameter<String> myTagFile = new PluginParameter.Builder<String>("tagFile", null, String.class).guiName("Tag File").required(false).inFile()
            .description("Input file with list of tags to process.  Any number of tab-delimited columns, but must have 1 called Tags.").build();

    public GetTagTaxaDistFromDBPlugin() {
        super(null, false);
    }

    public GetTagTaxaDistFromDBPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public GetTagTaxaDistFromDBPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

 
    @Override
    public DataSet processData(DataSet input) {
        TagDataWriter tdw = new TagDataSQLite(inputDB());
        try {
            BufferedWriter fileWriter = new BufferedWriter(new FileWriter(outputFile()));
            StringBuilder strB = new StringBuilder();
            TaxaList taxaList = tdw.getTaxaList(); // Used for printing taxon column headers
            strB.append("Tag");
            taxaList.stream().forEach(item -> { // column names are the taxon names
                strB.append("\t");
                strB.append(item.getName());
            });
            strB.append("\n");
            fileWriter.write(strB.toString());
            strB.setLength(0);

            Set<Tag> myTags = new HashSet<Tag>();
            if (tagFile() != null) {
                // Will provide data only for user-supplied tags
                BufferedReader tagFileReader = new BufferedReader(new FileReader(new File(tagFile())));
                String[] headers1 =  tagFileReader.readLine().split("\t");
                List<String> headers = new ArrayList<String>(Arrays.asList(headers1));
                headers.replaceAll(String::toUpperCase); // allow user to have tags/Tags/TAGS etc
                int tagColumn = headers.indexOf("TAGS");
                if (tagColumn == -1) {
                    tagFileReader.close();
                    throw new IllegalArgumentException("Missing Tags column in tag file: " + tagFile());
                }

                String line = null;
                // create list of Tag objects from the user provided file
                while ((line = tagFileReader.readLine()) != null) {
                    String tagSeq = line.split("\t")[tagColumn];
                    Tag tag = TagBuilder.instance(tagSeq).build();
                    myTags.add(tag);
                }
                tagFileReader.close();
            } else {
                // Data will be provided for all tags in the database
                myTags = tdw.getTags();
            }

            // This code expects the user tag sequences to match those found in the db
            int tagcount = 0;
            System.out.println("GetTagTaxaDist: number of tags to process: " + myTags.size());
            for (Tag myTag: myTags) {
                tagcount++;
                // get dist for each taxa
                TaxaDistribution tagTD = tdw.getTaxaDistribution(myTag);
                if (tagTD == null) {
                    System.out.println("GetTagTaxaDist: got null tagTD at tagcount " + tagcount);
                    return null;
                }
                int[] depths = tagTD.depths(); // gives us the depths for each taxon
                strB.append(myTag.sequence());
                for (int idx = 0; idx < depths.length; idx++) {
                    strB.append("\t"); 
                    strB.append(depths[idx]);  // add tag depth                     
                }
                strB.append("\n"); // end of line - start next tag
                fileWriter.write(strB.toString());
                strB.setLength(0);
            }

            fileWriter.close();
            ((TagDataSQLite)tdw).close();  
            myLogger.info("TagsTaxaDistToTabDelim: Finished writing TaxaDistribution \n");
        } catch (Exception exc) {
            myLogger.error("GetTagTaxaDistFromDBPlugin: caught error " + exc);
            exc.printStackTrace();
        }
        return null;       
    }

    public static void main(String[] args) {
        GeneratePluginCode.generate(GetTagTaxaDistFromDBPlugin.class);
    }
    @Override
    public ImageIcon getIcon() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getButtonName() {
        return ("GetTagTaxaDistFromDB");
    }

    @Override
    public String getToolTipText() {
        return ("From a given database, data on taxa distribution of tags is compiled and printed.");
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
    public GetTagTaxaDistFromDBPlugin inputDB(String value) {
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
    public GetTagTaxaDistFromDBPlugin outputFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
        return this;
    }

    /**
     * Input file with list of tags to process.  Any number
     * of tab-delimited columns, but must have 1 called Tags.
     *
     * @return Tag File
     */
    public String tagFile() {
        return myTagFile.value();
    }

    /**
     * Set Tag File. Input file with list of tags to process.
     *  Any number of tab-delimited columns, but must have
     * 1 called Tags.
     *
     * @param value Tag File
     *
     * @return this plugin
     */
    public GetTagTaxaDistFromDBPlugin tagFile(String value) {
        myTagFile = new PluginParameter<>(myTagFile, value);
        return this;
    }


}
