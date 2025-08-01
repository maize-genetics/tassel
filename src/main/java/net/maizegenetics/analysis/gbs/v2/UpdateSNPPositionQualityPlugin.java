/**
 * 
 */
package net.maizegenetics.analysis.gbs.v2;

import java.awt.Frame;
import java.nio.file.Path;
import java.nio.file.Paths;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListIOUtils;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TagDataWriter;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * This plugin takes as input:
 *    tab-delimited txt file with columns CHROM, POS, QUALITYSCORE 
 *    dbFile:  a GBSv2 database with snppositions recorded
 *  A PositionList of positions with quality scores is sent to the database
 *  where the snpposition table is updated with a qualityScore value for the
 *  specified chromosome and position.
 *  
 * @author lcj34
 *
 */
public class UpdateSNPPositionQualityPlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(UpdateSNPPositionQualityPlugin.class);

    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with SNP positions stored").build();
    private PluginParameter<String> myQSFile = new PluginParameter.Builder<String>("qsFile", null, String.class).guiName("Quality Score File").required(true).inFile()
            .description("tab-delimited txt file containing headers CHROM(String), POS(Integer) and QUALITYSCORE(Float) for filtering SNP positions from database").build();

    public UpdateSNPPositionQualityPlugin() {
        super(null, false);
    }

    public UpdateSNPPositionQualityPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public UpdateSNPPositionQualityPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        try {
            Path qsPath= Paths.get(qsFile()).toAbsolutePath();
            PositionList posListWithQS = PositionListIOUtils.readQualityScoreFile(qsPath.toString());
            if (posListWithQS == null) {
            	String errMsg = "Error: no data read from file  " + qsPath.toString();
            	myLogger.error(errMsg);
            	return null;
            }
            TagDataWriter tdw=new TagDataSQLite(inputDB());
            // Write the qualityScore positions to the snpposition table
            tdw.putSNPPositionQS(posListWithQS);      
            ((TagDataSQLite)tdw).close();  
            myLogger.info("UupdateSNPPositionQualityPlugin: Finished writing quality scores file to snpposition table.\n");
        } catch (Exception exc) {
            myLogger.error("Caught error adding quality scores to the database " + exc);
            exc.printStackTrace();
        }
        return null;
    }

    @Override
    public String getToolTipText() {
        return "Update SNP Position quality score";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Update SNP Position Quality";
    }

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
    public UpdateSNPPositionQualityPlugin inputDB(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }

    /**
     * Input quality score file to use for setting the qualityScore
     * in the database's snpposition table.
     *
     * @return QualityScore File
     */
    public String qsFile() {
        return myQSFile.value();
    }

    /**
     * Set qsFile. Input quality score file used for setting
     * the qualityScore field in the snpposition table.
     *
     * @param value Quality Score File
     *
     * @return this plugin
     */
    public UpdateSNPPositionQualityPlugin qsFile(String value) {
        myQSFile = new PluginParameter<>(myQSFile, value);
        return this;
    }
}
