/*
 *  FeatureListToPositionsPlugin
 * 
 *  Created on Nov 3, 2015
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.sql.Connection;
import java.sql.ResultSet;
import javax.swing.ImageIcon;
import static net.maizegenetics.analysis.data.GenomeAnnosDBQueryToPositionListPlugin.connectToDB;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class FeatureListToPositionsPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(FeatureListToPositionsPlugin.class);

    private PluginParameter<String> myConnConfigFile = new PluginParameter.Builder<>("cf", null, String.class)
            .required(true)
            .inFile()
            .guiName("DB Config File")
            .description("DB connection config file")
            .build();

    private PluginParameter<String> myFeatureType = new PluginParameter.Builder<>("featureType", "gene", String.class)
            .description("Feature Type.  Examples: three_prime_utr, gene, exon, five_prime_utr, transcript, cds")
            .build();

    private PluginParameter<String> myFeatureListFilename = new PluginParameter.Builder<>("featureListFilename", null, String.class)
            .description("File with list of features.  If this is not specified, all features of specified type are retrieved.")
            .inFile()
            .build();

    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<>("outputFile", null, String.class)
            .description("Output filename")
            .required(true)
            .outFile()
            .build();

    private PluginParameter<Integer> myNumberBasesPlusOrMinus = new PluginParameter.Builder<>("numberBasesPlusOrMinus", 0, Integer.class)
            .description("Number of bases plus or minus the range of each feature")
            .build();

    public FeatureListToPositionsPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        // String query = "select feature_name, chr_number, feature_range from feature where feature_type = 'gene' AND feature_name in ('AC233950.1_FG002', 'GRMZM5G880025');";
        StringBuilder builder = new StringBuilder();
        builder.append("select feature_name, chr_number, feature_range from feature where feature_type = '");
        builder.append(featureType());
        builder.append("'");

        if ((featureListFilename() != null) && (!featureListFilename().isEmpty())) {
            builder.append(" AND feature_name in (");
            try (BufferedReader reader = Utils.getBufferedReader(featureListFilename())) {
                boolean first = true;
                String line = reader.readLine();
                while (line != null) {
                    line = line.trim();
                    if (!line.isEmpty()) {
                        if (!first) {
                            builder.append(", ");
                        }
                        builder.append("'");
                        builder.append(line);
                        builder.append("'");
                        first = false;
                    }
                    line = reader.readLine();
                }
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("processData: Problem reading gene list file: " + featureListFilename());
            }
            builder.append(")");
        }

        builder.append(" order by chr_number, feature_range;");

        String query = builder.toString();
        myLogger.info("processData: query statement: " + query);

        Connection dbConnection = connectToDB(connConfigFile());
        if (dbConnection == null) {
            throw new IllegalStateException("processData: Problem connecting to database.");
        }

        try (ResultSet rs = dbConnection.createStatement().executeQuery(query)) {

            try (BufferedWriter writer = Utils.getBufferedWriter(Utils.addSuffixIfNeeded(outputFile(), ".bed"))) {

                while (rs.next()) {
                    String featureName = rs.getString("feature_name");
                    String range = rs.getString("feature_range");
                    writer.write(rs.getString("chr_number").trim());
                    writer.write("\t");
                    String[] tokens = range.split(",");
                    String startPosStr = tokens[0].trim().substring(1);
                    String endPosStr = tokens[1].trim().substring(0, tokens[1].length() - 1);
                    // minus one because bed files are 0-based
                    // start position from database and bed files are inclusive
                    int startPos = Math.max(0, Integer.parseInt(startPosStr) - numberBasesPlusOrMinus() - 1);
                    writer.write(String.valueOf(startPos));
                    writer.write("\t");
                    // minus one because bed files are 0-based
                    // end position from database and bed files are exclusive
                    int endPos = Integer.parseInt(endPosStr) + numberBasesPlusOrMinus() - 1;
                    writer.write(String.valueOf(endPos));
                    writer.write("\t");
                    writer.write(featureName);
                    writer.write("\t");
                    writer.write("0");
                    writer.write("\t");
                    writer.write("+");
                    writer.write("\t");
                    writer.write(startPosStr);
                    writer.write("\t");
                    writer.write(endPosStr);
                    writer.write("\n");
                }

            }

        } catch (Exception se) {
            myLogger.debug(se.getMessage(), se);
            throw new IllegalStateException("processData: Problem querying the database: " + se.getMessage());
        }

        return null;

    }
    
    /**
     * DB connection config file
     *
     * @return DB Config File
     */
    public String connConfigFile() {
        return myConnConfigFile.value();
    }

    /**
     * Set DB Config File. DB connection config file
     *
     * @param value DB Config File
     *
     * @return this plugin
     */
    public FeatureListToPositionsPlugin connConfigFile(String value) {
        myConnConfigFile = new PluginParameter<>(myConnConfigFile, value);
        return this;
    }

    /**
     * Feature Type
     *
     * @return Feature Type
     */
    public String featureType() {
        return myFeatureType.value();
    }

    /**
     * Set Feature Type. Feature Type
     *
     * @param value Feature Type
     *
     * @return this plugin
     */
    public FeatureListToPositionsPlugin featureType(String value) {
        myFeatureType = new PluginParameter<>(myFeatureType, value);
        return this;
    }

    /**
     * File with list of features
     *
     * @return Feature List Filename
     */
    public String featureListFilename() {
        return myFeatureListFilename.value();
    }

    /**
     * Set Feature List Filename. File with list of features
     *
     * @param value Feature List Filename
     *
     * @return this plugin
     */
    public FeatureListToPositionsPlugin featureListFilename(String value) {
        myFeatureListFilename = new PluginParameter<>(myFeatureListFilename, value);
        return this;
    }

    /**
     * Output filename
     *
     * @return Output File
     */
    public String outputFile() {
        return myOutputFile.value();
    }

    /**
     * Set Output File. Output filename
     *
     * @param value Output File
     *
     * @return this plugin
     */
    public FeatureListToPositionsPlugin outputFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
        return this;
    }

    /**
     * Number of bases plus or minus the range of each feature
     *
     * @return Number Bases Plus Or Minus
     */
    public Integer numberBasesPlusOrMinus() {
        return myNumberBasesPlusOrMinus.value();
    }

    /**
     * Set Number Bases Plus Or Minus. Number of bases plus or minus the range
     * of each feature
     *
     * @param value Number Bases Plus Or Minus
     *
     * @return this plugin
     */
    public FeatureListToPositionsPlugin numberBasesPlusOrMinus(Integer value) {
        myNumberBasesPlusOrMinus = new PluginParameter<>(myNumberBasesPlusOrMinus, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Feature List to Positions";
    }

    @Override
    public String getToolTipText() {
        return "Feature List to Positions";
    }

}
