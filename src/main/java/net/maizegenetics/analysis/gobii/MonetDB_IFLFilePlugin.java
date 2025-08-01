/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;

import javax.swing.ImageIcon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

/**
 * This class has methods to create GOBII intermediary files to be used when 
 * creating a GOBII monetdb dataset table.  The GOBII IFL scripts require
 * 3 files:  a matrix of variants, a list of maker_ids,
 * and a list of dnarun ids.  The matrix of variants is
 * created from the hmp.txt file at the same time the intermediary
 * files for the postgres marker and dnarun related tables are
 * created. (see MarkerDNARunFromHMP_IFLFilePlugin)
 * 
 * Once the IFL scripts have been run to populate these postgres
 * tables, the marker_id and dnarun_id values are created.
 * 
 * This script pulls the id values from the postgres DB to create the final
 * 2 intermediary files needed by the GOBII loadVariantMatrix.py script
 * 
 * See this link for details:
 *  http://cbsugobii05.tc.cornell.edu:6084/display/TD/MonetDB+IFL
 *  
 * The db config file should look like this:
 *   host=cbsudc01.tc.cornell.edu
 *   user=<user id>
 *   password=<your secret password>
 *   DB=gobii_maizeifltest  (or other postgres db you want to query)
 *   
 * The outputDir field should also contain the prefix for the file.
 * Previously this was the dataset.name.  But GOBII names their .h5
 * and monetdb table with DS_<datatset_id>.  We want to do the same to
 * be consistent.  THis will need to be queried from the db before
 * running this plugin.  If the dataset name is known, this is a simple
 * query.
 * 
 * AUgust2016 UPDATE:  The list of markers in the marker_id, and dnarun_ids in
 * dnarun_id file must be in proper order.  They must be in the order the markers
 * are stored in the monetdb table, ie in the order they are stored in the
 * variant file.  TO achieve this, the DB query orders the output by
 * marker_idx  (marker query) and dnarun_idx (dnarun query).  These idx values
 * were created sequentially when the marker and dnarun tables were created.
 *   
 * @author lcj34
 *  
 */
public class MonetDB_IFLFilePlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(MonetDB_IFLFilePlugin.class);

    private PluginParameter<String> dbConfigFile= new PluginParameter.Builder<>("dbConfigFile",null,String.class).guiName("dbConfigFile").required(true)
            .description("DB connection config file").build();
    private PluginParameter<String> datasetName= new PluginParameter.Builder<>("datasetName",null,String.class).guiName("dataset name").required(true)
            .description("Name of dataset whose marker and dnarun IDs are to be pulled").build();
    private PluginParameter<String> outputDir= new PluginParameter.Builder<>("outputDir",null,String.class).guiName("Path of output directory").required(true)
            .description("Full path name of output directory, must end with a /").build();

    public MonetDB_IFLFilePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    public MonetDB_IFLFilePlugin() {
        super(null, false);
    }

    @Override
    public ImageIcon getIcon() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getButtonName() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getToolTipText() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public DataSet processData(DataSet input) {
        // THis method will create the marker_id and dnarun_id files for
        // loading into monetdb.
        String dnarunFile = outputDir() + ".dnarun_id"; // the outputDir should include the DS_<dataset_id>
        String markerFile = outputDir() + ".marker_id";

        //  process the input data file
        try {
            BufferedWriter writerRunID = Utils.getBufferedWriter(dnarunFile);
            BufferedWriter writerMarkerID = Utils.getBufferedWriter(markerFile);

            long time=System.nanoTime();
            // Connect to db
            Connection dbConnection = GOBIIDbUtils.connectToDB(dbConfigFile());
            if (dbConnection == null) {
                throw new IllegalStateException("MonetDB_IFLFilePlugin:processData: Problem connecting to database.");
            }
            
            // create dnarun_id query:
            // select dnarun_id from dataset_dnarun, dataset where dataset.name = datasetName and dataset.dataset_id = dataset_dnarun.dataset_id;
            StringBuilder builder = new StringBuilder();
            builder.append("select dnarun_id from dataset_dnarun, dataset where dataset.name = '");
            builder.append(datasetName());
            builder.append("'");
            builder.append(" and dataset.dataset_id = dataset_dnarun.dataset_id order by dnarun_idx;"); // lcj - added "Order by" Aug 26,2016
            
            String query = builder.toString();
            myLogger.info("processData: query statement: " + query);
            
            System.out.println("MonetDB_IFLFilePlugin: execute query: " + query);
            dbConnection.setAutoCommit(false); // required for Cursor processing (fetchSize)
            Statement st = dbConnection.createStatement();
            st.setFetchSize(100000); // should return results in batches
            ResultSet rs = st.executeQuery(query);
           // ResultSet rs = dbConnection.createStatement().executeQuery(query);   
            // GOBII IFL scripts need a header in the dnarun file (but not in
            // the .variant or .marker files !!
            writerRunID.write("Header\n");
            while (rs.next()) {
                int dnarun_id = rs.getInt("dnarun_id");
                writerRunID.write(String.valueOf(dnarun_id));
                writerRunID.write("\n");
            }
            System.out.printf("TotalTime for dnarun_id query %g sec%n", (double) (System.nanoTime() - time) / 1e9);
            
            // create marker query:
            // select marker_id from dataset_marker, dataset where dataset.name = datasetName and dataset.dataset_id = dataset_marker.dataset_id
            builder = new StringBuilder();
            builder.append("select marker_id from dataset_marker, dataset where dataset.name = '");
            builder.append(datasetName());
            builder.append("'");
            builder.append(" and dataset.dataset_id = dataset_marker.dataset_id order by marker_idx;"); // lcj - added order by 8/26/16
//            builder.append("select marker_id from dataset_marker, dataset where dataset.name = '");
//            builder.append(datasetName());
//            builder.append("'");
//            builder.append(" and dataset.dataset_id = dataset_marker.dataset_id order by marker_id;"); // lcj - added order by 8/26/16
            
            query = builder.toString();
            myLogger.info("processData: query statement: " + query);
            
            System.out.println("MonetDB_IFLFilePlugin: execute query: " + query);
            st.setFetchSize(100000); // shouldn't need to set this again
            rs = st.executeQuery(query);
           // rs = dbConnection.createStatement().executeQuery(query);           
            while (rs.next()) {
                int marker_id = rs.getInt("marker_id");
                writerMarkerID.write(String.valueOf(marker_id));
                writerMarkerID.write("\n");
            }
            
            writerRunID.close();
            writerMarkerID.close();
            System.out.printf("TotalTime for marker_id query: %g sec%n", (double) (System.nanoTime() - time) / 1e9);

        } catch (Exception exc) {
            System.out.println("Monetdb_IFLFile:  caught exception processing writing files");
            exc.printStackTrace();
        }
        System.out.println("\nFiles written to " + dnarunFile + " and " + markerFile);
        return null;
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.

    //     public static void main(String[] args) {
    //         GeneratePluginCode.generate(MonetDB_IFLFilePlugin.class);
    //     }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    //     public <Type> runPlugin(DataSet input) {
    //         return (<Type>) performFunction(input).getData(0).getData();
    //     }

    /**
     * DB connection config file
     *
     * @return dbConfigFile
     */
    public String dbConfigFile() {
        return dbConfigFile.value();
    }

    /**
     * Set dbConfigFile. DB connection config file
     *
     * @param value dbConfigFile
     *
     * @return this plugin
     */
    public MonetDB_IFLFilePlugin dbConfigFile(String value) {
        dbConfigFile = new PluginParameter<>(dbConfigFile, value);
        return this;
    }

    /**
     * Name of dataset whose marker and dnarun IDs are to
     * be pulled
     *
     * @return dataset name
     */
    public String datasetName() {
        return datasetName.value();
    }

    /**
     * Set dataset name. Name of dataset whose marker and
     * dnarun IDs are to be pulled
     *
     * @param value dataset name
     *
     * @return this plugin
     */
    public MonetDB_IFLFilePlugin datasetName(String value) {
        datasetName = new PluginParameter<>(datasetName, value);
        return this;
    }

    /**
     * Full path name of output directory, must end with a
     * /
     *
     * @return Path of output directory
     */
    public String outputDir() {
        return outputDir.value();
    }

    /**
     * Set Path of output directory. Full path name of output
     * directory, must end with a /
     *
     * @param value Path of output directory
     *
     * @return this plugin
     */
    public MonetDB_IFLFilePlugin outputDir(String value) {
        outputDir = new PluginParameter<>(outputDir, value);
        return this;
    }
}
