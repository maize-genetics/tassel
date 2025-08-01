/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.sql.Connection;
import java.sql.ResultSet;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ImageIcon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

/**
 * This plugin should be run prior to creating the intermediate files for marker and
 * dnarun.
 * 
 * There are 3 purposes to this plugin's.  Using the mapping file created for the dataset:
 *
 *   1.  Identify duplicate/missing germplasm/dnasample entries, create intermediate
 *       file for germplasm and dnasmple tables, load any missing entries.  Duplicates are skipped.
 *   2.  Identify duplicate libraryPrepIds.  Write a list of duplicate libraryPrepIds, write
 *        to a file.  
 *   3.  Provide mapping data to load new marker/dnarun related tables.  Create intermediate
 *       files, load via GOBII IFL scripts
 *   
 * For the first 2 purposes, the database must be queried.  Missing entries
 * entries are defined as below:
 *   germplasm table:  From the db,Get list of distinct MGIDs (they should all be distinct).  use this
 *       list to compare to MGIDs in the file.  For any MGIDs that don't appear, create a
 *       line in the *.germplasm intermediate file used to add values.
 *   dnasample table:  From the db, Get a list of dnasample names.  These names are a string comprised
 *       of these components:  GID:plate:well.  From the input file, for each entry, create
 *       a concatenanted string of GID:plate:well.  compare to list from db.  For any names
 *       that don't appear, create a line in the *.dnasample intermediate file for loading.
 *       This file needs the "name" field to be a concatenation of GID:plate:well as this will be
 *       unique and GOBII dnasample.dupmap looks at only the name field.  Code can be MGID if
 *       we need that stored (which I think we do).  It takes "external code" column instead of
 *       germplasm_id as that maps to the external_code field in the germplasm table when GOBII IFL
 *       looks to find the germplasm_id from DB.
 *       This file also needs project_name, which comes from the mapping file.
 *   dnarun table:  From the db, Get a list of all dnasample.name fields.  These should be
 *       distinct library prep id.  Compare to libraryPrepIds from the mapping file.  IF there
 *       are duplicate, write to a file to show the biologist.
 *       
 * NOTES:  GOBII uses dnasample.name and dnasample.num to determine duplicates
 *    BL is not populating dnasample.num.  "num" has been removed from the dnasample.dupmap
 *    file when running this.  For some reason, with it present, but all values "null", the
 *    script believed the values were different and I ended up duplicating all dnasamples when
 *    sending the file through the GOBII scripts.  When I removed this line, the scripts only
 *    checked the "name" field and project id and it worked.
 *    
 * For step 3:  The intermediate files are created by the MarkerDNARunMGID_fromHMPIFIFIlePLugin.java.
 * Note the dnasample and germplasm entries must be loaded to the db before loading the marker/
 * dnarun intermediate files or the necssary db ids will not be found..  
 * 
 * @author lcj34
 *
 */
public class PreProcessGOBIIMappingFilePlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(PreProcessGOBIIMappingFilePlugin.class);
    
    private PluginParameter<String> dbConfigFile= new PluginParameter.Builder<>("dbConfigFile",null,String.class).guiName("dbConfigFile").required(true)
            .description("DB connection config file").build();
    private PluginParameter<String> datasetName= new PluginParameter.Builder<>("datasetName",null,String.class).guiName("dataset name").required(true)
            .description("Name of existing database dataset.  Will be used to pull dataset_id from the db.  This id is incorporated into the output file names.").build();
    private PluginParameter<String> mappingFile = new PluginParameter.Builder<>("mappingFile",null,String.class).guiName("mappingFile").required(true)
            .description("tab-delimited File containing columns: taxaColumn, name, source,MGID, GID,libraryID, plate_code, well, species, type, project").build();
    private PluginParameter<String> outputDir= new PluginParameter.Builder<>("outputDir",null,String.class).guiName("Path of output directory").required(true)
            .description("Full path name of output directory, must end with a /").build();
    
    public PreProcessGOBIIMappingFilePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    public PreProcessGOBIIMappingFilePlugin() {
        super(null, false);
    }
    
    @Override
    public DataSet processData(DataSet input) {
        // THis method will create germplasm and dnasample IFL intermediate files for any
        // data represented in the mapping file that does not already exist in the germplasm or
        // dnasample table.  We expect most if not all to be in there, but as time goes by, there
        // may be some that need to be added.
        
        // Additionally, this method will create an output file containing a list of library prep Ids
        // that already exist in the dnarun table as "name".  The user is responsible for checking this
        // list and  taking appropriate action.  The libraryPrepIds should be unique and should only
        // occur once in this db.  If we've already processed this, we want to know.  
        
        // QUSTION FOR CINTA:  did we say there may be a reason we run the same sample/libraryPrepId
        // multiple times for different analysis?  Or does that analysis result in a different library
        // prep ID, so this test is still valid ??
        
        //  process the input data file
        try {

            // Connect to db
            Connection dbConnection = GOBIIDbUtils.connectToDB(dbConfigFile());
            if (dbConnection == null) {
                throw new IllegalStateException("PreProcessGOBIIMappingFilePlugin: Problem connecting to database.");
            }
            // get the dataset id;
            StringBuilder sb = new StringBuilder();
            sb.append("select dataset_id from dataset where name = '");
            sb.append(datasetName());
            sb.append("';");
            ResultSet rs = dbConnection.createStatement().executeQuery(sb.toString());
            int datasetId = -1;
            while (rs.next()) {
                datasetId = rs.getInt("dataset_id");                
            }
            if (datasetId < 0) {
                System.out.println("Could not find datasetId from datasetName " + datasetName() + " please check name and try again !!");
                return null;
            }
            System.out.println("LCJ - processing dataset number " + datasetId);
            String germplasmFile = outputDir() + "DS_" + datasetId + ".germplasm"; // the outputDir should include the DS_<dataset_id>
            String dnasampleFile = outputDir() + "DS_" + datasetId + ".dnasample";
            String dupLibIDFile = outputDir() + "DS_" + datasetId + ".dup_libraryPrepIds";
            
            BufferedWriter writergp = Utils.getBufferedWriter(germplasmFile);
            BufferedWriter writerdna = Utils.getBufferedWriter(dnasampleFile);
            BufferedWriter writerlib = Utils.getBufferedWriter(dupLibIDFile);
            // create germplasm query: We want all of them to compare to those in our mapping file
            // Select external_code from germplasm.  WE are storing nothing in name, and GID in
            // the external code field.  Every run with a new GID gets loaded into germplasm
            StringBuilder builder = new StringBuilder();
            builder.append("select external_code from germplasm;");           
            String query = builder.toString();
            myLogger.info("processData: query statement: " + query);
            
            System.out.println("PreProcessGOBIIMappingFilePlugin: execute query: " + query);
            rs = dbConnection.createStatement().executeQuery(query);   
            // GOBII IFL scripts need a header in the dnarun file (but not in
            // the .variant or .marker files !!
            List<String> existingGermplasm = new ArrayList<String>();
            while (rs.next()) {
                existingGermplasm.add(rs.getString("external_code"));
            }
                  
            // NOTE: If there are no current entries in the germplasm list (e.g.
            // prior to populating the first dataset in the db) then all from the
            // mapping file will be added.
            
            // create dnasample query:
            // select name from dnasample: This will eventually be an extraction_id,
            // from the lab where it is processed, or a created one by Cinta (or me)
            // for now (july 11,2016) we have the name as GID:plate:wel
            builder = new StringBuilder();
            builder.append("select name from dnasample;");
            query = builder.toString();
            myLogger.info("processData: query statement: " + query);
            
            System.out.println("PreProcessGOBIIMappingFilePlugin: execute query: " + query);
            rs = dbConnection.createStatement().executeQuery(query); 
            List<String> existingDNASample = new ArrayList<String>();
            while (rs.next()) {
                existingDNASample.add(rs.getString("name")); 
            }
 
            // create dnarun query.
            builder = new StringBuilder();
            builder.append("select name from dnarun;");
            query = builder.toString();
            myLogger.info("processData: query statement: " + query);
            
            System.out.println("PreProcessGOBIIMappingFilePlugin: execute query: " + query);
            rs = dbConnection.createStatement().executeQuery(query); 
            List<String> existingLib = new ArrayList<String>();
            while (rs.next()) {
                existingLib.add(rs.getString("name")); 
            }           
            // AFTER we have the lists, go through the mapping file ONCE, reading
            // all the colums we need for each line, and comparing.  
            // Get the GID value.  If it is NOT on the existingGermplasm list, then
            // create a line for the germplasm intermediate file and add it.
            // Check the GID:plate:well - concatenate, check against the dnasample name,
            // add missing names to dnasample intermediate file
            // Finally, get the libraryPrepID;  check against the dnarun name.  if a
            // duplicate, write to the dup_libID file
            
            // Create string builders for the 3 files - we'll append data, then write
            StringBuilder germplasmSB = new StringBuilder();
            StringBuilder dnasampleSB = new StringBuilder();
            
            // write header lines:
            germplasmSB.append("name\texternal_code\tspecies_name\ttype_name\tcreated_by\tcreated_date\tmodified_by\tmodified_date\tstatus\tcode\n");
            writergp.write(germplasmSB.toString());
            
            // THis doesn't work without all columns.  In order to ignore "num" as a duplicate column,
            // the dnasample.dupmap file was changed:  the "num" entry was removed.
            dnasampleSB.append("name\tcode\tplatename\tnum\twell_row\twell_col\tproject_name\texternal_code\tcreated_by\tcreated_date\tmodified_by\tmodified_date\tstatus\n");
            writerdna.write(dnasampleSB.toString());
            
            // column names
            int taxaIdx=-1, nameIdx=-1, sourceIdx=-1, mgidIdx=-1, gidIdx=-1, libIdx=-1, plateIdx=-1;
            int wellIdx=-1, speciesIdx=-1, typeIdx=-1, projectIdx = -1, sampleNameIdx = -1;
            BufferedReader germplasmBR = Utils.getBufferedReader(mappingFile());
            String mline = germplasmBR.readLine(); // header line
 
            // write header line to duplicate library prep id. 
            writerlib.write(mline); // will the \n still be there?
            
            System.out.println("\nPreprocessGObii: getting header columns from mline: " + mline);
            String [] headers = mline.split("\\t");
            //parse headers
            if (mline.contains("TaxaColumn")) {
                int idx = 0;
                for (String header : headers) {
                    if (header.trim().toUpperCase().equals("TAXACOLUMN")) {
                        taxaIdx = idx;
                    } else if (header.trim().toUpperCase().equals("NAME")) {
                        nameIdx = idx;
                    }else if (header.trim().toUpperCase().equals("SOURCE")) {
                        sourceIdx = idx;
                    }else if (header.trim().toUpperCase().equals("MGID")) {
                        mgidIdx = idx;
                    }else if (header.trim().toUpperCase().equals("GID")) {
                        gidIdx = idx;
                    }else if (header.trim().toUpperCase().equals("LIBRARYID")) {
                        libIdx = idx;
                    }else if (header.trim().toUpperCase().equals("PLATE_CODE")) {
                        plateIdx = idx;
                    }else if (header.trim().toUpperCase().equals("WELL")) {
                        wellIdx = idx;
                    }else if (header.trim().toUpperCase().equals("SPECIES")) {
                        speciesIdx = idx;
                    }else if (header.trim().toUpperCase().equals("TYPE")) {
                        typeIdx = idx;
                    }else if (header.trim().toUpperCase().equals("PROJECT")) {
                        projectIdx = idx;
                    }else if (header.trim().toUpperCase().equals("SAMPLENAME")) {
                        sampleNameIdx = idx;
                    }                       
                    idx++;
                }
            } else {
                System.out.println("Mappingfile is missing header line !!!");
                return null;
            }
 
            // LCJ ***** FIx the use of SampleName whenvery ou get a good xcel sheet.  I HATE EXCeL !!
            if (taxaIdx == -1 || nameIdx == -1 || sourceIdx == -1 || mgidIdx == -1 || gidIdx == -1 || libIdx == -1 ||
                  //  plateIdx == -1 || wellIdx == -1 || speciesIdx == -1 || typeIdx == -1 || projectIdx == -1 ) {
                 plateIdx == -1 || wellIdx == -1 || speciesIdx == -1 || typeIdx == -1 || projectIdx == -1 || sampleNameIdx == -1) {
                System.out.println("\nMappingfile is missing required header line.  Expecting columns: TaxaColumn, name, source, MGID, GID, libraryID, plate_code, well, species, type, project, SampleName\n");
                return null;
            }
            System.out.println("PreprocessGobii: processing mapping file: " + mappingFile());
            // read all lines, check for missing germplasm/dnasample or duplice libprepIDs
            int dnaNotAdded = 0;
            int germplasmNotAdded = 0;
            int totalLines = 0;
            List<String> addedGermplasm = new ArrayList<String>();
            while ((mline = germplasmBR.readLine()) != null) {
                String[] data = mline.split("\\t");
                totalLines++;
                // Check germplasm
                // There shoudl only be entries in this file that have assigned GIDs. Sometimes there
                // is a problem, and Robert didn't get something merged.  However, I don't want these
                // in the file.  Should not need the check for data[gidIdx] != null
                if (data[gidIdx].trim() != null && 
                        !(existingGermplasm.contains(data[gidIdx].trim())) &&
                        !(addedGermplasm.contains(data[gidIdx].trim()))) {
                    // tab over fields not filled in, add values for others.  Name is skipped (first column)
                    // external_code (2nd column) is now GID
                    String gpentry = "\t" + data[gidIdx].trim() + "\t" + data[speciesIdx].trim() + "\t" + data[typeIdx].trim() + "\t\t\t\t\t1\t0\n";
                    writergp.write(gpentry);
                    addedGermplasm.add(data[gidIdx].trim()); // record that we already have this GID - don't want dups!
                } else {
                    System.out.println("LCJ - not adding " + data[nameIdx].trim() + " to germplasm file");
                    germplasmNotAdded++;
                }
                //System.out.println("\n");
                
                // LCJ **** Replace sampleName = sampleNameIdx line when excel file is fixed !!
                String sampleName = null;
                if (sampleNameIdx == -1 || data[sampleNameIdx].trim().equals("")) {                   
                    sampleName = data[gidIdx].trim() + ":" + data[plateIdx].trim() + ":" + data[wellIdx].trim();
                } else {
                    sampleName = data[sampleNameIdx].trim();
                }
                // These were 2 original values of sampleName.  We want the value from sampleNameIdx from
                // the mapping file.  BUt if the sampleName column doesn't exist, or the value is null
                // (Cinta had a file will some missing values) then use the gid with plate and well.  Could
                // be that some of those columns are null as well ....
                //String sampleName = data[sampleNameIdx].trim();    
               // String sampleName = data[gidIdx].trim() + ":" + data[plateIdx].trim() + ":" + data[wellIdx].trim();
                if (data[gidIdx] != null && !(existingDNASample.contains(sampleName))) {
                    // tab over fields not filled in, add values for others, 
                    // "dummycode" is stored for code.  format is:
                    // name,code,platename,num,well_row,well_col,project_name,external_code,created_by,created_date,modified_by,modified_date,status
                    String wellRow = data[wellIdx].trim();
                    String wellCol = "";
                    if (!(wellRow).equals("")){ // this field may be blank
                        wellRow = wellRow.substring(0,1);
                        wellCol = data[wellIdx].substring(1);
                    }
                    String dnaentry = sampleName + "\tdummycode\t" + data[plateIdx].trim() + "\t\t" + wellRow + "\t" + wellCol + "\t" 
                            + data[projectIdx].trim() + "\t" + data[gidIdx].trim() + "\t" + "6\t\t\t\t1\n";
                    writerdna.write(dnaentry);
                } else {
                    System.out.println("LCJ - not adding " + data[nameIdx].trim() + " to dnasample file");
                    dnaNotAdded++;
                }
                if (data[gidIdx].trim() != null && (existingLib.contains(data[libIdx].trim()))) {
                    // Found a duplicate library prep ID - record it
                    writerlib.write(mline);
                }
            }
            writergp.close();
            writerdna.close();
            writerlib.close();
            System.out.println("\nFiles written to " + germplasmFile + " and " + dnasampleFile);
            System.out.println("Total mapping file lines: " + totalLines + " Not added to germplasm:" 
                + germplasmNotAdded + ", not added to dnasample:" + dnaNotAdded + "\n");

        } catch (Exception exc) {
            System.out.println("PreProcessGOBIIMappingFilePlugin:  caught exception processing or writing files");
            exc.printStackTrace();
        }
     
        
        return null;
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
    
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(PreProcessGOBIIMappingFilePlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
//    public <Type> runPlugin(DataSet input) {
//        return (<Type>) performFunction(input).getData(0).getData();
//    }

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
    public PreProcessGOBIIMappingFilePlugin dbConfigFile(String value) {
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
    public PreProcessGOBIIMappingFilePlugin datasetName(String value) {
        datasetName = new PluginParameter<>(datasetName, value);
        return this;
    }

    /**
     * tab-delimited File containing columns: taxaColumn,
     * name, MGID, GID,libraryID, plate_code, well, species,
     * type, project_id, experiment_name, platform_name, reference_name
     * and dataset_name
     *
     * @return mappingFile
     */
    public String mappingFile() {
        return mappingFile.value();
    }

    /**
     * Set mappingFile. tab-delimited File containing columns:
     * taxaColumn, name, MGID, GID,libraryID, plate_code,
     * well, species, type, project_id, experiment_name, platform_name,
     * reference_name and dataset_name
     *
     * @param value mappingFile
     *
     * @return this plugin
     */
    public PreProcessGOBIIMappingFilePlugin mappingFile(String value) {
        mappingFile = new PluginParameter<>(mappingFile, value);
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
    public PreProcessGOBIIMappingFilePlugin outputDir(String value) {
        outputDir = new PluginParameter<>(outputDir, value);
        return this;
    }
}
