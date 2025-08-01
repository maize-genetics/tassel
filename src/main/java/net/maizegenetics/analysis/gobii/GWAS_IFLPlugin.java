/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import java.awt.Frame;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.ImageIcon;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;

/**
 * The tables populated from this plugin are described in TAS-1162
 * 
 * The plugin takes files of gwas data and adds them to the gwas_data
 * in a GOBII instance.  It is assume the gwas_method and gwas_experiment
 * tables associated with this data have already been populated.  These 
 * are all proprietary tables currently only in use by Buckler Lab.
 * 
 * TO speed up procesing, the experimentId and methodIds are hard-coded
 * in the juint that calls this plugin.  Those tables are generally small,
 * and if you have to go to it to get the name, you might as well just
 * input the ID anda save GOBII IFL processing time.
 * 
 * The values stored in the "values" column will be stored as "real" in the gwas_data table 
 * This is because there is 1 "value" field, which holds values for all data, of any type.  
 * The method table will provide specifics on how to interpret each statistic.
 * 
 * In addition to the .gz files of gwas data, a mapping file of 
 * phenotype names to IDs is created - data pulled from b4R table.
 * @author lcj34
 *
 */
public class GWAS_IFLPlugin extends AbstractPlugin 
{
    private PluginParameter<String> b4rConfigFile= new PluginParameter.Builder<>("b4rConfigFile",null,String.class).guiName("B4R Config File").required(true)
            .description("DB config file containing connection information to the B4r database").build();
    private PluginParameter<String> inputFile= new PluginParameter.Builder<>("inputFile",null,String.class).guiName("Input File").required(true)
            .description("Tab-delimited Input File containing a header line and entries, or Directory containing Tab-delimited files of gwas data.  \nIf parameter is a directory, each file must contain a header line, and the files must end with .txt or .txt.gz").build();
    private PluginParameter<String> methodIds= new PluginParameter.Builder<>("methodIds",null,String.class).guiName("Method IDs").required(true)
            .description("Method Id from the method table.  Must be in same order as the statNames, one for each statname.").build();
    private PluginParameter<String> expID= new PluginParameter.Builder<>("expId",null,String.class).guiName("GWAS Experiment  ID").required(true)
            .description("ID of the GWAS experiment to which this data belongs. This experiment must already existin the gwas_experiment table.").build();
    private PluginParameter<String> statNames = new PluginParameter.Builder<>("statNames",null,String.class).guiName("Statistic Names").required(true)
            .description("Comma separated list of column names from which to pull data.  \nThese names must match values in the statistics field of the gwas_method table for the specified method name.").build();
    private PluginParameter<String> outputDir= new PluginParameter.Builder<>("outputDir",null,String.class).guiName("Path of output directory").required(true)
            .description("Full path name of directory to which output files will be written.  If no user prefix, then end with /").build();
//   private PluginParameter<String> methodIds= new PluginParameter.Builder<>("methodIds",null,String.class).guiName("GWAS Method").required(true)
//           .description("Method id's for each used to create the statistics, e.g. Fast Association, etc").build();
//    private PluginParameter<String> traitMapFile = new PluginParameter.Builder<>("traitMapFile",null,String.class).guiName("Phenotype Trait Mapping File").required(true)
//            .description("tab-delimited File containing a Trait and an ID column to be used to map gwas traits with the B4R ID.").build();

    public GWAS_IFLPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    public GWAS_IFLPlugin() {
        super(null, false);
    }
    private static HashMap<Integer,BufferedWriter> traitWriters= null;
    private static HashMap<Integer,String> traitHM= null;
    private static HashMap<String,String> methodHM = null; // this is trait (from traitHM)/method_id
    private static int chrCol= -1;
    private static int markerCol= -1;
    private static int posCol= -1;
    private static int traitCol = -1;
    @Override
    public DataSet processData(DataSet input) {
        long totalTime = System.nanoTime();
        
        File dataFile = new File(inputFile());
        if (!dataFile.exists()) {
            System.out.println("ERROR - input file doesn't exit: " +  inputFile());
            return null;
        }
        // Create list of files
        List<Path> directoryFiles = new ArrayList<>();
        String inputFileGlob="glob:*{txt,txt.gz,}";
        if (dataFile.isDirectory()) {
            System.out.println("LCJ - input file is a directory");
            directoryFiles= DirectoryCrawler.listPaths(inputFileGlob, Paths.get(inputFile.value()).toAbsolutePath());
            Collections.sort(directoryFiles);
        } else {
            Path inputPath= Paths.get(inputFile()).toAbsolutePath();
            directoryFiles.add(inputPath);
        }
        System.out.println("LCJ - postProcessParamers: size of DirectoryFiles is " + directoryFiles.size());
        for (int idx = 0;idx<directoryFiles.size();idx++) {
            System.out.println("File " + idx +": " + directoryFiles.get(idx));
        }
        
        // get the trait/id mapping from the b4r instance
        Connection dbConnection = GOBIIDbUtils.connectToDB(b4rConfigFile());
        if (dbConnection == null) {
            throw new IllegalStateException("GWAS_IFLPlugin: Problem connecting to database.");
        }
        
        // Find the headers, create the writers
        BufferedReader filebr = Utils.getBufferedReader(directoryFiles.get(0).toString(), 1 << 22);       
        // first line should contain headers.  Find the columns for each data to process
        String headerline;
        try {
            headerline = filebr.readLine();
            filebr.close();
            boolean success = findColumns(headerline,statNames(),methodIds(), outputDir());           
            if (!success) {
                System.out.println("FindColumns failed - quitting");
                return null;
            }
        } catch (IOException ioe) {
            // TODO Auto-generated catch block
            ioe.printStackTrace();
            return null;
        }
               
        // get the phenotype ids;       
        String query = "select id,name from master.variable;";

        Map<String,Integer> traitIDMap = new HashMap<String,Integer>();
        getTraitIds( dbConnection, query, traitIDMap);
        System.out.println("LCJ - Number of entries from b4r trait map: " + traitIDMap.size());
        
        // Create header line for output file
        StringBuilder sb = new StringBuilder();
        // process the data
        try {
            
            // OK, this is more complicated.  For each line in Terry's file, we
            // need to create 3 db lines.  
            // We need to have multiple writers - one for each data column.  in each
            // file, we write all the same data 
            // (phenotype_id/chr/position/method_name/experiment_name/statistic_name/value)
            // but to each one we add the column type, e.g. p-value, or r2.  This goes
            // under "statistic_name"
            // Each one of these files is written to the gwas_data table, but IFL uses the
            // method_name and the value_name column to find the method_id from the gwas_method
            // table.
            // IFL will find expeirment_id from experiment_name column
            String headerLine = "phenotype_id\tmarker\tchr\tposition\texperiment_id\tmethod_id\tstatistic_name\tvalue\n";
            writeHeaderLineToFiles(headerLine);

            sb.setLength(0); // reset the buffer
            for (int idx = 0; idx < directoryFiles.size(); idx++) {               
                int totalLines = 0;
                long time=System.nanoTime();
                Path infile = directoryFiles.get(idx);
                String infileString = infile.toString();
                System.out.println("GWAS_IFLPlugin: processing file " + infileString);
                filebr = Utils.getBufferedReader(infileString, 1 << 22);
                // first line should contain headers.  Find the columns for each data to process
                // and header line must be consistent with all files in the directory
                String line = filebr.readLine(); // toss the header 
                
                // read the data, write to the files.  Each of these files will have the same
                // headers and all will be written into the IFL map table
                sb.setLength(0);
                int count = 0;
                while ((line=filebr.readLine()) != null) {
                    sb.setLength(0);
                    String[] tokens = line.split("\t");
                    // FInd phenotype id
                    int phenotypeId = traitIDMap.get(tokens[traitCol]);
                    sb.append(Integer.toString(phenotypeId));
                    sb.append("\t");
                    sb.append(tokens[markerCol]);
                    sb.append("\t");
                    sb.append(tokens[chrCol]);
                    sb.append("\t");
                    sb.append(tokens[posCol]);
                    sb.append("\t");
                    sb.append(expID());
                    sb.append("\t");
                    writeValues(tokens,sb.toString());
                }
                System.out.println("Process took " + (System.nanoTime() - time)/1e9 + " seconds for file " + infileString);

            }
            System.out.println("TOtalTime for all files: " + (System.nanoTime() - totalTime)/1e9);
            filebr.close();
            Shutdown();
        } catch (IOException ioe) {
            System.out.println("LCJ - caught exception readding or writing data file "  );
            ioe.printStackTrace();
        }
        return null;       
    }
    
    // Grab the trait and id's from the B4R phenotype table.  These will be used
    // to map traits from the input file to phenotype Ids to store in the gwas_data
    // table, field phenotype_id.
    private static void getTraitIds(Connection conn, String query, Map<String,Integer> traitIDMap) {       
        try {
            ResultSet rs = conn.createStatement().executeQuery(query);
            while (rs.next()) {
                int id = rs.getInt("id");
                String trait = rs.getString("name");
                traitIDMap.put(trait, id);
            }
            
        } catch (SQLException sqle) {
            System.out.println("getTraitIds barfed on query: " + query);
            sqle.printStackTrace();
            return;
        }
       return;
    }
    
    private static boolean findColumns( String headerLine, String statNames, String methodIds, String outputDir) {
        System.out.println("LCJ - header line: " + headerLine);
        String [] headers = headerLine.split("\t");
        int index = 0;
        for (String header : headers) {
            if (header.trim().toUpperCase().equals("CHR")) {
                chrCol = index;
            } else if (header.trim().toUpperCase().equals("POS")) {
                posCol = index;
            } else if (header.trim().toUpperCase().equals("TRAIT")) {
                traitCol = index;
            } else if (header.trim().toUpperCase().equals("MARKER")) {
                markerCol = index;
            }
            index++;
        }       
        if (chrCol == -1 || posCol == -1 || traitCol == -1 || markerCol == -1) {
            System.out.println("LCJ - didn't find chr or pos or trait column - quitting");
            return false;
        }
        //Get indices of column names for trait statistics
        int ncols = 0;
        if (statNames!=null && !statNames.equalsIgnoreCase("null") ) {
            traitHM= new HashMap<>(); 
            methodHM = new HashMap<>();
            traitWriters= new HashMap<>(); 
            int search= -1;
            
            String[] traitTokens = statNames.split(",");
            String[] methodTokens = methodIds.split(",");
            for (String str:traitTokens) {
                for (int col = 0; col < headers.length; col++) {
                    if (str.equalsIgnoreCase(headers[col])) {
                        search= col; 
                        break;
                    }
                }
                if (search<0) {
                    System.out.println("Cannot find column "+str);
                    return false;
                } else {
                    String strToStore = str;
                    if (str.equals("p")) strToStore = "pvalue"; // Terry had jus "p", which I don't like
                    traitHM.put(search, strToStore);
                    
                    int methodIdx = Arrays.asList(traitTokens).indexOf(str); // use "str" to find "p", not pvalue
                    methodHM.put(strToStore, methodTokens[methodIdx]);
                    ncols++;
                    System.out.println("Found "+str +" in column "+search);
                }                           
            }
        }
        // create the file writers
        if (traitHM!=null) {
            for (Integer idx:traitHM.keySet()) {
                String outFile = outputDir + traitHM.get(idx) + ".gwas_data"; // IFL wants table name as suffix
                BufferedWriter statWriter = Utils.getBufferedWriter(outFile);
                traitWriters.put(idx,statWriter);
            }
        }
        if (ncols<1) throw new IllegalStateException("No valid columns to read in!"); 
        return true;
    }
    
    // Each file gets same header line
    private void writeHeaderLineToFiles(String headerline) {
        try {
            if (traitHM!=null) {
                for (Integer idx:traitWriters.keySet()) {
                    // write the full line                   
                     traitWriters.get(idx).write(headerline);
                }
            }
        } catch (IOException ioe) {
            System.out.println("LCJ - error writing trait files");
            ioe.printStackTrace();
        }        
    }
    
    private void writeValues(String[] next, String initialLine) {
        // next contains the dataline. We pull only the trait we want from the column
        // that contains that trait.  traitHM contains a list of column.  The name
        // of the trait goes into the "statistic_name" column.  The value goes in last
        // header order:  
        //  String headerLine = "phenotype_id\tchr\tposition\texperiment_id\tmethod_id\tstatistic_name\tvalue\n";
        
        // NOTE:  statistic name is ignored, I kept in so I'd know the value
        try {
            if (traitHM!=null) {
                for (Integer idx:traitWriters.keySet()) {
                    // write the full line
                    String methodId = methodHM.get(traitHM.get(idx)) + "\t";
                    String lineToWrite = initialLine + methodId + traitHM.get(idx) + "\t" + next[idx.intValue()] + "\n";
                     traitWriters.get(idx).write(lineToWrite);
                }
            }
        } catch (IOException ioe) {
            System.out.println("LCJ - error writing trait files");
            ioe.printStackTrace();
        }
    }
       
    private void Shutdown() {
        try {           
            if (traitHM!=null) {
                for (Integer i:traitWriters.keySet()) {
                    traitWriters.get(i).close();
                }
            }

        } catch (Exception exc) {
            System.out.println("Problem with shutdown");
            exc.printStackTrace();
        }
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
    public static void main(String[] args) {
        GeneratePluginCode.generate(GWAS_IFLPlugin.class);
    }
    
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(GWAS_IFLPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
//    public <Type> runPlugin(DataSet input) {
//        return (<Type>) performFunction(input).getData(0).getData();
//    }

    /**
     * Tab-delimited Input File containing a header line and
     * entries, or Directory containing Tab-delimited files
     * of gwas data.  
     * If parameter is a directory, each file must contain
     * a header line, and the files must end with .txt or
     * .txt.gz
     *
     * @return Input File
     */
    public String inputFile() {
        return inputFile.value();
    }

    /**
     * Set Input File. Tab-delimited Input File containing
     * a header line and entries, or Directory containing
     * Tab-delimited files of gwas data.  
     * If parameter is a directory, each file must contain
     * a header line, and the files must end with .txt or
     * .txt.gz
     *
     * @param value Input File
     *
     * @return this plugin
     */
    public GWAS_IFLPlugin inputFile(String value) {
        inputFile = new PluginParameter<>(inputFile, value);
        return this;
    }

    /**
     * Name of GWAS method as it appears in the name field
     * of the gwas_method table
     *
     * @return Method Name
     */
    public String methodIds() {
        return methodIds.value();
    }

    /**
     * Set Method Name. Name of GWAS method as it appears
     * in the name field of the gwas_method table
     *
     * @param value Method Name
     *
     * @return this plugin
     */
    public GWAS_IFLPlugin methodIds(String value) {
        methodIds = new PluginParameter<>(methodIds, value);
        return this;
    }

    /**
     * Name of the GWAS experiment to which this data belongs.
     * This name must already exist in the name field of the
     * gwas_experiment table.
     *
     * @return GWAS Experiment  Name
     */
    public String expID() {
        return expID.value();
    }

    /**
     * Set GWAS Experiment  Name. Name of the GWAS experiment
     * to which this data belongs. This name must already
     * exist in the name field of the gwas_experiment table.
     *
     * @param value GWAS Experiment  Name
     *
     * @return this plugin
     */
    public GWAS_IFLPlugin expID(String value) {
        expID = new PluginParameter<>(expID, value);
        return this;
    }

    /**
     * Comma separated list of column names from which to
     * pull data.  
     * These names must match values in the statistics field
     * of the gwas_method table for the specified method name.
     *
     * @return Statistic Names
     */
    public String statNames() {
        return statNames.value();
    }

    /**
     * Set Statistic Names. Comma separated list of column
     * names from which to pull data.  
     * These names must match values in the statistics field
     * of the gwas_method table for the specified method name.
     *
     * @param value Statistic Names
     *
     * @return this plugin
     */
    public GWAS_IFLPlugin statNames(String value) {
        statNames = new PluginParameter<>(statNames, value);
        return this;
    }

    /**
     * Full path name of directory to which output files will
     * be written, must end with a /
     *
     * @return Path of output directory
     */
    public String outputDir() {
        return outputDir.value();
    }

    /**
     * Set Path of output directory. Full path name of directory
     * to which output files will be written, must end with
     * a /
     *
     * @param value Path of output directory
     *
     * @return this plugin
     */
    public GWAS_IFLPlugin outputDir(String value) {
        outputDir = new PluginParameter<>(outputDir, value);
        return this;
    }

//    /**
//     * Method used to create the statistics, e.g. Fast Association,
//     * etc
//     *
//     * @return GWAS Method
//     */
//    public String method() {
//        return method.value();
//    }
//
//    /**
//     * Set GWAS Method. Method used to create the statistics,
//     * e.g. Fast Association, etc
//     *
//     * @param value GWAS Method
//     *
//     * @return this plugin
//     */
//    public GWAS_IFLPlugin method(String value) {
//        method = new PluginParameter<>(method, value);
//        return this;
//    }
    
    /**
     * DB config file containing connection information to
     * the B4r database
     *
     * @return B4R Config File
     */
    public String b4rConfigFile() {
        return b4rConfigFile.value();
    }

    /**
     * Set B4R Config File. DB config file containing connection
     * information to the B4r database
     *
     * @param value B4R Config File
     *
     * @return this plugin
     */
    public GWAS_IFLPlugin b4rConfigFile(String value) {
        b4rConfigFile = new PluginParameter<>(b4rConfigFile, value);
        return this;
    }

//    /**
//     * tab-delimited File containing a Trait and an ID column
//     * to be used to map gwas traits with the B4R ID.
//     *
//     * @return Phenotype Trait Mapping File
//     */
//    public String traitMapFile() {
//        return traitMapFile.value();
//    }
//
//    /**
//     * Set Phenotype Trait Mapping File. tab-delimited File
//     * containing a Trait and an ID column to be used to map
//     * gwas traits with the B4R ID.
//     *
//     * @param value Phenotype Trait Mapping File
//     *
//     * @return this plugin
//     */
//    public GWAS_IFLPlugin traitMapFile(String value) {
//        traitMapFile = new PluginParameter<>(traitMapFile, value);
//        return this;
//    }


}
