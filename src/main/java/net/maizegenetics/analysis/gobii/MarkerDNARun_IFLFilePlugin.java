/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import java.awt.Frame;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
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
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GenomeSequence;
import net.maizegenetics.dna.map.GenomeSequenceBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;

/**
 * BEWARE:  Whenever you "pull" to update the GOBII projects, there could be
 * changes that effect these plugins.  Need to check if the .nmap or .dupmap
 * files have changed, and make corresponding changes here as required.
 * 
 * This class takes a hmp.txt file(s) or vcf file(s) with a mapping file and creates the
 * intermediate files for the marker, marker_linkage_group, dataset_marker,
 * dnarun, and dataset_dnarun tables.
 * 
 * The inputFile variable can be a file or a directory.  If it is a directory,
 * the code will look for all files with format *(hmp.txt,hmp.txt.gz,vcf,vcf.gz)
 * and process them.  It is assumed all files use the same taxa.
 * 
 * Because we assume it is all the same taxa, the dnarun and dataset_dnarun files
 * are created from the first *.hmp.txt file processed.  These are the intermediate
 * files that map this dnarun to a dataset, and contain one entry for each taxa
 * which contains the taxa name (in the name field), libraryPrepID (as the code field),
 * and ids into experiment and dnasample tables.
 * 
 * The "mapping file" needs to contain columnns for the following data:
 *   taxaname: as appears in the vcf/hmp file
 *   name: taxa name it maps to (do I need this?)
 *   MGID: MGID for this taxa name
 *   GID: GID for this dnarun
 *   libraryID: same as in dnasample file
 *   project_name: db will be queried to get project_id from project name.
 *      Needed by IFL get get dnasample_id
 *   experiment_name: name of experiment needed for dnarun table (IFL maps to id)
 *   platform_name: name of platform needed for marker table, (IFL maps to ID)
 *   reference_name: name of reference table (IFL maps to ID)
 *   dataset_name: needed for dataset_dnarun and dataset_marker tables (IFL Maps to ID)
 *   samplename: will be used for table dnasample.name field
 *   
 *   The mapping file needs an entry for all taxa that may appear in the data input file.
 *   It is ok if multiple taxa names appear with the same MGID/GID/etc.  These are
 *   synonyms.  We mostly aren't storing the names, just the MGID.  It must be 
 *   identified in the mapping file.
 *   
 * THe dataset id must be gotten from the database.  Check the dataset from the mapping
 * file, query the database to get the dataset_id.  GOBII creates the data_table and
 * data_file names from the GUI when it creates the ID.  It always names them DS_<dataset_id>.h5
 * and DS_<dataset_id> for the table.  We must do this by hand as we want to maintain
 * consistency.
 * 
 * The marker_linkage_group:  Their mapping now requires both marker_name and platform_id.
 * So Platform_name must also be a parameter.  The software will query the db to get the
 * platform_id from platform name. Both the marker and the marker_linkage_group intermediate 
 * files need the platform_name.  This could be moved to the mapping file, but currently
 * is an input parameter.
 * 
 * VCF file headers have these fields:
 *   #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  taxa1      taxa2 ...
 * HMP.txt file headers have these fields:
 *   rs#  alleles chrom   pos  strand  assembly#  center  protLSID  assayLSID panelLSID QCcode taxa1 ...
 *   
 * Class Gobii_IFLUtils is used to find chrom,pos, alt and strand values based on
 * file type of hmp or vcf.  THe type of file is determined by the file suffix (hmp.txt, 
 * hmp.txt.gz, vcf or vcf.gz)
 * 
 * July 6, 2016:  In addition, this method will create a file to be used
 * with PostProcessMarkerPlugin().  This file will contain the marker name, platformid
 * and alts array.  It may be used at a future date to find existing markers in the DB and 
 * update the alts array.  See PostProcessMarkerPlugin() for details.  Currently any allele 
 * in A/C/G/T that is not the reference will appear on the alt list.  This is per Ed who
 * says given a large enough population, each allele will appear as an alternate.
 * 
 * August 3:  BEcause we continue to change the data that makes up the sample name
 * (was GID:plate:well, now is extraction_id) I have added a column called "SampleName"
 * to the mapping file.  The software will take whatever is stored here and use it as
 * the dnasample name.  Biologists can then change at will without a need to change
 * the software
 * 
 * Problems with GOBII loaders:  The GOBII dnarun.nmap file now takes "num" instead of "platename"
 * EIther one can be a problem for BL as they are not required fields, and we often don't have
 * values for these columns.  The IFL scripts preprocessor_ifile.py does not check for IS NULL.
 * It merely checks if the input file and exsiting db column match.  You can't compare blank to
 * null in postgres.  BEcause of this, I changed our copy of the dnarun.nmap file and removed "num"
 * as a mapping criteria.  I have left it in this code to alert me if I do a "pull" on GOBII
 * and move over new scripts.  Uploading DS_X.dnarun should reulst in no entries, which will
 * hopefully remind me to chaneg the mapping script again.
 *   
 * @author lcj34
 *
 */
public class MarkerDNARun_IFLFilePlugin extends AbstractPlugin {

    //  String inputFile() = "/Users/lcj34/notes_files/gobiiANDBms/gobii_curator_training/Maize282_noComments.hmp.txt";
    //  String outputFileDir() = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/gobii_ifl_files/";
    //  String platform_name = "GBSv27";

    private PluginParameter<String> dbConfigFile= new PluginParameter.Builder<>("dbConfigFile",null,String.class).guiName("dbConfigFile").required(true)
            .description("DB connection config file").build();
    private PluginParameter<String> inputFile= new PluginParameter.Builder<>("inputFile",null,String.class).guiName("inputFile").required(true)
            .description("Full path of file or directory, each file including the header line. Files with format *.hmp.txt, *.hmp.txt.gz,*.vcf, *.vcf.gz will be processed.").build();
    private PluginParameter<String> outputFileDir= new PluginParameter.Builder<>("outputFileDir",null,String.class).guiName("outputFileDir").required(true)
            .description("Directory where created files will be written.  Should end with /").build();
    private PluginParameter<String> refFile= new PluginParameter.Builder<>("refFile",null,String.class).guiName("Reference File").required(true)
            .description("Species reference file used to determine ref allele at marker position").build();
    private PluginParameter<String> mappingFile = new PluginParameter.Builder<>("mappingFile",null,String.class).guiName("mappingFile").required(true)
            .description("tab-delimited File containing columns: TaxaColumn, name, source,MGID, GID,libraryID, plate_code, well, species, type and project").build();
    private PluginParameter<String> mapsetName = new PluginParameter.Builder<>("mapsetName",null,String.class).guiName("Mapset Name").required(true)
            .description("mapset name from the mapset table.  Used to identify correct linkage group, e.g. chrom 1 from agpv2 vs chrom 1 from agpv3").build();
    private PluginParameter<String> expName= new PluginParameter.Builder<>("expName",null,String.class).guiName("Experiment Name").required(true)
            .description("Name of experiment to which this data belongs.  Must match an experiment name from the db experiment table.").build();
    private PluginParameter<String> platformName= new PluginParameter.Builder<>("platformName",null,String.class).guiName("Platform Name").required(true)
            .description("THe platform on which this data set was run, e.g. GBSv27.  Must match a platform name from the platform db table").build();
    private PluginParameter<String> refName= new PluginParameter.Builder<>("refName",null,String.class).guiName("Reference Name").required(true)
            .description("Name of referenece, e.g agpv2.  Must match name from entry in reference table in db.").build();
    private PluginParameter<String> datasetName= new PluginParameter.Builder<>("datasetName",null,String.class).guiName("Dataset Name").required(true)
            .description("Name of dataset for this data.  Must match the name of one of the administered datasets in the db.").build();

    static int datasetId = -1;
    static int projectId = -1;
    static int platformId = -1;
    static int mapsetId = -1;
    public MarkerDNARun_IFLFilePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    public MarkerDNARun_IFLFilePlugin() {
        super(null, false);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
 
    }

    @Override
    protected void postProcessParameters() {

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

        // THis method creates tab-delimited text files to be used for loading the specified tables
        DataOutputStream writerMarker = null;
        DataOutputStream writerMarkerProp = null; // we'll ignore this for now until Cinta gives us properites file
        DataOutputStream writerMarkerLink = null;
        DataOutputStream writerDSMarker = null;
        DataOutputStream writerVariants = null;
        DataOutputStream writerMarkerAlts = null;
        
        int marker_idx = 0; // incr by 1 as add in sequence to the marker table.  THis is how GOBII populates the field.
                     
        long totalTime=System.nanoTime();
        // Check if inputFile is file or directory
        File dataFile = new File(inputFile());
        if (!dataFile.exists()) {
            System.out.println("ERROR - input file doesn't exit: " +  inputFile());
            return null;
        }
        // Create list of files
        List<Path> directoryFiles = new ArrayList<>();
        String inputFileGlob="glob:*{hmp.txt,hmp.txt.gz,vcf,vcf.gz}";
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
        
        // Get db connection.  needed to query for table ids based on name
        Connection dbConnection = GOBIIDbUtils.connectToDB(dbConfigFile());
        if (dbConnection == null) {
            throw new IllegalStateException("MarkerDNARun_IFLFilePlugin:process_data: Problem connecting to database.");
        }
        
        // get the dataset id;
        StringBuilder sb = new StringBuilder();
        sb.append("select dataset_id from dataset where name = '");
        sb.append(datasetName());
        sb.append("';");
        datasetId = getTableId( dbConnection, sb.toString(), "dataset_id");
        
        // get the platformId
        sb.setLength(0);
        sb.append("select platform_id from platform where name='");
        sb.append(platformName());
        sb.append("';");
        platformId = getTableId( dbConnection, sb.toString(), "platform_id");
        
        // get the mapsetId
        sb.setLength(0);
        sb.append("select mapset_id from mapset where name='");
        sb.append(mapsetName());
        sb.append("';");        
        mapsetId = getTableId( dbConnection, sb.toString(), "mapset_id");
        
        // Get taxa mapping object.  THis is used for mapping dnarun and dataset_dnarun values
        HashMap<String,HmpTaxaData> taxaDataMap = createTaxaMap(dbConnection, mappingFile());
        if (taxaDataMap == null) return null;
        System.out.println("MarkerDNARunMGID: finished creating taxaDataMap, size: " + taxaDataMap.size());

        // The "name" field is not required by the db.  We require it here in
        // order to find the correct dataset, AND to have a consistent identifier
        // for this group of files.  Constraints on the name are NO SPACES - must be
        // alpha-numeric with underscore as term separator  (my constraints for both GOBII
        // IFL scripts and for this file)

        // Grab from the input file
        // THese need to be named DS_<dataset_id>.<table name>
        String markerOutFile = outputFileDir() + "DS_" + datasetId + ".marker";
        String markerPropOutFile = outputFileDir() + "DS_" + datasetId + ".marker_prop";
        String markerLinkOutFile = outputFileDir() + "DS_" + datasetId + ".marker_linkage_group";
        String dsMarkerOutFile = outputFileDir() + "DS_" + datasetId + ".dataset_marker";
        String dnarunOutFile = outputFileDir() + "DS_" + datasetId + ".dnarun"; // passed to createDNARunFiles
        String dsdnarunOutFile = outputFileDir() + "DS_" + datasetId + ".dataset_dnarun"; // passed to createDNARunFiles
        String variantOutFile = outputFileDir() + "DS_" + datasetId + ".variant";
        String markerAltsFile = outputFileDir() + "DS_" + datasetId + "_markerAlts.txt";
        BufferedReader markerbr = null;

        // Create reference file 
        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(refFile());
        byte[] refChromBytes = null; // holds ref bytes on a per-chrom basis

        // Create string builders for the 8 files - we'll append data, then write
        StringBuilder markerSB = new StringBuilder();
        StringBuilder markerPropSB = new StringBuilder();
        StringBuilder markerLinkageSB = new StringBuilder();
        StringBuilder dsMarkerSB = new StringBuilder();
        StringBuilder variantsSB = new StringBuilder();
        StringBuilder markerAltsSB = new StringBuilder();
        

        try {
            // GOBII currently doesn't have an IFL file for dataset, so that should be loaded by hand
            // prior to running this script.  Or, could consider connecting and adding it in here
            // Currently I create the dataset manually before running this plugin
            writerMarker = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(markerOutFile)));
            writerMarkerProp = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(markerPropOutFile)));
            writerMarkerLink = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(markerLinkOutFile)));
            writerDSMarker = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(dsMarkerOutFile)));
            writerVariants = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(variantOutFile)));
            writerMarkerAlts = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(markerAltsFile)));
            
            // write header lines - these are used when creating foreign table by the IFL scripts
            markerSB.append("platform_name\tvariant_id\tname\tcode\tref\talts\tsequence\treference_name\tprimers\tprobsets\tstrand_name\tstatus\n");
            writerMarker.writeBytes(markerSB.toString());
            markerPropSB.append("marker_name\tprops\n");
            writerMarkerProp.writeBytes(markerPropSB.toString());
            markerLinkageSB.append("marker_name\tplatform_id\tstart\tstop\tlinkage_group_name\tmap_id\n");
            writerMarkerLink.writeBytes(markerLinkageSB.toString());
            dsMarkerSB.append("dataset_name\tmarker_name\tplatform_id\tcall_rate\tmaf\treproducibility\tscores\tmarker_idx\n");
            writerDSMarker.writeBytes(dsMarkerSB.toString());
            markerAltsSB.append("name\tplatform_id\talts\n");
            //writerMarkerAlts.writeBytes(markerAltsSB.toString());
            
            // reset lengths to 0 after writing string
            markerSB.setLength(0);
            markerPropSB.setLength(0);
            markerLinkageSB.setLength(0);
            dsMarkerSB.setLength(0);
            variantsSB.setLength(0);
            markerAltsSB.setLength(0);
            
            int prevChrom = -1;
            int curChrom = -1;
            
            int[] tabPos = new int[11]; // there are many values, we care about the first 10           
            // Process all input files.  We could have 1 hmp.txt file, or there could be
            // a directory of them (generally split by chromosome).  Process all files on the list
            System.out.println("LCJ - size of directoryFiles: " + directoryFiles);
            for (int idx = 0; idx < directoryFiles.size(); idx++) {
                int totalLines = 0;
                long time=System.nanoTime();
                Path infile = directoryFiles.get(idx);
                String infileString = infile.toString();
                System.out.println("\nMarkerDNARun_IFLFilePlugin: processing file " + infileString);
                markerbr = Utils.getBufferedReader(infileString, 1 << 22);
                // check each file on the directlyFiles list.  If it ends with .vcf or .vcf.gz process vcf
                // if file ends with hmp.txt or hmp.txt.gz process hapmap.  Anything else should have been
                // tossed.
                boolean isVCF = false;
                boolean wroteHeader = false;
                if (infileString.endsWith("vcf.gz") || infileString.endsWith("vcf")) {
                    isVCF = true;
                }
                String mline;
                // Process each line in the file
                while ( (mline=markerbr.readLine() )!= null) {
                    totalLines++;
                    if (mline.startsWith("##")) continue; // toss comments, assumes all comments are at top of file, followed by header
                    if (!wroteHeader) {
                        // after comments we get the header - write the dnarun and dataset_dnarun files
                        // dnarun IFL maps experiment name and dnasample name.
                        // format is: experiment_name, dnasample_name, name,code
                        boolean dnaSuccess = createDNARunFiles( dbConnection,mline, dnarunOutFile,dsdnarunOutFile, expName(),
                                projectId, datasetName(), taxaDataMap,isVCF);                                
                        if (!dnaSuccess) {
                            System.out.println("ERROR processing dnarun and dataset_dnarun tables - time to quit!");
                            writerMarker.close();
                            writerMarkerProp.close();
                            writerMarkerLink.close();
                            writerDSMarker.close(); 
                            writerVariants.close();
                            writerMarkerAlts.close();
                            return null;
                        } 
                        wroteHeader = true;
                        System.out.println("LCJ - found header at line " + totalLines);
                        continue;
                    }
                                     
                 // Get tab positions in the string
                    int len = mline.length();
                    int tabIndex = 0;
                    for (int i = 0; (tabIndex < 11) && (i < len); i++) {
                        if (mline.charAt(i) == '\t') {
                            tabPos[tabIndex++] = i;
                        }
                    }
                    // Get the chromosome 
                    curChrom = GOBII_IFLUtils.getChromFromLine(mline, isVCF, tabPos);
     
                    if (curChrom < 1 || curChrom > 10) continue; // skipping all but chroms 1-10
                    if (curChrom != prevChrom) {
                        // get reference for this chromosome
                        System.out.println("LCJ - processing chromosome : " + curChrom);
                        Chromosome newChrom = new Chromosome(Integer.toString(curChrom));
                        try {
                            refChromBytes = myRefSequence.chromosomeSequence(newChrom);
                        } catch (Exception exc) {
                            System.out.println("LCJ - no data for chrom " + curChrom + " continuing ...");
                            continue;
                        }
                        if (refChromBytes == null) {
                            System.out.println("LCJ - NO BYTES found for chrom " + curChrom);
                            continue;
                        }
                        prevChrom = curChrom;
                    }
                    String linkageGroupName = Integer.toString(curChrom); // chrom is 3rd column

                    // Data for marker intermediary table :
                    // default mappings are reference name, platform name, strand name
                    // format is: platform_name/variant_id/name/code/ref/alts/sequence/reference_name/primers/probsets/strand_name/status
     
                    //platformname - will be mapped to platform id in db                   
                    markerSB.append(platformName());
                    markerSB.append("\t");
                    markerSB.append("\t"); //skip variant field
                    int position = GOBII_IFLUtils.getPosFromLine(mline, isVCF, tabPos);
                    String markerName = GOBII_IFLUtils.getMarkerNameFromLine(mline,isVCF,tabPos,mapsetName());
                    if (markerName == null) {
                        System.out.println("LCJ - failure from call to Gobii_IFLUtils.getMarkerNameFromLine !!!");
                        writerMarker.close();
                        writerMarkerProp.close();
                        writerMarkerLink.close();
                        writerDSMarker.close(); 
                        writerVariants.close();
                        writerMarkerAlts.close();
                        return null;
                    }

                    markerSB.append(markerName); // name field
                    markerSB.append("\t");
                    markerSB.append("dummycode\t"); // code field
                    
                    byte[] oneAllele = new byte[1];
                    // Position from hmp.txt or vcf file is 1 based, nucleotideByteToString wants 0 based
                    oneAllele[0] = refChromBytes[position-1];
                    String ref = NucleotideAlignmentConstants.nucleotideBytetoString(oneAllele);
                    //String ref = aTokens[0]; // assumes ref is first one - believe GOBII processes this way
                    markerSB.append(ref);
                    markerSB.append("\t");
                    
                    // Alts are A/C/G/T minus the ref                   
                    String alts = GOBII_IFLUtils.getAltsForRef(ref);
                    markerSB.append(alts);
                    markerSB.append("\t");
                    markerSB.append("\t"); // no sequence, just tab over
                    markerSB.append(refName());
                    markerSB.append("\t");
                    markerSB.append("\t\t"); // skip over primers and probsets - we leave these fields blank
                    // find the strand - column 5 in hmp.txt file, not present in VCF
                    String strand = GOBII_IFLUtils.getStrandFromLine(mline, isVCF, tabPos);
                    markerSB.append(strand);
                    markerSB.append("\t");
                    markerSB.append("1\n"); // last item in table, the status - default to 1
                    
                    // add the markerAlts values
                    // Used when we handle deletions - need actual alts from the data file
                    String altsFromFile = GOBII_IFLUtils.getAltsFromLine(mline, ref, isVCF, tabPos);
                    markerAltsSB.append(markerName);
                    markerAltsSB.append("\t");
                    markerAltsSB.append(platformId);
                    markerAltsSB.append("\t");
                    markerAltsSB.append(altsFromFile);
                    markerAltsSB.append("\n");
                    
                    // add the marker_prop Table entries
                    // format is: marker_name/props
                    markerPropSB.append(markerName);
                    markerPropSB.append("\t");
                    
                    // NOTE:  Currently the markerProp table is NOT populated.  I leave the
                    // creation of it here as an example of how to create the JSONB, but never 
                    // run the IFL scripts with the created file.  
                    //
                    // What entries from the hmp file do we want?  SHould add more to marker
                    // from the headers.  GOBII handles this mapping on the loader form.
                    // It has a table with specific marker properties and lets user fill in value.
                    // This would be an issue for us - we'd have to know what properties to add.
                    // Perhaps I keep a table somewhere that indicates type of input file,
                    // the group (AGPv2, AGPv3, etc) and can set them.
                    
                    // I'm going to default to species and genome_build and source
                    // Note the cv entries here are hard-coded. These are not valid if that table changes.
                    // How to handle this?  Could make db queries from here, or perhaps we don't fill in the prop table yet?
                    
                    // NOTE - you need all the escaped " shown below for the postgres json entries.
                    // the actual file will look like this:
                    // S1_10045        "{""23"": ""2"",""24"": ""Zea Mays"",""25"": ""Maize282_GBSv27_raw_MAF02""}"
                    String propsString = "\"{\"\"23\"\": \"\"2\"\",\"\"24\"\": \"\"Zea Mays\"\",\"\"25\"\": \"\"Maize282_GBSv27_raw_MAF02\"\"}\"";
                    
                    markerPropSB.append(propsString); // last column, no tab needed
                    markerPropSB.append("\n");
                    
                    // add marker_linkage_group_entries
                    // format is:  marker_name, platformId,start,stop, linkage_group_name (IFL maps to linkage_group_id)
                    markerLinkageSB.append(markerName);
                    markerLinkageSB.append("\t");
                    
                    markerLinkageSB.append(platformId); // store the platformId
                    markerLinkageSB.append("\t");
                    
                    // the position is both the start and stop - we grabbed it above when finding the ref
                    markerLinkageSB.append(position); // no need to convert to string?
                    markerLinkageSB.append("\t");
                    markerLinkageSB.append(position);
                    markerLinkageSB.append("\t");
     
                    // add linkage group id - IFL maps to linkage_group_id based on name
                    markerLinkageSB.append(linkageGroupName); // last column, no tab needed
                    markerLinkageSB.append("\t");
                    markerLinkageSB.append(mapsetId); // IFL uses linkageGroupName and mapset_id  together to get linkage_group_id
                    markerLinkageSB.append("\n"); // end of line
                    
                    // Write the dataset_marker file
                    // format is: dataset_name,marker_name,call_rate,maf,reproducibility,scores,marker_idx
                    // by default, dataset_marker.nmap only maps the marker_name.  I added dataset_name to this conversion file
                    // only populating the first 3 columns
                    dsMarkerSB.append(datasetName());
                    //int dsID = 2; // if using dataset_id vs name column .  Shouldn't need id anymore
                    //dsMarkerSB.append(Integer.toString(dsID)); // switch back if I get mapping to work
                    dsMarkerSB.append("\t");
                    dsMarkerSB.append(markerName);
                    dsMarkerSB.append("\t");
                    dsMarkerSB.append(platformId);
                    //dsMarkerSB.append("\t\t\t\t\t\n"); // skip rest of columns
                    dsMarkerSB.append("\t\t\t\t\t");
                    dsMarkerSB.append(Integer.toString(marker_idx));       
                    dsMarkerSB.append("\n");  // end of line
                    
                    marker_idx++; // Marker_idx takes values 0-N, incr after each value is added
                    
                    // Create the monetdb variants file with taxa info - no header
                    // This file is one of 3 files used for creating and loading
                    // the monetdb table for this dataset.  The other 2 files are
                    // generated from postgres queries of table IDs, and cannot be
                    // created until the files created from this class have been
                    // run through the IFL scripts and loaded.  THis is because the monetdb
                    // table creator script needs the dnarun_id and marker_id which will
                    // be generated by postgres upon loading those tables.
                    
                    String variantLine = GOBII_IFLUtils.addMonetdbVariantData(ref,altsFromFile,mline, isVCF, tabPos);
                    if (variantLine != null) {
                        variantsSB.append(variantLine);
                    } else {
                        System.out.println("LCJ - failure from call to Gobii_IFLUtils.addMonetdbVariantData !!!");
                        writerMarker.close();
                        writerMarkerProp.close();
                        writerMarkerLink.close();
                        writerDSMarker.close(); 
                        writerVariants.close();
                        writerMarkerAlts.close();
                        return null;
                    }
                                        
                    // Write lines to all files
                    if (totalLines >1000) {
                        writerMarker.writeBytes(markerSB.toString());
                        writerMarkerProp.writeBytes(markerPropSB.toString());
                        writerMarkerLink.writeBytes(markerLinkageSB.toString());
                        writerDSMarker.writeBytes(dsMarkerSB.toString());
                        writerVariants.writeBytes(variantsSB.toString()); 
                        //writerMarkerAlts.writeBytes(markerAltsSB.toString());
                        
                        // Reset strings to null before processing next batch
                        markerSB.setLength(0);
                        markerPropSB.setLength(0);
                        markerLinkageSB.setLength(0);
                        dsMarkerSB.setLength(0);
                        variantsSB.setLength(0);
                        markerAltsSB.setLength(0);
                        
                        totalLines = 0;
                    }
                }
                if (totalLines > 0) {
                    writerMarker.writeBytes(markerSB.toString());
                    writerMarkerProp.writeBytes(markerPropSB.toString());
                    writerMarkerLink.writeBytes(markerLinkageSB.toString());
                    writerDSMarker.writeBytes(dsMarkerSB.toString());
                    writerVariants.writeBytes(variantsSB.toString()); 
                   // writerMarkerAlts.writeBytes(markerAltsSB.toString());
                    
                    // Reset strings to null before processing next batch
                    markerSB.setLength(0);
                    markerPropSB.setLength(0);
                    markerLinkageSB.setLength(0);
                    dsMarkerSB.setLength(0);
                    variantsSB.setLength(0);
                    markerAltsSB.setLength(0);
                }
                System.out.println("Process took " + (System.nanoTime() - time)/1e9 + " seconds for file " + infileString);
                 
            }
            writerMarker.close();
            writerMarkerProp.close();
            writerMarkerLink.close();
            writerDSMarker.close(); 
            writerVariants.close();
            writerMarkerAlts.close();
            
            if (markerbr != null) {
                markerbr.close();
            }          
        } catch (IOException ioe) {
            System.out.println("Caugh exception reading or writing hmp.txt files");
            ioe.printStackTrace();
        }
        System.out.println("Total time to process all files: " + (System.nanoTime() - totalTime)/1e9 + " seconds ");
        return null;
    }
    
    public HashMap<String,HmpTaxaData> createTaxaMap (Connection conn, String mappingFile) {
        HashMap<String,HmpTaxaData> taxaDataMap = new HashMap<>();
        // open the mapping file
        BufferedReader mappingbr = Utils.getBufferedReader(mappingFile);
        
        // read the mapping file, create hashmap with name as the key,
        // then an object of MGID, libraryID.  Project name, experiment,platform, reference and dataset
        // names should be the same for all in this group and are separate parameters
        // to this method.  This allows the same format for germplasm/dnasample as we have
        // for individual experiment/dataset runs.
        
        System.out.println("MarkerDNARunMGID: creating taxaDataMap from mapping file");
        try {
         // column names: storing idxes so column names may appear in any order.
            int taxaIdx=-1, nameIdx=-1, sourceIdx=-1, mgidIdx=-1, gidIdx=-1, libIdx=-1;
            int plateIdx=-1, wellIdx=-1, speciesIdx=-1, typeIdx=-1, projectIdx = -1, sampleIdx = -1;
            String mappingLine = mappingbr.readLine(); // header line
            String [] headers = mappingLine.split("\\t");
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
                    sampleIdx = idx;
                }
                idx++;
            }
            // LCJ _ return this when files are fixed to include sample name - WGS is not fixed!!
            if (taxaIdx == -1 || nameIdx == -1 || sourceIdx == -1 || mgidIdx == -1 || gidIdx == -1 || libIdx == -1 ||
                    plateIdx == -1 || wellIdx == -1 || speciesIdx == -1 || typeIdx == -1 || projectIdx == -1 || sampleIdx == -1) {
                   System.out.println("Mappingfile is missing required header line.  Expecting columns: TaxaColumn, name, source, MGID, GID, libraryID, plate_code, well, species, type, project, SampleName");
                   return null;
            }
//            if (taxaIdx == -1 || nameIdx == -1 || sourceIdx == -1 || mgidIdx == -1 || gidIdx == -1 || libIdx == -1 ||
//                    plateIdx == -1 || wellIdx == -1 || speciesIdx == -1 || typeIdx == -1 || projectIdx == -1 ) {
//                   System.out.println("Mappingfile is missing required header line.  Expecting columns: TaxaColumn, name, source, MGID, GID, libraryID, plate_code, well, species, type, project, SampleName");
//                   return null;
//            }
            boolean first = true;
            while ((mappingLine = mappingbr.readLine()) != null) {
                String[] data = mappingLine.split("\\t");
                if (first) {
                    // get the project id from project name - should be same for all taxa in file
                    // get the dataset id;
                    StringBuilder sb = new StringBuilder();
                    sb.append("select project_id from project where name = '");
                    sb.append(data[projectIdx].trim());
                    sb.append("';");
                    projectId = getTableId( conn, sb.toString(), "project_id");
                    first=false;
                }
                
                String taxa = data[taxaIdx].trim(); // taxa name to match from input data file
                String mgid = data[mgidIdx].trim();
                String libID = data[libIdx].trim();
                String plateName = data[plateIdx].trim();
                String well = data[wellIdx].trim();
                String gid = data[gidIdx].trim();
                String sampleName = null;                
                if (sampleIdx == -1 || data[sampleIdx].trim().equals("")) {                   
                    sampleName = data[gidIdx].trim() + ":" + data[plateIdx].trim() + ":" + data[wellIdx].trim();
                } else {
                    sampleName = data[sampleIdx].trim();
                }
//                if (sampleIdx == -1) {
//                    sampleName = data[gidIdx].trim() + ":" + data[plateIdx].trim() + ":" + data[wellIdx].trim();
//                } else {
//                    sampleName = data[sampleIdx].trim();
//                }
                
                HmpTaxaData taxaDataItem = new HmpTaxaData(mgid, gid, libID,plateName,well,sampleName);
                taxaDataMap.put(taxa, taxaDataItem);           
            }
        } catch (IOException ioe) {
            System.out.println("Caught exception reading mapping file");
            ioe.printStackTrace();
            return null;
        }
        return taxaDataMap;
    }
    
    private static int getTableId(Connection conn, String query, String column) {
        try {
            ResultSet rs = conn.createStatement().executeQuery(query);
            while (rs.next()) {
                int id = rs.getInt(column);
                return id;
            }
            
        } catch (SQLException sqle) {
            System.out.println("getTableId barfed on query: " + query);
            sqle.printStackTrace();
            return -1;
        }
       return -1;
    }
    private static boolean createDNARunFiles(Connection dbConnection, String hdrLine, String dnarunOutFile,
            String dsdnarunOutFile,String expName, int projectId, String dsName, 
            HashMap<String,HmpTaxaData> taxaDataMap, boolean isVCF) {
        // create the dnarun and dataset_dnarun IFL files
        System.out.println("LCJ - createDNARunFiles - begin ");
        int dnarun_idx = 0;
        boolean returnVal = true;
        try {
            // Get experiment id - needed together with name for dataset_dnarun to get dnarun ID
            StringBuilder builder = new StringBuilder();
            builder.append("select experiment_id from experiment where name= '");
            builder.append(expName);
            builder.append("';");
            
            int experiment_id = getTableId( dbConnection, builder.toString(), "experiment_id");
            if (experiment_id == -1) return false;
            
            DataOutputStream writerDNARun = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(dnarunOutFile)));
            DataOutputStream writerDSdnaRun = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(dsdnarunOutFile)));
            String[] hdrTokens = hdrLine.split("\t");
            
            System.out.println("LCJ - createDNARunFiles: size of hdrTokens: " + hdrTokens.length);
            StringBuilder dsDNArunSB = new StringBuilder();
            StringBuilder dnaRunSB = new StringBuilder();
            
            // create the header line dnasample_name,num,project_id (previously used platename vs num field - IFL changed)
            // NOTE:  BL doesn't use "num" and it has been removed from our copy of the loading scripts.  The null values
            // would not compare in postges (need IS NULL, GOBII IFL scripts hard to add for general query)
           // dnaRunSB.append("experiment_name\tdnasample_name\tnum\tproject_id\tname\tcode\n");
            dnaRunSB.append("experiment_name\tdnasample_name\tproject_id\tname\tcode\n");
            writerDNARun.writeBytes(dnaRunSB.toString());
            
            // create header line for dataset_dnarun
            dsDNArunSB.append("dataset_name\tdnarun_name\texperiment_id\tdnarun_idx\n");
            writerDSdnaRun.writeBytes(dsDNArunSB.toString());
            // The initial tokens are skipped
            int idx = 11; // for hmp, taxa start in 12th column (11 when 0-based)
            if (isVCF) idx= 9;  // for vcf, taxa starts at 10th column (9 when 0-based)
            
            for ( ; idx < hdrTokens.length; idx++) {
                //System.out.println("LCJ - createDNARUnFiles: processing taxa: " + hdrTokens[idx]);
                dsDNArunSB.setLength(0);
                dnaRunSB.setLength(0);
                
                HmpTaxaData taxaData =  taxaDataMap.get(hdrTokens[idx].trim());
                if (taxaData == null) {
                    System.out.println("LCJ - createDNARunFiles - NO DATA FOR taxa " + hdrTokens[idx]);
                    returnVal = false;
                    continue;
                }
                // create line for dnarun
                dnaRunSB.append(expName);
                dnaRunSB.append("\t");
                // add dnasample_name.  IFL gets dnasample_id from dnasmample name(extractionid, now "SampleName" in mapping file), num (null for us), project_id=project_id
                // Used to use plateName instead of num.  We don't always have either one.  HOpefully it just matches on the null entry
                dnaRunSB.append(taxaData.getDnasampleName()); // name
                dnaRunSB.append("\t");
                // LCJ - not appending getNum() as postgres does not consider "" to be equivalent to NULL
                // SO these fields don't match when running the IFL scripts and no data is loaded.
                //dnaRunSB.append(taxaData.getNum()); // should be "" - BL doesn't use this.  
                //dnaRunSB.append("\t"); // took out column and data now! don't need extra \t
                dnaRunSB.append(projectId);
                dnaRunSB.append("\t");
                // Programmer's meeting on 7/7/16: decided dnarun name is just library prep ID, NOT TAXA
                String dnarun_name = taxaData.getLibraryID();
                dnaRunSB.append(dnarun_name); // name field - store just the library prep ID
                dnaRunSB.append("\t");
               // dnaRunSB.append(taxaData.getMGID()); // MGID stored in "code" field
                dnaRunSB.append("dummycode");
                dnaRunSB.append("\n"); // that's it - end the line
                writerDNARun.writeBytes(dnaRunSB.toString());
                
                // create line for dataset_dnarun
                // format: datset_name, dnarun_name, experiment_id, dnarun_idx
                dsDNArunSB.append(dsName); // IFL maps to dataset_id (I added this mapping )
                dsDNArunSB.append("\t");
                dsDNArunSB.append(dnarun_name); // IFL maps name to dnarun_id
                dsDNArunSB.append("\t");
                dsDNArunSB.append(experiment_id);
                //dsDNArunSB.append("\t\n"); // skip the last field
                // add the dnarun_idx value
                dsDNArunSB.append("\t");
                dsDNArunSB.append(Integer.toString(dnarun_idx));
                dsDNArunSB.append("\n"); // end of line - add new line
                
                dnarun_idx++; // next dnarun_idx
                
                writerDSdnaRun.writeBytes(dsDNArunSB.toString());                
            }
            writerDNARun.close();
            writerDSdnaRun.close();
            
        } catch (IOException ioe) {
            System.out.println("LCJ - error processing IFL files for dnarun or dataset_dnarun table");
            ioe.printStackTrace();
            return false;
        }
        System.out.println("LCJ - successful creation of DNARun and dataset_dnarun files");
        //return true; // successful processing
//        if (returnVal == false) {
//            System.out.println("LCJ - hdrline has these taxa: ");
//            System.out.println(hdrLine);
//        }
        return returnVal;
        //return false;  // LCJ - RETURN LINE ABOVE !! _ This is just to force quit after creating dnarun tables.
    }
    public static void main(String[] args) {
        GeneratePluginCode.generate(MarkerDNARun_IFLFilePlugin.class);
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(MarkerDNARun_IFLFilePlugin.class);
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
    public MarkerDNARun_IFLFilePlugin dbConfigFile(String value) {
        dbConfigFile = new PluginParameter<>(dbConfigFile, value);
        return this;
    }
    /**
     * hmp.txt file including, including the header line,
     * which will be used to create marker related and dnarun
     * related intermediary files for GOBII loading
     *
     * @return inputFile
     */
    public String inputFile() {
        return inputFile.value();
    }

    /**
     * Set inputFile. hmp.txt file including, including the
     * header line, which will be used to create marker related
     * and dnarun related intermediary files for GOBII loading
     *
     * @param value inputFile
     *
     * @return this plugin
     */
    public MarkerDNARun_IFLFilePlugin inputFile(String value) {
        inputFile = new PluginParameter<>(inputFile, value);
        return this;
    }

    /**
     * Directory where created files will be written
     *
     * @return outputFileDir
     */
    public String outputFileDir() {
        return outputFileDir.value();
    }

    /**
     * Set outputFileDir. Directory where created files will
     * be written
     *
     * @param value outputFileDir
     *
     * @return this plugin
     */
    public MarkerDNARun_IFLFilePlugin outputFileDir(String value) {
        outputFileDir = new PluginParameter<>(outputFileDir, value);
        return this;
    }

    /**
     * Species reference file used to determine ref allele
     * at marker position
     *
     * @return Reference File
     */
    public String refFile() {
        return refFile.value();
    }

    /**
     * Set Reference File. Species reference file used to
     * determine ref allele at marker position
     *
     * @param value Reference File
     *
     * @return this plugin
     */
    public MarkerDNARun_IFLFilePlugin refFile(String value) {
        refFile = new PluginParameter<>(refFile, value);
        return this;
    }

    /**
     * tab-delimited File containing columns for taxaname,
     * name, MGID, libraryID, project_id, experiment_name,
     * platform_name, reference_name and dataset_name
     *
     * @return mappingFile
     */
    public String mappingFile() {
        return mappingFile.value();
    }

    /**
     * Set mappingFile. tab-delimited File containing columns
     * for taxaname, name, MGID, libraryID, project_id, experiment_name,
     * platform_name, reference_name and dataset_name
     *
     * @param value mappingFile
     *
     * @return this plugin
     */
    public MarkerDNARun_IFLFilePlugin mappingFile(String value) {
        mappingFile = new PluginParameter<>(mappingFile, value);
        return this;
    }

    /**
     * Integer identifying the mapset_id from the linkage group
     * table to use when mapping to marker_linkage_group.
     *
     * @return mapsetId
     */
    public String mapsetName() {
        return mapsetName.value();
    }

    /**
     * Set mapsetId
     * 
     * @param value mapsetId
     *
     * @return this plugin
     */
    public MarkerDNARun_IFLFilePlugin mapsetName(String value) {
        mapsetName = new PluginParameter<>(mapsetName, value);
        return this;
    }
    

    /**
     * Name of experiment to which this data belongs.  Must
     * match an experiment name from the db experiment table.
     *
     * @return Experiment Name
     */
    public String expName() {
        return expName.value();
    }

    /**
     * Set Experiment Name. Name of experiment to which this
     * data belongs.  Must match an experiment name from the
     * db experiment table.
     *
     * @param value Experiment Name
     *
     * @return this plugin
     */
    public MarkerDNARun_IFLFilePlugin expName(String value) {
        expName = new PluginParameter<>(expName, value);
        return this;
    }

    /**
     * THe platform on which this data set was run, e.g. GBSv27.
     *  Must match a platform name from the platform db table
     *
     * @return Platform Name
     */
    public String platformName() {
        return platformName.value();
    }

    /**
     * Set Platform Name. THe platform on which this data
     * set was run, e.g. GBSv27.  Must match a platform name
     * from the platform db table
     *
     * @param value Platform Name
     *
     * @return this plugin
     */
    public MarkerDNARun_IFLFilePlugin platformName(String value) {
        platformName = new PluginParameter<>(platformName, value);
        return this;
    }

    /**
     * Name of referenece, e.g agpv2.  Must match name from
     * entry in reference table in db.
     *
     * @return Reference Name
     */
    public String refName() {
        return refName.value();
    }

    /**
     * Set Reference Name. Name of referenece, e.g agpv2.
     *  Must match name from entry in reference table in db.
     *
     * @param value Reference Name
     *
     * @return this plugin
     */
    public MarkerDNARun_IFLFilePlugin refName(String value) {
        refName = new PluginParameter<>(refName, value);
        return this;
    }

    /**
     * Name of dataset for this data.  Must match one the
     * name of one of the administered datasets in the db.
     *
     * @return Dataset Name
     */
    public String datasetName() {
        return datasetName.value();
    }

    /**
     * Set Dataset Name. Name of dataset for this data.  Must
     * match one the name of one of the administered datasets
     * in the db.
     *
     * @param value Dataset Name
     *
     * @return this plugin
     */
    public MarkerDNARun_IFLFilePlugin datasetName(String value) {
        datasetName = new PluginParameter<>(datasetName, value);
        return this;
    }

    public static class HmpTaxaData  {
        //While these values are acutally ints, they will come in
        // as strings from the mapping file, and need to be stored
        // in the output file as strings, so they're declared as 
        // strings to save on conversion processing
        
        private  String MGID;
        private String libraryID;
        private String plateName;
        private String dnasampleName;
        private String num;
        
        HmpTaxaData (String MGID, String GID, String libraryID, String plateName, String well,String sampleName) {                       
            this.MGID = MGID;
            this.libraryID = libraryID;
            this.plateName = plateName;
            // THis will be changed to extraction_id.  That will be passed in
            // this.dnasampleName = GID + ":" + plateName + ":" + well;
            // because it keeps changing, we now have a SampleName column, the value of which is passed in here.
            this.dnasampleName = sampleName;
            this.num = ""; // Currently, BL doesn't use this field but is needed for GOBII IFL mapping file
        } 
        public String getMGID() {
            return this.MGID;
        }
        public String getLibraryID() {
            return this.libraryID;
        }
        public String getPlateName() {
            return this.plateName;
        }  
        public String getDnasampleName() {
            return this.dnasampleName;
        } 
        public String getNum() {
            return this.num;
        }
    }
}
