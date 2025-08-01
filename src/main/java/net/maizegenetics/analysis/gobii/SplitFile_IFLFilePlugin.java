/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import java.awt.Frame;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.nio.file.Path;
import java.nio.file.Paths;

import javax.swing.ImageIcon;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

/**
 * @author lcj34
 *
 */
public class SplitFile_IFLFilePlugin extends AbstractPlugin {
    private PluginParameter<String> inFile= new PluginParameter.Builder<>("inFile",null,String.class).guiName("Input FIle").required(true)
            .description("Name of file to split.  It must have 1 and only 1 period in the name to separate dataset name from table name").build();
    private PluginParameter<Integer> maxSize= new PluginParameter.Builder<>("maxSize",null,Integer.class).guiName("Maximum file size").required(true)
            .description("Maximum size of each file after splitting.").build();
    private PluginParameter<String> outputDir= new PluginParameter.Builder<>("outputDir",null,String.class).guiName("Path of output directory").required(true)
            .description("Full path name of output directory, must end with a /").build();
    
    public SplitFile_IFLFilePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    public SplitFile_IFLFilePlugin() {
        super(null, false);
    }
    
    @Override
    public DataSet processData(DataSet input) {
        BufferedReader br = Utils.getBufferedReader(inFile(), 1 << 22);
        //DataOutputStream writerMarker = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile)));
        DataOutputStream bw = null;
        StringBuilder sb = new StringBuilder();
        int fileCount = 1;
        Path filePath=Paths.get(inFile());
        String fileName = filePath.getFileName().toString();
        //String fileNameSubstr = fileName.substring(0,fileName.indexOf("."));
 
        
        String[] fileparts = fileName.split("\\."); // should be dataset_name.<tableName> - only 1 period in file name !!
        System.out.println("FileName is " + fileName + ", fileparts.length " + fileparts.length + "\n");
        if (fileparts.length != 2) {
            System.out.println("Input file must have 1 and only 1 period in the name to differentiate dataset name from table name");
            System.out.println("Inputfile example:  DS_4.marker");
            return null;
        }
        try {
            String line = br.readLine(); // should be the header
            String headerLine = line + "\n";
            int linecount = 0;
            int totalLines = 0;
            int maxBuffer = 10000; // max lines before we write out the buffer.
            String outfile = outputDir() + fileparts[0] + "_" + fileCount + "." + fileparts[1];
            // add header to file
            bw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile)));
            bw.writeBytes(headerLine); // already added \n to header line
            fileCount++;
            while ((line = br.readLine()) != null){
                sb.append(line);
                sb.append("\n"); // readline removes the carriage return
                linecount++;
                totalLines++;
                if (totalLines == maxSize() || linecount == maxBuffer) {
                    // write to file, clear out buffer
                    bw.writeBytes(sb.toString());
                    sb.setLength(0); // zero out the buffer                   
                    if (totalLines == maxSize()) {
                        // close the writer, get new one for next set
                        bw.close();
                        System.out.println("Wrote file " + outfile + " with lines " + totalLines);
                        outfile = outputDir() + fileparts[0] + "_" + fileCount + "." + fileparts[1]; // filecount is different
                        bw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile)));
                        bw.writeBytes(headerLine);
                        fileCount++; 
                        totalLines = 0;
                    } 
                    linecount = 0; // set linecount to 0, start filling buffer to write
                }               
            }
            if (linecount > 0) {
                // write remaining lines
                bw.writeBytes(sb.toString());
                System.out.println("Wrote remaining lines to file " + outfile + " with lines " + linecount);
            }
            bw.close();  // close out buffer, return
                    
        } catch (Exception exc) {
            System.out.println("Whoops - error reading/writing file " + inFile() + " or outfile " + fileCount);
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

//    public static void main(String[] args) {
//        GeneratePluginCode.generate(SplitFile_IFLFilePlugin.class);
//    }
    
    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
//    public <Type> runPlugin(DataSet input) {
//        return (<Type>) performFunction(input).getData(0).getData();
//    }

    /**
     * Name of file to split.  It must have 1 and only 1 period
     * in the name to separate dataset name from table name
     *
     * @return Input FIle
     */
    public String inFile() {
        return inFile.value();
    }

    /**
     * Set Input FIle. Name of file to split.  It must have
     * 1 and only 1 period in the name to separate dataset
     * name from table name
     *
     * @param value Input FIle
     *
     * @return this plugin
     */
    public SplitFile_IFLFilePlugin inFile(String value) {
        inFile = new PluginParameter<>(inFile, value);
        return this;
    }

    /**
     * Maximum size of each file after splitting.
     *
     * @return Maximum file size
     */
    public Integer maxSize() {
        return maxSize.value();
    }

    /**
     * Set Maximum file size. Maximum size of each file after
     * splitting.
     *
     * @param value Maximum file size
     *
     * @return this plugin
     */
    public SplitFile_IFLFilePlugin maxSize(Integer value) {
        maxSize = new PluginParameter<>(maxSize, value);
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
    public SplitFile_IFLFilePlugin outputDir(String value) {
        outputDir = new PluginParameter<>(outputDir, value);
        return this;
    }
}
