/**
 */
package net.maizegenetics.analysis.monetdb;

import java.awt.Frame;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GenomeSequence;
import net.maizegenetics.dna.map.GenomeSequenceBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;

/**
 * This plugin is copied from the ColumnsToBinarySNPOnlyTablePlugin (which was copied Kelly's
 * ColumnsToBinaryPlugin, which re-worked  various lynn methods) to create 
 * binary files for loading into the hmp321_snp table in the maizeFullGEnomeDB of the
 * Rare Alleles monetdb instance.  THe difference between this plugin and the 
 * ColumnsToBinarySNPOnlyTablePlugin is the latter creates entries only for identified hmp321 SNPs
 * This method (ColumnsToBinaryFullGenomePlugin) creates entries  for all positions listed in 
 * in the reference genome.
 * 
 * To be consistent with the existing columns in monetdb maize tables, the reference file must
 * be a link to a copy of the Zea_mays.AGPv3.20.dna.genome.fa file stored 
 * on andersonii in Research/Zea/Genotypes/Annotations/monetDB/refGenomeFiles.
 * 
 * The "inputFile" parameter can be either a single file containing data for all chromosomes or a 
 * directory of files split by chromosome.  If all chroms are in a single file, that file
 * must be sorted by chromosome and position, and must contain a header line that contains
 * the columns "chr" and "pos" along with user specified data columns.
 * 
 * If the "inputFile" parameter is a directory, the only *.txt files it holds must be files
 * intended for this processing.  These files must be split by chromosome and must be named
 * such that they will sort lexicographically in chromosome order.  For example:  files named
 * chr01.txt, chr02.txt ... chr09.txt, chr10.txt will sort from 1-10.  But files named
 * chr1.txt, chr2,txt ... chr9.txt, chr10.txt will not.  In the latter case, chr10.txt will
 * be processed before the other files.
 * 
 * This plugin may also be used to create binaries for the maizeChrom10DB.  In this case,
 * we still use the full reference genome, but only chrom 10 is processed.  The inputFile
 * paramaeter shoudl be just 1 file containing chromosome 10 data.
 * 
 * @author kelly
 * @author lcj34
 *
 */
public class ColumnsToBinaryFullGenomeTablePlugin extends AbstractPlugin {
    private PluginParameter<String> inputFile= new PluginParameter.Builder<>("inputFile",null,String.class).guiName("Input File").required(true)
            .description("Input File containing a header line and entries sorted by chromosome and position, or Directory containing Tab-delimited files split by chromosome and sorted by position.  \nIf parameter is a directory, each file must contain a header line, and the files must end with .txt and be named in a maaner that sorts lexicographically by chromosome").build();
    private PluginParameter<String> refFile= new PluginParameter.Builder<>("refFile",null,String.class).guiName("refFile").required(true)
            .description("A link to the species reference file. This must be the same reference file used by other columns in the monetdb table for your species.").build();
    private PluginParameter<String> outBase= new PluginParameter.Builder<>("outBase",null,String.class).guiName("outBase").required(true)
            .description("Output directory and base filename to hold the binary files. Will make directory if neccesary").build();
    private PluginParameter<String> colsFloat= new PluginParameter.Builder<>("colNamesFloat",null,String.class).guiName("Columns keep as Real (Float)").required(false)
            .description("Comma separated list of column names to generate real binaries for").build();
    private PluginParameter<String> colsInt= new PluginParameter.Builder<>("colNamesInt",null,String.class).guiName("Columns keep as Int").required(false)
            .description("Comma separated list of column names to generate int  binaries for").build();
    private PluginParameter<String> colsShort= new PluginParameter.Builder<>("colNamesShort",null,String.class).guiName("Columns keep as Short").required(false)
            .description("Comma separated list of column names to generate short binaries for").build();
    private PluginParameter<String> colsLong= new PluginParameter.Builder<>("colNamesLong",null,String.class).guiName("Columns keep as Long").required(false)
            .description("Comma separated list of column names to generate long binaries for").build();
    private PluginParameter<String> colsChar= new PluginParameter.Builder<>("colNamesChar",null,String.class).guiName("Columns keep as Char String").required(false)
            .description("Comma separated list of column names to generate character/text binaries for").build();
    private PluginParameter<String> colsByte= new PluginParameter.Builder<>("colNamesByte",null,String.class).guiName("Columns keep as Byte").required(false)
            .description("Comma separated list of column names to generate byte binaries for").build();
    private PluginParameter<String> colsAllele= new PluginParameter.Builder<>("colNamesAllele",null,String.class).guiName("Columns translate Alleles").required(false)
            .description("Comma separated list of column names to generate for single char alleles to be translated to 0-5 for A/C/G/T/+/-").build();
    private PluginParameter<String> colsLog10= new PluginParameter.Builder<>("colNamesLog10",null,String.class).guiName("Columns to Keep and transform -log10").required(false)
            .description("Comma separated list of column names to first transform using -log10 then generate binaries for").build();
    private PluginParameter<Boolean> range= new PluginParameter.Builder<>("range",false,Boolean.class).guiName("Range information?").required(false)
            .description("Columns for range data. If true, will look for 'start' and 'end' (inclusive exclusive) or 'first' 'last' (inclusive inclusive) instead of 'Pos' ").build();
    private PluginParameter<Boolean> negToZero= new PluginParameter.Builder<>("negToZero",false,Boolean.class).guiName("Negative floats to zero?").required(false)
            .description("Will set negative column values to zero (otherwise they the negative value is stored").build();
    private PluginParameter<Boolean> missToZero= new PluginParameter.Builder<>("missToZero",false,Boolean.class).guiName("Missing float values to zero?").required(false)
            .description("Will set missing  column values to zero (otherwise they are set to null").build();
    private PluginParameter<Boolean> oneBased= new PluginParameter.Builder<>("oneBased",true,Boolean.class).guiName("positions are 1 based?").required(false)
            .description("Will assume all positions and ranges are 1-based unless this value is set to false").build();
    
    private HashMap<Integer,String> log10HM= null;
    private HashMap<Integer,LittleEndianDataOutputStream> log10Writers= null;
    private HashMap<Integer,String> realHM= null;
    private HashMap<Integer,LittleEndianDataOutputStream> realWriters= null;
    private HashMap<Integer,String> shortHM= null;
    private HashMap<Integer,LittleEndianDataOutputStream> shortWriters= null;
    private HashMap<Integer,String> longHM= null;
    private HashMap<Integer,LittleEndianDataOutputStream> longWriters= null;
    ///PUT CHARACTER BACK (figure out what monetDB wants)
    private HashMap<Integer,String> charHM= null; //Strings written as Strings, because that's what MonetDB likes
    private HashMap<Integer,LittleEndianDataOutputStream> charWriters= null;
    private HashMap<Integer,String> intHM= null;
    private HashMap<Integer,LittleEndianDataOutputStream> intWriters= null;
    private HashMap<Integer,String> byteHM= null;
    private HashMap<Integer,LittleEndianDataOutputStream> byteWriters= null;
    private HashMap<Integer,String> alleleHM= null;
    private HashMap<Integer,LittleEndianDataOutputStream> alleleWriters= null;
    private int chrCol= -1;
    private String chrName= "CHR";
    private int posCol= -1;
    private int totalLines;
    private int totalZeroLines;
    private int totalNonZeroLines;
    private int linesForChr;
    private FileInputStream fis= null;
    private Scanner scanner= null;
    private boolean inclusive= false;
    private int startCol= -1;
    private int endCol= -1;
    
    private GenomeSequence myRefSequence = null;
    
    List<Path> infiles;
    
    public ColumnsToBinaryFullGenomeTablePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    public ColumnsToBinaryFullGenomeTablePlugin() {
        super(null, false);
    }
    
    @Override
    protected void preProcessParameters(DataSet input) {
    }
    
    @Override
    protected void postProcessParameters() {
 
        File out= new File(outBase.value());
        if (!out.getParentFile().exists()) out.getParentFile().mkdirs();
        
        // create list of directory files
        File dirList = new File(inputFile());
        if (!(dirList.exists())) {
            throw new IllegalStateException("Input file or directory not found !!") ;
        }
 
        if (dirList.isDirectory()) {
            String inputFileGlob="glob:*{.txt}";
            infiles= DirectoryCrawler.listPaths(inputFileGlob, Paths.get(inputFile.value()).toAbsolutePath());
            if (infiles.size() < 1) {
                throw new IllegalStateException("no .txt files found in input directory !!");
            }
        } else {
            infiles=new ArrayList<>();
            Path filePath= Paths.get(inputFile()).toAbsolutePath(); 
            infiles.add(filePath);
        }
        Collections.sort(infiles); // files should be chr01.txt ...chr10.txt so they will sort correctly
        // all chroms should have same file headers - get columns from first file
        // This is a requirement to have headers in the files.
        FindColumns(infiles.get(0).toString());
    }

    @Override
    public DataSet processData(DataSet input) {
        long time= System.nanoTime();
        Initialize();
        
        // Create reference genome
        myRefSequence = GenomeSequenceBuilder.instance(refFile());
        try {
            String currChr = null;
            Chromosome currChrom = null; // for getting reference sequence
            System.out.println("processData: number of infiles: " + infiles.size());
            for (Path filepath: infiles) {
                String inputFile = filepath.toString();
                int refChromPosSize = 0;
                int lastStart= -1; int lastEnd= -1;
                // Full genome lastLinePos defaults to 1, not 0, as our ref genome is 1-based.
                int lastLinePos = 1; //lastLinePos holds the first position not already written. 
                boolean first = true;
                fis = new FileInputStream(inputFile);
                System.out.println("Processing file: " + inputFile);
                scanner = new Scanner(fis);  
 
                while (scanner.hasNextLine()) {
                    String line = scanner.nextLine();
                    if ((line.isEmpty() || line.startsWith("#")) && first) continue; // skip comment lines
                    if (line.isEmpty() && !first) break;
                    String[] next = line.split("\t");
                    if (next[chrCol].toUpperCase().equals(chrName)) {continue;} // this is header line - skip it
                    if (first) { //initialize the reference position list                       
                        currChr = next[chrCol];
                        System.out.println("...working on chr "+currChr);
                        currChrom = new Chromosome(currChr);
                        refChromPosSize = myRefSequence.chromosomeSize(currChrom);
                        first= false;
                    }
 
                    //When the current line changes chromosome, finish writing end of last, reinitialize for new
                    // The code below is not hit when files are broken up by chromosome.  It IS hit if all
                    // chromosome data is contained in a single file.
                    if (!currChr.equals(next[chrCol])) {
                        // Write everything to null between the last line and the end of the chromosome
                        writeNull(lastLinePos,refChromPosSize);
                        System.out.println("Chrom within file changed, " + linesForChr+" lines output for chr "+currChr);
                        linesForChr= 0; lastLinePos = 1;
                        //if chromosomes not ordered 1-10, shut everything down and throw exception
                        if (Integer.parseInt(currChr)!=(Integer.parseInt(next[chrCol])-1)) {
                            Shutdown();
                            String message = "Chromosomes files must be named so they sort by chromosomes in ascending order !!\nYour chromosome " +
                               currChr + " came before chromosome " + next[chrCol];
                            throw new IllegalStateException(message);
                        }
                        currChr= next[chrCol]; lastStart= -1; lastEnd= -1;
                        currChrom = new Chromosome(currChr);
                        refChromPosSize = myRefSequence.chromosomeSize(currChrom);
                        System.out.println("...working on chr "+currChr);
                    }
                    int position= -1; //problems if in exponential notation (sometimes R does this) so do it as a try as double and cast if problems
                    try {
                        position= (posCol>-1)?Integer.parseInt(next[posCol]):Integer.parseInt(next[startCol]);
                        if (!oneBased()) { 
                            position++; // increment position so it matches the 1-based stored in reference files
                        }
                    } catch(Exception e) {
                        System.out.println(line); 
                        position= Double.valueOf(next[posCol]).intValue();
                        if (!oneBased()) position++; // increment if input file has 0-based positions
                    }
                    int currPosOnChr= position; // 1-based position
                    //Only take the first one, if position/chromosome is duplicated, as long as lines are (if range, just checks start)
                    if (lastStart==currPosOnChr) {
                        System.out.println("WARNING: Already recorded value for this position. Excluded\n"+line); continue;
                    }
                    if (range.value() && inclusive && lastEnd==currPosOnChr) { //if inclusive and lastEnd position same as start of next throw exception
                        throw new IllegalStateException("Ranges are overlapping and supposed to be inclusive!!!!!. If not inclusive, denote ranges by start, end");
                    }
                    // Write everything to null between the last line and the current line
                    // This only writes null if the start (lastLinePos) is less than the end (currPosOnChr)
                    // Because we are doing this all as 1-based, it must be set to 1 initially, not 0.
                    // SNP-only uses 0 as it is indexing an array, which is 0-based.
                    writeNull(lastLinePos,currPosOnChr-1); 
                    //Write out the current line for debugging
//                    if (currPosOnChr==refChromPosList.length-1 || linesForChr>currPosOnChr || (currPosOnChr>9584780 && currPosOnChr<9584800)) {
//                        //System.out.println(currPosOnChr+"\t"+(refChromPosList.length-1)+"\t"+totalLines);
//                        System.out.println(line);
//                    }
                    if (posCol>-1) { // process position-based data file
                        WriteValues(next); //Write out values for this line  ("next" is the line we're processing)
                        lastLinePos= currPosOnChr+1; lastStart= currPosOnChr;
                        totalLines++; linesForChr++; totalNonZeroLines++;
                    } else { // processing ranges
                        // Looking for the value from the refChromPosList array for this chromosome.
                        int end = Integer.parseInt(next[endCol]);
                        if (!oneBased()) end++;
                        end= (inclusive?end+1:end);
                        for (int i = currPosOnChr; i < end; i++) {
                            WriteValues(next); //Write out values for this line
                            totalLines++; linesForChr++; totalNonZeroLines++;
                        }
                        lastLinePos= end; //end, because even if inclusive has been transformed to exclusive
                        lastStart= currPosOnChr; //last is just for making sure that lines aren't duplicated, and checks start of range
                        lastEnd= end-1;
                    }
                } //  end while loop for processing this chrom file
                // Write everything to null between the last line and the end of the chromosome when it hits the end of the file
                // add +1 because writeNull writes until " < end", assuming 0-based, but this is 1-based.
                writeNull(lastLinePos,refChromPosSize); 
                System.out.println("File changed, " + linesForChr+"  lines output for chr "+currChr  
                        + " lastLinePos " + lastLinePos + " refChromPosSize:" + refChromPosSize);
                linesForChr = 0;
            } //  end for loop processing each file
 
            Shutdown(); // close out writers
        } catch (Exception e) {
            e.printStackTrace();
            Shutdown();
            return null;
        }
        System.out.println("Process took " + (System.nanoTime() - time)/1e9 + " seconds.");
        System.out.println("Wrote the files with totalReads: " + totalLines 
                + " totalZeroLines written: " + totalZeroLines + " totalNonZeroLines written: " + totalNonZeroLines);
        return null;
    } 
    
    
    // NOTE:  This method altered from the SNPOnly version, which used "end"
    // as "exclusive" (ie, it looped until posIdx was < end)
    private void writeNull(int start, int end) {
        try {
        for (int posIdx=start;posIdx <= end; posIdx++) {
            if (realHM!=null) {for (Integer i:realWriters.keySet()) {
                realWriters.get(i).writeFloat(missToZero.value()?0:-Float.MAX_VALUE);
            }}
            if (log10HM!=null) {for (Integer i:log10Writers.keySet()) {
                log10Writers.get(i).writeFloat(missToZero()?0:-Float.MAX_VALUE);
            }}
            if (shortHM!=null) {for (Integer i:shortWriters.keySet()) {
                shortWriters.get(i).writeShort(missToZero() ? 0:Short.MIN_VALUE);
            }}
            if (longHM!=null) {for (Integer i:longWriters.keySet()) {
                longWriters.get(i).writeLong(missToZero() ? 0 :Long.MIN_VALUE);
            }}
            // You ALWAYS need the \n, this worked in TE_LTRSuperFamily ... for MIchelle
            //String fam = sampLine.substring(tabPos[12] + 1, tabPos[13]) + "\n";
            //writerFam.writeChars(fam.getBytes()); 
            if (charHM!=null) {
                for (Integer i:charWriters.keySet()) {
                    String nullString = "NULL\n";
                  charWriters.get(i).writeChars(nullString.getBytes()); // LCJ - added NULL - see if it worksZZ
                }
            }
          
            if (intHM!=null) {for (Integer i:intWriters.keySet()) {
                intWriters.get(i).writeInt(missToZero() ? 0 : Integer.MIN_VALUE);
            }}
            if (byteHM!=null) {for (Integer i:byteWriters.keySet()) {
                byteWriters.get(i).writeByte(missToZero() ? 0:Byte.MIN_VALUE);
            }}
            if (alleleHM!=null) {for (Integer i:alleleWriters.keySet()) {
                alleleWriters.get(i).writeByte(Byte.MIN_VALUE); // for missing alleles write null
            }}
            totalLines++; linesForChr++; totalZeroLines++;
        }
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Problem writing nulls between "+start+" and "+end);
        }
    }
    
    private void WriteValues(String[] next) {
        try {
            if (realHM!=null) {for (Integer i:realWriters.keySet()) {
                double val= Double.MIN_VALUE; float storedVal;
                try {
                    val = Double.parseDouble(next[i.intValue()]);
                } catch (NumberFormatException nfe) {
                    val = Double.MIN_VALUE;
                }
                if (val==Double.MIN_VALUE) storedVal= missToZero.value()?0:-Float.MAX_VALUE;
                else if (negToZero.value() && val<=0) storedVal= 0;
                else storedVal= (float) val;
                realWriters.get(i).writeFloat(storedVal);
            }}
            if (log10HM!=null) {for (Integer i:log10Writers.keySet()) {
                double val= Double.MIN_VALUE; double storedVal;
                try {
                    val = Double.parseDouble(next[i.intValue()]);
                } catch (NumberFormatException nfe) {
                    val = Double.MIN_VALUE;
                }
                if (val==Double.MIN_VALUE) storedVal= missToZero.value()?0:-Float.MAX_VALUE;
                else if (negToZero.value() && val<=0) storedVal= 0;
                else storedVal = val==Double.MIN_VALUE?-Float.MAX_VALUE:(float)(-Math.log10(val));
                log10Writers.get(i).writeFloat((float)storedVal);
                //log10Writers.get(i).writeFloat(val==Double.MIN_VALUE?-Float.MAX_VALUE:(float)(-Math.log10(val)));
            }}
            if (shortHM!=null) {for (Integer i:shortWriters.keySet()) {
                short val= Short.MIN_VALUE; short storedVal;
                try {
                    val = Short.parseShort(next[i.intValue()]);
                } catch (NumberFormatException nfe) {
                    val = Short.MIN_VALUE;
                }
                if (val==Short.MIN_VALUE) storedVal= missToZero.value()?0:Short.MIN_VALUE;
                else if (negToZero.value() && val<=0) storedVal= 0;
                else storedVal = val;
                shortWriters.get(i).writeShort(storedVal);
                //shortWriters.get(i).writeShort(val);
            }}
            if (longHM!=null) {for (Integer i:longWriters.keySet()) {
                long val= Long.MIN_VALUE;long storedVal;
                try {
                    val = Long.parseLong(next[i.intValue()]);
                } catch (NumberFormatException nfe) {
                    val = Long.MIN_VALUE;
                }
                if (val==Long.MIN_VALUE) storedVal= missToZero.value()?0:Long.MIN_VALUE;
                else if (negToZero.value() && val<=0) storedVal= 0;
                else storedVal = val;
                longWriters.get(i).writeLong(storedVal);
                //longWriters.get(i).writeLong(val);
            }}
            // always want \n !!
            if (charHM!=null) {
                for (Integer i:charWriters.keySet()) {
                  String charData = next[i.intValue()] + "\n";
                  charWriters.get(i).writeChars(charData.getBytes());
                }
            }
            if (intHM!=null) {for (Integer i:intWriters.keySet()) {
                int val= Integer.MIN_VALUE;int storedVal;
                try {val = Integer.parseInt(next[i.intValue()]);
                } catch (NumberFormatException nfe) {val = Integer.MIN_VALUE;}
                if (val==Integer.MIN_VALUE) storedVal= missToZero.value()?0:Integer.MIN_VALUE;
                else if (negToZero.value() && val<=0) storedVal= 0;
                else storedVal = val;
                intWriters.get(i).writeInt(storedVal);
                //intWriters.get(i).writeInt(val);
            }}
            if (byteHM!=null) {
                for (Integer i:byteWriters.keySet()) {
                  byte val= Byte.MIN_VALUE;byte storedVal;                 
                  try {
                    val = (byte)(Integer.parseInt(next[i.intValue()]));                    
                  } catch (NumberFormatException nfe) {
                    val = Byte.MIN_VALUE;
                  }
                  if (val==Byte.MIN_VALUE) storedVal= missToZero.value()?0:Byte.MIN_VALUE;
                  else if (negToZero.value() && val<=0) storedVal= 0;
                  else storedVal = val;
                  byteWriters.get(i).writeByte(val);
                }
             }
            if (alleleHM!=null){
                // Not checking missToZero or negToZero - these are character alleles, missing is null
                // otherwise stored as is done in TASSEL as 0-5 for A/G/C/T/+/-
                for (Integer idx:alleleWriters.keySet()) {
                    byte val = Byte.MIN_VALUE;// monetdb stores Byte.MIN_VALUE as null
                    try {
                        String allele = next[idx.intValue()];
                        int alleleInt = ColumnsToBinaryUtils.getIntFromSeq(allele);
                        // getIntFromSeq allows for ins/del as 4/5. 
                        if (alleleInt >= 0 || alleleInt <= 5) val =(byte) alleleInt; 
                        else 
                            throw new IllegalStateException("Illegal allele char, must be A/C/G/T/+/- but found " + allele);
                        alleleWriters.get(idx).writeByte(val); // store null
                    } catch (Exception exc) {
                        exc.printStackTrace();
                        throw new IllegalStateException("Error parsing  allele column " + alleleWriters.get(idx));
                    }
                }
            }
                
        } catch (Exception e) {
            for (String i:next) {System.out.println(i);}
            e.printStackTrace();
        }
    }
    
    private void FindColumns(String filename) {
        // When broken by chromosome, column names shoudl be the
        // same for all files.  Send in any chromosome file and use
        // it to set the columns names
        String[] next= null;
        try { 
            FileInputStream fis = new FileInputStream(filename);
            Scanner scanner = new Scanner(fis); boolean first= true;
            while (scanner.hasNextLine()) {
                if (!first) break;
                String line = scanner.nextLine();
                if (line.isEmpty() || line.startsWith("#")) {continue;}
                next = line.split("\t");
                if (first) {
                    first= false;
                }
            }
            scanner.close();
            fis.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        //get index of chromosome
        int search= -1;
        for (int col = 0; col < next.length; col++) {
            if (next[col].toLowerCase().startsWith("chr")) {search= col; break;}
        }
        if (search<0) {throw new IllegalStateException("Couldn't find chromosome in first non-empty line of input file");
        } else {
            chrCol= search;
            //chrName= next[search]; // now defaulting this to "CHR"
            System.out.println("Chr found in column "+search);
        }
        //Get position
        search= -1;
        if (range.value()) { //if range, looks for first last (inclusive) or start end (exclusive)
            for (int col = 0; col < next.length; col++) {
                if (next[col].equalsIgnoreCase("start")) {
                    search= col; break;
                }
            }
            if (search<0) {
                inclusive= true;
                for (int col = 0; col < next.length; col++) {
                    if (next[col].equalsIgnoreCase("first")) {
                        search= col; break;
                    }
                }
                if (search<0) throw new IllegalStateException("Couldn't find first or start in first line of input file to indicate range");
                else startCol= search;
            } else startCol= search;
            search= -1;
            // look for "last" or "end" for end of range
            for (int col = 0; col < next.length; col++) {
                if (next[col].equalsIgnoreCase(inclusive?"last":"end")) {
                    search= col; break;
                }
            }
            if (search<0) {
                if (search<0) throw new IllegalStateException("Couldn't find "+(inclusive?"last":"end")+" in first line of input file to indicate range");
            } else endCol= search;
            System.out.println("Using columns "+startCol+" and "+endCol+" to indicate ends of ranges, "+(inclusive?"inclusive":"exclusive"));
        } else {
            for (int col = 0; col < next.length; col++) {
                if (next[col].toLowerCase().startsWith("pos")) {
                    search= col; 
                    break;
                }
            }
            if (search<0) {System.out.println("Cannot find column for pos");
            } else {posCol= search;System.out.println("Pos found in column "+search);}
            if (chrCol<0 || posCol<0) throw new IllegalStateException("Couldn't find chromosome or position in first line of input file");
        }
        int ncols= 0;
        //Get indices of column names to keep and transform with -log10, if not null. Don't use binary search, because not sorted and won't work
        if (colsLog10.value()!=null && !colsLog10.value().equalsIgnoreCase("null") && !colsLog10.value().isEmpty()) {
            log10HM= new HashMap<>(); search= -1;
            for (String s:colsLog10.value().split(",")) {
                search= -1;
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {search= col; break;}
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {
                    log10HM.put(search, s);
                    log10Writers= new HashMap<>();
                    ncols++;
                    System.out.println("Found "+s+" in column "+search+" as -log10 transformed float");
                }
            }
        }
        //Get indices of column names to keep as float
        if (colsFloat.value()!=null && !colsFloat.value().equalsIgnoreCase("null") &&!colsFloat.value().isEmpty()) {
        realHM= new HashMap<>(); realWriters= new HashMap<>(); search= -1;
            for (String s:colsFloat.value().split(",")) {
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {search= col; break;}
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {realHM.put(search, s);ncols++;System.out.println("Found "+s+" in column "+search+" as real (float)");}
            }
        }
        //Get indices of column names to keep as int
        if (colsInt.value()!=null && !colsInt.value().equalsIgnoreCase("null") && !colsInt.value().isEmpty()) {
        intHM= new HashMap<>(); intWriters= new HashMap<>(); search= -1;
            for (String s:colsInt.value().split(",")) {
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {search= col; break;}
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {intHM.put(search, s);ncols++;System.out.println("Found "+s+" in column "+search+" as int");}
            }
        }
        //Get indices of column names to keep as short
        if (colsShort.value()!=null && !colsShort.value().equalsIgnoreCase("null") && !colsShort.value().isEmpty()) {
        shortHM= new HashMap<>(); shortWriters= new HashMap<>(); search= -1;
            for (String s:colsShort.value().split(",")) {
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {search= col; break;}
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {shortHM.put(search, s);ncols++;System.out.println("Found "+s+" in column "+search+" as short");}
            }
        }
        //Get indices of column names to keep as long
        if (colsLong.value()!=null && !colsLong.value().equalsIgnoreCase("null") && !colsLong.value().isEmpty()) {
        longHM= new HashMap<>(); longWriters= new HashMap<>(); search= -1;
            for (String s:colsLong.value().split(",")) {
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {search= col; break;}
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {longHM.put(search, s);ncols++;System.out.println("Found "+s+" in column "+search+" as long");}
            }
        }
        //Get indices of column names to keep as char
        if (colsChar.value()!=null && !colsChar.value().equalsIgnoreCase("null") && !colsChar.value().isEmpty()) {
        charHM= new HashMap<>(); charWriters= new HashMap<>(); search= -1;
            for (String s:colsChar.value().split(",")) {
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {
                        search= col; break;
                    }
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {charHM.put(search, s);ncols++;System.out.println("Found "+s+" in column "+search+" as char");}
            }
        }  
        //Get indices of column names to keep as char
        if (colsByte.value()!=null && !colsByte.value().equalsIgnoreCase("null") && !colsByte.value().isEmpty()) {
        byteHM= new HashMap<>(); byteWriters= new HashMap<>(); search= -1;
            for (String s:colsByte.value().split(",")) {
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {
                        search= col; break;
                    }
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {byteHM.put(search, s);ncols++;System.out.println("Found "+s+" in column "+search+" as byte");}
            }
        } 
        //Get indices of column names to keep as alleles translated to 0-5
        if (colsAllele.value()!=null && !colsAllele.value().equalsIgnoreCase("null") && !colsAllele.value().isEmpty()) {
        alleleHM= new HashMap<>(); alleleWriters= new HashMap<>(); search= -1;
            for (String s:colsAllele.value().split(",")) {
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {
                        search= col; break;
                    }
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {alleleHM.put(search, s);ncols++;System.out.println("Found "+s+" in column "+search+" as allele");}
            }
        }  
         if (ncols<1) throw new IllegalStateException("No valid columns to read in!");   
    }
    
    private void Initialize() {
        // This method creates the .sql files for creating the table and copying into the table.
        // These values can be appended to the larger table.  Table name must be filled in by
        // user, but the skeleton of the sql command is there to feed into monetdb.
        totalLines = 0;
        totalZeroLines = 0;
        totalNonZeroLines = 0;
        linesForChr = 0;
        try {
            System.out.println("\nColumnsToBInaryFullGenomeTable: outBase is : " + outBase() + "\n\n");
            DataOutputStream outCopy = Utils.getDataOutputStream(outBase.value()+"copy.sql", 1040);
            DataOutputStream outCreate = Utils.getDataOutputStream(outBase.value()+"create.sql", 1040);
            outCopy.writeBytes("COPY binary into ... from ('snpChrFile.bin','snpPosFile.bin',"); boolean first= true;
            outCreate.writeBytes("CREATE TABLE ... (chr int,pos int,"); 
            String lastChar = outBase().substring(outBase().length()-1,outBase().length());
            String out = null;
            if (!(lastChar.equals("/"))) {
                out = new File(outBase.value()).getName();
            }
            //String outDir= new File(outBase.value()).getParent();
            if (realHM!=null) {
                for (Integer i:realHM.keySet()) {
                    String copyFile;
                    String createFile;
                    if (out != null) {
                        copyFile = out + realHM.get(i) + "_real.bin";
                        createFile = out +  realHM.get(i) + " real";
                    } else {
                        copyFile = realHM.get(i) + "_real.bin";
                        createFile = realHM.get(i) + " real";
                    }
                    String outFile = outBase() + realHM.get(i) + "_real.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+copyFile+"'"); first= false;
                    outCreate.writeBytes( createFile );
                   // outCreate.writeBytes(out + realHM.get(i)+" real"); 
                    realWriters.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile))));
                    //realWriters.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outDir+"/"+outFile))));
                }
            }
            if (intHM!=null) {
                for (Integer i:intHM.keySet()) {
                    String copyFile;
                    String createFile;
                    if (out != null) {
                        copyFile = out +  intHM.get(i) + "_int.bin";
                        createFile = out +  intHM.get(i) + " int";
                    } else {
                        copyFile = intHM.get(i) + "_int.bin";
                        createFile = intHM.get(i) + " int";
                    }
                    
                    String outFile = outBase() +  intHM.get(i) + "_int.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+copyFile+"'"); 
                    outCreate.writeBytes(createFile );
                    intWriters.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile))));
                }
            }
            if (byteHM!=null) { // tinyint is a byte
                for (Integer i:byteHM.keySet()) {
                    String copyFile;
                    String createFile;
                    if (out != null) {
                        copyFile = out + byteHM.get(i) + "_byte.bin";
                        createFile = out + byteHM.get(i) + " byte";
                    } else {
                        copyFile = byteHM.get(i) + "_byte.bin";
                        createFile = byteHM.get(i) + " byte";
                    }
                    String outFile = outBase() + byteHM.get(i) + "_byte.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+copyFile+"'"); 
                    outCreate.writeBytes(createFile); // monetdb doesn't have "byte", tinyint instead
                    byteWriters.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile))));
                }
            }
            if (alleleHM!=null) { // alleles are stored as bytes, here as tinyints
                for (Integer i:alleleHM.keySet()) {
                    String copyFile;
                    String createFile;
                    if (out != null) {
                        copyFile = out + alleleHM.get(i) + "_allelebyte.bin";
                        createFile = out + alleleHM.get(i) + " tinyint";
                    } else {
                        copyFile = alleleHM.get(i) + "_allelebyte.bin";
                        createFile = alleleHM.get(i) + " tinyint";
                    }
                    String outFile = outBase() + alleleHM.get(i) + "_allelebyte.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+copyFile+"'"); 
                    outCreate.writeBytes(createFile); // monetdb doesn't have "byte", tinyint instead
                    alleleWriters.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile))));
                }
            }
            if (log10HM!=null) {
                for (Integer i:log10HM.keySet()) {
                    String copyFile;
                    String createFile;
                    if (out != null) {
                        copyFile = out +  log10HM.get(i) + "_neglog10_real.bin";
                        createFile = out +  log10HM.get(i) + "_neglog10 real";
                    } else {
                        copyFile = log10HM.get(i) + "_neglog10_real.bin";
                        createFile = log10HM.get(i) + "_neglog10 real";
                    }
                    String outFile = outBase() +  log10HM.get(i) + "_neglog10_real.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+copyFile+"'");
                    outCreate.writeBytes(createFile); 
                    log10Writers.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile))));
                }
            }
            if (shortHM!=null) {
                for (Integer i:shortHM.keySet()) {
                    String copyFile;
                    String createFile;
                    if (out != null) {
                        copyFile = out + shortHM.get(i) + "_short.bin";
                        createFile = out + shortHM.get(i) + " short";
                    } else {
                        copyFile = shortHM.get(i) + "_short.bin";
                        createFile = shortHM.get(i) + " short";
                    }
                    String outFile = outBase() + shortHM.get(i) + "_short.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+copyFile+"'");
                    outCreate.writeBytes(createFile); 
                    shortWriters.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile))));
                }
            }
            if (longHM!=null) {
                for (Integer i:longHM.keySet()) {
                    String copyFile;
                    String createFile;
                    if (out != null) {
                        copyFile = out + "_"+longHM.get(i) + "_long.bin";
                        createFile = out + "_"+longHM.get(i) + " long";
                    } else {
                        copyFile = longHM.get(i) + "_long.bin";
                        createFile = longHM.get(i) + " long";
                    }
                    String outFile = outBase() + longHM.get(i) + "_long.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+copyFile+"'");
                    outCreate.writeBytes(createFile); 
                    longWriters.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile))));
                }
            }
            if (charHM!=null) {
                for (Integer i:charHM.keySet()) {
                    String copyFile;
                    String createFile;
                    if (out != null) {
                        copyFile = out + charHM.get(i) + "_char.bin";
                        createFile = out + charHM.get(i) + " char";
                    } else {
                        copyFile = charHM.get(i) + "_char.bin";
                        createFile = charHM.get(i) + " char";
                    }
                    String outFile = outBase() + charHM.get(i) + "_char.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+copyFile+"'");
                    outCreate.writeBytes(createFile); 
                    charWriters.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile))));
                }
            }
        outCopy.writeBytes(");"); outCopy.close();
        outCreate.writeBytes(");"); outCreate.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
    
    private void Shutdown() {
        try {
            scanner.close();
            fis.close();
            if (realHM!=null) {for (Integer i:realWriters.keySet()) {realWriters.get(i).close();}}
            if (intHM!=null) {for (Integer i:intWriters.keySet()) {intWriters.get(i).close();}}
            if (shortHM!=null) {for (Integer i:shortWriters.keySet()) {shortWriters.get(i).close();}}
            if (longHM!=null) {for (Integer i:longWriters.keySet()) {longWriters.get(i).close();}}
            if (charHM!=null) {for (Integer i:charWriters.keySet()) {charWriters.get(i).close();}}
            if (log10HM!=null) {for (Integer i:log10Writers.keySet()) {log10Writers.get(i).close();}}
            if (byteHM!=null) {for (Integer i:byteWriters.keySet()) {byteWriters.get(i).close();}}

        } catch (Exception e) {
            System.out.println("Problem with shutdown");
            e.printStackTrace();
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
    
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
     public static void main(String[] args) {
         GeneratePluginCode.generate(ColumnsToBinaryFullGenomeTablePlugin.class);
     }
     /**
      * Convenience method to run plugin with one return object.
      */
     // TODO: Replace <Type> with specific type.
//     public <Type> runPlugin(DataSet input) {
//         return (<Type>) performFunction(input).getData(0).getData();
//     }

     /**
      * File or a Directory containing Tab-delimited files with data
      * to add to the database. Files must be named chr01.txt,
      * chr02.txt etc!
      *
      * @return Input File 
      */
     public String inputFile() {
         return inputFile.value();
     }

     /**
      * Set Input Directory. Directory containing Tab-delimited
      * files with data to add to the database. Files must
      * be named chr01.txt, chr02.txt etc!
      *
      * @param value Input Directory or File
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin inputFile(String value) {
         inputFile = new PluginParameter<>(inputFile, value);
         return this;
     }

     /**
      * Link to maize reference file.
      *
      * @return refFile
      */
     public String refFile() {
         return refFile.value();
     }

     /**
      * Set refFile.
      *
      * @param value refFile
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin refFile(String value) {
         refFile = new PluginParameter<>(refFile, value);
         return this;
     }

     /**
      * Output directory and base filename to hold the binary
      * files. Will make directory if neccesary
      *
      * @return outBase
      */
     public String outBase() {
         return outBase.value();
     }

     /**
      * Set outBase. Output directory and base filename to
      * hold the binary files. Will make directory if neccesary
      *
      * @param value outBase
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin outBase(String value) {
         outBase = new PluginParameter<>(outBase, value);
         return this;
     }

     /**
      * Comma separated list of column names to generate real
      * binaries for
      *
      * @return Columns keep as Real (Float)
      */
     public String colsFloat() {
         return colsFloat.value();
     }

     /**
      * Set Columns keep as Real (Float). Comma separated list
      * of column names to generate real binaries for
      *
      * @param value Columns keep as Real (Float)
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin colsFloat(String value) {
         colsFloat = new PluginParameter<>(colsFloat, value);
         return this;
     }

     /**
      * Comma separated list of column names to generate int
      *  binaries for
      *
      * @return Columns keep as Int
      */
     public String colsInt() {
         return colsInt.value();
     }

     /**
      * Set Columns keep as Int. Comma separated list of column
      * names to generate int  binaries for
      *
      * @param value Columns keep as Int
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin colsInt(String value) {
         colsInt = new PluginParameter<>(colsInt, value);
         return this;
     }

     /**
      * Comma separated list of column names to generate short
      * binaries for
      *
      * @return Columns keep as Short
      */
     public String colsShort() {
         return colsShort.value();
     }

     /**
      * Set Columns keep as Short. Comma separated list of
      * column names to generate short binaries for
      *
      * @param value Columns keep as Short
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin colsShort(String value) {
         colsShort = new PluginParameter<>(colsShort, value);
         return this;
     }

     /**
      * Comma separated list of column names to generate long
      * binaries for
      *
      * @return Columns keep as Long
      */
     public String colsLong() {
         return colsLong.value();
     }

     /**
      * Set Columns keep as Long. Comma separated list of column
      * names to generate long binaries for
      *
      * @param value Columns keep as Long
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin colsLong(String value) {
         colsLong = new PluginParameter<>(colsLong, value);
         return this;
     }

     /**
      * Comma separated list of column names to generate byte
      * binaries for
      *
      * @return Columns keep as Byte
      */
     public String colsByte() {
         return colsByte.value();
     }

     /**
      * Set Columns keep as Byte. Comma separated list of column
      * names to generate byte binaries for
      *
      * @param value Columns keep as bhte
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin colByte(String value) {
         colsByte = new PluginParameter<>(colsByte, value);
         return this;
     }
     /**
      * Comma separated list of column names to generate byte
      * binaries for
      *
      * @return Columns keep as Byte
      */
     public String colsAllele() {
         return colsAllele.value();
     }

     /**
      * Set Columns translate allele values to 0-5. Comma separated list of column
      * names of single character alleles to be translated from A/C/G/T/+/- to
      * 0-5.
      *
      * @param value Columns keep for translated alleles stored as bytes
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin colsAllele(String value) {
         colsAllele = new PluginParameter<>(colsAllele, value);
         return this;
     }     
     /**
      * Comma separated list of column names to generate char
      * binaries for
      *
      * @return Columns keep as Char
      */
     public String colsChar() {
         return colsChar.value();
     }

     /**
      * Set Columns keep as Char. Comma separated list of column
      * names to generate char binaries for
      *
      * @param value Columns keep as Char
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin colsChar(String value) {
         colsChar = new PluginParameter<>(colsChar, value);
         return this;
     }

     /**
      * Comma separated list of column names to first transform
      * using -log10 then generate binaries for
      *
      * @return Columns to Keep and transform -log10
      */
     public String colsLog10() {
         return colsLog10.value();
     }

     /**
      * Set Columns to Keep and transform -log10. Comma separated
      * list of column names to first transform using -log10
      * then generate binaries for
      *
      * @param value Columns to Keep and transform -log10
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin colsLog10(String value) {
         colsLog10 = new PluginParameter<>(colsLog10, value);
         return this;
     }

     /**
      * Columns for range data. If true, will look for 'start'
      * and 'end' (inclusive exclusive) or 'first' 'last' (inclusive
      * inclusive) instead of 'Pos' 
      *
      * @return Range information?
      */
     public Boolean range() {
         return range.value();
     }

     /**
      * Set Range information?. Columns for range data. If
      * true, will look for 'start' and 'end' (inclusive exclusive)
      * or 'first' 'last' (inclusive inclusive) instead of
      * 'Pos' 
      *
      * @param value Range information?
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin range(Boolean value) {
         range = new PluginParameter<>(range, value);
         return this;
     }

     /**
      * Will set negative float denoted column values to zero
      * (otherwise they are set to Float.MIN
      *
      * @return Negative floats to zero?
      */
     public Boolean negToZero() {
         return negToZero.value();
     }

     /**
      * Set Negative floats to zero?. Will set negative float
      * denoted column values to zero (otherwise they are set
      * to Float.MIN
      *
      * @param value Negative floats to zero?
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin negToZero(Boolean value) {
         negToZero = new PluginParameter<>(negToZero, value);
         return this;
     }

     /**
      * Will set missing float denoted column values to zero
      * (otherwise they are set to Float.MIN
      *
      * @return Missing float values to zero?
      */
     public Boolean missToZero() {
         return missToZero.value();
     }

     /**
      * Set Missing float values to zero?. Will set missing
      * float denoted column values to zero (otherwise they
      * are set to Float.MIN
      *
      * @param value Missing float values to zero?
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin missToZero(Boolean value) {
         missToZero = new PluginParameter<>(missToZero, value);
         return this;
     }
     /**
      * Will assume positions and ranges are 1-based
      * unless otherwise set.
      * 
      * @return positions are 1-based?
      */
     public Boolean oneBased() {
         return oneBased.value();
     }

     /**
      * Set positions to 1-based unless this variable
      * is "false".  
      *
      * @param value 1-based positions?
      *
      * @return this plugin
      */
     public ColumnsToBinaryFullGenomeTablePlugin oneBased(Boolean value) {
         oneBased = new PluginParameter<>(oneBased, value);
         return this;
     }

}
