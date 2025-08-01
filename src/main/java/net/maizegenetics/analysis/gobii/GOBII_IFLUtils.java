
package net.maizegenetics.analysis.gobii;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.StringTokenizer;

import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.Utils;

/**
 * This class contains utility methods for pulling values out of
 * hmp or vcf files needed when creating the intermediate files
 * for loading into GOBII postgres and monetdb instances
 * 
 * @author lcj34
 *
 */
public class GOBII_IFLUtils {
    // Chromosome is in first column for VCF file, in 3rd column for hmp file
    // first tab occurs after the first column (tabPos[0])
    public static int getChromFromLine(String mline, boolean isVCF, int[] tabPos) {
        int chrom;
        if (isVCF) {
            chrom = Integer.parseInt(mline.substring(0,tabPos[0]));
        } else {
            // tabPos[1] + 1 is beginning of the 3rd column
            chrom = Integer.parseInt(mline.substring(tabPos[1]+1,tabPos[2]));
        }
        return chrom;
    }
    
    public static String getMarkerNameFromLine(String mline, boolean isVCF, int[] tabPos, String mapsetname){
        String name = null;
        if (isVCF) {
            name = mline.substring(tabPos[1]+1,tabPos[2]);
            if (name.equals(".")) { // "." is used in vcf for unknown. 
                // Marker name becomes S<chr>_<pos>, e.g. S10_20 for position 20 on chrom 10
                // Marker name becomes PZ.V.chrom.position, e.g. PZ.2.9.123132
                String pos = mline.substring(tabPos[0]+1,tabPos[1]);
                String chrom = mline.substring(0,tabPos[0]);
               // name = "S" + chrom + "_" + pos;
                String mapset = "";
                if (mapsetname.toUpperCase().equals("AGPV2")){
                    mapset = "2";
                } else if (mapsetname.toUpperCase().equals("AGPV3")) {
                    mapset = "3";
                } else if (mapsetname.toUpperCase().equals("AGPV4")) {
                    mapset = "4";
                } else {
                    System.out.println("WARNING: getMarkerNameFromLine - bad mapset name: " + mapsetname);
                }
                
                name = "PZ." + mapset + "." + chrom + "." + pos;
                //System.out.println("LCJ GOBII_IFLUtils.getMarkerNameFromLine- creatwed marker name : " + name);
            }
            // There COULD have multiple identifiers.  If so, they are separated
            // by colons with no white space.  Will this be encountered in our files??

        } else {
            name = mline.substring(0, tabPos[0]); // store rs# as name
            if (name == null) {
                System.out.println("WARNING: getMarkerNameFromLine: hmp name rs field is NULL");
            }
        }
        return name;
    }
    
    // position is in second column for VCF file, in 4th column for hmp file
    // first tab occurs after the first column (tabPos[0])
    public static int getPosFromLine(String mline, boolean isVCF, int[] tabPos) {
        int pos;
        if (isVCF) {
            pos = Integer.parseInt(mline.substring(tabPos[0]+1,tabPos[1]));
        } else {
            // tabPos[0] + 1 is beginning of the 2nd column
            pos = Integer.parseInt(mline.substring(tabPos[2]+1,tabPos[3]));
        }
        return pos;
    }
    
    // find the strand - column 5 in hmp.txt file.  Hapmap does not have
    // a strand field.
    public static String getStrandFromLine(String mline, boolean isVCF, int[] tabPos) {
        String strand = null;
        // for hmp, strand is column 5
        if (isVCF) {
            strand = "Unknown";
        } else {
            strand = mline.substring(tabPos[3]+1,tabPos[4]);
            if (strand.equals("+")) {
                strand = "Forward";
            } else if (strand.equals("-")) {
                strand = "Reverse";
            } else {
                strand = "Unknown";
            }
        }
        return strand;
    }

    public static String addMonetdbVariantData(String ref, String altsOrig, String mline, boolean isVCF, int[]tabPos){
           // boolean isVCF, int[]tabPos,String platformName){ // add this when do illumina
        StringBuilder variantsSB = new StringBuilder();
        if (isVCF) {
            // The VCF shows the allele call for each taxa as x/y
            // from the vcf format definition:
            // The allele values are 0 for the reference allele (what is in the REF field),
            // 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on.
            
            // strip the "{" off of alts
            String alts = altsOrig.substring(1,altsOrig.length()-1);
            // Create alleles array
            String[] altsTokens = alts.split("/"); // need to split to get rid of the /'s !!!
            char[] alleles = new char[altsTokens.length + 1]; // length of alleles string plus 1 for ref
            alleles[0] = ref.charAt(0); // should just be 1
            for (int idx = 1,altsIdx=0; idx < alleles.length; idx++,altsIdx++) {
                //alleles[idx]=alts.charAt(altsIdx); //
                alleles[idx] = altsTokens[altsIdx].charAt(0); // 
            }
            
            String taxaString = mline.substring(tabPos[8]+1); // skip non-taxa headers
            StringTokenizer taxaValues = new StringTokenizer(taxaString);
                       
            //System.out.println("LCJ - addMonetdbVariantData, ref: " + ref + ", altsOrig: " + altsOrig + " alts after modify: " + alts);
            boolean firstTaxa = true;
            while (taxaValues.hasMoreTokens()){
                if (!firstTaxa) {
                    variantsSB.append("\t"); // only append tab if we know there is a next taxa
                } else {
                    firstTaxa = false;
                }
                String nextTaxa = taxaValues.nextToken();
                //System.out.println("LCJ - nextTaxa from taxaVlues tokens is: " + nextTaxa);
 
                if (nextTaxa.equals(".") || nextTaxa.equals("./.")) {
                    byte unknown = (byte)0xFF;
                    variantsSB.append(NucleotideAlignmentConstants.getNucleotideIUPAC(unknown));
                    
                } else {
                    int end = nextTaxa.indexOf(":");
                    if (end != -1) {
                        // nothing beyond the x/x, ie there is no :
                        // Robert has removed this data in some of the files.
                        nextTaxa = nextTaxa.substring(0,end); // grab just the first part, e.g. 0/0 or 0/1 or ./. etc
                    }
                    
                    // get index into list alleles array
                    int a1 = nextTaxa.charAt(0) - '0';
                    int a2 = nextTaxa.charAt(2) - '0';
                    
                    if (a1 < 0 || a2 < 0) {
                        // THis is same as above, when there is . or ./.
                        // It is seen in maize chrom 10 where ref is N and allt is ./1 (vs ./. above)
                        // this makes a1 be negative.  MarkerDNARun_IFLFilePlugin only allows for 1 ref,
                        // and that is the ref from the genome fasta file.  I notice in vcf format 4.3, you
                        // can have multiple values for the ref, e.g. Ref=GTC, ALTs = G,GTCT.   We many
                        // need to handle this.  FOr now, this is catching the "." in the ref spot.
                        // SHould be ok to continue processing with these notices.:w!
                        
                        System.out.println("LCJ - OOPS negative values!! set to N  for nextTaxa " + nextTaxa);
                        System.out.println("  mline: " + mline.substring(0, tabPos[8]));
                        variantsSB.append(NucleotideAlignmentConstants.getNucleotideIUPAC((byte)(0xFF)));
                    } else {
                        // get bytes, set IUPAC values
                        byte highbyte = NucleotideAlignmentConstants.getNucleotideAlleleByte(alleles[a1]);
                        byte lowbyte = NucleotideAlignmentConstants.getNucleotideAlleleByte(alleles[a2]);
                        String iupacByte = NucleotideAlignmentConstants.getNucleotideIUPAC((byte)((highbyte << 4) | lowbyte));                       
                        variantsSB.append(iupacByte);   
                    }
                }                   
            }
            variantsSB.append("\n");
            //System.out.println("LCJ - variants line: " + variantsSB.toString());
                    
        } else {
            // hmp - just need to grab the values
            String taxaString = mline.substring(tabPos[10]+1); // skip non-taxa headers
            // This gets all the remaining values on the line
            StringTokenizer taxaValues = new StringTokenizer(taxaString);
            boolean firstTaxa = true;
            while (taxaValues.hasMoreTokens()){
                if (!firstTaxa) {
                    variantsSB.append("\t"); // only append tab if we know there is a next taxa
                } else {
                    firstTaxa = false;
                }
                // The illumina code below needs testing!
//                if (platformName.toUpperCase().equals("ILLUMINA")) {
//                    // Illumina 50K has double letter values, need to be translated
//                    // to IUPAC values. 
//                    String nextTaxa = taxaValues.nextToken();
//                    byte highbyte = NucleotideAlignmentConstants.getNucleotideAlleleByte(nextTaxa.charAt(0));
//                    byte lowbyte = NucleotideAlignmentConstants.getNucleotideAlleleByte(nextTaxa.charAt(1));
//                    String iupacByte = NucleotideAlignmentConstants.getNucleotideIUPAC((byte)((highbyte << 4) | lowbyte));
//                    variantsSB.append(iupacByte);  
//                } else {
//                    variantsSB.append(taxaValues.nextToken());
//                }
                // COMMENT LINE BELOW when add support for Illumina
                variantsSB.append(taxaValues.nextToken());
            }
            variantsSB.append("\n");
        }
        return variantsSB.toString();
    }
    public static String getAltsForRef(String ref) {
        StringBuilder alts = new StringBuilder();
        String[] aTokens = {"A","C","G","T"};
        boolean first = true;
        alts.append("{");
        for (String allele: aTokens) {
            if (!(allele.equals(ref))){
                if (!first) {
                    alts.append("/");                               
                }
                alts.append(allele);
                first = false;                           
            }
        }
        alts.append("}");
        return alts.toString();
    }
    // get alts from hapmap or vcf file.  This method loses some of its usefulness.
    // We now initially set the alts to whichever of A/C/G/T that is not the reference.
    // Eventually we need to handle deletions.  Will leave this code in place now as
    // the actual Alts from the vcf/hmp file are recorded to a separate file for 
    // use in updating the db entries.  That may be done when we handle deletions.
    public static String getAltsFromLine(String mline, String ref, boolean isVCF, int[] tabPos) {
        StringBuilder alts = new StringBuilder();
        if (isVCF) {
            // look at ref and alts column.  Robert has probably identified th ref
            // correctly, but if the vcf was created by TASSEL, this is not necessarily so.
            String vcfRef = mline.substring(tabPos[2]+1, tabPos[3]);
            String alleles = mline.substring(tabPos[3]+1, tabPos[4]);
            String[] aTokens = alleles.split(",");
                        
            alts.append("{");
            boolean first = true;
            for (String allele: aTokens) {
                if (!(allele.equals(ref))){
                    if (!first) {
                        alts.append("/");                               
                    }
                    // VCF can have values inside angle brackets,
                    // e.g. <INS> or <DEL>.  Must look for this pattern
                    if (allele.startsWith("<")) {
                        // Strip off the <>
                        String val = allele.substring(1,allele.length()-1);
                        if (val.equals("INS")) {
                            alts.append("+");
                        } else if (val.equals("DEL")) {
                            alts.append("-");
                        } // other values within <> ???
                    } else {
                        alts.append(allele);
                    }                   
                    first = false;                           
                }
            }
            if (!(vcfRef.equals(ref))) {
                // The passed in ref allele is based on a fasta reference file.
                // The VCF file also declares a ref.  We use the passed in ref as 
                // the ref.  If the VCF declared ref doesn't match this, then add
                // it as an alt.
                alts.append(vcfRef);
            }
            alts.append("}"); 
        } else {
            // get alts.  Pull the alleles from column 2 of the hmp.txt file
            // Alleles are separated by "/".  Any allele that does not match
            // the reference allele gets stored in the alts array.
            String alleles = mline.substring(tabPos[0]+1,tabPos[1]);
            String[] aTokens = alleles.split("/");
                           
            alts.append("{");
            boolean first = true;
            for (String allele: aTokens) {
                if (!(allele.equals(ref))){
                    if (!first) {
                        alts.append("/");                               
                    }
                    alts.append(allele);
                    first = false;                           
                }
            }
            alts.append("}");
        }
        return alts.toString();
    }
    
    /**
     * This method is created to break up very large files that GOBII can't handle.
     * For example:  one dataset has 83M+ lines.  Takes GOBII overnight just to load
     * the marker_linkage_group or dataset_marker data.  Yaw suggests we break it into
     * files of size 10M.  So pass in the file, give an output directory, and split
     * the file into smaller files.  Each smaller file must retain the header, and it
     * must end with the same table name (e.g. DS_4.marker_linkage_gropu must still 
     * 
     * @param infile
     * @param outdir
     * @param maxSize
     */
    public static void splitIFLFile(String infile, String outdir, int maxSize) {
        BufferedReader br = Utils.getBufferedReader(infile, 1 << 22);
        //DataOutputStream writerMarker = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile)));
        DataOutputStream bw = null;
        StringBuilder sb = new StringBuilder();
        int fileCount = 1;
        Path filePath=Paths.get(infile);
        String fileName = filePath.getFileName().toString();
        //String fileNameSubstr = fileName.substring(0,fileName.indexOf("."));
 
        
        String[] fileparts = fileName.split("\\."); // should be dataset_name.<tableName> - only 1 period in file name !!
        System.out.println("LCJ - fileName is " + fileName + ", fileparts.length " + fileparts.length + "\n");
        if (fileparts.length != 2) {
            System.out.println("LCJ - input file must have 1 and only 1 period in the name to differentiate dataset name from table name");
            System.out.println("Inputfile example:  DS_4.marker");
        }
        try {
            String line = br.readLine(); // should be the header
            String headerLine = line + "\n";
            int linecount = 0;
            int totalLines = 0;
            int maxBuffer = 10000; // max lines before we write out the buffer.
            String outfile = outdir + fileparts[0] + "_" + fileCount + "." + fileparts[1];
            // add header to file
            bw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile)));
            bw.writeBytes(headerLine); // already added \n to header line
            fileCount++;
            while ((line = br.readLine()) != null){
                sb.append(line);
                sb.append("\n"); // readline removes the carriage return
                linecount++;
                totalLines++;
                if (totalLines == maxSize || linecount == maxBuffer) {
                    // write to file, clear out buffer
                    bw.writeBytes(sb.toString());
                    sb.setLength(0); // zero out the buffer                   
                    if (totalLines == maxSize) {
                        System.out.println("LCJ - linecout " + linecount + " is max size, close out file " + outfile);
                        // close the writer, get new one for next set
                        bw.close();
                        System.out.println("LCJ - wrote file " + outfile + " with lines " + linecount);
                        outfile = outdir + fileparts[0] + "_" + fileCount + "." + fileparts[1]; // filecount is different
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
            }
            bw.close();  // close out buffer, return
                    
        } catch (Exception exc) {
            System.out.println("LCJ - whoops - error reading/writing file " + infile + " or outfile " + fileCount);
        }
    }
}
