/**
 * 
 */
package net.maizegenetics.util;

import org.junit.Test;

/**
 * This test calls ConvertOldFastqToModernFormatPlugin to
 * transform fastq files from an old illumina format which
 * contained multiple taxa with barcode added to each read sequence
 * into the newer format consisting of 1 fastq file per taxon
 * and no barcode in the read sequences.
 * 
 * @author lcj34
 *
 */
public class ConvertOldFastqToModernFormatTest {

    
    @Test
    public void testConvertFiles() {
        System.out.println("Testing ConvertOldFastqToMOdernFormat ...");
        
        
        String inputFile = "/Volumes/Samsung_T1/repgen/gbs_junit_fastq";
        String keyFile = "/Users/lcj34/git/tassel-5-test/dataFiles/GBS/Pipeline_Testing_key.txt";
        
        String outputDir = "/Volumes/Samsung_T1/repgen/convertGBS_testjunit";
        String projectName = "chrom9-10junit";
        String enzyme = "ApeKI";
        new ConvertOldFastqToModernFormatPlugin()
        .inputFile(inputFile)
        .keyFile(keyFile)
        .outputDir(outputDir)
        .projectName(projectName)
        .enzyme(enzyme)
        .performFunction(null);
    }

}
