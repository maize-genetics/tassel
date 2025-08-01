/**
 * 
 */
package net.maizegenetics.analysis.monetdb;

import org.junit.Test;

/**
 * @author lcj34
 *
 */
public class ColumnsToBinarySNPOnlyTablePluginTest {
    @Test
    public void testColumnsToBinarySNPOnlyTest() throws Exception {
       // String inputFile = "/Volumes/Samsung_T1/machineLearningFiles/kelly_dbTestCode/test.FirstLast.txt";
        // In test file Ames_DTA.txt, the start/end are inclusive/exclusive, which is Kelly's default for start/end 
       
       // String inputFile = "/Users/lcj34/notes_files/machineLearningDB/monetdb_dir/kelly_scripts/Ames_DTA.txt";
        //String inputFile = "/Users/lcj34/notes_files/machineLearningDB/monetdb_dir/kelly_scripts/Ames_DTA_byChrom";
        String inputFile = "/Volumes/Samsung_T1/machineLearningFiles/kelly_snpOnly/fstMappingPops/MAFAndCov.popInfo.txt";
        //String refDir = "/Volumes/Samsung_T1/machineLearningFiles/kelly_dbTestCode";
        String refDir = "/Volumes/Samsung_T1/machineLearningFiles/refGenomeFiles";
        //String outBase = "/Volumes/Samsung_T1/machineLearningFiles/kelly_dbTestCode/lcjTestResults/kelly_junit";
        //
        String outBase = "/Users/lcj34/notes_files/machineLearningDB/monetdb_dir/kelly_scripts/kelly_output/maf_andCov_output/";
        String floatColNames = "medianR,medianS,stdR,stdS,meansS,wilPS30,wilP30Bon";
        //String intColNames = "N";
        String intColNames = "Flint_NMaj,Flint_NNonMiss,Ames_NMaj,Ames_NNonMiss,Dent_NMaj,Dent_NNonMiss,TurkeyPen_NMaj,TurkeyPen_NNonMiss,CNNAM_NMaj,CNNAM_NNonMiss,NAM_NMaj,NAM_NNonMiss";
        String alleleColNames = "Flint_MajAll,Ames_MajAll,Dent_MajAll,TurkeyPen_MajAll,CNNAM_MajAll,NAM_MajAll";
        
        
        System.out.println("LCJ - begine testColumnsToBinaryTest ...");
        new ColumnsToBinarySNPOnlyTablePlugin()
        .inputFile(inputFile)
        .refDir(refDir)
        .outBase(outBase)
        //.colsFloat(floatColNames)
        .colsInt(intColNames)
        .colsAllele(alleleColNames)
        //.colsByte(byteColNames)
        //.range(true)  // remove if using position specific file vs ranges
        //.colsShort(shortColNames)
        .performFunction(null);
        
//        new ColumnsToBinaryPlugin()
//        .inputFile(inputFile)
//        .refDir(refDir)
//        .outBase(outBase)
//        .colsFloat(floatColNames)
//        .colsInt(intColNames)
//        .range(true)  // remove if using position specific file vs ranges
//        //.colsShort(shortColNames)
//        .performFunction(null);
//        
        //Columns are now created, need to put into monetdb
        // before loading, must ftp them to cbsudc01, which I haven't figured out
        // yet how to do programmatically      
        
        // Connect to DB, write script 
        
        System.out.println("\nFInished creating columns via testColumnsToBinaryTest !");
    }

    
 
}
