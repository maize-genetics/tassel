/**
 * 
 */
package net.maizegenetics.analysis.monetdb;

import org.junit.Test;

/**
 * @author lcj34
 *
 */
public class ColumnsToBinaryFullGenomeTablePluginTest {
    @Test
    public void testColumnsToBinaryFullGenomeTest() throws Exception {
       // String inputFile = "/Volumes/Samsung_T1/machineLearningFiles/kelly_dbTestCode/test.FirstLast.txt";
        // In test file Ames_DTA.txt, the start/end are inclusive/exclusive, which is Kelly's default for start/end 
     
        String inputFile = "/Volumes/Samsung_T1/machineLearningFiles/eli_xoDensity/filesForScript";
        String refFile = "/Volumes/Samsung_T1/machineLearningFiles/refGenomeFiles/Zea_mays.AGPv3.29.dna.genome.fa";
        String outBase = "/Volumes/Samsung_T1/machineLearningFiles/monetdb_binaryFiles/eliXODensity_binaries/from_plugin/eli_xoDensity";
        String floatColNames = "mean,corrected";
        
        System.out.println("LCJ - begin testColumnsToBinaryTest ...");
        new ColumnsToBinaryFullGenomeTablePlugin()
        .inputFile(inputFile)
        .refFile(refFile)
        .outBase(outBase)
        .colsFloat(floatColNames)
        //.colsInt(intColNames)
        .range(true)  // remove if using position specific file vs ranges
        //.colsShort(shortColNames)
        .performFunction(null);
        
        //Columns are now created, need to put into monetdb temp table
        // The values were then compared with the values stored in the full genome table
        // that were created via Lynn's original method.  Spot checking showed transitions
        // between chromosomes was correct, transitions from null to values was correct.
        
        System.out.println("\nFInished creating columns via testColumnsToBinaryTest !");
    }
}
