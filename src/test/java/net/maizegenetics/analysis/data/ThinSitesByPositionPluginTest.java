/**
 * 
 */
package net.maizegenetics.analysis.data;

import java.nio.file.Files;
import java.nio.file.Paths;
import static org.junit.Assert.*;

import org.junit.Test;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.TaxaList;

/**
 * The data for this test was created as below:
 *  -  inputfile:  The the h5 file created by running ProductionSNPCallerPluginV2Test
 *         after EvaluateSNPCallQualityOfPipelineTest was run using .01 MAF pipeline
 *         on Nov 3, 2015 (file may be different if created at a different time).  THis
 *         file has been stored in the ExpectedResults/GenotypeTableFiles directory.
 *         We could consider storing the .vcf output file from ProductionSNPCallerPluginV2
 *         but it is much larger.  
 *  - expectedResults100:  Using VCF tools and the VCF version of the file created above, ran
 *         this command:
 *      vcftools --vcf  maize_chr0_10.vcf --thin 100 --recode --recode-INFO-all --out thin
 *  - expected Results40000: created using the same command as above, but with
 *     --thin 40000 
 *     
 *  The test verifies the code from this plugin creates the same results as vcf tools
 *  running the commands above.
 *  
 *  
 * @author lcj34
 *
 */
public class ThinSitesByPositionPluginTest {
    String actualResults100 = GeneralConstants.TEMP_DIR + "thin100.recode.vcf";
    String actualResults40000 = GeneralConstants.TEMP_DIR + "thin40000.recode.vcf";
    String actualResultsLARGE = GeneralConstants.TEMP_DIR + "thinVariable.recode.vcf";
    String expectedResults100 = GeneralConstants.DATA_DIR + "GenotypeTableTests/maize_chr9_10thin100.recode.vcf";
    String expectedResults40000 = GeneralConstants.DATA_DIR + "GenotypeTableTests/maize_chr9_10thin40000.recode.vcf";
 
    // String inputFile = GeneralConstants.DATA_DIR + "GenotypeTableTests/maize_chr9_10.vcf";
    // See comment above for how the input file below was created.  I stored the .h5
    // version of the input file rather than the .vcf version becuase it was much smaller
    // to put into git.
    String inputFile = GeneralConstants.DATA_DIR + "GenotypeTableTests/maize_chr9_10.h5";
    
    @Test
    public void ThinSitesByPosition100() throws Exception {       
        System.out.println("Starting test ThinSitesByPositionPluginTest");
        
        //GenotypeTable result = ImportUtils.readGuessFormat(inputFile);
        GenotypeTable result = GenotypeTableBuilder.getInstance(inputFile);
        //GenotypeTable result = (GenotypeTable)ImportUtils.readDataSet(inputFile);
        ThinSitesByPositionPlugin thinSitesByPosition = new ThinSitesByPositionPlugin(null, false);
        thinSitesByPosition.outfile(actualResults100);
        thinSitesByPosition.minDist(100);
        
        // Since becoming GUI compatible, this line now calls
        // the plugin
        DataSet myDS = thinSitesByPosition.performFunction(new DataSet(new Datum("inputFile", result, null),null));
        
//      new ThinSitesByPositionPlugin()
//      .infile(inputFile)
//      .outfile(actualResults100)
//      .minDist(100)
//      .performFunction(null); 
        
        //Read the genotype files - create genotype tables
        GenotypeTable actualGenos100 = ImportUtils.readFromVCF(actualResults100, null);
        GenotypeTable expectedGenos100 = ImportUtils.readFromVCF(expectedResults100,null);

        // Get taxa list from genotype table
        TaxaList actualTaxaList100 = actualGenos100.taxa();
        int nActualTaxa = actualTaxaList100.numberOfTaxa();
        
        TaxaList expectedTaxaList100 = expectedGenos100.taxa();
        int nExpectedTaxa = expectedTaxaList100.numberOfTaxa();
        
        System.out.println("\nExpected number of taxa: " + nExpectedTaxa 
                + ", actual number of taxa: " + nActualTaxa);
        assertEquals(nExpectedTaxa, nActualTaxa); // LCJ - uncomment, jus wnat to see results of next assert
        
        PositionList actualPL100 = actualGenos100.positions();
        int nActualPL100 = actualPL100.size();
        
        PositionList expectedPL100 = expectedGenos100.positions();
        int nExpectedPL100 = expectedPL100.size();
        
        System.out.println("Expected number of positions: " + nExpectedPL100 
                + ", actual number of positions: " + nActualPL100);
        assertEquals(nExpectedPL100, nActualPL100);

        //Clean up
        Files.deleteIfExists(Paths.get(actualResults100));
    }  
    
    @Test
    public void ThinSitesByPosition40000() throws Exception {       
        System.out.println("Starting test ThinSitesByPosition40000");
        
        GenotypeTable result = GenotypeTableBuilder.getInstance(inputFile);
       // GenotypeTable result = ImportUtils.readGuessFormat(inputFile);
        ThinSitesByPositionPlugin thinSitesByPosition = new ThinSitesByPositionPlugin(null, false);
        thinSitesByPosition.outfile(actualResults40000);
        thinSitesByPosition.minDist(40000);
        
        DataSet myDS = thinSitesByPosition.performFunction(new DataSet(new Datum("inputFile", result, null),null));
        
//        new ThinSitesByPositionPlugin()
//        .infile(inputFile)
//        .outfile(actualResults40000)
//        .minDist(40000)
//        .performFunction(null); 
        
        //Read the genotype files - create genotype tables
        GenotypeTable actualGenos40000 = ImportUtils.readFromVCF(actualResults40000,null);
        GenotypeTable expectedGenos40000 = ImportUtils.readFromVCF(expectedResults40000,null);
  
        // Get taxa list from genotype table
        TaxaList actualTaxaList40000 = actualGenos40000.taxa();
        int nActualTaxa = actualTaxaList40000.numberOfTaxa();
        
        TaxaList expectedTaxaList40000 = expectedGenos40000.taxa();
        int nExpectedTaxa = expectedTaxaList40000.numberOfTaxa();
        
        System.out.println("\nExpected number of taxa: " + nExpectedTaxa 
                + ", actual number of taxa: " + nActualTaxa);
        assertEquals(nExpectedTaxa, nActualTaxa);
        
        PositionList actualPL40000 = actualGenos40000.positions();
        int nActualPL40000 = actualPL40000.size();
        
        PositionList expectedPL40000 = expectedGenos40000.positions();
        int nExpectedPL40000 = expectedPL40000.size();
        
        System.out.println("Expected number of positions: " + nExpectedPL40000 
                + ", actual number of positions: " + nActualPL40000);
        assertEquals(nExpectedPL40000, nActualPL40000);

        //Clean up
        Files.delete(Paths.get(actualResults40000));
    }   
    
    @Test
    public void ThinSitesByPositionLARGE() throws Exception {       
        System.out.println("Starting test ThinSitesByPositionLARGE");
        
        GenotypeTable result = GenotypeTableBuilder.getInstance(inputFile);
        //GenotypeTable result = ImportUtils.readGuessFormat(inputFile);
        ThinSitesByPositionPlugin thinSitesByPosition = new ThinSitesByPositionPlugin(null, false);
        thinSitesByPosition.outfile(actualResultsLARGE);
        thinSitesByPosition.minDist(2147483600);
        
        DataSet myDS = thinSitesByPosition.performFunction(new DataSet(new Datum("inputFile", result, null),null));

//        new ThinSitesByPositionPlugin()
//        .infile(inputFile)
//        .outfile(actualResults)
//        .minDist(2147483600)
//        .performFunction(null); 
        
        //Read the genotype files - create genotype tables
        GenotypeTable actualGenos= ImportUtils.readFromVCF(actualResultsLARGE,null);
        PositionList actualPL = actualGenos.positions();
        int nActualPL = actualPL.size();
        
        // We gave it a really large number.  There are just 2 chromosomes.
        // We expect this will return just the 1st position for each chromosomes,
        // resulting in number of positions = 2.
        
        System.out.println("\nNumber of positions found: " + nActualPL);
        assertEquals(nActualPL,2);

        //Clean up
        Files.deleteIfExists(Paths.get(actualResultsLARGE));
    }        
}
