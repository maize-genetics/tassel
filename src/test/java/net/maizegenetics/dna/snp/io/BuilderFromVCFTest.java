package net.maizegenetics.dna.snp.io;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;

import org.junit.Ignore;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.io.BuilderFromVCF;
import net.maizegenetics.util.BuilderFromVCFUtil;

import org.apache.commons.io.FileUtils;

public class BuilderFromVCFTest {

	//String directoryName = "/Users/zrm22/Desktop/VCFTestData/VCF_files/TestFiles/";
	String directoryName = GeneralConstants.DATA_DIR+"CandidateTests/VCFFiles/";
	
	@Rule
	public ExpectedException thrown = ExpectedException.none();

	
	//TESTS TODO
	//-Verify that the information is the same during import/export
	//-Verify information is same during sorting
	//-Verify GenotypeTables are correctly generated
	
	//Currently Ignoring tests as the file paths will not work unless set up correctly.
	
	//Test to make sure it can load correct VCF files
	//SuccessCondition: BuilderFromVCF generates a non-null GenotypeTable
	@Test
	public void standardVCFTest() {
		String fileName = directoryName+"correctVCF1.vcf";
		BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
		GenotypeTable table = builder.buildAndSortInMemory();
		assertNotNull("Standard VCF Test Table should not be null",table);
	}

	//Test to make sure it can handle ./. as the GT field
	//SuccessCondition: BuilderFromVCF generates a non-null GenotypeTable
	@Test
	public void missingAlleleTest() {
		String fileName = directoryName+"missingAlleleVCF1.vcf";
		BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
		GenotypeTable table = builder.buildAndSortInMemory();
		assertNotNull("Standard VCF Test Table should not be null",table);
	}
	
	//Test to make sure it catches when a haploid is encountered
	//SuccessCondition: BuilderFromVCF throws an exception with an appropriate error message
	//TODO fix this as we can now load haploid entries
	@Ignore
	@Test
	public void tas_509Test() {
	    String fileName = directoryName+"tas_509VCF1.vcf";
	    BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
	    thrown.expect(IllegalStateException.class);
	    thrown.expectMessage("java.lang.IllegalStateException: Error Processing VCF block: Found haploid information for the element: .:.:.:.:..\nExpected a diploid entry.");
	    GenotypeTable table = builder.buildAndSortInMemory();
	    /*
	    try{
		GenotypeTable table = builder.buildAndSortInMemory();
		fail("Should have thrown an exception at the string \".:.:.:.:.\"");
	    }catch(IllegalStateException e) {
		    assertEquals(e.getMessage(),"java.lang.IllegalStateException: Error Processing VCF block: Found haploid information for the element: .:.:.:.:..\nExpected a diploid entry.");
		}
		*/
	}
	
	//Test for Missing FORMAT Tag
	@Test
	public void formatTest() {
	    String fileName = directoryName+"MissingFormat.vcf";
	    BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
	    thrown.expect(IllegalStateException.class);
	    thrown.expectMessage("java.lang.IllegalStateException: Error Processing VCF Block: Missing FORMAT tag.");
	    GenotypeTable table = builder.buildAndSortInMemory();
	    /*
            try{
                GenotypeTable table = builder.buildAndSortInMemory();
                fail("Should have thrown an exception when FORMAT is hit");
            }catch(IllegalStateException e) {
                    assertEquals(e.getMessage(),"java.lang.IllegalStateException: Error Processing VCF Block: Missing FORMAT tag.");
            }
            */
	}
	
	//Test for Incorrectly ordered FORMAT Tag
	@Test
	public void formatOrderedTag() {
	    String fileName = directoryName+"GTFieldIncorrect.vcf";
            BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
            thrown.expect(IllegalStateException.class);
            thrown.expectMessage("java.lang.IllegalStateException: Error Processing VCF Block: GT field is not in first position of FORMAT.");
            GenotypeTable table = builder.buildAndSortInMemory();
            /*
            try{
                GenotypeTable table = builder.buildAndSortInMemory();
                fail("Should have thrown an exception when FORMAT is hit and GT is not first");
            }catch(IllegalStateException e) {
                    assertEquals(e.getMessage(),"java.lang.IllegalStateException: Error Processing VCF Block: GT field is not in first position of FORMAT.");
            }
            */
	}
	
	//Test to make sure it catches when a file is unsorted by position
	//SuccessCondition: BuilderFromVCF throws an IllegalStateException
	@Test
	public void unSortedTest() {
	    String fileName = directoryName+"unsortedVCF1.vcf";
	    BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
	    thrown.expect(IllegalStateException.class);
            thrown.expectMessage("BuilderFromVCF: Ordering incorrect. VCF file must be ordered by position. Please first use SortGenotypeFilePlugin to correctly order the file.");
            GenotypeTable table = builder.build();
            
            /*
	    try{
	        GenotypeTable table = builder.build();
	        fail("Should have thrown an exception as the file is unsorted by position");
	    }catch(IllegalStateException e) {
	        assertEquals(e.getMessage(),"BuilderFromVCF: Ordering incorrect. VCF file must be ordered by position. Please first use SortGenotypeFilePlugin to correctly order the file.");
	    }
	    */
	}
	
	@Test
	public void spacedHeaderTest() {
	    String fileName = directoryName+"correctHeaderVCF1.vcf";
	    BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
	    GenotypeTable table = builder.buildAndSortInMemory();
	    assertEquals("Failed Header with Spaces Test",table.numberOfTaxa(),1);
	}
	
	@Test
	//Basically round-trips a VCF file with INFO Tags
	public void infoTagPerSiteTest() {
	    //Read in file with INFO Tags
	    String fileName = directoryName+"infoTagVCF1.vcf";
	    BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
	    //Build GenotypeTable
	    GenotypeTable table = builder.buildAndSortInMemory();
	    
	    //Export GenotypeTable to VCF
	    String outputFileName = GeneralConstants.TEMP_DIR+"infoTagVCFExport.vcf";
            ExportUtils.writeToVCF(table, outputFileName, true);
            
            //Read in Exported and compare to a Input VCF file
            try {
                boolean filesMatch = FileUtils.contentEquals(new File(fileName), new File(outputFileName));              
                assertTrue("Exported File does not match Original File.  INFO tags are not consistent",filesMatch);
                
                //Check to see if Exported file matches expected Result
                String outputFileNameCorrect = GeneralConstants.EXPECTED_RESULTS_DIR+"infoTagVCFExpected.vcf";
                filesMatch = FileUtils.contentEquals(new File(outputFileName), new File(outputFileNameCorrect));
                assertTrue("Exported File does not match Expected File. INFO tags are not consistent",filesMatch);
                
                //Could check to see if we remove info that it will not match, but this is trivial as it is the inverse of the previous test
            }
            catch(Exception e) {
                System.out.println(e);
            }
	}
	//Test to make sure we can roundtrip a vcf file with indels at certain positions
	@Test
	public void indelTest() {
	    String fileName = directoryName + "indelTestVCFInput.vcf";
	    BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
	    
	    //Build GenotypeTable
            GenotypeTable table = builder.buildAndSortInMemory();
            
	    //Export GenotypeTable to VCF
	    String outputFileName = GeneralConstants.TEMP_DIR + "indelTextVCFExport.vcf";
	    ExportUtils.writeToVCF(table, outputFileName, true);   
	
	    //Read in Exported and compare to a Input VCF file
            try {
                boolean filesMatch = FileUtils.contentEquals(new File(fileName), new File(outputFileName));              
                assertTrue("Exported File does not match Original File.  Indel Processing is incomplete",filesMatch);
            }
            catch(Exception e) {
                System.out.println(e);
            }
	}
	
	@Test
	public void refNotSeenTest() {
	    String fileName = directoryName + "refNotSeenInput.vcf";
	    BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
            
            //Build GenotypeTable
            GenotypeTable table = builder.buildAndSortInMemory();
            
            //Export GenotypeTable to VCF
            String outputFileName = GeneralConstants.TEMP_DIR + "refNotSeenVCFExport.vcf";
            ExportUtils.writeToVCF(table, outputFileName, true);   
        
            //Read in Exported and compare to a Input VCF file
            try {
                boolean filesMatch = FileUtils.contentEquals(new File(fileName), new File(outputFileName));              
                assertTrue("Exported File does not match Original File.  Ref not represented but is swapped for major allele",filesMatch);
            }
            catch(Exception e) {
                System.out.println(e);
            }
	}
	
	@Ignore
	@Test
	public void singleLineIndelTest() {
	    //Test createSmallCallTable
	    ArrayList<String> variants = new ArrayList<String>();
	   

	    variants.add("T");
	    variants.add("TGATATGTAGCGGTTCC");
	    variants.add("TGATATCTAGCGGTTCC");
	    variants.add("G");
	    variants.add("TAATATGTACCTGTTCC");
	    variants.add("TGATATATAGCGGTTCC");
	    variants.add("TGATATGTAGCAGTTCC");
	    char[][] smallCallTable = BuilderFromVCFUtil.createSmallCallTable(variants);
	    
	    for(char[] variant : smallCallTable) {
	        System.out.println(Arrays.toString(variant));
	    }
	}
	
	@Ignore
	@Test
	public void kellyTest() {
	    //String fileName = "/Users/zrm22/Desktop/KellyTests/merged_flt_c1.imputed.vcf";
	    String fileName = "/Users/zrm22/Desktop/KellyTests/TestVCFFileNewExpanded.vcf";
            
	    String outName = "/Users/zrm22/Desktop/KellyTests/output.h5";
	    BuilderFromVCF.getBuilder(fileName).convertToHDF5(outName).keepDepth().build();
	}
	
	@Test
	public void testBadFormatedDataLines() {
	    String fileName = directoryName + "badDataLineFormatted.vcf";
            BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
            
            thrown.expect(IllegalStateException.class);
            //thrown.expectMessage("java.lang.IllegalStateException: Error Processing VCF Block: Incorrect Formatting of a Data Row.\nThe Line which causes the error is this:\nThis Line should trigger error.");
            GenotypeTable table = builder.buildAndSortInMemory();
	}
	
	@Test
        //Basically round-trips a VCF file there should be a change with the INFO tag
        public void infoTagDPChangingPerSiteTest() {
            //Read in file with INFO Tags
            String fileName = directoryName+"infoTagDPChangingVCF1.vcf";
            BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName).keepDepth();
            //Build GenotypeTable
            GenotypeTable table = builder.buildAndSortInMemory();
            
            //Export GenotypeTable to VCF
            String outputFileName = GeneralConstants.TEMP_DIR+"infoTagDPChangingVCFExport.vcf";
            ExportUtils.writeToVCF(table, outputFileName, true);
            
            //Read in Exported and compare to a Input VCF file
            try {
                
                //Check to see if Exported file matches expected Result
                String outputFileNameCorrect = GeneralConstants.EXPECTED_RESULTS_DIR+"infoTagDFChangingVCFExpected.vcf";
                boolean filesMatch = FileUtils.contentEquals(new File(outputFileName), new File(outputFileNameCorrect));
                assertTrue("Exported File does not match Expected File. INFO tags are not consistent",filesMatch);
                
                //Could check to see if we remove info that it will not match, but this is trivial as it is the inverse of the previous test
            }
            catch(Exception e) {
                System.out.println(e);
            }
        }
}
