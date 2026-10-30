package net.maizegenetics.dna.snp.io;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.AlignmentTestingUtils;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.io.BuilderFromVCF;
import net.maizegenetics.util.BuilderFromVCFUtil;

import org.apache.commons.io.FileUtils;

public class BuilderFromVCFTest {

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
	
	//Test that a haploid-style missing entry (".:.:.:.:.") loads as a missing genotype
	//SuccessCondition: BuilderFromVCF builds the table and the entry is unknown
	@Test
	public void tas_509Test() {
	    String fileName = directoryName+"tas_509VCF1.vcf";
	    BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
	    GenotypeTable table = builder.buildAndSortInMemory();
	    assertEquals(2, table.numberOfSites());
	    assertEquals("T", table.genotypeAsString(0, 0));
	    assertEquals(GenotypeTable.UNKNOWN_DIPLOID_ALLELE, table.genotype(0, 1));
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

	    String[] expected = {"T----------------", "TGATATGTAGCGGTTCC", "TGATATCTAGCGGTTCC", "G----------------",
	            "TAATATGTACCTGTTCC", "TGATATATAGCGGTTCC", "TGATATGTAGCAGTTCC"};
	    assertEquals(expected.length, smallCallTable.length);
	    for(int i = 0; i < expected.length; i++) {
	        assertEquals(expected[i], new String(smallCallTable[i]));
	    }
	}

	//Test that converting a VCF with depth straight to HDF5 keeps the genotypes
	//SuccessCondition: the HDF5 file has the same genotypes as the in-memory build
	@Test
	public void convertToHDF5KeepDepthGenotypesTest() throws Exception {
	    String fileName = TutorialConstants.TUTORIAL_DIR + "mdp_genotype_w_depth.vcf";
	    String outName = GeneralConstants.TEMP_DIR + "vcfWithDepthToHDF5Genotypes.h5";
	    Files.deleteIfExists(Paths.get(outName));

	    GenotypeTable expected = BuilderFromVCF.getBuilder(fileName).keepDepth().build();
	    BuilderFromVCF.getBuilder(fileName).convertToHDF5(outName).keepDepth().build();
	    AlignmentTestingUtils.alignmentsEqual(expected, GenotypeTableBuilder.getInstance(outName));
	}

	//Test that converting a VCF with depth straight to HDF5 keeps the allele depths
	//SuccessCondition: the HDF5 file has the same depths as the in-memory build
	@Test
	public void convertToHDF5KeepDepthTest() throws Exception {
	    String fileName = TutorialConstants.TUTORIAL_DIR + "mdp_genotype_w_depth.vcf";
	    String outName = GeneralConstants.TEMP_DIR + "vcfWithDepthToHDF5.h5";
	    Files.deleteIfExists(Paths.get(outName));

	    GenotypeTable expected = BuilderFromVCF.getBuilder(fileName).keepDepth().build();
	    BuilderFromVCF.getBuilder(fileName).convertToHDF5(outName).keepDepth().build();
	    GenotypeTable fromHDF5 = GenotypeTableBuilder.getInstance(outName);

	    assertTrue("Expected depth in the VCF fixture", expected.hasDepth());
	    assertTrue("HDF5 file lost the allele depths", fromHDF5.hasDepth());
	    for(int t = 0; t < expected.numberOfTaxa(); t++) {
	        for(int s = 0; s < expected.numberOfSites(); s++) {
	            assertArrayEquals("Depth differs at taxon " + t + " site " + s,
	                    expected.depthForAlleles(t, s), fromHDF5.depthForAlleles(t, s));
	        }
	    }
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
