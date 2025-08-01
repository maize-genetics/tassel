/**
 * 
 */
package net.maizegenetics.analysis.gbs.v2;

import junit.framework.Assert;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.tag.TagData;
import net.maizegenetics.dna.tag.TagDataSQLite;

import org.junit.Test;

/**
 * @author lcj34
 *
 */
public class SNPCutPosTagVerificationPluginTest {

	// Test assumes the Evaluate Pipeline was run prior to 
	// running this test - GBSv2.db has data.
	@Test
	public void SNPCutPosTagVerificationTest() throws Exception{
		System.out.println("Running SNPCutPosTagVerificationTest");

		System.out.println("\nSNPCutPosTagVerificationTest:  creating cut positions file\n");
		TagData tagData =new TagDataSQLite(GBSConstants.GBS_GBS2DB_FILE);
		//TagData tagData =new TagDataSQLite("/Users/lcj34/Documents/MissouriWorkshopData/GBSv2/expectedResults/GBSv2.db");
		PositionList dbCutPositions=tagData.getTagCutPositions(true);
		int numCutFiles = dbCutPositions.size() < 5 ? dbCutPositions.size() : 5;
		
		for (int idx = 0; idx < numCutFiles; idx++) {
			Position cPosition = dbCutPositions.get(idx);
			int aCutPosition = cPosition.getPosition();
			byte strand = cPosition.getStrand();
			String sChromName = cPosition.getChromosome().getName();
			String cutOutFile = GBSConstants.GBS_TEMP_DIR + "/cutPositionData" + idx +".txt";
			new SNPCutPosTagVerificationPlugin()
			.inputDB(GBSConstants.GBS_GBS2DB_FILE)
			//.inputDB("/Users/lcj34/Documents/MissouriWorkshopData/GBSv2/expectedResults/GBSv2.db")
			.cutOrSnpPosition(aCutPosition)
			.chrom(sChromName)
			.positionType("cut")
			.strand(strand)
			.outputFile(cutOutFile)
			//.outputFile("/Users/lcj34/notes_files/gbsv2/qi_ref/debug_March14/cutPositionData" + idx + ".txt")
			//.outputFile("/Users/lcj34/notes_files/gbsv2/debug_problems/tas758_debugPluginPullCutPositionTags/tmpSNPFiles/snpPositionData" + idx + ".txt")
			.performFunction(null);
		}				

		System.out.println("\nSNPCutPosTagVerificationTest:  creating SNP positions file\n");
		PositionList dbSNPPositions = tagData.getSNPPositions();
		int numFiles = dbSNPPositions.size() < 5 ? dbSNPPositions.size() : 5;
		for (int idx = 0; idx < numFiles; idx ++) {
			Position sPosition = dbSNPPositions.get(idx);
			int aSNPPosition = sPosition.getPosition();
			byte strand = sPosition.getStrand();
			String sChromName = sPosition.getChromosome().getName();
			String snpOutFile = GBSConstants.GBS_TEMP_DIR + "/snpPositionData" + idx +".txt";
			new SNPCutPosTagVerificationPlugin()
			.inputDB(GBSConstants.GBS_GBS2DB_FILE)
			//.inputDB("/Users/lcj34/Documents/MissouriWorkshopData/GBSv2/expectedResults/GBSv2.db")
			.cutOrSnpPosition(aSNPPosition)
			.chrom(sChromName)
			.positionType("snp")
			.strand(strand)
			.outputFile(snpOutFile)
			//.outputFile("/Users/lcj34/notes_files/gbsv2/qi_ref/debug_March14/snpPositionData" + idx + ".txt").
			//.outputFile("/Users/lcj34/notes_files/gbsv2/debug_problems/tas758_debugPluginPullCutPositionTags/tmpSNPFiles/snpPositionData" + idx + ".txt")
			.performFunction(null);
		}

		((TagDataSQLite)tagData).close();
		System.out.println("SNPCutPosTagVerificationTest finished successfully !!!");
	}
	
	       @Test
	        public void SNPCutPosTagSinglePositionTest() throws Exception{
	                System.out.println("Running SNPCutPosTagSinglePositionTest");

	                String cutOutFile = GBSConstants.GBS_TEMP_DIR + "/cutData66648.txt";
                        new SNPCutPosTagVerificationPlugin()
                        .inputDB(GBSConstants.GBS_GBS2DB_FILE)
                        .cutOrSnpPosition(66648)
                        .chrom("9")
                        .positionType("cut")
                        .strand((byte)0)
                        .outputFile(cutOutFile)
                        //.outputFile("/Users/lcj34/notes_files/gbsv2/debug_problems/tas758_debugPluginPullCutPositionTags/tmpSNPFiles/snpPositionData" + idx + ".txt")
                        .performFunction(null);

                        String snpOutFile = GBSConstants.GBS_TEMP_DIR + "/snpData66696.txt";
                        new SNPCutPosTagVerificationPlugin()
                        .inputDB(GBSConstants.GBS_GBS2DB_FILE)
                        .cutOrSnpPosition(66696)
                        .chrom("9")
                        .positionType("snp")
                        .strand((byte)1)
                        .outputFile(snpOutFile)
                        //.outputFile("/Users/lcj34/notes_files/gbsv2/debug_problems/tas758_debugPluginPullCutPositionTags/tmpSNPFiles/snpPositionData" + idx + ".txt")
                        .performFunction(null);

	               
	                System.out.println("SNPCutPosTagSinglePositionTest finished successfully !!!");
	        }
}
