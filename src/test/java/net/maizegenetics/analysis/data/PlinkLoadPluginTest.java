/**
 * Tests importing/exporting of PLINK files by comparing
 * to hapmap files.  
 * - import the PLINK files, convert the data to a GenotypeTable
 * - export the GenotypeTable as a hapmap file
 * - import the new hapmap file, convert again to GenotypeTable 
 * - export this table as PLINK files
 * - import the new PLINK files
 * - comparisons along the way to verify the data remains consistent
 */
package net.maizegenetics.analysis.data;

import static org.junit.Assert.*;

import java.util.Iterator;
import java.util.List;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.AlignmentTestingUtils;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import org.junit.Test;

/**
 * @author lcj34
 *
 */
public class PlinkLoadPluginTest {

	static final String plinkPedFile=GeneralConstants.DATA_DIR+"CandidateTests/mdp_genotype.plk.ped";
	static final String plinkMapFile=GeneralConstants.DATA_DIR+"CandidateTests/mdp_genotype.plk.map";

	@Test
	public void testPlinkLoadPlugin() throws Exception {
		System.out.println("Running testPlinkLoadTestPlugin");
		PlinkLoadPlugin plugin = new PlinkLoadPlugin(null,false);
		plugin.pedFile(plinkPedFile);
		plugin.mapFile(plinkMapFile);
		DataSet dataPlink = plugin.performFunction(null);

		List<Datum> inputsPL = dataPlink.getDataOfType(GenotypeTable.class);
		
		if (inputsPL == null  || inputsPL.size() == 0) {
			fail("ERROR:PlinkLoadPluginTest: inputsPL are NULL");
			return;
		}
		GenotypeTable genotypePL = (GenotypeTable) inputsPL.get(0).getData();

		// Export PLINK table as hapmap
		ExportUtils.writeToHapmap(genotypePL, GeneralConstants.TEMP_DIR+"plinkToHapmap");

		// Import the new hapmap file 
		FileLoadPlugin flp = new FileLoadPlugin(null, false);
		String filename = GeneralConstants.TEMP_DIR + "plinkToHapmap.hmp.txt";
		flp.setOpenFiles(new String[]{filename});
		DataSet dataHapmap = flp.performFunction(null);

		// Convert to a GenotypeTable for comparison
		List<Datum> inputsHM = dataHapmap.getDataOfType(GenotypeTable.class);
		if (inputsHM == null  || inputsHM.size() == 0) {
			fail(" ERROR PlinkLoadPluginTest: HM inputs are NULL");
			return;
		}
		GenotypeTable genotypeHM = (GenotypeTable) inputsHM.get(0).getData();

		// Compare PLINK/HapMap values for equality
		assertEquals("Plink DataSet and hapmap Dataset sizes are different ", dataPlink.getSize(), dataHapmap.getSize());		         
		comparePositions(genotypePL, genotypeHM);
		AlignmentTestingUtils.alignmentsEqual(genotypePL, genotypeHM);

		// Export the hapmap file as plink, import new plink files  
		ExportUtils.writeToPlink(genotypeHM, GeneralConstants.TEMP_DIR+"hapmapToPlink", '\t');	          

		plugin.pedFile(GeneralConstants.TEMP_DIR+"hapmapToPlink.plk.ped");
		plugin.mapFile(GeneralConstants.TEMP_DIR+"hapmapToPlink.plk.map");
		DataSet dataPlink2 = plugin.performFunction(null);

		List<Datum> inputsPL2 = dataPlink2.getDataOfType(GenotypeTable.class);
		if (inputsPL == null  || inputsPL.size() == 0) {
			fail("ERROR PlinkLoadPluginTest: inputsPL2 are NULL");
			return;
		}
		GenotypeTable genotypePL2 = (GenotypeTable) inputsPL2.get(0).getData();

		// Compare genotype tables from original and new PLINK files for equality 
		assertEquals("Plink DataSet and hapmap Dataset sizes are different ", dataPlink.getSize(), dataPlink2.getSize());		   
		comparePositions(genotypePL, genotypePL2);
		AlignmentTestingUtils.alignmentsEqual(genotypePL,  genotypePL2);
	}

	public void comparePositions(GenotypeTable genotype1, GenotypeTable genotype2) {
		// compare the Positions 
		PositionList onePos = genotype1.positions();
		PositionList twoPos = genotype2.positions();
		assertEquals("Expected and Actual PositionLists are not the same size", onePos.size(), twoPos.size());
		Iterator<Position> ExpIter = onePos.iterator();
		Iterator<Position> ActIter = twoPos.iterator();
		while (ExpIter.hasNext()) {
			Position expPos = ExpIter.next();
			Position actPos = ActIter.next();
			assertEquals("Expected and Actual Positions are not the same for postion "+expPos.getSNPID()+" vs "+actPos.getSNPID(), expPos.compareTo(actPos), 0);
		}
	}
}
