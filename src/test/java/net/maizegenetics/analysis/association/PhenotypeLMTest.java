package net.maizegenetics.analysis.association;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportTestUtils;
import net.maizegenetics.util.TableReportUtils;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class PhenotypeLMTest {
	private PhenotypeLM myModel;
	private static float tol = 1e-4f;
	
	@Before
	public void setUp() throws Exception {
		FileLoadPlugin flp = new FileLoadPlugin(null, false);
		String filename = GeneralConstants.DATA_DIR + "CandidateTests/earht_with_rep.txt";
		flp.setOpenFiles(new String[]{filename});
		DataSet ds = flp.performFunction(null);
		
		myModel = new PhenotypeLM(ds.getData(0));
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testBlues() {
		float[] expectedValues = new float[] {94.0315456025016f, 64.7254997617167f, 45.6519525059068f, 29.7524935013531f, 62.5555279926411f, 38.8168399291231f, 
				67.7808311612335f, 70.8414557486813f, 45.2450936651001f, 63.8476743870665f, 66.9296408321963f, 137.361806914807f, 83.4337287433577f, 78.2158230657591f, 
				66.1571330200281f, 58.7176743870664f, 55.9297855278656f, 43.3691463840685f, 68.9777451837748f, 66.3652248691787f, 77.0025258082017f, 52.3935891603886f, 
				93.7872183320187f, 72.7176743870665f, 45.7176743870664f, 51.9557860580037f, 85.6218663548891f, 62.9402081302404f, 67.2176743870665f, 62.8193914236268f}; //values calculate in R
		Phenotype myBlues = myModel.blues();
		NumericAttribute blueAttribute = (NumericAttribute) myBlues.attribute(1);
		float[] observedValues = blueAttribute.floatValues();
		assertArrayEquals("", expectedValues, observedValues, tol);
	}

	@Test
	public void testReport() {
		TableReport myReport = myModel.report();
		
		String expectedResults = GeneralConstants.EXPECTED_RESULTS_DIR + "phenotypeLMreport.txt";

		//to save the results to use for validation, uncomment the following line
//		TableReportUtils.saveDelimitedTableReport(myReport, "\t", new File(expectedResults));

		TableReport expectedSiteReport = TableReportUtils.readDelimitedTableReport(expectedResults, "\t");
		
		File tempFile; 
		try {
			tempFile = File.createTempFile("test", ".txt");
			TableReportUtils.saveDelimitedTableReport(myReport, "\t", tempFile);
			TableReport observedSiteReport = TableReportUtils.readDelimitedTableReport(tempFile.getPath(), "\t");
			TableReportTestUtils.compareTableReports(expectedSiteReport, observedSiteReport);
		} catch(IOException e) {
			e.printStackTrace();
		}
	}

}
