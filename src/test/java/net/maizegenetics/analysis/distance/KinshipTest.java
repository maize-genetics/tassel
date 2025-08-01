package net.maizegenetics.analysis.distance;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportUtils;

import org.junit.Test;

public class KinshipTest {

    /**
     * Test of getDistanceMatrix method, of class DistanceMatrixPlugin.
     */
    @Test
    public void testKinship() throws IOException {
        System.out.println("Testing Endelman Kinship Matrix ...");

        File input = new File(TutorialConstants.HAPMAP_FILENAME);
        System.out.println("   Input File: " + input.getCanonicalPath());
        if (!input.exists()) {
            fail("Input File: " + TutorialConstants.HAPMAP_FILENAME + " doesn't exist.");
        }
        GenotypeTable inputAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
        System.out.println("   Input Alignment: Taxa Count: " + inputAlign.numberOfTaxa() + "  Site Count: " + inputAlign.numberOfSites());

        GenotypeTable filteredGeno = GenotypeTableUtils.removeSitesBasedOnFreqIgnoreMissing(inputAlign, 0.01, 1, 250);
        //remove the marker with 4 alleles
        filteredGeno = FilterGenotypeTable.getInstanceRemoveSiteNames(filteredGeno, new String[]{"PZB00063.1"});

        //export this alignment for validation (run once)
        //ExportUtils.writeToHapmap(filteredGeno, "/Users/pbradbury/Documents/projects/tassel/validation/kinshipTestData.hmp.txt");
        System.out.println("Testing Endelman method.");
        DistanceMatrix kin = EndelmanDistanceMatrix.getInstance(filteredGeno);

        String saveFilename = GeneralConstants.EXPECTED_RESULTS_DIR + "KinshipTestDistanceMatrix.txt";

        //File expectedEndelmanResults = new File(saveFilename);
        //export the kinship matrix for validation (run once)
        //TableReportUtils.saveDelimitedTableReport(kin.getDm(), "\t", expectedEndelmanResults);
        //import the table report for comparing with the expected results (which have been validated in R)
        TableReport expectedResults = TableReportUtils.readDelimitedTableReport(saveFilename, "\t");
        TableReport actualResults = kin;

        long nrowsActual = actualResults.getRowCount();
        long nrowsExpected = expectedResults.getRowCount();
        assertEquals("Expected number of rows != actual number of rows", nrowsActual, nrowsExpected);
        int ncolsActual = actualResults.getColumnCount();
        int ncolsExpected = expectedResults.getColumnCount();
        assertEquals("Expected number of columns != actual number of columns", ncolsActual, ncolsExpected);
        double tol = 1e-5;
        for (long r = 0; r < nrowsExpected; r++) {
            for (int c = 1; c < ncolsExpected; c++) {
                double expectedValue = Double.parseDouble(expectedResults.getValueAt(r, c).toString());
                double actualValue = Double.parseDouble(actualResults.getValueAt(r, c).toString());
                double absdiff = Math.abs(expectedValue - actualValue);
                assertTrue(String.format("Expected and actual values different at row %d, column %d.", r, c), absdiff < tol);
            }
        }

    }

    @Test
    public void testEndelman() throws IOException {
        System.out.println("Testing Endelman Kinship Matrix ...");

        File input = new File(TutorialConstants.HAPMAP_FILENAME);
        System.out.println("   Input File: " + input.getCanonicalPath());
        if (!input.exists()) {
            fail("Input File: " + TutorialConstants.HAPMAP_FILENAME + " doesn't exist.");
        }
        GenotypeTable inputAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
        System.out.println("   Input Alignment: Taxa Count: " + inputAlign.numberOfTaxa() + "  Site Count: " + inputAlign.numberOfSites());

        GenotypeTable filteredGeno = GenotypeTableUtils.removeSitesBasedOnFreqIgnoreMissing(inputAlign, 0.01, 1, 250);
        //remove the marker with 4 alleles
        filteredGeno = FilterGenotypeTable.getInstanceRemoveSiteNames(filteredGeno, new String[]{"PZB00063.1"});

        //export this alignment for validation (run once)
        //ExportUtils.writeToHapmap(filteredGeno, "/Users/pbradbury/Documents/projects/tassel/validation/kinshipTestData.hmp.txt");
        System.out.println("Testing Endelman method.");
        DistanceMatrix kin = EndelmanDistanceMatrix.getInstance(filteredGeno);

        String saveFilename = GeneralConstants.EXPECTED_RESULTS_DIR + "KinshipTestDistanceMatrix.txt";
        
        //File expectedEndelmanResults = new File(saveFilename);
        //export the kinship matrix for validation (run once)
        //TableReportUtils.saveDelimitedTableReport(kin.getDm(), "\t", expectedEndelmanResults);
        //import the table report for comparing with the expected results (which have been validated in R)
        TableReport expectedResults = TableReportUtils.readDelimitedTableReport(saveFilename, "\t");
        TableReport actualResults = kin;

        long nrowsActual = actualResults.getRowCount();
        long nrowsExpected = expectedResults.getRowCount();
        assertEquals("Expected number of rows != actual number of rows", nrowsActual, nrowsExpected);
        int ncolsActual = actualResults.getColumnCount();
        int ncolsExpected = expectedResults.getColumnCount();
        assertEquals("Expected number of columns != actual number of columns", ncolsActual, ncolsExpected);
        double tol = 1e-5;
        for (long r = 0; r < nrowsExpected; r++) {
            for (int c = 1; c < ncolsExpected; c++) {
                double expectedValue = Double.parseDouble(expectedResults.getValueAt(r, c).toString());
                double actualValue = Double.parseDouble(actualResults.getValueAt(r, c).toString());
                double absdiff = Math.abs(expectedValue - actualValue);
                assertTrue(String.format("Expected and actual values different at row %d, column %d.", r, c), absdiff < tol);
            }
        }
    }
}
