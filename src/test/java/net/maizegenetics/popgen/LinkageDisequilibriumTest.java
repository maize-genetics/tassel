/*
 * LinkageDisequilibriumTest
 */
package net.maizegenetics.popgen;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportUtils;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.util.TableReportTestUtils;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

/**
 * TODO: Need to test r2, D', P with missing data and hets.
 *
 * @author terry
 */
public class LinkageDisequilibriumTest {

    public LinkageDisequilibriumTest() {
    }

    @Test
    public void testR2CalcsForHaplotypes() {
        System.out.println("Testing LinkageDisequilibrium Haplotypes Calcs()...");
        String[] seqs = {"AAA", "AAA", "AAA", "AAA", "AAA", "AAA", "AAA", "AAA", "AAA", "AAA", "AAA", "AAA",
            "CCC", "CCC", "CCC", "CCC", "CCC", "CCC", "CCC", "CCC", "CCC", "CCC", "CCC", "CCC", "CCC", "CCA", "NCA"};
        TaxaListBuilder tLB = new TaxaListBuilder();
        for (int i = 0; i < seqs.length; i++) {
            tLB.add(new Taxon("T" + i));
        }
        PositionListBuilder pLB = new PositionListBuilder();
        for (int i = 0; i < seqs[0].length(); i++) {
            pLB.add(new GeneralPosition.Builder(Chromosome.UNKNOWN, i).build());
        }
        GenotypeCallTable g = GenotypeCallTableBuilder.getInstance(seqs.length, seqs[0].length()).setBases(seqs).build();
        GenotypeTable theTBA = GenotypeTableBuilder.getInstance(g, pLB.build(), tLB.build());
        LinkageDisequilibrium allLD = new LinkageDisequilibrium(theTBA, 0, LinkageDisequilibrium.testDesign.All,
                0, null, false, 0, null, LinkageDisequilibrium.HetTreatment.Homozygous);
        allLD.run();
        double result = allLD.getRSqr(0, 1);
        System.out.println(result);
        assertEquals(1, result, 0.00001);
        result = allLD.getRSqr(0, 2);
        System.out.println(result);
        assertEquals(0.85714, result, 0.001);
        int[] f = {60, 10, 20, 10};
        result = LinkageDisequilibrium.calculateRSqr(f[0], f[1], f[2], f[3], 20);
        assertEquals(0.0476, result, 0.001);
        result = LinkageDisequilibrium.calculateDPrime(f[0], f[1], f[2], f[3], 20);
        assertEquals(0.2857, result, 0.001);
        result = LinkageDisequilibrium.calculateDPrime(f[1], f[2], f[3], f[0], 20);
        assertEquals(0.2857, result, 0.001);
        result = LinkageDisequilibrium.calculateRSqr(f[0], f[1], f[2], f[3], 200);  //too few to calculate
        System.out.println(result);
        assertEquals(Double.NaN, result, 0.001);
        seqs[1] = "CAA";
    }

    /**
     * Test class LinkageDisequilibrium.
     */
    @Test
    public void testLD() throws IOException {
        testLDforTreatment(LinkageDisequilibrium.HetTreatment.Haplotype);
        testLDforTreatment(LinkageDisequilibrium.HetTreatment.Homozygous);
    }

    private void testLDforTreatment(LinkageDisequilibrium.HetTreatment treatment) throws IOException {
        System.out.println("Testing LD...  " + treatment);

        File input = new File(TutorialConstants.HAPMAP_CHR_9_10_FILENAME);
        System.out.println("   Input File: " + input.getCanonicalPath());
        if (!input.exists()) {
            fail("Input File: " + TutorialConstants.HAPMAP_CHR_9_10_FILENAME + " doesn't exist.");
        }
        GenotypeTable inputAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_CHR_9_10_FILENAME, null);
        System.out.println("   Input Alignment: Taxa Count: " + inputAlign.numberOfTaxa() + "  Site Count: " + inputAlign.numberOfSites());

        System.out.println("   Testing LD All Sites by All Sites...");

        int windowSize = -1;
        LinkageDisequilibrium.testDesign ldType = LinkageDisequilibrium.testDesign.All;
        int testSite = -1;
        boolean isAccumulateResults = false;
        int numAccumulateIntervals = -1;
        int[] possibleSiteList = null;
        LinkageDisequilibrium allLD = new LinkageDisequilibrium(inputAlign, windowSize, ldType, testSite, null,
                isAccumulateResults, numAccumulateIntervals, possibleSiteList, treatment);
        allLD.run();

        int numSites = allLD.getSiteCount();
        assertEquals("LD Result Should have same number sites as Input Alignment", numSites, inputAlign.numberOfSites());

        for (int i = 0; i < numSites; i++) {
            for (int j = i; j < numSites; j++) {
                assertEquals("LD RSqr Differ: i: " + i + " j: " + j, allLD.getRSqr(i, j), allLD.getRSqr(j, i), 0.0);
                assertEquals("LD DPrime Differ: i: " + i + " j: " + j, allLD.getDPrime(i, j), allLD.getDPrime(j, i), 0.0);
                assertEquals("LD PVal Differ: i: " + i + " j: " + j, allLD.getPVal(i, j), allLD.getPVal(j, i), 0.0);
                assertEquals("LD SampleSize Differ: i: " + i + " j: " + j, allLD.getSampleSize(i, j), allLD.getSampleSize(j, i), 0.0);
            }
        }

        // Compare against previously calculated All by All results.
        TableReport expected = TableReportUtils.readDelimitedTableReport(GeneralConstants.EXPECTED_RESULTS_DIR + "LDResults_" + treatment + ".txt.gz", "\t");
        TableReportTestUtils.compareTableReportValues(expected, allLD, 0.01);

        // Compare Sliding Window against corrosponding All by All results
        System.out.println("   Testing LD Sliding Window...");

        windowSize = 50;
        ldType = LinkageDisequilibrium.testDesign.SlidingWindow;
        testSite = -1;
        isAccumulateResults = false;
        numAccumulateIntervals = -1;
        possibleSiteList = null;
        LinkageDisequilibrium slidingWinLD = new LinkageDisequilibrium(inputAlign, windowSize, ldType, testSite, null,
                isAccumulateResults, numAccumulateIntervals, possibleSiteList, treatment);
        slidingWinLD.run();

        long rowCount = slidingWinLD.getRowCount();
        int columnCount = expected.getColumnCount();
        assertEquals("Sliding Window LD Result Should have same number Columns as Expected", columnCount, slidingWinLD.getColumnCount());

        for (long i = 0; i < rowCount; i++) {
            Object[] actualRow = slidingWinLD.getRow(i);
            int site1 = (Integer) actualRow[2];
            int site2 = (Integer) actualRow[8];
            assertEquals("Sliding Window LD RSqr Differ from Full LD: site1: " + site1 + " site2: " + site2, allLD.getRSqr(site1, site2), slidingWinLD.getRSqr(site1, site2), 0.0);
            assertEquals("Sliding Window LD DPrime Differ from Full LD: site1: " + site1 + " site2: " + site2, allLD.getDPrime(site1, site2), slidingWinLD.getDPrime(site1, site2), 0.0);
            //TODO: Fix this assertion to work with new re-sized array capability of FisherExact.java (LCJ) 
            //assertEquals("Sliding Window LD PVal Differ from Full LD: site1: " + site1 + " site2: " + site2, allLD.getPVal(site1, site2), slidingWinLD.getPVal(site1, site2), 0.0);
            assertEquals("Sliding Window LD SampleSize Differ from Full LD: site1: " + site1 + " site2: " + site2, allLD.getSampleSize(site1, site2), slidingWinLD.getSampleSize(site1, site2), 0.0);
        }
        slidingWinLD = null;

        // Compare Site By All against corrosponding All by All results
        System.out.println("   Testing LD Site By All...");

        windowSize = -1;
        ldType = LinkageDisequilibrium.testDesign.SiteByAll;
        testSite = 3;
        isAccumulateResults = false;
        numAccumulateIntervals = -1;
        possibleSiteList = null;
        LinkageDisequilibrium siteByAllLD = new LinkageDisequilibrium(inputAlign, windowSize, ldType, testSite, null,
                isAccumulateResults, numAccumulateIntervals, possibleSiteList, treatment);
        siteByAllLD.run();

        rowCount = siteByAllLD.getRowCount();
        columnCount = expected.getColumnCount();
        assertEquals("Site By All LD Result Should have same number Columns as Expected", columnCount, siteByAllLD.getColumnCount());

        for (long i = 0; i < rowCount; i++) {
            Object[] actualRow = siteByAllLD.getRow(i);
            int site1 = (Integer) actualRow[2];
            int site2 = (Integer) actualRow[8];
            assertEquals("Site By All LD RSqr Differ from Full LD: site1: " + site1 + " site2: " + site2, allLD.getRSqr(site1, site2), siteByAllLD.getRSqr(site1, site2), 0.0);
            assertEquals("Site By All LD DPrime Differ from Full LD: site1: " + site1 + " site2: " + site2, allLD.getDPrime(site1, site2), siteByAllLD.getDPrime(site1, site2), 0.0);
            //TODO: Fix this assertion to match with new resizing array ability in FisherExact.java (LCJ)
            // assertEquals("Site By All LD PVal Differ from Full LD: site1: " + site1 + " site2: " + site2, allLD.getPVal(site1, site2), siteByAllLD.getPVal(site1, site2), 0.0);
            assertEquals("Site By All LD SampleSize Differ from Full LD: site1: " + site1 + " site2: " + site2, allLD.getSampleSize(site1, site2), siteByAllLD.getSampleSize(site1, site2), 0.0);
        }
        siteByAllLD = null;

        // Compare Site List against corrosponding All by All results
        System.out.println("   Testing LD Site List...");

        windowSize = -1;
        ldType = LinkageDisequilibrium.testDesign.SiteList;
        testSite = -1;
        isAccumulateResults = false;
        numAccumulateIntervals = -1;
        possibleSiteList = new int[]{3, 5, 49, 25, 10, 123, 35, 1, 36};
        LinkageDisequilibrium siteListLD = new LinkageDisequilibrium(inputAlign, windowSize, ldType, testSite, null,
                isAccumulateResults, numAccumulateIntervals, possibleSiteList, treatment);
        siteListLD.run();

        rowCount = siteListLD.getRowCount();
        columnCount = expected.getColumnCount();
        assertEquals("Site List LD Result Should have same number Columns as Expected", columnCount, siteListLD.getColumnCount());

        for (long i = 0; i < rowCount; i++) {
            Object[] actualRow = siteListLD.getRow(i);
            int site1 = (Integer) actualRow[2];
            int site2 = (Integer) actualRow[8];
            assertEquals("Site List LD RSqr Differ from Full LD: site1: " + site1 + " site2: " + site2, allLD.getRSqr(site1, site2), siteListLD.getRSqr(site1, site2), 0.0);
            assertEquals("Site List LD DPrime Differ from Full LD: site1: " + site1 + " site2: " + site2, allLD.getDPrime(site1, site2), siteListLD.getDPrime(site1, site2), 0.0);
            // LCJ 3: assertEquals("Site List LD PVal Differ from Full LD: site1: " + site1 + " site2: " + site2, allLD.getPVal(site1, site2), siteListLD.getPVal(site1, site2), 0.0);
            assertEquals("Site List LD SampleSize Differ from Full LD: site1: " + site1 + " site2: " + site2, allLD.getSampleSize(site1, site2), siteListLD.getSampleSize(site1, site2), 0.0);
        }
        siteListLD = null;
    }
}
