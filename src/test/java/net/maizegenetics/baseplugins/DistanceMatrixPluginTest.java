/*
 * DistanceMatrixPluginTest
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportTestUtils;
import net.maizegenetics.util.TableReportUtils;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import org.junit.*;

import java.io.File;
import java.io.IOException;
import net.maizegenetics.analysis.distance.DistanceMatrixPlugin;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

/**
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class DistanceMatrixPluginTest {

    /**
     * Test of getDistanceMatrix method, of class DistanceMatrixPlugin.
     */
    @Test
    public void testGetDistanceMatrixAlignment() throws IOException {
        System.out.println("Testing Distance Matrix with Filtering Taxa...");

        String expectedResultsFilename = GeneralConstants.EXPECTED_RESULTS_DIR + "DistanceMatrixPluginTest.hmp.txt";
        GenotypeTable expectedAlign = ImportUtils.readFromHapmap(expectedResultsFilename, null);

        File input = new File(TutorialConstants.HAPMAP_FILENAME);
        System.out.println("   Input File: " + input.getCanonicalPath());
        if (!input.exists()) {
            fail("Input File: " + TutorialConstants.HAPMAP_FILENAME + " doesn't exist.");
        }
        GenotypeTable inputAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
        System.out.println("   Input Alignment: Taxa Count: " + inputAlign.numberOfTaxa() + "  Site Count: " + inputAlign.numberOfSites());

        TaxaList origTaxa = inputAlign.taxa();
        Taxon[] newIds = new Taxon[origTaxa.numberOfTaxa() / 3 + 1];
        int count = 0;
        for (int i = 0; i < origTaxa.numberOfTaxa(); i = i + 3) {
            newIds[count++] = origTaxa.get(i);
        }
        TaxaList taxa = new TaxaListBuilder().addAll(newIds).build();
        GenotypeTable filterTaxa = FilterGenotypeTable.getInstance(inputAlign, taxa);

        DistanceMatrix matrix1 = DistanceMatrixPlugin.getDistanceMatrix(expectedAlign);
        DistanceMatrix matrix2 = DistanceMatrixPlugin.getDistanceMatrix(filterTaxa);

        TableReport expected = TableReportUtils.readDelimitedTableReport(GeneralConstants.EXPECTED_RESULTS_DIR + "DistanceMatrixPluginTest.txt", "\t");

        TableReportTestUtils.compareTableReports(expected, matrix1, 0.0000001);
        TableReportTestUtils.compareTableReports(expected, matrix2, 0.0000001);

    }

    @Test
    public void testComputeHetBitDistances_3args() {
        System.out.println("Testing computeHetBitDistances()...");
        int taxon1 = 0;
        int taxon2 = 1;
        String[] seqs = {"AAA", "AAA"};
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
        double result = IBSDistanceMatrix.computeHetBitDistances(theTBA, taxon1, taxon2)[0];
        assertEquals(0.0, result, 0.00001);

        seqs[1] = "CAA";
        g = GenotypeCallTableBuilder.getInstance(seqs.length, seqs[0].length()).setBases(seqs).build();
        theTBA = GenotypeTableBuilder.getInstance(g, pLB.build(), tLB.build());
        result = IBSDistanceMatrix.computeHetBitDistances(theTBA, taxon1, taxon2)[0];
        assertEquals(0.3333333333, result, 0.00001);

        seqs[1] = "CGT";
        g = GenotypeCallTableBuilder.getInstance(seqs.length, seqs[0].length()).setBases(seqs).build();
        theTBA = GenotypeTableBuilder.getInstance(g, pLB.build(), tLB.build());
        result = IBSDistanceMatrix.computeHetBitDistances(theTBA, taxon1, taxon2)[0];
        assertEquals(1.0, result, 0.00001);

        seqs[1] = "AAR";
        g = GenotypeCallTableBuilder.getInstance(seqs.length, seqs[0].length()).setBases(seqs).build();
        theTBA = GenotypeTableBuilder.getInstance(g, pLB.build(), tLB.build());
        result = IBSDistanceMatrix.computeHetBitDistances(theTBA, taxon1, taxon2)[0];
        assertEquals(1.0 / 6.0, result, 0.00001);

        seqs[0] = "RMR";
        seqs[1] = "RMR";
        g = GenotypeCallTableBuilder.getInstance(seqs.length, seqs[0].length()).setBases(seqs).build();
        theTBA = GenotypeTableBuilder.getInstance(g, pLB.build(), tLB.build());
        result = IBSDistanceMatrix.computeHetBitDistances(theTBA, taxon1, taxon2)[0];
        assertEquals(0.5, result, 0.00001);
    }
}
