/*
 *  AlleleFreqCacheTest
 * 
 *  Created on Sep 22, 2015
 */
package net.maizegenetics.dna.snp.genotypecall;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.util.Utils;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Terry Casstevens
 */
public class AlleleFreqCacheTest {

    public AlleleFreqCacheTest() {
    }

    /**
     * Test of getAllelesSortedByFrequency method, of class AlleleFreqCache.
     */
    @Test
    public void testGetAllelesSortedByFrequency() throws IOException {
        GenotypeTable inputAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
        String expectedFile = GeneralConstants.EXPECTED_RESULTS_DIR + "AlleleFreqForTutorialHapmap.txt";
        try (BufferedReader reader = Utils.getBufferedReader(expectedFile)) {
            int numSites = inputAlign.numberOfSites();
            for (int s = 0; s < numSites; s++) {
                String[] expectedAlleles = reader.readLine().trim().split("\t");
                int[][] observedAlleles = inputAlign.allelesSortedByFrequency(s);
                assertEquals("Site: " + s + " Expected number alleles: " + (expectedAlleles.length / 2) + " Observed: " + observedAlleles[0].length, (expectedAlleles.length / 2), observedAlleles[0].length);
                for (int i = 0; i < observedAlleles[0].length; i++) {
                    assertEquals("Site: " + s + " Expected Allele: " + expectedAlleles[i * 2] + " Observed: " + observedAlleles[0][i], Integer.parseInt(expectedAlleles[i * 2]), observedAlleles[0][i]);
                    assertEquals("Site: " + s + " Expected Allele Count: " + expectedAlleles[i * 2 + 1] + " Observed: " + observedAlleles[1][i], Integer.parseInt(expectedAlleles[i * 2 + 1]), observedAlleles[1][i]);
                }
            }
        }
    }

    public static void main(String[] args) {
        GenotypeTable inputAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
        String outputFile = GeneralConstants.EXPECTED_RESULTS_DIR + "AlleleFreqForTutorialHapmap.txt";
        try (BufferedWriter writer = Utils.getBufferedWriter(outputFile)) {
            int numSites = inputAlign.numberOfSites();
            for (int s = 0; s < numSites; s++) {
                int[][] alleles = inputAlign.allelesSortedByFrequency(s);
                for (int i = 0; i < alleles[0].length; i++) {
                    if (i != 0) {
                        writer.write("\t");
                    }
                    writer.write(String.valueOf(alleles[0][i]));
                    writer.write("\t");
                    writer.write(String.valueOf(alleles[1][i]));
                }
                writer.write("\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
