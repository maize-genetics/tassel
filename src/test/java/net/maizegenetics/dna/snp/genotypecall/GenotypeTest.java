/*
 *  GenotypeTest
 */
package net.maizegenetics.dna.snp.genotypecall;

import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import junit.framework.Assert;

import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;

import org.junit.*;

import static org.junit.Assert.assertEquals;

/**
 *
 * @author Terry Casstevens
 */
public class GenotypeTest {

    public GenotypeTest() {
    }

    @BeforeClass
    public static void setUpClass() {
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    @Test
    public void testInstanceCopy() {
        int numTaxa = 10000;
        int numSites = 100000;
        SuperByteMatrix temp = createByteMatrixRandomNucleotideData(numTaxa, numSites);
        GenotypeCallTable genotype = new NucleotideGenotypeCallTable(temp, false);
        System.out.println("genotype class: " + genotype.getClass().getName());
        long previous = System.nanoTime();
        GenotypeCallTable genotypeCopy = GenotypeCallTableBuilder.getInstanceCopy(genotype).build();
        double timeToCopy = (double) (System.nanoTime() - previous) / 1_000_000_000.0;
        System.out.println("GenotypeTest: testInstanceCopy: taxa: " + numTaxa + " sites: " + numSites + " time: " + timeToCopy);
        compareGenotypeCallTables(genotype, genotypeCopy);
    }

    /**
     * Tests the thread safety of Genotypes
     */
    @Test
    public void threadSafeTest() {
        System.out.println("GenotypeTest: threadSafeTest...");
        int numTaxa = 1000;
        int numSites = 100000;
        SuperByteMatrix temp = createByteMatrixRandomNucleotideData(numTaxa, numSites);
        SuperByteMatrix temp1 = SuperByteMatrixBuilder.getInstanceCopy(temp);
        GenotypeCallTable genotype = new NucleotideGenotypeCallTable(temp, false);
        GenotypeCallTable genotype1 = new ByteGenotypeCallTable(temp1, false, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
        ExecutorService pool = Executors.newFixedThreadPool(5);
        pool.execute(new Thread(new VerifyThreadSafe(genotype, genotype1, 0, -1)));
        pool.execute(new Thread(new VerifyThreadSafe(genotype, genotype1, numSites / 10, -1)));
        pool.execute(new Thread(new VerifyThreadSafe(genotype, genotype1, 0, -1)));
        pool.execute(new Thread(new VerifyThreadSafe(genotype, genotype1, 0, numSites / 2)));
        pool.execute(new Thread(new VerifyThreadSafe(genotype, genotype1, 0, -1)));
        try {
            pool.shutdown();
            if (!pool.awaitTermination(60, TimeUnit.SECONDS)) {
                throw new IllegalStateException("GenotypeTest: threadSafeTest: processing threads timed out.");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static class VerifyThreadSafe implements Runnable {

        private final GenotypeCallTable myGenotype;
        private final GenotypeCallTable myGenotype1;
        private final int myStart;
        private final int myEnd;

        public VerifyThreadSafe(GenotypeCallTable genotype, GenotypeCallTable genotype1, int startSite, int endSite) {
            myGenotype = genotype;
            myGenotype1 = genotype1;
            myStart = startSite;
            myEnd = endSite;
        }

        @Override
        public void run() {
            verifyResults(myGenotype, myGenotype1, myStart, myEnd);
        }
    }

    private static void verifyResults(GenotypeCallTable genotype, GenotypeCallTable genotype1, int startSite, int endSite) {

        int numSites = genotype.numberOfSites();
        int numTaxa = genotype.numberOfTaxa();
        System.out.println("GenotypeTest: verifyResults: sites: " + numSites + "  taxa: " + numTaxa + "  start: " + startSite + "  end: " + endSite);
        if (endSite < 0) {
            endSite = numSites;
        }

        long previous;

        long otherTime = 0;
        long firstTime = 0;
        for (int s = startSite; s < endSite; s++) {
            previous = System.nanoTime();
            int[][] alleleFreq = genotype.allelesSortedByFrequency(s);
            firstTime += System.nanoTime() - previous;
            previous = System.nanoTime();
            int[][] alleleFreq1 = genotype1.allelesSortedByFrequency(s);
            otherTime += System.nanoTime() - previous;
            int numAlleles = alleleFreq[0].length;
            for (int a = 0; a < numAlleles; a++) {
                if (alleleFreq[0][a] != alleleFreq1[0][a]) {
                    System.out.println("GenotypeTest: verifyResults: First Genotype: " + genotype.getClass().getName());
                    printAlleles(alleleFreq);
                    System.out.println("GenotypeTest: verifyResults: Second Genotype: " + genotype1.getClass().getName());
                    printAlleles(alleleFreq1);
                    assertEquals("GenotypeTest: verifyResults: Not equal: site: " + s + " allele: " + a, alleleFreq[0][a], alleleFreq1[0][a]);
                }

                if (alleleFreq[1][a] != alleleFreq1[1][a]) {
                    System.out.println("GenotypeTest: verifyResults: First Genotype: " + genotype.getClass().getName());
                    printAlleles(alleleFreq);
                    System.out.println("GenotypeTest: verifyResults: Second Genotype: " + genotype1.getClass().getName());
                    printAlleles(alleleFreq1);
                    assertEquals("GenotypeTest: verifyResults: Not equal: site: " + s + " allele: " + a, alleleFreq[1][a], alleleFreq1[1][a]);
                }
            }
        }

        System.out.println("GenotypeTest: verifyResults: First Time: " + genotype.getClass().getName() + ": " + ((double) firstTime / 1_000_000_000.0) + "  start: " + startSite + "  end: " + endSite);
        System.out.println("GenotypeTest: verifyResults: Second Time: " + genotype1.getClass().getName() + ": " + ((double) otherTime / 1_000_000_000.0) + "  start: " + startSite + "  end: " + endSite);
    }

    private static void printAlleles(int[][] alleles) {
        int numAlleles = alleles[0].length;
        for (int a = 0; a < numAlleles; a++) {
            System.out.println(alleles[0][a] + ": " + alleles[1][a]);
        }
    }

    @Ignore
    @Test
    public void iterateAllBases() {
        int numTaxa = 10000;
        int numSites = 100000;
        SuperByteMatrix temp = createByteMatrixRandomNucleotideData(numTaxa, numSites);
        GenotypeCallTable genotype = new NucleotideGenotypeCallTable(temp, false);

        long previous = System.nanoTime();
        for (int t = 0; t < numTaxa; t++) {
            for (int s = 0; s < numSites; s++) {
                genotype.genotype(t, s);
            }
        }
        double timeToIterate = (double) (System.nanoTime() - previous) / 1_000_000_000.0;
        System.out.println("GenotypeTest: iterateAllBases: taxa: " + numTaxa + " sites: " + numSites + " time: " + timeToIterate);
        Assert.assertTrue("GenotypeTest: iterateAllBases: " + timeToIterate + " sec. longer than expected: 0.05", timeToIterate < 0.05);

        previous = System.nanoTime();
        genotype.transposeData(false);
        timeToIterate = (double) (System.nanoTime() - previous) / 1_000_000_000.0;
        System.out.println("GenotypeTest: time to transpose: taxa: " + numTaxa + " sites: " + numSites + " time: " + timeToIterate);
        Assert.assertTrue("GenotypeTest: time to transpose: " + timeToIterate + " sec. longer than expected: 1.0", timeToIterate < 1.0);

        previous = System.nanoTime();
        for (int s = 0; s < numSites; s++) {
            for (int t = 0; t < numTaxa; t++) {
                genotype.genotype(t, s);
            }
        }
        timeToIterate = (double) (System.nanoTime() - previous) / 1_000_000_000.0;
        System.out.println("GenotypeTest: iterateAllBases transposed: taxa: " + numTaxa + " sites: " + numSites + " time: " + timeToIterate);
        Assert.assertTrue("GenotypeTest: iterateAllBases transposed: " + timeToIterate + " sec. longer than expected: 7.0", timeToIterate < 7.0);
    }

    @Ignore
    @Test
    public void getMajorAlleleForAllSites() {
        int numTaxa = 10000;
        int numSites = 100000;
        SuperByteMatrix temp = createByteMatrixRandomNucleotideData(numTaxa, numSites);
        GenotypeCallTable genotype = new NucleotideGenotypeCallTable(temp, false);
        long previous = System.nanoTime();
        genotype.majorAlleleForAllSites();
        double timeToIterate = (double) (System.nanoTime() - previous) / 1_000_000_000.0;
        System.out.println("GenotypeTest: getMajorAlleleForAllSites: " + numTaxa + " sites: " + numSites + " time: " + timeToIterate);
        Assert.assertTrue("GenotypeTest: getMajorAlleleForAllSites: " + timeToIterate + " sec. longer than expected: 10.0", timeToIterate < 10.0);
    }

    @Test
    public void testSetBaseRangeForTaxon() {
        int numTaxa = 10000;
        int numSites = 100000;
        SuperByteMatrix temp = createByteMatrixRandomNucleotideData(numTaxa, numSites);
        GenotypeCallTable genotype = new NucleotideGenotypeCallTable(temp, false);
        long previous = System.nanoTime();
        GenotypeCallTableBuilder builder = GenotypeCallTableBuilder.getInstance(numTaxa, numSites);
        for (int t = 0; t < numTaxa; t++) {
            byte[] current = genotype.genotypeAllSites(t);
            builder.setBaseRangeForTaxon(t, 0, current);
        }
        double timeLapsed = (double) (System.nanoTime() - previous) / 1_000_000_000.0;
        System.out.println("GenotypeTest: testSetBaseRangeForTaxon: " + numTaxa + " sites: " + numSites + " time: " + timeLapsed);
        GenotypeCallTable newGenotype = builder.build();
        compareGenotypeCallTables(genotype, newGenotype);
        Assert.assertTrue("GenotypeTest: testSetBaseRangeForTaxon: " + timeLapsed + " sec. longer than expected: 5.0", timeLapsed < 5.0);
    }

    public static void compareGenotypeCallTables(GenotypeCallTable table1, GenotypeCallTable table2) {
        int numTaxa = table1.numberOfTaxa();
        int numSites = table1.numberOfSites();
        assertEquals("Number of Taxa Differ", numTaxa, table2.numberOfTaxa());
        assertEquals("Number of Sites Differ", numSites, table2.numberOfSites());
        for (int t = 0; t < numTaxa; t++) {
            for (int s = 0; s < numSites; s++) {
                assertEquals("Genotype Call differ", table1.genotype(t, s), table2.genotype(t, s));
            }
        }
    }

    /**
     * Sets genotype bases with randomly set values using only valid nucleotide
     * encodings.
     *
     * @param numTaxa number of taxa
     * @param numSites number of sites
     *
     * @return byte matrix
     */
    public static SuperByteMatrix createByteMatrixRandomNucleotideData(int numTaxa, int numSites) {
        SuperByteMatrix matrix = SuperByteMatrixBuilder.getInstance(numTaxa, numSites);
        Random random = new Random();
        for (int t = 0; t < numTaxa; t++) {
            for (int s = 0; s < numSites; s++) {
                int value = random.nextInt(127);
                if ((value >>> 4) > 5) {
                    value |= 0xF0;
                }
                if ((value & 0xF) > 5) {
                    value |= 0xF;
                }
                matrix.set(t, s, (byte) value);
            }
        }
        return matrix;
    }

    public static void setRandomNucleotideData(GenotypeCallTableBuilder matrix) {
        int numTaxa = matrix.getTaxaCount();
        int numSites = matrix.getSiteCount();
        Random random = new Random();
        for (int t = 0; t < numTaxa; t++) {
            for (int s = 0; s < numSites; s++) {
                int value = random.nextInt(127);
                if ((value >>> 4) > 5) {
                    value |= 0xF0;
                }
                if ((value & 0xF) > 5) {
                    value |= 0xF;
                }
                matrix.setBase(t, s, (byte) value);
            }
        }
    }
}
