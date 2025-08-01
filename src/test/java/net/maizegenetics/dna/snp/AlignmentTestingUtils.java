/*
 * AlignmentTestingUtils
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.score.AlleleDepth;
import net.maizegenetics.dna.snp.score.AlleleDepthBuilder;
import net.maizegenetics.util.BitSet;

import java.util.Random;

import static net.maizegenetics.dna.snp.GenotypeTable.UNKNOWN_ALLELE;
import static org.junit.Assert.assertEquals;

/**
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class AlignmentTestingUtils {

    private AlignmentTestingUtils() {
        // utility
    }

    public static void alignmentsEqual(GenotypeTable align1, GenotypeTable align2) {
        alignmentsEqual(align1, align2, false);
    }

    public static void alignmentsEqual(GenotypeTable align1, GenotypeTable align2, boolean allowDifferentTaxaSort) {

        System.out.println("AlignmentTestingUtils: alignmentsEqual: Alignment 1: " + align1.getClass().getName() + ", " + align1.genotypeMatrix().getClass().getName() + ", " + align1.taxa().getClass().getName() + ", " + align1.positions().getClass().getName());
        System.out.println("AlignmentTestingUtils: alignmentsEqual: Alignment 2: " + align2.getClass().getName() + ", " + align2.genotypeMatrix().getClass().getName() + ", " + align2.taxa().getClass().getName() + ", " + align2.positions().getClass().getName());

        int siteCount = align1.numberOfSites();
        assertEquals("Expected Site Count: " + siteCount + " Actual: " + align2.numberOfSites(), siteCount, align2.numberOfSites());

        int taxaCount = align1.numberOfTaxa();
        assertEquals("Expected Taxa Count: " + taxaCount + " Actual: " + align2.numberOfTaxa(), taxaCount, align2.numberOfTaxa());

        int[] t1to2 = new int[align2.numberOfTaxa()];
        boolean sameTaxaSort = true;
        for (int t1 = 0; t1 < align2.numberOfTaxa(); t1++) {
            if (allowDifferentTaxaSort) {
                t1to2[t1] = align2.taxa().indexOf(align1.taxa().get(t1));
            } else {
                t1to2[t1] = t1;
            }
            if (t1to2[t1] != t1) {
                sameTaxaSort = false;
            }
        }
        System.out.println("Sort of taxa in same order:" + sameTaxaSort);

        for (int t = 0; t < taxaCount; t++) {
            assertEquals("Taxa Name: ", align1.taxaName(t), align2.taxaName(t1to2[t]));
            assertEquals("Taxa Object: ", align1.taxa().get(t), align2.taxa().get(t1to2[t]));
        }

        int[] physicalPositions = align1.physicalPositions();
        assertEquals("Number of Physical Positions: ", siteCount, physicalPositions.length);

        Chromosome[] loci1 = align1.chromosomes();
        Chromosome[] loci2 = align2.chromosomes();
        assertEquals("Number of Loci: ", loci1.length, loci2.length);
        int[] loci1Offsets = align1.chromosomesOffsets();
        int[] loci2Offsets = align2.chromosomesOffsets();
        assertEquals("Number of Loci Offsets: ", loci1Offsets.length, loci2Offsets.length);
        for (int i = 0; i < loci1.length; i++) {
            int end = align1.numberOfSites();
            if (i < loci1.length - 1) {
                end = loci1Offsets[i + 1];
            }
            int expectedLength = end - loci1Offsets[i];
            assertEquals("Site Length of Chromosome: ", expectedLength, align1.chromosomeSiteCount(loci1[i]));
            assertEquals("Site Length of Chromosome: ", align1.chromosomeSiteCount(loci1[i]), align2.chromosomeSiteCount(loci2[i]));
        }

        for (int t = 0; t < taxaCount; t++) {

            for (int s = 0; s < siteCount; s++) {

                byte first = align1.genotype(t, s);
                byte second = align2.genotype(t1to2[t], s);
                if ((first == (byte) 0xFE) || (first == (byte) 0xEF) || (first == (byte) 0xEE)) {
                    System.out.println("Z at taxa: " + t + " site: " + s + " value: " + first);
                } else {
                    if (first != second) {
                        byte reversed = (byte) ((first << 4) | (first >>> 4));
                        assertEquals("Base at taxa: " + t + " site: " + s + " doesn't match.", reversed, second);
                    }

                    byte[] firstArray = align1.genotypeArray(t, s);
                    byte[] secondArray = align2.genotypeArray(t1to2[t], s);
                    if (((firstArray[0] != secondArray[0]) || (firstArray[1] != secondArray[1]))
                            && ((firstArray[1] != secondArray[0]) || (firstArray[0] != secondArray[1]))) {
                        throw new IllegalStateException("Base Array at taxa: " + t + " site: " + s + " doesn't match.  first array: " + firstArray[0] + ":" + firstArray[1] + " second array: " + secondArray[0] + ":" + secondArray[1]);
                    }

                    String firstStr = align1.genotypeAsString(t, s);
                    String secondStr = align2.genotypeAsString(t1to2[t], s);
                    assertEquals("Base as String at taxa: " + t + " site: " + s + " doesn't match.", firstStr, secondStr);

                    assertEquals("Is Heterozygous: taxa: " + t + " site: " + s, align1.isHeterozygous(t, s), align2.isHeterozygous(t1to2[t], s));
                }

            }

        }

        for (int s = 0; s < siteCount; s++) {

            // allelesSortedByFrequency()
            int[][] align1SortedFreq = align1.allelesSortedByFrequency(s);
            int[][] align2SortedFreq = align2.allelesSortedByFrequency(s);
            assertEquals("getAllelesSortedByFrequency() length differ at Site: " + s, align1SortedFreq[0].length, align2SortedFreq[0].length);
            for (int i = 0, n = align1SortedFreq[0].length; i < n; i++) {
                assertEquals("getAllelesSortedByFrequency() values differ at Site: " + s, align1SortedFreq[0][i], align2SortedFreq[0][i]);
                assertEquals("getAllelesSortedByFrequency() counts differ at Site: " + s, align1SortedFreq[1][i], align2SortedFreq[1][i]);
            }

            assertEquals("Physical Position at Site: " + s, align1.chromosomalPosition(s), align2.chromosomalPosition(s));
            String snpID1 = align1.siteName(s);
            String snpID2 = align2.siteName(s);
            assertEquals("SNP ID at Site: " + s, snpID1, snpID2);
            Chromosome Chromosome1 = align1.chromosome(s);
            Chromosome Chromosome2 = align2.chromosome(s);
            assertEquals("Chromosome at Site: " + s, Chromosome1.getName(), Chromosome2.getName());
            assertEquals("Site at Physical Position: " + s, align1.siteOfPhysicalPosition(physicalPositions[s], Chromosome1), align2.siteOfPhysicalPosition(physicalPositions[s], Chromosome2));
            assertEquals("Site at Physical Position with SNP ID: " + s, align1.siteOfPhysicalPosition(physicalPositions[s], Chromosome1, snpID1), align2.siteOfPhysicalPosition(physicalPositions[s], Chromosome2, snpID2));
            assertEquals("Major Allele at Site: " + s, align1.majorAllele(s), align2.majorAllele(s));
            assertEquals("Minor Allele at Site: " + s, align1.minorAllele(s), align2.minorAllele(s));
            assertEquals("Major Allele Freq. at Site: " + s, align1.majorAlleleFrequency(s), align2.majorAlleleFrequency(s), 0.0001);
            assertEquals("Minor Allele Freq. at Site: " + s, align1.minorAlleleFrequency(s), align2.minorAlleleFrequency(s), 0.0001);
            assertEquals("Major Allele Count at Site: " + s, align1.majorAlleleCount(s), align2.majorAlleleCount(s), 0.0001);
            assertEquals("Minor Allele Count at Site: " + s, align1.minorAlleleCount(s), align2.minorAlleleCount(s), 0.0001);
            assertEquals("Heterozygous Count at Site: " + s, align1.heterozygousCount(s), align2.heterozygousCount(s));

            // alleles()
            byte[] align1Alleles = align1.alleles(s);
            byte[] align2Alleles = align2.alleles(s);
            assertEquals("getAlleles() length differ at Site: " + s, align1Alleles.length, align2Alleles.length);
            for (int i = 0, n = align1Alleles.length; i < n; i++) {
                assertEquals("getAlleles() values differ at Site: " + s, align1Alleles[i], align2Alleles[i]);
            }

        }

        for (int t = 0; t < taxaCount; t++) {
            byte[] align1Row = align1.genotypeAllSites(t);
            byte[] align2Row = align2.genotypeAllSites(t1to2[t]);
            for (int s = 0; s < siteCount; s++) {
                byte first = align1Row[s];
                byte second = align2Row[s];
                if (first != second) {
                    byte reversed = (byte) ((first << 4) | (first >>> 4));
                    assertEquals("genotypeAllSites: taxon: " + t + " site: " + s + " doesn't match.", reversed, second);
                }
            }
        }

        for (int s = 0; s < siteCount; s++) {
            byte[] align1Row = align1.genotypeAllTaxa(s);
            byte[] align2Row = align2.genotypeAllTaxa(s);
            for (int t = 0; t < taxaCount; t++) {
                byte first = align1Row[t];
                byte second = align2Row[t1to2[t]];
                if (first != second) {
                    byte reversed = (byte) ((first << 4) | (first >>> 4));
                    assertEquals("genotypeAllTaxa: taxon: " + t + " site: " + s + " doesn't match.", reversed, second);
                }
            }
        }

        byte[] align1Range = align1.genotypeRange(0, 0, siteCount);
        byte[] align2Range = align2.genotypeRange(t1to2[0], 0, siteCount);
        for (int i = 0; i < siteCount; i++) {
            byte first = align1Range[i];
            byte second = align2Range[i];
            if (first != second) {
                byte reversed = (byte) ((first << 4) | (first >>> 4));
                assertEquals("getBaseRange: taxon: 0 site: " + i + " doesn't match.", reversed, second);
            }
        }

        depthEqual(align1, align2);

    }

    public static void alignmentsTaxaBitSetsEqual(GenotypeTable align1, GenotypeTable align2) {

        int taxaCount = align1.numberOfTaxa();
        assertEquals("Expected Taxa Count: " + taxaCount + " Actual: " + align2.numberOfTaxa(), taxaCount, align2.numberOfTaxa());

        for (int t = 0; t < taxaCount; t++) {
            for (WHICH_ALLELE a : new WHICH_ALLELE[]{WHICH_ALLELE.Major, WHICH_ALLELE.Minor}) {
                BitSet bits1 = align1.allelePresenceForAllSites(t, a);
                BitSet bits2 = align2.allelePresenceForAllSites(t, a);
                assertEquals("Bit Cardinality Not Equal: Taxa: " + t + " Allele: " + a, bits1.cardinality(), bits2.cardinality());
                int numLongs = bits1.getNumWords();
                assertEquals("Number of Bit Words Not Equal: Taxa: " + t + " Allele: " + a, numLongs, bits2.getNumWords());
                for (int w = 0; w < numLongs; w++) {
                    assertEquals("Bits Not Equal: Taxa: " + t + " Allele: " + a, bits1.getBit(w), bits2.getBit(w));
                }
            }
        }

    }

    public static void alignmentsSiteBitSetsEqual(GenotypeTable align1, GenotypeTable align2) {

        int siteCount = align1.numberOfSites();
        assertEquals("Expected Site Count: " + siteCount + " Actual: " + align2.numberOfSites(), siteCount, align2.numberOfSites());

        for (int s = 0; s < siteCount; s++) {
            for (WHICH_ALLELE a : new WHICH_ALLELE[]{WHICH_ALLELE.Major, WHICH_ALLELE.Minor}) {
                BitSet bits1 = align1.allelePresenceForAllTaxa(s, a);
                BitSet bits2 = align2.allelePresenceForAllTaxa(s, a);
                assertEquals("Bit Cardinality Not Equal: Site: " + s + " Allele: " + a, bits1.cardinality(), bits2.cardinality());
                int numLongs = bits1.getNumWords();
                assertEquals("Number of Bit Words Not Equal: Site: " + s + " Allele: " + a, numLongs, bits2.getNumWords());
                for (int w = 0; w < numLongs; w++) {
                    assertEquals("Bits Not Equal: Site: " + s + " Allele: " + a, bits1.getBit(w), bits2.getBit(w));
                }
            }
        }

    }

    public static void alignmentsScopeEqual(GenotypeTable.ALLELE_SORT_TYPE scope, GenotypeTable align1, GenotypeTable align2) {

        int siteCount = align1.numberOfSites();
        assertEquals("Expected Site Count: " + siteCount + " Actual: " + align2.numberOfSites(), siteCount, align2.numberOfSites());

        // allelePresenceForAllTaxaBySortType
        for (int s = 0; s < siteCount; s++) {
            for (WHICH_ALLELE a : new WHICH_ALLELE[]{WHICH_ALLELE.Major, WHICH_ALLELE.Minor, WHICH_ALLELE.Unknown}) {
                BitSet bits1 = align1.allelePresenceForAllTaxa(s, a);
                BitSet bits2 = align2.allelePresenceForAllTaxa(s, a);
                assertEquals("Bit Cardinality Not Equal: Site: " + s + " Allele: " + a, bits1.cardinality(), bits2.cardinality());
                int numLongs = bits1.getNumWords();
                assertEquals("Number of Bit Words Not Equal: Site: " + s + " Allele: " + a, numLongs, bits2.getNumWords());
                for (int w = 0; w < numLongs; w++) {
                    assertEquals("Bits Not Equal: Site: " + s + " Allele: " + a, bits1.getBit(w), bits2.getBit(w));
                }
            }
        }

        // allelesBySortType
        for (int s = 0; s < siteCount; s++) {
            byte[] align1Alleles = align1.allelesBySortType(scope, s);
            byte[] align2Alleles = align2.allelesBySortType(scope, s);
            assertEquals("getAllelesByScope() length differ at Site: " + s, align1Alleles.length, align2Alleles.length);
            for (int i = 0, n = align1Alleles.length; i < n; i++) {
                assertEquals("getAllelesByScope() values differ at Site: " + s, align1Alleles[i], align2Alleles[i]);
            }
        }

    }

    public static void depthEqual(GenotypeTable genotype1, GenotypeTable genotype2) {

        AlleleDepth depth1 = genotype1.depth();
        AlleleDepth depth2 = genotype2.depth();

        if (depth1 == null || depth2 == null) {
            return;
        }

        int numSites = depth1.numSites();
        int numTaxa = depth1.numTaxa();

        assertEquals("Depth number sites not equal: ", numSites, depth2.numSites());
        assertEquals("Depth number taxa not equal: ", numTaxa, depth2.numTaxa());

        for (int t = 0; t < numTaxa; t++) {
            for (int s = 0; s < numSites; s++) {
                assertEquals("Depth values not equal: ", depth1.depth(t, s), depth2.depth(t, s));
            }
        }

    }

    public static GenotypeTable createRandomDepthForGenotypeTable(GenotypeTable a, int randomSeed) {
        AlleleDepth randomAd = createRandomDepth(a, randomSeed);
        GenotypeTable aWithDepth = GenotypeTableBuilder.getInstance(a.genotypeMatrix(), a.positions(), a.taxa(), randomAd);
        return aWithDepth;
    }

    public static AlleleDepth createRandomDepth(GenotypeTable gt, int randomSeed) {
        AlleleDepthBuilder adb = AlleleDepthBuilder.getInstance(gt.numberOfTaxa(), gt.numberOfSites(), gt.taxa());
        Random r = new Random(randomSeed);
        for (int t = 0; t < gt.numberOfTaxa(); t++) {
            byte[][] depths = new byte[6][gt.numberOfSites()];
            for (int s = 0; s < gt.numberOfSites(); s++) {
                byte[] a = gt.genotypeArray(t, s);
                if (a[0] != UNKNOWN_ALLELE) {
                    depths[a[0]][s] = (byte) (r.nextInt(4) + 1);
                }
                if (a[1] != UNKNOWN_ALLELE) {
                    depths[a[1]][s] = (byte) (r.nextInt(4) + 1);
                }
            }
            adb.addTaxon(t, depths);
        }
        return adb.build();
    }

}
