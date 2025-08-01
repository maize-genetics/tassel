package net.maizegenetics.analysis.distance;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.dna.snp.GenotypeTable.GENOTYPE_TABLE_COMPONENT;
import net.maizegenetics.taxa.distance.DistanceMatrix;

import java.util.function.BiConsumer;
import java.util.stream.IntStream;

/**
 * @author Terry Casstevens
 * @author Zhiwu Zhang
 * @author Peter Bradbury
 *
 * @deprecated This utility has been replaced. See KinshipPlugin. Please discuss
 * with Terry if something needs to be salvaged from this class.
 */
public class Kinship {

    //scale the numeric matrix produced by the transform function or from probabilities which code phenotypes as {1,0.5,0}
    public static double MATRIX_MULTIPLIER = 2;

    public static enum KINSHIP_TYPE {

        Endelman
    };

    private Kinship() {
        // utility to create kinship matrix
    }

    /**
     * @deprecated Replaced by {@link EndelmanDistanceMatrix#getInstance(net.maizegenetics.dna.snp.GenotypeTable, int, net.maizegenetics.util.ProgressListener)
     * }
     */
    private static DistanceMatrix createKinship(GenotypeTable mar, KINSHIP_TYPE kinshipType, GENOTYPE_TABLE_COMPONENT dataType) {
        System.out.println("Starting Kinship.buildFromMarker().");
        long start = System.currentTimeMillis();
        DistanceMatrix result;
        if (kinshipType == KINSHIP_TYPE.Endelman) {
            result = buildFromMarker(mar, kinshipType, dataType);
        } else {
            throw new IllegalArgumentException("Kinship: createKinship: unknown kinship algorithm: " + kinshipType);
        }
        System.out.printf("Built Kinship in %d millisec.\n", System.currentTimeMillis() - start);
        return result;
    }

    private static DistanceMatrix buildFromMarker(GenotypeTable mar, KINSHIP_TYPE kinshipType, GENOTYPE_TABLE_COMPONENT dataType) {
        if (dataType == GENOTYPE_TABLE_COMPONENT.Genotype) {
            return calculateKinshipFromMarkers(mar);
        } else if (dataType == GENOTYPE_TABLE_COMPONENT.ReferenceProbability) {
            return calculateRelationshipKinshipFromReferenceProbability(mar);
        } else {
            throw new IllegalArgumentException("The supplied data type is not currently supported by the Kinship method: " + dataType);
        }

    }

    /**
     * Calculates a kinship matrix from genotypes using the method described in
     * Endelman and Jannink (2012) G3 2:1407-1413. It is best to impute missing
     * data before calculating. However, if data is missing it is replaced by
     * the allele average at that site.
     *
     */
    private static DistanceMatrix calculateKinshipFromMarkers(GenotypeTable mar) {
        //mar is the input genotype table
        byte missingAllele = GenotypeTable.UNKNOWN_ALLELE;

        // from Endelman and Jannink. 2012. G3 2:1407ff
        // A = WW'/[2*sumk(pk*qk)]
        // where W = centered genotype matrix (centered on marker mean value, marker coded as 2,1,0)
        // where marker is multi-allelic, leave one allele out to keep markers independent
        int ntaxa = mar.numberOfTaxa();
        int nsites = mar.numberOfSites();
        double[][] distance = new double[ntaxa][ntaxa];
        double sumpi = 0;

        //calculate WW' by summing ww' for each allele, where w is a column vector of centered allele counts {2,1,0}
        for (int s = 0; s < nsites; s++) {
            int[][] alleleFreq = mar.allelesSortedByFrequency(s);
            int nalleles = alleleFreq[0].length;
            int totalAlleleCount = mar.totalGametesNonMissingForSite(s);

            for (int a = 0; a < nalleles - 1; a++) {
                double pi = ((double) alleleFreq[1][a]) / ((double) totalAlleleCount);
                double pix2 = 2 * pi;
                sumpi += pi * (1 - pi);
                double[] scores = new double[ntaxa];

                for (int t = 0; t < ntaxa; t++) {
                    byte[] geno = GenotypeTableUtils.getDiploidValues(mar.genotype(t, s));
                    double thisScore = 0;
                    if (geno[0] != missingAllele) {
                        if (geno[0] == alleleFreq[0][a]) {
                            thisScore++;
                        }
                        if (geno[1] == alleleFreq[0][a]) {
                            thisScore++;
                        }
                        thisScore -= pix2;
                    }
                    scores[t] = thisScore;
                }

                for (int r = 0; r < ntaxa; r++) {
                    double rowval = scores[r];
                    distance[r][r] += rowval * rowval;
                    for (int c = r + 1; c < ntaxa; c++) {
                        distance[r][c] += rowval * scores[c];
                    }
                }
            }
        }

        double sumpk = 2 * sumpi;

        for (int r = 0; r < ntaxa; r++) {
            distance[r][r] /= sumpk;
            for (int c = r + 1; c < ntaxa; c++) {
                distance[r][c] = distance[c][r] = distance[r][c] / sumpk;
            }
        }

        return new DistanceMatrix(distance, mar.taxa());
    }

    /**
     * Calculates a kinship matrix from genotypes using the method described in
     * Endelman and Jannink (2012) G3 2:1407-1413. It is best to impute missing
     * data before calculating. However, if data is missing it is replaced by
     * the allele average at that site.
     *
     */
    private static DistanceMatrix calculateKinshipFromMarkersV2(GenotypeTable mar) {
        //mar is the input genotype table
        byte missingAllele = GenotypeTable.UNKNOWN_ALLELE;

        // from Endelman and Jannink. 2012. G3 2:1407ff
        // A = WW'/[2*sumk(pk*qk)]
        // where W = centered genotype matrix (centered on marker mean value, marker coded as 2,1,0)
        // where marker is multi-allelic, leave one allele out to keep markers independent
        int ntaxa = mar.numberOfTaxa();
        int nsites = mar.numberOfSites();
        double[][] distance = new double[ntaxa][ntaxa];
        double sumpq = IntStream.range(0, nsites).parallel().mapToDouble(s -> {
            int[][] alleleFreq = mar.allelesSortedByFrequency(s);
            int nalleles = alleleFreq[0].length;
            double totalAlleleCount = mar.totalGametesNonMissingForSite(s);
            double pq = 0;
            for (int a = 0; a < nalleles - 1; a++) {
                double p = alleleFreq[1][a] / totalAlleleCount;
                pq += p * (1 - p);
            }
            return pq;
        }).sum();

        //calculate WW' by summing ww' for each allele, where w is a column vector of centered allele counts {2,1,0}
        BiConsumer<double[][], Integer> siteDistance = (dist, site) -> {
            int s = site.intValue();
            int[][] alleleFreq = mar.allelesSortedByFrequency(s);
            int nalleles = alleleFreq[0].length;
            int totalAlleleCount = mar.totalGametesNonMissingForSite(s);

            for (int a = 0; a < nalleles - 1; a++) {
                double pi = ((double) alleleFreq[1][a]) / ((double) totalAlleleCount);
                double pix2 = 2 * pi;
                double[] scores = new double[ntaxa];

                for (int t = 0; t < ntaxa; t++) {
                    byte[] geno = GenotypeTableUtils.getDiploidValues(mar.genotype(t, s));
                    double thisScore = 0;
                    if (geno[0] != missingAllele) {
                        if (geno[0] == alleleFreq[0][a]) {
                            thisScore++;
                        }
                        if (geno[1] == alleleFreq[0][a]) {
                            thisScore++;
                        }
                        thisScore -= pix2;
                    }
                    scores[t] = thisScore;
                }

                for (int r = 0; r < ntaxa; r++) {
                    double rowval = scores[r];
                    dist[r][r] += rowval * rowval;
                    for (int c = r + 1; c < ntaxa; c++) {
                        dist[r][c] += rowval * scores[c];
                    }
                }
            }
        };

        BiConsumer<double[][], double[][]> mergeDistance = (d1, d2) -> {
            for (int r = 0; r < ntaxa; r++) {
                double[] row1 = d1[r];
                double[] row2 = d2[r];
                for (int c = 0; c < ntaxa; c++) {
                    row1[c] += row2[c];
                }
            }
        };

        distance = IntStream.range(0, nsites).boxed().parallel().collect(() -> new double[ntaxa][ntaxa], siteDistance, mergeDistance);

        double sumpk = 2 * sumpq;

        for (int r = 0; r < ntaxa; r++) {
            distance[r][r] /= sumpk;
            for (int c = r + 1; c < ntaxa; c++) {
                distance[r][c] = distance[c][r] = distance[r][c] / sumpk;
            }
        }

        return new DistanceMatrix(distance, mar.taxa());
    }

    private static DistanceMatrix calculateRelationshipKinshipFromReferenceProbability(GenotypeTable mar) {
        ReferenceProbability referenceP = mar.referenceProbability();

        //calculate the column averages and sumpq, center W
        int nrow = referenceP.numTaxa();
        int ncol = referenceP.numSites();
        double[][] W = new double[nrow][ncol];

        for (int r = 0; r < nrow; r++) {
            for (int c = 0; c < ncol; c++) {
                W[r][c] = referenceP.value(r, c) * MATRIX_MULTIPLIER;
            }
        }

        double sumpq = 0;
        for (int c = 0; c < ncol; c++) {
            double colTotal = 0;
            int colCount = 0;
            for (int r = 0; r < nrow; r++) {
                if (!Double.isNaN(W[r][c])) {
                    colTotal += W[r][c];
                    colCount++;
                }
            }

            double pi = colTotal / colCount / 2.0;
            double pix2 = pi * 2;
            sumpq += pi * (1 - pi);
            for (int r = 0; r < nrow; r++) {
                if (Double.isNaN(W[r][c])) {
                    W[r][c] = 0;
                } else {
                    W[r][c] -= pix2;
                }
            }
        }

        DoubleMatrix WWt = DoubleMatrixFactory.DEFAULT.make(W).tcrossproduct();

        double[][] scaledIBS = new double[nrow][nrow];
        for (int r = 0; r < nrow; r++) {
            for (int c = 0; c < nrow; c++) {
                scaledIBS[r][c] = WWt.get(r, c) / sumpq / 2;
            }
        }
        return new DistanceMatrix(scaledIBS, mar.taxa());
    }

    private static DistanceMatrix calculateRelationshipKinshipFromAlleleProbabilities() {
        //TODO implement
        return null;
    }

}
