package net.maizegenetics.analysis.numericaltransform;

import org.junit.Test;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Random;

public class kNearestNeighborsTest {

    public double[][] data() {
        double[][] genotype = new double[10][10];
        double[] options = new double[]{1.0, 0.5, 0.0, -0.0 / 0.0};
        Random ranGen = new Random();
        ranGen.setSeed(123456789);
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                int randValue = ranGen.nextInt(4);
                genotype[i][j] = options[randValue];
            }

        }
        return genotype;
    }

    @Test
    public void testEuclidean() {
        kNearestNeighborsTest knntest = new kNearestNeighborsTest();
        double[][] genotype = knntest.data();
        int k = 3;
        boolean isManhattan = false;
        boolean isCosine = false;
        double[][] imputedData = kNearestNeighbors.impute(genotype, k, isManhattan, isCosine);

        ArrayList<Double> imputedValues = new ArrayList<>();

        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                if (Double.isNaN(genotype[i][j])) {
                    imputedValues.add(imputedData[i][j]);
                }

            }
        }
        double[] expectedResults = new double[]{0.5, 0.5, 0.5, 0.0, 0.5, 0.25, 0.0, 0.25, 0.75, 0.5, 0.75, 0.75, 0.0, 0.25, 0.0, 1.0, 0.5};
        double[] actualResults = new double[17];
        for (int l = 0; l < 17; l++) {
            actualResults[l] = imputedValues.get(l);
        }
        assertArrayEquals(expectedResults, actualResults, 0.001);
    }

    @Test
    public void testManhattan() {
        kNearestNeighborsTest knnTest = new kNearestNeighborsTest();
        double[][] genotype = knnTest.data();
        int k = 3;
        boolean isManhattan = true;
        boolean isCosine = false;
        double[][] imputedData = kNearestNeighbors.impute(genotype, k, isManhattan, isCosine);

        ArrayList<Double> imputedValues = new ArrayList<>();

        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                if (Double.isNaN(genotype[i][j])) {
                    imputedValues.add(imputedData[i][j]);
                }

            }
        }
        double[] expectedResults = new double[]{0.5, 0.5, 0.25, 1.0, 0.5, 0.25, 0.0, 0.25, 0.25, 0.5, 0.75, 0.75, 0.0, 0.25, 0.5, 0.42857142857142855, 0.5};
        double[] actualResults = new double[17];
        for (int l = 0; l < 17; l++) {
            actualResults[l] = imputedValues.get(l);
        }
        assertArrayEquals(expectedResults, actualResults, 0.001);

    }

    @Test
    public void testCosine() {
        kNearestNeighborsTest knnTest = new kNearestNeighborsTest();
        double[][] genotype = knnTest.data();
        int k = 3;
        boolean isManhattan = false;
        boolean isCosine = true;
        double[][] imputedData = kNearestNeighbors.impute(genotype, k, isManhattan, isCosine);

        ArrayList<Double> imputedValues = new ArrayList<>();

        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                if (Double.isNaN(genotype[i][j])) {
                    imputedValues.add(imputedData[i][j]);
                }

            }
        }
        double[] expectedResults = new double[]{0.25, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.5, 0.25, 0.0, 0.0, 0.25, 0.5, 0.5, 0.25};
        double[] actualResults = new double[17];
        for (int l = 0; l < 17; l++) {
            actualResults[l] = imputedValues.get(l);
        }
        assertArrayEquals(expectedResults, actualResults, 0.001);

    }

}
