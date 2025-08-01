package net.maizegenetics.analysis.numericaltransform;

import org.junit.Test;

import java.util.ArrayList;
import java.util.Random;

import static org.junit.Assert.*;

public class ImputationByMeanTest {

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
    public void test() {
        ImputationByMeanTest meanTest = new ImputationByMeanTest();
        double[][] genotype = meanTest.data();
        double[][] imputedData = ImputationByMean.impute(genotype);
        ArrayList<Double> imputedValues = new ArrayList<>();

        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                if (Double.isNaN(genotype[i][j])) {
                    imputedValues.add(imputedData[i][j]);
                }

            }
        }
        double[] expectedResults = new double[]{0.625,0.6875,0.3889,0.4286,0.6875,0.2143,0.375,0.3889,0.375,0.4286,0.4286,0.4286,0.2143,0.4286,0.625,0.4286,0.2143};
        double[] actualResults = new double[17];
        for (int l = 0; l < 17; l++) {
            actualResults[l] = imputedValues.get(l);
        }
        assertArrayEquals(expectedResults, actualResults, 0.001);

    }

}
