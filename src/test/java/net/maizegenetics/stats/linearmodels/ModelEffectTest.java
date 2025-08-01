package net.maizegenetics.stats.linearmodels;

import static org.junit.Assert.*;

import java.util.Arrays;

import org.junit.Test;

public class ModelEffectTest {

    @Test
    public void testIntLevels() {
        int[] origLevels = new int[] {0,1,2,3,4,0,0,1,1,2,2,3,3,4,4};
        int[] reset = ModelEffectUtils.getIntegerLevels(origLevels);
        assertArrayEquals(origLevels, reset);
        
        int[] modifiedLevels = new int[] {0,1,2,4,6,0,0,1,1,2,2,4,4,6,6};
        reset = ModelEffectUtils.getIntegerLevels(modifiedLevels);
        assertArrayEquals(origLevels, reset);
    }
    
    @Test
    public void testSubsetCovariate() {
        double[] origLevels = new double[] {0,1,2,3,4,0,0,1,1,2,2,3,3,4,4};
        int[] subset = new int[]{0,2,2,4,6,7,8,10,14};
        CovariateModelEffect cme = new CovariateModelEffect(origLevels);
        int nsubset = subset.length;
        double[] subsetLevels = Arrays.stream(subset).mapToDouble(i -> origLevels[i]).toArray();
        CovariateModelEffect subcme = (CovariateModelEffect) cme.getSubSample(subset);
        assertArrayEquals(subsetLevels, subcme.getCovariate(), 1e-10);
    }
    
    @Test
    public void testSubsetFactor() {
        String[] factor = new String[]{"A","A","B","B","C","C","C","C","D","D","D","D","E","E","E"};
        FactorModelEffect fme = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(factor), false);
        int[] subset = new int[]{0,2,2,4,6,7,8,10,14};
        FactorModelEffect subfme = (FactorModelEffect) fme.getSubSample(subset);
        int[] expectedLevels = new int[]{0,1,1,2,2,2,3,3,4};
        assertArrayEquals(expectedLevels, subfme.getLevels());
    }

}
