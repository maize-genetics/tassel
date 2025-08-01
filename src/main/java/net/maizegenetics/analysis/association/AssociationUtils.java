package net.maizegenetics.analysis.association;

import java.util.Arrays;
import java.util.List;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffectUtils;
import net.maizegenetics.util.BitSet;

public class AssociationUtils {
	
    public static DoubleMatrix createFixedEffectsArray(List<PhenotypeAttribute> factorList, List<PhenotypeAttribute> covariateList, BitSet missing, int numberOfObservations) {
    	int numberOfFactors;
    	int numberOfCovariates;
    	if (factorList == null) numberOfFactors = 0;
    	else numberOfFactors = factorList.size();
    	if (covariateList == null) numberOfCovariates = 0;
    	else numberOfCovariates = covariateList.size();
    	int numberOfEffects = 1 + numberOfFactors + numberOfCovariates;
    	DoubleMatrix[][] theMatrices = new DoubleMatrix[1][numberOfEffects];
    	
    	//the mean
    	int count = 0;
    	theMatrices[0][count++] = DoubleMatrixFactory.DEFAULT.make(numberOfObservations, 1, 1.0);

    	for (PhenotypeAttribute factorAttribute : factorList) {
    		CategoricalAttribute attr = (CategoricalAttribute) factorAttribute;
    		String[] nonMissingFactorLevels = AssociationUtils.getNonMissingValues(attr.allLabels(), missing);
    		int[] levels = ModelEffectUtils.getIntegerLevels(nonMissingFactorLevels, null);
    		FactorModelEffect fme = new FactorModelEffect(levels, true);
    		theMatrices[0][count++] = fme.getX();
    	}

    	for (PhenotypeAttribute covariateAttribute : covariateList) {
    		NumericAttribute attr = (NumericAttribute) covariateAttribute;
    		double[] nonMissingValues = AssociationUtils.getNonMissingDoubles(attr.floatValues(), missing);
    		theMatrices[0][count++] = DoubleMatrixFactory.DEFAULT.make(numberOfObservations, 1, nonMissingValues);
    	}
    	
    	if (theMatrices[0].length == 1) return theMatrices[0][0];
    	return DoubleMatrixFactory.DEFAULT.compose(theMatrices);
    }

	public static double[] getNonMissingDoubles(double[] allData, BitSet missing) {
		int originalLength = allData.length;
		int resultLength = originalLength - (int) missing.cardinality();
		double[] result = new double[resultLength];
		int resultCount = 0;
		for (int i = 0; i < originalLength;i++) {
			if (!missing.fastGet(i)) result[resultCount++] = allData[i];
		}
		return result;
	}

	public static double[] getNonMissingDoubles(float[] allData, BitSet missing) {
		int originalLength = allData.length;
		int resultLength = originalLength - (int) missing.cardinality();
		double[] result = new double[resultLength];
		int resultCount = 0;
		for (int i = 0; i < originalLength;i++) {
			if (!missing.fastGet(i)) result[resultCount++] = allData[i];
		}
		return result;
	}

	public static byte[] getNonMissingBytes(byte[] allData, BitSet missing) {
		int originalLength = allData.length;
		int resultLength = originalLength - (int) missing.cardinality();
		byte[] result = new byte[resultLength];
		int resultCount = 0;
		for (int i = 0; i < originalLength;i++) {
			if (!missing.fastGet(i)) result[resultCount++] = allData[i];
		}
		return result;
	}

	public static <T> T[] getNonMissingValues(T[] allData, BitSet missing) {
		int originalLength = allData.length;
		int resultLength = originalLength - (int) missing.cardinality();
		T[] result = Arrays.copyOf(allData, resultLength);
		int resultCount = 0;
		for (int i = 0; i < originalLength;i++) {
			if (!missing.fastGet(i)) result[resultCount++] = allData[i];
		}
		return result;
	}

	public static DoubleMatrix getNonMissingValues(DoubleMatrix allData, BitSet missing) {
		int originalLength = allData.numberOfRows();
		int numberMissing = (int) missing.cardinality();
		if (numberMissing == 0) return allData;
		if (originalLength > 1) {
			int[] select = new int[originalLength - numberMissing];
			int notMissingCount = 0;
			for (int i = 0; i < originalLength; i++) {
				if (!missing.fastGet(i)) select[notMissingCount++] = i;
			}
			return allData.getSelection(select, null);
		} else {
			originalLength = allData.numberOfColumns();
			int[] select = new int[originalLength - numberMissing];
			int notMissingCount = 0;
			for (int i = 0; i < originalLength; i++) {
				if (!missing.fastGet(i)) select[notMissingCount++] = i;
			}
			return allData.getSelection(null, select);
		}
	}

	public static boolean isMonomorphic(float[] probs) {
		float first = Float.NaN;
		int n = probs.length;
		for (int i = 0; i < n; i++) {
			if (!Float.isNaN(probs[i])) {
				if (Float.isNaN(first)) first = probs[i];
				else if (probs[i] != first) return false;
			}
		}
		return true;
	}
	
	public static double[] convertFloatArrayToDouble(float[] floatArray) {
		int n = floatArray.length;
		double[] doubleArray = new double[n];
		for (int i = 0; i < n; i++) doubleArray[i] = floatArray[i];
		return doubleArray;
	}
	
	public static float[] convertDoubleArrayToFloat(double[] doubleArray) {
		int n = doubleArray.length;
		float[] floatArray = new float[n];
		for (int i = 0; i < n; i++) floatArray[i] = (float) doubleArray[i];
		return floatArray;
	}
	
}
