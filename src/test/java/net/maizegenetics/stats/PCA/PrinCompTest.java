package net.maizegenetics.stats.PCA;

import static org.junit.Assert.*;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.stats.PCA.PrinComp.PC_TYPE;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class PrinCompTest {
	DoubleMatrix myData;
	PrinComp myPCA;
	double tol = 1e-12;
	double[] testEigenval;
	double[][] testEigenvec;
	double[][] testPC;
	
	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void test() {	
		useDataSet1();
		myPCA = new PrinComp(myData, PC_TYPE.cov);
		System.out.printf("DoubleMatrixFactory type is %s\n", DoubleMatrixFactory.DEFAULT.getType());
		System.out.println("Eigenvalues:");
		for (double val:myPCA.getEigenValues()) System.out.print(val + " ");
		System.out.println();
		
		System.out.println("Eigenvectors:");
		System.out.println(myPCA.getEigenVectors());
		
		System.out.println("PC's:");
		System.out.println(myPCA.getPrincipalComponents());
		
		assertTrue("Testing eigenvalues for test set 1...", testEigenvalues());
		assertTrue("Testing eigenvectors for test set 1...", testEigenvectors());
		assertTrue("Testing PC's for test set 1...", testPC());
		
		useDataSet2();
		myPCA = new PrinComp(myData, PC_TYPE.corr);
		System.out.println("Eigenvalues:");
		for (double val:myPCA.getEigenValues()) System.out.print(val + " ");
		System.out.println();
		
		System.out.println("Eigenvectors:");
		System.out.println(myPCA.getEigenVectors());
		
		System.out.println("PC's:");
		System.out.println(myPCA.getPrincipalComponents());
		
		assertTrue("Testing eigenvalues for test set 2...", testEigenvalues());
		assertTrue("Testing eigenvectors  for test set 2...", testEigenvectors());
		assertTrue("Testing PC's for test set 2...", testPC());
	}

	public void useDataSet1() {
		double[][] data = new double[][]{{1,2,3,2,3,4},{3,2,3,4,5,5},{4,3,4,5,6,5},{4,8,7,3,4,5}};
		myData = DoubleMatrixFactory.DEFAULT.make(data);
		
		//test values from R: prcomp()
		testEigenval = new double[]{1.265282e+01, 4.716365e+00, 4.747914e-02, 2.501067e-31};
		
		testEigenvec = new double[][]{
			{-0.26040315,  0.49115810,  0.2994106, -0.34767622},
			{-0.80318313, -0.13582600,  0.1136607, -0.36391773},
			{-0.53113010, -0.03333491, -0.4290802,  0.64855531},
			{0.01164988,  0.59364919, -0.2433303,  0.03151932},
			{0.01164988,  0.59364919, -0.2433303,  0.03151932},
			{-0.06867629,  0.18537818,  0.7713011,  0.56927516}
		};
		testPC = new double[][]{
				{2.6068470, -2.6229333, -0.10986215,  3.885781e-16},
				{2.0639639,  0.9193579,  0.28693917, -3.053113e-16},
				{0.4925473,  2.4286535, -0.21573033, -2.775558e-17},
				{-5.1633582, -0.7250780,  0.03865331,  2.775558e-17}
		};
	}
	
	public void useDataSet2() { //use corr with this data set
		double[][] data = new double[][]{
				{9, 5, 7, 9, 4, 7, 1},
				{2, 1, 1, 1, 1, 7, 6},
				{7, 6, 9, 8, 4, 2, 7},
				{2, 6, 8, 7, 9, 8, 4},
				{3, 7, 4, 7, 4, 4, 4}
		};
		myData = DoubleMatrixFactory.DEFAULT.make(data);
		
		//test values from R: prcomp()
		testEigenval = new double[]{3.537200e+00, 1.688439e+00, 1.278101e+00, 4.962601e-01, 3.970544e-32};
		
		testEigenvec = new double[][] {
				{0.3123845, -0.19660085, -0.64877339,  0.3224638, -0.11369759},
				{0.4622523, -0.07987140,  0.26331346, -0.5401578,  0.45748534},
				{0.4798065, -0.06589826,  0.12565918,  0.5645608,  0.51010183},
				{0.5247024, -0.01732956, -0.12095154, -0.1184232, -0.43041095},
				{0.3413618,  0.37243242,  0.51229863,  0.1914193, -0.50807637},
				{-0.1426675,  0.72445284, -0.06051857,  0.2737018,  0.26207596},
				{-0.2152440, -0.53552586,  0.46176532,  0.3999245, -0.07414278}	
		};
		testPC = new double[][]{
				{1.2309943, 0.83516558, -1.7306678,  0.08630054, -1.249001e-16},
				{-3.2826861, 0.01445715, -0.2166971,  0.24060123, -1.228184e-15},
				{1.0824049, -1.94998024,  0.2254010,  0.53510557, -2.053913e-15},
				{0.8133263, 1.45798609,  1.3790800,  0.36330020, -3.042705e-15},
				{0.1559606, -0.35762858,  0.3428839, -1.22530754,  7.088080e-15}
		};
	}
	
	public boolean testEigenvalues() {
		int n = testEigenval.length;
		double[] calculatedEigenval = myPCA.getEigenValues();
		boolean pass = true;
		for (int i = 0; i < n; i++) {
			if (Math.abs(testEigenval[i] - calculatedEigenval[i]) > Math.max(tol, testEigenval[i] * 1e-6)) {
				pass = false;
				break;
			}
		}
		return pass;
	}
	
	public boolean testEigenvectors() {
		int ncol = testEigenvec[0].length;
		int nrow = testEigenvec.length;
		double[] calculatedEigenval = myPCA.getEigenValues();
		DoubleMatrix calculatedEigenvec = myPCA.getEigenVectors();
		
		boolean pass = true;
		for (int c = 0; c < ncol; c++) {
			//if the eigenvalue ~= 0, the eigenvectors may not be the same but that does not matter
			if (Math.abs(calculatedEigenval[c]) > tol && pass) {
				for (int r = 0; r < nrow; r++) {
					if (Math.abs(testEigenvec[r][c] - calculatedEigenvec.get(r, c)) > Math.max(tol, Math.abs(testEigenvec[r][c]) * 1e-6)) {
						pass = false;
						break;
					}
				}
				
				//calculated eigenvector may be of opposite sign
				if (!pass) {
					pass = true;
					for (int r = 0; r < nrow; r++) {
						if (Math.abs(testEigenvec[r][c] + calculatedEigenvec.get(r, c)) > Math.max(tol, Math.abs(testEigenvec[r][c]) * 1e-6)) {
							pass = false;
							break;
						}
					}
				}
			}
		}
		
		return pass;
	}
	
	public boolean testPC() {
		int ncol = testPC[0].length;
		int nrow = testPC.length;
		double[] calculatedEigenval = myPCA.getEigenValues();
		DoubleMatrix calculatedPC = myPCA.getPrincipalComponents();
		
		boolean pass = true;
		for (int c = 0; c < ncol; c++) {
			//if the eigenvalue ~= 0, the eigenvectors may not be the same but that does not matter
			if (Math.abs(calculatedEigenval[c]) > tol && pass) {
				for (int r = 0; r < nrow; r++) {
					if (Math.abs(testPC[r][c] - calculatedPC.get(r, c)) > Math.max(tol, Math.abs(testPC[r][c]) * 1e-6)) {
						pass = false;
						break;
					}
				}
				
				//calculated eigenvector may be of opposite sign
				if (!pass) {
					pass = true;
					for (int r = 0; r < nrow; r++) {
						if (Math.abs(testPC[r][c] + calculatedPC.get(r, c)) > Math.max(tol, Math.abs(testPC[r][c]) * 1e-6)) {
							pass = false;
							break;
						}
					}
				}
			}
		}
		
		return pass;
	}
}
