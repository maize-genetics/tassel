package net.maizegenetics.matrixalgebra.Matrix;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;

import net.maizegenetics.matrixalgebra.decomposition.EigenvalueDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.SingularValueDecomposition;

import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

public class DoubleMatrixTest {
	DoubleMatrix A, B;
	double mintol = Double.MIN_NORMAL * 100; //maximum difference tolerated between two doubles
	double relaxedTol = 1e-12;
	
	@Before
	public void setUp() throws Exception {
		
	}

	@Test
	public void testEJML() {
		System.out.println("Testing ejml");
		DoubleMatrixFactory.setDefault(DoubleMatrixFactory.FactoryType.ejml);
		intializeTestMatrices();
		runTests();
	}
	
	@Test
	public void testColt() {
		System.out.println("Testing Colt");
		DoubleMatrixFactory.setDefault(DoubleMatrixFactory.FactoryType.colt);
		intializeTestMatrices();
		runTests();
	}

	@Ignore
	@Test
	public void testBlas() {
		System.out.println("Testing BLAS");
		DoubleMatrixFactory.setDefault(DoubleMatrixFactory.FactoryType.blas);
		intializeTestMatrices();
		runTests();
	}
	
	public void runTests() {
		testGetsSets();
		testTranspose();
		testMult();
		testCrossProduct();
		testSVD();
		testEigenvalues();
		testInverse();
		testConcatenate();
		testTcrossproduct();
		testXtXGM();
		testSelection();
		testPlusMinus();
		testScalarAddMult();
		testRow();
		testColumn();
	}
	
	public void intializeTestMatrices() {
		A = DoubleMatrixFactory.DEFAULT.make(4, 5,new double[]{0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9}, true);
		B = DoubleMatrixFactory.DEFAULT.make(5, 3, new double[]{0,1,2,3,4,5,6,7,0,1,2,3,4,5,6}, false);
	}
	
	public void testGetsSets() {
		assertEquals(0, A.get(0, 0), mintol);
		assertEquals(5, A.get(1, 1), mintol);
		assertEquals(9, A.get(3, 4), mintol);
		assertEquals(0, A.getChecked(0, 0), mintol);
		assertEquals(0, A.getChecked(2, 2), mintol);
		
		assertEquals(1, B.get(0, 1), mintol);
		assertEquals(7, B.get(2, 1), mintol);
		assertEquals(3, B.get(3, 2), mintol);
		
		double tmp = A.get(3, 4);
		A.set(3, 4, 20);
		assertEquals(20, A.get(3,4), mintol);
		A.set(3, 4, tmp);
		
		tmp = A.get(2, 1);
		A.set(2, 1, 20);
		assertEquals(20, A.get(2, 1), mintol);
		A.set(2, 1, tmp);
	}
	
	public void testTranspose() {
		DoubleMatrix AT = A.transpose();
		assertEquals(0, AT.get(0, 0), mintol);
		assertEquals(5, AT.get(1, 1), mintol);
		assertEquals(9, AT.get(4, 3), mintol);
		
		DoubleMatrix BT = B.transpose();
		assertEquals(1, BT.get(1, 0), mintol);
		assertEquals(7, BT.get(1, 2), mintol);
		assertEquals(3, BT.get(2, 3), mintol);
	}

	public void testMult() {
		DoubleMatrix AB = DoubleMatrixFactory.DEFAULT.make(new double[][]{{86,106,62},{100,125,78},{54,74,94},{68,93,110}});
		DoubleMatrix C = A.mult(B);
		compareMatrices(AB, C, false);
		
		C = A.multadd(B, null, 1, 0, false, false);
		compareMatrices(AB, C, false);
	}
	
	public void testCrossProduct() {
		DoubleMatrix ATA = DoubleMatrixFactory.DEFAULT.make(new double[][]{{14,38,12,26,50},{38,126,84,82,170},{12,84,146,48,120},{26,82,48,54,110},{50,170,120,110,230}});
		DoubleMatrix C = A.crossproduct();
		compareMatrices(ATA, C, false);
		
		DoubleMatrix D = DoubleMatrixFactory.DEFAULT.make(new double[][]{{0,4,8},{1,5,9},{2,6,0},{3,7,1}});
		DoubleMatrix ATD = DoubleMatrixFactory.DEFAULT.make(new double[][]{{14,38,12},{38,126,84},{12,84,146},{26,82,48},{50,170,120}});
		C = A.crossproduct(D);
		compareMatrices(ATD, C, false);
	}
	
	public void testSVD() {
		SingularValueDecomposition svd = B.getSingularValueDecomposition();
		double[][] uArray = new double[][]{{-0.119594996694745,-0.243028865682223,-0.759364212490157},{-0.485856848446252,-0.318246381010037,0.112527169127685},{-0.5671276459022,0.806169962745293,-0.168689728260741},{-0.241682280611914,-0.268101370791494,-0.468733751950877},{-0.607944132363421,-0.343318886119308,0.403157629666964}};
		double[][] vArray = new double[][]{{-0.535872570953034,0.39036677868278,0.748635001720474},{-0.6797103902171,0.326517808497542,-0.656795178244238},{-0.500833678111909,-0.860813509904899,0.0903644179647362}};
		double[] sArray = new double[]{14.058930498008,5.74050826183117,0.626927546076496};
		double[] svector = svd.getSingularValues();
		DoubleMatrix uMatrix = svd.getU(false);
		DoubleMatrix vMatrix = svd.getV(false);
		
		//compare singular values
		int[] index = indexSort(svector, false);
		for (int i = 0; i < 3; i++) assertEquals(sArray[i], svector[index[i]], relaxedTol);
		
		//compare u and v
		assertEquals(5, uMatrix.numberOfRows());
		assertEquals(3, uMatrix.numberOfColumns());
		assertEquals(3, vMatrix.numberOfRows());
		assertEquals(3, vMatrix.numberOfColumns());
		
		for (int r = 0; r < 5; r++) {
			for (int c = 0; c < 3; c++) {
				assertEquals(Math.abs(uArray[r][c]), Math.abs(uMatrix.get(r, index[c])), relaxedTol);
			}
		}
		
		for (int r = 0; r < 3; r++) {
			for (int c = 0; c < 3; c++) {
				assertEquals(Math.abs(vArray[r][c]), Math.abs(vMatrix.get(r, index[c])), relaxedTol);
			}
		}

	}
	
	public void testEigenvalues() {
		double[] eigenval = new double[]{197.653526747819,32.9534351041519,0.393038148029496};
		double[][] eigenvec = new double[][]{{-0.535872570953034,-0.39036677868278,0.748635001720475},{-0.6797103902171,-0.326517808497542,-0.656795178244237},{-0.500833678111908,0.860813509904899,0.0903644179647361}};
		DoubleMatrix BTB = B.crossproduct();
		EigenvalueDecomposition ed = BTB.getEigenvalueDecomposition();
		
		//test eigenvalues
		double[] ev = ed.getEigenvalues();
		int[] index = indexSort(ev, false);
		for (int i = 0; i < 3; i++) assertEquals(eigenval[i], ev[index[i]], relaxedTol);
		
		//test eigenvectors
		DoubleMatrix vectors = ed.getEigenvectors();
		for (int r = 0; r < 3; r++) {
			for (int c = 0; c < 3; c++) {
				assertEquals(Math.abs(eigenvec[r][c]), Math.abs(vectors.get(r, index[c])), relaxedTol);
			}
		}

	}
	
	public void testInverse() {
		DoubleMatrix BTB = B.crossproduct();
		DoubleMatrix inv = BTB.inverse();
		DoubleMatrix expected = DoubleMatrixFactory.DEFAULT.make(new double[][]{{1.43203124999999,-1.2453125,0.16328125},{-1.2453125,1.103125,-0.1578125},{0.163281250000001,-0.157812500000001,0.0445312500000002}});
		compareMatrices(expected, inv, true);
	}
	
	public void testConcatenate() {
		double[] v1 = new double[]{1,1,1,1,1,1};
		double[] v2 = new double[]{2,2,2,2,2,2};
		DoubleMatrix m1 = DoubleMatrixFactory.DEFAULT.make(1,6,v1,true);
		DoubleMatrix m2 = DoubleMatrixFactory.DEFAULT.make(1,6,v2,true);
		DoubleMatrix m3 = m1.concatenate(m2, true);
		assertEquals(6, m3.numberOfColumns());
		assertEquals(2, m3.numberOfRows());
		assertEquals(1, m3.get(0, 3), mintol);
		assertEquals(2, m3.get(1, 3), mintol);

		m1 = DoubleMatrixFactory.DEFAULT.make(6,1,v1,true);
		m2 = DoubleMatrixFactory.DEFAULT.make(6,1,v2,true);
		m3 = m1.concatenate(m2, false);
		assertEquals(2, m3.numberOfColumns());
		assertEquals(6, m3.numberOfRows());
		assertEquals(1, m3.get(3, 0), mintol);
		assertEquals(2, m3.get(3, 1), mintol);
	}
	
	public void testTcrossproduct() {
		DoubleMatrix AAT = DoubleMatrixFactory.DEFAULT.make(new double[][]{{120,140,80,100},{140,165,100,125},{80,100,120,140},{100,125,140,165}});
		DoubleMatrix C = A.tcrossproduct();
		compareMatrices(AAT, C, false);
		
		DoubleMatrix D = DoubleMatrixFactory.DEFAULT.make(new double[][]{{0,2,4,6,8},{1,3,5,7,9}});
		DoubleMatrix ADT = DoubleMatrixFactory.DEFAULT.make(new double[][]{{100,120},{120,145},{100,120},{120,145}});
		C = A.tcrossproduct(D);
		compareMatrices(ADT, C, false);
	}
	
	public void testXtXGM() {
		DoubleMatrix XTX = DoubleMatrixFactory.DEFAULT.make(new double[][]{{62,76,42},{76,95,58},{42,58,74}});
		DoubleMatrix G = DoubleMatrixFactory.DEFAULT.make(new double[][]{{1.43203124999999,-1.2453125,0.16328125},{-1.2453125,1.103125,-0.1578125},{0.163281250000001,-0.157812500000001,0.0445312500000002}});
		DoubleMatrix M = DoubleMatrixFactory.DEFAULT.make(new double[][]{{0.350000000000004,-0.050000000000001,-2.66453525910038e-15,-0.449999999999997,0.149999999999997},{-0.0499999999999984,0.649999999999996,-6.55031584528842e-15,-0.15,-0.450000000000006},{-7.32747196252603e-15,-6.32827124036339e-15,-1.50990331349021e-14,-6.99440505513849e-15,-5.99520433297585e-15},{-0.449999999999997,-0.150000000000002,-4.44089209850063e-15,0.650000000000001,-0.0500000000000038},{0.150000000000001,-0.450000000000005,-7.105427357601e-15,-0.0500000000000013,0.349999999999993}});
		
		DoubleMatrix[] xtxgm = B.getXtXGM();
		
		compareMatrices(XTX, xtxgm[0], true);
		compareMatrices(G, xtxgm[1], true);
		compareMatrices(M, xtxgm[2], true);
	}
	
	public void testSelection() {
		DoubleMatrix S = A.getSelection(new int[]{0,1,2}, new int[]{1,2,3});
		DoubleMatrix C = DoubleMatrixFactory.DEFAULT.make(new double[][]{{4,8,2},{5,9,3},{6,0,4}});
		compareMatrices(C, S, false);
		
		S = A.getSelection(new int[]{1,3}, new int[]{2,3,4});
		C = DoubleMatrixFactory.DEFAULT.make(new double[][]{{9,3,7},{1,5,9}});
		compareMatrices(C, S, false);

		S = A.getSelection(new int[]{0,0,1}, new int[]{3,1,0});
		C = DoubleMatrixFactory.DEFAULT.make(new double[][]{{2,4,0},{2,4,0},{3,5,1}});
		compareMatrices(C, S, false);
	}
	
	public void testPlusMinus() {
		DoubleMatrix A2 = DoubleMatrixFactory.DEFAULT.make(new double[][]{{0,8,16,4,12},{2,10,18,6,14},{4,12,0,8,16},{6,14,2,10,18}});
		
		DoubleMatrix C = A.plus(A);
		compareMatrices(A2, C, false);
		
		DoubleMatrix D = C.minus(A);
		compareMatrices(A, D, false);
		
		C = A.copy();
		C.plusEquals(A);
		compareMatrices(A2, C, false);
		
		C.minusEquals(A);
		compareMatrices(A, C, false);
	}
	
	public void testScalarAddMult() {
		DoubleMatrix A2 = DoubleMatrixFactory.DEFAULT.make(new double[][]{{0,8,16,4,12},{2,10,18,6,14},{4,12,0,8,16},{6,14,2,10,18}});
		DoubleMatrix Aplus3 = DoubleMatrixFactory.DEFAULT.make(new double[][]{{3,7,11,5,9},{4,8,12,6,10},{5,9,3,7,11},{6,10,4,8,12}});
		DoubleMatrix Aminus2 = DoubleMatrixFactory.DEFAULT.make(new double[][]{{-2,2,6,0,4},{-1,3,7,1,5},{0,4,-2,2,6},{1,5,-1,3,7}});
		
		DoubleMatrix C = A.scalarAdd(3);
		compareMatrices(Aplus3, C, false);
		
		C = A.scalarAdd(-2);
		compareMatrices(Aminus2, C, false);
		
		C = A.scalarMult(2);
		compareMatrices(A2, C, false);
		
		C = A.copy();
		C.scalarAddEquals(3);
		compareMatrices(Aplus3, C, false);
		
		C = A.copy();
		C.scalarAddEquals(-2);
		compareMatrices(Aminus2, C, false);
		
		C = A.copy();
		C.scalarMultEquals(2);
		compareMatrices(A2, C, false);
		
	}
	
	public void testRow() {
		double[][] testMatrix = new double[][]{{0,8,16,4,12},{2,10,18,6,14},{4,12,0,8,16},{6,14,2,10,18}};
		DoubleMatrix A = DoubleMatrixFactory.DEFAULT.make(testMatrix);
		for (int i = 0; i < 4; i++) {
			System.out.println("i = " + i);
			DoubleMatrix B = A.row(i);
			double[] thisRow = new double[5];
			for (int j = 0; j < 5; j++) {
				thisRow[j] = B.get(j, 0);
			}
			assertArrayEquals(testMatrix[i], thisRow, mintol);
		}
	}
	
	public void testColumn() {
		double[][] testMatrix = new double[][]{{0,8,16,4,12},{2,10,18,6,14},{4,12,0,8,16},{6,14,2,10,18}};
		DoubleMatrix A = DoubleMatrixFactory.DEFAULT.make(testMatrix);
		for (int i = 0; i < 5; i++) {
			System.out.println("i = " + i);
			DoubleMatrix B = A.column(i);
			double[] thisCol = new double[4];
			double[] testCol = new double[4];
			for (int j = 0; j < 4; j++) {
				thisCol[j] = B.get(j, 0);
				testCol[j] = testMatrix[j][i];
			}
			assertArrayEquals(testCol, thisCol, mintol);
		}
	}
	
	private void compareMatrices(DoubleMatrix m1, DoubleMatrix m2, boolean relaxed) {
		double delta;
		if (relaxed) delta = relaxedTol;
		else delta = mintol;
		assertEquals(m1.numberOfRows(), m2.numberOfRows());
		assertEquals(m1.numberOfColumns(), m2.numberOfColumns());
		for (int r = 0; r < m1.numberOfRows(); r++) {
			for (int c = 0; c < m1.numberOfColumns(); c++) {
				assertEquals(m1.get(r, c), m2.get(r, c), delta);
			}
		}
	}
	
	private void printDoubleArray(double[] array) {
		for (double d : array) System.out.format("%3.5f ", d);
		System.out.println();
	}
	
	private int[] indexSort(double[] array, boolean ascending) {
		class indexpair implements Comparable<indexpair> {
			double val;
			int ndx;
			indexpair(double val, int ndx) {
				this.val = val;
				this.ndx = ndx;
			}
			@Override
			public int compareTo(indexpair o) {
				if (val > o.val) return 1;
				if (val < o.val) return -1;
				return 0;
			}
			
		}
		
		LinkedList<indexpair> pairlist = new LinkedList<indexpair>();
		int n = array.length;
		for (int i = 0; i < n; i++) pairlist.add(new indexpair(array[i], i));
		Collections.sort(pairlist);
		int[] index = new int[n];
		if (ascending) {
			for (int i = 0; i < n; i++) {
				index[i] = pairlist.get(i).ndx;
			}
		} else {
			for (int i = 0; i < n; i++) {
				index[i] = pairlist.get(n - i - 1).ndx;
			}
		}

		return index;
	}
}
