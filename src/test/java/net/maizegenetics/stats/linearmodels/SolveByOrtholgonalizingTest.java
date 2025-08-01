package net.maizegenetics.stats.linearmodels;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

public class SolveByOrtholgonalizingTest {
	static double[] x1, x2, x3, y;
	
	@BeforeClass
	public static void beforeTest() {
		x1 = new double[]{1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2};
		x2 = new double[]{1.151,-0.037,0.155,1.192,1.411,0.268,-0.618,2.097,-0.65,1.3,1.277,3.544,0.946,3.564,1.205,3.267,2.013,4.555,2.939,2.571};
		x3 = new double[]{1.352,-2.207,-0.118,3.145,-1.14,5.457,1.589,0.303,-4.265,-0.879,-3.544,3.845,3.185,3.825,2.599,0.973,1.758,3.392,1.424,-1.636};
		y = new double[]{3.384,1.175,0.381,-3.952,0.856,7.744,-2.734,5.256,2.854,-1.167,0.762,9.54,4.031,16.509,3.838,0.832,5.839,7.596,10.753,2.363};
	}
	
	@Test
	public void testVectors() {
		List<double[]> covList = new ArrayList<>();
		List<double[]> dataList = new ArrayList<>();
		covList.add(x1);
		dataList.add(y);
		SolveByOrthogonalizing sbo = SolveByOrthogonalizing.getInstanceFromVectors(covList, dataList);
		SolveByOrthogonalizing.Marker myMarker = sbo.solveForR(null, x3);
		assertEquals(0.1642605, myMarker.vector1()[0], 1e-5);
		double F = SolveByOrthogonalizing.calculateFfromR2(myMarker.vector1()[0], 1, 17);
		assertEquals(3.341267, F, 1e-5);
		
		covList.add(x2);
		sbo = SolveByOrthogonalizing.getInstanceFromVectors(covList, dataList);
		myMarker = sbo.solveForR(null, x3);
		assertEquals(0.1121142, myMarker.vector1()[0], 1e-5);
		F = SolveByOrthogonalizing.calculateFfromR2(myMarker.vector1()[0], 1, 16);
		assertEquals(2.020335, F, 1e-5);
		
		covList = new ArrayList<>();
		sbo = SolveByOrthogonalizing.getInstanceFromVectors(covList, dataList);
		myMarker = sbo.solveForR(null, x3);
		assertEquals(0.214881, myMarker.vector1()[0], 1e-5);
		F = SolveByOrthogonalizing.calculateFfromR2(myMarker.vector1()[0], 1, 18);
		assertEquals(4.92646, F, 1e-5);

		covList = new ArrayList<>();
		sbo = SolveByOrthogonalizing.getInstanceFromVectors(covList, dataList);
		myMarker = sbo.solveForR(null, x2, x3);
		assertEquals(0.4105654, myMarker.vector1()[0], 1e-5);
		F = SolveByOrthogonalizing.calculateFfromR2(myMarker.vector1()[0], 2, 17);
		assertEquals(5.9206, F, 1e-5);

	}
	
	@Test
	public void testModels() {
		List<ModelEffect> covList = new ArrayList<>();
		List<double[]> dataList = new ArrayList<>();
		ModelEffect me1 = new CovariateModelEffect(x1);
		me1.setID("x1");
		covList.add(me1);
		dataList.add(y);
		SolveByOrthogonalizing sbo = SolveByOrthogonalizing.getInstanceFromModel(covList, dataList);
		SolveByOrthogonalizing.Marker myMarker = sbo.solveForR(null, x3);
		assertEquals(0.1642605, myMarker.vector1()[0], 1e-5);
		
		ModelEffect me2 = new CovariateModelEffect(x2);
		me2.setID("x2");
		covList.add(me2);
		dataList.add(y);
		sbo = SolveByOrthogonalizing.getInstanceFromModel(covList, dataList);
		myMarker = sbo.solveForR(null, x3);
		assertEquals(0.1121142, myMarker.vector1()[0], 1e-5);
		
		covList = new ArrayList<>();
		covList.add(me1);
		sbo = SolveByOrthogonalizing.getInstanceFromModel(covList, dataList);
		myMarker = sbo.solveForR(null, x2, x3);
		assertEquals(0.2362625, myMarker.vector1()[0], 1e-5);
		double F = SolveByOrthogonalizing.calculateFfromR2(myMarker.vector1()[0], 2, 16);
		assertEquals(2.474803, F, 1e-5);
	}

}
