/**
 * 
 */
package net.maizegenetics.stats.statistics;
import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * @author lcj34
 *
 */
public class FisherExactTest {

	@Test
	public void TestFisherExact() {

		int[][] argInts = new int[15][4];
		argInts[0] = new int[]{2, 3, 6, 4};
		argInts[1] = new int[]{2, 1, 3, 0};
		argInts[2] = new int[]{3, 0, 2, 1};
		argInts[3] = new int[]{1, 2, 0, 3};
		argInts[4] = new int[]{3, 1, 1, 3};
		argInts[5] = new int[]{1, 3, 3, 1};
		argInts[6] = new int[]{0, 1, 1, 0};
		argInts[7] = new int[]{1, 0, 0, 1};
		argInts[8] = new int[]{11, 0, 0, 6};
		argInts[9] = new int[]{10, 1, 1, 5};
		argInts[10] = new int[]{5, 6, 6, 0};
		argInts[11] = new int[]{9, 2, 2, 4};
		argInts[12] = new int[]{6, 5, 5, 1};
		argInts[13] = new int[]{8, 3, 3, 3};
		argInts[14] = new int[]{7, 4, 4, 2};

		// Store pValues from R using values above.
		// Values were obtained by running R-Studio for each of the
		// argInts values above as the x values (e.g. x = c(2,3,6,4) )
		// followed by the commands:
		//    y = matrix(x,2,2, byrow=TRUE)
		//    ft = fisher.test(y)
		//    ft$p  ( to get 2- tailed p-value)
		//    fisher.test(y, alternative='less') - for the leftTailedR value
		//    fisher.test(y, alternative='greater') - for the rightTailed R value

		double[] cumulativeR = new double[15];  // What is equivalent in R's Fisher Test?
		double[] twoTailedR = new double[15];
		twoTailedR[0] = 0.6083916;
		twoTailedR[1] = 1;
		twoTailedR[2] = 1;
		twoTailedR[3] = 1;
		twoTailedR[4] = 0.4857143;
		twoTailedR[5] = 0.4857143;
		twoTailedR[6] = 1;
		twoTailedR[7] = 1;
		twoTailedR[8] = 8.080155e-05;
		twoTailedR[9] = 0.005413704;
		twoTailedR[10] = 0.04274402;
		twoTailedR[11] = 0.1094053;
		twoTailedR[12] = 0.3333872;
		twoTailedR[13] = 0.6000323;
		twoTailedR[14] = 1;
		
		double[] rightTailedR = new double[15];
		rightTailedR[0] = 0.8998;
		rightTailedR[1] = 1;
		rightTailedR[2] = 0.5;
		rightTailedR[3] = 0.5;
		rightTailedR[4] = 0.2429;
		rightTailedR[5] = 0.9857;
		rightTailedR[6] = 1;
		rightTailedR[7] = 0.5;
		rightTailedR[8] = 8.08e-05;
		rightTailedR[9] = 0.005414;
		rightTailedR[10] = 1;
		rightTailedR[11] = 0.07207;
		rightTailedR[12] = 0.9627;
		rightTailedR[13] = 0.3387;
		rightTailedR[14] = 0.7387;

		double[] leftTailedR = new double[15];  
		leftTailedR[0] = 0.4266;
		leftTailedR[1] = 0.5;
		leftTailedR[2] = 1;
		leftTailedR[3] = 1;
		leftTailedR[4] = 0.9857;
		leftTailedR[5] = 0.2429;
		leftTailedR[6] = 0.5;
		leftTailedR[7] = 1;
		leftTailedR[8] = 1;
		leftTailedR[9] = 0.9999;
		leftTailedR[10] = 0.03733;
		leftTailedR[11] = 0.9946;
		leftTailedR[12] = 0.2613;
		leftTailedR[13] = 0.9279;
		leftTailedR[14] = 0.6613;		

		FisherExact fe = FisherExact.getInstance(100);

		final float PRECISION_LEVEL = 0.001f;
		for (int i = 0; i < argInts.length; i++) {
			System.out.println("\na=" + argInts[i][0] + " b=" + argInts[i][1] + " c=" + argInts[i][2] + " d=" + argInts[i][3]);
			System.out.print("*****Original algorithm: ");
			double cumulativeP = fe.getCumlativeP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
			System.out.println("\tcumulativeP = " + cumulativeP);

			System.out.print("*****Left Tailed: ");
			double leftTailedP = fe.getLeftTailedP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
			System.out.println("\tleftTailedP = " + leftTailedP + " leftTailedR = " + leftTailedR[i]);
			assertTrue(Math.abs(leftTailedP - leftTailedR[i]) < PRECISION_LEVEL);

			System.out.print("*****Right Tailed: ");
			double rightTailedP = fe.getRightTailedP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
			System.out.println("\trightTailedP = " + rightTailedP + " rightTailedR = " + rightTailedR[i]);
			assertTrue(Math.abs(rightTailedP - rightTailedR[i]) < PRECISION_LEVEL);

			System.out.print("*****Two Tailed: ");
			double twoTailedP = fe.getTwoTailedP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
			System.out.println("\ttwoTailedP = " + twoTailedP + " twoTailedR = " + twoTailedR[i]);
			assertTrue(Math.abs(twoTailedP - twoTailedR[i]) < PRECISION_LEVEL);
		}
	}
}
