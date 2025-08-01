/*
 *  SuperByteMatrixTest
 */
package net.maizegenetics.util;

import java.util.Random;
import junit.framework.Assert;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author terry
 */
public class SuperByteMatrixTest {

    public SuperByteMatrixTest() {
    }

    @BeforeClass
    public static void setUpClass() {
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    @Test
    public void reorderColumnsTest() {

        System.out.println("SuperByteMatrixTest: reorderColumnsTest...");

        int numRows = 1_000;
        int numColumns = 10_000;

        long previous = System.nanoTime();
        SuperByteMatrix matrix = SuperByteMatrixBuilder.getInstance(numRows, numColumns);
        fillWithRandomData(matrix);
        long current = System.nanoTime();
        System.out.println("   matrix: " + matrix.getClass().getName());
        System.out.println("   create matrix time: " + ((double) (current - previous) / 1_000_000_000.0));

        previous = System.nanoTime();
        SuperByteMatrix copyMatrix = SuperByteMatrixBuilder.getInstanceTranspose(matrix);
        current = System.nanoTime();
        System.out.println("   copyMatrix: " + copyMatrix.getClass().getName());
        System.out.println("   copy matrix time: " + ((double) (current - previous) / 1_000_000_000.0));

        int[] reorderArray = createReorderArray(numColumns);

        previous = System.nanoTime();
        copyMatrix.reorderColumns(reorderArray);
        current = System.nanoTime();
        System.out.println("   reorder matrix time: " + ((double) (current - previous) / 1_000_000_000.0));

        previous = System.nanoTime();
        compareTwoMatricesReorderedColumns(matrix, copyMatrix, reorderArray);
        current = System.nanoTime();
        System.out.println("   compare reorder columns matrix time: " + ((double) (current - previous) / 1_000_000_000.0));

    }

    @Test
    public void reorderRowsTest() {

        System.out.println("SuperByteMatrixTest: reorderRowsTest...");

        int numRows = 1_000;
        int numColumns = 10_000;

        long previous = System.nanoTime();
        SuperByteMatrix matrix = SuperByteMatrixBuilder.getInstance(numRows, numColumns);
        fillWithRandomData(matrix);
        long current = System.nanoTime();
        System.out.println("   matrix: " + matrix.getClass().getName());
        System.out.println("   create matrix time: " + ((double) (current - previous) / 1_000_000_000.0));

        previous = System.nanoTime();
        SuperByteMatrix copyMatrix = SuperByteMatrixBuilder.getInstanceCopy(matrix);
        current = System.nanoTime();
        System.out.println("   copyMatrix: " + copyMatrix.getClass().getName());
        System.out.println("   copy matrix time: " + ((double) (current - previous) / 1_000_000_000.0));

        int[] reorderArray = createReorderArray(numRows);

        previous = System.nanoTime();
        copyMatrix.reorderRows(reorderArray);
        current = System.nanoTime();
        System.out.println("   reorder matrix time: " + ((double) (current - previous) / 1_000_000_000.0));

        previous = System.nanoTime();
        compareTwoMatricesReorderedRows(matrix, copyMatrix, reorderArray);
        current = System.nanoTime();
        System.out.println("   compare reorder rows matrix time: " + ((double) (current - previous) / 1_000_000_000.0));

    }

    private static int[] createReorderArray(int size) {
        int[] result = new int[size];
        for (int i = 0; i < size; i++) {
            result[i] = i;
        }
        Random random = new Random();
        int numSwaps = 1_000_000;
        for (int i = 0; i < numSwaps; i++) {
            int first = random.nextInt(size);
            int second = random.nextInt(size);
            int temp = result[first];
            result[first] = result[second];
            result[second] = temp;
        }
        return result;
    }

    private static void fillWithRandomData(SuperByteMatrix matrix) {
        int numRows = matrix.getNumRows();
        int numColumns = matrix.getNumColumns();
        Random random = new Random();
        if (matrix instanceof SuperByteMatrixTranspose) {
            for (int s = 0; s < numColumns; s++) {
                for (int t = 0; t < numRows; t++) {
                    matrix.set(t, s, (byte) random.nextInt(127));
                }
            }
        } else {
            for (int t = 0; t < numRows; t++) {
                for (int s = 0; s < numColumns; s++) {
                    matrix.set(t, s, (byte) random.nextInt(127));
                }
            }
        }
    }

    private static void compareTwoMatricesReorderedRows(SuperByteMatrix matrix1, SuperByteMatrix matrix2, int[] newIndices) {
        int numRows = matrix1.getNumRows();
        int numColumns = matrix1.getNumColumns();
        if (matrix1 instanceof SuperByteMatrixTranspose) {
            for (int c = 0; c < numColumns; c++) {
                for (int r = 0; r < numRows; r++) {
                    Assert.assertTrue("SuperByteMatrixTest: compareTwoMatricesReorderedRows: Not Equal row: " + r + "  column: " + c + " newRow: " + newIndices[r], (matrix1.get(newIndices[r], c) == matrix2.get(r, c)));
                }
            }
        } else {
            for (int r = 0; r < numRows; r++) {
                for (int c = 0; c < numColumns; c++) {
                    Assert.assertTrue("SuperByteMatrixTest: compareTwoMatricesReorderedRows: Not Equal row: " + r + "  column: " + c + " newRow: " + newIndices[r], (matrix1.get(newIndices[r], c) == matrix2.get(r, c)));
                }
            }
        }
    }

    private static void compareTwoMatricesReorderedColumns(SuperByteMatrix matrix1, SuperByteMatrix matrix2, int[] newIndices) {
        int numRows = matrix1.getNumRows();
        int numColumns = matrix1.getNumColumns();
        if (matrix1 instanceof SuperByteMatrixTranspose) {
            for (int c = 0; c < numColumns; c++) {
                for (int r = 0; r < numRows; r++) {
                    Assert.assertTrue("SuperByteMatrixTest: compareTwoMatricesReorderedColumns: Not Equal row: " + r + "  column: " + c + " newColumn: " + newIndices[c], (matrix1.get(r, newIndices[c]) == matrix2.get(r, c)));
                }
            }
        } else {
            for (int r = 0; r < numRows; r++) {
                for (int c = 0; c < numColumns; c++) {
                    Assert.assertTrue("SuperByteMatrixTest: compareTwoMatricesReorderedColumns: Not Equal row: " + r + "  column: " + c + " newColumn: " + newIndices[c], (matrix1.get(r, newIndices[c]) == matrix2.get(r, c)));
                }
            }
        }
    }
}