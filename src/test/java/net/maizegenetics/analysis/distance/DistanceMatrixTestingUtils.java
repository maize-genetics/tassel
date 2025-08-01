/*
 *  DistanceMatrixTestingUtils
 * 
 *  Created on Dec 18, 2015
 */
package net.maizegenetics.analysis.distance;

import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.DistanceMatrixWithCounts;
import static org.junit.Assert.assertEquals;

/**
 *
 * @author Terry Casstevens
 */
public class DistanceMatrixTestingUtils {

    private DistanceMatrixTestingUtils() {
        // utility
    }

    public static void compare(DataSet first, DataSet second, double delta) {
        compare((DistanceMatrix) first.getDataOfType(DistanceMatrix.class).get(0).getData(), (DistanceMatrix) second.getDataOfType(DistanceMatrix.class).get(0).getData(), delta);
    }

    public static void compare(DistanceMatrix first, DistanceMatrix second, double delta) {

        int numColumns = first.getColumnCount() - 1;
        int numRows = (int) first.getRowCount();
        assertEquals("Columns not equal: ", numColumns, second.getColumnCount() - 1);
        assertEquals("Rows not equal: ", numRows, second.getRowCount());

        for (int r = 0; r < numRows; r++) {
            for (int c = 0; c < numColumns; c++) {
                assertEquals("Row: " + r + "  Column: " + c + " difference greater than delta: " + delta, first.getDistance(r, c), second.getDistance(r, c), delta);
            }
        }

        if ((first instanceof DistanceMatrixWithCounts) && (second instanceof DistanceMatrixWithCounts)) {
            for (int r = 0; r < numRows; r++) {
                for (int c = 0; c < numColumns; c++) {
                    assertEquals("Row: " + r + "  Column: " + c + " difference greater than delta: " + delta, ((DistanceMatrixWithCounts) first).getCount(r, c), ((DistanceMatrixWithCounts) second).getCount(r, c));
                }
            }
        }

    }

}
