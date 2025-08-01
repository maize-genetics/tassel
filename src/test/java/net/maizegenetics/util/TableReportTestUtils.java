/*
 * TableReportTestUtils
 */
package net.maizegenetics.util;

import static org.junit.Assert.*;

/**
 *
 * @author Terry Casstevens
 */
public class TableReportTestUtils {

    private TableReportTestUtils() {
        // utility
    }

    public static void compareTableReports(TableReport report1, TableReport report2) {

        long numRows = report1.getRowCount();
        assertEquals("Mismatched Row Count: ", numRows, report2.getRowCount());

        int numColumns = report1.getColumnCount();
        assertEquals("Mismatched Column Count: ", numColumns, report2.getColumnCount());

        for (long r = 0; r < numRows; r++) {
            for (int c = 0; c < numColumns; c++) {
                assertEquals("Different Values in Distance Matrix: Row: " + r + " Column: " + c, report1.getValueAt(r, c).toString(), report2.getValueAt(r, c).toString());
            }
        }
    }

    /**
     * Use this function to compare TableReports that might hold values as Float
     * or Double.
     * <p>
     * @param expected expected TableReport
     * @param observed observed TableReport
     * @param delta the maximum difference between two Double or Float values
     * for which they will be considered equal
     */
    public static void compareTableReports(TableReport expected, TableReport observed, double delta) {
        long numRows = expected.getRowCount();
        assertEquals("Mismatched Row Count: ", numRows, observed.getRowCount());

        int numColumns = expected.getColumnCount();
        assertEquals("Mismatched Column Count: ", numColumns, observed.getColumnCount());

        for (long r = 0; r < numRows; r++) {
            for (int c = 0; c < numColumns; c++) {
                Object expectedValue = expected.getValueAt(r, c);
                Object observedValue = observed.getValueAt(r, c);

                if (observedValue instanceof Double) {
                    double d1 = -1.0;
                    if (expectedValue instanceof Double) {
                        d1 = (Double) expectedValue;
                    } else if (expectedValue instanceof Float) {
                        d1 = ((Float) expectedValue).doubleValue();
                    } else {
                        d1 = Double.parseDouble(expectedValue.toString());
                    }
                    double d2 = (Double) observedValue;
                    assertEquals("Observed is Double: Different Values in Distance Matrix: Row: " + r + " Column: " + c, d1, d2, delta);
                } else if (observedValue instanceof Float) {
                    double d1 = -1.0;
                    if (expectedValue instanceof Double) {
                        d1 = (Double) expectedValue;
                    } else if (expectedValue instanceof Float) {
                        d1 = ((Float) expectedValue).doubleValue();
                    } else {
                        d1 = Double.parseDouble(expectedValue.toString());
                    }
                    double d2 = ((Float) observedValue).doubleValue();
                    assertEquals("Observed is Float: Different Values in Distance Matrix: Row: " + r + " Column: " + c, d1, d2, delta);
                } else if (expectedValue instanceof Double) {
                    double d1 = (Double) expectedValue;
                    double d2 = -1.0;
                    if (observedValue instanceof Double) {
                        d2 = (Double) observedValue;
                    } else if (observedValue instanceof Float) {
                        d2 = ((Float) observedValue).doubleValue();
                    } else {
                        d2 = Double.parseDouble(observedValue.toString());
                    }
                    assertEquals("Expected is Double: Different Values in Distance Matrix: Row: " + r + " Column: " + c, d1, d2, delta);
                } else if (expectedValue instanceof Float) {
                    double d1 = ((Float) expectedValue).doubleValue();
                    double d2 = -1.0;
                    if (observedValue instanceof Double) {
                        d2 = (Double) observedValue;
                    } else if (observedValue instanceof Float) {
                        d2 = ((Float) observedValue).doubleValue();
                    } else {
                        d2 = Double.parseDouble(observedValue.toString());
                    }
                    assertEquals("Expected is Float: Different Values in Distance Matrix: Row: " + r + " Column: " + c, d1, d2, delta);
                } else {
                    assertEquals("Comparing as Strings: Different Values in Distance Matrix: Row: " + r + " Column: " + c, expectedValue.toString(), observedValue.toString());
                }
            }
        }
    }

    /**
     * Use this function to compare TableReports that might hold values as Float
     * or Double.
     * <p>
     * @param expected expected TableReport
     * @param observed observed TableReport
     * @param delta the maximum difference between two Double or Float values
     * for which they will be considered equal Before comparing the two values,
     * both are divided by the first value. This provides for appropriate
     * comparison of values that are very different from 1, such as p-values.
     */
    public static void compareTableReportValues(TableReport expected, TableReport observed, double delta) {
        long numRows = expected.getRowCount();
        assertEquals("Mismatched Row Count: ", numRows, observed.getRowCount());

        int numColumns = expected.getColumnCount();
        assertEquals("Mismatched Column Count: ", numColumns, observed.getColumnCount());

        for (long r = 0; r < numRows; r++) {
            for (int c = 0; c < numColumns; c++) {
                Object expectedValue = expected.getValueAt(r, c);
                Object observedValue = observed.getValueAt(r, c);

                if (observedValue instanceof Double) {
                    double d1 = -1.0;
                    if (expectedValue instanceof Double) {
                        d1 = (Double) expectedValue;
                    } else if (expectedValue instanceof Float) {
                        d1 = ((Float) expectedValue).doubleValue();
                    } else {
                        d1 = Double.parseDouble(expectedValue.toString());
                    }
                    double d2 = (Double) observedValue;
                    double d3 = d1 / d1;
                    double d4 = d2 / d1;
                    assertEquals("Observed is Double: Different Values in Distance Matrix: Row: " + r + " Column: " + c, d3, d4, delta);
                } else if (observedValue instanceof Float) {
                    double d1 = -1.0;
                    if (expectedValue instanceof Double) {
                        d1 = (Double) expectedValue;
                    } else if (expectedValue instanceof Float) {
                        d1 = ((Float) expectedValue).doubleValue();
                    } else {
                        d1 = Double.parseDouble(expectedValue.toString());
                    }
                    double d2 = ((Float) observedValue).doubleValue();
                    double d3 = d1 / d1;
                    double d4 = d2 / d1;
                    assertEquals("Observed is Float: Different Values in Distance Matrix: Row: " + r + " Column: " + c, d3, d4, delta);
                } else if (expectedValue instanceof Double) {
                    double d1 = (Double) expectedValue;
                    double d2 = -1.0;
                    if (observedValue instanceof Double) {
                        d2 = (Double) observedValue;
                    } else if (observedValue instanceof Float) {
                        d2 = ((Float) observedValue).doubleValue();
                    } else {
                        d2 = Double.parseDouble(observedValue.toString());
                    }
                    double d3 = d1 / d1;
                    double d4 = d2 / d1;
                    assertEquals("Expected is Double: Different Values in Distance Matrix: Row: " + r + " Column: " + c, d3, d4, delta);
                } else if (expectedValue instanceof Float) {
                    double d1 = ((Float) expectedValue).doubleValue();
                    double d2 = -1.0;
                    if (observedValue instanceof Double) {
                        d2 = (Double) observedValue;
                    } else if (observedValue instanceof Float) {
                        d2 = ((Float) observedValue).doubleValue();
                    } else {
                        d2 = Double.parseDouble(observedValue.toString());
                    }
                    double d3 = d1 / d1;
                    double d4 = d2 / d1;
                    assertEquals("Expected is Float: Different Values in Distance Matrix: Row: " + r + " Column: " + c, d3, d4, delta);
                } else {
                    assertEquals("Comparing as Strings: Different Values in Distance Matrix: Row: " + r + " Column: " + c, expectedValue.toString(), observedValue.toString());
                }
            }
        }
    }

}
