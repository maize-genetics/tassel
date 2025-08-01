/*
 * TableReportManhattanDataset
 */
package net.maizegenetics.analysis.chart;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import net.maizegenetics.analysis.association.AssociationConstants;
import net.maizegenetics.util.TableReport;
import org.jfree.data.Range;
import org.jfree.data.xy.DefaultTableXYDataset;

/**
 *
 * @author yz79
 * @author Terry Casstevens
 */
public class TableReportManhattanDataset extends DefaultTableXYDataset {

    private String[] seriesNames;
    private int[] seriesOffsets;
    private String myTrait;
    private int numberYAxes;
    private String[] myChromNames;
    private Object[] myColumnNames;
    private double[] myPValues;
    private double[] myLogPValues;
    private double myMaxLogValue;
    private long[] myPositions;
    private String[] myMarkers;
    private HashMap<Double, Integer> myLookupTable;
    private int myPValueColumnIndex = -1;
    private int myChromColumnIndex = -1;
    private int myPositionColumnIndex = -1;
    private int myMarkerColumnIndex = -1;
    private int myTraitColumnIndex = -1;
    private final int myNumRows;
    private final int myStartIndex;
    private final int myEndIndex;
    private final List<Long> myActualPositions = new ArrayList<>();

    public TableReportManhattanDataset(TableReport theTable, int start, int end) {
        numberYAxes = 0;
        myStartIndex = start;
        myEndIndex = end;
        myNumRows = myEndIndex - start;
        setTableReport(theTable);
    }

    @Override
    public Range getDomainBounds(boolean includeInterval) {
        return new Range(getDomainLowerBound(includeInterval), getDomainUpperBound(includeInterval));
    }

    @Override
    public double getDomainUpperBound(boolean includeInterval) {
        return myPositions[myNumRows - 1];
    }

    @Override
    public double getDomainLowerBound(boolean includeInterval) {
        return myPositions[0];
    }

    public Range getRangeBounds() {
        return new Range(0.0, myMaxLogValue);
    }

    @Override
    public int getItemCount(int parm1) {
        return seriesOffsets[parm1 + 1] - seriesOffsets[parm1];
    }

    @Override
    public Number getX(int series, int item) {
        return myPositions[seriesOffsets[series] + item];
    }

    @Override
    public int getSeriesCount() {
        return numberYAxes;
    }

    @Override
    public Number getY(int series, int item) {
        return myLogPValues[seriesOffsets[series] + item];
    }

    public String getSeriesName(int series) {
        return seriesNames[series];
    }

    @Override
    public String getSeriesKey(int series) {
        return seriesNames[series];
    }

    public String getXName() {
        return "Position";
    }

    public List<Long> getActualPostions() {
        return myActualPositions;
    }

    private void setTraitColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals(AssociationConstants.STATS_HEADER_TRAIT)) {
                myTraitColumnIndex = i;
                return;
            }
        }
    }

    private void setPValueColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals(AssociationConstants.STATS_HEADER_P_VALUE)) {
                myPValueColumnIndex = i;
                return;
            }
        }
        throw new IllegalArgumentException("No P-values in selected data");
    }

    private void setChromColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals(AssociationConstants.STATS_HEADER_CHR)) {
                myChromColumnIndex = i;
                return;
            }
        }
        throw new IllegalArgumentException("No Chromosome names in selected data");
    }

    private void setMarkerColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals(AssociationConstants.STATS_HEADER_MARKER)) {
                myMarkerColumnIndex = i;
                return;
            }
        }
    }

    private void setPositionColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals(AssociationConstants.STATS_HEADER_POSITION)) {
                myPositionColumnIndex = i;
                return;
            }
        }
        throw new IllegalArgumentException("No positions in selected data");
    }

    private void setMarkers(TableReport myTableReport) {
        for (int i = 0; i < myMarkers.length; i++) {
            myMarkers[i] = ((String) myTableReport.getValueAt(myStartIndex + i, myMarkerColumnIndex));
        }
    }

    private void setPValues(TableReport myTableReport) {
        Object temp = myTableReport.getValueAt(myStartIndex, myPValueColumnIndex);
        if (temp instanceof Double) {
            for (int i = 0; i < myPValues.length; i++) {
                myPValues[i] = ((Double) myTableReport.getValueAt(myStartIndex + i, myPValueColumnIndex));
                if (myPValues[i] == 0) {
                    myPValues[i] = Double.MIN_VALUE;
                }
                myLookupTable.put(myPValues[i], i);
            }
        } else if (temp instanceof String) {
            for (int i = 0; i < myPValues.length; i++) {
                myPValues[i] = Double.parseDouble((String) myTableReport.getValueAt(myStartIndex + i, myPValueColumnIndex));
                if (myPValues[i] == 0) {
                    myPValues[i] = Double.MIN_VALUE;
                }
                myLookupTable.put(myPValues[i], i);
            }
        } else {
            throw new IllegalStateException("TableReportManhattanDataset: setPValues: Unknown data type of P values: " + temp.getClass().getName());
        }
    }

    private void setPositions(TableReport myTableReport) {

        long offset = 0;
        int currentPosition = 0;
        int previousPosition = 0;
        String currentChrom = "";
        List<Integer> offsets = new ArrayList<>();
        // GLM positions formatted as int
        boolean isNewChromosome = true;
        List<String> chromosomeNames = new ArrayList<>();
        for (int i = 0; i < myNumRows; i++) {

            currentPosition = Integer.valueOf((myTableReport.getValueAt(myStartIndex + i, myPositionColumnIndex)).toString());

            myChromNames[i] = myTableReport.getValueAt(myStartIndex + i, myChromColumnIndex).toString();
            if (!currentChrom.equals(myChromNames[i])) {
                numberYAxes++;
                currentChrom = myChromNames[i];
                chromosomeNames.add(currentChrom);
                previousPosition = 0;
                offsets.add(i);
                isNewChromosome = true;
            }

            myPositions[i] = currentPosition - previousPosition + offset;
            if (isNewChromosome) {
                myActualPositions.add(myPositions[i] - currentPosition);
                isNewChromosome = false;
            }
            previousPosition = currentPosition;
            offset = myPositions[i];

        }

        offsets.add(myNumRows);
        seriesOffsets = new int[numberYAxes + 1];
        for (int i = 0; i < numberYAxes + 1; i++) {
            seriesOffsets[i] = offsets.get(i);
        }

        seriesNames = new String[numberYAxes];
        chromosomeNames.toArray(seriesNames);

    }

    private void setTrait(TableReport table) {
        myTrait = (String) table.getValueAt(myStartIndex, myTraitColumnIndex);
    }

    public String getMarker(int series, int item) {
        return myMarkers[seriesOffsets[series] + item];
    }

    public String[] getChroms() {
        return myChromNames;
    }

    public String getTrait() {
        return myTrait;
    }

    private void setLogPValues() {
        for (int i = 0; i < myLogPValues.length; i++) {
            myLogPValues[i] = -Math.log10(myPValues[i]);
            if (myLogPValues[i] > myMaxLogValue) {
                myMaxLogValue = myLogPValues[i];
            }
        }
    }

    public void setTableReport(TableReport theTable) {
        myColumnNames = theTable.getTableColumnNames();
        setPValueColumnIndex();
        setChromColumnIndex();
        setPositionColumnIndex();
        setMarkerColumnIndex();
        setTraitColumnIndex();
        myPValues = new double[myNumRows];
        myLogPValues = new double[myNumRows];
        myChromNames = new String[myNumRows];
        myPositions = new long[myNumRows];
        myMarkers = new String[myNumRows];
        myLookupTable = new HashMap<>(myNumRows);
        setPValues(theTable);
        setLogPValues();
        setPositions(theTable);
        setMarkers(theTable);
        setTrait(theTable);
    }
}
