/*
 * TableReportNoPagingTableModel.java
 *
 */
package net.maizegenetics.gui;

import javax.swing.table.AbstractTableModel;

import net.maizegenetics.util.TableReport;

/**
 *
 * @author Terry Casstevens
 */
public class TableReportNoPagingTableModel extends AbstractTableModel {

    private TableReport myTable = null;
    private Object[] myColumnHeadings = null;

    public TableReportNoPagingTableModel(TableReport table) {
        if (table == null) {
            throw new IllegalArgumentException("TableReportNoPagingTableModel: init: table can not be null.");
        }

        myTable = table;
        myColumnHeadings = myTable.getTableColumnNames();

    }

    // Return values appropriate for the visible table part
    @Override
    public int getRowCount() {
        return (int) Math.min((long) Integer.MAX_VALUE, myTable.getRowCount());
    }

    @Override
    public int getColumnCount() {
        return myTable.getColumnCount();
    }

    @Override
    public Object getValueAt(int row, int col) {
        return myTable.getValueAt(row, col);
    }

    @Override
    public String getColumnName(int col) {
        return myColumnHeadings[col].toString();
    }

    /**
     * Resets the table backing this matrix table model to an empty table.
     */
    public void resetTable() {
    }

    public Object getColumnObject(int columnIndex) {
        return myColumnHeadings[columnIndex];
    }

    /**
     * Always returns false.
     */
    @Override
    public boolean isCellEditable(int rowIndex, int columnIndex) {
        return false;
    }

    /**
     * No operation.
     */
    @Override
    public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
        // NO OPERATION
    }

    public void fireTableChanged() {
        fireTableStructureChanged();
    }

    @Override
    public Class<?> getColumnClass(int columnIndex) {
        Object value = getValueAt(0, columnIndex);
        if (value == null) {
            return String.class;
        } else {
            return value.getClass();
        }
    }
}
