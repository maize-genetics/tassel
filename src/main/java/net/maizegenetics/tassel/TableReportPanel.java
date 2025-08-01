package net.maizegenetics.tassel;

import net.maizegenetics.gui.TableReportNoPagingTableModel;
import net.maizegenetics.util.TableReport;

import net.maizegenetics.util.DoubleFormat;

import javax.swing.*;
import javax.swing.table.*;
import java.awt.*;
import java.util.Map;
import java.util.WeakHashMap;

import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListTableReport;
import net.maizegenetics.dna.map.TOPMInterface;
import net.maizegenetics.dna.map.TOPMTableReport;
import net.maizegenetics.dna.snp.FilterList;
import net.maizegenetics.dna.snp.FilterTableReport;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListTableReport;

/**
 */
public class TableReportPanel extends JPanel {

    private static final Map<Object, TableReportPanel> INSTANCES = new WeakHashMap<>();

    private JTable myDataTable;
    private int myTaxaColumnIndex = -1;

    private TableReportPanel(TableReport theTableSource) {

        TableModel theModel = null;

        if (theTableSource instanceof TableModel) {
            theModel = (TableModel) theTableSource;
        } else {
            theModel = new TableReportNoPagingTableModel(theTableSource);
        }

        myDataTable = new JTable(theModel);

        myDataTable.setDefaultRenderer(Double.class, new DefaultTableCellRenderer() {
            @Override
            public void setValue(Object value) {
                if (value == null) {
                    setText("NaN");
                } else if (value.getClass() == Integer.class) {
                    setText(String.valueOf(value));
                } else if (value.getClass() == Double.class) {
                    setText(DoubleFormat.format((Double) value));
                } else {
                    setText(value.toString());
                }
            }
        });

        ((DefaultTableCellRenderer) myDataTable.getDefaultRenderer(String.class)).setHorizontalAlignment(JLabel.RIGHT);
        ((DefaultTableCellRenderer) myDataTable.getDefaultRenderer(Double.class)).setHorizontalAlignment(JLabel.RIGHT);
        ((DefaultTableCellRenderer) myDataTable.getDefaultRenderer(Integer.class)).setHorizontalAlignment(JLabel.RIGHT);

        myDataTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        TableRowSorter<TableModel> theSorter = new TableRowSorter<TableModel>(theModel);
        myDataTable.setRowSorter(theSorter);
        JScrollPane jsp = new JScrollPane(myDataTable);

        String[][] rowtable = getRowHeadings(theTableSource);
        if (myTaxaColumnIndex != -1) {

            JTable rowNameTable = new JTable(rowtable, new String[]{"Taxa"}) {
                public boolean isCellEditable(int row, int column) {
                    return false;
                }
            };
            rowNameTable.setRowSorter(theSorter);
            jsp.setRowHeaderView(rowNameTable);
            int numTaxaToCheck = Math.min(500, rowtable.length);
            int maxChars = 0;
            for (int i = 0; i < numTaxaToCheck; i++) {
                if (rowtable[i][0].length() > maxChars) {
                    maxChars = rowtable[i][0].length();
                }
            }
            int preferedWidth = Math.max(100, maxChars * 9);
            jsp.getRowHeader().setPreferredSize(new Dimension(preferedWidth, myDataTable.getHeight()));
            jsp.setCorner(JScrollPane.UPPER_LEFT_CORNER, rowNameTable.getTableHeader());

            TableColumn taxaColumn = myDataTable.getColumnModel().getColumn(myTaxaColumnIndex);
            myDataTable.getColumnModel().removeColumn(taxaColumn);
        }

        setLayout(new BorderLayout());
        add(jsp, BorderLayout.CENTER);

        setVisible(true);

    }

    public static TableReportPanel getInstance(TableReport theTableSource) {
        TableReportPanel result = INSTANCES.get(theTableSource);
        if (result == null) {
            result = new TableReportPanel(theTableSource);
            INSTANCES.put(theTableSource, result);
        }
        return result;
    }

    public static TableReportPanel getInstance(FilterList filter) {
        TableReportPanel result = INSTANCES.get(filter);
        if (result == null) {
            result = new TableReportPanel(new FilterTableReport(filter));
            INSTANCES.put(filter, result);
        }
        return result;
    }

    public static TableReportPanel getInstance(TOPMInterface topm) {
        TableReportPanel result = INSTANCES.get(topm);
        if (result == null) {
            result = new TableReportPanel(new TOPMTableReport(topm));
            INSTANCES.put(topm, result);
        }
        return result;
    }

    public static TableReportPanel getInstance(TaxaList taxaList) {
        TableReportPanel result = INSTANCES.get(taxaList);
        if (result == null) {
            result = new TableReportPanel(new TaxaListTableReport(taxaList));
            INSTANCES.put(taxaList, result);
        }
        return result;
    }

    public static TableReportPanel getInstance(PositionList positionList) {
        TableReportPanel result = INSTANCES.get(positionList);
        if (result == null) {
            result = new TableReportPanel(new PositionListTableReport(positionList));
            INSTANCES.put(positionList, result);
        }
        return result;
    }

    private String[][] getRowHeadings(TableReport report) {

        Object[] headings = report.getTableColumnNames();

        if (headings.length < 2) {
            return null;
        }

        for (int i = 0, n = headings.length; i < n; i++) {
            if (headings[i].toString().equalsIgnoreCase("taxa")) {
                myTaxaColumnIndex = i;
                break;
            }
        }

        if (myTaxaColumnIndex == -1) {
            return null;
        }

        int numRows = (int) Math.min((long) Integer.MAX_VALUE, report.getRowCount());
        String[][] result = new String[numRows][1];
        for (int i = 0, n = numRows; i < n; i++) {
            result[i][0] = report.getValueAt(i, myTaxaColumnIndex).toString();
        }

        return result;

    }
}
