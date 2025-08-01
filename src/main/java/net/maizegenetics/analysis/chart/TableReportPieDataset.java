package net.maizegenetics.analysis.chart;

import net.maizegenetics.util.TableReport;
import org.jfree.data.general.DefaultPieDataset;

import java.util.Vector;

/**
 * <p>Title: </p>
 * <p>Description: This will find the categories and frequency of categories from a Tablereport</p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: USDA-ARS</p>
 * @author Ed Buckler
 * @version 1.0
 */

public class TableReportPieDataset extends DefaultPieDataset {
  String seriesNames;

  public TableReportPieDataset(TableReport theTable, int seriesCategory) {
    setTableReport(theTable, seriesCategory);
  }

  public boolean setTableReport(TableReport theTable, int seriesCategory) {
    Vector theCategories = new Vector();
    Object[] theSN = theTable.getTableColumnNames();
    seriesNames=theSN[seriesCategory].toString();
    for (int i = 0; i < theTable.getRowCount(); i++) {
        Object current = theTable.getValueAt(i, seriesCategory);
      if(theCategories.contains(current)==false) {
        theCategories.add(current);
      }
    }
    int[] catCount=new int[theCategories.size()];
    for (int i = 0; i < theTable.getRowCount(); i++) {
        Object current = theTable.getValueAt(i, seriesCategory);
      int cat=theCategories.indexOf(current);
      catCount[cat]++;
    }
    for (int i = 0; i < theCategories.size(); i++) {
      this.setValue(theCategories.get(i).toString(), catCount[i]);
    }
      return true;
  }

}
