package net.maizegenetics.analysis.chart;

import net.maizegenetics.util.TableReport;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;

import java.util.ArrayList;
import java.util.Vector;

/**
 * <p>Title: </p>
 * <p>Description: This will find the categories and estimate the box and wiskers for a Tablereport</p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: </p>
 * @author Ed Buckler
 * @version 1.0
 */

public class TableReportBoxWhiskerCatDataset extends DefaultBoxAndWhiskerCategoryDataset {
  String[] seriesNames;

  public TableReportBoxWhiskerCatDataset(TableReport theTable, int seriesCategory, int[] seriesY) {
    setTableReport(theTable, seriesCategory, seriesY);
  }

  public boolean setTableReport(TableReport theTable, int seriesCategory, int[] seriesY) {
    int numRows = (int) theTable.getRowCount();
    Vector theCategories = new Vector();
    seriesNames = new String[seriesY.length];
    Object[] theSN = theTable.getTableColumnNames();
    for(int x=0; x<seriesY.length; x++){seriesNames[x]=theSN[seriesY[x]].toString();}
    for (int i = 0; i < numRows; i++) {
      Object current = theTable.getValueAt(i, seriesCategory);
      if(theCategories.contains(current)==false) {
        theCategories.add(current);
      }
    }
    ArrayList[][] catData=new ArrayList[theCategories.size()][seriesY.length];
    for (int i = 0; i < theCategories.size(); i++) {
       for(int x=0; x<seriesY.length; x++){catData[i][x]=new ArrayList();}
    }
    for (int i = 0; i < numRows; i++) {
      int cat=theCategories.indexOf(theTable.getValueAt(i, seriesCategory));
      Double d;
      for(int x=0; x<seriesY.length; x++){
        try {d=new Double(theTable.getValueAt(i,seriesY[x]).toString());
            if(d.isNaN()==false) catData[cat][x].add(d);}
        catch (NumberFormatException ex) {}
      }
    }
    for (int i = 0; i < theCategories.size(); i++) {
      for (int x = 0; x < seriesY.length; x++) {
        //This throw errors when were are zero, but the graph is fine
        this.add(catData[i][x], seriesNames[x], theCategories.get(i).toString());
      }
    }
      return true;
  }

}
