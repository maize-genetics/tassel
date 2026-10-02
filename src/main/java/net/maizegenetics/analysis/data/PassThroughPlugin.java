/*
 * PassThroughPlugin.java
 *
 * Created on June 9, 2007
 *
 */

package net.maizegenetics.analysis.data;


import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginEvent;

import javax.swing.*;
import net.maizegenetics.plugindef.GeneratedGuiBoilerplate;

/**
 *
 * @author terryc
 */
public class PassThroughPlugin extends AbstractPlugin{
    
    
    /**
     * Creates a new instance of PassThroughPlugin
     */
    public PassThroughPlugin() {
        super(null, false);
    }
    
    public DataSet performFunction(DataSet dataSet) {
        DataSet result = new DataSet(dataSet.getDataSet(), this);
        fireDataSetReturned(result);
        return result;
    }
    
    @GeneratedGuiBoilerplate
    public String getToolTipText() {
        return "";
    }
    
    @GeneratedGuiBoilerplate
    public ImageIcon getIcon() {
        return null;
    }
    
    @GeneratedGuiBoilerplate
    public String getButtonName() {
        return "";
    }
    
    public void dataSetReturned(PluginEvent event) {
        DataSet input = (DataSet) event.getSource();
        performFunction(input);
    }
    
}
