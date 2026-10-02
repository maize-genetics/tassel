/*
 *  MemoryUsagePlugin
 * 
 *  Created on April 18, 2016
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.Sizeof;
import net.maizegenetics.plugindef.GeneratedGuiBoilerplate;

/**
 *
 * @author Terry Casstevens
 */
public class MemoryUsagePlugin extends AbstractPlugin {

    public MemoryUsagePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        Sizeof.printMemoryUse();
        return input;
    }

    @Override
    @GeneratedGuiBoilerplate
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    @GeneratedGuiBoilerplate
    public String getButtonName() {
        return "Memory Usage";
    }

    @Override
    @GeneratedGuiBoilerplate
    public String getToolTipText() {
        return "Memory Usage";
    }

}
