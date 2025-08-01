/*
 * ListPlugins
 */
package net.maizegenetics.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import java.awt.Frame;
import javax.swing.ImageIcon;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Plugin;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;

/**
 *
 * @author Terry Casstevens
 */
public class ListPlugins extends AbstractPlugin {

    private PluginParameter<Boolean> myFull = new PluginParameter.Builder<>("full", false, Boolean.class).build();
    private PluginParameter<Boolean> myUsage = new PluginParameter.Builder<>("usage", false, Boolean.class).build();

    public ListPlugins(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);

        LogManager.getLogger("net.maizegenetics").atLevel(Level.OFF);
        LogManager.getLogger("net.maizegenetics.plugindef").atLevel(Level.INFO);
    }

    @Override
    public DataSet processData(DataSet input) {

        Set<String> classes = Utils.getTasselClasses();
        List<String> temp = new ArrayList<>();
        for (String current : classes) {
            if (Plugin.isPlugin(current)) {
                if (full() || usage()) {
                    temp.add(current);
                } else {
                    temp.add(Utils.getBasename(current));
                }
            }
        }

        Collections.sort(temp);

        for (String current : temp) {
            System.out.println(current);
            if (usage()) {
                System.out.println("--------------------------------------");
                Plugin plugin = Plugin.getPluginInstance(current, null);
                if (plugin != null) {
                    System.out.println(plugin.getUsage());
                }
                System.out.println("\n");
            }
        }

        return null;

    }

    public boolean full() {
        return myFull.value();
    }

    public boolean usage() {
        return myUsage.value();
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "List Plugins";
    }

    @Override
    public String getToolTipText() {
        return "List Plugins";
    }
}
