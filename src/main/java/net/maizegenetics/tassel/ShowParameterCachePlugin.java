package net.maizegenetics.tassel;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.ParameterCache;
import net.maizegenetics.util.TableReportBuilder;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;

/**
 * This plugin creates a table report of the key / value pairs stored in the @{@link ParameterCache}
 *
 * @author Terry Casstevens
 * Created September 08, 2018
 */
public class ShowParameterCachePlugin extends AbstractPlugin {

    public ShowParameterCachePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        Enumeration<?> keys = ParameterCache.keys();
        if (keys == null) {
            throw new IllegalArgumentException("ShowParameterCachePlugin: processData: there is no parameter cache. Use File -> Preferences to set.");
        }

        List<String> keyList = new ArrayList<>();
        while (keys.hasMoreElements()) {
            keyList.add(keys.nextElement().toString());
        }

        keyList.sort((o1, o2) -> {
            int o1Periods = countPeriods(o1);
            int o2Periods = countPeriods(o2);
            if (o1Periods < o2Periods) {
                return -1;
            } else if (o1Periods > o2Periods) {
                return 1;
            }
            return o1.compareTo(o2);
        });

        TableReportBuilder builder = TableReportBuilder.getInstance("Parameter Cache", new String[]{"plugin", "parameter", "value"});
        for (String key : keyList) {
            int index = key.lastIndexOf('.');
            String plugin = "";
            String parameter = null;
            if (index == -1) {
                parameter = key;
            } else {
                plugin = key.substring(0, index);
                parameter = key.substring(index + 1);
            }
            builder.add(new Object[]{plugin, parameter, ParameterCache.value(key).get()});
        }

        return new DataSet(new Datum("Parameter Cache", builder.build(), null), this);

    }

    private int countPeriods(String str) {
        int result = 0;
        for (int i = 0; i < str.length(); i++) {
            if (str.charAt(i) == '.') {
                result++;
            }
        }
        return result;
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = TasselLogging.class.getResource("/net/maizegenetics/analysis/images/properties.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Show Parameter Cache";
    }

    @Override
    public String getToolTipText() {
        return "Show Parameter Cache";
    }


}
