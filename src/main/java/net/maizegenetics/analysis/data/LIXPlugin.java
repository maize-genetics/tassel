/*
 *  LIXPlugin
 * 
 *  Created on Sep 9, 2015
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import javax.swing.ImageIcon;
import net.maizegenetics.dna.snp.io.LineIndexBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;

/**
 *
 * @author Terry Casstevens
 */
public class LIXPlugin extends AbstractPlugin {

    private PluginParameter<String> myCreateIndex = new PluginParameter.Builder<String>("createIndex", null, String.class)
            .inFile()
            .description("Create Index for given file.")
            .build();

    public LIXPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        String genotypeFilename = createIndex();
        if ((genotypeFilename != null) && (!genotypeFilename.isEmpty())) {
            LineIndexBuilder.buildHapmapIndex(genotypeFilename);
        }

        return null;

    }

    /**
     * Create Index for given file.
     *
     * @return Create Index
     */
    public String createIndex() {
        return myCreateIndex.value();
    }

    /**
     * Set Create Index. Create Index for given file.
     *
     * @param value Create Index
     *
     * @return this plugin
     */
    public LIXPlugin createIndex(String value) {
        myCreateIndex = new PluginParameter<>(myCreateIndex, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "LIX";
    }

    @Override
    public String getToolTipText() {
        return "Line Index Plugin";
    }

}
