/*
 * SequenceDiversityPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.analysis.popgen;

import com.google.common.collect.Range;
import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import javax.swing.*;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.SimpleTableReport;

/**
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public class SequenceDiversityPlugin extends AbstractPlugin {

    private PluginParameter<Integer> myStartSite = new PluginParameter.Builder<>("startSite", 0, Integer.class)
            .description("Start Site")
            .range(Range.atLeast(0))
            .build();

    private PluginParameter<Integer> myEndSite = new PluginParameter.Builder<>("endSite", null, Integer.class)
            .description("End Site")
            .range(Range.atLeast(0))
            .build();

    private PluginParameter<Boolean> myIsSlidingWindowAnalysis = new PluginParameter.Builder<>("slidingWindowAnalysis", false, Boolean.class)
            .description("Whether this is sliding window analysis")
            .build();

    private PluginParameter<Integer> myStepSize = new PluginParameter.Builder<>("stepSize", 100, Integer.class)
            .description("Step Size")
            .range(Range.atLeast(0))
            .dependentOnParameter(myIsSlidingWindowAnalysis)
            .build();

    private PluginParameter<Integer> myWindowSize = new PluginParameter.Builder<>("windowSize", 500, Integer.class)
            .description("Window Size")
            .range(Range.atLeast(0))
            .dependentOnParameter(myIsSlidingWindowAnalysis)
            .build();

    private GenotypeTable myGenotypeTable = null;

    public SequenceDiversityPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
        if (alignInList.size() != 1) {
            throw new IllegalArgumentException("SequenceDiversityPlugin: Please select one Genotype Table");
        }
        myGenotypeTable = (GenotypeTable) alignInList.get(0).getData();
        if (isInteractive()) {
            setParameter(myEndSite, myGenotypeTable.numberOfSites() - 1);
        }
    }

    @Override
    protected void postProcessParameters() {
        if (endSite() == null) {
            setParameter(myEndSite, myGenotypeTable.numberOfSites() - 1);
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        Datum current = input.getDataOfType(GenotypeTable.class).get(0);
        PolymorphismDistribution pda = new PolymorphismDistribution();
        DiversityAnalyses theDA = new DiversityAnalyses(myGenotypeTable, isSlidingWindowAnalysis(),
                startSite(), endSite(), windowSize(), stepSize(), pda);
        List<Datum> results = new ArrayList<>();
        results.add(new Datum("Diversity:" + current.getName(), theDA, "Diversity Analysis"));
        results.add(new Datum("PolyDist:" + current.getName(), new SimpleTableReport(pda.getTableTitle(), pda.getTableColumnNames(), pda.getTableData()), "Polymorphism Distribution"));
        return new DataSet(results, this);
    }

    /**
     * Start Site
     *
     * @return Start Site
     */
    public Integer startSite() {
        return myStartSite.value();
    }

    /**
     * Set Start Site. Start Site
     *
     * @param value Start Site
     *
     * @return this plugin
     */
    public SequenceDiversityPlugin startSite(Integer value) {
        myStartSite = new PluginParameter<>(myStartSite, value);
        return this;
    }

    /**
     * End Site
     *
     * @return End Site
     */
    public Integer endSite() {
        return myEndSite.value();
    }

    /**
     * Set End Site. End Site
     *
     * @param value End Site
     *
     * @return this plugin
     */
    public SequenceDiversityPlugin endSite(Integer value) {
        myEndSite = new PluginParameter<>(myEndSite, value);
        return this;
    }

    /**
     * Whether this is sliding window analysis
     *
     * @return Sliding Window Analysis
     */
    public Boolean isSlidingWindowAnalysis() {
        return myIsSlidingWindowAnalysis.value();
    }

    /**
     * Set Sliding Window Analysis. Whether this is sliding window analysis
     *
     * @param value Sliding Window Analysis
     *
     * @return this plugin
     */
    public SequenceDiversityPlugin isSlidingWindowAnalysis(Boolean value) {
        myIsSlidingWindowAnalysis = new PluginParameter<>(myIsSlidingWindowAnalysis, value);
        return this;
    }

    /**
     * Window Size
     *
     * @return Window Size
     */
    public Integer windowSize() {
        return myWindowSize.value();
    }

    /**
     * Set Window Size. Window Size
     *
     * @param value Window Size
     *
     * @return this plugin
     */
    public SequenceDiversityPlugin windowSize(Integer value) {
        myWindowSize = new PluginParameter<>(myWindowSize, value);
        return this;
    }

    /**
     * Step Size
     *
     * @return Step Size
     */
    public Integer stepSize() {
        return myStepSize.value();
    }

    /**
     * Set Step Size. Step Size
     *
     * @param value Step Size
     *
     * @return this plugin
     */
    public SequenceDiversityPlugin stepSize(Integer value) {
        myStepSize = new PluginParameter<>(myStepSize, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = SequenceDiversityPlugin.class.getResource("/net/maizegenetics/analysis/images/Diversity.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Diversity";
    }

    @Override
    public String getToolTipText() {
        return "Basic description of diversity";
    }
}
