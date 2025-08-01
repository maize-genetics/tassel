/*
 *  TasselPipelineStepsDialog
 * 
 *  Created on Mar 23, 2015
 */
package net.maizegenetics.pipeline;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dialog;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Insets;
import java.awt.Rectangle;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import javax.swing.AbstractAction;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;
import javax.swing.ScrollPaneConstants;
import javax.swing.text.MutableAttributeSet;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyledDocument;
import net.maizegenetics.analysis.association.EqtlAssociationPlugin;
import net.maizegenetics.analysis.association.FixedEffectLMPlugin;
import net.maizegenetics.analysis.association.MLMPlugin;
import net.maizegenetics.analysis.association.RidgeRegressionEmmaPlugin;
import net.maizegenetics.analysis.association.WeightedMLMPlugin;
import net.maizegenetics.analysis.chart.ChartDisplayPlugin;
import net.maizegenetics.analysis.chart.ManhattanDisplayPlugin;
import net.maizegenetics.analysis.chart.QQDisplayPlugin;
import net.maizegenetics.analysis.chart.TableDisplayPlugin;
import net.maizegenetics.analysis.data.CombineDataSetsPlugin;
import net.maizegenetics.analysis.data.ExportPlugin;
import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.data.GenotypeSummaryPlugin;
import net.maizegenetics.analysis.data.GetPositionListPlugin;
import net.maizegenetics.analysis.data.GetTaxaListPlugin;
import net.maizegenetics.analysis.data.HetsToUnknownPlugin;
import net.maizegenetics.analysis.data.IntersectionAlignmentPlugin;
import net.maizegenetics.analysis.data.MergeGenotypeTablesPlugin;
import net.maizegenetics.analysis.data.PrincipalComponentsPlugin;
import net.maizegenetics.analysis.data.SeparatePlugin;
import net.maizegenetics.analysis.data.SortGenotypeFilePlugin;
import net.maizegenetics.analysis.data.SynonymizerPlugin;
import net.maizegenetics.analysis.data.UnionAlignmentPlugin;
import net.maizegenetics.analysis.distance.DistanceMatrixPlugin;
import net.maizegenetics.analysis.distance.KinshipPlugin;
import net.maizegenetics.analysis.distance.MultiDimensionalScalingPlugin;
import net.maizegenetics.analysis.filter.FilterAlignmentPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.analysis.filter.FilterSiteNamePlugin;
import net.maizegenetics.analysis.filter.FilterTaxaBuilderPlugin;
import net.maizegenetics.analysis.filter.FilterTaxaPropertiesPlugin;
import net.maizegenetics.analysis.filter.FilterTraitsPlugin;
import net.maizegenetics.analysis.imputation.FILLINFindHaplotypesPlugin;
import net.maizegenetics.analysis.imputation.FILLINImputationPlugin;
import net.maizegenetics.analysis.imputation.FSFHapImputationPlugin;
import net.maizegenetics.analysis.imputation.RemoveIndelsForBeaglePlugin;
import net.maizegenetics.analysis.modelfitter.StepwiseOLSModelFitterPlugin;
import net.maizegenetics.analysis.numericaltransform.ImputationPlugin;
import net.maizegenetics.analysis.numericaltransform.NumericalGenotypePlugin;
import net.maizegenetics.analysis.numericaltransform.TransformDataPlugin;
import net.maizegenetics.analysis.popgen.LinkageDiseqDisplayPlugin;
import net.maizegenetics.analysis.popgen.LinkageDisequilibriumPlugin;
import net.maizegenetics.analysis.popgen.SequenceDiversityPlugin;
import net.maizegenetics.analysis.tree.ArchaeopteryxPlugin;
import net.maizegenetics.analysis.tree.CreateTreePlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Plugin;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.PluginListener;
import net.maizegenetics.progress.ProgressPanel;
import net.maizegenetics.tassel.TASSELMainFrame;
import net.maizegenetics.util.Utils;

/**
 *
 * @author Terry Casstevens
 */
public class TasselPipelineStepsDialog extends JDialog implements PluginListener {

    private static final Map<Class, String> MENU_LOCATIONS = new HashMap<>();

    static {
        MENU_LOCATIONS.put(FileLoadPlugin.class, "Data");
        MENU_LOCATIONS.put(ExportPlugin.class, "Data");
        MENU_LOCATIONS.put(GetTaxaListPlugin.class, "Data");
        MENU_LOCATIONS.put(GetPositionListPlugin.class, "Data");
        MENU_LOCATIONS.put(SortGenotypeFilePlugin.class, "Data");
        MENU_LOCATIONS.put(SynonymizerPlugin.class, "Data");
        MENU_LOCATIONS.put(IntersectionAlignmentPlugin.class, "Data");
        MENU_LOCATIONS.put(UnionAlignmentPlugin.class, "Data");
        MENU_LOCATIONS.put(MergeGenotypeTablesPlugin.class, "Data");
        MENU_LOCATIONS.put(SeparatePlugin.class, "Data");
        MENU_LOCATIONS.put(HetsToUnknownPlugin.class, "Data");
        MENU_LOCATIONS.put(TransformDataPlugin.class, "Data");
        MENU_LOCATIONS.put(NumericalGenotypePlugin.class, "Data");

        MENU_LOCATIONS.put(FILLINFindHaplotypesPlugin.class, "Impute");
        MENU_LOCATIONS.put(FILLINImputationPlugin.class, "Impute");
        MENU_LOCATIONS.put(FSFHapImputationPlugin.class, "Impute");
        MENU_LOCATIONS.put(ImputationPlugin.class, "Impute");
        MENU_LOCATIONS.put(RemoveIndelsForBeaglePlugin.class, "Impute");

        MENU_LOCATIONS.put(FilterSiteBuilderPlugin.class, "Filter");
        MENU_LOCATIONS.put(FilterAlignmentPlugin.class, "Filter");
        MENU_LOCATIONS.put(FilterSiteNamePlugin.class, "Filter");
        MENU_LOCATIONS.put(FilterTaxaBuilderPlugin.class, "Filter");
        MENU_LOCATIONS.put(FilterTaxaPropertiesPlugin.class, "Filter");
        MENU_LOCATIONS.put(FilterTraitsPlugin.class, "Filter");

        MENU_LOCATIONS.put(SequenceDiversityPlugin.class, "Analysis");
        MENU_LOCATIONS.put(LinkageDisequilibriumPlugin.class, "Analysis");
        MENU_LOCATIONS.put(DistanceMatrixPlugin.class, "Analysis");
        MENU_LOCATIONS.put(CreateTreePlugin.class, "Analysis");
        MENU_LOCATIONS.put(KinshipPlugin.class, "Analysis");
        MENU_LOCATIONS.put(PrincipalComponentsPlugin.class, "Analysis");
        MENU_LOCATIONS.put(MultiDimensionalScalingPlugin.class, "Analysis");
        MENU_LOCATIONS.put(FixedEffectLMPlugin.class, "Analysis");
        MENU_LOCATIONS.put(MLMPlugin.class, "Analysis");
        MENU_LOCATIONS.put(WeightedMLMPlugin.class, "Analysis");
        MENU_LOCATIONS.put(RidgeRegressionEmmaPlugin.class, "Analysis");
        MENU_LOCATIONS.put(GenotypeSummaryPlugin.class, "Analysis");
        MENU_LOCATIONS.put(StepwiseOLSModelFitterPlugin.class, "Analysis");
        MENU_LOCATIONS.put(EqtlAssociationPlugin.class, "Analysis");

        MENU_LOCATIONS.put(TableDisplayPlugin.class, "Results");
        MENU_LOCATIONS.put(ArchaeopteryxPlugin.class, "Results");
        MENU_LOCATIONS.put(LinkageDiseqDisplayPlugin.class, "Results");
        MENU_LOCATIONS.put(ChartDisplayPlugin.class, "Results");
        MENU_LOCATIONS.put(QQDisplayPlugin.class, "Results");
        MENU_LOCATIONS.put(ManhattanDisplayPlugin.class, "Results");

        MENU_LOCATIONS.put(CombineDataSetsPlugin.class, "Select Multiple Datasets while holding Ctrl or Command");
    }

    private static final Color CURRENT_STEP_COLOR = new Color(0x2c712b);

    private final List<Plugin> myPluginOrder = new ArrayList<>();
    private final List<String> myPluginDescriptions = new ArrayList<>();
    private final Map<Plugin, JTextPane> myTextAreas = new LinkedHashMap<>();
    private final JPanel myMainPane = new JPanel();
    private String myDescription = null;
    private String myCitation = null;
    private final ProgressPanel myProgressPanel;

    public TasselPipelineStepsDialog(TASSELMainFrame mainFrame, String name) {
        super((Window) null, "Tassel Workflow: " + name, Dialog.ModalityType.MODELESS);
        setLayout(new BorderLayout());
        myMainPane.setLayout(new BoxLayout(myMainPane, BoxLayout.Y_AXIS));
        JScrollPane scroll = new JScrollPane(myMainPane);
        scroll.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        getContentPane().add(scroll, BorderLayout.CENTER);
        JButton ok = new JButton(new AbstractAction("Ok") {

            @Override
            public void actionPerformed(ActionEvent e) {
                setVisible(false);
            }
        });
        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new FlowLayout());
        buttonPanel.add(ok);
        getContentPane().add(buttonPanel, BorderLayout.SOUTH);
        setPreferredSize(new Dimension(400, 650));
        setResizable(false);

        myProgressPanel = mainFrame.getProgressPanel();
    }

    public void setOverallDescription(String description) {
        myDescription = description;
    }

    public void setCitation(String citation) {
        myCitation = citation;
    }

    public void addPlugin(Plugin plugin, String description) {
        myPluginOrder.add(plugin);
        myPluginDescriptions.add(description);
        plugin.addListener(this);
    }

    public void showDialog() {

        if (myDescription != null) {
            JTextPane text = new JTextPane();
            text.setContentType("text/html");
            text.setText("<h3>" + myDescription + "</h3>");

            text.setMargin(new Insets(5, 10, 3, 10));
            text.setEditable(false);
            myMainPane.add(text);
        }

        if (myCitation != null) {
            JTextPane text = new JTextPane();
            text.setContentType("text/html");
            text.setText("<h3>" + myCitation + "</h3>");

            text.setMargin(new Insets(3, 10, 5, 10));
            text.setEditable(false);
            myMainPane.add(text);
        }

        boolean firstCreated = true;
        for (int i = 0; i < myPluginOrder.size(); i++) {
            myMainPane.add(Box.createRigidArea(new Dimension(20, 2)));
            Plugin plugin = myPluginOrder.get(i);
            Class clazz = plugin.getClass();
            StringBuilder builder = new StringBuilder();
            builder.append("<html>");
            String menu = MENU_LOCATIONS.get(clazz);
            builder.append("<h3>");
            if (menu == null) {
                builder.append(Utils.getBasename(clazz.getName()));
            } else {
                builder.append(menu);
                builder.append(" -> ");
                builder.append(plugin.getButtonName());
            }
            builder.append("</h3>\n");
            if ((myPluginDescriptions.get(i) == null) || (myPluginDescriptions.get(i).length() == 0)) {

            } else {
                builder.append(myPluginDescriptions.get(i));
            }
            builder.append("</html>");
            JTextPane text = new JTextPane();
            text.setContentType("text/html");
            text.setText(builder.toString());

            text.setMargin(new Insets(10, 10, 10, 10));
            text.setEditable(false);
            //JScrollPane textScroll = new JScrollPane();
            //textScroll.setViewportView(text);
            //textScroll.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
            myMainPane.add(text);
            myTextAreas.put(plugin, text);

            if (firstCreated) {
                setTextColor(plugin, CURRENT_STEP_COLOR);
                firstCreated = false;
            }
        }
        pack();
        setVisible(true);
    }

    @Override
    public synchronized void dataSetReturned(PluginEvent event) {
        DataSet dataSet = (DataSet) event.getSource();
        Plugin plugin = dataSet.getCreator();
        System.out.println("data set returned from: " + plugin.getClass().getName());
        setTextColor(plugin, Color.gray);
        int index = myPluginOrder.indexOf(plugin);
        if (index < myPluginOrder.size() - 1) {
            setTextColor(myPluginOrder.get(index + 1), CURRENT_STEP_COLOR);
            myProgressPanel.addPlugin(myPluginOrder.get(index + 1));
            if (index < myPluginOrder.size() - 2) {
                Rectangle rec = myTextAreas.get(myPluginOrder.get(index + 2)).getBounds();
                myMainPane.scrollRectToVisible(rec);
            }
        }
    }

    private void setTextColor(Plugin plugin, Color color) {
        JTextPane pane = myTextAreas.get(plugin);
        MutableAttributeSet attrs = pane.getInputAttributes();
        StyleConstants.setForeground(attrs, color);
        StyledDocument doc = pane.getStyledDocument();
        doc.setCharacterAttributes(0, doc.getLength(), attrs, true);
    }

    @Override
    public void progress(PluginEvent event) {
        // do nothing
    }

}
