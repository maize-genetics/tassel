/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license, version 2 and without
 * any warranty or technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General
 * public license.
 *
 */
//Title:      TASSELMainFrame
//Version:
//Copyright:  Copyright (c) 1997
//Author:     Ed Buckler
//Company:    NCSU
package net.maizegenetics.tassel;

import net.maizegenetics.analysis.association.FastMultithreadedAssociationPlugin;
import net.maizegenetics.analysis.association.FixedEffectLMPlugin;
import net.maizegenetics.analysis.association.GenomicSelectionPlugin;
import net.maizegenetics.analysis.association.MLMPlugin;
import net.maizegenetics.analysis.association.MultivariateStepwisePlugin;
import net.maizegenetics.analysis.association.WeightedMLMPlugin;
import net.maizegenetics.analysis.chart.ChartDisplayPlugin;
import net.maizegenetics.analysis.chart.ManhattanDisplayPlugin;
import net.maizegenetics.analysis.chart.QQDisplayPlugin;
import net.maizegenetics.analysis.chart.TableDisplayPlugin;
import net.maizegenetics.analysis.data.*;
import net.maizegenetics.analysis.distance.AMatrixPlugin;
import net.maizegenetics.analysis.distance.AddDistanceMatrixPlugin;
import net.maizegenetics.analysis.distance.DistanceMatrixPlugin;
import net.maizegenetics.analysis.distance.HMatrixPlugin;
import net.maizegenetics.analysis.distance.KinshipPlugin;
import net.maizegenetics.analysis.distance.MultiDimensionalScalingPlugin;
import net.maizegenetics.analysis.distance.RemoveNaNFromDistanceMatrixPlugin;
import net.maizegenetics.analysis.distance.SubtractDistanceMatrixPlugin;
import net.maizegenetics.analysis.distance.VCAPScanPlugin;
import net.maizegenetics.analysis.filter.FilterAlignmentPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.analysis.filter.FilterSiteNamePlugin;
import net.maizegenetics.analysis.filter.FilterTaxaBuilderPlugin;
import net.maizegenetics.analysis.filter.FilterTaxaPropertiesPlugin;
import net.maizegenetics.analysis.filter.FilterTraitsPlugin;
import net.maizegenetics.analysis.gbs.*;
import net.maizegenetics.analysis.gbs.v2.DiscoverySNPCallerPluginV2;
import net.maizegenetics.analysis.gbs.v2.GBSSeqToTagDBPlugin;
import net.maizegenetics.analysis.gbs.v2.GetTagSequenceFromDBPlugin;
import net.maizegenetics.analysis.gbs.v2.ProductionSNPCallerPluginV2;
import net.maizegenetics.analysis.gbs.v2.SAMToGBSdbPlugin;
import net.maizegenetics.analysis.gbs.v2.SNPCutPosTagVerificationPlugin;
import net.maizegenetics.analysis.gbs.v2.SNPQualityProfilerPlugin;
import net.maizegenetics.analysis.gbs.v2.TagExportToFastqPlugin;
import net.maizegenetics.analysis.gbs.v2.UpdateSNPPositionQualityPlugin;
import net.maizegenetics.analysis.imputation.ClusterGenotypesPlugin;
import net.maizegenetics.analysis.imputation.FILLINFindHaplotypesPlugin;
import net.maizegenetics.analysis.imputation.FILLINImputationPlugin;
import net.maizegenetics.analysis.imputation.FSFHapImputationPlugin;
import net.maizegenetics.analysis.imputation.ImputationAccuracyPlugin;
import net.maizegenetics.analysis.imputation.LDKNNiImputationPlugin;
import net.maizegenetics.analysis.imputation.RemoveIndelsForBeaglePlugin;
import net.maizegenetics.analysis.modelfitter.StepwiseOLSModelFitterPlugin;
import net.maizegenetics.analysis.numericaltransform.ImputationPlugin;
import net.maizegenetics.analysis.numericaltransform.NumericalGenotypePlugin;
import net.maizegenetics.analysis.numericaltransform.TransformDataPlugin;
import net.maizegenetics.analysis.popgen.LinkageDiseqDisplayPlugin;
import net.maizegenetics.analysis.popgen.LinkageDisequilibriumPlugin;
import net.maizegenetics.analysis.popgen.SequenceDiversityPlugin;
import net.maizegenetics.analysis.tree.*;
import net.maizegenetics.analysis.workflow.WorkflowPlugin;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.gui.PrintHeapAction;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Plugin;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.ThreadedPluginListener;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.progress.ProgressPanel;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import javax.swing.*;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowEvent;
import java.awt.font.TextAttribute;
import java.io.File;
import java.io.InputStream;
import java.net.URI;
import java.net.URL;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * TASSELMainFrame
 */
public class TASSELMainFrame extends JFrame implements ActionListener {

    private static final Logger myLogger = LogManager.getLogger(TASSELMainFrame.class);
    public static final String version = "5.2.96";
    public static final String versionDate = "March 31, 2025";
    private DataTreePanel myDataTreePanel;
    //a variable to control when the progress bar was last updated
    private JFileChooser filerSave = new JFileChooser();
    private JFileChooser filerOpen = new JFileChooser();
    private JScrollPane reportPanelScrollPane = new JScrollPane();
    private JTextArea reportPanelTextArea = new JTextArea();
    JScrollPane mainPanelScrollPane = new JScrollPane();
    JPanel mainDisplayPanel = new JPanel();
    private JTextArea mainPanelTextArea = new JTextArea();
    private JTextField myStatusTextField = new JTextField();
    private final ProgressPanel myProgressPanel = ProgressPanel.getInstance();
    private HashMap<JMenuItem, Plugin> myMenuItemHash = new HashMap<JMenuItem, Plugin>();
    private PreferencesDialog myPreferences = null;

    public TASSELMainFrame() {
        try {
            loadSettings();
            myDataTreePanel = new DataTreePanel(this);
            myDataTreePanel.setToolTipText("Data Tree Panel");
            addMenuBar();
            initializeMyFrame();

            this.setTitle("TASSEL (Trait Analysis by aSSociation, Evolution, and Linkage) " + version);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    //Component initialization
    private void initializeMyFrame() throws Exception {
        getContentPane().setLayout(new BorderLayout());

        int xDim = TasselPrefs.getXDim();
        int yDim = TasselPrefs.getYDim();
        if ((xDim < 50) || (yDim < 50)) {
            Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
            setSize(new Dimension(screenSize.width * 19 / 20, screenSize.height * 19 / 20));
        } else {
            setSize(xDim, yDim);
        }
        setTitle("TASSEL (Trait Analysis by aSSociation, Evolution, and Linkage)");
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                TasselLogging.closeInstance();
                TasselPrefs.putXDim(getWidth());
                TasselPrefs.putYDim(getHeight());
                System.exit(0);
            }
        });

        filerSave.setDialogType(JFileChooser.SAVE_DIALOG);

        JSplitPane dataTreeReportPanelsSplitPanel = new JSplitPane();
        dataTreeReportPanelsSplitPanel.setOrientation(JSplitPane.VERTICAL_SPLIT);

        reportPanelTextArea.setEditable(false);
        reportPanelTextArea.setToolTipText("Report Panel");
        reportPanelTextArea.setLineWrap(true);
        reportPanelTextArea.setWrapStyleWord(true);
        mainPanelTextArea.setDoubleBuffered(true);
        mainPanelTextArea.setEditable(false);
        mainPanelTextArea.setFont(new java.awt.Font("Monospaced", 0, 12));
        mainPanelTextArea.setToolTipText("Main Panel");

        myStatusTextField.setBackground(Color.lightGray);
        myStatusTextField.setBorder(null);

        JSplitPane dataTreeReportMainPanelsSplitPanel = new JSplitPane();

        getContentPane().add(dataTreeReportMainPanelsSplitPanel, BorderLayout.CENTER);

        dataTreeReportMainPanelsSplitPanel.add(dataTreeReportPanelsSplitPanel, JSplitPane.LEFT);

        dataTreeReportPanelsSplitPanel.add(myDataTreePanel, JSplitPane.TOP);

        JSplitPane reportProgressSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);

        JPanel reportPanel = new JPanel(new BorderLayout());

        reportProgressSplitPane.add(reportPanel, JSplitPane.TOP);

        reportPanel.add(reportPanelScrollPane, BorderLayout.CENTER);

        reportPanelScrollPane.getViewport().add(reportPanelTextArea, null);

        reportProgressSplitPane.add(new JScrollPane(myProgressPanel), JSplitPane.BOTTOM);

        dataTreeReportPanelsSplitPanel.add(reportProgressSplitPane, JSplitPane.BOTTOM);

        dataTreeReportMainPanelsSplitPanel.add(mainDisplayPanel, JSplitPane.RIGHT);

        mainDisplayPanel.setLayout(new BorderLayout());
        mainPanelScrollPane.getViewport().add(mainPanelTextArea, null);
        mainDisplayPanel.add(mainPanelScrollPane, BorderLayout.CENTER);

        mainPanelScrollPane.getViewport().add(mainPanelTextArea, null);

        getContentPane().add(myStatusTextField, BorderLayout.SOUTH);

        dataTreeReportMainPanelsSplitPanel.setDividerLocation(getSize().width / 4);

        dataTreeReportPanelsSplitPanel.setDividerLocation((int) (getSize().height / 3.5));

        reportProgressSplitPane.setDividerLocation((int) (getSize().height / 3.5));
    }

    public void loadPreferences() {
        myPreferences.processData(null);
    }

    private void addMenuBar() {

        JMenuBar jMenuBar = new JMenuBar();

        jMenuBar.add(getFileMenu());
        jMenuBar.add(getDataMenu());
        jMenuBar.add(getImputeMenu());
        jMenuBar.add(getFiltersMenu());
        jMenuBar.add(getAnalysisMenu());
        jMenuBar.add(getResultsMenu());
        jMenuBar.add(getGBSv2Menu());
        jMenuBar.add(getGBSMenu());
        jMenuBar.add(getWorkflowMenu());
        addThirdPartyMenus(jMenuBar);
        jMenuBar.add(Box.createHorizontalGlue());
        jMenuBar.add(getHelpMenu());

        this.setJMenuBar(jMenuBar);

        updatePluginsWithGlobalConfigParameters();

    }

    /**
     * This adds plugins to the TASSEL menu specified by files named tassel_menu.xml. This can be used in jars that
     * implement TASSEL plugins to dynamically add a menu.
     *
     * @param jMenuBar TASSEL Menu Bar
     */
    private void addThirdPartyMenus(JMenuBar jMenuBar) {

        for (String current : Utils.getFullyQualifiedResourceNames("tassel_menu.xml")) {
            myLogger.info("TASSELMainFrame: addThirdPartyMenus: adding: " + current);
            addThirdPartyMenus(jMenuBar, current);
        }

    }

    private void addThirdPartyMenus(JMenuBar jMenuBar, String filename) {

        try (InputStream input = TASSELMainFrame.class.getResourceAsStream(filename)) {
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(input);
            doc.getDocumentElement().normalize();
            Element rootElement = doc.getDocumentElement();
            if (!(rootElement.getNodeName().equalsIgnoreCase("TasselMenus"))) {
                throw new IllegalArgumentException("TASSELMainFrame: addThirdPartyMenus: Root Node must be TasselMenus: " + rootElement.getNodeName());
            }
            NodeList children = rootElement.getChildNodes();
            for (int i = 0; i < children.getLength(); i++) {
                Node current = children.item(i);
                JMenu menu = createSubMenu(current);
                if (menu != null) {
                    jMenuBar.add(menu);
                }
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("TASSELMainFrame: addThirdPartyMenus: Problem reading XML file: " + filename + "\n" + e.getMessage());
        }

    }

    private JMenu createSubMenu(Node rootElement) {

        try {

            if (!rootElement.getNodeName().trim().equalsIgnoreCase("TasselMenu")) {
                return null;
            }

            JMenu result = new JMenu("Menu");

            NodeList children = rootElement.getChildNodes();
            for (int i = 0; i < children.getLength(); i++) {
                Node current = children.item(i);
                String elementName = current.getNodeName().trim();
                if ((current.getNodeType() == Node.ELEMENT_NODE) && (elementName.equalsIgnoreCase("name"))) {
                    result.setText(current.getTextContent().trim());
                } else if ((current.getNodeType() == Node.ELEMENT_NODE) && (elementName.equalsIgnoreCase("TasselMenu"))) {
                    JMenu menu = createSubMenu(current);
                    if (menu != null) {
                        result.add(menu);
                    }
                } else if ((current.getNodeType() == Node.ELEMENT_NODE) && (elementName.equalsIgnoreCase("plugin"))) {
                    try {

                        Plugin plugin = null;
                        String possibleClassName = current.getTextContent().trim();
                        List<String> matches = Utils.getFullyQualifiedClassNames(possibleClassName);
                        for (String match : matches) {
                            plugin = Plugin.getPluginInstance(match, this, true);
                            if (plugin != null) {
                                break;
                            }
                        }

                        if (plugin == null) {
                            plugin = Plugin.getPluginInstance(possibleClassName, this, true);
                        }

                        if (plugin != null) {
                            result.add(createMenuItem(plugin, false));
                        }

                    } catch (UnsupportedOperationException usoe) {
                        throw new IllegalArgumentException("TASSELMainFrame: addSubMenu: this plugin is not self-described: " + current);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TASSELMainFrame: addSubMenu: Unknown parameter: " + current);
                    }
                } else if ((current.getNodeType() == Node.ELEMENT_NODE) && (elementName.equalsIgnoreCase("separator"))) {
                    result.addSeparator();
                }

            }

            return result;

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("TASSELMainFrame: addSubMenu: Problem converting XML to Args: " + e.getMessage());
        }

    }

    //Help | About action performed
    private void helpAbout_actionPerformed(ActionEvent e) {
        AboutBox dlg = new AboutBox(this);
        Dimension dlgSize = dlg.getPreferredSize();
        Dimension frmSize = getSize();
        Point loc = getLocation();
        dlg.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x, (frmSize.height - dlgSize.height) / 2 + loc.y);
        dlg.setModal(true);
        dlg.setVisible(true);
    }

    public void sendMessage(String text) {
        myStatusTextField.setForeground(Color.BLACK);
        myStatusTextField.setText(text);
    }

    public void sendErrorMessage(String text) {
        myStatusTextField.setForeground(Color.RED);
        myStatusTextField.setText(text);
    }

    public void setMainText(String text) {
        mainPanelTextArea.setText(text);
        mainPanelTextArea.setCaretPosition(0);
    }

    public void setNoteText(String text) {
        reportPanelTextArea.setText(text);
        reportPanelTextArea.setCaretPosition(0);
    }

    private void loadSettings() {
        filerOpen.setCurrentDirectory(new File(TasselPrefs.getOpenDir()));
        filerSave.setCurrentDirectory(new File(TasselPrefs.getSaveDir()));
    }

    public void addDataSet(DataSet theDataSet, String defaultNode) {
        myDataTreePanel.addDataSet(theDataSet, defaultNode);
    }

    public void updateMainDisplayPanel(JPanel panel) {
        mainDisplayPanel.removeAll();
        mainDisplayPanel.add(panel, BorderLayout.CENTER);
        mainDisplayPanel.repaint();
        mainDisplayPanel.validate();
    }

    public DataTreePanel getDataTreePanel() {
        return myDataTreePanel;
    }

    public ProgressPanel getProgressPanel() {
        return myProgressPanel;
    }

    private JMenuItem createMenuItem(Plugin theTP) {
        return createMenuItem(theTP, -1);
    }

    private JMenuItem createMenuItem(Plugin theTP, boolean adjustForIcon) {
        return createMenuItem(theTP, -1, adjustForIcon);
    }

    private JMenuItem createMenuItem(Plugin theTP, int mnemonic) {
        return createMenuItem(theTP, mnemonic, true);
    }

    private static final int ICON_WIDTH_PLUS_GAP = 30;

    private JMenuItem createMenuItem(Plugin theTP, int mnemonic, boolean adjustForIcon) {
        ImageIcon icon = theTP.getIcon();
        JMenuItem menuItem = new JMenuItem(theTP.getButtonName(), icon);
        if (mnemonic != -1) {
            menuItem.setMnemonic(mnemonic);
        }
        if (adjustForIcon) {
            int pixels = ICON_WIDTH_PLUS_GAP;
            if (icon != null) {
                pixels -= icon.getIconWidth();
                pixels /= 2;
            }
            menuItem.setIconTextGap(pixels);
        }
        menuItem.setBackground(Color.white);
        menuItem.setMargin(new Insets(2, 2, 2, 2));
        menuItem.setToolTipText(theTP.getToolTipText());
        menuItem.addActionListener(this);
        theTP.addListener(myDataTreePanel);
        myMenuItemHash.put(menuItem, theTP);
        return menuItem;
    }

    private JMenuItem createMenuItem(Action action, int mnemonic) {

        Icon icon = (Icon) action.getValue("SmallIcon");

        JMenuItem menuItem = new JMenuItem(action.getValue("Name").toString(), icon);
        if (mnemonic != -1) {
            menuItem.setMnemonic(mnemonic);
        }
        int pixels = ICON_WIDTH_PLUS_GAP;
        if (icon != null) {
            pixels -= icon.getIconWidth();
            pixels /= 2;
        }
        menuItem.setIconTextGap(pixels);
        menuItem.setBackground(Color.white);
        menuItem.setMargin(new Insets(2, 2, 2, 2));
        menuItem.addActionListener(action);
        return menuItem;

    }

    private JMenuItem createMenuItem(String label, boolean adjustForIcon) {

        JMenuItem menuItem = new JMenuItem(label);
        if (adjustForIcon) {
            menuItem.setIconTextGap(ICON_WIDTH_PLUS_GAP);
        }
        menuItem.setBackground(Color.white);
        menuItem.setMargin(new Insets(2, 2, 2, 2));
        menuItem.setEnabled(false);
        return menuItem;

    }

    void updatePluginsWithGlobalConfigParameters() {

        for (Plugin plugin : myMenuItemHash.values()) {
            ((AbstractPlugin) plugin).setConfigParameters();
        }

    }

    private JMenu getFiltersMenu() {
        JMenu result = new JMenu("Filter");
        result.setMnemonic(KeyEvent.VK_F);
        result.add(createMenuItem(new FilterSiteBuilderPlugin(this, true)));
        result.add(createMenuItem(new FilterTaxaBuilderPlugin(this, true)));
        result.add(createMenuItem(new FilterTraitsPlugin(this, true)));
        result.addSeparator();
        result.add(createMenuItem("Not Recommended", false));
        result.add(createMenuItem(new FilterAlignmentPlugin(this, true)));
        result.add(createMenuItem(new FilterSiteNamePlugin(this, true)));
        result.add(createMenuItem(new FilterTaxaPropertiesPlugin(this, true)));
        return result;
    }

    private JMenu getDataMenu() {

        JMenu result = new JMenu("Data");
        result.setMnemonic(KeyEvent.VK_D);

        result.add(createMenuItem(new GetTaxaListPlugin(this, true)));
        result.add(createMenuItem(new GetPositionListPlugin(this, true)));
        result.add(createMenuItem(new SortTaxaAlphabeticallyPlugin(this, true)));
        result.add(createMenuItem(new SortGenotypeFilePlugin(this, true)));
        result.add(createMenuItem(new SynonymizerPlugin(this, true)));
        result.add(createMenuItem(new IntersectionAlignmentPlugin(this, true)));
        result.add(createMenuItem(new UnionAlignmentPlugin(this, true)));
        result.add(createMenuItem(new MergeGenotypeTablesPlugin(this, true)));
        result.add(createMenuItem(new SeparatePlugin(this, true)));
        result.add(createMenuItem(new HetsToUnknownPlugin(this, true)));
        result.add(createMenuItem(new IndelsToUnknownPlugin(this, true)));
        result.add(createMenuItem(new SetLowDepthGenosToMissingPlugin(this, true)));
        result.add(createMenuItem(new TransformDataPlugin(this, true)));
        result.add(createMenuItem(new NumericalGenotypePlugin(this, true)));
        result.add(createMenuItem(new GenosToABHPlugin(this, true)));
        result.add(createMenuItem(new ThinSitesByPositionPlugin(this, true)));
        result.add(createMenuItem(new ClusterGenotypesPlugin(this, true)));
        result.add(createMenuItem(new MaskGenotypePlugin(this, true)));
        result.add(createMenuItem(new FindInversionsPlugin(this, true)));
        result.add(createMenuItem(new CreateHybridGenotypesPlugin(this, true)));
        result.add(createMenuItem(new GenotypeSummaryPlugin(this, true)));

        return result;
    }

    private JMenu getImputeMenu() {

        JMenu result = new JMenu("Impute");
        result.setMnemonic(KeyEvent.VK_I);

        result.add(createMenuItem(new FILLINFindHaplotypesPlugin(this, true), false));
        result.add(createMenuItem(new FILLINImputationPlugin(this, true), false));
        result.add(createMenuItem(new FSFHapImputationPlugin(this, true), false));
        result.add(createMenuItem(new ImputationPlugin(this, true), false));
        result.add(createMenuItem(new RemoveIndelsForBeaglePlugin(this, true), false));
        result.add(createMenuItem(new LDKNNiImputationPlugin(this, true), false));
        result.add(createMenuItem(new ImputationAccuracyPlugin(this, true), false));
        return result;
    }

    private JMenu getAnalysisMenu() {

        JMenu result = new JMenu("Analysis");
        result.setMnemonic(KeyEvent.VK_A);

        result.add(getPopGenMenu());
        result.add(getDistanceMenu());
        result.add(getAssociationMenu());

        return result;

    }

    private JMenu getPopGenMenu() {

        JMenu result = new JMenu("Diversity");
        result.setMnemonic(KeyEvent.VK_D);

        result.add(createMenuItem(new SequenceDiversityPlugin(this, true)));
        result.add(createMenuItem(new LinkageDisequilibriumPlugin(this, true)));
        return result;

    }

    private JMenu getDistanceMenu() {

        JMenu result = new JMenu("Relatedness");
        result.setMnemonic(KeyEvent.VK_R);

        result.add(createMenuItem(new DistanceMatrixPlugin(this, true)));
        result.add(createMenuItem(new KinshipPlugin(this, true)));
        result.add(createMenuItem(new CreateTreePlugin(this, true)));
        result.add(createMenuItem(new MergeTreesPlugin(this, true)));
        result.add(createMenuItem(new SubsetTreePlugin(this, true)));
        result.add(createMenuItem(new RemoveBranchLengthsPlugin(this, true)));
        result.add(createMenuItem(new ConvertTreeNamesPlugin(this, true)));
        result.add(createMenuItem(new AMatrixPlugin(this, true)));
        result.add(createMenuItem(new HMatrixPlugin(this, true)));
        result.add(createMenuItem(new MultiDimensionalScalingPlugin(this, true)));
        result.add(createMenuItem(new PrincipalComponentsPlugin(this, true)));
        result.addSeparator();
        result.add(createMenuItem(new RemoveNaNFromDistanceMatrixPlugin(this, true)));
        result.add(createMenuItem(new SubtractDistanceMatrixPlugin(this, true)));
        result.add(createMenuItem(new AddDistanceMatrixPlugin(this, true)));
        result.add(createMenuItem(new IdentityRecognitionPlugin(this, true)));
        return result;

    }

    private JMenu getAssociationMenu() {

        JMenu result = new JMenu("Association");
        result.setMnemonic(KeyEvent.VK_A);

        result.add(createMenuItem(new FixedEffectLMPlugin(this, true)));
        result.add(createMenuItem(new MLMPlugin(this, true)));
        result.add(createMenuItem(new WeightedMLMPlugin(this, true)));
        result.add(createMenuItem(new GenomicSelectionPlugin(this, true)));
        result.add(createMenuItem(new StepwiseOLSModelFitterPlugin(this, true)));
        result.add(createMenuItem(new FastMultithreadedAssociationPlugin(this, true)));
        result.add(createMenuItem(new MultivariateStepwisePlugin(this, true)));
        result.add(createMenuItem(new VCAPScanPlugin(this, true)));
        return result;

    }

    private JMenu getResultsMenu() {

        JMenu result = new JMenu("Results");
        result.setMnemonic(KeyEvent.VK_R);

        result.add(createMenuItem(new TableDisplayPlugin(this, true)));
        result.add(createMenuItem(new ArchaeopteryxPlugin(this, true)));
        result.add(createMenuItem(new LinkageDiseqDisplayPlugin(this, true)));
        result.add(createMenuItem(new ChartDisplayPlugin(this, true)));
        result.add(createMenuItem(new QQDisplayPlugin(this, true)));
        result.add(createMenuItem(new ManhattanDisplayPlugin(this, true)));
        return result;

    }

    private JMenu getFileMenu() {

        JMenu fileMenu = new JMenu();
        fileMenu.setText("File");

        FileLoadPlugin autoGuessPlugin = new FileLoadPlugin(this, true, true);
        autoGuessPlugin.addListener(myDataTreePanel);
        fileMenu.add(createMenuItem(new AbstractAction("Open", autoGuessPlugin.getIcon()) {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    autoGuessPlugin.processData(null);
                } catch (Exception ex) {
                    myLogger.debug(ex.getMessage(), ex);
                    DialogUtils.showError(ex.getMessage() + "\n", autoGuessPlugin.getParentFrame());
                }
            }
        }, KeyEvent.VK_O));

        fileMenu.add(createMenuItem(new FileLoadPlugin(this, true)));
        fileMenu.add(createMenuItem(new ExportPlugin(this, true)));

        URL deleteImageURL = TASSELMainFrame.class.getResource("/net/maizegenetics/analysis/images/trash.gif");
        Icon deleteIcon = null;
        if (deleteImageURL != null) {
            deleteIcon = new ImageIcon(deleteImageURL);
        }
        fileMenu.add(createMenuItem(new AbstractAction("Delete Dataset", deleteIcon) {
            @Override
            public void actionPerformed(ActionEvent e) {
                myDataTreePanel.deleteSelectedNodes();
            }
        }, -1));

        myPreferences = new PreferencesDialog(this, true);
        fileMenu.add(createMenuItem(myPreferences, true));
        fileMenu.add(createMenuItem(new ShowParameterCachePlugin(this, true), true));

        fileMenu.addSeparator();

        JMenuItem exitMenuItem = new JMenuItem("Exit");
        exitMenuItem.addActionListener((ActionEvent e) -> {
            TasselLogging.closeInstance();
            TasselPrefs.putXDim(getWidth());
            TasselPrefs.putYDim(getHeight());
            System.exit(0);
        });
        fileMenu.add(exitMenuItem);

        return fileMenu;

    }

    private JMenu getGBSv2Menu() {

        JMenu result = new JMenu("GBSv2");
        result.setMnemonic(KeyEvent.VK_G);

        result.add(createMenuItem(new GBSSeqToTagDBPlugin(this, true), false));
        result.add(createMenuItem(new TagExportToFastqPlugin(this, true), false));
        result.add(createMenuItem("Align to Reference", false));
        result.add(createMenuItem(new SAMToGBSdbPlugin(this, true), false));
        result.add(createMenuItem(new DiscoverySNPCallerPluginV2(this, true), false));
        result.add(createMenuItem(new SNPQualityProfilerPlugin(this, true), false));
        result.add(createMenuItem(new UpdateSNPPositionQualityPlugin(this, true), false));
        result.addSeparator();
        result.add(createMenuItem(new ProductionSNPCallerPluginV2(this, true), false));
        result.addSeparator();
        result.add(createMenuItem(new GetTagSequenceFromDBPlugin(this, true), false));
        result.addSeparator();
        result.add(createMenuItem(new SNPCutPosTagVerificationPlugin(this, true), false));

        return result;
    }

    private JMenu getGBSMenu() {

        JMenu result = new JMenu("GBS");
        Map attributes = result.getFont().getAttributes();
        attributes.put(TextAttribute.STRIKETHROUGH, TextAttribute.STRIKETHROUGH_ON);
        result.setFont(new Font(attributes));

        addMenuItemDeprecated(result, createMenuItem(new BinaryToTextPlugin(this, true), false));
        addMenuItemDeprecated(result, createMenuItem(new FastqToTagCountPlugin(this, true), false));
        addMenuItemDeprecated(result, createMenuItem(new MergeMultipleTagCountPlugin(this, true), false));
        result.addSeparator();
        addMenuItemDeprecated(result, getGBSReferenceMenu());
        result.addSeparator();
        addMenuItemDeprecated(result, getUNEAKMenu());
        result.addSeparator();
        addMenuItemDeprecated(result, createMenuItem(new SeqToTBTHDF5Plugin(this, true), false));
        addMenuItemDeprecated(result, createMenuItem(new ModifyTBTHDF5Plugin(this, true), false));
        addMenuItemDeprecated(result, createMenuItem(new DiscoverySNPCallerPlugin(this, true), false));
        addMenuItemDeprecated(result, createMenuItem(new ProductionSNPCallerPlugin(this, true), false));

        return result;
    }

    private JMenu getGBSReferenceMenu() {
        JMenu result = new JMenu("Reference Genome");
        addMenuItemDeprecated(result, createMenuItem(new TagCountToFastqPlugin(this, true), false));
        addMenuItemDeprecated(result, createMenuItem("Align to Reference", false));
        addMenuItemDeprecated(result, createMenuItem(new SAMConverterPlugin(this, true), false));
        return result;
    }

    private JMenu getUNEAKMenu() {
        JMenu result = new JMenu("UNEAK (No Reference)");
        addMenuItemDeprecated(result, createMenuItem(new UTagCountToTagPairPlugin(this, true), false));
        addMenuItemDeprecated(result, createMenuItem(new UTagPairToTOPMPlugin(this, true), false));
        return result;
    }

    private void addMenuItemDeprecated(JMenu menu, JMenuItem item) {
        Map attributes = item.getFont().getAttributes();
        attributes.put(TextAttribute.STRIKETHROUGH, TextAttribute.STRIKETHROUGH_ON);
        item.setFont(new Font(attributes));
        menu.add(item);
    }

    private void addMenuItemDeprecated(JMenu menu, JMenu item) {
        Map attributes = item.getFont().getAttributes();
        attributes.put(TextAttribute.STRIKETHROUGH, TextAttribute.STRIKETHROUGH_ON);
        item.setFont(new Font(attributes));
        menu.add(item);
    }

    private JMenu getWorkflowMenu() {

        JMenu result = new JMenu("Workflow");

        List<WorkflowPlugin> workflows = WorkflowPlugin.getInstances(this);

        for (Plugin current : workflows) {
            result.add(createMenuItem(current));
        }

        return result;
    }

    private JMenu getHelpMenu() {
        JMenu helpMenu = new JMenu();
        helpMenu.setMnemonic(KeyEvent.VK_H);
        helpMenu.setText("Help");

        URL infoImageURL = TASSELMainFrame.class.getResource("/net/maizegenetics/analysis/images/info.gif");
        Icon infoIcon = null;
        if (infoImageURL != null) {
            infoIcon = new ImageIcon(infoImageURL);
        }
        helpMenu.add(createMenuItem(new AbstractAction("Help Manual", infoIcon) {

            @Override
            public void actionPerformed(ActionEvent e) {
                final String html = "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Home";
                try {
                    Desktop desktop = Desktop.getDesktop();
                    URI uri = new URI(html);
                    desktop.browse(uri);
                } catch (Exception ex) {
                    myLogger.warn("Problem showing Tassel Wiki URl in Browser.", ex);
                }
            }
        }, -1));

        URL aboutImageURL = TASSELMainFrame.class.getResource("/net/maizegenetics/analysis/images/Tassel_Logo16.png");
        Icon aboutIcon = null;
        if (aboutImageURL != null) {
            aboutIcon = new ImageIcon(aboutImageURL);
        }
        helpMenu.add(createMenuItem(new AbstractAction("About", aboutIcon) {

            @Override
            public void actionPerformed(ActionEvent e) {
                helpAbout_actionPerformed(e);
            }
        }, -1));

        helpMenu.add(createMenuItem(PrintHeapAction.getInstance(this), -1));
        helpMenu.add(createMenuItem(TasselLogging.getInstance(this)));

        return helpMenu;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        JMenuItem theMenuItem = (JMenuItem) e.getSource();
        Plugin theTP = this.myMenuItemHash.get(theMenuItem);
        PluginEvent event = new PluginEvent(myDataTreePanel.getSelectedTasselDataSet());
        ProgressPanel progressPanel = getProgressPanel();
        progressPanel.addPlugin(theTP);
        ThreadedPluginListener thread = new ThreadedPluginListener(theTP, event);
        thread.start();
    }
}
