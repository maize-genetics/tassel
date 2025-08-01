/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General
 * public license.
 *
 */
package net.maizegenetics.tassel;

import net.maizegenetics.analysis.association.EqtlAssociationPlugin;
import net.maizegenetics.analysis.association.FastMultithreadedAssociationPlugin;
import net.maizegenetics.analysis.association.FixedEffectLMPlugin;
import net.maizegenetics.analysis.association.MLMPlugin;
import net.maizegenetics.analysis.data.GenotypeSummaryPlugin;
import net.maizegenetics.analysis.popgen.LinkageDisequilibriumPlugin;
import net.maizegenetics.analysis.popgen.SequenceDiversityPlugin;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.TOPMInterface;
import net.maizegenetics.dna.snp.FilterList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.gui.GenotypeTableMask;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.Plugin;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.PluginListener;
import net.maizegenetics.taxa.IdentifierSynonymizer;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.tree.SimpleTree;
import net.maizegenetics.taxa.tree.Tree;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.HDF5TableReport;
import net.maizegenetics.util.TableReport;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import javax.swing.event.TreeModelEvent;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.MutableTreeNode;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;
import java.awt.*;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.Serializable;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class DataTreePanel extends JPanel implements PluginListener {

    private static final Logger myLogger = LogManager.getLogger(DataTreePanel.class);
    public static final String NODE_TYPE_DATA = "Data";
    public static final String NODE_TYPE_RESULT = "Result";
    public static final String NODE_TYPE_SEQUENCE = "Sequence";
    public static final String NODE_TYPE_POLYMORPHISMS = "Polymorphisms";
    public static final String NODE_TYPE_NUMERICAL = "Numerical";
    public static final String NODE_TYPE_HDF5_SCHEMA = "HDF5 Schema";
    public static final String NODE_TYPE_MATRIX = "Matrix";
    public static final String NODE_TYPE_TREE = "Tree";
    public static final String NODE_TYPE_LISTS = "Lists";
    public static final String NODE_TYPE_FILTERS = "Filters";
    public static final String NODE_TYPE_FUSIONS = "Fusions";
    public static final String NODE_TYPE_SYNONYMS = "Synonyms";
    public static final String NODE_TYPE_DIVERSITY = "Diversity";
    public static final String NODE_TYPE_SNP_ASSAYS = "SNP Assays";
    public static final String NODE_TYPE_LD = "LD";
    public static final String NODE_TYPE_ASSOCIATIONS = "Association";
    public static final String NODE_TYPE_VARIANCES = "Variances";
    public static final String NODE_TYPE_SYNONYMIZER = "Synonymizer";
    public static final String NODE_TYPE_STEPWISE = "Stepwise";
    public static final String NODE_TYPE_TOPM = "TOPM";
    public static final String NODE_TYPE_GENO_SUMMARY = "Genotype Summary";
    public static final String NODE_TYPE_DEFAULT = NODE_TYPE_DATA;
    private static final List<String> NODE_TYPE_DATA_CHILDREN = new ArrayList<>();

    static {
        NODE_TYPE_DATA_CHILDREN.add(NODE_TYPE_SEQUENCE);
        NODE_TYPE_DATA_CHILDREN.add(NODE_TYPE_POLYMORPHISMS);
        NODE_TYPE_DATA_CHILDREN.add(NODE_TYPE_NUMERICAL);
        NODE_TYPE_DATA_CHILDREN.add(NODE_TYPE_HDF5_SCHEMA);
        NODE_TYPE_DATA_CHILDREN.add(NODE_TYPE_MATRIX);
        NODE_TYPE_DATA_CHILDREN.add(NODE_TYPE_TREE);
        NODE_TYPE_DATA_CHILDREN.add(NODE_TYPE_FUSIONS);
        NODE_TYPE_DATA_CHILDREN.add(NODE_TYPE_TOPM);
        NODE_TYPE_DATA_CHILDREN.add(NODE_TYPE_LISTS);
        NODE_TYPE_DATA_CHILDREN.add(NODE_TYPE_FILTERS);
    }
    //Possible line styles...
    //"Angled", "Horizontal", and "None" (the default).
    private String myLineStyle = "Angled";
    private JTree myTree;
    private DefaultTreeModel myTreeModel;
    private TASSELMainFrame myTASSELMainFrame;
    private HashMap<String, DefaultMutableTreeNode> myNodeHash = new HashMap<>();
    private LinkedHashMap myDataSetList = new LinkedHashMap();
    private Datum myLastBookSelected;
    private DefaultMutableTreeNode myTopNode;
    private DefaultMutableTreeNode myDataNode;
    private DefaultMutableTreeNode myResultNode;
    private DefaultMutableTreeNode myLastNode = null;
    private TreeSelectionListener myTreeSelectionListener = null;

    public DataTreePanel(TASSELMainFrame theQAF) {
        super();

        myTASSELMainFrame = theQAF;

        myTopNode = new DefaultMutableTreeNode("top");
        createNodes();

        myTreeModel = new DefaultTreeModel(myTopNode);

        myTree = new JTree(myTreeModel);
        myTree.setRootVisible(false);
        myTree.setEditable(true);
        myTreeModel.addTreeModelListener(new MyTreeModelListener(this));
        myTree.getSelectionModel().setSelectionMode(TreeSelectionModel.DISCONTIGUOUS_TREE_SELECTION);
        myTree.setLargeModel(true);
        //Listen for when the selection changes.
        initSelectionListener();
        initKeyStrokeListener();

        myTree.putClientProperty("JTree.lineStyle", myLineStyle);

        myTree.setCellRenderer(new DefaultTreeCellRenderer() {
            public Component getTreeCellRendererComponent(JTree pTree,
                    Object pValue, boolean pIsSelected, boolean pIsExpanded,
                    boolean pIsLeaf, int pRow, boolean pHasFocus) {
                DefaultMutableTreeNode node = (DefaultMutableTreeNode) pValue;
                Datum nodeInfo = (Datum) node.getUserObject();
                Component result = super.getTreeCellRendererComponent(pTree, pValue, pIsSelected,
                        pIsExpanded, pIsLeaf, pRow, pHasFocus);
                Object data = nodeInfo.getData();
                if (data instanceof GenotypeTableMask) {
                    result.setForeground(((GenotypeTableMask) nodeInfo.getData()).getColor());
                }

                return result;
            }
        });

        try {
            jbInit();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public DataTreePanel getThis() {
        return this;
    }

    private void processMyKeyStroke(KeyEvent e) {

        // only respond if the KeyEvent is a key released event.
        if (e.getID() == KeyEvent.KEY_RELEASED) {
            if (KeyEvent.VK_DELETE == e.getKeyCode()) {

                TreePath[] currentSelection = myTree.getSelectionPaths();
                String pluralQuestion = "Are you sure you want to delete these " + currentSelection.length + " nodes?";
                String question = "Are you sure you want to delete this node?";

                if (currentSelection.length > 1) {
                    question = pluralQuestion;
                }
                int userResponse = JOptionPane.showOptionDialog(myTASSELMainFrame, question, "Node Deletion",
                        JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);

                if (userResponse == JOptionPane.YES_OPTION) {
                    deleteSelectedNodes();
                }
            }
        }
    }

    private void initKeyStrokeListener() {

        myTree.addKeyListener(new KeyAdapter() {
            /**
             * Invoked when a key has been released.
             */
            public void keyReleased(KeyEvent e) {
                super.keyReleased(e);
                processMyKeyStroke(e);
            }
        });
    }

    private void initSelectionListener() {

        myTreeSelectionListener = (new TreeSelectionListener() {
            @Override
            public void valueChanged(TreeSelectionEvent e) {
                try {
                    DefaultMutableTreeNode node = null;
                    if ((e != null) && (e.getSource() instanceof DefaultMutableTreeNode)) {
                        node = (DefaultMutableTreeNode) e.getSource();
                    } else {
                        node = (DefaultMutableTreeNode) myTree.getLastSelectedPathComponent();
                    }

                    if (node == null) {
                        return;
                    }
                    if (myLastNode == node) {
                        return;
                    }
                    myLastNode = node;
                    Object nodeInfo = node.getUserObject();
                    myLastBookSelected = (Datum) nodeInfo;

                    Datum book = (Datum) nodeInfo;
                    myLogger.info("initSelectionListener: node type: " + book.getDataType());
                    StringBuilder builder = new StringBuilder();
                    if (book.getData() instanceof GenotypeTable) {
                        GenotypeTable a = (GenotypeTable) book.getData();
                        builder.append("Number of taxa: ");
                        builder.append(a.numberOfTaxa());
                        builder.append("\n");
                        builder.append("Number of sites: ");
                        builder.append(a.numberOfSites());
                        builder.append("\n");
                        Chromosome[] loci = a.chromosomes();
                        PositionList positions = null;
                        boolean first = true;
                        int numLoci = 0;
                        if (loci != null) {
                            numLoci = loci.length;
                            positions = a.positions();
                        }
                        for (int i = 0; i < numLoci; i++) {
                            String name = loci[i].getName();
                            if ((name == null) || (name.length() == 0)) {
                                continue;
                            }
                            if (first) {
                                builder.append("\nChromosomes...\n\n");
                                first = false;
                            }
                            builder.append(name);
                            int[] firstLast = a.firstLastSiteOfChromosome(loci[i]);
                            builder.append(": ");
                            builder.append(firstLast[1] - firstLast[0] + 1);
                            builder.append(" sites:\n");
                            builder.append(firstLast[0]);
                            builder.append(" (");
                            builder.append(positions.get(firstLast[0]).getPosition());
                            builder.append(") - ");
                            builder.append(firstLast[1]);
                            builder.append(" (");
                            builder.append(positions.get(firstLast[1]).getPosition());
                            builder.append(")");
                            builder.append("\n\n");
                        }
                    }
                    if (book.getData() instanceof GenotypeTable) {
                        try {
                            GenotypeTable a = (GenotypeTable) book.getData();
                            if ((a.hasGenotype()) && (NucleotideAlignmentConstants.isNucleotideEncodings(a.alleleDefinitions()))) {
                                builder.append("\n");
                                builder.append("Nucleotide Codes\n(Derived from IUPAC)...\n");
                                builder.append("A     A:A\n");
                                builder.append("C     C:C\n");
                                builder.append("G     G:G\n");
                                builder.append("T     T:T\n");
                                builder.append("R     A:G\n");
                                builder.append("Y     C:T\n");
                                builder.append("S     C:G\n");
                                builder.append("W     A:T\n");
                                builder.append("K     G:T\n");
                                builder.append("M     A:C\n");
                                builder.append("+     +:+ (insertion)\n");
                                builder.append("0     +:-\n");
                                builder.append("-     -:- (deletion)\n");
                                builder.append("N     Unknown\n\n");
                            }
                        } catch (Exception ex) {
                            myLogger.debug(ex.getMessage(), ex);
                        }

                    }
                    if (book.getData() instanceof TableReport) {
                        TableReport tr = (TableReport) book.getData();
                        long numColumns = tr.getColumnCount();
                        long numRows = tr.getRowCount();
                        builder.append("Table Title: ");
                        builder.append(tr.getTableTitle());
                        builder.append("\n");
                        builder.append("Number of columns: ");
                        builder.append(numColumns);
                        builder.append("\n");
                        builder.append("Number of rows: ");
                        builder.append(numRows);
                        builder.append("\n");
                        builder.append("Matrix size (excludes row headers): ");
                        builder.append((numColumns - 1) * numRows);
                        builder.append("\n");
                    }
                    if (book.getData() instanceof TOPMInterface) {
                        TOPMInterface topm = (TOPMInterface) book.getData();
                        builder.append("Number of Tags: ");
                        builder.append(topm.getSize());
                        builder.append("\n");
                    }
                    if (book.getData() instanceof DistanceMatrix) {
                        DistanceMatrix matrix = (DistanceMatrix) book.getData();
                        GeneralAnnotation annotations = matrix.annotations();
                        if ((annotations != null) && (annotations.numAnnotations() != 0)) {
                            builder.append("\nAnnotations...\n");
                            for (Map.Entry<String, String> entry : annotations.getAllAnnotationEntries()) {
                                builder.append(entry.getKey());
                                builder.append(": ");
                                builder.append(entry.getValue());
                                builder.append("\n");
                            }
                            builder.append("\n");
                        }
                    }
                    if (book.getData() instanceof TaxaList) {
                        TaxaList taxa = (TaxaList) book.getData();
                        builder.append("Number of Taxa: ");
                        builder.append(taxa.numberOfTaxa());
                        builder.append("\n\n");
                    }

                    String comment = book.getComment();
                    if ((comment != null) && (comment.length() != 0)) {
                        builder.append(comment);
                    }
                    myTASSELMainFrame.setNoteText(builder.toString());
                    if (book.getData() != null) {
                        myTASSELMainFrame.mainDisplayPanel.removeAll();

                        if (book.getData() instanceof TableReport) {
                            long size = ((TableReport) book.getData()).getElementCount();
                            myLogger.info("initSelectionListener: Table Report Size: " + size);
                            if (size == 0) {
                                JPanel blankPanel = new JPanel();
                                blankPanel.setLayout(new BorderLayout());
                                blankPanel.add(new JLabel("     Nothing to Display"), BorderLayout.CENTER);
                                myTASSELMainFrame.mainDisplayPanel.add(blankPanel, BorderLayout.CENTER);
                            } else {
                                TableReportPanel theATP = TableReportPanel.getInstance((TableReport) book.getData());
                                myTASSELMainFrame.mainDisplayPanel.add(theATP, BorderLayout.CENTER);
                            }
                        } else if (book.getData() instanceof FilterList) {
                            TableReportPanel theATP = TableReportPanel.getInstance((FilterList) book.getData());
                            myTASSELMainFrame.mainDisplayPanel.add(theATP, BorderLayout.CENTER);
                        } else if (book.getData() instanceof TOPMInterface) {
                            int size = ((TOPMInterface) book.getData()).getTagCount();
                            myLogger.info("initSelectionListener: TOPM Tag Count: " + size);
                            if (size == 0) {
                                JPanel blankPanel = new JPanel();
                                blankPanel.setLayout(new BorderLayout());
                                blankPanel.add(new JLabel("     No Tags"), BorderLayout.CENTER);
                                myTASSELMainFrame.mainDisplayPanel.add(blankPanel, BorderLayout.CENTER);
                            } else {
                                String name = book.getName();
                                if (name.endsWith("(Graphical)")) {
                                    SeqViewerPanel seqViewer = SeqViewerPanel.getInstance((TOPMInterface) book.getData(), getThis());
                                    myTASSELMainFrame.mainDisplayPanel.add(seqViewer, BorderLayout.CENTER);
                                } else {
                                    TableReportPanel theATP = TableReportPanel.getInstance((TOPMInterface) book.getData());
                                    myTASSELMainFrame.mainDisplayPanel.add(theATP, BorderLayout.CENTER);
                                }

                            }
                        } else if (book.getData() instanceof TaxaList) {
                            int size = ((TaxaList) book.getData()).numberOfTaxa();
                            myLogger.info("initSelectionListener: Number of Taxa: " + size);
                            if (size == 0) {
                                JPanel blankPanel = new JPanel();
                                blankPanel.setLayout(new BorderLayout());
                                blankPanel.add(new JLabel("     No Taxa"), BorderLayout.CENTER);
                                myTASSELMainFrame.mainDisplayPanel.add(blankPanel, BorderLayout.CENTER);
                            } else {
                                TableReportPanel theATP = TableReportPanel.getInstance((TaxaList) book.getData());
                                myTASSELMainFrame.mainDisplayPanel.add(theATP, BorderLayout.CENTER);
                            }
                        } else if (book.getData() instanceof PositionList) {
                            int size = ((PositionList) book.getData()).numberOfSites();
                            myLogger.info("initSelectionListener: Number of Positions: " + size);
                            if (size == 0) {
                                JPanel blankPanel = new JPanel();
                                blankPanel.setLayout(new BorderLayout());
                                blankPanel.add(new JLabel("     No Positions"), BorderLayout.CENTER);
                                myTASSELMainFrame.mainDisplayPanel.add(blankPanel, BorderLayout.CENTER);
                            } else {
                                TableReportPanel theATP = TableReportPanel.getInstance((PositionList) book.getData());
                                myTASSELMainFrame.mainDisplayPanel.add(theATP, BorderLayout.CENTER);
                            }
                        } else if (book.getData() instanceof GenotypeTable) {
                            GenotypeTable align = (GenotypeTable) book.getData();
                            List masks = new ArrayList();
                            for (int i = 0, n = node.getChildCount(); i < n; i++) {
                                try {
                                    DefaultMutableTreeNode currentNode = (DefaultMutableTreeNode) node.getChildAt(i);
                                    Datum currentDatum = (Datum) currentNode.getUserObject();
                                    Object currentMask = currentDatum.getData();
                                    if ((currentMask != null) && (currentMask instanceof GenotypeTableMask)) {
                                        masks.add(currentMask);
                                    }
                                } catch (Exception ex) {
                                    ex.printStackTrace();
                                }
                            }

                            SeqViewerPanel seqViewer = null;
                            if (masks.size() != 0) {
                                GenotypeTableMask[] masksArray = new GenotypeTableMask[masks.size()];
                                masks.toArray(masksArray);
                                seqViewer = SeqViewerPanel.getInstance(align, masksArray, getThis());
                            } else {
                                seqViewer = SeqViewerPanel.getInstance(align, getThis());
                            }
                            myTASSELMainFrame.mainDisplayPanel.add(seqViewer, BorderLayout.CENTER);
                        } else if (book.getData() instanceof GenotypeTableMask) {
                            GenotypeTableMask mask = (GenotypeTableMask) book.getData();
                            if (mask.getColor() != null) {
                                Color tempColor = JColorChooser.showDialog(myTASSELMainFrame, "Select Mask Color...", mask.getColor());
                                //myTree.setSelectionPath(new TreePath(parentNode.getPath()));
                                if (tempColor != null) {
                                    mask.setColor(tempColor);
                                    book.setName(mask.toString());
                                }
                            }
                            myLastNode = null;
                            DefaultMutableTreeNode parentNode = (DefaultMutableTreeNode) node.getParent();
                            TreeSelectionEvent event = new TreeSelectionEvent(parentNode, null, false, null, null);
                            myTreeSelectionListener.valueChanged(event);
                        } else if (book.getData() instanceof SimpleTree) {
                            myTASSELMainFrame.mainDisplayPanel.add(myTASSELMainFrame.mainPanelScrollPane, BorderLayout.CENTER);
                            SimpleTree tree = (SimpleTree) book.getData();
                            StringWriter writer = null;
                            try {
                                if (tree.getExternalNodeCount() < 300) {
                                    writer = new StringWriter();
                                    tree.report(writer);
                                    myTASSELMainFrame.setMainText(writer.toString());
                                } else {
                                    myTASSELMainFrame.setMainText(tree.toString());
                                }
                            } catch (Exception ex) {
                                myLogger.debug(ex.getMessage(), ex);
                            } finally {
                                try {
                                    writer.close();
                                } catch (Exception ex) {
                                    // do nothing
                                }
                            }
                        } else {
                            myTASSELMainFrame.mainDisplayPanel.add(myTASSELMainFrame.mainPanelScrollPane, BorderLayout.CENTER);
                            String s = book.getData().toString();
                            if (s.length() > 2000000) {
                                s = s.substring(1, 2000000) + "\n Truncated view.  Too much to display.  Save it to a file.";
                            }
                            myTASSELMainFrame.setMainText(s);
                        }

                        myTASSELMainFrame.sendMessage(book.getData().getClass().toString());

                        myTASSELMainFrame.mainDisplayPanel.revalidate();
                        myTASSELMainFrame.mainDisplayPanel.repaint();

                    }

                } catch (Throwable e1) {

                    String userMessage = "TASSEL has experienced an error.  ";
                    if (e1 instanceof java.lang.OutOfMemoryError) {
                        userMessage = "You have used up all of the memory allocated to the Java Virtual Machine.  \n"
                                + "It is recommneded that you adjust your heap settings and possibly add more memory to the computer\n"
                                + "Additionally, some operations are not recommended on a full dataset, i.e., select SNPs *before* determining LD";
                    }
                    JOptionPane.showMessageDialog(myTASSELMainFrame, userMessage, "Fatal Error", JOptionPane.ERROR_MESSAGE);
                    e1.printStackTrace();
                }
            }
        });

        myTree.addTreeSelectionListener(myTreeSelectionListener);

    }

    private void createNodes() {
        myDataNode = new DefaultMutableTreeNode(new Datum(NODE_TYPE_DATA, "Node on data tree", "Holds the basic data structures"));
        myTopNode.add(myDataNode);
        myNodeHash.put(NODE_TYPE_DATA, myDataNode);
        myResultNode = new DefaultMutableTreeNode(new Datum(NODE_TYPE_RESULT, "Node on data tree", "Holds the basic result structures"));
        myTopNode.add(myResultNode);
        myNodeHash.put(NODE_TYPE_RESULT, myResultNode);
    }

    private synchronized DefaultMutableTreeNode getTreeNode(String nodeString) {
        DefaultMutableTreeNode result = myNodeHash.get(nodeString);
        if (result == null) {
            if (NODE_TYPE_DATA_CHILDREN.contains(nodeString)) {
                result = new DefaultMutableTreeNode(new Datum(nodeString, nodeString, nodeString));
                int childIndex = getInsertLocation(myDataNode, nodeString);
                myTreeModel.insertNodeInto(result, myDataNode, childIndex);
                myNodeHash.put(nodeString, result);
            } else {
                result = new DefaultMutableTreeNode(new Datum(nodeString, nodeString, nodeString));
                int childIndex = getInsertLocation(myResultNode, nodeString);
                myTreeModel.insertNodeInto(result, myResultNode, childIndex);
                myNodeHash.put(nodeString, result);
            }
        }
        return result;
    }

    private int getInsertLocation(DefaultMutableTreeNode node, String newNode) {
        Enumeration<TreeNode> children = node.children();
        int result = 0;
        while (children.hasMoreElements()) {
            DefaultMutableTreeNode currentNode = (DefaultMutableTreeNode) children.nextElement();
            Datum current = (Datum) currentNode.getUserObject();
            if (current.getName().compareTo(newNode) > 0) {
                break;
            }
            result++;
        }
        return result;
    }

    private void jbInit() throws Exception {
        JScrollPane jScrollPane1 = new JScrollPane(myTree);
        jScrollPane1.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        setLayout(new BorderLayout());
        add(jScrollPane1, BorderLayout.CENTER);
    }

    private void deleteDatumFromList(Datum datum) {

        Iterator itr = myDataSetList.keySet().iterator();
        while (itr.hasNext()) {
            Datum current = (Datum) itr.next();
            if (current == datum) {
                itr.remove();
                break;
            }
        }

    }

    public void addDataSet(DataSet theDataSet, String defaultNode) {

        Plugin theCreator = theDataSet.getCreator();
        for (int i = 0; i < theDataSet.getSize(); i++) {
            Datum d = theDataSet.getData(i);
            if ((theCreator instanceof MLMPlugin)
                    || (theCreator instanceof FixedEffectLMPlugin)
                    || (theCreator instanceof FastMultithreadedAssociationPlugin)
                    || (theCreator instanceof EqtlAssociationPlugin)) {
                addDatum(NODE_TYPE_ASSOCIATIONS, d);
                continue;
            }

            if (theCreator instanceof GenotypeSummaryPlugin) {
                addDatum(NODE_TYPE_GENO_SUMMARY, d);
                continue;
            }

            if (theCreator instanceof SequenceDiversityPlugin) {
                addDatum(NODE_TYPE_DIVERSITY, d);
                continue;
            }

            if (theCreator instanceof LinkageDisequilibriumPlugin) {
                addDatum(NODE_TYPE_LD, d);
                continue;
            }

            if (d.getData() instanceof GenotypeTable) {
                addDatum(NODE_TYPE_SEQUENCE, d);
                continue;
            }

            if (d.getData() instanceof GenotypeTableMask) {
                addDatum(d);
                continue;
            }

            if (d.getData() instanceof IdentifierSynonymizer) {
                addDatum(NODE_TYPE_SYNONYMIZER, d);
                continue;
            }

            if (d.getData() instanceof Phenotype) {
                addDatum(NODE_TYPE_NUMERICAL, d);
                continue;
            }

            if (d.getData() instanceof DistanceMatrix) {
                addDatum(NODE_TYPE_MATRIX, d);
                continue;
            }

            if (d.getData() instanceof HDF5TableReport) {
                addDatum(NODE_TYPE_HDF5_SCHEMA, d);
                continue;
            }

            if (d.getData() instanceof TableReport) {
                addDatum(NODE_TYPE_NUMERICAL, d);
                continue;
            }

            if (d.getData() instanceof FilterList) {
                addDatum(NODE_TYPE_FILTERS, d);
                continue;
            }

            if (d.getData() instanceof Tree) {
                addDatum(NODE_TYPE_TREE, d);
                continue;
            }

            if (d.getData() instanceof TaxaList) {
                addDatum(NODE_TYPE_LISTS, d);
                continue;
            }

            if (d.getData() instanceof PositionList) {
                addDatum(NODE_TYPE_LISTS, d);
                continue;
            }

            if (d.getData() instanceof TOPMInterface) {
                addDatum(NODE_TYPE_TOPM, new Datum(d.getName() + " (Text)", d.getData(), null));
                addDatum(NODE_TYPE_TOPM, new Datum(d.getName() + " (Graphical)", d.getData(), null));
                continue;
            }

            if (defaultNode == null) {
                addDatum(NODE_TYPE_DEFAULT, d);
                continue;
            }

            addDatum(defaultNode, d);
        }

    }

    public synchronized void addDatum(String dataParent, Datum theDatum) {
        if (theDatum.getData() instanceof GenotypeTableMask) {
            addDatum(theDatum);
            return;
        }

        DefaultMutableTreeNode temp = findOnTree(theDatum);
        if (temp != null) {
            myTree.scrollPathToVisible(new TreePath(temp.getPath()));
            myTree.validate();
            myTree.repaint();
            return;
        }

        DefaultMutableTreeNode parentNode = getTreeNode(dataParent);
        DefaultMutableTreeNode childNode = new DefaultMutableTreeNode(theDatum);
        myNodeHash.put(theDatum.getName(), childNode);
        myTreeModel.insertNodeInto(childNode, parentNode, parentNode.getChildCount());
        myTree.scrollPathToVisible(new TreePath(childNode.getPath()));
        if (theDatum instanceof Serializable) {
            myDataSetList.put(theDatum, dataParent);
        }
    }

    public synchronized void addDatum(Datum theDatum) {
        GenotypeTableMask mask = null;
        try {
            mask = (GenotypeTableMask) theDatum.getData();
        } catch (Exception e) {
            throw new IllegalArgumentException("DataTreePanel: addDatum: input must be AlignmentMask.");
        }
        GenotypeTable align = mask.getAlignment();
        Iterator itr = myNodeHash.values().iterator();
        while (itr.hasNext()) {
            DefaultMutableTreeNode parentNode = (DefaultMutableTreeNode) itr.next();
            Object current = ((Datum) parentNode.getUserObject()).getData();
            if ((current instanceof GenotypeTable) && (current == align)) {
                DefaultMutableTreeNode childNode = new DefaultMutableTreeNode(theDatum);
                myTreeModel.insertNodeInto(childNode, parentNode, parentNode.getChildCount());
                myTree.scrollPathToVisible(new TreePath(childNode.getPath()));
                myLastNode = null;
                myTree.setSelectionPath(new TreePath(parentNode.getPath()));
                myTree.scrollPathToVisible(new TreePath(parentNode.getPath()));
                myTreeSelectionListener.valueChanged(null);
                if (theDatum instanceof Serializable) {
                    myDataSetList.put(theDatum, "NA");
                }
                break;
            }
        }
    }

    public void setSelectionPath(Datum theDatum) {
        DefaultMutableTreeNode temp = findOnTree(theDatum);
        if (temp != null) {
            myTree.scrollPathToVisible(new TreePath(temp.getPath()));
            myTree.setSelectionPath(new TreePath(temp.getPath()));
            myTree.validate();
            myTree.repaint();
        }
    }

    private DefaultMutableTreeNode findOnTree(Datum theDatum) {

        if (theDatum == null) {
            return null;
        }

        Iterator itr = myNodeHash.values().iterator();
        while (itr.hasNext()) {
            DefaultMutableTreeNode parentNode = (DefaultMutableTreeNode) itr.next();
            Datum current = (Datum) parentNode.getUserObject();
            if ((current != null) && (current.equals(theDatum))) {
                return parentNode;
            }
        }
        return null;

    }

    public DataSet getSelectedTasselDataSet() {
        int n = myTree.getSelectionCount();
        String parentNodeName = null;
        ArrayList<Datum> data = new ArrayList<Datum>();
        TreePath[] tp = myTree.getSelectionPaths();
        for (int i = 0; i < n; i++) {
            DefaultMutableTreeNode node = (DefaultMutableTreeNode) tp[i].getLastPathComponent();
            if (node == null) {
                return null;
            }
            //if the parents don't share the same name then leave blank

            if (i == 0) {
                parentNodeName = node.getParent().toString();
            } else if (!node.getParent().toString().equals(parentNodeName)) {
                parentNodeName = null;
            }

            data.add((Datum) node.getUserObject());
        }

        return new DataSet(data, null);
    }

    public void deleteSelectedNodes() {

        TreePath[] currentSelection = myTree.getSelectionPaths();
        if (currentSelection != null) {
            for (int i = 0; i < currentSelection.length; i++) {
                myLogger.info("Start Deleting at Selection: " + currentSelection[i]);
                if (currentSelection[i] != null) {
                    DefaultMutableTreeNode currentNode = (DefaultMutableTreeNode) (currentSelection[i].getLastPathComponent());
                    myLogger.info("Deleting Node : " + currentNode);
                    deleteNode(currentNode);
                }

            }
        }

    }

    private void deleteNode(DefaultMutableTreeNode currentNode) {

        int size = currentNode.getChildCount();
        for (int i = size - 1; i >= 0; i--) {
            DefaultMutableTreeNode current = (DefaultMutableTreeNode) currentNode.getChildAt(i);
            deleteNode(current);
        }

        if ((currentNode == myDataNode) || (currentNode == myResultNode)) {
            return;
        }

        MutableTreeNode parent = (MutableTreeNode) (currentNode.getParent());
        if ((parent != null) && (parent != myDataNode) && (parent != myResultNode)) {
            myTreeModel.removeNodeFromParent(currentNode);
            Iterator itr = myNodeHash.keySet().iterator();
            while (itr.hasNext()) {
                Object key = itr.next();
                if (myNodeHash.get(key) == currentNode) {
                    itr.remove();
                }

            }

            try {
                Datum datum = (Datum) currentNode.getUserObject();
                deleteDatumFromList(datum);
            } catch (Exception ex) {
                // do nothing
            }

            try {
                if (((Datum) currentNode.getUserObject()).getData() instanceof GenotypeTableMask) {
                    myLastNode = null;
                    myTree.setSelectionPath(new TreePath(((DefaultMutableTreeNode) parent).getPath()));
                    myTree.scrollPathToVisible(new TreePath(((DefaultMutableTreeNode) parent).getPath()));
                    myTreeSelectionListener.valueChanged(null);
                }
            } catch (Exception e) {
                // do nothing
            }

            try {
                if (((Datum) currentNode.getUserObject()).getData() instanceof GenotypeTable) {
                    SeqViewerPanel.removeInstance((GenotypeTable) ((Datum) currentNode.getUserObject()).getData());
                }
            } catch (Exception e) {
                // do nothing
            }
        }

    }

    public Map getDataList() {
        return Collections.unmodifiableMap(myDataSetList);
    }

    public JTree getTree() {
        return myTree;
    }

    /**
     * Returns Tassel data set after complete.
     *
     * @param event event
     */
    public void dataSetReturned(PluginEvent event) {
        DataSet tds = (DataSet) event.getSource();
        if (tds != null) {
            addDataSet(tds, DataTreePanel.NODE_TYPE_DEFAULT);
            setSelectionPath(tds.getData(0));
        }
    }

    /**
     * Returns progress of execution.
     *
     * @param event event
     */
    public void progress(PluginEvent event) {
        //do nothing
    }

    class MyTreeModelListener implements javax.swing.event.TreeModelListener {

        DataTreePanel theDTP;

        public MyTreeModelListener(DataTreePanel theDTP) {
            this.theDTP = theDTP;
        }

        public void treeNodesChanged(TreeModelEvent e) {
            DefaultMutableTreeNode node;
            node = (DefaultMutableTreeNode) (e.getTreePath().getLastPathComponent());

            /*
             * If the event lists children, then the changed
             * node is the child of the node we've already
             * gotten.  Otherwise, the changed node and the
             * specified node are the same.
             */
            try {
                int index = e.getChildIndices()[0];
                node = (DefaultMutableTreeNode) (node.getChildAt(index));
            } catch (NullPointerException exc) {
            }
            String s = (String) node.getUserObject();
            theDTP.myLastBookSelected.setName(s);
            node.setUserObject(theDTP.myLastBookSelected);
            System.out.println("The user has finished editing the node.");
            System.out.println("New value: " + node.getUserObject());
        }

        public void treeNodesInserted(TreeModelEvent e) {
        }

        public void treeNodesRemoved(TreeModelEvent e) {
        }

        public void treeStructureChanged(TreeModelEvent e) {
        }
    }
}
