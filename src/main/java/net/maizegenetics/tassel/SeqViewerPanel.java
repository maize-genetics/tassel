/*
 * SeqViewerPanel
 */
package net.maizegenetics.tassel;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import java.net.URL;

import java.util.Map;
import java.util.WeakHashMap;

import javax.swing.AbstractAction;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.event.MouseInputListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.plaf.basic.BasicTableUI;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import net.maizegenetics.dna.map.TOPMInterface;
import net.maizegenetics.dna.snp.FilterGenotypeTable;

import net.maizegenetics.gui.AlignmentTableCellRenderer;
import net.maizegenetics.gui.AlignmentTableModel;
import net.maizegenetics.gui.RowHeaderRenderer;
import net.maizegenetics.gui.TableRowHeaderListModel;
import net.maizegenetics.gui.VerticalLabelUI;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.TOPMGenotypeTable;
import net.maizegenetics.gui.GenotypeTableMask;
import net.maizegenetics.gui.GenotypeTableMaskGeneticDistance;
import net.maizegenetics.gui.GenotypeTableMaskReference;
import net.maizegenetics.gui.TOPMGenotypeTableCellRenderer;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.Tuple;

/**
 *
 * @author Terry Casstevens
 */
public class SeqViewerPanel extends JPanel implements ComponentListener, TableModelListener {

    private static final Map<Object, Object[]> INSTANCES = new WeakHashMap<>();
    private static final int ROW_HEADER_WIDTH = 150;
    private static final int SCROLL_BAR_WIDTH = 25;
    private static final int SLIDER_TEXT_WIDTH = 125;
    private JSlider mySlider;
    private JButton myLeftButton;
    private JButton myRightButton;
    private final JTable myTable;
    private final AlignmentTableModel myTableModel;
    private final GenotypeTable myAlignment;
    private final JScrollPane myScrollPane;
    private JPanel mySliderPane;
    private JLabel searchLabel;
    private JTextField searchField;
    private JButton searchButton;
    private int start;
    private int end;
    private int startPos;
    private int endPos;
    private int siteCount;
    private final DataTreePanel myDataTreePanel;
    private final AlignmentTableCellRenderer myTableCellRenderer;
    private final JComboBox<AlignmentTableCellRenderer.RENDERING_TYPE> myHighlightingComboBox;

    private SeqViewerPanel(GenotypeTable alignment, GenotypeTableMask[] masks, DataTreePanel dataTreePanel) {
        this(alignment, masks, dataTreePanel, -1);
    }

    private SeqViewerPanel(TOPMGenotypeTable alignment, DataTreePanel dataTreePanel) {
        this(alignment, null, dataTreePanel, 0);
    }

    private SeqViewerPanel(GenotypeTable alignment, GenotypeTableMask[] masks, DataTreePanel dataTreePanel, int sliderPosition) {

        setLayout(new BorderLayout());
        myAlignment = alignment;
        myDataTreePanel = dataTreePanel;

        myTableModel = new AlignmentTableModel(alignment);

        if (myAlignment instanceof TOPMGenotypeTable) {
            myTableCellRenderer = new TOPMGenotypeTableCellRenderer(myTableModel, (TOPMGenotypeTable) myAlignment);
        } else {
            myTableCellRenderer = new AlignmentTableCellRenderer(myTableModel, myAlignment, masks);
        }

        myHighlightingComboBox = new JComboBox<>(myTableCellRenderer.getRenderingTypes());
        myHighlightingComboBox.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                AlignmentTableCellRenderer.RENDERING_TYPE type = (AlignmentTableCellRenderer.RENDERING_TYPE) myHighlightingComboBox.getSelectedItem();
                myTableCellRenderer.setRenderingType(type);
                myTableModel.setHorizontalPageSize(calculateNumColumns());
                myTableModel.fireTableChanged();
            }
        });
        myHighlightingComboBox.setSelectedItem(myTableCellRenderer.getRenderingType());

        siteCount = myAlignment.numberOfSites();
        start = 0;
        end = siteCount - 1;
        startPos = myAlignment.chromosomalPosition(0);
        endPos = myAlignment.chromosomalPosition(end);

        mySlider = new JSlider();
        mySlider.addChangeListener(myTableModel);

        if (sliderPosition == -1) {
            myTableModel.adjustPositionToCenter();
        } else {
            myTableModel.adjustPositionToSite(sliderPosition);
        }
        myTableModel.addTableModelListener(this);
        myTable = new JTable(myTableModel);
        myTable.setUI(new MyTableUI());
        myTable.setDefaultRenderer(myTable.getColumnClass(0), myTableCellRenderer);
        myTable.setAutoResizeMode(JTable.AUTO_RESIZE_SUBSEQUENT_COLUMNS);
        myTable.getTableHeader().setReorderingAllowed(false);
        JList rowHeaders = new JList(new TableRowHeaderListModel(myTableModel.getRowHeaders())) {
            public String getToolTipText(MouseEvent evt) {

                int index = locationToIndex(evt.getPoint());
                Taxon id = (Taxon) getModel().getElementAt(index);

                return id.getName();

            }
        };

        rowHeaders.setFixedCellWidth(ROW_HEADER_WIDTH);
        rowHeaders.setFixedCellHeight(myTable.getRowHeight());
        rowHeaders.setCellRenderer(new RowHeaderRenderer(myTable));

        myScrollPane = new JScrollPane(myTable);
        myScrollPane.addComponentListener(this);
        myScrollPane.setRowHeaderView(rowHeaders);

        add(getControls(), BorderLayout.NORTH);
        add(myScrollPane, BorderLayout.CENTER);

        boolean multipleAlignments = myAlignment.numChromosomes() > 1;

        if (multipleAlignments) {
            myTableModel.setColumnNameType(AlignmentTableModel.COLUMN_NAME_TYPE.siteNumber);
            updateSliderSiteNumbers();
        } else {
            myTableModel.setColumnNameType(AlignmentTableModel.COLUMN_NAME_TYPE.physicalPosition);
            updateSliderPhysicalPositions();
        }

        initMenu();
    }

    public static SeqViewerPanel getInstance(GenotypeTable alignment, GenotypeTableMask[] masks, DataTreePanel dataTreePanel) {
        Object[] instance = (Object[]) INSTANCES.get(alignment);
        SeqViewerPanel result = null;
        if (instance == null) {
            result = new SeqViewerPanel(alignment, masks, dataTreePanel);
            saveInstance(result, alignment, masks);
        } else {
            result = (SeqViewerPanel) instance[0];
            result.setMasks(masks);
            saveInstance(result, alignment, masks);
        }
        return result;
    }

    public static SeqViewerPanel getInstance(TOPMInterface topm, DataTreePanel dataTreePanel) {
        Object[] instance = (Object[]) INSTANCES.get(topm);
        SeqViewerPanel result = null;
        if (instance == null) {
            result = new SeqViewerPanel(new TOPMGenotypeTable(topm), dataTreePanel);
            saveInstance(result, topm, null);
        } else {
            result = (SeqViewerPanel) instance[0];
        }
        return result;
    }

    private static void saveInstance(SeqViewerPanel panel, Object alignment, GenotypeTableMask[] masks) {
        int arraySize = 1;
        if (masks != null) {
            arraySize += masks.length;
        }
        Object[] instance = new Object[arraySize];
        instance[0] = panel;
        if (masks != null) {
            for (int i = 0; i < masks.length; i++) {
                instance[i + 1] = masks[i];
            }
        }
        INSTANCES.put(alignment, instance);
    }

    public static void removeInstance(GenotypeTable alignment) {
        try {
            INSTANCES.remove(alignment);
        } catch (Exception e) {
            // do nothing
        }
    }

    public static SeqViewerPanel getInstance(GenotypeTable alignment, DataTreePanel dataTreePanel) {
        return getInstance(alignment, null, dataTreePanel);
    }

    public void setMasks(GenotypeTableMask[] masks) {
        myTableCellRenderer.setMasks(masks);
    }

    private void initMenu() {

        JPopupMenu menu = new JPopupMenu();
        menu.setInvoker(this);

        JMenuItem useAsReference = new JMenuItem("Use this Taxa as Reference");
        useAsReference.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                int index = myTable.getSelectedRow();
                GenotypeTableMaskReference mask = GenotypeTableMaskReference.getInstanceCompareReference(myAlignment, index);
                myHighlightingComboBox.setSelectedItem(AlignmentTableCellRenderer.RENDERING_TYPE.ReferenceMasks);
                myDataTreePanel.addDatum(new Datum(mask.toString(), mask, null));
            }
        });
        menu.add(useAsReference);

        JMenuItem useAsGeneticDistance = new JMenuItem("Use this Taxa for Genetic Distance");
        useAsGeneticDistance.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                int index = myTable.getSelectedRow();
                Tuple<GenotypeTable, double[]> sortedByGeneticDistance = FilterGenotypeTable.getInstanceTaxaOrderedByGeneticDistance(myAlignment, index);
                DataSet dataSet = new DataSet(new Datum("Genetic Distance with " + myAlignment.taxaName(index), sortedByGeneticDistance.x, null), null);
                myDataTreePanel.addDataSet(dataSet, null);
                GenotypeTableMaskGeneticDistance mask = GenotypeTableMaskGeneticDistance.getInstanceCompareReference(sortedByGeneticDistance.x, 0, sortedByGeneticDistance.y);
                SeqViewerPanel newViewer = getInstance(sortedByGeneticDistance.x, new GenotypeTableMask[]{mask}, myDataTreePanel);
                newViewer.myHighlightingComboBox.setSelectedItem(AlignmentTableCellRenderer.RENDERING_TYPE.GeneticDistanceMasks);
                myDataTreePanel.addDatum(new Datum(mask.toString(), mask, null));
            }
        });
        menu.add(useAsGeneticDistance);

        myTable.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                checkPopup(e);
            }

            @Override
            public void mouseClicked(MouseEvent e) {
                checkPopup(e);
            }

            @Override
            public void mouseReleased(MouseEvent e) {
                checkPopup(e);
            }

            private void checkPopup(MouseEvent e) {
                if (e.isPopupTrigger()) {
                    menu.show(myTable, e.getX(), e.getY());
                }
            }
        });

    }

    private JPanel getControls() {

        JPanel result = new JPanel();
        result.setLayout(new BoxLayout(result, BoxLayout.Y_AXIS));

        RadioListener radioListener = new RadioListener();

        ButtonGroup buttonGroup = new ButtonGroup();

        boolean multipleAlignments = myAlignment.numChromosomes() > 1;

        JRadioButton physicalPosition = null;
        if (!multipleAlignments) {
            physicalPosition = new JRadioButton("Physical Positions");
            physicalPosition.setActionCommand(AlignmentTableModel.COLUMN_NAME_TYPE.physicalPosition.toString());
            physicalPosition.addActionListener(radioListener);
        }

        JRadioButton siteNumber = new JRadioButton("Site Numbers");
        siteNumber.setActionCommand(AlignmentTableModel.COLUMN_NAME_TYPE.siteNumber.toString());
        siteNumber.addActionListener(radioListener);

        JRadioButton locus = new JRadioButton("Locus");
        locus.setActionCommand(AlignmentTableModel.COLUMN_NAME_TYPE.locus.toString());
        locus.addActionListener(radioListener);

        JRadioButton siteName = new JRadioButton("Site Name");
        siteName.setActionCommand(AlignmentTableModel.COLUMN_NAME_TYPE.siteName.toString());
        siteName.addActionListener(radioListener);

        JRadioButton alleles = null;
        if (myAlignment.hasGenotype()) {
            alleles = new JRadioButton("Alleles");
            alleles.setActionCommand(AlignmentTableModel.COLUMN_NAME_TYPE.alleles.toString());
            alleles.addActionListener(radioListener);
        }

        searchLabel = new JLabel();
        searchLabel.setPreferredSize(new Dimension(70, 25));
        searchField = new JTextField();
        searchField.setPreferredSize(new Dimension(225, 25));
        searchField.setText("(Enter physical position)");
        searchButton = new JButton("Search");

        searchButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                searchButton_actionPerformed(e);
            }
        });

        JPanel selectColumnHeadings = new JPanel(new FlowLayout(FlowLayout.LEFT));

        FlowLayout layout = new FlowLayout(FlowLayout.CENTER, 0, 0) {
            @Override
            public Dimension preferredLayoutSize(Container target) {

                int maxWidth = (target.getSize().width == 0) ? Integer.MAX_VALUE : target.getSize().width;

                int resultWidth = 0;
                int resultHeight = 0;

                int currentWidth = 0;
                int currentHeight = 0;
                Component[] components = target.getComponents();
                for (int i = 0; i < components.length; i++) {

                    if (components[i].isVisible()) {
                        Dimension currentSize = components[i].getPreferredSize();

                        if (currentWidth + currentSize.width > maxWidth) {
                            resultWidth = Math.max(resultWidth, currentWidth);
                            resultHeight += currentHeight;
                            currentWidth = 0;
                            currentHeight = 0;
                        }

                        currentWidth += currentSize.width;
                        currentHeight = Math.max(currentHeight, currentSize.height);
                    }

                }

                resultWidth = Math.max(resultWidth, currentWidth);
                resultHeight += currentHeight;

                return new Dimension(resultWidth, resultHeight);
            }
        };
        selectColumnHeadings.setLayout(layout);

        if (!multipleAlignments) {
            buttonGroup.add(physicalPosition);
        }
        buttonGroup.add(siteNumber);
        buttonGroup.add(locus);
        buttonGroup.add(siteName);
        if (myAlignment.hasGenotype()) {
            buttonGroup.add(alleles);
        }

        if (multipleAlignments) {
            buttonGroup.setSelected(siteNumber.getModel(), true);
        } else {
            buttonGroup.setSelected(physicalPosition.getModel(), true);
        }

        if (!multipleAlignments) {
            selectColumnHeadings.add(physicalPosition);
        }
        selectColumnHeadings.add(siteNumber);
        selectColumnHeadings.add(locus);
        selectColumnHeadings.add(siteName);
        if (myAlignment.hasGenotype()) {
            selectColumnHeadings.add(alleles);
        }

        JLabel blankSpace = new JLabel();
        blankSpace.setPreferredSize(new Dimension(25, 25));
        selectColumnHeadings.add(blankSpace);

        selectColumnHeadings.add(myHighlightingComboBox);

        JLabel blankSpace2 = new JLabel();
        blankSpace2.setPreferredSize(new Dimension(25, 25));
        selectColumnHeadings.add(blankSpace2);

        selectColumnHeadings.add(searchLabel);
        selectColumnHeadings.add(searchField);
        selectColumnHeadings.add(searchButton);

        result.add(selectColumnHeadings);

        result.add(getSliderPane());

        return result;

    }

    public AlignmentTableCellRenderer.RENDERING_TYPE getCellRenderingType() {
        return (AlignmentTableCellRenderer.RENDERING_TYPE) myHighlightingComboBox.getSelectedItem();
    }

    private void searchButton_actionPerformed(ActionEvent e) {
        try {
            int searchValue = (int) Double.parseDouble(searchField.getText().trim());
            if (myTableModel.getColumnNameType().equals(AlignmentTableModel.COLUMN_NAME_TYPE.physicalPosition)) {
                if (searchValue > endPos) {
                    JOptionPane.showMessageDialog(this.getParent(), "Physical position must be between " + startPos + " and " + endPos + ".");
                } else if (searchValue < startPos) {
                    JOptionPane.showMessageDialog(this.getParent(), "Physical position must be between " + startPos + " and " + endPos + ".");
                } else {
                    mySlider.setValue(searchValue);
                }
            } else if (myTableModel.getColumnNameType().equals(AlignmentTableModel.COLUMN_NAME_TYPE.siteNumber)) {
                if (searchValue > end) {
                    JOptionPane.showMessageDialog(this.getParent(), "Site number must be between " + start + " and " + end + ".");
                } else if (searchValue < start) {
                    JOptionPane.showMessageDialog(this.getParent(), "Site number must be between " + start + " and " + end + ".");
                } else {
                    mySlider.setValue(searchValue);
                }
            }
        } catch (Exception ex) {
            String positionType = "";
            if (myTableModel.getColumnNameType().equals(AlignmentTableModel.COLUMN_NAME_TYPE.physicalPosition)) {
                positionType = "physical position";
            } else if (myTableModel.getColumnNameType().equals(AlignmentTableModel.COLUMN_NAME_TYPE.siteNumber)) {
                positionType = "site number";
            }
            JOptionPane.showMessageDialog(this.getParent(), "Invalid " + positionType + ": " + searchField.getText().trim());
        }
    }

    private void hideSearchFunction() {
        searchField.setVisible(false);
        searchButton.setVisible(false);
    }

    private void showSearchFunction() {
        searchField.setVisible(true);
        searchButton.setVisible(true);
    }

    private void changeTextPhysicalPosition() {
        searchField.setText("(Enter physical position)");
    }

    private void changeTextSiteNumber() {
        searchField.setText("(Enter site number)");
    }

    private JPanel getSliderPane() {

        mySliderPane = new JPanel(new BorderLayout());

        URL imageURL = SeqViewerPanel.class.getResource("left.gif");
        ImageIcon imageIcon = null;
        if (imageURL != null) {
            imageIcon = new ImageIcon(imageURL);
        }

        myLeftButton = new JButton(imageIcon);
        myLeftButton.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                if (myTableModel.isPhysicalPosition()) {
                    int newSite = myTableModel.getHorizontalCenter() - myTableModel.getHorizontalPageSize() * 3 / 4;
                    newSite = Math.max(0, newSite);
                    mySlider.setValue(myAlignment.chromosomalPosition(newSite));
                } else {
                    int newValue = mySlider.getValue() - myTableModel.getHorizontalPageSize() * 3 / 4;
                    newValue = Math.max(mySlider.getMinimum(), newValue);
                    mySlider.setValue(newValue);
                }
            }
        });
        mySliderPane.add(myLeftButton, BorderLayout.WEST);

        mySliderPane.add(mySlider, BorderLayout.CENTER);

        imageURL = SeqViewerPanel.class.getResource("right.gif");
        imageIcon = null;
        if (imageURL != null) {
            imageIcon = new ImageIcon(imageURL);
        }

        myRightButton = new JButton(imageIcon);
        myRightButton.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                if (myTableModel.isPhysicalPosition()) {
                    int newSite = myTableModel.getHorizontalCenter() + myTableModel.getHorizontalPageSize() * 3 / 4;
                    newSite = Math.min(myAlignment.numberOfSites() - 1, newSite);
                    mySlider.setValue(myAlignment.chromosomalPosition(newSite));
                } else {
                    int newValue = mySlider.getValue() + myTableModel.getHorizontalPageSize() * 3 / 4;
                    newValue = Math.min(mySlider.getMaximum(), newValue);
                    mySlider.setValue(newValue);
                }
            }
        });
        mySliderPane.add(myRightButton, BorderLayout.EAST);

        return mySliderPane;

    }

    public int getSliderPosition() {
        return mySlider.getValue();
    }

    public int getSliderPositionAsSite() {
        return myTableModel.getHorizontalCenter();
    }

    public GenotypeTable getAlignment() {
        return myAlignment;
    }

    private void updateSliderPhysicalPositions() {

        int min = myAlignment.chromosomalPosition(0);
        int max = myAlignment.chromosomalPosition(myAlignment.numberOfSites() - 1);
        int tableSize = max - min + 1;
        int center = myTableModel.getHorizontalCenter();
        mySlider.setMinimum(min);
        mySlider.setMaximum(max);
        mySlider.setPaintTicks(true);
        mySlider.setPaintLabels(true);
        int scrollWidth = myScrollPane.getWidth();
        int numTicks = 0;
        if (scrollWidth <= 0) {
            numTicks = 10;
        } else {
            numTicks = scrollWidth / SLIDER_TEXT_WIDTH;
        }
        int spacing = tableSize / numTicks;
        if (spacing < 1) {
            spacing = 1;
        }
        mySlider.setLabelTable(mySlider.createStandardLabels(spacing));
        mySlider.setMajorTickSpacing(spacing);
        mySlider.setValue(myAlignment.chromosomalPosition(center));
        mySlider.validate();
        mySliderPane.validate();
        mySliderPane.repaint();
        repaint();

    }

    private void updateSliderSiteNumbers() {

        int tableSize = myTableModel.getRealColumnCount();
        int center = myTableModel.getHorizontalCenter();
        mySlider.setMinimum(0);
        mySlider.setMaximum(tableSize - 1);
        mySlider.setPaintTicks(true);
        mySlider.setPaintLabels(true);
        int scrollWidth = myScrollPane.getWidth();
        int numTicks = 0;
        if (scrollWidth <= 0) {
            numTicks = 10;
        } else {
            numTicks = scrollWidth / SLIDER_TEXT_WIDTH;
        }
        int spacing = tableSize / numTicks;
        if (spacing < 1) {
            spacing = 1;
        }
        mySlider.setLabelTable(mySlider.createStandardLabels(spacing));
        mySlider.setMajorTickSpacing(spacing);
        mySlider.setValue(center);
        mySlider.validate();
        mySliderPane.validate();
        mySliderPane.repaint();
        repaint();

    }

    public void componentResized(ComponentEvent e) {
        myTableModel.setHorizontalPageSize(calculateNumColumns());

        if (myTableModel.isPhysicalPosition()) {
            updateSliderPhysicalPositions();
        } else if (!myTableModel.isPhysicalPosition()) {
            updateSliderSiteNumbers();
        }
    }

    private int calculateNumColumns() {
        if (myScrollPane == null) {
            return 10;
        }
        Dimension size = myScrollPane.getSize();
        int width = (int) size.getWidth() - ROW_HEADER_WIDTH - SCROLL_BAR_WIDTH;
        return width / myTableModel.getColumnWidth();
    }

    public void componentMoved(ComponentEvent e) {
        // do nothing
    }

    public void componentShown(ComponentEvent e) {
        // do nothing
    }

    public void componentHidden(ComponentEvent e) {
        // do nothing
    }

    public void tableChanged(TableModelEvent e) {
        for (int c = 0; c < myTable.getColumnCount(); c++) {
            TableColumn col = myTable.getColumnModel().getColumn(c);
            col.setHeaderRenderer(new TableCellRenderer() {
                @Override
                public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
                    JLabel label = new JLabel(value.toString());
                    label.setBorder(BorderFactory.createEtchedBorder());
                    label.setUI(VerticalLabelUI.getInstance());
                    return label;
                }
            });
        }

        if (myTableModel.isPhysicalPosition()) {
            updateSliderPhysicalPositions();
        } else if (!myTableModel.isPhysicalPosition()) {
            updateSliderSiteNumbers();
        }
    }

    class RadioListener implements ActionListener {

        public void actionPerformed(ActionEvent e) {
            if (e.getActionCommand().equals(AlignmentTableModel.COLUMN_NAME_TYPE.physicalPosition.toString())) {
                myTableModel.setColumnNameType(AlignmentTableModel.COLUMN_NAME_TYPE.physicalPosition);
                updateSliderPhysicalPositions();
                showSearchFunction();
                changeTextPhysicalPosition();
            } else if (e.getActionCommand().equals(AlignmentTableModel.COLUMN_NAME_TYPE.siteNumber.toString())) {
                myTableModel.setColumnNameType(AlignmentTableModel.COLUMN_NAME_TYPE.siteNumber);
                updateSliderSiteNumbers();
                showSearchFunction();
                changeTextSiteNumber();
            } else if (e.getActionCommand().equals(AlignmentTableModel.COLUMN_NAME_TYPE.locus.toString())) {
                myTableModel.setColumnNameType(AlignmentTableModel.COLUMN_NAME_TYPE.locus);
                hideSearchFunction();
            } else if (e.getActionCommand().equals(AlignmentTableModel.COLUMN_NAME_TYPE.alleles.toString())) {
                myTableModel.setColumnNameType(AlignmentTableModel.COLUMN_NAME_TYPE.alleles);
                hideSearchFunction();
            } else if (e.getActionCommand().equals(AlignmentTableModel.COLUMN_NAME_TYPE.siteName.toString())) {
                myTableModel.setColumnNameType(AlignmentTableModel.COLUMN_NAME_TYPE.siteName);
                hideSearchFunction();
            }
        }
    }

    public class MyTableUI extends BasicTableUI {

        private boolean justClicked = false;

        protected MouseInputListener createMouseInputListener() {
            return (new MyMouseInputListener());
        }

        public class MyMouseInputListener extends BasicTableUI.MouseInputHandler {

            public void mousePressed(MouseEvent e) {

                JTable table = (JTable) e.getSource();

                if (!table.isRowSelected(table.rowAtPoint(e.getPoint()))) {

                    int mod = e.getModifiers();

                    if ((mod & MouseEvent.BUTTON3_MASK) == MouseEvent.BUTTON3_MASK) {
                        mod = mod ^ MouseEvent.BUTTON3_MASK;
                        mod = mod | MouseEvent.BUTTON1_MASK;
                        e = new MouseEvent((Component) e.getSource(), e.getID(), e.getWhen(), mod, e.getX(), e.getY(), e.getClickCount(), true);
                    }

                    super.mousePressed(e);

                    justClicked = true;

                }

            }

            public void mouseClicked(MouseEvent e) {

                if (!justClicked) {
                    super.mousePressed(e);
                } else {
                    justClicked = false;
                }

            }

            public void mouseDragged(MouseEvent e) {
            }
        }
    }
}
