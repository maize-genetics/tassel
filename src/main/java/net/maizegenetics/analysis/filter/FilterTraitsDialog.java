package net.maizegenetics.analysis.filter;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;

import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;

import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;

public class FilterTraitsDialog extends JDialog implements ActionListener, TableModelListener {

    private Phenotype myPhenotype;
    private TraitTableModel traitModel;
    private JScrollPane jsp;
    private JTable traitTable;
    private boolean clickedOK = false;
    public final static String CMD_EXCLUDE = "exclude";
    public final static String CMD_INCLUDE = "include";
    public final static String CMD_EXCLUDE_ALL = "excludeall";
    public final static String CMD_INCLUDE_ALL = "includeall";
    public final static String CMD_OK = "ok";
    public final static String CMD_CANCEL = "cancel";
    public final static String CMD_CHANGE_DATA = "change2data";
    public final static String CMD_CHANGE_COV = "change2cov";

    public FilterTraitsDialog(Frame parent, Phenotype aPhenotype) {
        super(parent);
        myPhenotype = aPhenotype;
        traitModel = new TraitTableModel(myPhenotype);
        init();
    }

    private void init() {
        setTitle("Filter Traits / Modify Trait Properties");
        setSize(new Dimension(600, 800));
        setModal(true);
        Container contentPane = getContentPane();
        contentPane.setLayout(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();

        //add the table
        traitModel = new TraitTableModel(myPhenotype);
        traitTable = new JTable(traitModel);
        jsp = new JScrollPane(traitTable);

        //create a combo box editor for type
        JComboBox<ATTRIBUTE_TYPE> comboType = new JComboBox();
        comboType.addItem(ATTRIBUTE_TYPE.data);
//        comboType.addItem(ATTRIBUTE_TYPE.factor);
        comboType.addItem(ATTRIBUTE_TYPE.covariate);
        int typeColNumber = traitModel.getTypeColumnNumber();

        TableColumn typeColumn = traitTable.getColumnModel().getColumn(typeColNumber);
        typeColumn.setCellEditor(new DefaultCellEditor(comboType));

        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridheight = 1;
        gbc.gridwidth = 2;
        contentPane.add(jsp, gbc);

        //add buttons
        JButton btnExclude = new JButton("Exclude Selected");
        btnExclude.setActionCommand(CMD_EXCLUDE);
        btnExclude.addActionListener(this);
        JButton btnInclude = new JButton("Include Selected");
        btnInclude.setActionCommand(CMD_INCLUDE);
        btnInclude.addActionListener(this);
        JButton btnExcludeAll = new JButton("Exclude All");
        btnExcludeAll.setActionCommand(CMD_EXCLUDE_ALL);
        btnExcludeAll.addActionListener(this);
        JButton btnIncludeAll = new JButton("Include All");
        btnIncludeAll.setActionCommand(CMD_INCLUDE_ALL);
        btnIncludeAll.addActionListener(this);
        
        JButton btnChangeData = new JButton("Change Selected Type to Data");
        btnChangeData.setActionCommand(CMD_CHANGE_DATA);
        btnChangeData.addActionListener(this);
        JButton btnChangeCov = new JButton("Change Selected Type to Covariate");
        btnChangeCov.setActionCommand(CMD_CHANGE_COV);
        btnChangeCov.addActionListener(this);
        
        JButton btnOK = new JButton("OK");
        btnOK.setActionCommand(CMD_OK);
        btnOK.addActionListener(this);
        JButton btnCancel = new JButton("Cancel");
        btnCancel.setActionCommand(CMD_CANCEL);
        btnCancel.addActionListener(this);

        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.gridwidth = 1;
        gbc.anchor = GridBagConstraints.EAST;
        gbc.weightx = 1;
        gbc.insets = new Insets(5, 5, 5, 5); //top, left, bottom, right
        contentPane.add(btnExclude, gbc);

        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.insets = new Insets(5, 5, 5, 5); //top, left, bottom, right
        gbc.anchor = GridBagConstraints.WEST;
        contentPane.add(btnInclude, gbc);

        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.insets = new Insets(5, 5, 5, 5); //top, left, bottom, right
        gbc.anchor = GridBagConstraints.EAST;
        contentPane.add(btnExcludeAll, gbc);

        gbc.gridx = 1;
        gbc.gridy = 2;
        gbc.insets = new Insets(5, 5, 5, 5); //top, left, bottom, right
        gbc.anchor = GridBagConstraints.WEST;
        contentPane.add(btnIncludeAll, gbc);

        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.gridwidth = 2;
        gbc.insets = new Insets(15, 5, 5, 5); //top, left, bottom, right
        gbc.anchor = GridBagConstraints.CENTER;
        contentPane.add(btnChangeData, gbc);
        gbc.gridx = 0;
        gbc.gridy = 4;
        gbc.insets = new Insets(5, 5, 5, 5); //top, left, bottom, right
        gbc.anchor = GridBagConstraints.CENTER;
        contentPane.add(btnChangeCov, gbc);
        
        gbc.gridx = 0;
        gbc.gridy = 5;
        gbc.gridwidth = 1;
        gbc.insets = new Insets(5, 5, 5, 5); //top, left, bottom, right
        gbc.anchor = GridBagConstraints.EAST;
        contentPane.add(btnOK, gbc);

        gbc.gridx = 1;
        gbc.gridy = 5;
        gbc.insets = new Insets(5, 5, 5, 5); //top, left, bottom, right
        gbc.anchor = GridBagConstraints.WEST;
        contentPane.add(btnCancel, gbc);
    }

    /**
     * @return an index of the traits to be included in the filter phenotype
     */
    public int[] getIncludedTraits() {
        if (!clickedOK) {
            return null;
        }
        Boolean[] isIncluded = traitModel.include;
        int nOrig = isIncluded.length;
        int[] includedTrait = new int[nOrig];
        
        int includedCount = 0;
        for (int i = 0; i < nOrig; i++) {
            if (isIncluded[i]) {
                includedTrait[includedCount++] = i;
            }
        }
        return Arrays.copyOf(includedTrait, includedCount);
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        if (e.getActionCommand().equals(CMD_EXCLUDE)) {
            traitModel.excludeSome(traitTable.getSelectedRows());
        } else if (e.getActionCommand().equals(CMD_INCLUDE)) {
            traitModel.includeSome(traitTable.getSelectedRows());
        } else if (e.getActionCommand().equals(CMD_EXCLUDE_ALL)) {
            traitModel.excludeAll();
        } else if (e.getActionCommand().equals(CMD_INCLUDE_ALL)) {
            traitModel.includeAll();
        } else if (e.getActionCommand().equals(CMD_OK)) {
            setVisible(false);
            clickedOK = true;
        } else if (e.getActionCommand().equals(CMD_CANCEL)) {
            setVisible(false);
            clickedOK = false;
        } else if (e.getActionCommand().equals(CMD_CHANGE_DATA)) {
        	changeSelectedType(ATTRIBUTE_TYPE.data);
        } else if (e.getActionCommand().equals(CMD_CHANGE_COV)) {
        	changeSelectedType(ATTRIBUTE_TYPE.covariate);
        }
        
    }

    private void changeSelectedType(ATTRIBUTE_TYPE toNewType) {
    	for (int i : traitTable.getSelectedRows()) {
    		traitModel.setValueAt(toNewType, i, traitModel.getTypeColumnNumber());
    	}
    	traitModel.fireTableDataChanged();
    }
    
    @Override
    public void tableChanged(TableModelEvent e) {
        System.out.println("table model changed");
    }

    public boolean getClickedOK() {
        return clickedOK;
    }
    
    public  Map<PhenotypeAttribute, ATTRIBUTE_TYPE> getTypeChangeMap() {
    	Map<PhenotypeAttribute, ATTRIBUTE_TYPE> typeChangeMap = new HashMap<>();
    	for (int trait : getIncludedTraits()) {
    		if (traitModel.types[trait] != myPhenotype.attributeType(trait)) {
    			typeChangeMap.put(myPhenotype.attribute(trait), traitModel.types[trait]);
    		}
    	}
    	return typeChangeMap;
    }
}

class TraitTableModel extends AbstractTableModel {

//    ArrayList<Trait> traitList;
    List<PhenotypeAttribute> attributeList;
    Phenotype myPhenotype;
    Boolean[] include;
    ATTRIBUTE_TYPE[] types;
    String[] colName;
    int numberOfTraits;
    int numberOfColumns;
    int typeColumnNumber;

    TraitTableModel(Phenotype aPhenotype) {
        super();
        myPhenotype = aPhenotype;
        attributeList = myPhenotype.attributeListCopy();
        numberOfTraits = attributeList.size();
        numberOfColumns = 3;
        setColumnNames();
        include = new Boolean[numberOfTraits];
        types = new ATTRIBUTE_TYPE[numberOfTraits];
        for (int i = 0; i < numberOfTraits; i++) {
            include[i] = true;
            types[i] = myPhenotype.attributeType(i);
        }
    }

    private void setColumnNames() {
        colName = new String[numberOfColumns];
        int col = 0;
        colName[col++] = "Trait";
        colName[col++] = "Type";
        colName[col++] = "Include";
        typeColumnNumber = 1;
    }

    @Override
    public int getColumnCount() {
        return colName.length;
    }

    @Override
    public int getRowCount() {
        return include.length;
    }

    @Override
    public Object getValueAt(int row, int col) {
    	if (col == 0) return myPhenotype.attribute(row).name();
    	if (col == 1) return types[row];
    	if (col == 2) return include[row];
    	else return "";
    }

    @Override
    public Class<?> getColumnClass(int columnIndex) {
        if (columnIndex < 2) {
            return String.class;
        } else {
            return Boolean.class;
        }
    }

    @Override
    public String getColumnName(int column) {
        return colName[column];
    }

    @Override
    public boolean isCellEditable(int rowIndex, int columnIndex) {
    	if (myPhenotype.attributeType(rowIndex) == ATTRIBUTE_TYPE.taxa) return false;
    	else if (myPhenotype.attributeType(rowIndex) == ATTRIBUTE_TYPE.factor) return false;
    	else if (columnIndex > 0)  return true;
        return false;
    }

    @Override
    public void setValueAt(Object value, int rowIndex, int columnIndex) {
        if (columnIndex == 2) {
            include[rowIndex] = (Boolean) value;
        }
        if (columnIndex == 1 && value instanceof ATTRIBUTE_TYPE) {
            types[rowIndex] = (ATTRIBUTE_TYPE) value;
        }
    }

    public void excludeAll() {
        for (int i = 0; i < numberOfTraits; i++) {
            if (types[i] != ATTRIBUTE_TYPE.taxa) include[i] = Boolean.FALSE;
        }
        fireTableDataChanged();
    }

    public void includeAll() {
        for (int i = 0; i < numberOfTraits; i++) {
            include[i] = Boolean.TRUE;
        }
        fireTableDataChanged();
    }

    public void excludeSome(int[] index) {
        for (int i : index) {
        	if (types[i] != ATTRIBUTE_TYPE.taxa) include[i] = Boolean.FALSE;
        }
        fireTableDataChanged();
    }

    public void includeSome(int[] index) {
        for (int i : index) {
            include[i] = Boolean.TRUE;
        }
        fireTableDataChanged();
    }

    public int getTypeColumnNumber() {
        return typeColumnNumber;
    }
}
