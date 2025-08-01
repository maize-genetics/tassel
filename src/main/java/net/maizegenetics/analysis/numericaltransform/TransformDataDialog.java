package net.maizegenetics.analysis.numericaltransform;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextField;

import net.maizegenetics.analysis.numericaltransform.TransformDataPlugin.BASE;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.util.OpenBitSet;

public class TransformDataDialog extends JDialog implements ActionListener {
	private final List<NumericAttribute> dataAttributes;
	private final List<CategoricalAttribute> factorAttributes;
	private JList<NumericAttribute> lstTraits;
	private JList<CategoricalAttribute> lstFactors;
	private JCheckBox chkPower;
	private JCheckBox chkLog;
	private JRadioButton radNatural;
	private JRadioButton radBase2;
	private JRadioButton radBase10;
	private JTextField txtPower;
	private JCheckBox chkStandard;
	private JButton btnOk;
	private JButton btnCancel;
	private boolean wasCancelled = true;
	private final boolean hasFactors;

	public TransformDataDialog(Frame parentFrame, List<NumericAttribute> dataAttributes, List<CategoricalAttribute> factorAttributes) {
		super(parentFrame);
		this.dataAttributes = dataAttributes;
		this.factorAttributes = factorAttributes;
		if (factorAttributes == null || factorAttributes.size() == 0) hasFactors = false;
		else hasFactors = true;
		init();
	}

	//for testing dialog
	public static void main(String[] args) {
		List<NumericAttribute> dataAttributes = new ArrayList<>();
		dataAttributes.add(new NumericAttribute("trait1", new float[]{1f,2f,3f}, new OpenBitSet(3)));
		dataAttributes.add(new NumericAttribute("trait2", new float[]{1f,2f,3f}, new OpenBitSet(3)));
		dataAttributes.add(new NumericAttribute("trait3", new float[]{1f,2f,3f}, new OpenBitSet(3)));
		dataAttributes.add(new NumericAttribute("trait4", new float[]{1f,2f,3f}, new OpenBitSet(3)));
		dataAttributes.add(new NumericAttribute("trait5", new float[]{1f,2f,3f}, new OpenBitSet(3)));
		dataAttributes.add(new NumericAttribute("trait6", new float[]{1f,2f,3f}, new OpenBitSet(3)));
		dataAttributes.add(new NumericAttribute("trait7", new float[]{1f,2f,3f}, new OpenBitSet(3)));

		List<CategoricalAttribute> factorAttributes = new ArrayList<>();
		factorAttributes.add(new CategoricalAttribute("factor1", new String[]{"one","one","two"}));

		TransformDataDialog tdd = new TransformDataDialog(null, dataAttributes, null);
		tdd.setVisible(true);
		
		System.out.printf("power transformation: %b; exponent = %1.1f\n", tdd.powerTransformation(), tdd.exponent());
		
		System.out.printf("log transformation: %b; base = %s\n", tdd.logTransformation(), tdd.base().name());
		System.out.println("Traits to transform:");
		for (NumericAttribute na : tdd.traitsToTransform()) System.out.println(na.name());
		System.out.println("Factors to standardize:");
		for (CategoricalAttribute ca : tdd.factorsForStandardizing()) System.out.println(ca.name());
		if (tdd.standardize()) System.out.println("Standardize data.");
		else System.out.println("Do not standardize data.");
		if(tdd.wasCancelled) System.out.println("was cancelled.");
		else System.out.println("was not cancelled.");
		System.exit(0);
	}

	private void init() {
		setTitle("Transform Data");
		setModalityType(ModalityType.APPLICATION_MODAL);

		lstTraits = new JList<>(dataAttributes.toArray(new NumericAttribute[0]));
		if (hasFactors) lstFactors = new JList<>(factorAttributes.toArray(new CategoricalAttribute[0]));

		chkPower = new JCheckBox("Power", false);
		chkPower.setActionCommand("power");
		chkPower.addActionListener(this);
		txtPower = new JTextField("Enter exponent");
		txtPower.setEnabled(false);

		chkLog = new JCheckBox("Log", false);
		chkLog.setActionCommand("log");
		chkLog.addActionListener(this);
		radNatural = new JRadioButton("Natural", true);
		radNatural.setEnabled(false);
		radBase2 = new JRadioButton("Base 2", false);
		radBase2.setEnabled(false);
		radBase10 = new JRadioButton("Base 10", false);
		radBase10.setEnabled(false);
		ButtonGroup logBaseGroup = new ButtonGroup();
		logBaseGroup.add(radNatural);
		logBaseGroup.add(radBase2);
		logBaseGroup.add(radBase10);

		chkStandard = new JCheckBox("Standardize");
		btnOk = new JButton("OK");
		btnOk.setActionCommand("ok");
		btnOk.addActionListener(this);
		btnCancel = new JButton("Cancel");
		btnCancel.setActionCommand("cancel");
		btnCancel.addActionListener(this);

		JPanel mainPanel = new JPanel(new GridBagLayout());
		
		GridBagConstraints gbc = new GridBagConstraints();
		gbc.gridx = 0;
		gbc.gridy = 0;
		if (hasFactors) gbc.gridwidth = 1;
		else gbc.gridwidth = 2;
		gbc.anchor = GridBagConstraints.CENTER;
		if (!hasFactors) gbc.insets = new Insets(15, 60, 2, 60);
		else gbc.insets = new Insets(15, 30, 2, 30); //top, left, bottom, right
		gbc.weightx = 1;
		gbc.weighty = 1;
		JLabel lblTraits = new JLabel("Select Traits to Transform:");
		Font labelFont = lblTraits.getFont().deriveFont(Font.BOLD, 14f);
		lblTraits.setFont(labelFont);
		mainPanel.add(lblTraits, gbc);

		if (hasFactors) {
			gbc.gridx++;
			JLabel lblFactors = new JLabel("Standardize within Selected Factor:");
			lblFactors.setFont(labelFont);
			mainPanel.add(lblFactors, gbc);
		}

		gbc.gridx = 0;
		gbc.gridy++;
		JScrollPane traitScroller = new JScrollPane(lstTraits);
		traitScroller.setPreferredSize(new Dimension(150,100));
		gbc.insets = new Insets(5, 30, 20, 30);
		mainPanel.add(traitScroller, gbc);

		if (hasFactors) {
			gbc.gridx++;
			JScrollPane factorScroller = new JScrollPane(lstFactors);
			factorScroller.setPreferredSize(new Dimension(150,100));
			gbc.insets = new Insets(2, 30, 20, 30);
			mainPanel.add(factorScroller, gbc);
		}

		gbc.gridx = 0;
		gbc.gridy++;
		gbc.gridwidth = 2;
		gbc.insets = new Insets(5, 2, 5, 2); //top, left, bottom, right
		JLabel lblMethod = new JLabel("Select Transformation Method:");
		lblMethod.setFont(labelFont);
		
		mainPanel.add(lblMethod, gbc);

		gbc.gridy++;
		gbc.gridwidth = 1;
		gbc.insets = new Insets(2, 2, 2, 2); //top, left, bottom, right
		gbc.anchor = GridBagConstraints.EAST;
		mainPanel.add(chkPower,gbc);

		gbc.gridx++;
		gbc.anchor = GridBagConstraints.WEST;
		mainPanel.add(txtPower,gbc);

		gbc.gridy++;
		gbc.gridx = 0;
		gbc.gridwidth = 1;
		gbc.anchor = GridBagConstraints.EAST;
		gbc.insets = new Insets(20, 2, 2, 2); //top, left, bottom, right
		mainPanel.add(chkLog,gbc);

		gbc.gridx++;
		gbc.anchor = GridBagConstraints.WEST;
		mainPanel.add(radNatural,gbc);

		gbc.gridy++;
		gbc.insets = new Insets(2, 2, 2, 2); //top, left, bottom, right
		mainPanel.add(radBase2,gbc);

		gbc.gridy++;
		mainPanel.add(radBase10,gbc);

		gbc.gridy++;
		gbc.gridx = 0;
		gbc.gridwidth = 2;
		gbc.anchor = GridBagConstraints.CENTER;
		gbc.insets = new Insets(20, 2, 2, 2); //top, left, bottom, right
		mainPanel.add(chkStandard,gbc);

		gbc.gridy++;
		gbc.gridwidth = 1;
		gbc.anchor = GridBagConstraints.EAST;
		gbc.insets = new Insets(20, 2, 30, 2); //top, left, bottom, right
		mainPanel.add(btnOk,gbc);

		gbc.gridx++;
		gbc.anchor = GridBagConstraints.WEST;
		mainPanel.add(btnCancel,gbc);

		getContentPane().add(mainPanel);
		pack();
		setLocationRelativeTo(getParent());
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("power")) {
			if (chkPower.isSelected()) txtPower.setEnabled(true);
			else txtPower.setEnabled(false);
		} else if (e.getActionCommand().equals("log")) {
			if (chkLog.isSelected()) {
				radNatural.setEnabled(true);
				radBase2.setEnabled(true);
				radBase10.setEnabled(true);
			} else {
				radNatural.setEnabled(false);
				radBase2.setEnabled(false);
				radBase10.setEnabled(false);
			}
		} else if (e.getActionCommand().equals("ok")) {
			wasCancelled = false;
			setVisible(false);
		} else if (e.getActionCommand().equals("cancel")) {
			wasCancelled = true;
			setVisible(false);
		}
		
	}
	
	//getters
	public List<NumericAttribute> traitsToTransform() {
		if (wasCancelled) return new ArrayList<>();
		return lstTraits.getSelectedValuesList();
	}
	
	public List<CategoricalAttribute> factorsForStandardizing() {
		if (hasFactors && !wasCancelled) return lstFactors.getSelectedValuesList();
		else return new ArrayList<CategoricalAttribute>();
	}
	
	public boolean powerTransformation() {
		return chkPower.isSelected();
	}
	
	public double exponent() {
		try {
			return Double.parseDouble(txtPower.getText());
		} catch (NumberFormatException nfe) {
			return Double.NaN;
		}
	}
	
	public boolean logTransformation() {
		return chkLog.isSelected();
	}
	
	public BASE base() {
		if (radBase2.isSelected()) return BASE.base_2;
		if (radBase10.isSelected()) return BASE.base_10;
		else return BASE.natural;
	}

	public boolean standardize() {
		return chkStandard.isSelected();
	}
	
	public boolean wasCancelled() {
		return wasCancelled;
	}
}
