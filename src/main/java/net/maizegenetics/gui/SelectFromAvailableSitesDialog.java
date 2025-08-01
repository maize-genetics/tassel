package net.maizegenetics.gui;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import net.maizegenetics.util.BitSet;

public class SelectFromAvailableSitesDialog extends SelectFromAvailableDialog {
	
	private Frame parentFrame;
	private List<int[]> myListOfSelectedIndices = new ArrayList<>();
	private List<int[]> myListOfCovariateIndices = new ArrayList<>();
	private List<int[]> myListOfFactorIndices = new ArrayList<>();
	
	public SelectFromAvailableSitesDialog(Frame frame, String title, AbstractAvailableListModel availableModel) {
		super(frame, title, availableModel);
		parentFrame = frame;
		myIsCanceled = true;
	}
	
	@Override
	protected JPanel getBottomPanel() {
        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);

        JPanel firstRow = new JPanel(new FlowLayout());

        firstRow.add(getSitesToCovariateButton());
        firstRow.add(Box.createRigidArea(new Dimension(15, 1)));
        firstRow.add(getSitesToFactorsButton());
        firstRow.add(Box.createRigidArea(new Dimension(15, 1)));
        firstRow.add(getRemoveButton());
        
        JPanel secondRow = new JPanel(new FlowLayout());

        secondRow.add(getSelectedButton());
        secondRow.add(Box.createRigidArea(new Dimension(15, 1)));
        secondRow.add(getUnselectedButton());
        secondRow.add(Box.createRigidArea(new Dimension(15, 1)));
        secondRow.add(getOkButton());
        secondRow.add(Box.createRigidArea(new Dimension(15, 1)));
        secondRow.add(getCancelButton());

        result.add(firstRow);
        result.add(Box.createRigidArea(new Dimension(1, 10)));
        result.add(secondRow);
        result.add(Box.createRigidArea(new Dimension(1, 10)));
        return result;
	}
	
    protected JButton getSitesToCovariateButton() {
        JButton okButton = new JButton("Selected Sites => Covariates");
        okButton.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
            	String msg = "Selected sites will be converted to covariates.";
            	int response = JOptionPane.showConfirmDialog(parentFrame, msg, "Sites to Covariates", JOptionPane.OK_CANCEL_OPTION, JOptionPane.INFORMATION_MESSAGE);
            	if (response == JOptionPane.OK_OPTION) myListOfCovariateIndices.add(mySelectedListModel.getBitSet().getIndicesOfSetBits());
            }
        });
        return okButton;
    }

    protected JButton getSitesToFactorsButton() {
        JButton okButton = new JButton("Selected Sites => Factors");
        okButton.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
            	String msg = "Selected sites will be converted to factors.";
            	int response = JOptionPane.showConfirmDialog(parentFrame, msg, "Sites to Factors", JOptionPane.OK_CANCEL_OPTION, JOptionPane.INFORMATION_MESSAGE);
            	if (response == JOptionPane.OK_OPTION) myListOfFactorIndices.add(mySelectedListModel.getBitSet().getIndicesOfSetBits());
            }
        });
        return okButton;
    }

    protected JButton getOkButton() {
        JButton okButton = new JButton("OK");
        okButton.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
            	setVisible(false);
            	myIsCanceled = false;
            }
        });
        return okButton;
    }

	@Override
	protected JButton getSelectedButton() {
        JButton okButton = new JButton("Capture Selected");
        okButton.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
            	String msg = "A new genotype data set will be created using only the selected sites.";
            	int response = JOptionPane.showConfirmDialog(parentFrame, msg, "Sites to Factors", JOptionPane.OK_CANCEL_OPTION, JOptionPane.INFORMATION_MESSAGE);
            	if (response == JOptionPane.OK_OPTION) myListOfSelectedIndices.add(mySelectedListModel.getBitSet().getIndicesOfSetBits());
            }
        });
        return okButton;
	}

	@Override
	protected JButton getUnselectedButton() {
        JButton okButton = new JButton("Capture Unselected");
        okButton.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
            	String msg = "A new genotype data set with the selected sites excluded.";
            	int response = JOptionPane.showConfirmDialog(parentFrame, msg, "Sites to Factors", JOptionPane.OK_CANCEL_OPTION, JOptionPane.INFORMATION_MESSAGE);
            	if (response == JOptionPane.OK_OPTION) {
                    BitSet temp = mySelectedListModel.getBitSet();
                    temp.flip(0, myAvailableListModel.getRealSize());
            		myListOfSelectedIndices.add(temp.getIndicesOfSetBits());
            	}
            }
        });
        return okButton;
	}
    
    public List<int[]> listOfSelectedIndices() {
    	return myListOfSelectedIndices;
    }
    
    public List<int[]> listOfCovariateIndices() {
    	return myListOfCovariateIndices;
    }

    public List<int[]> listOfFactorIndices() {
    	return myListOfFactorIndices;
    }

}
