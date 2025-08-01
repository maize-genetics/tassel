package net.maizegenetics.analysis.association;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.gui.ReportDestinationDialog;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.distance.DistanceMatrix;

public class MLMPlugin extends AbstractPlugin {
    
    public boolean isUseP3D() {
		return useP3D;
	}

	public void setUseP3D(boolean useP3D) {
		this.useP3D = useP3D;
	}

	public boolean isUseGenotype() {
		return useGenotype;
	}

	public void setUseGenotype(boolean useGenotype) {
		this.useGenotype = useGenotype;
	}

	public boolean isUseRefProb() {
		return useRefProb;
	}

	public void setUseRefProb(boolean useRefProb) {
		this.useRefProb = useRefProb;
	}

	private static final Logger myLogger = LogManager.getLogger(MLMPlugin.class);
    protected DistanceMatrix kinshipMatrix;
    protected boolean analyzeByColumn;
    protected boolean useP3D = true;
    protected CompressionType compressionType = CompressionType.None;
    protected double compression = 1;
    private boolean writeOutputToFile = false;
    private String outputName = null;
    private boolean filterOutput = false;
    private double maxp = 1;
    private boolean useGenotype = true;
    private boolean useRefProb = false;
    private boolean useAlleleProb = false;

    public enum CompressionType {
        Optimum, Custom, None
    };

    public MLMPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);

    }

    @Override
    public DataSet processData(DataSet input) {
        boolean hasGenotypePhenotype = true;
        try {

            java.util.List<Datum> alignInList = input.getDataOfType(GenotypePhenotype.class);
            if (alignInList.size() == 0) {
                hasGenotypePhenotype = false;
                alignInList = input.getDataOfType(Phenotype.class);
            }
            java.util.List<Datum> kinshipList = input.getDataOfType(DistanceMatrix.class);

            if (alignInList.size() != 1) {
                String message = "Invalid selection. Please select one dataset with marker and trait data.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    myLogger.error("performFunction: " + message);
                }
                return null;
            }
            
            if (kinshipList.size() != 1) {
                String message = "Please select exactly one kinship matrix.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    myLogger.error("performFunction: " + message);
                }
                return null;
            }

            kinshipMatrix = (DistanceMatrix) kinshipList.get(0).getData();
            Iterator<Datum> itr = alignInList.iterator();

            if (isInteractive()) {
                MLMOptionDialog theOptions;
                if (hasGenotypePhenotype) {
                    GenotypePhenotype gp = (GenotypePhenotype) alignInList.get(0).getData();
                    theOptions = new MLMOptionDialog(getParentFrame(), hasDataTypes(gp));
                } else {
                    theOptions = new MLMOptionDialog(getParentFrame(), new boolean[]{false, false, false});
                }

                if (theOptions.runClicked) {
                    useP3D = theOptions.useP3D();
                    compressionType = theOptions.getCompressionType();
                    compression = theOptions.getCompressionLevel();
                    theOptions.dispose();

                    // give the user the option of sending the output to a file
                    ReportDestinationDialog rdd = new ReportDestinationDialog();
                    rdd.setLocationRelativeTo(getParentFrame());
                    rdd.setVisible(true);
                    if (!rdd.isOkayChecked()) {
                        return null;
                    }
                    writeOutputToFile = rdd.wasUseFileChecked();
                    if (writeOutputToFile) {
                        outputName = rdd.getOutputFileName();
                    }
                    filterOutput = rdd.wasRestrictOutputChecked();
                    if (filterOutput) {
                        maxp = rdd.getMaxP();
                    }

                } else {
                    theOptions.dispose();
                    return null;
                }

            } else {         //non-interactive stuff
            }

            List<Datum> myResults = new ArrayList<Datum>();
            while (itr.hasNext()) {

            	Datum current = itr.next();
            	CompressedMLMusingDoubleMatrix theAnalysis;
            	
            	if (hasGenotypePhenotype) {
                    GenotypeTable myGenotype = ((GenotypePhenotype) current.getData()).genotypeTable();
                    useGenotype = myGenotype.hasGenotype();
                    if (!useGenotype) useRefProb = myGenotype.hasReferenceProbablity();
            	} else {
            	    useGenotype = false;
            	    useRefProb = false;
            	}

            	if (useP3D) {
            		if (compressionType.equals(CompressionType.Optimum)) {
            			theAnalysis = new CompressedMLMusingDoubleMatrix(this, current, kinshipMatrix, true, true, Double.NaN);
            		} else if (compressionType.equals(CompressionType.Custom)) {
            			theAnalysis = new CompressedMLMusingDoubleMatrix(this, current, kinshipMatrix, true, true, compression);
            		} else {
            			theAnalysis = new CompressedMLMusingDoubleMatrix(this, current, kinshipMatrix, false, true, Double.NaN);
            		}
            		theAnalysis.useGenotypeCalls(useGenotype);
            		theAnalysis.useReferenceProbability(useRefProb);
            		theAnalysis.useAlleleProbabilities(useAlleleProb);

            	} else {
            		if (compressionType.equals(CompressionType.Optimum)) {
            			theAnalysis = new CompressedMLMusingDoubleMatrix(this, current, kinshipMatrix, true, false, Double.NaN);
            		} else if (compressionType.equals(CompressionType.Custom)) {
            			theAnalysis = new CompressedMLMusingDoubleMatrix(this, current, kinshipMatrix, true, false, compression);
            		} else {
            			theAnalysis = new CompressedMLMusingDoubleMatrix(this, current, kinshipMatrix, false, false, Double.NaN);
            		}

            		theAnalysis.useGenotypeCalls(useGenotype);
            		theAnalysis.useReferenceProbability(useRefProb);
            		theAnalysis.useAlleleProbabilities(useAlleleProb);
            		
            	}

            	myResults.addAll(theAnalysis.solve());
           }

            if (myResults.size() > 0 && !writeOutputToFile) {
            	fireDataSetReturned(new DataSet(myResults, this));
            	return new DataSet(myResults, this);
            }
            else return null;

        } finally {
            fireProgress(100);
        }

    }

    private boolean[] hasDataTypes(GenotypePhenotype gp) {
    	boolean[] hasTypes = new boolean[]{false, false, false};
    	if (gp.genotypeTable().hasGenotype()) hasTypes[0] = true;
    	if (gp.genotypeTable().hasReference()) hasTypes[1] = true;
    	if (gp.genotypeTable().hasAlleleProbabilities()) hasTypes[2] = true;
    	return hasTypes;
    }
    
    public ImageIcon getIcon() {
        URL imageURL = MLMPlugin.class.getResource("/net/maizegenetics/analysis/images/Mix.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    public String getButtonName() {
        return "MLM";
    }

    public String getToolTipText() {
        return "Association analysis using mixed model";
    }

    //a few stub functions to avoid producing errors in existing pipeline code
    public void setAnalyzeByColumn(boolean analyzeByColumn) {
        this.analyzeByColumn = analyzeByColumn;
    }

    public void setMaximumNumOfIteration(int max) {/*does nothing*/

    }

    public void setFinalIterMarker(boolean myFinalIterMarker) {/*does nothing*/

    }

    public void addFactors(int[] factors) {/*does nothing*/

    }

    public void setColumnTypes(String[] types) {/*does nothing*/

    }

    public void addFactors(String[] names) {/*does nothing*/

    }

    public void updateProgress(int progress) {
        if (progress < 0) {
            progress = 0;
        } else if (progress > 100) {
            progress = 100;
        }
        fireProgress(progress);
    }

    public void setVarCompEst(String value) {
        if (value.equalsIgnoreCase("P3D")) {
            useP3D = true;
        } else if (value.equalsIgnoreCase("EachMarker")) {
            useP3D = false;
        } else {
            throw new IllegalArgumentException("MLMPlugin: setVarCompEst: don't know how to handle value: " + value);
        }
    }

    public void setCompressionType(CompressionType type) {
        compressionType = type;
    }

    public boolean isWriteOutputToFile() {
        return writeOutputToFile;
    }

    public void setWriteOutputToFile(boolean writeOutputToFile) {
        this.writeOutputToFile = writeOutputToFile;
    }

    public String getOutputName() {
        return outputName;
    }

    public void setOutputName(String outputName) {
        this.outputName = outputName;
        this.writeOutputToFile = true;
    }

    public boolean isFilterOutput() {
        return filterOutput;
    }

    public void setFilterOutput(boolean filterOutput) {
        this.filterOutput = filterOutput;
    }

    public double getMaxp() {
        return maxp;
    }

    public void setMaxp(double maxp) {
        this.maxp = maxp;
        this.filterOutput = true;
    }

    public double getCustomCompression() {
        return compression;
    }

    public void setCustomCompression(double value) {
        compression = value;
    }
    
    public void useGenotypeCalls() {
        useGenotype = true;
        useRefProb = false;
        useAlleleProb = false;
    }
    
    public void useReferenceProbability() {
        useGenotype = false;
        useRefProb = true;
        useAlleleProb = false;
    }
    
    public void useAlleleProbabilities() {
        useGenotype = false;
        useRefProb = false;
        useAlleleProb = true;
    }

}

class MLMOptionDialog extends JDialog implements ActionListener {

    JRadioButton btnOptimum, btnCustom, btnNoCompression, btnEachMarker, btnP3D;
    ButtonGroup bgCompress, bgVariance, bgType;
    JTextField txtCustom;
    boolean runClicked = false;
    boolean useP3D = true;
    MLMPlugin.CompressionType compressionType = MLMPlugin.CompressionType.None;
    boolean useDiscrete = true;
    boolean useRefprob = false;
    boolean useAlleleprob = false;
    
    MLMOptionDialog(Frame parentFrame, boolean[] hasTypes) {
        super(parentFrame, true);
        final Frame pframe = parentFrame;
        setTitle("MLM Options");
        setSize(new Dimension(350, 300));
        setLocationRelativeTo(pframe);
        Container theContentPane = getContentPane();
        theContentPane.setLayout(new BorderLayout());
        JPanel compressionPanel = new JPanel(new GridBagLayout());
        compressionPanel.setBorder(BorderFactory.createTitledBorder("Compression Level"));

        //the method radio buttons
        btnOptimum = new JRadioButton("Optimum Level", false);
        btnOptimum.setActionCommand("Optimum");
        btnOptimum.addActionListener(this);
        btnCustom = new JRadioButton("Custom Level:", false);
        btnCustom.setActionCommand("Custom");
        btnCustom.addActionListener(this);
        btnNoCompression = new JRadioButton("No Compression", true);
        btnNoCompression.setActionCommand("None");
        btnNoCompression.addActionListener(this);
        btnEachMarker = new JRadioButton("Re-estimate after each marker", false);
        btnEachMarker.setActionCommand("Eachmarker");
        btnEachMarker.addActionListener(this);
        btnP3D = new JRadioButton("P3D (estimate once)", true);
        btnP3D.setActionCommand("P3D");
        btnP3D.addActionListener(this);

        bgCompress = new ButtonGroup();
        bgCompress.add(btnOptimum);
        bgCompress.add(btnCustom);
        bgCompress.add(btnNoCompression);
        bgVariance = new ButtonGroup();
        bgVariance.add(btnEachMarker);
        bgVariance.add(btnP3D);

        txtCustom = new JTextField(5);
        Insets inset1 = new Insets(5, 15, 5, 5);
        Insets inset2 = new Insets(5, 5, 5, 5);
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridwidth = 2;
        gbc.gridy = 0;
        gbc.weightx = 0;
        gbc.weighty = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = inset1; //top, left, bottom, right
        compressionPanel.add(btnOptimum, gbc);

        gbc.gridy++;
        gbc.gridwidth = 1;
        compressionPanel.add(btnCustom, gbc);
        gbc.gridx++;
        gbc.insets = inset2;
        compressionPanel.add(txtCustom, gbc);

        gbc.gridy++;
        gbc.gridx = 0;
        gbc.gridwidth = 2;
        gbc.insets = inset1;
        compressionPanel.add(btnNoCompression, gbc);

        theContentPane.add(compressionPanel, BorderLayout.NORTH);

        JPanel variancePanel = new JPanel(new GridBagLayout());
        variancePanel.setBorder(BorderFactory.createTitledBorder("Variance Component Estimation"));
        gbc.gridy = 0;
        variancePanel.add(btnP3D, gbc);
        gbc.gridy++;
        variancePanel.add(btnEachMarker, gbc);
        theContentPane.add(variancePanel, BorderLayout.CENTER);

        //panel for choosing a data type, if there is more than one
        int numberOfTypes = 0;
        for (boolean b : hasTypes) if (b) numberOfTypes++;
        if (numberOfTypes > 1) {
        	bgType = new ButtonGroup();
        	JPanel typePanel = new JPanel(new GridBagLayout());
        	typePanel.setBorder(BorderFactory.createTitledBorder("Choose genotype data"));
        	gbc.gridy = 0;
        	boolean initialValue = true;
        	
        	if (hasTypes[0]) {
        		JRadioButton btnDiscrete = new JRadioButton("Discrete type, e.g. SNP", initialValue);
                btnDiscrete.setActionCommand("discrete");
                btnDiscrete.addActionListener(this);
                initialValue = false;
                bgType.add(btnDiscrete);
                typePanel.add(btnDiscrete);
                gbc.gridy++;
        	}
        	if (hasTypes[1]) {
        		JRadioButton btnRef = new JRadioButton("Numeric genotype", initialValue);
        		btnRef.setActionCommand("reference");
        		btnRef.addActionListener(this);
                initialValue = false;
                bgType.add(btnRef);
                typePanel.add(btnRef);
                gbc.gridy++;
        	}
        	if (hasTypes[2]) {
        		JRadioButton btnAllele = new JRadioButton("Allele probabilities", initialValue);
        		btnAllele.setActionCommand("allele");
        		btnAllele.addActionListener(this);
                initialValue = false;
                bgType.add(btnAllele);
                typePanel.add(btnAllele);
                gbc.gridy++;
        	}
 
        	
        }
        
        
        //the help me button
        JButton btnHelpme = new JButton("Help Me Choose");
        final String msg = "For faster analysis, impute marker values before running MLM and use P3D.\n"
                + "With imputed marker values (no missing data), P3D will be very fast but compression will actually increase execution time somewhat.\n"
                + "However, because compression will improve the overall model fit, it should still be used.\n"
                + "If there is missing marker data, compression with P3D will probably be faster than P3D alone.\n"
                + "With small to moderate sized data sets any of the methods should give reasonable performance. \n"
                + "With large data sets consider imputing marker data then using EMMA with both compression and P3D";
        btnHelpme.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent arg0) {
                JOptionPane.showMessageDialog(pframe, msg, "EMMA methods", JOptionPane.INFORMATION_MESSAGE);
            }
        });


        //the run and cancel buttons
        JButton btnRun = new JButton("Run");
        JButton btnCancel = new JButton("Cancel");
        btnRun.setActionCommand("run");
        btnRun.addActionListener(this);
        btnCancel.setActionCommand("cancel");
        btnCancel.addActionListener(this);

        Box buttonBox = Box.createHorizontalBox();
        buttonBox.add(Box.createGlue());
        buttonBox.add(btnRun);
        buttonBox.add(Box.createHorizontalStrut(50));
        buttonBox.add(btnCancel);
        buttonBox.add(Box.createHorizontalStrut(50));
        buttonBox.add(btnHelpme);
        buttonBox.add(Box.createGlue());
        theContentPane.add(buttonBox, BorderLayout.SOUTH);
        this.pack();
        setLocationRelativeTo(getParent());
        this.setVisible(true);

    }

    public boolean useP3D() {
        return useP3D;
    }

    public MLMPlugin.CompressionType getCompressionType() {
        return compressionType;
    }

    public boolean useDiscrete() { return useDiscrete; }
    public boolean useRefProb() { return useRefprob; }
    public boolean useAlleleProb() { return useAlleleprob; }
    
    double getCompressionLevel() {
        double comp;
        try {
            comp = Double.parseDouble(txtCustom.getText());
        } catch (Exception e) {
            comp = Double.NaN;
        }
        return comp;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        if (e.getActionCommand().equals("run")) {
            runClicked = true;
            this.setVisible(false);
        } else if (e.getActionCommand().equals("cancel")) {
            this.setVisible(false);
        } else if (e.getActionCommand().equals("Optimum")) {
            compressionType = MLMPlugin.CompressionType.Optimum;
        } else if (e.getActionCommand().equals("Custom")) {
            compressionType = MLMPlugin.CompressionType.Custom;
        } else if (e.getActionCommand().equals("None")) {
            compressionType = MLMPlugin.CompressionType.None;
        } else if (e.getActionCommand().equals("Eachmarker")) {
            useP3D = false;
        } else if (e.getActionCommand().equals("P3D")) {
            useP3D = true;
        } else if (e.getActionCommand().equals("discrete")) {
            useDiscrete = true;
            useRefprob = false;
            useAlleleprob = false;
        } else if (e.getActionCommand().equals("reference")) {
            useDiscrete = false;
            useRefprob = true;
            useAlleleprob = false;
        } else if (e.getActionCommand().equals("allele")) {
            useDiscrete = false;
            useRefprob = false;
            useAlleleprob = true;
        }
    }
    
    public static void main(String[] args) {
    	MLMOptionDialog mod = new MLMOptionDialog(null, new boolean[]{true, true, true});
    	mod.dispose();
    }
}