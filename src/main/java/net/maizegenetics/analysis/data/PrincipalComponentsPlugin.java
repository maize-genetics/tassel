package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

import javax.swing.ImageIcon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.google.common.collect.Range;

import net.maizegenetics.analysis.numericaltransform.ImputationPlugin;
import net.maizegenetics.analysis.numericaltransform.NumericalGenotypePlugin;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.stats.PCA.PrinComp;
import net.maizegenetics.stats.PCA.PrinComp.PC_TYPE;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.SimpleTableReport;

public class PrincipalComponentsPlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(PrincipalComponentsPlugin.class);
    public static enum PCA_LIMIT {number_of_components, min_eigenvalue, total_variance};
    
    private PluginParameter<Boolean> useCovariance = new PluginParameter.Builder<>("covariance", true, Boolean.class)
    		.description("If the box is checked, then the analysis will do an eigenvalue decomposition of the covariance matrix. "
    				+ "If the box is unchecked, it will use a correlation matrix. Using the covariance matrix "
    				+ "is recommended for genotypes while the correlation matrix is often used for phenotypes.")
    		.guiName("covariance (alternative = correlation)")
    		.build();
    private PluginParameter<PCA_LIMIT> limitBy = new PluginParameter.Builder<>("limitBy", PCA_LIMIT.number_of_components, PCA_LIMIT.class)
    		.description("This parameter determines the type of value that will be used to limit the number of principal components (axes) returned. "
    				+ "The possible choices are number_of_components, min_eigenvalue, and total_variance.")
    		.guiName("limit number of components by")
    		.build();
    private PluginParameter<Integer> numberOfComponents = new PluginParameter.Builder<>("ncomponents", 5, Integer.class)
    		.description("The analysis will return this many principal components up to the number of taxa.")
    		.guiName("number of components")
    		.dependentOnParameter(limitBy, PCA_LIMIT.number_of_components)
    		.build();
    private PluginParameter<Double> minEigenval = new PluginParameter.Builder<>("minEigenval", 0.0, Double.class)
    		.description("All principal components with an eigenvalue greater than or equal to this value will be returned.")
    		.guiName("minimum eigenvalue")
    		.dependentOnParameter(limitBy, PCA_LIMIT.min_eigenvalue)
    		.build();
    private PluginParameter<Double> totalVar = new PluginParameter.Builder<>("totalVar", 0.5, Double.class)
    		.description("The first principal components that together explain this proportion of the total variance will be returned.")
    		.range(Range.closed(0.0, 1.0))
    		.guiName("total variance")
    		.dependentOnParameter(limitBy, PCA_LIMIT.total_variance)
    		.build();
    private PluginParameter<Boolean> reportEigenvalues = new PluginParameter.Builder<>("reportEigenvalues", true, Boolean.class)
    		.description("Returns a list of eigenvalues sorted high to low.")
    		.guiName("Return Eigenvalues")
    		.build();
    private PluginParameter<Boolean> reportEigenvectors = new PluginParameter.Builder<>("reportEigenvectors", true, Boolean.class)
    		.description("Returns the eigenvectors calculated from a Singular Value Decomposition of the data. The resulting table can be quite large if the number of variants and taxa are big.")
    		.guiName("Return Eigenvectors")
    		.build();
    
	public PrincipalComponentsPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
	}
	
	public DataSet processData(DataSet input){
		List<Datum> myResults = new ArrayList<>();
		List<Datum> myData = input.getDataOfType(new Class[]{Phenotype.class, GenotypeTable.class});
		for (Datum aDatum : myData) {
			if (aDatum.getData() instanceof Phenotype) {
				//check for missing values, throw an IllegalArgumentException if there is any missing data
				Phenotype myPhenotype = (Phenotype) aDatum.getData();
				if (areAnyPhenotypesMissing(myPhenotype.dataAttributeStream())) {
					StringBuilder msgBuilder = new StringBuilder();
					msgBuilder.append("There are missing values in ")
						.append(aDatum.getName())
						.append(". PCA will not be run.");
					throw new IllegalArgumentException(msgBuilder.toString());
				}

				//create the matrix, rows are observations, columns are data attributes
				List<PhenotypeAttribute> dataAttributes = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data);
				int nAttributes = dataAttributes.size();
				int nobs = myPhenotype.numberOfObservations();
				DoubleMatrix dataMatrix = DoubleMatrixFactory.DEFAULT.make(nobs, nAttributes);
				int colCount = 0;
				for (PhenotypeAttribute attr: dataAttributes) {
					float[] colData = ((NumericAttribute) attr).floatValues();
					for (int i = 0; i < nobs; i++) dataMatrix.set(i, colCount, colData[i]);
					colCount++;
				}

				//run PCA
				PC_TYPE pctype;
				if (useCovariance.value()) pctype = PC_TYPE.cov;
				else pctype = PC_TYPE.corr;
				PrinComp pca = new PrinComp(dataMatrix, pctype);
				
				//get results
				myResults.addAll(addResultsToDatumList(pca, myPhenotype.taxaAttribute(), dataAttributes, aDatum.getName()));
				
			} else {
				GenotypeTable myGenotype = (GenotypeTable) aDatum.getData();
				
				//is there a reference probability? If not, create one and impute missing values
				if (!myGenotype.hasReferenceProbablity()) {
					myGenotype = NumericalGenotypePlugin.setAlternateMinorAllelesToMinor(myGenotype);
					DataSet myDataset = new DataSet(new Datum("name", myGenotype, "comment"), this);
					ImputationPlugin myImputer = new ImputationPlugin(null, false);
					myImputer.by_mean(true);
					DataSet imputedDataset = myImputer.performFunction(myDataset);
					myGenotype = (GenotypeTable) imputedDataset.getData(0).getData();
				} else if (areAnyGenotypesMissingInReferenceProbability(myGenotype)) {
					StringBuilder msgBuilder = new StringBuilder();
					msgBuilder.append("There are missing values in ")
						.append(aDatum.getName())
						.append(". PCA will not be run.");
					throw new IllegalArgumentException(msgBuilder.toString());
				}
				
				//create the matrix, rows are taxa, columns are sites
				int ntaxa = myGenotype.numberOfTaxa();
				int nsites = myGenotype.numberOfSites();
				DoubleMatrix dataMatrix = DoubleMatrixFactory.DEFAULT.make(ntaxa, nsites);
				for (int t = 0; t < ntaxa; t++) {
					for (int s = 0; s < nsites; s++) {
						dataMatrix.set(t, s, myGenotype.referenceProbability(t, s));
					}
				}

				//run PCA
				PC_TYPE pctype;
				if (useCovariance.value()) pctype = PC_TYPE.cov;
				else pctype = PC_TYPE.corr;
				PrinComp pca = new PrinComp(dataMatrix, pctype);

				//get results
				myResults.addAll(addResultsToDatumList(pca, myGenotype, aDatum.getName()));
				
			}
		}

		return new DataSet(myResults, this);
	}
	
	private boolean areAnyPhenotypesMissing(Stream<NumericAttribute> attributes) {
		Optional<NumericAttribute> na = attributes.filter(a -> a.missing().cardinality() > 0).findAny();
		return na.isPresent();
	}
	
	private boolean areAnyGenotypesMissingInReferenceProbability(GenotypeTable myGenotype) {
		int ntaxa = myGenotype.numberOfTaxa();
		int nsites = myGenotype.numberOfSites();
		ReferenceProbability refprob = myGenotype.referenceProbability();
		for (int s = 0; s < nsites; s++) {
			for (int t = 0; t < ntaxa; t++) {
				if (Float.isNaN(refprob.value(t, s))) return true;
			}
		}
		return false;
	}
	
	private List<Datum> addResultsToDatumList(PrinComp pca, TaxaAttribute myTaxa, List<PhenotypeAttribute> data, String datasetName) {
		List<Datum> results = new ArrayList<>();
		
		//determine how many pc's to return
		int numberOfPCs;
		double[] eigenvalues = pca.getEigenValues();
		int nvalues = eigenvalues.length;
		double[] cumulativeEigenvalues = Arrays.copyOf(eigenvalues, nvalues);
		for (int i = 1; i < nvalues; i++) {
			cumulativeEigenvalues[i] += cumulativeEigenvalues[i - 1];
		}
		
		if (limitBy.value() == PCA_LIMIT.number_of_components) {
			numberOfPCs = Math.min(numberOfComponents.value(), nvalues);
		} else if (limitBy.value() == PCA_LIMIT.total_variance) {
			double limit = totalVar.value() * cumulativeEigenvalues[nvalues - 1];
			int ndx = Arrays.binarySearch(cumulativeEigenvalues, limit);
			if (ndx < -1) numberOfPCs = - ndx;
			else numberOfPCs = ndx + 1;
			numberOfPCs = Math.min(numberOfPCs, nvalues);
		} else {    //min_eigenvalue
			int ndx = Arrays.binarySearch(eigenvalues, minEigenval.value());
			if (ndx < -1) numberOfPCs = - ndx;
			else numberOfPCs = ndx + 1;
			numberOfPCs = Math.min(numberOfPCs, nvalues);
		}
		
		//create a Phenotype with the requested number of PCs
		DoubleMatrix pcs = pca.getPrincipalComponents();
		List<PhenotypeAttribute> attributes = new ArrayList<>();
		List<ATTRIBUTE_TYPE> types = new ArrayList<>();
		attributes.add(myTaxa);
		types.add(ATTRIBUTE_TYPE.taxa);
		int ntaxa = myTaxa.size();
		for (int i = 0; i < numberOfPCs; i++) {
			String pcname = "PC" + (i + 1);
			float[] pcvalue = new float[ntaxa];
			for (int t = 0; t < ntaxa; t++) {
				pcvalue[t] = (float) pcs.get(t, i);
			}
			NumericAttribute na = new NumericAttribute(pcname, pcvalue, new OpenBitSet(ntaxa));
			attributes.add(na);
			types.add(ATTRIBUTE_TYPE.covariate);
		}
		Phenotype pcPhenotype = new PhenotypeBuilder().fromAttributeList(attributes, types).build().get(0);
		StringBuilder nameBuilder = new StringBuilder();
		nameBuilder.append("PC_").append(datasetName);
		StringBuilder commentBuilder = new StringBuilder("\nPrincipalComponents stored as covariates.\n");
		commentBuilder.append("calculated from ").append(datasetName);
		
		results.add(new Datum(nameBuilder.toString(), pcPhenotype, commentBuilder.toString()));
		
		//create a tableReport with eigenvalues, if requested
		if (reportEigenvalues.value()) {
			String name = "Proportion of Variance Explained";
			String[] columnNames = new String[]{"PC","eigenvalue","proportion of total","cumulative proportion"};
			int nEigenvalues = eigenvalues.length;
			Object[][] tableData = new Object[nEigenvalues][4];
			double sumvalues = cumulativeEigenvalues[nEigenvalues - 1];
			for (int i = 0; i < nEigenvalues; i++) {
				tableData[i][0] = String.format("%d",i+1);
				tableData[i][1] = Double.valueOf(eigenvalues[i]);
				tableData[i][2] = new Double(eigenvalues[i]/sumvalues);
				tableData[i][3] = new Double(cumulativeEigenvalues[i]/sumvalues);
			}
			
			nameBuilder = new StringBuilder();
			nameBuilder.append("Eigenvalues_").append(datasetName);
			commentBuilder = new StringBuilder("\nEigenvalues and proportion of variance explained by PCs.\n");
			commentBuilder.append("calculated from ").append(datasetName);
			SimpleTableReport str = new SimpleTableReport(name, columnNames, tableData);
			results.add(new Datum(nameBuilder.toString(), str, commentBuilder.toString()));
		}
		
		//create a tableReport with eigenvectors, if requested
		DoubleMatrix eigenvectors = pca.getEigenVectors();
		if (reportEigenvectors.value()) {
			String name = "Eigenvectors";
			int ncol = numberOfPCs + 1;
			int nrows = data.size();
			String[] columnNames = new String[ncol];
			columnNames[0] = "Trait";
			for (int c = 1; c < ncol; c++) columnNames[c] = String.format("Eigenvector%d",c);
			Object[][] tableData = new Object[nrows][ncol];
			for (int r = 0; r < nrows; r++) {
				tableData[r][0] = data.get(r).name();
				for (int c = 1; c < ncol; c++) tableData[r][c] = new Double(eigenvectors.get(r, c - 1));
			}
			
			nameBuilder = new StringBuilder();
			nameBuilder.append("Eigenvectors_").append(datasetName);
			commentBuilder = new StringBuilder("\nEigenvectors for requested PCs.\n");
			commentBuilder.append("calculated from ").append(datasetName);
			SimpleTableReport str = new SimpleTableReport(name, columnNames, tableData);
			results.add(new Datum(nameBuilder.toString(), str, commentBuilder.toString()));
			
		}
		
		return results;
	}
	
	private List<Datum> addResultsToDatumList(PrinComp pca, GenotypeTable myGenotype, String datasetName) {
		List<Datum> results = new ArrayList<>();
		
		//determine how many pc's to return
		int numberOfPCs;
		double[] eigenvalues = pca.getEigenValues();
		int nvalues = eigenvalues.length;
		double[] cumulativeEigenvalues = Arrays.copyOf(eigenvalues, nvalues);
		for (int i = 1; i < nvalues; i++) {
			cumulativeEigenvalues[i] += cumulativeEigenvalues[i - 1];
		}
		
		if (limitBy.value() == PCA_LIMIT.number_of_components) {
			numberOfPCs = Math.min(numberOfComponents.value(), nvalues);
		} else if (limitBy.value() == PCA_LIMIT.total_variance) {
			double limit = totalVar.value() * cumulativeEigenvalues[nvalues - 1];
			int ndx = Arrays.binarySearch(cumulativeEigenvalues, limit);
			if (ndx < -1) numberOfPCs = - ndx;
			else numberOfPCs = ndx + 1;
			numberOfPCs = Math.min(numberOfPCs, nvalues);
		} else {    //min_eigenvalue
			int ndx = Arrays.binarySearch(eigenvalues, minEigenval.value());
			if (ndx < -1) numberOfPCs = - ndx;
			else numberOfPCs = ndx + 1;
			numberOfPCs = Math.min(numberOfPCs, nvalues);
		}
		
		//create a Phenotype with the requested number of PCs
		DoubleMatrix pcs = pca.getPrincipalComponents();
		List<PhenotypeAttribute> attributes = new ArrayList<>();
		List<ATTRIBUTE_TYPE> types = new ArrayList<>();
		attributes.add(new TaxaAttribute(myGenotype.taxa()));
		types.add(ATTRIBUTE_TYPE.taxa);
		int ntaxa = myGenotype.numberOfTaxa();
		for (int i = 0; i < numberOfPCs; i++) {
			String pcname = "PC" + (i + 1);
			float[] pcvalue = new float[ntaxa];
			for (int t = 0; t < ntaxa; t++) {
				pcvalue[t] = (float) pcs.get(t, i);
			}
			NumericAttribute na = new NumericAttribute(pcname, pcvalue, new OpenBitSet(ntaxa));
			attributes.add(na);
			types.add(ATTRIBUTE_TYPE.covariate);
		}
		Phenotype pcPhenotype = new PhenotypeBuilder().fromAttributeList(attributes, types).build().get(0);
		StringBuilder nameBuilder = new StringBuilder();
		nameBuilder.append("PC_").append(datasetName);
		StringBuilder commentBuilder = new StringBuilder("\nPrincipalComponents stored as covariates.\n");
		commentBuilder.append("calculated from ").append(datasetName);
		
		results.add(new Datum(nameBuilder.toString(), pcPhenotype, commentBuilder.toString()));
		
		//create a tableReport with eigenvalues, if requested
		if (reportEigenvalues.value()) {
			String name = "Proportion of Variance Explained";
			String[] columnNames = new String[]{"PC","eigenvalue","proportion of total","cumulative proportion"};
			int nEigenvalues = eigenvalues.length;
			Object[][] tableData = new Object[nEigenvalues][4];
			double sumvalues = cumulativeEigenvalues[nEigenvalues - 1];
			for (int i = 0; i < nEigenvalues; i++) {
				tableData[i][0] = String.format("%d",i);
				tableData[i][1] = new Double(eigenvalues[i]);
				tableData[i][2] = new Double(eigenvalues[i]/sumvalues);
				tableData[i][3] = new Double(cumulativeEigenvalues[i]/sumvalues);
			}
			
			nameBuilder = new StringBuilder();
			nameBuilder.append("Eigenvalues_").append(datasetName);
			commentBuilder = new StringBuilder("\nEigenvalues and proportion of variance explained by PCs.\n");
			commentBuilder.append("calculated from ").append(datasetName);
			SimpleTableReport str = new SimpleTableReport(name, columnNames, tableData);
			results.add(new Datum(nameBuilder.toString(), str, commentBuilder.toString()));
		}
		
		//create a tableReport with eigenvectors, if requested
		DoubleMatrix eigenvectors = pca.getEigenVectors();
		if (reportEigenvectors.value()) {
			String name = "Eigenvectors";
			int ncol = numberOfPCs + 1;
			int nrows = myGenotype.numberOfSites();
			String[] columnNames = new String[ncol];
			columnNames[0] = "Trait";
			for (int c = 1; c < ncol; c++) columnNames[c] = String.format("Eigenvector%d",c);
			Object[][] tableData = new Object[nrows][ncol];
			for (int r = 0; r < nrows; r++) {
				tableData[r][0] = myGenotype.positions().siteName(r);
				for (int c = 1; c < ncol; c++) tableData[r][c] = new Double(eigenvectors.get(r, c - 1));
			}
			
			nameBuilder = new StringBuilder();
			nameBuilder.append("Eigenvectors_").append(datasetName);
			commentBuilder = new StringBuilder("\nEigenvectors for requested PCs.\n");
			commentBuilder.append("calculated from ").append(datasetName);
			SimpleTableReport str = new SimpleTableReport(name, columnNames, tableData);
			results.add(new Datum(nameBuilder.toString(), str, commentBuilder.toString()));
			
		}
		
		return results;
	}

	@Override
	public String pluginDescription() {
		return "This plugin performs principal components analysis and returns the requested number of PC axes (components), and, optionally, the eigenvalues and eigenvectors. "
				+ "It can take as input either phenotype data or ReferenceProbability from a GenotypeTable.";
	}

	@Override
	public ImageIcon getIcon() {
        URL imageURL = FileLoadPlugin.class.getResource("/net/maizegenetics/analysis/images/pca.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
	}

	@Override
	public String getButtonName() {
		return "PCA";
	}

	@Override
	public String getToolTipText() {
		return "Performs principal components analysis";
	}

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(PrincipalComponentsPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    public Phenotype runPlugin(DataSet input) {
        return (Phenotype) performFunction(input).getData(0).getData();
    }

    /**
     * If the box is checked, then the analysis will do an
     * eigenvalue decomposition of the covariance matrix.
     * If the box is unchecked, it will use a correlation
     * matrix. Using the covariance matrix is recommended
     * for genotypes while the correlation matrix is often
     * used for phenotypes.
     *
     * @return covariance (alternative = correlation)
     */
    public Boolean covariance() {
        return useCovariance.value();
    }

    /**
     * Set covariance (alternative = correlation). If the
     * box is checked, then the analysis will do an eigenvalue
     * decomposition of the covariance matrix. If the box
     * is unchecked, it will use a correlation matrix. Using
     * the covariance matrix is recommended for genotypes
     * while the correlation matrix is often used for phenotypes.
     *
     * @param value covariance (alternative = correlation)
     *
     * @return this plugin
     */
    public PrincipalComponentsPlugin covariance(Boolean value) {
        useCovariance = new PluginParameter<>(useCovariance, value);
        return this;
    }

    /**
     * This parameter determines the type of value that will
     * be used to limit the number of principal components
     * (axes) returned. The possible choices are number_of_components,
     * min_eigenvalue, and total_variance.
     *
     * @return limit number of components by
     */
    public PCA_LIMIT limitNumberOfComponentsBy() {
        return limitBy.value();
    }

    /**
     * Set limit number of components by. This parameter determines
     * the type of value that will be used to limit the number
     * of principal components (axes) returned. The possible
     * choices are number_of_components, min_eigenvalue, and
     * total_variance.
     *
     * @param value limit number of components by
     *
     * @return this plugin
     */
    public PrincipalComponentsPlugin limitNumberOfComponentsBy(PCA_LIMIT value) {
        limitBy = new PluginParameter<>(limitBy, value);
        return this;
    }

    /**
     * The analysis will return this many principal components
     * up to the number of taxa.
     *
     * @return number of components
     */
    public Integer numberOfComponents() {
        return numberOfComponents.value();
    }

    /**
     * Set number of components. The analysis will return
     * this many principal components up to the number of
     * taxa.
     *
     * @param value number of components
     *
     * @return this plugin
     */
    public PrincipalComponentsPlugin numberOfComponents(Integer value) {
        numberOfComponents = new PluginParameter<>(numberOfComponents, value);
        return this;
    }

    /**
     * All principal components with an eigenvalue greater
     * than or equal to this value will be returned.
     *
     * @return minimum eigenvalue
     */
    public Double minimumEigenvalue() {
        return minEigenval.value();
    }

    /**
     * Set minimum eigenvalue. All principal components with
     * an eigenvalue greater than or equal to this value will
     * be returned.
     *
     * @param value minimum eigenvalue
     *
     * @return this plugin
     */
    public PrincipalComponentsPlugin minimumEigenvalue(Double value) {
        minEigenval = new PluginParameter<>(minEigenval, value);
        return this;
    }

    /**
     * The first principal components that together explain
     * this proportion of the total variance will be returned.
     *
     * @return total variance
     */
    public Double totalVariance() {
        return totalVar.value();
    }

    /**
     * Set total variance. The first principal components
     * that together explain this proportion of the total
     * variance will be returned.
     *
     * @param value total variance
     *
     * @return this plugin
     */
    public PrincipalComponentsPlugin totalVariance(Double value) {
        totalVar = new PluginParameter<>(totalVar, value);
        return this;
    }

    /**
     * Returns a list of eigenvalues sorted high to low.
     *
     * @return Return Eigenvalues
     */
    public Boolean returnEigenvalues() {
        return reportEigenvalues.value();
    }

    /**
     * Set Return Eigenvalues. Returns a list of eigenvalues
     * sorted high to low.
     *
     * @param value Return Eigenvalues
     *
     * @return this plugin
     */
    public PrincipalComponentsPlugin returnEigenvalues(Boolean value) {
        reportEigenvalues = new PluginParameter<>(reportEigenvalues, value);
        return this;
    }

    /**
     * Returns the eigenvectors calculated from a Singular
     * Value Decomposition of the data. The resulting table
     * can be quite large if the number of variants and taxa
     * are big.
     *
     * @return Return Eigenvectors
     */
    public Boolean returnEigenvectors() {
        return reportEigenvectors.value();
    }

    /**
     * Set Return Eigenvectors. Returns the eigenvectors calculated
     * from a Singular Value Decomposition of the data. The
     * resulting table can be quite large if the number of
     * variants and taxa are big.
     *
     * @param value Return Eigenvectors
     *
     * @return this plugin
     */
    public PrincipalComponentsPlugin returnEigenvectors(Boolean value) {
        reportEigenvectors = new PluginParameter<>(reportEigenvectors, value);
        return this;
    }


}
