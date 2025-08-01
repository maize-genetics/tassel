package net.maizegenetics.analysis.modelfitter;

import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.LinearModelUtils;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffectUtils;
import net.maizegenetics.stats.linearmodels.NestedCovariateModelEffect;
import net.maizegenetics.stats.linearmodels.PartitionedLinearModel;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.SimpleTableReport;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.analysis.association.AssociationUtils;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.stats.linearmodels.BasicShuffler;
import net.maizegenetics.taxa.Taxon;

public class StepwiseOLSModelFitter {
	private final GenotypePhenotype myData;

	//settable parameters
	private double[] enterlimits = null;
	private double[] exitlimits = null;
	private double enterlimit = 1e-5;
	private double exitlimit = 2e-5;
	private int maxNumberOfMarkers = 1000;
	private FactorModelEffect nestingEffect;
	private int nestingFactorIndex;
	private boolean isNested;
    private boolean calculateVIF = true;
    private double VIFTolerance = 0.001;
    private VIF_TYPE VIFType = VIF_TYPE.average; // options are min, average, and max
    public enum VIF_TYPE {min, average, max};
	private ArrayList<String> nestingFactorNames;

	//global variables used by the analysis
	private ArrayList<ModelEffect> currentModel;
	private int currentPhenotypeIndex;
	private int numberOfBaseModelEffects;
	private double[] y;
	private OpenBitSet missing;
	int numberNotMissing;
	int totalNumber;
	private LinkedList<Object[]> resultRowsAnova = new LinkedList<Object[]>();
	private LinkedList<Object[]> resultRowsAnovaWithCI = new LinkedList<Object[]>();        
	private LinkedList<Object[]> rowsSiteEffectTable = new LinkedList<Object[]>();
	private LinkedList<Object[]> rowsSiteEffectTableWithCI = new LinkedList<Object[]>();
    private ArrayList<Integer> excludedSNPs = new ArrayList<>();
	private MODEL_TYPE modelType = MODEL_TYPE.mbic;
	private String datasetName;
	private double globalbestbic = Double.MAX_VALUE; //Rationale: the optimal model has minimum BIC. Thus, initialize bestbic with the largest possible double java can handle.
	private double globalbestmbic = Double.MAX_VALUE;
	private double globalbestaic = Double.MAX_VALUE;
	public enum MODEL_TYPE {pvalue, bic, mbic, aic};
	private boolean test;
    private int[][] theUpperAndLowerBound;

    private int numberOfPermutations = 0;
    private double alpha = 0.05;
    private double[] minPvalues = null;

	private final String[] anovaReportHeader = new String[]{"Trait", "Name","Locus","Position","df","SS","MS", "F", "pr>F", "BIC", "mBIC", "AIC", "Model Rsq"};
    private final String[] anovaReportWithCIHeader = new String[]{"Trait", "Name","Locus","Position","df","SS","MS", "F", "pr>F", "BIC", "mBIC", "AIC", "Model Rsq", "SuppLeft", "SuppRight"};
    
	private final GenotypeTable myGenotype;
	private final Phenotype myPhenotype;
	List<PhenotypeAttribute> dataAttributeList;
	List<PhenotypeAttribute> factorAttributeList;
	List<PhenotypeAttribute> covariateAttributeList;

	public void setModelType(MODEL_TYPE modelType){
		this.modelType = modelType;
	}

	public StepwiseOLSModelFitter(GenotypePhenotype genoPheno, String datasetName) {
		myData = genoPheno;       
		myGenotype = myData.genotypeTable();
        myPhenotype = myData.phenotype();

		this.datasetName = datasetName;
	}

	public DataSet runAnalysis() {
		//numbers of various model components
		dataAttributeList = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data);
		factorAttributeList = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor);
		covariateAttributeList = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.covariate);

		int numberOfPhenotypes = dataAttributeList.size();

		//cycle through the phenotypes
		//notation: 
		//X is the design matrix without the markers, rows of X will be deleted if marker data is missing
		//Xm is the design matrix with markers
		//y is the data
		currentPhenotypeIndex = -1;
		for (PhenotypeAttribute attr : dataAttributeList) {
			NumericAttribute currentTrait = (NumericAttribute) attr;
			currentPhenotypeIndex++;
			if (enterlimits != null) enterlimit = enterlimits[currentPhenotypeIndex];
			if (exitlimits != null) exitlimit = exitlimits[currentPhenotypeIndex];
			globalbestbic = Double.MAX_VALUE; 
			globalbestmbic = Double.MAX_VALUE;
			globalbestaic = Double.MAX_VALUE;
			
			//get phenotype data
			float[] phenotypeData = currentTrait.floatValues();

			//keep track of missing rows
			missing = new OpenBitSet(currentTrait.missing());
			for (PhenotypeAttribute factor : factorAttributeList) missing.or(factor.missing());
			for (PhenotypeAttribute cov : covariateAttributeList) missing.or(cov.missing());

			//remove missing values from the arrays
			totalNumber = phenotypeData.length;
			numberNotMissing = totalNumber - (int) missing.cardinality();

			y = AssociationUtils.getNonMissingDoubles(phenotypeData, missing);
			fitModel();
			scanAndFindCI();    
		}

		return null;
	}

	public void fitModel() {
		//build the base model
		currentModel = new ArrayList<ModelEffect>();
		int numberOfTaxa = y.length;
		int[] mean = new int[numberOfTaxa];

		FactorModelEffect meanEffect = new FactorModelEffect(mean, false);
		meanEffect.setID("mean");
		currentModel.add(meanEffect);

		for (PhenotypeAttribute factor:factorAttributeList) {
			ArrayList<String> ids = new ArrayList<String>();
			CategoricalAttribute ca = (CategoricalAttribute) factor;
			String[] factorLabels = AssociationUtils.getNonMissingValues(ca.allLabels(), missing);
			int[] levels = ModelEffectUtils.getIntegerLevels(factorLabels, ids);
			FactorModelEffect fme = new FactorModelEffect(levels, true, new Object[]{ca.name(), ids});
			if (isNested && myPhenotype.attributeIndexForName(factor.name()) == nestingFactorIndex) {
				nestingEffect = fme;
				nestingFactorNames = ids; 
			}
			currentModel.add(fme);
		}
		
		//add the covariate effects
		for (PhenotypeAttribute cov:covariateAttributeList) {
			NumericAttribute na = (NumericAttribute) cov;
			double[] covData = AssociationUtils.getNonMissingDoubles(na.floatValues(), missing);
			CovariateModelEffect cme = new CovariateModelEffect(covData, cov.name());
			currentModel.add(cme);
		}
		numberOfBaseModelEffects = currentModel.size();
		
		//run permutation test
		System.out.println("-----Number of Permutations = " + numberOfPermutations + "------------");
		if(numberOfPermutations > 0) runPermutationTestNoMissingData(); 
		System.out.println("--------------the Enter limit is " + enterlimit);
		System.out.println("--------------the Exit limit is " + exitlimit);  			

		while(forwardStep()){
			if((calculateVIF)&&(currentModel.size()> (numberOfBaseModelEffects+1))){
				ArrayList<Double> theVIFValues = removeCollinearMarkers(true);
			}
			while(backwardStep());
		}
		if((calculateVIF)&&(currentModel.size()> (numberOfBaseModelEffects+1))){
			//Begin code for testing out the model
			ArrayList<SNP> theSNPs = new ArrayList<SNP>();
			for(int i = numberOfBaseModelEffects; i < currentModel.size();i++){
				SNP thisSNP = (SNP) currentModel.get(i).getID();
				theSNPs.add(thisSNP);
			}

			System.out.println("Here are all of the SNPs in the model prior to identifying a collinear effect: ");
			System.out.println(theSNPs.toString());
			//End code for testing out the model

			ArrayList<Double> theVIFValues = removeCollinearMarkers(true);
		}

		SweepFastLinearModel theLinearModel = new SweepFastLinearModel(currentModel, y);
		appendAnovaResults(theLinearModel);
		appendSiteEffectEstimates(theLinearModel);

	}

	public ArrayList<Double> removeCollinearMarkers(boolean removeCollinearMarker){
		//Obtain the design matrix of the current model (everything except for the "BaseModelEffects")        
		DoubleMatrix theDesignMatrix = null;
		for(int i=numberOfBaseModelEffects; i < currentModel.size(); i++){//I am starting the i loop at 1 because I do not want to have the intercept in the model
			DoubleMatrix partialDesignMatrix = currentModel.get(i).getX();
			if(i==numberOfBaseModelEffects){
				theDesignMatrix = partialDesignMatrix;
			}else{
				theDesignMatrix = theDesignMatrix.concatenate(partialDesignMatrix, false);
			}
		}

		//Scale and center the design matrix
		DoubleMatrix theCenteredDesignmatrix = centerCols(theDesignMatrix);         
		DoubleMatrix theCenteredAndScaledDesignMatrix = scaleCenteredMatrix(theCenteredDesignmatrix);
		double invSqrtSampleSizeMinusOne = 1/Math.sqrt(((double)theCenteredAndScaledDesignMatrix.numberOfRows())-1);
		theCenteredAndScaledDesignMatrix = theCenteredAndScaledDesignMatrix.scalarMult(invSqrtSampleSizeMinusOne);

		//Obtain the correlation matrix of the explanatory variables, r_{xx}
		DoubleMatrix theCorrelationMatrix = theCenteredAndScaledDesignMatrix.crossproduct(); 


		//Invert r_{xx}, and call it r_{xx}^{-1}. The diagonal elements of r_{xx}^{-1} are the VIF values
		DoubleMatrix theInvertedCorrelationMatrix = theCorrelationMatrix.inverse();
		ArrayList<Double> theVIFForEachExplanatoryVariable = new ArrayList<Double>();

		//For loop through the marker effects //NOTE: Turn these two for loops into separate methods
		int counter = 0;
		double theMinimumValue = Double.MAX_VALUE;
		int indexOfCollinearEffect = 0;
		ArrayList<Integer> indiciesOfCollinearExplanatoryVariables = new ArrayList<Integer>();
		ModelEffect theCollinearEffect = null;
		boolean InfiniteVIFsPresent = false;
		for(int i=numberOfBaseModelEffects; i < currentModel.size(); i++){
			int numberOfExplanatoryVariables = currentModel.get(i).getX().numberOfColumns();

			//For loop through all explanatory variables for a marker effect
			double theSumVIFValues = 0; 
			double theMinimumToleranceValueWithinX = Double.MAX_VALUE;
			double theMaximumToleranceValueWithinX = Double.MIN_VALUE;
			for(int j = 0; j < numberOfExplanatoryVariables; j++){
				//If one of the diagonoal elements of theInvertedCorrelationMatrix (i.e. the VIFs)
				// is infinity, this means that the most recent term added to the model has has a VIF
				// of infinity. The following if statment will remove this term from the model.
				double theVIF = theInvertedCorrelationMatrix.get((counter+j), (counter+j));
				if(Double.isNaN(theVIF)){
					theMinimumValue = 0;
					theCollinearEffect = currentModel.get((currentModel.size()-1));
					InfiniteVIFsPresent = true;
					break;
				}
				double theInvertedVIF = 1/theInvertedCorrelationMatrix.get((counter+j), (counter+j));
				if((j == 0) | (j == numberOfExplanatoryVariables)) System.out.println("theInvertedVIF is: " + theInvertedVIF);
				if(theInvertedVIF < theMinimumToleranceValueWithinX)theMinimumToleranceValueWithinX = theInvertedVIF;
				if(theInvertedVIF > theMaximumToleranceValueWithinX)theMaximumToleranceValueWithinX = theInvertedVIF;
				theSumVIFValues += theInvertedVIF;
			}

			if(InfiniteVIFsPresent) break;
			if(VIFType == VIF_TYPE.average){
				//Take the average (1/VIF) value of the tested marker
				double theAverageVIFValue = theSumVIFValues/numberOfExplanatoryVariables;
				if(theAverageVIFValue < theMinimumValue){
					theMinimumValue = theAverageVIFValue;
					theCollinearEffect = currentModel.get(i);
					indexOfCollinearEffect = i-numberOfBaseModelEffects; //i-numberOfBaseModelEffects because java indexes things starting at zero
				}
				theVIFForEachExplanatoryVariable.add(theAverageVIFValue);
			}else if(VIFType == VIF_TYPE.min){
				if(theMinimumToleranceValueWithinX < theMinimumValue){
					theMinimumValue = theMinimumToleranceValueWithinX;
					theCollinearEffect = currentModel.get(i);
					indexOfCollinearEffect = i-numberOfBaseModelEffects;//i-1 because java indexes things starting at zero
				} 
				theVIFForEachExplanatoryVariable.add(theMinimumToleranceValueWithinX);
			}else if(VIFType == VIF_TYPE.max){
				if(theMaximumToleranceValueWithinX < theMinimumValue){
					theMinimumValue = theMaximumToleranceValueWithinX;
					theCollinearEffect = currentModel.get(i);
					indexOfCollinearEffect = i-numberOfBaseModelEffects;//i-1 because java indexes things starting at zero
				} 
				theVIFForEachExplanatoryVariable.add(theMaximumToleranceValueWithinX);

			}

			counter = counter + numberOfExplanatoryVariables;
			System.out.println("indexOfCollinearEffect: " + indexOfCollinearEffect);
			System.out.println("theVIFForEachExplanatoryVariable");
			System.out.println(theVIFForEachExplanatoryVariable.toString());

			//End for loop through the marker effects
		}
		//If this minimum exceeds the tolerence threshold, remove the corresponding marker from the model
		if(removeCollinearMarker){
			if(theMinimumValue < VIFTolerance){
				ArrayList<SNP> theSNPs = new ArrayList<SNP>();
				for(int i = numberOfBaseModelEffects; i < currentModel.size();i++){
					SNP thisSNP = (SNP) currentModel.get(i).getID();
					theSNPs.add(thisSNP);
				}
				System.out.println("Here are all of the SNPs in the model prior to identifying a collinear effect: ");
				System.out.println(theSNPs.toString());

				currentModel.remove(theCollinearEffect);
				excludedSNPs.add(((SNP)theCollinearEffect.getID()).index);

				SNP snpRemoved = (SNP) theCollinearEffect.getID(); 
				if(!InfiniteVIFsPresent)theVIFForEachExplanatoryVariable.remove(indexOfCollinearEffect);
				System.out.println("------------------"+snpRemoved +
						" is causing multicollinearity in the model. It has been removed from the model------------------------------");
			}
		}

		return theVIFForEachExplanatoryVariable;

	}    
	
	//scaleCenteredMatrix was written by in PrinComp by (I think) Peter Bradbury

	private DoubleMatrix centerCols(DoubleMatrix data) {
		int nrows = data.numberOfRows();
		int ncols = data.numberOfColumns();
		DoubleMatrix dm = data.copy();
		for (int c = 0; c < ncols; c++) {
			double colmean = dm.columnSum(c) / nrows;
			for (int r = 0; r < nrows; r++) {
				dm.set(r, c, dm.get(r, c) - colmean);
			}
		}

		return dm;
	}
	
	private DoubleMatrix scaleCenteredMatrix(DoubleMatrix data) {
		int nrows = data.numberOfRows();
		int ncols = data.numberOfColumns();
		for (int c = 0; c < ncols; c++) {
			double sumsq = 0;
			for (int r = 0; r < nrows; r++) {
				double val = data.get(r, c);
				sumsq += val * val;
			}
			double stdDev = Math.sqrt(sumsq/(nrows - 1));
			for (int r = 0; r < nrows; r++) {
				double val = data.get(r, c);
				data.set(r, c, val/stdDev);
			}
		}
		return data;
	}

    public void scanAndFindCI() {
        //rescan the model to refine position and find CI's
        //rescan requires a map

        //for each snp in the model:
        //1. add an adjacent snp to left of the original
        //2. solve the new model
        //3. if the marginal p of the model snp is sig stop and repeat on right of snp.
        //4. find the max ss of any point with the CI interval (not including the original marker)
        //5. if the max ss is greater than the original marker ss, replace the original and go to 1.
        //report steps as process proceeds
        SweepFastLinearModel sflm;
        int numberOfTerms = currentModel.size();
        int firstMarkerIndex = 1;
        firstMarkerIndex += factorAttributeList.size();
        firstMarkerIndex += covariateAttributeList.size();
        
        System.out.println("firstMarkerIndex " + firstMarkerIndex);
        int[] upperbound = new int[numberOfTerms - firstMarkerIndex];
        int[] lowerbound = new int[numberOfTerms - firstMarkerIndex];
        theUpperAndLowerBound = new int[numberOfTerms - firstMarkerIndex][2];
        for (int t = firstMarkerIndex; t < numberOfTerms; t++){ //0 is the mean, 1 is populations
            //find the CI bounds
            lowerbound[t - firstMarkerIndex] = scanASide(true, t);
            System.out.println("lowerbound[t - firstMarkerIndex]");
            System.out.println(lowerbound[t - firstMarkerIndex]);
            upperbound[t - firstMarkerIndex] = scanASide(false, t);
            System.out.println("upperbound[t - firstMarkerIndex]");
            System.out.println(upperbound[t - firstMarkerIndex]);
            theUpperAndLowerBound[t-firstMarkerIndex][0] = lowerbound[t-firstMarkerIndex];
            theUpperAndLowerBound[t-firstMarkerIndex][1] = upperbound[t-firstMarkerIndex];
             
            //find the marker with highest ss in the CI
            SNP bestsnp = null;
            double bestss = 0;
            ModelEffect besteffect = null;
            ModelEffect currentme = currentModel.remove(t);
            sflm = new SweepFastLinearModel(currentModel, y);
            PartitionedLinearModel plm = new PartitionedLinearModel(currentModel, sflm);
            
            			//NEW CODE: create the appropriate marker effect
            for (int m = lowerbound[t - firstMarkerIndex]; m <= upperbound[t - firstMarkerIndex]; m++) {
                ModelEffect markerEffect = null;
                SNP testsnp = new SNP(myGenotype.siteName(m), myGenotype.chromosome(m), myGenotype.chromosomalPosition(m), m);
                
                //if the Genotype has a reference allele use that, else if genotype use that else nothing, allele probability not implemented
                //for genotype use an additive model for now
                if (myGenotype.hasReferenceProbablity()) {
                	double[] cov = new double[numberNotMissing];
                	int count = 0;
                	float[] siteReferenceProb = myData.referenceProb(m);
                	for (int i = 0; i < totalNumber; i++) {
                		if (!missing.fastGet(i)) cov[count++] = siteReferenceProb[i];
                	}

                	if (isNested) {
                		CovariateModelEffect cme = new CovariateModelEffect(cov);
                		markerEffect = new NestedCovariateModelEffect(cme, nestingEffect);
                		markerEffect.setID(testsnp);
                	} else {
                		markerEffect = new CovariateModelEffect(cov, testsnp);
                	}
 
                } else if (myGenotype.hasGenotype()) {
                	byte minor = myGenotype.minorAllele(m);
                	double[] cov = new double[numberNotMissing];
                	byte[] siteGeno = AssociationUtils.getNonMissingBytes( myData.genotypeAllTaxa(m), missing);
                	for (int i = 0; i < numberNotMissing; i++) {
                		byte[] diploidAlleles = GenotypeTableUtils.getDiploidValues(siteGeno[i]);
                		if (diploidAlleles[0] == minor) cov[i] += 0.5;
                		if (diploidAlleles[1] == minor) cov[i] += 0.5;
                	}
                	
                	if (isNested) {
                		CovariateModelEffect cme = new CovariateModelEffect(cov);
                		markerEffect = new NestedCovariateModelEffect(cme, nestingEffect);
                		markerEffect.setID(testsnp);
                	} else {
                		markerEffect = new CovariateModelEffect(cov, testsnp);
                	}
               	
                } else {
                	throw new IllegalArgumentException("No genotypes or reference probabilities in the data set.");
                }
                
                plm.testNewModelEffect(markerEffect);
                double modelss = plm.getModelSS();
                if (modelss > bestss) {
                    bestss = modelss;
                    bestsnp = testsnp;
                    besteffect = markerEffect;
                }
            }
            
            currentModel.add(t, besteffect);
            ////
            //did the best marker change?
            boolean markerchanged = false;
            SNP currentSnp = (SNP) currentme.getID();
            if (currentSnp.index != bestsnp.index) markerchanged = true;

            //if this marker is different than the one tested, reset the CI
            if (markerchanged) {
                lowerbound[t - firstMarkerIndex] = scanASide(true, t);
                upperbound[t - firstMarkerIndex] = scanASide(false, t);
                theUpperAndLowerBound[t-firstMarkerIndex][0] = lowerbound[t-firstMarkerIndex];
                theUpperAndLowerBound[t-firstMarkerIndex][1] = upperbound[t-firstMarkerIndex];
            }
            System.out.println("upperAndLowerBound[t-firstMarkerIndex][0]");
            System.out.println(theUpperAndLowerBound[t-firstMarkerIndex][0]); 
            System.out.println("upperAndLowerBound[t-firstMarkerIndex][1]");
            System.out.println(theUpperAndLowerBound[t-firstMarkerIndex][1]);
        }
        
        //report the results including CI's
        //writeModelWithCI2File(lowerbound, upperbound); Turn this into a global matrix of upper and lower bounds
        SweepFastLinearModel theLinearModel = new SweepFastLinearModel(currentModel, y);
        appendAnovaResultsWithCI(theLinearModel);
        appendSiteEffectEstimatesWithCI(theLinearModel);

 
    }        
        
  
    
    private int scanASide(boolean left, int whichModelTerm) {

    	double alpha = 0.05;
    	int minIndex = 0;
    	int maxIndex = myGenotype.numberOfSites() - 1;
    	int incr;
    	if (left) {
    		incr = -1;
    	} else {
    		incr = 1;
    	}

    	SNP modelsnp = (SNP) currentModel.get(whichModelTerm).getID();
    	int markerIndex = modelsnp.index;
    	Chromosome chr = modelsnp.locus;
    	//int chrAsInt = chr.getChromosomeNumber();
    	boolean boundfound = false;
    	int testIndex = markerIndex;
    	int lastterm = currentModel.size();

    	///
    	do{
    		testIndex += incr;
    		if(testIndex < minIndex || testIndex > maxIndex){
    			testIndex -= incr;
    			boundfound = true;
    			break;
    		}
    		ModelEffect markerEffect = null;
    		SNP snp = new SNP(myGenotype.siteName(testIndex), myGenotype.chromosome(testIndex), myGenotype.chromosomalPosition(testIndex), testIndex);
    		Chromosome chrOfTestedSnp = snp.locus;
    		//int chrOfTestedSnpAsInt = chrOfTestedSnp.getChromosomeNumber();
    		if (snp == null || !chrOfTestedSnp.equals(chr)) {
    		//if (snp == null || chrOfTestedSnpAsInt != chrAsInt) {
    			testIndex -= incr;
    			boundfound = true;
    		} else {
    			//create the appropriate marker effect

    			
                //if the Genotype has a reference allele use that, else if genotype use that else nothing, allele probability not implemented
                //for genotype use an additive model for now
                if (myGenotype.hasReferenceProbablity()) {
                	double[] cov = new double[numberNotMissing];
                	int count = 0;
                	float[] referenceProb = myData.referenceProb(testIndex);
                	for (int i = 0; i < totalNumber; i++) {
                		if (!missing.fastGet(i)) cov[count++] = referenceProb[i];
                	}

                	if (isNested) {
                		CovariateModelEffect cme = new CovariateModelEffect(cov);
                		markerEffect = new NestedCovariateModelEffect(cme, nestingEffect);
                		markerEffect.setID(snp);
                	} else {
                		markerEffect = new CovariateModelEffect(cov, snp);
                	}
 
                } else if (myGenotype.hasGenotype()) {
                	byte minor = myGenotype.minorAllele(testIndex);
                	double[] cov = new double[numberNotMissing];
                	byte[] siteGeno = AssociationUtils.getNonMissingBytes(myData.genotypeAllTaxa(testIndex), missing);
                	for (int i = 0; i < numberNotMissing; i++) {
                		byte[] diploidAlleles = GenotypeTableUtils.getDiploidValues(siteGeno[i]);
                		if (diploidAlleles[0] == minor) cov[i] += 0.5;
                		if (diploidAlleles[1] == minor) cov[i] += 0.5;
                	}
                	
                	if (isNested) {
                		CovariateModelEffect cme = new CovariateModelEffect(cov);
                		markerEffect = new NestedCovariateModelEffect(cme, nestingEffect);
                		markerEffect.setID(snp);
                	} else {
                		markerEffect = new CovariateModelEffect(cov, snp);
                	}
               	
                } else {
                	throw new IllegalArgumentException("No genotypes or reference probabilities in the data set.");
                }

    			///
    			currentModel.add(markerEffect);
    			SweepFastLinearModel sflm = new SweepFastLinearModel(currentModel, y);
    			double[] snpssdf = sflm.getMarginalSSdf(whichModelTerm);
    			double[] errorssdf = sflm.getResidualSSdf();
    			double F = snpssdf[0] / snpssdf[1] / errorssdf[0] * errorssdf[1];
    			double p;
    			try {
    				p = LinearModelUtils.Ftest(F, snpssdf[1], errorssdf[1]);
    			} catch(Exception e) {
    				p = 1;
    			}     

    			if(p < alpha){
    				boundfound = true;
    			}
    			currentModel.remove(lastterm);
    		}    
    	}  while (!boundfound && testIndex > minIndex && testIndex < maxIndex);      
    	////      
    	return testIndex;
    }      
        
	public boolean forwardStep() {
		double bestss = 0;
		double bestbic = Double.MAX_VALUE;
		double bestmbic = Double.MAX_VALUE;
		double bestaic = Double.MAX_VALUE;
		ModelEffect besteffect = null;

		SweepFastLinearModel sflm = new SweepFastLinearModel(currentModel, y);
		PartitionedLinearModel plm = new PartitionedLinearModel(currentModel, sflm);
		int numberOfSites = myGenotype.numberOfSites();

		System.out.println("We are in forwardStep()");

		double[] temperrorssdf = sflm.getResidualSSdf();
		System.out.println("The value of errorss before the loop in forwardStep() is: " + temperrorssdf[0]);
		for (int s = 0; s < numberOfSites; s++) if (!excludedSNPs.contains(s)) {
			//create the appropriate marker effect
			ModelEffect markerEffect = null;
			SNP snp = new SNP(myGenotype.siteName(s), myGenotype.chromosome(s), myGenotype.chromosomalPosition(s), s);

            //if the Genotype has a reference allele use that, else if genotype use that else nothing, allele probability not implemented
            //for genotype use an additive model for now
            if (myGenotype.hasReferenceProbablity()) {
            	double[] cov = new double[numberNotMissing];
            	int count = 0;
            	float[] referenceProb = myData.referenceProb(s);
            	for (int i = 0; i < totalNumber; i++) {
            		if (!missing.fastGet(i)) cov[count++] = referenceProb[i];
            	}

            	if (isNested) {
            		CovariateModelEffect cme = new CovariateModelEffect(cov);
            		markerEffect = new NestedCovariateModelEffect(cme, nestingEffect);
            		markerEffect.setID(snp);
            	} else {
            		markerEffect = new CovariateModelEffect(cov, snp);
            	}

            } else if (myGenotype.hasGenotype()) {
            	byte minor = myGenotype.minorAllele(s);
            	double[] cov = new double[numberNotMissing];
            	byte[] siteGeno = AssociationUtils.getNonMissingBytes(myData.genotypeAllTaxa(s), missing);
            	for (int i = 0; i < numberNotMissing; i++) {
            		byte[] diploidAlleles = GenotypeTableUtils.getDiploidValues(siteGeno[i]);
            		if (diploidAlleles[0] == minor) cov[i] += 0.5;
            		if (diploidAlleles[1] == minor) cov[i] += 0.5;
            	}
            	
            	if (isNested) {
            		CovariateModelEffect cme = new CovariateModelEffect(cov);
            		markerEffect = new NestedCovariateModelEffect(cme, nestingEffect);
            		markerEffect.setID(snp);
            	} else {
            		markerEffect = new CovariateModelEffect(cov, snp);
            	}
           	
            } else {
            	throw new IllegalArgumentException("No genotypes or reference probabilities in the data set.");
            }

			plm.testNewModelEffect(markerEffect);
			double modelss = plm.getModelSS();

			currentModel.add(markerEffect); //Temporary; Remove if wrong
			SweepFastLinearModel sflm2 = new SweepFastLinearModel(currentModel, y);
			
			int n = numberNotMissing;
			//Calculate the BIC
			double [] errorss = sflm2.getResidualSSdf();
			double [] modeldf = sflm2.getFullModelSSdf();
			double pForBIC  = modeldf[1]; 
			double bic = (n*Math.log(errorss[0]/n)) + (pForBIC*Math.log(n)) ;
			double aic = (n*Math.log(errorss[0]/n)) + (2*pForBIC) ;
			if(s==2) System.out.println("The value of errorss in forwardStep() is: " + errorss[0]);

			//Calculate the mBIC
			int numberOfTwoWayInteractions = (numberOfSites*(numberOfSites-1))/2;
			double pForMbic = modeldf[1];
			double qForMbic = 0; //This is going to be the number of two-way interaction terms. For now, this is not implemented, and thus I am setting it equal to zer
			double mbic;
			if(qForMbic == 0){
				mbic = (n*Math.log(errorss[0])) + ((pForMbic+qForMbic)*Math.log(n)) + (2*pForMbic*(Math.log((numberOfSites/2.2)-1)));
			}else{
				mbic = (n*Math.log(errorss[0])) + ((pForMbic+qForMbic)*Math.log(n)) + (2*pForMbic*(Math.log((numberOfSites/2.2)-1))) + 
						(2*qForMbic*(Math.log((numberOfTwoWayInteractions/2.2)-1) )); 
			}


			switch(modelType){
			case pvalue:
				test = modelss > bestss;
				break;
			case bic:
				test = bic < bestbic;
				break;
			case mbic:
				test = mbic < bestmbic;
				break;
			case aic:
				test = aic < bestaic;
				break;
			}
			if (test) {
				bestmbic = mbic;
				bestbic = bic;
				bestaic = aic;
				bestss = modelss;
				besteffect = markerEffect;
			}
			currentModel.remove(markerEffect); //Temporary; Remove if wrong
		}

		//if the p-value for the select SNP is less than the enter limit, add it to the model and recalculate the model solution
		plm.testNewModelEffect(besteffect);
		double[] Fp = plm.getFp();
		if(modelType == MODEL_TYPE.pvalue){
			if ( Fp[1] < enterlimit) {
				currentModel.add(besteffect);
				if (currentModel.size() == maxNumberOfMarkers + numberOfBaseModelEffects) return false;
				return true;
			} else { 
				return false; 
			}     
		} else if(modelType == MODEL_TYPE.bic){
			if(bestbic < globalbestbic){
				globalbestbic = bestbic;
				currentModel.add(besteffect);
				if (currentModel.size() == maxNumberOfMarkers + numberOfBaseModelEffects) return false;
				return true;
			} else{
				return false;
			} 
		} else if(modelType == MODEL_TYPE.mbic){
			if(bestmbic < globalbestmbic){
				globalbestmbic = bestmbic;
				//System.out.println("***********The value of globalbestmbic at the end of forwardStep() is " + globalbestmbic);
				currentModel.add(besteffect);
				if (currentModel.size() == maxNumberOfMarkers + numberOfBaseModelEffects) return false;
				return true;
			} else{
				return false;
			}
		} else if(modelType == MODEL_TYPE.aic){
			if(bestaic < globalbestaic){
				globalbestaic = bestaic;
				//System.out.println("***********The value of globalbestmbic at the end of forwardStep() is " + globalbestmbic);
				currentModel.add(besteffect);
				if (currentModel.size() == maxNumberOfMarkers + numberOfBaseModelEffects) return false;
				return true;
			} else{
				return false;
			} 
		}  else{
			return false;
		}  


	}

	public boolean backwardStep() {
		int numberOfTerms = currentModel.size();
		if (numberOfTerms <= numberOfBaseModelEffects) return false;
		double bestbic = Double.MAX_VALUE;
		double bestmbic = Double.MAX_VALUE;
		double bestaic = Double.MAX_VALUE;

		int n = y.length;
		int numberOfSites = myGenotype.numberOfSites();

		SweepFastLinearModel sflm0 = new SweepFastLinearModel(currentModel, y);

		//find the model term (snps only) with the largest p-value
		double maxp = 0;
		double minF= -1;
		int maxterm = 0;
		ModelEffect worstMarkerEffect = null;
		double[] errorssdf = sflm0.getResidualSSdf();

		for (int t = numberOfBaseModelEffects; t < numberOfTerms; t++) {
			double bic;
			double mbic;
			double aic;
			ModelEffect markerEffect = null;
			double[] termssdf = sflm0.getIncrementalSSdf(t);
			double F = termssdf[0]/termssdf[1]/errorssdf[0]*errorssdf[1];
			double p;
			if(modelType != MODEL_TYPE.pvalue){

				ModelEffect meTest1 = currentModel.remove(t);
				SNP snpRemoved = (SNP) meTest1.getID(); 
				System.out.println("CORRECT The SNP just removed is: " + snpRemoved);

				SweepFastLinearModel sflm = new SweepFastLinearModel(currentModel, y); 
				//Calculate the BIC
				double [] errorss = sflm.getResidualSSdf();
				double [] modeldf = sflm.getFullModelSSdf();
				double pForBIC  = modeldf[1]; 
				bic = (n*Math.log(errorss[0]/n)) + (pForBIC*Math.log(n)) ; 
				aic = (n*Math.log(errorss[0]/n)) + (2*pForBIC) ;

				//Calculate the mBIC
				double numberOfTwoWayInteractions = (double) (numberOfSites*(numberOfSites-1))/2;			
				double pForMbic = modeldf[1];
				double qForMbic = 0; //This is going to be the number of two-way interaction terms. For now, this is not implemented, and thus I am setting it equal to zer
				if(qForMbic == 0){
					mbic = (n*Math.log(errorss[0])) + ((pForMbic+qForMbic)*Math.log(n)) + (2*pForMbic*(Math.log((numberOfSites/2.2)-1)));
				}else{
					mbic = (n*Math.log(errorss[0])) + ((pForMbic+qForMbic)*Math.log(n)) + (2*pForMbic*(Math.log((numberOfSites/2.2)-1))) + 
							(2*qForMbic*(Math.log((numberOfTwoWayInteractions/2.2)-1) )); 
				}


				currentModel.add(t, meTest1);

				sflm = new SweepFastLinearModel(currentModel, y); 


			} else{
				bic = Double.MAX_VALUE;
				mbic = Double.MAX_VALUE;
				aic = Double.MAX_VALUE;
			}

			try{
				p = LinearModelUtils.Ftest(F, termssdf[1], errorssdf[1]);  
			} catch (Exception e){
				p = Double.NaN;
			}


			switch(modelType){
			case pvalue:
				try {
					if (p > maxp) {
						maxterm = t;
						maxp = p;
						minF = F; 
					}
				} catch(Exception e){

				}
				break;
			case bic:
				if(bic < bestbic){
					bestbic = bic;
					bestaic = aic;
					bestmbic = mbic;
					maxterm = t;  //Actually, this should be "minterm" becasue we are finding the model with the minimum BIC, but I am keeping this as "maxterm" for simpilicity of the calculations
					worstMarkerEffect = markerEffect;
					maxp = p;
					minF = F; 
				}
				break;
			case mbic:
				if(mbic < bestmbic){
					bestbic = bic;
					bestaic = aic;
					bestmbic = mbic;
					maxterm = t;
					worstMarkerEffect = markerEffect;
					maxp = p;
					minF = F;
				}
				break;
			case aic:
				if(aic < bestaic){
					bestbic = bic;
					bestaic = aic;
					bestmbic = mbic;
					maxterm = t;
					worstMarkerEffect = markerEffect;
					maxp = p;
					minF = F;
				}
				break;
			} 

		}


		//if that maxp is >= exitlimit, then remove maxterm from the model, recalculate the model, and return true;
		switch(modelType){
		case pvalue:
			test = maxp >= exitlimit;
			break;
		case bic:
			test = bestbic < globalbestbic;
			break;
		case mbic:
			test = bestmbic < globalbestmbic;
			break;
		case aic:
			test = bestaic < globalbestaic;
			break;
		}


		if ((test)&&(maxterm != 0)) {
			ModelEffect me = currentModel.remove(maxterm);
			globalbestbic = bestbic;
			globalbestmbic = bestmbic;
			globalbestaic = bestaic;
			return true;
		}

		return false;
	}

	public void runPermutationTestNoMissingData() {
		ArrayList<double[]> permutedData = new ArrayList<double[]>();
		DoubleMatrix PvalueVectorAcrossMarkers = null;

		int indexOfThreshold = (int) (alpha*numberOfPermutations);

		int numberOfSites = myGenotype.numberOfSites();
		DoubleMatrix[][] Xmatrices;
		System.out.println("-----------------Running permutations...----------------");

		int numberOfObs = y.length;
		double totalSS = 0;
		for (int i = 0; i < numberOfObs; i++) totalSS += y[i] * y[i];

		ArrayList<double[]> baseModelSSdf = new ArrayList<double[]>();
		double[] permTotalSS = new double[numberOfPermutations];

		int[] mean = new int[numberOfObs];
		FactorModelEffect meanME = new FactorModelEffect(mean, false, "mean");
		int columnNumber = 2;
		columnNumber += covariateAttributeList.size();
		columnNumber += factorAttributeList.size();

		Xmatrices = new DoubleMatrix[1][columnNumber];//Intercept+family term+                

		Xmatrices[0][0] = meanME.getX();
		int effectCount = 1;
		
		for (PhenotypeAttribute attr : factorAttributeList) {
			CategoricalAttribute ca = (CategoricalAttribute) attr;
			ArrayList<String> ids = new ArrayList<String>();
			String[] labels = AssociationUtils.getNonMissingValues(ca.allLabels(), missing);
			int[] levels = ModelEffectUtils.getIntegerLevels(labels, ids);
			FactorModelEffect fme = new FactorModelEffect(levels, true, new Object[]{attr.name(), ids});
			if (isNested && myPhenotype.indexOfAttribute(attr) == nestingFactorIndex) {
				nestingEffect = fme;
				nestingFactorNames = ids; 
			}
			Xmatrices[0][effectCount++] = fme.getX();
		}

		for (PhenotypeAttribute attr : covariateAttributeList) {
			NumericAttribute na = (NumericAttribute) attr;
			double[] cov = AssociationUtils.getNonMissingDoubles(na.floatValues(), missing);
			CovariateModelEffect cme = new CovariateModelEffect(cov, na.name());
			Xmatrices[0][effectCount++] = cme.getX();
		}

		SweepFastLinearModel baseSflm = new SweepFastLinearModel(currentModel, y);


		DoubleMatrix resAsDoubleMatrix = baseSflm.getResiduals();
		DoubleMatrix predAsDoubleMatrix = baseSflm.getPredictedValues();
		double[] res = new double[numberOfObs];
		double[] pred = new double[numberOfObs];
		for(int i = 0; i < numberOfObs; i++){
			res[i] = resAsDoubleMatrix.get(i,0);
			pred[i] = predAsDoubleMatrix.get(i,0);
		}

		//permute data
		double[] minP = new double[numberOfPermutations];

		DoubleMatrix minPAsDoubleMatrix = null;
		for (int p = 0; p < numberOfPermutations; p++) {
			minP[p] = 1;
			double[] pdata = new double[numberOfObs];
			System.arraycopy(res, 0, pdata, 0, numberOfObs);
			BasicShuffler.shuffle(pdata);

			for(int i = 0; i < numberOfObs; i++){
				pdata[i] = pdata[i] + pred[i];
			} //Add the predicted value back to the residual
			totalSS = 0;
			for (int i = 0; i < numberOfObs; i++) totalSS += pdata[i] * pdata[i];
			permTotalSS[p] = totalSS;
			permutedData.add(pdata); 

			SweepFastLinearModel sflm = new SweepFastLinearModel(currentModel, pdata);
			baseModelSSdf.add(sflm.getFullModelSSdf());
		}


		for (int m = 0; m < numberOfSites; m++) {

			//create the appropriate marker effect
			ModelEffect markerEffect = null;
			SNP snp = new SNP(myGenotype.siteName(m), myGenotype.chromosome(m), myGenotype.chromosomalPosition(m), m);

            //if the Genotype has a reference allele use that, else if genotype use that else nothing, allele probability not implemented
            //for genotype use an additive model for now
            if (myGenotype.hasReferenceProbablity()) {
            	double[] cov = new double[numberNotMissing];
            	int count = 0;
            	float[] referenceProb = myData.referenceProb(m);
            	for (int i = 0; i < totalNumber; i++) {
            		if (!missing.fastGet(i)) cov[count++] = referenceProb[i];
            	}

            	if (isNested) {
            		CovariateModelEffect cme = new CovariateModelEffect(cov);
            		markerEffect = new NestedCovariateModelEffect(cme, nestingEffect);
            		markerEffect.setID(snp);
            	} else {
            		markerEffect = new CovariateModelEffect(cov, snp);
            	}

            } else if (myGenotype.hasGenotype()) {
            	byte minor = myGenotype.minorAllele(m);
            	double[] cov = new double[numberNotMissing];
            	byte[] siteGeno = AssociationUtils.getNonMissingBytes(myData.genotypeAllTaxa(m), missing);
            	for (int i = 0; i < numberNotMissing; i++) {
            		byte[] diploidAlleles = GenotypeTableUtils.getDiploidValues(siteGeno[i]);
            		if (diploidAlleles[0] == minor) cov[i] += 0.5;
            		if (diploidAlleles[1] == minor) cov[i] += 0.5;
            	}
            	
            	if (isNested) {
            		CovariateModelEffect cme = new CovariateModelEffect(cov);
            		markerEffect = new NestedCovariateModelEffect(cme, nestingEffect);
            		markerEffect.setID(snp);
            	} else {
            		markerEffect = new CovariateModelEffect(cov, snp);
            	}
           	
            } else {
            	throw new IllegalArgumentException("No genotypes or reference probabilities in the data set.");
            }
			
			Xmatrices[0][columnNumber-1] = markerEffect.getX(); //columnNumber-1 because java
			//starts counting things from 0
			DoubleMatrix X = DoubleMatrixFactory.DEFAULT.compose(Xmatrices);

			currentModel.add(markerEffect);
			SweepFastLinearModel sflm = new SweepFastLinearModel(currentModel, y);         
			double[] modelSSdf = sflm.getFullModelSSdf();
			DoubleMatrix G = sflm.getInverseOfXtX();

			for (int p = 0; p < numberOfPermutations; p++) {
				double[] pdata = permutedData.get(p);
				DoubleMatrix yperm = DoubleMatrixFactory.DEFAULT.make(numberOfObs, 1, pdata);
				totalSS = permTotalSS[p];

				DoubleMatrix Xty = X.crossproduct(yperm);
				double[] reducedSSdf = baseModelSSdf.get(p);
				double fullSS = Xty.crossproduct(G.mult(Xty)).get(0, 0);
				double fulldf = modelSSdf[1];
				double markerSS = fullSS - reducedSSdf[0];
				double markerdf = fulldf - reducedSSdf[1];
				double errorSS = totalSS - fullSS;
				double errordf = numberOfObs - fulldf;
				double F = markerSS / markerdf / errorSS * errordf;
				double probF;
				try {
					probF = LinearModelUtils.Ftest(F, markerdf, errordf);
				} catch(Exception e) {
					probF = 1;
				}
				minP[p] = probF;
				minPAsDoubleMatrix = DoubleMatrixFactory.DEFAULT.make(minP.length, 1, minP);
			}
			if(m == 0){
				PvalueVectorAcrossMarkers = minPAsDoubleMatrix;
			}else{
				PvalueVectorAcrossMarkers = PvalueVectorAcrossMarkers.concatenate(minPAsDoubleMatrix, false);//Change this to a DoubleMatrix; transpose it
				//for the next step                            
			}

			currentModel.remove(markerEffect);
		}

		System.out.println(PvalueVectorAcrossMarkers.toString());
		for(int p=0; p<numberOfPermutations; p++){
			double minFromOnePermutation = 1;
			for(int m=0; m < numberOfSites; m++){
				//System.out.println("Marker:" + m + ":Permutation:"+ p +":"+ PvalueVectorAcrossMarkers.get(p,m));
				minFromOnePermutation = Math.min(PvalueVectorAcrossMarkers.get(p,m), minFromOnePermutation);
			}
			//obtain the minimum P-value across each permutation
			//System.out.println("--------------p is " + p);
			minPvalues[p] = minFromOnePermutation;
			//System.out.println(" minPvalues[p] " +  minPvalues[p]);
			//PvalueVectorAcrossMarkers.
		}

		//sort the P-values from smallest to largest
		Arrays.sort(minPvalues);

		//select the (alpha*nperm)th smallest P-value
		enterlimit = minPvalues[indexOfThreshold];
		exitlimit = 2*enterlimit;

		//Turn this into a new method. Put it into the plugger. Have a logical statement to run this 
		//only if permutations are selected.
		//Create a table that has the permuted P-values

		System.out.println("--------------the Enter limit is " + enterlimit);
		System.out.println("--------------the Exit limit is " + exitlimit); 
		System.out.println("--------------indexOfThreshold is " + indexOfThreshold); 
	}

 	public TableReport getPermutationReport() {
		if (numberOfPermutations == 0) return null;
                
                LinkedList<Object[]> pValueReportTable = new LinkedList<Object[]>();;
                //Object[] pValueReportRow = new Object[numberOfPermutations];
                for(int i=0; i < numberOfPermutations; i++){
                    Object[] pValueReportRow = new Object[1];
                    pValueReportRow[0] = new Double(minPvalues[i]);          
                    //System.out.println("minPvalues[i]" + minPvalues[i]);
                    pValueReportTable.add(pValueReportRow);
                }
                
                //pValueReportTable.add(pValueReportRow);
               	Object[][] table = new Object[pValueReportTable.size()][];
		pValueReportTable.toArray(table);                
		//columns are trait, snp name, locus, position, factor value, estimate
		String reportName = numberOfPermutations + " Permutations for "  + datasetName;
		String[] reportHeader = new String[]{"P-value"};

		return new SimpleTableReport(reportName, reportHeader, table);
	}       
        
        
	public LinkedList<Object[]> createReportRowsFromCurrentModel(SweepFastLinearModel sflm) {
		String traitname = dataAttributeList.get(currentPhenotypeIndex).name();
		int ncol = anovaReportHeader.length;
		LinkedList<Object[]> reportTable = new LinkedList<Object[]>();
		double[] residualSSdf = sflm.getResidualSSdf();
		int n = y.length;
		int numberOfSites = myGenotype.numberOfSites();


		//Calcualte the BIC
		double [] errorss = sflm.getResidualSSdf();
		double [] modeldf = sflm.getFullModelSSdf();
		double pForBIC  = modeldf[1];
		double bic = (n*Math.log(errorss[0]/n)) + (pForBIC*Math.log(n)) ; 
		double aic = (n*Math.log(errorss[0]/n)) + (2*pForBIC) ;

		System.out.println("error ss is " + errorss[0]);
		System.out.println("pForBIC is " + pForBIC);
		System.out.println("n is " + n);


		//Calculate the mBIC
		double numberOfTwoWayInteractions = (double) (numberOfSites*(numberOfSites-1))/2;			
		double pForMbic = modeldf[1];
		double qForMbic = 0; //This is going to be the number of two-way interaction terms. For now, this is not implemented, and thus I am setting it equal to zer
		double mbic;
		if(qForMbic == 0){
			mbic = (n*Math.log(errorss[0])) + ((pForMbic+qForMbic)*Math.log(n)) + (2*pForMbic*(Math.log((numberOfSites/2.2)-1)));
		}else{
			mbic = (n*Math.log(errorss[0])) + ((pForMbic+qForMbic)*Math.log(n)) + (2*pForMbic*(Math.log((numberOfSites/2.2)-1))) + 
					(2*qForMbic*(Math.log((numberOfTwoWayInteractions/2.2)-1) )); 
		}
                
		int effectPtr = 0;
		for (ModelEffect me : currentModel) {
			Object[] reportRow = new Object[ncol];
			int ptr = 0;
			reportRow[ptr++] = traitname;

			if (me.getID() instanceof SNP) {
				SNP snp = (SNP) me.getID();
				reportRow[ptr++] = snp.name;
				reportRow[ptr++] = snp.locus.getName();
				reportRow[ptr++] = Integer.toString(snp.position);
			} else {
				if (me.getID() instanceof String) reportRow[ptr++] = me.getID().toString();
				else if (me instanceof FactorModelEffect) reportRow[ptr++] = ((Object[]) me.getID())[0].toString();
				else reportRow[ptr++] = me.getID().toString();
				reportRow[ptr++] = "--";
				reportRow[ptr++] = "--";
			}
			double[] effectSSdf = sflm.getMarginalSSdf(effectPtr);
			double ms = effectSSdf[0] / effectSSdf[1];
			double Fval = ms / residualSSdf[0] * residualSSdf[1];
			double pval;
			try {
				pval = LinearModelUtils.Ftest(Fval, effectSSdf[1], residualSSdf[1]);
			} catch (Exception e) {
				pval = Double.NaN;
			}
			reportRow[ptr++] = new Integer((int) effectSSdf[1]);
			reportRow[ptr++] = new Double(effectSSdf[0]);
			reportRow[ptr++] = new Double(ms);
			reportRow[ptr++] = new Double(Fval);
			reportRow[ptr++] = new Double(pval);
			reportRow[ptr++] = new Double(bic);
			reportRow[ptr++] = new Double(mbic);
			reportRow[ptr++] = new Double(aic);
			double[] modelSSdf = sflm.getModelcfmSSdf();
			reportRow[ptr++] = new Double(modelSSdf[0]/(modelSSdf[0] + residualSSdf[0]));
			reportTable.add(reportRow);
			effectPtr++;
		}
		int ptr = 0;
		Object[] reportRow = new Object[ncol];
		reportRow[ptr++] = traitname;
		reportRow[ptr++] = "Error";
		reportRow[ptr++] = "--";
		reportRow[ptr++] = "--";
		reportRow[ptr++] = new Integer((int) residualSSdf[1]);
		reportRow[ptr++] = new Double(residualSSdf[0]);
		reportRow[ptr++] = new Double(residualSSdf[0]/residualSSdf[1]);
		reportRow[ptr++] = new Double(Double.NaN);
		reportRow[ptr++] = new Double(Double.NaN);
		reportTable.add(reportRow);

		return reportTable;
	}

	public LinkedList<Object[]> createReportRowsFromCurrentModelAfterScanCI(SweepFastLinearModel sflm) {
		String traitname = dataAttributeList.get(currentPhenotypeIndex).name();
		int ncol = anovaReportWithCIHeader.length;
		LinkedList<Object[]> reportTable = new LinkedList<Object[]>();
		double[] residualSSdf = sflm.getResidualSSdf();
		int n = y.length;
		int numberOfSites = myGenotype.numberOfSites();
		int firstMarkerIndex = 1;

		firstMarkerIndex += factorAttributeList.size();
		firstMarkerIndex += covariateAttributeList.size();
		
		//Calcualte the BIC
		double [] errorss = sflm.getResidualSSdf();
		double [] modeldf = sflm.getFullModelSSdf();
		double pForBIC  = modeldf[1];
		double bic = (n*Math.log(errorss[0]/n)) + (pForBIC*Math.log(n)) ; 
		double aic = (n*Math.log(errorss[0]/n)) + (2*pForBIC) ;

		System.out.println("error ss is " + errorss[0]);
		System.out.println("pForBIC is " + pForBIC);
		System.out.println("n is " + n);


		//Calculate the mBIC
		double numberOfTwoWayInteractions = (double) (numberOfSites*(numberOfSites-1))/2;			
		double pForMbic = modeldf[1];
		double qForMbic = 0; //This is going to be the number of two-way interaction terms. For now, this is not implemented, and thus I am setting it equal to zer
		double mbic;
		if(qForMbic == 0){
			mbic = (n*Math.log(errorss[0])) + ((pForMbic+qForMbic)*Math.log(n)) + (2*pForMbic*(Math.log((numberOfSites/2.2)-1)));
		}else{
			mbic = (n*Math.log(errorss[0])) + ((pForMbic+qForMbic)*Math.log(n)) + (2*pForMbic*(Math.log((numberOfSites/2.2)-1))) + 
					(2*qForMbic*(Math.log((numberOfTwoWayInteractions/2.2)-1) )); 
		}

		int effectPtr = 0;
		for (ModelEffect me : currentModel) {
			Object[] reportRow = new Object[ncol];
			int ptr = 0;
			reportRow[ptr++] = traitname;

			if (me.getID() instanceof SNP) {
				SNP snp = (SNP) me.getID();
				reportRow[ptr++] = snp.name;
				reportRow[ptr++] = snp.locus.getName();
				reportRow[ptr++] = Integer.toString(snp.position);
			} else {
				if (me.getID() instanceof String) reportRow[ptr++] = me.getID().toString();
				else if (me instanceof FactorModelEffect) reportRow[ptr++] = ((Object[]) me.getID())[0].toString();
				else reportRow[ptr++] = me.getID().toString();
				reportRow[ptr++] = "--";
				reportRow[ptr++] = "--";
			}
			double[] effectSSdf = sflm.getMarginalSSdf(effectPtr);
			double ms = effectSSdf[0] / effectSSdf[1];
			double Fval = ms / residualSSdf[0] * residualSSdf[1];
			double pval;
			try {
				pval = LinearModelUtils.Ftest(Fval, effectSSdf[1], residualSSdf[1]);
			} catch (Exception e) {
				pval = Double.NaN;
			}
			reportRow[ptr++] = new Integer((int) effectSSdf[1]);
			reportRow[ptr++] = new Double(effectSSdf[0]);
			reportRow[ptr++] = new Double(ms);
			reportRow[ptr++] = new Double(Fval);
			reportRow[ptr++] = new Double(pval);
			reportRow[ptr++] = new Double(bic);
			reportRow[ptr++] = new Double(mbic);
			reportRow[ptr++] = new Double(aic);
			double[] modelSSdf = sflm.getModelcfmSSdf();
			reportRow[ptr++] = new Double(modelSSdf[0]/(modelSSdf[0] + residualSSdf[0]));

			//System.out.println("The value of effectPtr is "+effectPtr);

			if (me.getID() instanceof SNP) {
				reportRow[ptr++] = new String(myGenotype.siteName(theUpperAndLowerBound[effectPtr-firstMarkerIndex][0]));
				reportRow[ptr++] = new String(myGenotype.siteName(theUpperAndLowerBound[effectPtr-firstMarkerIndex][1]));

			} else {
				reportRow[ptr++] = "--";
				reportRow[ptr++] = "--";
			}
			reportTable.add(reportRow);
			effectPtr++;
		}
		int ptr = 0;
		Object[] reportRow = new Object[ncol];
		reportRow[ptr++] = traitname;
		reportRow[ptr++] = "Error";
		reportRow[ptr++] = "--";
		reportRow[ptr++] = "--";
		reportRow[ptr++] = new Integer((int) residualSSdf[1]);
		reportRow[ptr++] = new Double(residualSSdf[0]);
		reportRow[ptr++] = new Double(residualSSdf[0]/residualSSdf[1]);
		reportRow[ptr++] = new Double(Double.NaN);
		reportRow[ptr++] = new Double(Double.NaN);
		reportTable.add(reportRow);

		return reportTable;
	}

	public Datum createReportFromCurrentModel(SweepFastLinearModel sflm) {
		String traitname = dataAttributeList.get(currentPhenotypeIndex).name();
		LinkedList<Object[]> reportTable = createReportRowsFromCurrentModel(sflm);
		Object[][] table = new Object[reportTable.size()][];
		reportTable.toArray(table);
		String reportName = "ANOVA for " + traitname + ", " + datasetName;
		TableReport tr = new SimpleTableReport(reportName, anovaReportHeader, table);
		return new Datum(reportName, tr, "");
	}

	public Datum createReportFromCurrentModelWithCI(SweepFastLinearModel sflm) {
		String traitname = dataAttributeList.get(currentPhenotypeIndex).name();
		LinkedList<Object[]> reportTable = createReportRowsFromCurrentModelAfterScanCI(sflm);
		Object[][] table = new Object[reportTable.size()][];
		reportTable.toArray(table);
		String reportName = "ANOVA With CI for " + traitname + ", " + datasetName;
		TableReport tr = new SimpleTableReport(reportName, anovaReportWithCIHeader, table);
		return new Datum(reportName, tr, "");
	}
        
	public void appendAnovaResults(SweepFastLinearModel sflm) {
		resultRowsAnova.addAll(createReportRowsFromCurrentModel(sflm));
	}

        public void appendAnovaResultsWithCI(SweepFastLinearModel sflm) {
            resultRowsAnovaWithCI.addAll(createReportRowsFromCurrentModelAfterScanCI(sflm));
	}

	public TableReport getAnovaReport() {
		String reportName = "ANOVA table for " + datasetName;
		Object[][] table = new Object[resultRowsAnova.size()][];
		resultRowsAnova.toArray(table);
		return new SimpleTableReport(reportName, anovaReportHeader, table);
	}

	public TableReport getAnovaReportWithCI() {
		String reportName = "ANOVA table with CI scan for " + datasetName;
		Object[][] table = new Object[resultRowsAnovaWithCI.size()][];
		resultRowsAnovaWithCI.toArray(table);
		return new SimpleTableReport(reportName, anovaReportWithCIHeader, table);
	}
        
	public TableReport getMarkerEffectReport() {
		if (rowsSiteEffectTable.size() == 0) return null;
		//columns are trait, snp name, locus, position, factor value, estimate
		String reportName = "Marker effects for " + datasetName;
		String[] reportHeader = new String[]{"Trait","Snp","Locus","Position","Within","Estimate"};
		Object[][] table = new Object[rowsSiteEffectTable.size()][];
		rowsSiteEffectTable.toArray(table);
		return new SimpleTableReport(reportName, reportHeader, table);
	}

	public TableReport getMarkerEffectReportWithCI() {
		if (rowsSiteEffectTableWithCI.size() == 0) return null;
		//columns are trait, snp name, locus, position, factor value, estimate
		String reportName = "Marker effects for " + datasetName;
		String[] reportHeader = new String[]{"Trait","Snp","Locus","Position","Within","Estimate"};
		Object[][] table = new Object[rowsSiteEffectTableWithCI.size()][];
		rowsSiteEffectTableWithCI.toArray(table);
		return new SimpleTableReport(reportName, reportHeader, table);
	}

	public void appendSiteEffectEstimates(SweepFastLinearModel sflm) {

		//columns are trait, snp name, locus, position, factor value, estimate
		int nBaseModelParameters = 0;
		String traitname = dataAttributeList.get(currentPhenotypeIndex).name();
		for (int i = 0; i < numberOfBaseModelEffects; i++) {
			nBaseModelParameters += currentModel.get(i).getEffectSize();
		}

		double[] beta = sflm.getBeta();
		int parameterIndex = nBaseModelParameters;
		for (int s = numberOfBaseModelEffects; s < currentModel.size(); s++) {
			ModelEffect me = currentModel.get(s);
			if (me instanceof NestedCovariateModelEffect) {
				NestedCovariateModelEffect ncme = (NestedCovariateModelEffect) me;
				FactorModelEffect fme = ncme.getFactorModelEffect();
				SNP snp = (SNP) ncme.getID();
				int n = nestingFactorNames.size();
				for (int i = 0; i < n; i++) {
					Object[] rowValues = new Object[6];
					int ptr = 0;
					rowValues[ptr++] = traitname;
					rowValues[ptr++] = snp.name;
					rowValues[ptr++] = snp.locus.getName();
					rowValues[ptr++] = new Integer(snp.position);
					rowValues[ptr++] = nestingFactorNames.get(i);
					rowValues[ptr++] = new Double(beta[parameterIndex++]);
					rowsSiteEffectTable.add(rowValues);
				}
			} else if (me instanceof CovariateModelEffect) {
				SNP snp = (SNP) me.getID();
				Object[] rowValues = new Object[6];
				int ptr = 0;
				rowValues[ptr++] = traitname;
				rowValues[ptr++] = snp.name;
				rowValues[ptr++] = snp.locus.getName();
				rowValues[ptr++] = new Integer(snp.position);
				rowValues[ptr++] = "NA";
				rowValues[ptr++] = new Double(beta[parameterIndex++]);
				rowsSiteEffectTable.add(rowValues);
			} else if (me instanceof FactorModelEffect) {

			}
		}
	}

        public void appendSiteEffectEstimatesWithCI(SweepFastLinearModel sflm) {

		//columns are trait, snp name, locus, position, factor value, estimate
		int nBaseModelParameters = 0;
		String traitname = dataAttributeList.get(currentPhenotypeIndex).name();
		for (int i = 0; i < numberOfBaseModelEffects; i++) {
			nBaseModelParameters += currentModel.get(i).getEffectSize();
		}

		double[] beta = sflm.getBeta();
		int parameterIndex = nBaseModelParameters;
		for (int s = numberOfBaseModelEffects; s < currentModel.size(); s++) {
			ModelEffect me = currentModel.get(s);
			if (me instanceof NestedCovariateModelEffect) {
				NestedCovariateModelEffect ncme = (NestedCovariateModelEffect) me;
				FactorModelEffect fme = ncme.getFactorModelEffect();
				SNP snp = (SNP) ncme.getID();
				int n = nestingFactorNames.size();
				for (int i = 0; i < n; i++) {
					Object[] rowValues = new Object[6];
					int ptr = 0;
					rowValues[ptr++] = traitname;
					rowValues[ptr++] = snp.name;
					rowValues[ptr++] = snp.locus.getName();
					rowValues[ptr++] = new Integer(snp.position);
					rowValues[ptr++] = nestingFactorNames.get(i);
					rowValues[ptr++] = new Double(beta[parameterIndex++]);
					rowsSiteEffectTableWithCI.add(rowValues);
				}
			} else if (me instanceof CovariateModelEffect) {
				SNP snp = (SNP) me.getID();
				Object[] rowValues = new Object[6];
				int ptr = 0;
				rowValues[ptr++] = traitname;
				rowValues[ptr++] = snp.name;
				rowValues[ptr++] = snp.locus.getName();
				rowValues[ptr++] = new Integer(snp.position);
				rowValues[ptr++] = "NA";
				rowValues[ptr++] = new Double(beta[parameterIndex++]);
				rowsSiteEffectTableWithCI.add(rowValues);
			} else if (me instanceof FactorModelEffect) {

			}
		}
	}

        public Phenotype generateChromosomeResidualsFromCurrentModel() {
            List<PhenotypeAttribute> attributes = new ArrayList<>();
            List<ATTRIBUTE_TYPE> types = new ArrayList<>();

            Taxon[] allTaxa = myPhenotype.taxaAttribute().allTaxa();
            Taxon[] nonmissingTaxa = AssociationUtils.getNonMissingValues(allTaxa, missing);
            attributes.add(new TaxaAttribute(Arrays.asList(nonmissingTaxa)));
            types.add(ATTRIBUTE_TYPE.taxa);
            
            //this next step will include family in the return Phenotype
            //any covariates will not be included
            for (PhenotypeAttribute factor : factorAttributeList) {
                String[] values = ((CategoricalAttribute) factor).allLabels();
                attributes.add(new CategoricalAttribute(factor.name(), AssociationUtils.getNonMissingValues(values, missing)));
                types.add(ATTRIBUTE_TYPE.factor);
            }
            
            for (int c = 0; c <10; c++) {
                String chrname = String.format("c%d",c + 1);
                
                //create a model without this chromosome
                List<ModelEffect> chrModel = currentModel.stream()
                        .filter(notInChr(c + 1))
                        .collect(Collectors.toList());
                
                SweepFastLinearModel sflm = new SweepFastLinearModel(chrModel, y);
                
                //add the residuals to the Phenotype
                DoubleMatrix resid = sflm.getResiduals();
                float[] data = AssociationUtils.convertDoubleArrayToFloat(resid.to1DArray());
                attributes.add(new NumericAttribute(chrname, data, new OpenBitSet(data.length)));
                types.add(ATTRIBUTE_TYPE.data);
            }
            
            return new PhenotypeBuilder().fromAttributeList(attributes, types).build().get(0);
        }
        
        private Predicate<ModelEffect> notInChr(int chr) {
            return me -> (!(me.getID() instanceof SNP)) || ((me.getID() instanceof SNP) && ( ((SNP) me.getID()).locus.getChromosomeNumber() != chr ));
        }
        
	public void setEnterlimits(double[] enterlimits) {
		this.enterlimits = enterlimits;
	}

	public void setExitlimits(double[] exitlimits) {
		this.exitlimits = exitlimits;
	}

	public void setEnterlimit(double enterlimit) {
		this.enterlimit = enterlimit;
	}

	public void setExitlimit(double exitlimit) {
		this.exitlimit = exitlimit;
	}

	public void setMaxNumberOfMarkers(int maxNumberOfMarkers) {
		this.maxNumberOfMarkers = maxNumberOfMarkers;
	}

	public void setNestingEffectIndex(int nestingFactorIndex) {
		this.nestingFactorIndex = nestingFactorIndex;
	}

	public void setNested(boolean isNested) {
		this.isNested = isNested;
	}
        
	public void setNumberOfPermutations(int numberOfPermutations) {
		this.numberOfPermutations = numberOfPermutations;
                minPvalues = new double[this.numberOfPermutations];
}
 
	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}
}
