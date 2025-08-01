package net.maizegenetics.analysis.association;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.score.SiteScore;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.LinearModelUtils;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffectUtils;
import net.maizegenetics.stats.linearmodels.SolveByOrthogonalizing;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportBuilder;

public abstract class AbstractFixedEffectLM implements FixedEffectLM {
	protected static Logger myLogger = LogManager.getLogger(AbstractFixedEffectLM.class);
	
	protected final Datum myDatum;
	protected final GenotypePhenotype myGenoPheno;
	protected final int numberOfObservations;
	protected final int numberOfSites;
	protected final List<PhenotypeAttribute> myDataAttributes;
	protected final List<PhenotypeAttribute> myFactorAttributes;
	protected final List<PhenotypeAttribute> myCovariateAttributes;
	protected TableReportBuilder siteReportBuilder;
	protected TableReportBuilder alleleReportBuilder;
	protected int numberOfSiteReportColumns;
	protected int numberOfAlleleReportColumns;
	protected float[] allData;
	protected int myCurrentSite;
	protected int myCurrentSiteMinimumClassSize;
	protected double[] siteData;
	protected OpenBitSet missingObsForSite;
	protected String currentTraitName;
	protected boolean areTaxaReplicated;
	protected boolean saveToFile = false;
	protected String siteReportFilename;
	protected String alleleReportFilename;
	protected double maxP = 1.0;
	protected FixedEffectLMPlugin myParentPlugin;
	protected boolean appendAddDomEffects = false;
	
	//filtering criteria
	protected int minClassSize = 0;
	protected boolean biallelicOnly = false;
	protected boolean outputSiteStats = false;
	protected String siteStatsFile = null;
	
	//fields used for permutation testing
	protected boolean permute = false;
	protected int numberOfPermutations = 0;
	protected double[] minP = null;
	protected List<DoubleMatrix> permutedData;
	protected double[] baseErrorSSdf;
	protected double[] totalcfmSSdf;
	protected double[] markerSSdf;
	protected double[] errorSSdf;
	protected List<Object[]> siteTableReportRows;
	protected int markerpvalueColumn;
	protected int permpvalueColumn;
	protected ArrayList<ModelEffect> myModel;
	protected DoubleMatrix G;
	protected ArrayList<ModelEffect> myBaseModel;
	protected int numberOfBaseEffects;
	protected int taxaEffectNumber;
	protected int randomSeed;
	protected boolean useRandomSeed = false;
	protected Random rand = null;


    protected static final Map<SiteScore.SITE_SCORE_TYPE, String> typeNameMap;
    static {
    	typeNameMap = new HashMap<SiteScore.SITE_SCORE_TYPE, String>();
    	typeNameMap.put(SiteScore.SITE_SCORE_TYPE.ProbA, "A");
    	typeNameMap.put(SiteScore.SITE_SCORE_TYPE.ProbC, "C");
    	typeNameMap.put(SiteScore.SITE_SCORE_TYPE.ProbG, "G");
    	typeNameMap.put(SiteScore.SITE_SCORE_TYPE.ProbT, "T");
    	typeNameMap.put(SiteScore.SITE_SCORE_TYPE.ProbGap, "-");
    	typeNameMap.put(SiteScore.SITE_SCORE_TYPE.ProbInsertion, "+");
    }
    
	public AbstractFixedEffectLM(Datum dataset, FixedEffectLMPlugin parentPlugin) {
		myDatum = dataset;
		myParentPlugin = parentPlugin;
		
		myGenoPheno = (GenotypePhenotype) myDatum.getData();
		numberOfObservations = myGenoPheno.phenotype().numberOfObservations();
		numberOfSites = myGenoPheno.genotypeTable().numberOfSites();
		
		myDataAttributes = myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.data);
		myFactorAttributes = myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.factor);
		myCovariateAttributes = myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.covariate);
		siteTableReportRows = new ArrayList<Object[]>();
		testTaxaReplication();
	}

	@Override
	public void initializeReportBuilders() {
		String tableName = "GLM Statistics - " + myDatum.getName();
		String[] columnNames = siteReportColumnNames();
		numberOfSiteReportColumns = columnNames.length;
		if (saveToFile) siteReportBuilder = TableReportBuilder.getInstance(tableName, columnNames, siteReportFilename);
		else siteReportBuilder = TableReportBuilder.getInstance(tableName, columnNames);
		 
		tableName = "GLM Genotype Effects - " + myDatum.getName();
		columnNames = alleleReportColumnNames();
		numberOfAlleleReportColumns = columnNames.length;
		if (saveToFile) alleleReportBuilder = TableReportBuilder.getInstance(tableName, columnNames, alleleReportFilename);
		else alleleReportBuilder = TableReportBuilder.getInstance(tableName, columnNames);
	}
	
	@Override
	public void solve() {
		//loop through data attributes
		//	loop through sites
		initializeReportBuilders();
		int numberOfAttributes = myDataAttributes.size();
		int numberOfTestsTotal = numberOfAttributes * numberOfSites;
		int numberOfTestsCalculated = 0;
		int updateInterval = Math.max(1, numberOfTestsTotal / 100);
		
//		long start = System.currentTimeMillis();
		for (PhenotypeAttribute dataAttribute:myDataAttributes) {
			currentTraitName = dataAttribute.name();
			OpenBitSet missingObs = new OpenBitSet(dataAttribute.missing());
			for (PhenotypeAttribute attr:myFactorAttributes) missingObs.or(attr.missing());
			for (PhenotypeAttribute attr:myCovariateAttributes) missingObs.or(attr.missing());
			allData = (float[]) dataAttribute.allValues();
			if (permute) {
				missingObsForSite = missingObs;
				createPermutedData();
			}
			for (int s = 0; s < numberOfSites; s++) {
				//updata missing obs for this site
				myCurrentSite = s;
				getGenotypeAndUpdateMissing(missingObs);
				boolean keepSite = applySiteFilters();
				if (!keepSite) continue;
				siteData = AssociationUtils.getNonMissingDoubles(allData, missingObsForSite);
				myBaseModel = baseModel();
				numberOfBaseEffects = myBaseModel.size();
				analyzeSite();
				if (permute) updateMinP(missingObs);
				
				numberOfTestsCalculated++;
				if (numberOfTestsCalculated % updateInterval == 0) {
				        double percentTested = 100.0 * ((double) numberOfTestsCalculated) / numberOfTestsTotal;
				        percentTested = Math.min(percentTested, 100);
					if (myParentPlugin != null) myParentPlugin.updateProgress((int) percentTested);
				}
			}
//			System.out.printf("Sites analyzed in %d ms\n", System.currentTimeMillis() - start);
			if (permute) updateReportsWithPermutationP();
		}
		if (saveToFile) {
			siteReportBuilder.build();
			alleleReportBuilder.build();
		}
	}
	
	private boolean applySiteFilters() {
		//does the site pass the filter for biallelic sites
		//start with the sites to be analyzed
		if (!myGenoPheno.genotypeTable().hasGenotype()) return true;
		byte[] siteGeno = myGenoPheno.genotypeAllTaxa(myCurrentSite);
		int nsites = siteGeno.length;
		Map<Byte,Integer> genoCountMap = new HashMap<>();
		for (int s = 0; s < nsites; s++) {
			if (!missingObsForSite.get(s)) {
				Integer genoCount = genoCountMap.get(siteGeno[s]);
				if (genoCount == null) {
					genoCountMap.put(siteGeno[s], 1);
				} else {
					genoCountMap.put(siteGeno[s], genoCount + 1);
				}
			}
		}
		
		boolean keepSite = true;
		if (biallelicOnly) {
			keepSite = false;
			//the site is biallelic if genoCount = 2 or if genoCount == 3 and one of the genotypes is heterozygous
			if (genoCountMap.size() == 2) keepSite = true;
			else if (genoCountMap.size() == 3) {
				int hetCount = 0;
				for (Byte genoval : genoCountMap.keySet()) {
					if (GenotypeTableUtils.isHeterozygous(genoval)) hetCount++;
				}
				if (hetCount == 1) keepSite = true;
			}
		}
		
		//apply minClassSizeFilter
		if (keepSite && minClassSize > 0) {
			int numberBigEnough = 0;
			int numberTooSmall = 0;
			for (Integer ival : genoCountMap.values()) {
				if (ival < minClassSize) numberTooSmall++;
				else numberBigEnough++;
			}
			
			//if there is only one class that has enough taxa, eliminate the site
			if (numberBigEnough < 2) keepSite = false;
			//if the minimum class size is too small and there are more than two classes set that class to missing
			else if (numberTooSmall > 0) {
				for (Byte Bval : genoCountMap.keySet()) {
					int classSize = genoCountMap.get(Bval);
					if (classSize < minClassSize) {
						byte classValue = Bval.byteValue();
						for (int s = 0; s < nsites; s++) {
							if (siteGeno[s] == classValue) missingObsForSite.set(s);
						}
					}
				}
				getGenotypeAfterUpdatingMissing();
			}
		}
		
		//calculate the minimum class size
		//if two classes min class size = the smaller of the two class counts
		//if three classes return second largest site count
		List<Integer> classSizes = new ArrayList<>(genoCountMap.values());
		Collections.sort(classSizes);
		int nclasses = classSizes.size();
		if (nclasses > 1) myCurrentSiteMinimumClassSize = classSizes.get(nclasses - 2);
		else myCurrentSiteMinimumClassSize = 0;
		
		return keepSite;
	}
	
	@Override
	public TableReport siteReport() {
		saveToFile = true;
		return siteReportBuilder.build();
	}
	
	@Override
	public TableReport alleleReport() {
		saveToFile = true;
		return alleleReportBuilder.build();
	}
	
	@Override
	public List<Datum> datumList() {
		List<Datum> dataList = new ArrayList<Datum>();
		StringBuilder comment = new StringBuilder();
		comment.append("GLM Output\nStatistical Tests for individual variants.\n");
		comment.append("Input data: " + myDatum.getName()).append("\n");
		dataList.add(new Datum("GLM_Stats_" + myDatum.getName(), siteReport(), comment.toString()));
		comment = new StringBuilder();
		comment.append("GLM Output\nGenotype Effect Estimates\n");
		comment.append("Input data: " + myDatum.getName()).append("\n");
		dataList.add(new Datum("GLM_Genotypes_" + myDatum.getName(), alleleReport(), comment.toString()));
		return dataList;
	}
	
	@Override
	public void permutationTest(boolean permute, int nperm) {
		this.permute = permute;
		numberOfPermutations = nperm;
	}

	/**
	 * @param missingObsBeforeSite	a BitSet with bits set for observations missing in model covariates and data
	 */
	protected abstract void getGenotypeAndUpdateMissing(BitSet missingObsBeforeSite);
	
	/**
	 * 	updates the genotype after missingObsForSite has changed
	 */
	protected abstract void getGenotypeAfterUpdatingMissing();
	
	/**
	 * @param siteNumber		a site number
	 * This method tests the significance of this site and estimates the allele effects then appends the results to the site and allele reports.
	 */
	protected abstract void analyzeSite();
	
	protected String[] siteReportColumnNames() {
		markerpvalueColumn = 5;
		permpvalueColumn = 6;
		if (appendAddDomEffects && !permute) return new String[] {AssociationConstants.STATS_HEADER_TRAIT,AssociationConstants.STATS_HEADER_MARKER,AssociationConstants.STATS_HEADER_CHR,AssociationConstants.STATS_HEADER_POSITION,"marker_F",AssociationConstants.STATS_HEADER_P_VALUE,"marker_Rsq","add_F","add_p","dom_F","dom_p", "marker_df","marker_MS","error_df","error_MS","model_df","model_MS","minorObs", "addEffect", "domEffect"};
		if (appendAddDomEffects && permute) return new String[]{AssociationConstants.STATS_HEADER_TRAIT,AssociationConstants.STATS_HEADER_MARKER,AssociationConstants.STATS_HEADER_CHR,AssociationConstants.STATS_HEADER_POSITION,"marker_F",AssociationConstants.STATS_HEADER_P_VALUE,"perm_p","marker_Rsq","add_F","add_p","dom_F","dom_p", "marker_df","marker_MS","error_df","error_MS","model_df","model_MS","minorObs", "addEffect", "domEffect"};
		if (permute) return new String[]{AssociationConstants.STATS_HEADER_TRAIT,AssociationConstants.STATS_HEADER_MARKER,AssociationConstants.STATS_HEADER_CHR,AssociationConstants.STATS_HEADER_POSITION,"marker_F",AssociationConstants.STATS_HEADER_P_VALUE,"perm_p","marker_Rsq","add_F","add_p","dom_F","dom_p", "marker_df","marker_MS","error_df","error_MS","model_df","model_MS","minorObs"};
		return new String[] {AssociationConstants.STATS_HEADER_TRAIT,AssociationConstants.STATS_HEADER_MARKER,AssociationConstants.STATS_HEADER_CHR,AssociationConstants.STATS_HEADER_POSITION,"marker_F",AssociationConstants.STATS_HEADER_P_VALUE,"marker_Rsq","add_F","add_p","dom_F","dom_p", "marker_df","marker_MS","error_df","error_MS","model_df","model_MS","minorObs"};
	}
	
	protected String[] alleleReportColumnNames() {
		return new String[]{AssociationConstants.STATS_HEADER_TRAIT,AssociationConstants.STATS_HEADER_MARKER,AssociationConstants.STATS_HEADER_CHR,AssociationConstants.STATS_HEADER_POSITION,"Obs","Allele","Estimate"};
	}
	
	protected ArrayList<ModelEffect> baseModel() {
		int numberOfNonmissingObs = numberOfObservations - (int) missingObsForSite.cardinality();
		ArrayList<ModelEffect> modelEffects = new ArrayList<ModelEffect>();
		FactorModelEffect meanEffect = new FactorModelEffect(new int[numberOfNonmissingObs], false);
		meanEffect.setID("mean");
		modelEffects.add(meanEffect);
		
		//add factors to model
		for (PhenotypeAttribute attr:myFactorAttributes) {
			String[] factorLabels = AssociationUtils.getNonMissingValues((String[]) attr.allValues(), missingObsForSite);
			FactorModelEffect fme = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(factorLabels), true, attr.name());
			modelEffects.add(fme);
		}

		//add covariates to model
		for (PhenotypeAttribute attr:myCovariateAttributes) {
			double[] values = AssociationUtils.getNonMissingDoubles((float[]) attr.allValues(), missingObsForSite);
			CovariateModelEffect cme = new CovariateModelEffect(values, attr.name());
			modelEffects.add(cme);
		}

		return modelEffects;
	}
	
	protected void createPermutedData() {
		permutedData = new LinkedList<>();
		double[] y = AssociationUtils.getNonMissingDoubles(allData, missingObsForSite);
		SweepFastLinearModel sflm = new SweepFastLinearModel(baseModel(), y);
		DoubleMatrix residuals = sflm.getResiduals();
		DoubleMatrix predicted = sflm.getPredictedValues();
		baseErrorSSdf = sflm.getResidualSSdf();
		totalcfmSSdf = new double[2];
		double[] modelSSdf = sflm.getModelcfmSSdf();
		totalcfmSSdf[0] = baseErrorSSdf[0] + modelSSdf[0];
		totalcfmSSdf[1] = baseErrorSSdf[1] + modelSSdf[1];
		
		if (rand == null) {
			if (useRandomSeed) rand = new Random(randomSeed);
			else rand = new Random();
		}
		
		for (int p = 0; p < numberOfPermutations; p++) {
			LinearModelUtils.shuffle(residuals, rand);
			DoubleMatrix permdm = predicted.plus(residuals);
			permutedData.add(permdm);
		}
		
		minP = new double[numberOfPermutations];
		Arrays.fill(minP, 1.0);
	}

	protected void updateReportsWithPermutationP() {
		Arrays.sort(minP);
		for (Object[] row : siteTableReportRows) {
			double pval = ((Double) row[markerpvalueColumn]);
			int ndx = Arrays.binarySearch(minP, pval);
			
			if (ndx < 0) ndx = -(ndx + 1);
			double permPval = (double) (ndx + 1) / (double) numberOfPermutations;
			if (permPval > 1) permPval = 1;
			row[permpvalueColumn] = new Double(permPval);
			
		}
	}

	protected ModelEffect taxaEffect() {
        Taxon[] myTaxa = myGenoPheno.phenotype().taxaAttribute().allTaxa();
        Taxon[] myNonMissingTaxa = AssociationUtils.getNonMissingValues(myTaxa, missingObsForSite);
        int[] taxaLevels = ModelEffectUtils.getIntegerLevels(myNonMissingTaxa);
        FactorModelEffect taxaEffect = new FactorModelEffect(taxaLevels, true, "Taxon");
        return taxaEffect;
	}
	
//	protected void testTaxaReplication() {
//		int numberOfObservations = myGenoPheno.phenotype().numberOfObservations();
//		int numberOfTaxa = myGenoPheno.genotypeTable().numberOfTaxa();
//		if (numberOfTaxa < numberOfObservations) areTaxaReplicated = true;
//		else areTaxaReplicated = false;
//	}
	
	protected void testTaxaReplication() {
		areTaxaReplicated = false;
		int numberOfObservations = myGenoPheno.phenotype().numberOfObservations();
		int numberOfTaxa = myGenoPheno.genotypeTable().numberOfTaxa();
		if (numberOfTaxa < numberOfObservations) {
//			String msg = "Taxa are duplicated in the phenotype data set. Tassel version 5 will not run GLM when that is the case.";
			String msg = "There are more phenotype observations than taxa with genotypes. Either some taxa have multiple phenotypes or some taxa do not have genotypes. Tassel version 5 will not run GLM when that is the case. Be sure to use an intersect join to merge genotypes and phenotypes.";
//			myLogger.error(msg);
			throw new RuntimeException(msg);
		}
	}
	
	protected void updateMinP(BitSet missingObsBeforeSite) {
		boolean useFastMethod = false;
		int numberOfObsTotal = allData.length;
		int numberOfMissingBeforeSite = (int) missingObsBeforeSite.cardinality();
		int sizeOfPermutedData = permutedData.get(0).numberOfRows();
		BitSet newMissing = new OpenBitSet(sizeOfPermutedData);
		int permutedDataIndex = -1;
		for (int i = 0; i < numberOfObsTotal; i++) {
			if (!missingObsBeforeSite.fastGet(i)) {
				permutedDataIndex++;
				if (missingObsForSite.fastGet(i)) newMissing.fastSet(permutedDataIndex);
			}
		}
		
		if (areTaxaReplicated) {
			//if taxa are replicated
			int iter = 0;
			for (DoubleMatrix pdata : permutedData) {
				double[] y = AssociationUtils.getNonMissingDoubles(pdata.to1DArray(), newMissing);
		        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, y);
		        markerSSdf = sflm.getIncrementalSSdf(numberOfBaseEffects);
		        errorSSdf = sflm.getIncrementalSSdf(taxaEffectNumber);
		        double F = markerSSdf[0] / markerSSdf[1] / errorSSdf[0] * errorSSdf[1];
		        double p;
		        try {
		        	p = LinearModelUtils.Ftest(F, markerSSdf[1], errorSSdf[1]);
		        	if (minP[iter] > p) minP[iter] = p;
		        } catch (Exception e) {
		        	//do nothing
		        }
		        iter++;
			}
		} else if (useFastMethod) {
			int numberOfModelEffects = myModel.size();
			List<ModelEffect> thisBaseModel = new ArrayList<>(myModel);
			ModelEffect markerEffect = thisBaseModel.remove(myModel.size() - 1);
			List<double[]> permutedArrays = permutedData.stream()
					.map(dm -> dm.to1DArray())
					.map(da -> AssociationUtils.getNonMissingDoubles(da, newMissing))
					.collect(Collectors.toList());
			SolveByOrthogonalizing sbo = SolveByOrthogonalizing.getInstanceFromModel(thisBaseModel, permutedArrays);
			DoubleMatrix X = markerEffect.getX();
			SolveByOrthogonalizing.Marker markerRValues = null;
			if (X.numberOfColumns() == 1) {
				markerRValues = sbo.solveForR(null, X.to1DArray());
			} else if (X.numberOfColumns() == 2) {
				markerRValues = sbo.solveForR(null, X.column(0).to1DArray(), X.column(1).to1DArray());
			}
			
			if (markerRValues != null) {
				int n = X.numberOfRows();
				for (int iter = 0; iter < numberOfPermutations; iter++) {
					if (minP[iter] > markerRValues.vector2()[iter]) minP[iter] = markerRValues.vector2()[iter];
				}
			}
			
		} else {
			int iter = 0;
			int numberOfModelEffects = myModel.size();
			for (DoubleMatrix pdata : permutedData) {
				double[] y = AssociationUtils.getNonMissingDoubles(pdata.to1DArray(), newMissing);
		        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, y);
		        markerSSdf = sflm.getIncrementalSSdf(numberOfModelEffects - 1);
		        errorSSdf = sflm.getResidualSSdf();
		        double F = markerSSdf[0] / markerSSdf[1] / errorSSdf[0] * errorSSdf[1];
		        double p;
		        try {
		        	p = LinearModelUtils.Ftest(F, markerSSdf[1], errorSSdf[1]);
		        	if (minP[iter] > p) minP[iter] = p;
		        } catch (Exception e) {
		        	//do nothing
		        }
				
				iter++;
			}
		}
		
	}
	
	@Override
	public void maxP(double maxP) {
		this.maxP = maxP;
	}

	@Override
	public void siteReportFilepath(String savefile) {
		saveToFile = true;
		siteReportFilename = savefile;
		
	}

	@Override
	public void alleleReportFilepath(String savefile) {
		saveToFile = true;
		alleleReportFilename = savefile;
	}

	@Override
	public void biallelicOnly(boolean biallelic) {
		biallelicOnly = biallelic;
	}

	@Override
	public void minimumClassSize(int minsize) {
		minClassSize = minsize;
	}

	@Override
	public void saveSiteStats(boolean siteStats) {
		outputSiteStats = siteStats;
	}

	@Override
	public void siteStatsFile(String filename) {
		siteStatsFile = filename;
	}

	@Override
	public void appendAddDom(boolean append) {
		appendAddDomEffects = append;
	}

	/**
	 * This method is used mainly for testing in order to generate reproducible permutation results. 
	 * If the seed is not set, the current time is used to initialize the random number generator.
	 * @param seed	the seed used to initialize the random number generator used by permutation
	 */
	public void setRandomSeed(int seed) {
		randomSeed = seed;
		useRandomSeed = true;
	}
}
