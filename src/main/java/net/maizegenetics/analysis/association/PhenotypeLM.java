package net.maizegenetics.analysis.association;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.LinearModelUtils;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffectUtils;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportBuilder;

public class PhenotypeLM {
	protected Datum myDatum;
	protected Phenotype myPhenotype;
	protected Phenotype myBlues;
	protected TableReportBuilder reportBuilder;
	private TableReport report;
	boolean areTaxaReplicated;
	protected TaxaList myTaxaList;
	protected int numberOfObservations;
	protected List<PhenotypeAttribute> myDataAttributes;
	protected List<PhenotypeAttribute> myFactorAttributes;
	protected List<PhenotypeAttribute> myCovariateAttributes;
	protected Taxon[] myTaxa;
	protected ArrayList<Taxon> taxaInModel;
	protected Map<Taxon, Integer> taxaInModelMap;
	
	public PhenotypeLM(Datum phenotypeOnly) {
		myDatum = phenotypeOnly;
		Object data = myDatum.getData();
		if (data instanceof Phenotype) myPhenotype = (Phenotype) data;
		else if (data instanceof GenotypePhenotype) myPhenotype = ((GenotypePhenotype) data).phenotype();
		else return;
		
		initialize();
		solve();
	}
	
	public Phenotype blues() {
		return myBlues;
	}
	
	public TableReport report() {
		if (report == null) report = reportBuilder.build();
		return report;
	}
	
	public List<Datum> datumList() {
		ArrayList<Datum> returnDatum = new ArrayList<Datum>();
		String bluesName = "BLUEs_" + myDatum.getName();
		StringBuilder sb = new StringBuilder("Best Linear Unbiased Estimates\n");
		sb.append("From ").append(myDatum.getName());
		sb.append("\nNumber of Taxa = ").append(myTaxaList.size());
		sb.append("\nNumber of Traits = ").append(myBlues.numberOfAttributes() - 1);
		String bluesComment = sb.toString();
		String reportName = "Phenotype_ANOVA_" + myDatum.getName();
		sb = new StringBuilder("Taxa Statistical Tests");
		sb.append("From ").append(myDatum.getName());
		String reportComment = sb.toString();
		returnDatum.add(new Datum(bluesName, myBlues, bluesComment));
		returnDatum.add(new Datum(reportName, report(), reportComment));
		return returnDatum;
	}
	
	protected void initialize() {
		myTaxaList = myPhenotype.taxa();
		myTaxa = myPhenotype.taxaAttribute().allTaxa();

		String[] columns = new String[]{"Trait", "F", "p", "taxaDF", "taxaMS", "errorDF", "errorMS", "modelDF", "modelMS"};
		String name = "Phenotype analysis for " + myDatum.getName();
		reportBuilder = TableReportBuilder.getInstance(name, columns);
		numberOfObservations = myPhenotype.numberOfObservations();
		
		myDataAttributes = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data);
		myFactorAttributes = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor);
		myCovariateAttributes = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.covariate);
		
	}
	
	protected void solve() {
		
		OpenBitSet missingModelObs = new OpenBitSet(numberOfObservations);
		for (PhenotypeAttribute attr:myFactorAttributes) missingModelObs.or(attr.missing());
		for (PhenotypeAttribute attr:myCovariateAttributes) missingModelObs.or(attr.missing());
		List<PhenotypeAttribute> blueAttributes = new ArrayList<PhenotypeAttribute>();
		List<ATTRIBUTE_TYPE> blueAttributeTypes = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		blueAttributes.add(new TaxaAttribute(myTaxaList, myPhenotype.taxaAttribute().name()));
		blueAttributeTypes.add(ATTRIBUTE_TYPE.taxa);
		
		//iterate through phenotypes, create blues, add to report
		for (PhenotypeAttribute dataAttribute : myDataAttributes) {
			String currentTraitName = dataAttribute.name();
			
			//determine which observations have missing data
			OpenBitSet missingObs = new OpenBitSet(dataAttribute.missing());
			for (PhenotypeAttribute attr:myFactorAttributes) missingObs.or(attr.missing());
			for (PhenotypeAttribute attr:myCovariateAttributes) missingObs.or(attr.missing());
			
			//get y
			float[] allData = (float[]) dataAttribute.allValues();
			double[] y = AssociationUtils.getNonMissingDoubles(allData, missingObs);
			
			//build the model
			ArrayList<ModelEffect> myModel = model(missingObs, y.length);
			
			SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, y);
			
			//save BLUEs
			double[] beta = sflm.getBeta();
			
            // calculate the BLUEs for taxa
            // first calculate the average of the level estimates for all other factors (including the mean)
            double overallMean = beta[0]; //the mean
            int nEffects = myModel.size();
            int start = 1;
            for (int i = 1; i < nEffects - 1; i++) {
                ModelEffect me = myModel.get(i);
                if (me instanceof FactorModelEffect) {
                    FactorModelEffect fme = (FactorModelEffect) me;
                    int nLevels = fme.getNumberOfLevels();
                    int nEstimates;
                    if (fme.getRestricted()) {
                        nEstimates = nLevels - 1;
                    } else {
                        nEstimates = nLevels;
                    }
                    double factorMean = 0;
                    for (int j = 0; j < nEstimates; j++) {
                        factorMean += beta[j + start];
                    }
                    factorMean /= nLevels;
                    overallMean += factorMean;
                    start += nEstimates;
                } else {
                    start += me.getNumberOfLevels();
                }
            }

            float[] attrValues = new float[myTaxaList.size()];
            OpenBitSet missing = new OpenBitSet(myTaxaList.size());
            int lastTaxon = taxaInModel.size() - 1;
            for (int t = 0; t < myTaxaList.size(); t++) {
            	Integer ndx = taxaInModelMap.get(myTaxaList.get(t));
            	if (ndx == null) {
            		attrValues[t] = Float.NaN;
            		missing.fastSet(t);
            	} else if (ndx == lastTaxon) {
            		attrValues[t] = (float) overallMean; 
            	} else {
            		attrValues[t] = (float) (beta[ndx + start] + overallMean) ;
            	}
            }
            
            blueAttributes.add(new NumericAttribute(dataAttribute.name(), attrValues, missing));
            blueAttributeTypes.add(ATTRIBUTE_TYPE.data);
            
            //add row to report
            //record {"Trait", "F", "p", "taxaDF", "taxaMS", "errorDF", "errorMS", "modelDF", "modelMS"}
            double[] taxaSSdf = sflm.getIncrementalSSdf(nEffects - 1);
            double[] modelSSdf = sflm.getModelcfmSSdf();
            double[] errorSSdf = sflm.getResidualSSdf();
            
            double F, p;
            F = taxaSSdf[0] / taxaSSdf[1] / errorSSdf[0] * errorSSdf[1];
            try {
                p = LinearModelUtils.Ftest(F, taxaSSdf[1], errorSSdf[1]);
            } catch (Exception e) {
                p = Double.NaN;
            }

            Object[] result = new Object[9];
            result[0] = currentTraitName;
            result[1] = new Double(F);
            result[2] = new Double(p);
            result[3] = new Double(taxaSSdf[1]);
            result[4] = new Double(taxaSSdf[0] / taxaSSdf[1]);
            result[5] = new Double(errorSSdf[1]);
            result[6] = new Double(errorSSdf[0] / errorSSdf[1]);
            result[7] = new Double(modelSSdf[1]);
            result[8] = new Double(modelSSdf[0] / modelSSdf[1]);
            reportBuilder.add(result);

		}
		
		myBlues = new PhenotypeBuilder().fromAttributeList(blueAttributes, blueAttributeTypes).build().get(0);
	}
	
	protected void testTaxaReplication() {
		if (myTaxaList.numberOfTaxa() < myPhenotype.numberOfObservations()) areTaxaReplicated = true;
		else areTaxaReplicated = false;
	}
	
	protected ArrayList<ModelEffect> model(BitSet missingObs, int numberOfNonmissingObs) {
		ArrayList<ModelEffect> modelEffects = new ArrayList<ModelEffect>();
		FactorModelEffect meanEffect = new FactorModelEffect(new int[numberOfNonmissingObs], false, "mean");
		modelEffects.add(meanEffect);
		
		//add factors to model
		for (PhenotypeAttribute attr:myFactorAttributes) {
			String[] factorLabels = AssociationUtils.getNonMissingValues((String[]) attr.allValues(), missingObs);
			FactorModelEffect fme = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(factorLabels), true, attr.name());
			modelEffects.add(fme);
		}

		//add covariates to model
		for (PhenotypeAttribute attr:myCovariateAttributes) {
			NumericAttribute numericAttr = (NumericAttribute) attr;
			double[] values = AssociationUtils.getNonMissingDoubles(numericAttr.floatValues(), missingObs);
			CovariateModelEffect cme = new CovariateModelEffect(values, attr.name());
			modelEffects.add(cme);
		}

		//add Taxa to model
		Taxon[] taxa = AssociationUtils.getNonMissingValues(myTaxa, missingObs);
		taxaInModel = new ArrayList<Taxon>();
		int[] taxaLevels = ModelEffectUtils.getIntegerLevels(taxa, taxaInModel);
		FactorModelEffect taxaEffect = new FactorModelEffect(taxaLevels, true);
		modelEffects.add(taxaEffect);
		
		taxaInModelMap = new HashMap<>();
		int count = 0;
		for (Taxon taxon : taxaInModel) taxaInModelMap.put(taxon, count++);
		
		return modelEffects;
	}
	
}
