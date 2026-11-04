package net.maizegenetics.analysis.association;

import java.util.ArrayList;


import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.LinearModelUtils;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffectUtils;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

public class DiscreteSitesFELM extends AbstractFixedEffectLM {
	String[] siteGenotypes;

	public DiscreteSitesFELM(Datum dataset, FixedEffectLMPlugin parentPlugin) {
		super(dataset, parentPlugin);
	}

	@Override
	protected void analyzeSite() {
		myModel = new ArrayList<ModelEffect>(myBaseModel);
		GenotypeTable myGenotype = myGenoPheno.genotypeTable();
		
		//solve the full model
		//add the marker to the model
        String[] siteGenotypes = AssociationUtils.getNonMissingValues(myGenoPheno.getStringGenotype(myCurrentSite), missingObsForSite);
        ArrayList<String> markerIds = new ArrayList<>();
        int[] markerLevels = ModelEffectUtils.getIntegerLevelsSortedByGenotype(siteGenotypes, markerIds);
        String siteName = myGenotype.siteName(myCurrentSite);
        FactorModelEffect markerEffect = new FactorModelEffect(markerLevels, true, siteName);
        int firstMarkerParm = myModel.size();
        myModel.add(markerEffect);
        
        //add taxa:marker effect at end
        taxaEffectNumber = -1;
        if (areTaxaReplicated) {
        	taxaEffectNumber = myModel.size();
        	myModel.add(taxaEffect());
        }

        //calculate model
        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, siteData);
        double[] modelSSdf = sflm.getModelcfmSSdf();
        double[] residSSdf = sflm.getResidualSSdf();
        double[] totalSSdf = new double[2];
        totalSSdf[0] = modelSSdf[0] + residSSdf[0];
        totalSSdf[1] = modelSSdf[1] + residSSdf[1];
        
        markerSSdf = sflm.getIncrementalSSdf(numberOfBaseEffects);
        if (areTaxaReplicated) errorSSdf = sflm.getIncrementalSSdf(taxaEffectNumber);
        else errorSSdf = residSSdf;
        double F = markerSSdf[0] / markerSSdf[1] / errorSSdf[0] * errorSSdf[1];
        double p;
        try {
            p = LinearModelUtils.Ftest(F, markerSSdf[1], errorSSdf[1]);
        } catch (Exception e) {
            p = Double.NaN;
        }
        double[] beta = sflm.getBeta();
		G = null;
		if (permute) G = sflm.getInverseOfXtX();
		
		double domEffect = Double.NaN, addEffect = Double.NaN;
		if (appendAddDomEffects) {
			//test for bi-allelic and hets present
			if (myGenotype.alleles(myCurrentSite).length == 2) {
				if (markerIds.size() == 2) {
					addEffect = Math.abs(beta[firstMarkerParm] / 2);
				} else if (markerIds.size() == 3) {
					addEffect = Math.abs((beta[firstMarkerParm] -  beta[firstMarkerParm + 1])/ 2);
					domEffect = -(beta[firstMarkerParm] + beta[firstMarkerParm + 1])/ 2;
				}
			}
		}
		
		//calculate additive and dominance F and p
		ArrayList<ModelEffect> myAdditiveModel = new ArrayList<>(myBaseModel);
		byte[] alleles = myGenotype.alleles(myCurrentSite);
		int nAlleles = alleles.length;
		double Fadd, padd, Fdom, pdom;
		if (nAlleles < 2) {
			Fadd = Double.NaN;
			padd = Double.NaN;
			Fdom = Double.NaN;
			pdom = Double.NaN;
		} else {
			for (int a = 0; a < nAlleles - 1; a++) {
				String stringAllele = NucleotideAlignmentConstants.getHaplotypeNucleotide(alleles[a]);
				myAdditiveModel.add(new CovariateModelEffect(ModelEffectUtils.getNumericCodingForAdditiveModel(siteGenotypes, stringAllele)));
			}
			
	        //add taxa:marker effect at end
	        int addModelTaxaEffectNumber = -1;
	        if (areTaxaReplicated) {
	        	addModelTaxaEffectNumber = myAdditiveModel.size();
	        	myAdditiveModel.add(taxaEffect());
	        }
			
			SweepFastLinearModel sflmAdd = new SweepFastLinearModel(myAdditiveModel, siteData);
			double[] additiveErrorSSdf;
	        if (areTaxaReplicated) additiveErrorSSdf = sflmAdd.getIncrementalSSdf(addModelTaxaEffectNumber);
	        else additiveErrorSSdf = sflmAdd.getResidualSSdf();

			double[] addTermSSdf = new double[2];
			for (int a = 0; a < nAlleles - 1; a++) {
				double[] thisTermSSdf = sflmAdd.getIncrementalSSdf(a + numberOfBaseEffects);
				addTermSSdf[0] += thisTermSSdf[0];
				addTermSSdf[1] += thisTermSSdf[1];
			}
			
	        Fadd = addTermSSdf[0] / addTermSSdf[1] / additiveErrorSSdf[0] * additiveErrorSSdf[1];
	        
	        if (Double.isFinite(Fadd)) {
		        try {
		        	padd = LinearModelUtils.Ftest(Fadd, addTermSSdf[1], additiveErrorSSdf[1]);
		        } catch (Exception e) {
		        	padd = Double.NaN;
		        }
	        } else {
	        	padd = Double.NaN;
	        	Fadd = Double.NaN;
	        }

			//dominance term, F = (reduction in SS Error from full model compared to additive model)/ reduction in df / error MS for full model
	        double[] domTermSSdf = new double[2];
	        domTermSSdf[0] = markerSSdf[0] - addTermSSdf[0];
	        domTermSSdf[1] = markerSSdf[1] - addTermSSdf[1];
	        Fdom = domTermSSdf[0] / domTermSSdf[1] / errorSSdf[0] * errorSSdf[1];
	        if (Double.isFinite(Fdom)) {
		        try {
		        	pdom = LinearModelUtils.Ftest(Fdom, domTermSSdf[1], errorSSdf[1]);
		        } catch (Exception e) {
		        	pdom = Double.NaN;
		        }
	        } else {
	        	pdom = Double.NaN;
	        	Fdom = Double.NaN;
	        }

		}

		//add results if p <= maxP
        //add results to site report
        // column names {"Trait","Marker","Chr","Position","marker_F","marker_p","perm_p","marker_Rsq","add_F","add_p","dom_F","dom_p", "marker_df","marker_MS","error_df","error_MS","model_df","model_MS" }
        if (maxP == 1.0 || p <= maxP) {
    		Object[] rowData = new Object[numberOfSiteReportColumns];
            int columnCount = 0;
            rowData[columnCount++] = currentTraitName;
            rowData[columnCount++] = siteName;
            rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomeName(myCurrentSite);
            rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomalPosition(myCurrentSite);
            rowData[columnCount++] = Double.valueOf(F);
            rowData[columnCount++] = Double.valueOf(p);
            if (permute) rowData[columnCount++] = "";
            rowData[columnCount++] = Double.valueOf(markerSSdf[0]/totalSSdf[0]);
            rowData[columnCount++] = Double.valueOf(Fadd);
            rowData[columnCount++] = Double.valueOf(padd);
            rowData[columnCount++] = Double.valueOf(Fdom);
            rowData[columnCount++] = Double.valueOf(pdom);
            rowData[columnCount++] = Double.valueOf(markerSSdf[1]);
            rowData[columnCount++] = Double.valueOf(markerSSdf[0]/markerSSdf[1]);
            rowData[columnCount++] = Double.valueOf(errorSSdf[1]);
            rowData[columnCount++] = Double.valueOf(errorSSdf[0]/errorSSdf[1]);
            rowData[columnCount++] = Double.valueOf(modelSSdf[1]);
            rowData[columnCount++] = Double.valueOf(modelSSdf[0]/modelSSdf[1]);
            rowData[columnCount++] = Integer.valueOf(myCurrentSiteMinimumClassSize);
            if (appendAddDomEffects) {
            	rowData[columnCount++] = Double.valueOf(addEffect);
            	rowData[columnCount++] = Double.valueOf(domEffect);
            }
            siteReportBuilder.add(rowData);
            if (permute) siteTableReportRows.add(rowData);
            
            //add results to allele report if nAlleles > 1
            if (nAlleles > 1) {
            	int numberOfAlleles = markerIds.size();
            	int[] alleleCounts = markerEffect.getLevelCounts();
            	int firstEstimateIndex = beta.length - numberOfAlleles + 1;
            	for (int a = 0; a < numberOfAlleles; a++) {
            		rowData = new Object[numberOfAlleleReportColumns];
            		columnCount = 0;
            		rowData[columnCount++] = currentTraitName;
            		rowData[columnCount++] = siteName;
            		rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomeName(myCurrentSite);
            		rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomalPosition(myCurrentSite);
            		rowData[columnCount++] = Integer.valueOf(alleleCounts[a]);
            		rowData[columnCount++] = markerIds.get(a);
            		if (a < numberOfAlleles - 1) rowData[columnCount++] = Double.valueOf(beta[firstEstimateIndex + a]);
            		else rowData[columnCount++] = Double.valueOf(0.0);
            		alleleReportBuilder.add(rowData);
            	}
            }
        }
	}

	@Override
	protected void getGenotypeAndUpdateMissing(BitSet missingObsBeforeSite) {
		String[] allSiteGenotypes = myGenoPheno.getStringGenotype(myCurrentSite);
		int n = allSiteGenotypes.length;
		missingObsForSite = new OpenBitSet(missingObsBeforeSite);
		for (int i = 0; i < n; i++) {
			if (allSiteGenotypes[i].contains("N")) missingObsForSite.fastSet(i);
		}
		siteGenotypes = AssociationUtils.getNonMissingValues(allSiteGenotypes, missingObsForSite);
	}

	@Override
	protected void getGenotypeAfterUpdatingMissing() {
		siteGenotypes = AssociationUtils.getNonMissingValues(myGenoPheno.getStringGenotype(myCurrentSite), missingObsForSite);
	}

	
}
