package net.maizegenetics.analysis.association;

import java.util.ArrayList;

import net.maizegenetics.dna.snp.score.SiteScore.SITE_SCORE_TYPE;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.LinearModelUtils;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

public class ReferenceProbabilityFELM extends AbstractFixedEffectLM {
    double[] myProbabilities;

    public ReferenceProbabilityFELM(Datum data, FixedEffectLMPlugin parentPlugin) {
        super(data, parentPlugin);
    }

    @Override
    protected void analyzeSite() {
        myModel = new ArrayList<ModelEffect>(myBaseModel);
        String siteName = myGenoPheno.genotypeTable().siteName(myCurrentSite);

        myModel.add(new CovariateModelEffect(myProbabilities));

        if (areTaxaReplicated)
            myModel.add(taxaEffect());

        //solve the model
        SweepFastLinearModel markerModel = new SweepFastLinearModel(myModel, siteData);

        //calculate model
        double[] modelSSdf = markerModel.getModelcfmSSdf();
        markerSSdf = markerModel.getIncrementalSSdf(numberOfBaseEffects);
        if (areTaxaReplicated)
            errorSSdf = markerModel.getIncrementalSSdf(numberOfBaseEffects + 1);
        else
            errorSSdf = markerModel.getResidualSSdf();

        double rsq = markerSSdf[0] / (modelSSdf[0] + markerModel.getResidualSSdf()[0]);

        double F = markerSSdf[0] / markerSSdf[1] / errorSSdf[0] * errorSSdf[1];
        double p;
        try {
            p = LinearModelUtils.Ftest(F, markerSSdf[1], errorSSdf[1]);
        } catch (Exception e) {
            p = Double.NaN;
        }
        double[] beta = markerModel.getBeta();
        if (permute)
            G = markerModel.getInverseOfXtX();

        //add results to site report
        //{"Trait","Marker","Chr","Pos","marker_F","p","marker_Rsq","marker_df","marker_MS","error_df","error_MS","model_df","model_MS" }
        if (maxP == 1.0 || p <= maxP) {
            Object[] rowData = new Object[numberOfSiteReportColumns];
            int columnCount = 0;
            rowData[columnCount++] = currentTraitName;
            rowData[columnCount++] = siteName;
            rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomeName(myCurrentSite);
            rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomalPosition(myCurrentSite);
            rowData[columnCount++] = new Double(F);
            rowData[columnCount++] = new Double(p);
            if (permute)
                rowData[columnCount++] = "";
            rowData[columnCount++] = new Double(rsq);
            rowData[columnCount++] = new Double(markerSSdf[1]);
            rowData[columnCount++] = new Double(markerSSdf[0] / markerSSdf[1]);
            rowData[columnCount++] = new Double(errorSSdf[1]);
            rowData[columnCount++] = new Double(errorSSdf[0] / errorSSdf[1]);
            rowData[columnCount++] = new Double(modelSSdf[1]);
            rowData[columnCount++] = new Double(modelSSdf[0] / modelSSdf[1]);
            siteReportBuilder.add(rowData);
            if (permute)
                siteTableReportRows.add(rowData);

            //add results to allele report
            //{"Trait","Marker","Chr","Position","Estimate"}
            int estimateIndex = beta.length - 1;
            rowData = new Object[numberOfAlleleReportColumns];
            columnCount = 0;
            rowData[columnCount++] = currentTraitName;
            rowData[columnCount++] = siteName;
            rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomeName(myCurrentSite);
            rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomalPosition(myCurrentSite);
            rowData[columnCount++] = beta[estimateIndex];
            alleleReportBuilder.add(rowData);
        }
    }

    @Override
    protected void getGenotypeAndUpdateMissing(BitSet missingObsBeforeSite) {
        float[] allSiteProbs = myGenoPheno.referenceProb(myCurrentSite);

        int n = allSiteProbs.length;
        missingObsForSite = new OpenBitSet(missingObsBeforeSite);
        for (int i = 0; i < n; i++) {
            if (Float.isNaN(allSiteProbs[i]))
                missingObsForSite.fastSet(i);
        }
        myProbabilities = AssociationUtils.getNonMissingDoubles(allSiteProbs, missingObsForSite);
    }

    @Override
	protected void getGenotypeAfterUpdatingMissing() {
    	myProbabilities = AssociationUtils.getNonMissingDoubles(myGenoPheno.referenceProb(myCurrentSite), missingObsForSite);
	}

	@Override
    protected String[] siteReportColumnNames() {
        markerpvalueColumn = 5;
        permpvalueColumn = 6;
        if (permute)
            return new String[] { AssociationConstants.STATS_HEADER_TRAIT,
                    AssociationConstants.STATS_HEADER_MARKER,
                    AssociationConstants.STATS_HEADER_CHR,
                    AssociationConstants.STATS_HEADER_POSITION, "marker_F",
                    AssociationConstants.STATS_HEADER_P_VALUE,
                    "perm_p", "marker_Rsq", "marker_df", "marker_MS", "error_df", "error_MS",
                    "model_df", "model_MS" };
        return new String[] { AssociationConstants.STATS_HEADER_TRAIT,
                AssociationConstants.STATS_HEADER_MARKER, AssociationConstants.STATS_HEADER_CHR,
                AssociationConstants.STATS_HEADER_POSITION, "marker_F",
                AssociationConstants.STATS_HEADER_P_VALUE, "marker_Rsq", "marker_df", "marker_MS",
                "error_df", "error_MS", "model_df", "model_MS" };

    }

    @Override
    protected String[] alleleReportColumnNames() {
        // TODO Auto-generated method stub
        return new String[] { AssociationConstants.STATS_HEADER_TRAIT,
                AssociationConstants.STATS_HEADER_MARKER,
                AssociationConstants.STATS_HEADER_CHR, AssociationConstants.STATS_HEADER_POSITION,
                "Estimate" };
    }

}
