package net.maizegenetics.analysis.imputation;

import java.util.Arrays;
import java.util.List;

import net.maizegenetics.analysis.imputation.EmissionProbability;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;

public class CrossProgenyEmissionMatrix extends EmissionProbability {
	/*
	 * p(obs|state) matrix for each site 
	 * use depth to calculate p(obs) for each site/taxon
	 * 
	 * 4 state model, states are:
	 *  0 = h00/h10 [parent 0, haplotype 0 / parent 1, haplotype 0]
	 *  1 = h00/h11
	 *  2 = h01/h10
	 *  3 = h01/h11
	 *  obs are nucleotide genotypes
	 */
	
	private static final byte N = GenotypeTable.UNKNOWN_ALLELE;
	private static final byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
	
	private byte[][][] parentHaplotypes; //indexes are parent, chromosome, site
	private GenotypeTable progenyGenotypes;  //the original genotypes for the progeny of interest, not all taxa
	private double[][] parentHaplotypeProbabilities; //first index is in order h00, h01, h10, h11
	//TODO replace obsGivenState with method that accounts for depth
	private double[][] obsGivenState = new double[3][3]; //rows are imputed (assumed true), columns are observed
	private int myTaxonIndex;
	private int[] mySiteIndices;
	private boolean useProbabilityTables;
	
	/**
	 * @param hap	parent haplotypes
	 */
	public CrossProgenyEmissionMatrix(byte[][][] hap, GenotypeTable originalGenotypes, int taxonIndex, int[] siteIndices) {
		parentHaplotypes = hap;
		progenyGenotypes = originalGenotypes;
		useProbabilityTables = false;
		myTaxonIndex = taxonIndex;
		mySiteIndices = siteIndices;
	}
	
	public CrossProgenyEmissionMatrix(double[][] hapProbs, GenotypeTable originalGenotypes, double[][] stateByObs, int taxonIndex, int[] siteIndices) {
		parentHaplotypeProbabilities = hapProbs;
		progenyGenotypes = originalGenotypes;
		obsGivenState = stateByObs;
		useProbabilityTables = true;
		myTaxonIndex = taxonIndex;
		mySiteIndices = siteIndices;
	}
	
	@Override
	public double getProbObsGivenState(int state, int obs) {
		return Double.NaN;
	}

	@Override
	public double getProbObsGivenState(int state, int obs, int site) {
		if (useProbabilityTables) return probObsStateFromPhasedProgeny(state, obs, mySiteIndices[site]);
		return probObsStateUsingParentHaps(state, obs, mySiteIndices[site]);
	}
	
	private double probObsStateUsingParentHaps(int state, int obs, int site) {
		//obs = byte representing a diploid genotype as encoded by NucleotideAlignmentConstants
		
		//are the parent haplotypes known at this site
		//if not then base p(obs) on originalGenotypes
		byte obsgeno = (byte) obs;
		byte h00 = parentHaplotypes[0][0][site];
		byte h01 = parentHaplotypes[0][1][site];
		byte h10 = parentHaplotypes[1][0][site];
		byte h11 = parentHaplotypes[1][1][site];
		
//		if (h00 == N || h01 == N || h10 == N || h11 == N) return probObsAtSite(site, obs);
		if (h00 == N || h01 == N || h10 == N || h11 == N) return 0.25;
		
		byte h00h10 = GenotypeTableUtils.getUnphasedDiploidValue(h00, h10);
		byte h00h11 = GenotypeTableUtils.getUnphasedDiploidValue(h00, h11);
		byte h01h10 = GenotypeTableUtils.getUnphasedDiploidValue(h01, h10);
		byte h01h11 = GenotypeTableUtils.getUnphasedDiploidValue(h01, h11);
		
		//all parent haplotypes are known
		switch (state) {
		case 0:
			//state = h00/h10
			if (GenotypeTableUtils.isHeterozygous(h00h10)) {
				if (GenotypeTableUtils.isEqual(h00h10, obsgeno)) return 1 - probHetAsHom(site);
				else if (!GenotypeTableUtils.isHeterozygous(obsgeno))return probHetAsHom(site);
				else return 0.0001;
			} else {
				if (obsgeno == h00h10) return 0.998;
				else if (GenotypeTableUtils.isHeterozygous(obsgeno)) return 0.001;
				else return 0.001;
			}
		case 1:
			//state = h00/h11
			if (GenotypeTableUtils.isHeterozygous(h00h11)) {
				if (GenotypeTableUtils.isEqual(h00h11, obsgeno)) return 1 - probHetAsHom(site);
				else if (!GenotypeTableUtils.isHeterozygous(obsgeno)) return probHetAsHom(site);
				else return 0.0001;
			} else {
				if (obsgeno == h00h11) return 0.998;
				else if (GenotypeTableUtils.isHeterozygous(obsgeno)) return 0.001;
				else return 0.001;
			}
		case 2:
			//state = h01/h10
			if (GenotypeTableUtils.isHeterozygous(h01h10)) {
				if (GenotypeTableUtils.isEqual(h01h10, obsgeno)) return 1 - probHetAsHom(site);
				else if (!GenotypeTableUtils.isHeterozygous(obsgeno))return probHetAsHom(site);
				else return 0.0001;
			} else {
				if (obsgeno == h01h10) return 0.998;
				else if (GenotypeTableUtils.isHeterozygous(obsgeno)) return 0.001;
				else return 0.001;
			}
		case 3:
			//state = h01/h11
			if (GenotypeTableUtils.isHeterozygous(h01h11)) {
				if (GenotypeTableUtils.isEqual(h01h11, obsgeno)) return 1 - probHetAsHom(site);
				else if (!GenotypeTableUtils.isHeterozygous(obsgeno))return probHetAsHom(site);
				else return 0.0001;
			} else {
				if (obsgeno == h01h11) return 0.998;
				else if (GenotypeTableUtils.isHeterozygous(obsgeno)) return 0.001;
				else return 0.001;
			}
		default:
			throw new IllegalArgumentException("Illegal state as input to probObsStateUsingParentHaps().");
		}
	}
	
	private double probObsStateFromPhasedProgeny(int state, int obs, int site) {
		//obs = byte representing a diploid genotype as endoded by NucleotideAlignmentConstants
		
		byte obsgeno = (byte) obs;
		byte[] obsAlleles = GenotypeTableUtils.getDiploidValues(obsgeno);
		byte major = progenyGenotypes.majorAllele(site);
		byte minor = progenyGenotypes.minorAllele(site);
		
		//otherwise
		//0 = h0/h0, 1 = h0/h1, 2 = h1/h1
		//note if the code reaches this point, h0 != h1
		
		//P(obs|state=0) = P(obs|state=a) * P(state=a|state=0) + P(obs|state=h) * P(state=h|state=0) + P(obs|state=b) * P(state=b|state=0)
		double pHomMajor;
		double pHet;
		double pHomMinor;
		double pMajorh00 = parentHaplotypeProbabilities[0][site];
		double pMajorh01 = parentHaplotypeProbabilities[1][site];
		double pMajorh10 = parentHaplotypeProbabilities[2][site];
		double pMajorh11 = parentHaplotypeProbabilities[3][site];
		switch (state) {
		case 0:
			//state = h00/h10
			//following are probabilities given the state
			pHomMajor = pMajorh00 * pMajorh10;
			pHet = pMajorh00 *(1 - pMajorh10) + ( 1 - pMajorh00) * pMajorh10;
			pHomMinor = (1 - pMajorh00) *(1 - pMajorh10);
			break;
		case 1:
			//state = h00/h11
			pHomMajor = pMajorh00 * pMajorh11;
			pHet = pMajorh00 *(1 - pMajorh11) + ( 1 - pMajorh00) * pMajorh11;
			pHomMinor = (1 - pMajorh00) *(1 - pMajorh11);
			break;
		case 2:
			//state = h01/h10
			pHomMajor = pMajorh01 * pMajorh10;
			pHet = pMajorh01 *(1 - pMajorh10) + ( 1 - pMajorh01) * pMajorh10;
			pHomMinor = (1 - pMajorh01) *(1 - pMajorh10);
			break;
		case 3:
			//state = h01/h11
			pHomMajor = pMajorh01 * pMajorh11;
			pHet = pMajorh01 *(1 - pMajorh11) + ( 1 - pMajorh01) * pMajorh11;
			pHomMinor = (1 - pMajorh01) *(1 - pMajorh11);
			break;
		default:
			throw new IllegalArgumentException("Illegal state as input to probObsStateFromPhasedProgeny().");
		}
		
		if (obsgeno == GenotypeTableUtils.getDiploidValue(major, major)) {
			return obsGivenState[0][0] * pHomMajor + obsGivenState[1][0] * pHet + obsGivenState[2][0] * pHomMinor;
		} else if (GenotypeTableUtils.isEqual(obsgeno, GenotypeTableUtils.getDiploidValue(major, minor))) {
			return obsGivenState[0][1] * pHomMajor + obsGivenState[1][1] * pHet + obsGivenState[2][1] * pHomMinor;
		} else if (obsgeno == GenotypeTableUtils.getDiploidValue(minor, minor)){
			return obsGivenState[0][2] * pHomMajor + obsGivenState[1][2] * pHet + obsGivenState[2][2] * pHomMinor;
		} return Double.NaN;

	}
	
	private double probObsAtSite(int site, int obs) {
		double dtotal = progenyGenotypes.totalNonMissingForSite(site);
		byte obsgeno = (byte) obs;
		double pval;
		if (GenotypeTableUtils.isHeterozygous(obsgeno)) {
			int obscount = 0;
			for (int t = 0; t < progenyGenotypes.numberOfTaxa(); t++) {
				if (GenotypeTableUtils.isEqual(obsgeno, progenyGenotypes.genotype(t,site))) obscount++;
			}
			pval = obscount / dtotal;
		} else {
			int obscount = 0;
			for (int t = 0; t < progenyGenotypes.numberOfTaxa(); t++) {
				if (progenyGenotypes.genotype(t,site) == obsgeno) obscount++;
			}
			pval = obscount / dtotal;
		}
		
		return pval;
	}
		
	private double probHetAsHom(int site) {
		int totalDepth = progenyGenotypes.depth().depth(myTaxonIndex, site);
		if (totalDepth == 1) return 0.99;
		return Math.pow(0.5, totalDepth - 1);
	}
	
}
