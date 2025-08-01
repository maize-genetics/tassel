package net.maizegenetics.analysis.imputation;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.util.*;
import java.util.stream.Collectors;

import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.data.FileLoadPlugin.TasselFileType;
import net.maizegenetics.analysis.imputation.EmissionProbability;
import net.maizegenetics.analysis.imputation.TransitionProbability;
import net.maizegenetics.analysis.imputation.ViterbiAlgorithm;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.Tuple;

public class ImputeCrossProgeny {
	private Logger myLogger = LogManager.getLogger(ImputeCrossProgeny.class);
	
	private static final byte N = GenotypeTable.UNKNOWN_ALLELE;
	private static final byte missingState = 4;
	private static final byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
	private static final byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
	private static final byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("C");
	private static final byte GG = NucleotideAlignmentConstants.getNucleotideDiploidByte("G");
	private static final byte TT = NucleotideAlignmentConstants.getNucleotideDiploidByte("T");
	private static final byte[] stateNucleotide = new byte[]{AA, CC, GG, TT,NN};
	private static final byte[] ACGT = new byte[]{AA, CC, GG, TT};

	private List<String[]> plotList;
	private Map<String, byte[][]> haplotypeMap;
	private GenotypeTable myGenotype;
	
	//inputs
	private String parentCallFilename;
	
	//outputs
	private String phasedParentOutFilename;
	private String imputedGenotypeOutFilename;
	private String parentcallOutFilename;
	
	
	public void setPlotList(List<String[]> plotList) {
		this.plotList = plotList;
	}

	public void imputeAll() {
		String now = DateFormat.getDateInstance().format(new Date().getTime());

		int[] chrstart = myGenotype.chromosomesOffsets();
		int numberOfChromosomes = chrstart.length;
		int[] chrend = new int [10];
		System.arraycopy(chrstart, 1, chrend, 0, 9);
		chrend[9] = myGenotype.numberOfSites();

		Map<String, byte[]> phasedProgeny = new HashMap<>();
		Map<String, String[]> parentMap = new HashMap<>();
		
		//impute all progeny in that are in Genotype
		for (String[] plot : plotList) {
			int ndx = myGenotype.taxa().indexOf(plot[0]);
			if (ndx == -1) continue;
			byte[][] hap0 = haplotypeMap.get(plot[1]);
			byte[][] hap1 = haplotypeMap.get(plot[2]);
			if (hap0 != null && hap1 !=null && notMissingHap(hap0) && notMissingHap(hap1)) {
				phasedProgeny.put(plot[0], imputeCrossFromParents(plot[0], hap0, hap1));
				parentMap.put(plot[0], new String[]{plot[1], plot[2]});
			}
		}

		//fill gaps
		for (byte[] states : phasedProgeny.values()) fillgaps(states);
		
		//store the results in genotype tables
		//A,C,M for the parent haplotypes
		//nucleotides for the imputed data
		
		List<Taxon> taxa = phasedProgeny.keySet().stream().map(p -> new Taxon(p)).collect(Collectors.toList());
		Collections.sort(taxa);
		GenotypeTableBuilder parentCalls = GenotypeTableBuilder.getTaxaIncremental(myGenotype.positions());
		int nsites = myGenotype.numberOfSites();
		taxa.stream().forEach(t -> {
			byte[] states = phasedProgeny.get(t);
			byte[] geno = new byte[nsites];
			for (int s = 0; s < nsites; s++) {
				geno[s] = stateNucleotide[states[s]];
			}
			parentCalls.addTaxon(t, geno);
		});
		
		ExportUtils.writeToHapmap(parentCalls.build(), parentcallOutFilename);
		
		//fill in monomorphic sites where missing to deal with uncalled IBD
		GenotypeTableBuilder imputedGeno = GenotypeTableBuilder.getTaxaIncremental(myGenotype.positions());;
		taxa.stream().forEach(t -> {
			byte[] states = phasedProgeny.get(t);
			byte[] geno = new byte[nsites];
			String[] parents = parentMap.get(t.getName());
			byte[][] hap0 = haplotypeMap.get(parents[0]);
			byte[][] hap1 = haplotypeMap.get(parents[1]);
			for (int s = 0; s < nsites; s++) {
				if (hap0[0][s] == N || hap0[1][s] == N || hap1[0][s] == N || hap1[1][s] == N) {
					geno[s] = NN;
				} else if (states[s] == 0) {
					geno[s] = GenotypeTableUtils.getDiploidValue(hap0[0][s], hap1[0][s]);
				} else if (states[s] == 1) {
					geno[s] = GenotypeTableUtils.getDiploidValue(hap0[0][s], hap1[1][s]);
				} else if (states[s] == 2) {
					geno[s] = GenotypeTableUtils.getDiploidValue(hap0[1][s], hap1[0][s]);
				} else if (states[s] == 3) {
					geno[s] = GenotypeTableUtils.getDiploidValue(hap0[1][s], hap1[1][s]);
				} else {
					if (hap0[0][s] == hap0[1][s] && hap0[0][s] == hap1[0][s] && hap0[0][s] == hap1[1][s]) 
						geno[s] = GenotypeTableUtils.getDiploidValue(hap0[0][s], hap0[0][s]);
					else geno[s] = NN;
				}
				
			}
			imputedGeno.addTaxon(t, geno);
		});
		
		ExportUtils.writeToHapmap(imputedGeno.build(), imputedGenotypeOutFilename);

	}

	boolean notMissingHap(byte[][] hap) {
		int[] nm = new int[2];
		for (int i = 0; i < 2; i++) {
			for (byte b : hap[i]) if (b != N) nm[i]++;
		}
		if (nm[0] > 100 && nm[1] > 0) return true;
		return false;
	}

	public byte[] imputeCrossFromParents(String progeny, byte[][] hap0, byte[][] hap1) {
		//impute by chromosome
		int taxonIndex = myGenotype.taxa().indexOf(progeny);
		int nsites = myGenotype.numberOfSites();
		byte[] progenyGenotype = new byte[nsites];
		Arrays.fill(progenyGenotype, missingState);
		
		int[] chrstart = myGenotype.chromosomesOffsets();
		int numberOfChromosomes = chrstart.length;
		int[] chrend = new int[numberOfChromosomes];
		System.arraycopy(chrstart, 1, chrend, 0, numberOfChromosomes - 1);
		chrend[numberOfChromosomes - 1] = myGenotype.numberOfSites();
		
		for (int c = 0; c < numberOfChromosomes; c++) {
			System.out.printf("Imputing cross progeny %s, chr %d\n", progeny, c+1);
			//filter
			Chromosome chr = new Chromosome(Integer.toString(c + 1));
			
			//create the TransitionProbability, which will be the same for all taxa
			TransitionProbability tp = new TransitionProbability();
			
			//get the taxon genotype
			byte[] geno = myGenotype.genotypeRange(taxonIndex, chrstart[c], chrend[c]);

			//remove missing values
			//remove sites that are invariant in parents
			//consider removing values that are with 64 bp of each other
			int ngeno = geno.length;
			OpenBitSet isMissing = new OpenBitSet(ngeno);
			byte[] nonMissingGenotypes = new byte[ngeno];
			int[] nonMissingPositions = new int[ngeno];
			int[] nonMissingHaplotypeIndices = new int[ngeno];
			int numberNotMissing = 0;

			int prevpos = -100;
			int numberInconsistent = 0;
			for (int s = 0; s < ngeno; s++) {
				int siteIndex = chrstart[c] + s;
				int pos = myGenotype.positions().get(siteIndex).getPosition();
				int originalSite = siteIndex;
				if (geno[s] == NN || hap0[0][originalSite] == N || hap0[1][originalSite] == N ||hap1[0][originalSite] == N ||hap1[1][originalSite] == N || pos - prevpos < 64) {
					isMissing.fastSet(s);
				} else if (hap0[0][originalSite] == hap0[1][originalSite] && hap0[0][originalSite] == hap1[0][originalSite] && hap0[0][originalSite] == hap1[1][originalSite]) {
					isMissing.fastSet(s);
					if (geno[s] != GenotypeTableUtils.getDiploidValue(hap0[0][originalSite], hap0[0][originalSite])) numberInconsistent++;
//				} else if (hap0[0][originalSite] == hap0[1][originalSite] && hap1[0][originalSite] == hap1[1][originalSite] && hap0[0][originalSite] != hap1[1][originalSite]) {
				} else if (hap0[0][originalSite] == hap0[1][originalSite] && hap1[0][originalSite] == hap1[1][originalSite]) {
					isMissing.fastSet(s);  //this site is  not informative
				} else {
					nonMissingPositions[numberNotMissing] = pos;
					nonMissingGenotypes[numberNotMissing] = geno[s];
					nonMissingHaplotypeIndices[numberNotMissing] = originalSite;
					prevpos = pos;
					numberNotMissing++;
			}
			}
			System.out.printf("For %s, %d of %d sites inconsistent with parents\n", 
					progeny, numberInconsistent, numberNotMissing);
			nonMissingGenotypes = Arrays.copyOf(nonMissingGenotypes, numberNotMissing);
			nonMissingPositions = Arrays.copyOf(nonMissingPositions, numberNotMissing);
			nonMissingHaplotypeIndices = Arrays.copyOf(nonMissingHaplotypeIndices, numberNotMissing);

			//debug
			if (numberNotMissing < 100) {
				System.out.printf("Number not missing = %d for %s, chromosome %s\n", numberNotMissing, progeny, chr.getName());
				continue;
			}

			//update transition prob with avg segment length and (average) transition probability
			tp.setPositions(nonMissingPositions);

			//initialize the transition matrix
			double avgTP = 1 / ((double) (numberNotMissing - 1));  //on average 1 xo per chromosome
			double avgDO = avgTP * avgTP;
			double[][] transition = new double[][] {
					{1 - avgTP - avgDO, 0.5*avgTP, 0.5*avgTP, avgDO},
					{0.5*avgTP, 1 - avgTP - avgDO, avgDO, 0.5*avgTP},
					{0.5*avgTP, avgDO, 1 - avgTP - avgDO, 0.5*avgTP},
					{avgDO, 0.5*avgTP, 0.5*avgTP, 1 - avgTP - avgDO}
			};
			tp.setTransitionProbability(transition);
			tp.setAverageSegmentLength((nonMissingPositions[numberNotMissing - 1] - nonMissingPositions[0]) / (numberNotMissing - 1));

			//create an EmissionProbability
			EmissionProbability ep = new CrossProgenyEmissionMatrix(new byte[][][]{hap0, hap1}, myGenotype, taxonIndex, nonMissingHaplotypeIndices);

			//run Viterbi which returns a byte[] representing the state at each position
			ViterbiAlgorithm va = new ViterbiAlgorithm(nonMissingGenotypes, tp, ep, new double[]{0.25, 0.25, 0.25, 0.25});
			va.calculate();

			//report number of crossovers per chromosome and taxon
			byte[] states = va.getMostProbableStateSequence();

			for (int i = 0; i < numberNotMissing; i++) {
				progenyGenotype[nonMissingHaplotypeIndices[i]] = states[i];
			}
			
		}
		
		return progenyGenotype;
	}

	public void improveImputedProgenyStates() {
		int[] chrstart = myGenotype.chromosomesOffsets();
		int numberOfChromosomes = chrstart.length;
		int[] chrend = new int[numberOfChromosomes];
		for (int i = 1; i < numberOfChromosomes; i++) chrend[i - 1] = chrstart[i];
		chrend[numberOfChromosomes - 1] = myGenotype.numberOfSites();
		
		FileLoadPlugin flp = new FileLoadPlugin(null, false);
		flp.setTheFileType(TasselFileType.Unknown);
		flp.setOpenFiles(new File[]{new File(parentCallFilename)});
		GenotypeTable progenyStateGeno = (GenotypeTable) flp.performFunction(null).getData(0).getData();
		int nsites = progenyStateGeno.numberOfSites();
		
		Map<Byte, Byte> stateMap = new HashMap<>();
		stateMap.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("A"), new Byte((byte) 0));
		stateMap.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("C"), new Byte((byte) 1));
		stateMap.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("G"), new Byte((byte) 2));
		stateMap.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("T"), new Byte((byte) 3));
		stateMap.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("NN"), new Byte((byte) 4));

		Map<String, String[]> plotInfoMap = new HashMap<>();
		for (String[] plot:plotList) {
			plotInfoMap.put(plot[0], plot);
		}

		Map<String, String[]> parentMap = new HashMap<>();
		long[][] stateObsCounts = new long[3][3];
		
		myLogger.info("Loading progeny states and getting stateByObservation counts.");
		long start = System.currentTimeMillis();
		List<Taxon> myProgenyTaxaList = progenyStateGeno.taxa();
		int ntaxa = myProgenyTaxaList.size();
		byte[][] states = new byte[ntaxa][nsites];
		for (int s = 0; s < nsites; s++) {
			for (int t = 0; t < ntaxa; t++) {
				states[t][s] = stateMap.get(progenyStateGeno.genotype(t, s));
			}
		}
		
		myLogger.info(String.format("Progeny states loaded in %d ms.", System.currentTimeMillis() - start));
		start = System.currentTimeMillis();
		Map<String, byte[]> phasedProgenyMap = new HashMap<>();
		for (int t = 0; t < ntaxa; t++) {
			String taxonName = myProgenyTaxaList.get(t).getName();
			phasedProgenyMap.put(taxonName, states[t]);
			String[] plot = plotInfoMap.get(taxonName);
			byte[][] hap0 = haplotypeMap.get(plot[1]);
			byte[][] hap1 = haplotypeMap.get(plot[2]);
			updateStateObsCounts(stateObsCounts, hap0, hap1, states[t], taxonName);
			parentMap.put(taxonName, new String[]{plot[1], plot[2]});
		}
		myLogger.info(String.format("stateObsCounts updated in %d ms.", System.currentTimeMillis() - start));
		
		//rephase the parents based on imputed progeny states
		start = System.currentTimeMillis();
		myLogger.info("Rephasing parents");
		RephaseParents rephaser = new RephaseParents(myGenotype, phasedProgenyMap, plotList, haplotypeMap);
		Map<String, double[][]> haplotypeProbabilities = rephaser.rephaseUsingAlleleDepth(phasedParentOutFilename);
		myLogger.info(String.format("Time to rephase parents = %d ms.", System.currentTimeMillis() - start));
		
		//calculate probability of an observation given the state
		//the states are homozygous major, heterozygous, homozygous minor
		myLogger.info("Calculating stateGivenObs probability");
		double[][] stateGivenObs = new double[3][3];
		for (int i = 0; i < 3; i++) {
			double rowsum = (double) Arrays.stream(stateObsCounts[i]).sum();
			for (int col = 0; col < 3; col++) stateGivenObs[i][col] = stateObsCounts[i][col] / rowsum;
		}
		
		myLogger.debug("stateGivenObs probabilities");
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) sb.append(String.format("%1.4f ", stateGivenObs[i][j]));
		myLogger.debug(sb.toString());
		
		//reimpute the progeny
		start = System.currentTimeMillis();
		myLogger.info("Rephasing progeny.");
		List<Tuple<String,byte[]>> resultList = plotList.stream()
				.filter(p -> haplotypeProbabilities.get(p[1]) != null && haplotypeProbabilities.get(p[2]) != null)
				.map(p -> {
			double[][] hapProb0 = haplotypeProbabilities.get(p[1]);
			double[][] hapProb1 = haplotypeProbabilities.get(p[2]);
			byte[] progenyStates = imputeCrossFromParentsUsingProbabilities(p[0], hapProb0, hapProb1, stateGivenObs);
			return new Tuple<String, byte[]>(p[0], progenyStates);
		}).collect(Collectors.toList());
		Map<String, byte[]> rephasedProgeny = new HashMap<>();
		for (Tuple<String,byte[]> result : resultList) rephasedProgeny.put(result.x, result.y);
		myLogger.info(String.format("Progeny rephase in %d ms.", System.currentTimeMillis() - start)); 
		
		
		//fill gaps
		for (byte[] state : rephasedProgeny.values()) fillgaps(state);

		//output the results
		//store the results in genotype tables
		//A,C,M for the parent haplotypes
		//nucleotides for the imputed data

		GenotypeTableBuilder parentCalls = GenotypeTableBuilder.getTaxaIncremental(myGenotype.positions());
		myProgenyTaxaList.stream().forEach(t -> {
			byte[] state = rephasedProgeny.get(t.getName());
			byte[] geno = new byte[nsites];
			for (int s = 0; s < nsites; s++) {
				geno[s] = stateNucleotide[state[s]];
			}
			parentCalls.addTaxon(t, geno);
		});
		
		ExportUtils.writeToHapmap(parentCalls.build(), parentcallOutFilename);
		
		GenotypeTableBuilder imputedGeno = GenotypeTableBuilder.getTaxaIncremental(myGenotype.positions());
		double[] limit = new double[]{0.001, 0.999};
		myProgenyTaxaList.stream().forEach(t -> {
			byte[] state = rephasedProgeny.get(t.getName());
			byte[] geno = new byte[nsites];
			Arrays.fill(geno, NN);
			String[] parents = parentMap.get(t.getName());
			double[][] hapProb0 = haplotypeProbabilities.get(parents[0]);
			double[][] hapProb1 = haplotypeProbabilities.get(parents[1]);
			for (int s = 0; s < nsites; s++) {
				double prob[] = new double[]{Double.NaN, Double.NaN};
				switch(state[s]) {
				case 0:
					prob[0] = hapProb0[0][s];
					prob[1] = hapProb1[0][s];
					break;
				case 1:
					prob[0] = hapProb0[0][s];
					prob[1] = hapProb1[1][s];
					break;
				case 2:
					prob[0] = hapProb0[1][s];
					prob[1] = hapProb1[0][s];
					break;
				case 3:
					prob[0] = hapProb0[1][s];
					prob[1] = hapProb1[1][s];
					break;
				}

				if (Double.isNaN(prob[0]) || Double.isNaN(prob[1])) {
					geno[s] = NN;
				} else if ((prob[0] > limit[0] && prob[0] < limit[1]) || (prob[1] > limit[0] && prob[1] < limit[1])) {
					geno[s] = NN;
				} else {
					byte[] haplotype = new byte[2];
					if (prob[0] > 0.5) haplotype[0] = myGenotype.majorAllele(s);
					else haplotype[0] = myGenotype.minorAllele(s);
					if (prob[1] > 0.5) haplotype[1] = myGenotype.majorAllele(s);
					else haplotype[1] = myGenotype.minorAllele(s);
					geno[s] = GenotypeTableUtils.getDiploidValue(haplotype[0], haplotype[1]);
				}
			}
			imputedGeno.addTaxon(t, geno);
		});
		
		ExportUtils.writeToHapmap(imputedGeno.build(), imputedGenotypeOutFilename);
		
	}

	public byte[] imputeCrossFromParentsUsingProbabilities(String progeny, double[][] hapProb0, double[][] hapProb1, double[][] probObsGivenState) {
		int[] chrstart = myGenotype.chromosomesOffsets();
		int numberOfChromosomes = chrstart.length;
		int[] chrend = new int[numberOfChromosomes];
		System.arraycopy(chrstart, 1, chrend, 0, numberOfChromosomes - 1);
		chrend[numberOfChromosomes - 1] = myGenotype.numberOfSites();

		//debug
		if (hapProb0 == null) myLogger.debug(String.format("hapProb0 is null for progeny %s", progeny));
		if (hapProb1 == null) myLogger.debug(String.format("hapProb1 is null for progeny %s", progeny));
		
		//hapProb is the probability, for each chromosome, each parent, each position that the chromosome carries the major allele
		//impute by chromosome
		int nsites = myGenotype.numberOfSites();
		byte[] progenyGenotype = new byte[nsites];
		Arrays.fill(progenyGenotype, missingState);
		
		for (int c = 0; c < 10; c++) {
//			System.out.printf("Imputing cross progeny %s, chr %d\n", progeny, c+1);
			//filter
			Chromosome chr = new Chromosome(Integer.toString(c + 1));
			
			//create the TransitionProbability, which will be the same for all taxa
			TransitionProbability tp = new TransitionProbability();
			
			//get the taxon genotype
			int taxonIndex = myGenotype.taxa().indexOf(progeny);
			byte[] geno = myGenotype.genotypeRange(taxonIndex, chrstart[c], chrend[c]);

			
			//remove missing values
			//remove sites that are invariant in parents
			//consider removing values that are with 64 bp of each other
			int ngeno = geno.length;
			OpenBitSet isMissing = new OpenBitSet(ngeno);
			byte[] nonMissingGenotypes = new byte[ngeno];
			int[] nonMissingPositions = new int[ngeno];
			int[] nonMissingHaplotypeIndices = new int[ngeno];
			int numberNotMissing = 0;

			int prevpos = -100;
			int numberInconsistent = 0;
			int nNotNN = 0;
			int[] nNotNan = new int[4];
			for (int s = 0; s < ngeno; s++) {
				int siteIndex = chrstart[c] + s;
				int pos = myGenotype.positions().get(siteIndex).getPosition();
				int originalSite = siteIndex;
				
				if (geno[s] == NN || Double.isNaN(hapProb0[0][originalSite]) || Double.isNaN(hapProb0[1][originalSite]) || Double.isNaN(hapProb1[0][originalSite]) || Double.isNaN(hapProb1[1][originalSite]) || pos - prevpos < 64) {
					isMissing.fastSet(s);
				} else {
					nonMissingPositions[numberNotMissing] = pos;
					nonMissingGenotypes[numberNotMissing] = geno[s];
					nonMissingHaplotypeIndices[numberNotMissing] = originalSite;
					prevpos = pos;
					numberNotMissing++;
				}
			}
			
			if (numberNotMissing < 20) {
				myLogger.debug(String.format("numberNotMissing = %d for %s chr %d", numberNotMissing, progeny, c+1));
				myLogger.debug(String.format("notNN, notNaN(00,01,10,11): %d, %d, %d, %d, %d", nNotNN, nNotNan[0], nNotNan[1], nNotNan[2], nNotNan[3]));
				if (numberNotMissing < 5) continue;
			}
			
			nonMissingGenotypes = Arrays.copyOf(nonMissingGenotypes, numberNotMissing);
			nonMissingPositions = Arrays.copyOf(nonMissingPositions, numberNotMissing);
			nonMissingHaplotypeIndices = Arrays.copyOf(nonMissingHaplotypeIndices, numberNotMissing);

			//update transition prob with avg segment length and (average) transition probability
			tp.setPositions(nonMissingPositions);

			//initialize the transition matrix
			double avgTP = 1 / ((double) (numberNotMissing - 1));  //on average 1 xo per chromosome
			double avgDO = avgTP * avgTP;
			double[][] transition = new double[][] {
					{1 - avgTP - avgDO, 0.5*avgTP, 0.5*avgTP, avgDO},
					{0.5*avgTP, 1 - avgTP - avgDO, avgDO, 0.5*avgTP},
					{0.5*avgTP, avgDO, 1 - avgTP - avgDO, 0.5*avgTP},
					{avgDO, 0.5*avgTP, 0.5*avgTP, 1 - avgTP - avgDO}
			};
			tp.setTransitionProbability(transition);
			tp.setAverageSegmentLength((nonMissingPositions[numberNotMissing - 1] - nonMissingPositions[0]) / (numberNotMissing - 1));

			//create an EmissionProbability
			double[][] haplotypeProbs = new double[4][];
			haplotypeProbs[0] = hapProb0[0];
			haplotypeProbs[1] = hapProb0[1];
			haplotypeProbs[2] = hapProb1[0];
			haplotypeProbs[3] = hapProb1[1];
			
			//need double[][] stateByObs for emission probability
			//this should be calculated from progeny states and the original genotypes
			EmissionProbability ep = new CrossProgenyEmissionMatrix(haplotypeProbs, myGenotype, probObsGivenState, taxonIndex, nonMissingHaplotypeIndices);

			//run Viterbi which returns a byte[] representing the state at each position
			ViterbiAlgorithm va = new ViterbiAlgorithm(nonMissingGenotypes, tp, ep, new double[]{0.25, 0.25, 0.25, 0.25});
			va.calculate();

			//report number of crossovers per chromosome and taxon
			byte[] states = va.getMostProbableStateSequence();

//			int n = states.length;
//			int transitionCount = 0;
//			for (int i = 1; i < n; i++) {
//				if (states[i] != states[i - 1]) transitionCount++;
//			}

//			try (BufferedWriter bw = Files.newBufferedWriter(crossoverCountsPath, StandardOpenOption.CREATE,StandardOpenOption.APPEND)) {
//				bw.write(String.format("%s\t%d\t%d\n", progeny, c + 1, transitionCount));
//			} catch(IOException e) {
//				e.printStackTrace();
//			}
			
			for (int i = 0; i < numberNotMissing; i++) {
				progenyGenotype[nonMissingHaplotypeIndices[i]] = states[i];
			}
			
		}
		
		return progenyGenotype;
	}

	public void updateStateObsCounts(long[][] counts, byte[][] hap0, byte[][] hap1, byte[] states, String progenyName) {
		int ndx = myGenotype.taxa().indexOf(progenyName);
		int nsites = myGenotype.numberOfSites();
		int[][] stateMap = new int[][]{{0,0},{0,1},{1,0},{1,1}};
		
		//set state and obs index to
		//0 = homozygous major, 1 = het, 2 = homozygous minor
		for (int s = 0; s < nsites; s++) {
			int stateIndex = -1; 
			int obsIndex = -1;
			byte obs = myGenotype.genotype(ndx, s);
			if (obs == NN || states[s] == 4) continue; //missing data, no update
			
			byte major = myGenotype.majorAllele(s);
			byte minor = myGenotype.minorAllele(s);
			
			int[] myState = stateMap[states[s]];
			byte[] stateAlleles = new byte[]{hap0[myState[0]][s], hap1[myState[1]][s]};
			
			if (stateAlleles[0] == major) {
				if (stateAlleles[1] == major) stateIndex = 0;
				else if (stateAlleles[1] == minor) stateIndex = 1;
			} else if (stateAlleles[0] == minor) {
				if (stateAlleles[1] == major) stateIndex = 1;
				else if (stateAlleles[1] == minor) stateIndex = 2;
			}
			
			if (GenotypeTableUtils.isHeterozygous(obs)) obsIndex = 1;
			else {
				byte obsAllele = GenotypeTableUtils.getDiploidValues(obs)[0];
				if (obsAllele == major) obsIndex = 0;
				else if (obsAllele == minor) obsIndex = 2;
			}
			
			
			if (stateIndex > -1 && obsIndex > -1) {
				counts[stateIndex][obsIndex]++;
			}
		}
	}

	public void fillgaps(byte[] states) {
		int[] chrstart = myGenotype.chromosomesOffsets();
		int[] chrend = new int[10];
		System.arraycopy(chrstart, 1, chrend, 0, 9);
		chrend[9] = myGenotype.numberOfSites();
		for (int c = 0; c < 10; c++) {
			int segStart = chrstart[c];
			byte segNuke = states[segStart];
			for (int s = chrstart[c] + 1; s < chrend[c]; s++) {
				if (segNuke == missingState) {
					if (states[s] != missingState) {
						segStart = s;
						segNuke = states[s];
					}
				} else {
					if (states[s] == segNuke) {
						for (int i = segStart + 1; i < s; i++) states[i] = segNuke;
						segStart = s;
					} else if (states[s] != missingState) {
						segStart = s;
						segNuke = states[s];
					}
				}
			}
		}
	}

	public void setParentCallInputFilename(String parentCallFilename) {
		this.parentCallFilename = parentCallFilename;
	}

	public void setMyGenotype(GenotypeTable myGenotype) {
		this.myGenotype = myGenotype;
	}
	
	public void setParentage(String parentageFilename) {
    	plotList = new ArrayList<>();
		try (BufferedReader br = Files.newBufferedReader(Paths.get(parentageFilename))) {
			br.readLine();
			String input;
			while ((input = br.readLine()) != null) {
				String[] data = input.split("\t");
				if (data.length > 3) {
					plotList.add(data);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	public void setHaplotypeMap(String parentHaplotypeFilename) {
    	try {
    		FileInputStream fis = new FileInputStream(parentHaplotypeFilename);
            ObjectInputStream ois = new ObjectInputStream(fis);
            haplotypeMap = (Map<String, byte[][]>) ois.readObject();
            ois.close();
		} catch (IOException | ClassNotFoundException e) {
			haplotypeMap = null;
			e.printStackTrace();
		}
	}

	public void setImputedGenotypeOutFilename(String imputedGenotypeOutFilename) {
		this.imputedGenotypeOutFilename = imputedGenotypeOutFilename;
	}

	public void setParentcallOutFilename(String parentcallOutFilename) {
		this.parentcallOutFilename = parentcallOutFilename;
	}

	public void setPhasedParentOutFilename(String phasedParentOutFilename) {
		this.phasedParentOutFilename = phasedParentOutFilename;
	}

	public void writeBreakpointFile(GenotypeTable imputedGenotypes, String outFilename) {
		//call setParentage() to import the parentage file

		//get a list of parents
		//index the parent haplotypes

		//first get a set of parent haplotypes
		Set<String> haplotypeSet = new TreeSet<>();
		for (String[] record : plotList) {
			haplotypeSet.add(record[1] + "_0");
			haplotypeSet.add(record[1] + "_1");
			haplotypeSet.add(record[2] + "_0");
			haplotypeSet.add(record[2] + "_1");
		}

		//index the parent haplotypes
		Map<String, Integer> haplotypeMap = new HashMap<>();
		int ndx = 0;
		for (String haplotype : haplotypeSet) haplotypeMap.put(haplotype, ndx++);

		Map<String, int[]> progenyMap = new HashMap<>();
		for (String[] record : plotList) {
			int[] progenyParents = new int[4];
			progenyParents[0] = haplotypeMap.get(record[1] + "_0");
			progenyParents[1] = haplotypeMap.get(record[1] + "_1");
			progenyParents[2] = haplotypeMap.get(record[2] + "_0");
			progenyParents[3] = haplotypeMap.get(record[2] + "_1");
			progenyMap.put(record[0], progenyParents);
		}

		//write breakpoint file header
		//write index
		try (PrintWriter pwbp = new PrintWriter(outFilename)) {
			int numberOfParents = haplotypeMap.size();
			int numberOfProgeny = imputedGenotypes.numberOfTaxa();
			pwbp.printf("%d\t%d%n", numberOfParents, numberOfProgeny);
			pwbp.println("#indexed parents");
			for (Map.Entry<String, Integer> entry : haplotypeMap.entrySet()) pwbp.printf("%d\t%s%n", entry.getValue(), entry.getKey());
			pwbp.println("#progeny breakpoints");

			//for each of the progeny, create a list of breakpoints using the parental haplotype indices
			PositionList myPositions = imputedGenotypes.positions();
			List<String> bpRecords = new ArrayList<>();
			for (int taxonIndex = 0; taxonIndex < imputedGenotypes.numberOfTaxa(); taxonIndex++) {
				Position hapStartPos = null;
				Position hapEndPos = null;
				String taxonName = imputedGenotypes.taxaName(taxonIndex);
				StringBuilder bpRecordBuilder = new StringBuilder(taxonName);
				int[] myParents = progenyMap.get(taxonName);
				byte previousGeno = GenotypeTable.RARE_DIPLOID_ALLELE;
				for (int siteIndex = 0; siteIndex < imputedGenotypes.numberOfSites(); siteIndex++) {
					byte currentGeno = imputedGenotypes.genotype(taxonIndex, siteIndex);
					Position currentPos = imputedGenotypes.positions().get(siteIndex);

					//if hapStartPos = null and currentGeno in [ACGT]
					//   set hapStartPos = hapStartEnd = current, previousGeno = current
					//else if current chromosome != hapEndPos chromosome, output a breakpoint record
					//  and if currentGeno in [ACGT], set hapEndPos = hapStartPos = current and previousGeno= current
					//  else set hapStartPos = null, hapEndPos = null
					//else if currentGeno = previousGeno set hapEndPos = current
					//else if currentGeno in [ACGT]
					//   output breakpoint, set hapStartPos = hapEndPos = current, previousGeno = current
					//else do nothing

					if (hapStartPos == null && inACGT(currentGeno)) {
						hapStartPos = hapEndPos = currentPos;
						previousGeno = currentGeno;
					} else if (!currentPos.getChromosome().equals(hapEndPos.getChromosome())) {
						appendBreakpointRecord(bpRecordBuilder, hapStartPos, hapEndPos, myParents, previousGeno);
						if (inACGT(currentGeno)) {
							hapEndPos = hapStartPos = currentPos;
							previousGeno = currentGeno;
						} else hapEndPos = hapStartPos = null;
					} else if (currentGeno == previousGeno) {
						hapEndPos = currentPos;
					} else if (inACGT(currentGeno)) {
						appendBreakpointRecord(bpRecordBuilder, hapStartPos, hapEndPos, myParents, previousGeno);
						hapEndPos = hapStartPos = currentPos;
						previousGeno = currentGeno;
					}
				}
				//write to file
				pwbp.println(bpRecordBuilder.toString());
			}

		} catch(FileNotFoundException e ) {
			throw new IllegalArgumentException("Cannot open " + outFilename + " for writing.");
		}

	}

	private boolean inACGT(byte diploidGenotype) {
		for (byte geno : ACGT) if (diploidGenotype == geno) return true;
		return false;
	}

	private void appendBreakpointRecord(StringBuilder bp, Position start, Position end, int[] parents, byte geno) {
		if (!inACGT(geno) || start == null || end == null) return;
		bp.append(" ").append(start.getChromosome().getName()).append(":")
				.append(start.getPosition()).append(":").append(end.getPosition()).append(":");

//		from the parent calls:
//		A = mom, chr 0: dad, chr 0
//		C = mom, chr 0: dad, chr 1
//		G = mom, chr 1: dad, chr 0
//		T = mom, chr 1: dad, chr 1

		if (geno == AA) bp.append(parents[0]).append(":").append(parents[2]);
		if (geno == CC) bp.append(parents[0]).append(":").append(parents[3]);
		if (geno == GG) bp.append(parents[1]).append(":").append(parents[2]);
		if (geno == TT) bp.append(parents[1]).append(":").append(parents[3]);
	}



}
