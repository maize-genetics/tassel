package net.maizegenetics.analysis.imputation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.data.FileLoadPlugin.TasselFileType;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;


public class RephaseParents {
	private static Logger myLogger = LogManager.getLogger(RephaseParents.class);
	private static final byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
	private static final byte N = GenotypeTable.UNKNOWN_ALLELE;
	private static final byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
	private static final byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
	private static final byte GG = NucleotideAlignmentConstants.getNucleotideDiploidByte("GG");
	private static final byte TT = NucleotideAlignmentConstants.getNucleotideDiploidByte("TT");
	private static final byte missingState = (byte) 4;
	
	GenotypeTable origGeno;
	Map<String, byte[]> progenyStates;
	List<String[]> plotList;
//	Map<String, String[]> plotMap;
	Map<String, byte[][]> rephasedParents;
	Map<String, byte[][]> startingParents = null;
	Map<String, double[][]> parentHaplotypeProbabilities = null; //probability that the haplotype equals the major allele
	
	//parameters
	int minFamilySize = 10;
	int minDepth = 7;
	String outputFilename;
	
	
	
	public RephaseParents() {
		
	}
	
	/**
	 * @param originalGenotypes		the GenotypeTable containing the original data being imputed
	 * @param phasedProgeny			the previously imputed progeny states (0 - 4, 4 = missing)
	 * @param plotList				the parentage information
	 * @param parentHapmap			the parent haplotypes assigned in the previous round
	 */
	public RephaseParents(GenotypeTable originalGenotypes, Map<String, byte[]> phasedProgeny, List<String[]> plotList, Map<String, byte[][]> parentHapmap) {
		origGeno = originalGenotypes;
		progenyStates = phasedProgeny;
		this.plotList = plotList;
		startingParents = parentHapmap;
		
	}
	
	public RephaseParents(GenotypeTable originalGenotypes, String phasedProgeny, String parentage, String parentHaps) {
		
		origGeno = originalGenotypes;
		GenotypeTable progenyStatesTable = (GenotypeTable) FileLoadPlugin.runPlugin(phasedProgeny);
		progenyStates = progenyStates(progenyStatesTable);
		myLogger.info(String.format("progeny states loaded: %s", phasedProgeny));
		
		try {
			plotList = Files.lines(Paths.get(parentage)).skip(1).map(in -> in.split("\t")).collect(Collectors.toList());
			myLogger.info(String.format("plotList has %d entries", plotList.size()));
		} catch (IOException e) {
			throw new RuntimeException("Unable to read " + parentage, e);
		}

		startingParents = ImputationUtils.restorePhasedHaplotypes(Paths.get(parentHaps));
		myLogger.info(String.format("Starting parent haplotypes loaded: %s", parentHaps));
	}

	Map<String, double[][]> rephaseUsingAlleleDepth() {
		return rephaseUsingAlleleDepth(null);
	}
	
	Map<String, double[][]> rephaseUsingAlleleDepth(String saveFilename) {
		Map<String, double[][]> nextHaplotypeProbs = new HashMap<>();

		//can phase each parent using all progeny with a phased other parent
		//for each parent create a list of all plots for that parent for which the other parent is phased
		Map<String, List<String[]>> parentPlotMap = new HashMap<>();
		for (String[] plot : plotList) {

			if (startingParents.get(plot[2]) != null) {
				List<String[]> parentPlotList = parentPlotMap.get(plot[1]);
				if (parentPlotList == null) {
					parentPlotList = new ArrayList<>();
					parentPlotMap.put(plot[1], parentPlotList);
				}
				parentPlotList.add(plot);
			}
			
			if (!plot[2].equals(plot[1]) && startingParents.get(plot[1]) != null) {
				List<String[]> parentPlotList = parentPlotMap.get(plot[2]);
				if (parentPlotList == null) {
					parentPlotList = new ArrayList<>();
					parentPlotMap.put(plot[2], parentPlotList);
				}
				parentPlotList.add(plot);
			}

		}
		
		//for each parent haplotype, at each site, calculate the probability that the haplotype is the major allele
		for (String parent : parentPlotMap.keySet()) {
			if (startingParents.get(parent) != null) {
				double[][] probs = rephasePreviouslyPhased(parent, parentPlotMap.get(parent));
				nextHaplotypeProbs.put(parent, probs);
			}
		}
		
		parentHaplotypeProbabilities = nextHaplotypeProbs;
		
		//save the rephased haplotypes to a file
//		String savename = "/Users/pbradbury/Documents/projects/teosinte/haplotypes/C2_rephased_parents.bin";
		if (saveFilename != null && saveFilename.length() > 1) {
			try {
				FileOutputStream fos = new FileOutputStream(new File(saveFilename));
				ObjectOutputStream oos = new ObjectOutputStream(fos);
				oos.writeObject(parentHaplotypeProbabilities);
				oos.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		return parentHaplotypeProbabilities;
	}
	
	double[][] rephasePreviouslyPhased(String parent, List<String[]> plotList) {							
		//calculate P(hap=major) for each haplotype
		//P(hap=a | obs) = P(hap=a | obs,other=a)*P(other=a) + P(hap=a | obs,other=b)*P(other=b)
		//	~ P(obs | hap=a, other=a)*P(other=a) + P(obs | hap=a, other=b)*P(other=b)
		
		//for each site
		//P(hap=a)=P(hap=b)=0.5
		//for P(other=a), P(other=b) = allele freq of a and b in total population
		//for each parent haplotype
		//group progeny into those with other parent carrying the major allele (a) vs. the minor allele (b)
		//add up the allele depths for each group. These are the obs.
		//calculate P(obs|...) from the binomial distribution
		//calculate P(hap=major) and store the result
		
		double err = 0.01;
		int nsites = origGeno.numberOfSites();
		double[][] haplotypeProbability = new double[2][nsites];
		for (int i = 0; i < 2; i++) Arrays.fill(haplotypeProbability[i], Double.NaN);
		
		for (int s = 0; s < nsites; s++) {
			byte major = origGeno.majorAllele(s);
			byte minor = origGeno.minorAllele(s);
			if (minor == N) {
				for (int i = 0; i < 2; i++) haplotypeProbability[i][s] = 1;
			} else {
				double majorFreq = origGeno.majorAlleleFrequency(s);
				double minorFreq = origGeno.minorAlleleFrequency(s);

				//sum allele depths across all plots and by progeny state
				int[] totalAlleleDepths = new int[6];
				int[][] stateAlleleDepths = new int[4][6];
				int[][] parentStateAlleleDepths = new int[2][6];
				int[] switchState = new int[]{0,2,1,3};
				for (String[] plot : plotList) {
					int taxonNdx = origGeno.taxa().indexOf(plot[0]);
					int[] tempDepths = origGeno.depthForAlleles(taxonNdx, s);
					byte geno = origGeno.genotype(taxonNdx, s);
					int state = progenyStates.get(plot[0])[s];
					if (state < 4) {
						int parentState;
						int myState;
						if (parent.equals(plot[1])) {
							if (state == 0 || state == 1) parentState = 0;
							else parentState = 1;
							myState = state;
						} else {
							if (state == 0 || state == 2) parentState = 0;
							else parentState = 1;
							myState = switchState[state];
						}
						for (int i = 0; i < 6; i++) {
							totalAlleleDepths[i] += tempDepths[i];
							stateAlleleDepths[myState][i] += tempDepths[i];
							parentStateAlleleDepths[parentState][i] += tempDepths[i];
						}

					}
				}
				
				if (parentStateAlleleDepths[0][major] < 5 || parentStateAlleleDepths[0][major] < 5) {
					//skip, too little data
				} else {
					//for each parent haplotype
					//P(hap=a | obs) = P(obs | hap=a, other=a)*P(hap=a)*P(other=a) + P(obs | hap=a, other=b)*P(hap=a)*P(other=b)
					//then divide by sum of P(hap=a|obs) + P(hap=b|obs)
					//use origGeno allele freq for p(other=a) and p(hap=a)
					//consider using something other than allele frequency for this
					//perhaps using previously computed p's
					
					//parent haplotype 0, pMajor. Only consider the observations for h00 genotypes. h01 does not need to be considered
					//P(obs | h00=major) = [P(obs | h00=major, h10=major, h11=major) * P(h10=major) *P(h11=major)
					//                   + P(obs | h00=major, h10=major, h11=minor) * P(h10=major) * P(h11=minor)
					//                   + P(obs | h00=major, h10=minor, h11=major) * P(h10=minor) * P(h11=major)
					//                   + P(obs | h00=major, h10=minor, h11=minor) * P(h10=minor) * P(h11=minor))] * P(h00=major)

					//for h00 = major
					double ph10major = new BinomialDistribution(stateAlleleDepths[0][major] + stateAlleleDepths[0][minor], err).probability(stateAlleleDepths[0][minor]);  //h10 = major
					double ph10minor = new BinomialDistribution(stateAlleleDepths[0][major] + stateAlleleDepths[0][minor], 0.5).probability(stateAlleleDepths[0][minor]);  //h10 = minor
					double ph11major = new BinomialDistribution(stateAlleleDepths[1][major] + stateAlleleDepths[1][minor], err).probability(stateAlleleDepths[1][minor]);  //h11 = major
					double ph11minor = new BinomialDistribution(stateAlleleDepths[1][major] + stateAlleleDepths[1][minor], 0.5).probability(stateAlleleDepths[1][minor]);  //h11 = minor
					double ph00major = ph10major * ph11major * majorFreq * majorFreq 
							+ ph10major * ph11minor * majorFreq * minorFreq 
							+ ph10minor * ph11major * minorFreq * majorFreq 
							+ ph10minor * ph11minor * minorFreq * minorFreq;
					ph00major *= majorFreq;
					
					//for h00 = minor
					ph10major = new BinomialDistribution(stateAlleleDepths[0][major] + stateAlleleDepths[0][minor], 0.5).probability(stateAlleleDepths[0][major]);  //h10 = major
					ph10minor = new BinomialDistribution(stateAlleleDepths[0][major] + stateAlleleDepths[0][minor], err).probability(stateAlleleDepths[0][major]);  //h10 = minor
					ph11major = new BinomialDistribution(stateAlleleDepths[1][major] + stateAlleleDepths[1][minor], 0.5).probability(stateAlleleDepths[1][major]);  //h11 = major
					ph11minor = new BinomialDistribution(stateAlleleDepths[1][major] + stateAlleleDepths[1][minor], err).probability(stateAlleleDepths[1][major]);  //h11 = minor
					double ph00minor = ph10major * ph11major * majorFreq * majorFreq 
							+ ph10major * ph11minor * majorFreq * minorFreq 
							+ ph10minor * ph11major * minorFreq * majorFreq 
							+ ph10minor * ph11minor * minorFreq * minorFreq;
					ph00minor *= minorFreq;
					
					haplotypeProbability[0][s] = ph00major / (ph00major + ph00minor);
					
					//for h01 = major
					ph10major = new BinomialDistribution(stateAlleleDepths[2][major] + stateAlleleDepths[2][minor], err).probability(stateAlleleDepths[2][minor]);  //h10 = major
					ph10minor = new BinomialDistribution(stateAlleleDepths[2][major] + stateAlleleDepths[2][minor], 0.5).probability(stateAlleleDepths[2][minor]);  //h10 = minor
					ph11major = new BinomialDistribution(stateAlleleDepths[3][major] + stateAlleleDepths[3][minor], err).probability(stateAlleleDepths[3][minor]);  //h11 = major
					ph11minor = new BinomialDistribution(stateAlleleDepths[3][major] + stateAlleleDepths[3][minor], 0.5).probability(stateAlleleDepths[3][minor]);  //h11 = minor
					double ph01major = ph10major * ph11major * majorFreq * majorFreq 
							+ ph10major * ph11minor * majorFreq * minorFreq 
							+ ph10minor * ph11major * minorFreq * majorFreq 
							+ ph10minor * ph11minor * minorFreq * minorFreq;
					ph01major *= majorFreq;
					
					//for h01 = minor
					ph10major = new BinomialDistribution(stateAlleleDepths[2][major] + stateAlleleDepths[2][minor], 0.5).probability(stateAlleleDepths[2][major]);  //h10 = major
					ph10minor = new BinomialDistribution(stateAlleleDepths[2][major] + stateAlleleDepths[2][minor], err).probability(stateAlleleDepths[2][major]);  //h10 = minor
					ph11major = new BinomialDistribution(stateAlleleDepths[3][major] + stateAlleleDepths[3][minor], 0.5).probability(stateAlleleDepths[3][major]);  //h11 = major
					ph11minor = new BinomialDistribution(stateAlleleDepths[3][major] + stateAlleleDepths[3][minor], err).probability(stateAlleleDepths[3][major]);  //h11 = minor
					double ph01minor = ph10major * ph11major * majorFreq * majorFreq 
							+ ph10major * ph11minor * majorFreq * minorFreq 
							+ ph10minor * ph11major * minorFreq * majorFreq 
							+ ph10minor * ph11minor * minorFreq * minorFreq;
					ph01minor *= minorFreq;
					
					haplotypeProbability[1][s] = ph01major / (ph01major + ph01minor);

				}
				
			}
		}
		return haplotypeProbability;
	}
	
	Map<String, byte[][]> rephaseUsingCrossProgeny() {
		
		rephasedParents = new HashMap<>();
		int[] firstParentChr = new int[]{0,0,1,1};
		int[] secondParentChr = new int[]{0,1,0,1};
		
		//for each parent of an outcross, get list of outcross progeny plots, store in parentPlotMap
		Map<String, List<String[]>> parentPlotMap = new HashMap<>();
		for (String[] plot : plotList) {
			if (plot[3].equals("outcross")) {
				List<String[]> parentPlotList = parentPlotMap.get(plot[1]);
				if (parentPlotList == null) {
					parentPlotList = new ArrayList<>();
					parentPlotMap.put(plot[1], parentPlotList);
				}
				parentPlotList.add(plot);
				if (!plot[2].equals(plot[1])) {
					parentPlotList = parentPlotMap.get(plot[2]);
					if (parentPlotList == null) {
						parentPlotList = new ArrayList<>();
						parentPlotMap.put(plot[2], parentPlotList);
					}
					parentPlotList.add(plot);
				}
			}
		}
		
		int nsites = origGeno.numberOfSites();
		
		//iterate through parents
		for (String parent : parentPlotMap.keySet()) {
			//debug
			System.out.printf("Rephasing %s\n", parent);
			
			//haplotypeList contains a list of haplotypes inferred for each of the progeny of this parent
			List<byte[][]> haplotypeList = new ArrayList<>();
			List<String[]> parentPlotList = parentPlotMap.get(parent);
			if (parentPlotList == null) {
				System.out.printf("parentPlotList null for %s\n", parent);
				continue;
			}
			if (parentPlotList.size() < minFamilySize) {
				System.out.printf("parentPlotList has %d plots for %s\n", parentPlotList.size(), parent);
				continue;
			}
			
			//for each progeny (plot[0] is the progeny name)
			for (String[] plot : parentPlotList) {
				//is parent the plot[1] or plot[2]?
				String otherParent;
				boolean isFirstParent;
				if (plot[1].equals(parent)) {
					otherParent = plot[2];
					isFirstParent = true;
				}
				else {
					otherParent = plot[1];
					isFirstParent = false;
				}
				
				//use all other parents, it may help to use only self parents but whatever
				int progenyIndex = origGeno.taxa().indexOf(plot[0]);
				byte[][] parentHap = new byte[2][nsites];
				for (int i = 0; i < 2; i++) Arrays.fill(parentHap[i], N);
				
				byte[] myStates = progenyStates.get(plot[0]);
				byte[][] otherParentHap = startingParents.get(otherParent);
				if (otherParentHap == null || myStates == null) continue;
				
				for (int s = 0;s < nsites; s++) {
					if (myStates[s] == missingState) continue;
					
					byte myGenotype = origGeno.genotype(progenyIndex, s);
					if (myGenotype == NN) continue;
					int myChr, otherChr;
					if (isFirstParent) {
						myChr = firstParentChr[myStates[s]];
						otherChr = secondParentChr[myStates[s]];
					}
					else {
						otherChr = firstParentChr[myStates[s]];
						myChr = secondParentChr[myStates[s]];
					}
					byte otherAllele = otherParentHap[otherChr][s];
					if (otherAllele == N) continue;
					int mydepth = origGeno.depth().depth(progenyIndex, s);
					if (GenotypeTableUtils.isHeterozygous(myGenotype)) {
						byte[] alleles = GenotypeTableUtils.getDiploidValues(myGenotype);
						if (otherAllele == alleles[0])  parentHap[myChr][s] = alleles[1];
						if (otherAllele == alleles[1]) parentHap[myChr][s] = alleles[0];
					} else { //progeny has homozygous genotype
						byte[] alleles = GenotypeTableUtils.getDiploidValues(myGenotype);
						if (alleles[0] != otherAllele && mydepth < 9) parentHap[myChr][s] = alleles[0];
						else if (mydepth >= minDepth && otherAllele == alleles[0]) parentHap[myChr][s] = alleles[0];
					} 
				}
				
				haplotypeList.add(parentHap);
			}
			
			//get majority allele at each site and major allele frequency
			byte[][] newhaps = new byte[2][nsites];
			Arrays.fill(newhaps[0], N);
			Arrays.fill(newhaps[1], N);

			for (int s = 0; s < nsites; s++) {
				byte major = origGeno.majorAllele(s);
				byte minor = origGeno.minorAllele(s);
				int[][] alleleCount = new int[2][6];
				for (byte[][] hap : haplotypeList) {
					for (int i = 0; i < 2; i++) {
						byte val = hap[i][s];
						if (val < 6) alleleCount[i][val]++;
					}
				}
				
				for (int i = 0; i < 2; i++) {
					int[] order = countSortOrder(alleleCount[i], true);
					
					if (alleleCount[i][order[0]] > 2 * alleleCount[i][order[1]]) {
						byte myAllele = (byte) order[0];
						if (myAllele == major || myAllele == minor) newhaps[i][s] = myAllele;
					}
				}
				
			}
			
			rephasedParents.put(parent, newhaps);
		}
		
		return rephasedParents;
	}

	public void setMinDepth(int mindepth) {
		minDepth = mindepth;
	}
	
	public static int[] countSortOrder(int[] counts, boolean descending) {
		int n = counts.length;
		List<Integer> order = IntStream.range(0, n).boxed().collect(Collectors.toList());
		
		if (descending) {
			Collections.sort(order, (a,b) -> {
				if (counts[a] > counts[b]) return -1;
				if (counts[a] < counts[b]) return 1;
				return 0;
			});
		} else {
			Collections.sort(order, (a,b) -> {
				if (counts[a] > counts[b]) return 1;
				if (counts[a] < counts[b]) return -1;
				return 0;
			});
		}
		return order.stream().mapToInt(I -> I.intValue()).toArray();
	}

	public static Map<String, byte[]> progenyStates(GenotypeTable gt) {
		Map<String, byte[]> outputMap = new HashMap<>();
		int nsites = gt.numberOfSites();
		int ntaxa = gt.numberOfTaxa();
		for (int t = 0; t < ntaxa; t++) {
			byte[] states = new byte[nsites];
			for (int s = 0; s < nsites; s++) {
				byte val = GenotypeTableUtils.getDiploidValues(gt.genotype(t, s))[0];
				if (val > -1 && val < 3) states[s] = val;
				else states[s] = missingState;
			}
			outputMap.put(gt.taxaName(t), states);
		}
		return outputMap;
	}

}
