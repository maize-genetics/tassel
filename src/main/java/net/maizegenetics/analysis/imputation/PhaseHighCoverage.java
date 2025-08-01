package net.maizegenetics.analysis.imputation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.lang3.SerializationUtils;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.stat.inference.ChiSquareTest;

import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.data.FileLoadPlugin.TasselFileType;
import net.maizegenetics.analysis.distance.MultiDimensionalScalingPlugin;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.stats.PCA.ClassicMds;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.tree.Tree;
import net.maizegenetics.taxa.tree.TreeClusters;
import net.maizegenetics.taxa.tree.UPGMATree;
import net.maizegenetics.util.LoggingUtils;

public class PhaseHighCoverage {
	//this class uses sites with depth of 5 to phase parents
	//can phase parents using one progeny but will not know where xo's are
	private static final byte N = GenotypeTable.UNKNOWN_ALLELE;
	private static final byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;

	private Path parentagePath;
	private Path genopath;
	
	private Path outhapsSelfCross;
	private Path outhapsCross;
	private Path monomorphs;
	private Path hapsSelf;
	private Path outhapsCrossHighCover;
	private Path outhapsAllProgeny;
	
	private GenotypeTable myGenotypeTable;

	private static ChiSquareTest chisqTest = new ChiSquareTest();
	private int monoMultiplier = 100;
	
	Map<String, byte[][]> selfHaps;

	public PhaseHighCoverage(GenotypeTable genotype) {
		myGenotypeTable = genotype;
	}
	
	public List<String[]> loadPlotInfo() {
    	List<String[]> plotList = new ArrayList<>();
    	List<String> taxaNames = myGenotypeTable.taxa().stream().map(t -> t.getName()).collect(Collectors.toList());
		try (BufferedReader br = Files.newBufferedReader(parentagePath)) {
			br.readLine();
			String input;
			while ((input = br.readLine()) != null) {
				String[] data = input.split("\t");
				if (data.length > 3) {
					//check to make sure that myGenotypeTable contains the plot name
					if (taxaNames.contains(data[0]))
						plotList.add(data);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return plotList;
	}

	public void phaseParentsUsingAllAvailableProgeny(double minEigenRatio, Path savepath) {
		outhapsCrossHighCover = savepath;
		phaseParentsUsingAllAvailableProgeny(minEigenRatio);
	}
	
	public void phaseParentsUsingAllAvailableProgeny(double minEigenRatio) {
		System.out.println("Phasing parents using method phaseParentsUsingAllAvailableProgeny().");
		System.out.println("That in turn uses phaseParentUsingSelfAndCrossProgeny()");
		System.out.println("-------------------------------------------------------");
		int minNumberPhasedSites = 1000;
		Map<String, byte[][]> phasedParents = new HashMap<>();
//		myGenotypeTable = loadC2TeoGenotype();
		
		List<String[]> plotInfo = loadPlotInfo();
		
		TreeSet<String> parentSet = new TreeSet<>();
		for (String[] plot : plotInfo) {
			parentSet.add(plot[1]);
			parentSet.add(plot[2]);
		}
		
		for (String parent : parentSet) {
			System.out.printf("Phasing %s\n", parent);
			byte[][] phase = phaseParentUsingSelfAndCrossProgeny(parent, myGenotypeTable, minEigenRatio, plotInfo);
			if (phase == null) {
				System.out.println("Too few phased haplotypes, skipping.");
				System.out.println();
			} else {
				int nsites = phase[0].length;
				int phasedSiteCount = 0;
				for (int s = 0; s < nsites; s++) {
					if (phase[0][s] != N) {
						phasedSiteCount++;
					}
					
				}
				System.out.printf("%d sites phased for %s\n", phasedSiteCount, parent);
				if (phasedSiteCount < minNumberPhasedSites) {
					System.out.println("Too few sites phased, skipping.");
				} else {
					phasedParents.put(parent, phase);
				}
			}

		}
		
		SelfedHaplotypeFinder.serializePhasedHaplotypes(phasedParents, outhapsCrossHighCover);
		System.out.println("Finished phasing and storing haplotypes.");

	}
	
	public byte[][] phaseParentUsingSelfAndCrossProgeny(String parent, GenotypeTable myGeno, double minEigenRatio, List<String[]> plotInfo) {
		//set the min chisquare statistic
		double alpha = 0.0001;
		double chisqLimit = new ChiSquaredDistribution(1).inverseCumulativeProbability(1 - alpha);
		
		List<byte[]> phasedHaplotypes = new ArrayList<>();
		int parentIndex = myGeno.taxa().indexOf(parent);
		byte[] parentGeno = myGeno.genotypeAllSites(parentIndex);
		
		byte[][] phasedParent = new byte[2][];
		for (int i = 0; i < 2; i++) {
			phasedParent[i] = new byte[myGeno.numberOfSites()];
			Arrays.fill(phasedParent[i], N);
		}
		
		for (String[] plot : plotInfo) {
			if (plot[1].equals(parent) || plot[2].equals(parent)) {
				String otherParent;
				if (plot[1].equals(parent)) otherParent = plot[2];
				else otherParent = plot[1];
				System.out.println("phasing " + plot[0]);
				byte[][] haps = phaseParentUsingOneProgeny(parent, otherParent, plot[0], myGeno);
				phasedHaplotypes.add(haps[0]);
				phasedHaplotypes.add(haps[1]);
			}
		}
		
		if (phasedHaplotypes.size() < 10) return null;

		//need to decide how to handle IBD segments (minimize xo's?)
		//initially, only polymorphic sites are being phased, which will eliminate IBD segments
		//adding monomorphic sites back in will capture those
		//have to do this by chromosome
		int[] chrstart = myGeno.chromosomesOffsets();
		int[] chrend = new int[10];
		System.arraycopy(chrstart, 1, chrend, 0, 9);
		chrend[9] = myGeno.numberOfSites();
		int window = 40; //number of sites with sufficient coverage & polymorphic
		int minWindow = 20;
		int minPresent = 4;
		List<int[]> monomorphicSites = new ArrayList<>();
		for (int c = 0; c < 10; c++) { 
			int s = chrstart[c];
			List<Integer> prevHapIndices1 = null;
			List<Integer> prevHapIndices2 = null;
			int[] previousSiteIndex = null;
			boolean isPreviousHapValid = false;
			while (s < chrend[c]) {
				//index of polymorphic sites
				int[] siteIndex = new int[window];
				int indexCount = 0;
				
				//this while segment finds the next window polymorphic loci
				//stores the sites in siteIndex. indexCount is the number found, which is < window only at the end of the chromosome
				int startS = s;
				while (s < chrend[c] && indexCount < window) {
					int[] alleleCounts = countAllelesAtSite(phasedHaplotypes, s);
					int npresent = Arrays.stream(alleleCounts).sum();
					if (npresent > minPresent) {
						
						//generate sorted Index
						List<Integer> index = IntStream.range(0, 6).boxed().collect(Collectors.toList());
						Collections.sort(index, (a,b) -> alleleCounts[a] >= alleleCounts[b] ? -1 : 1);
						
						//test for polymorphism (minor allele count * 4 > major allele count)
						int mult = 10; //was 4 for phase 1, but changing to 10 in phase2 teosinte
						if (alleleCounts[index.get(1)] > 1 && alleleCounts[index.get(1)] * mult > alleleCounts[index.get(0)]) {
							siteIndex[indexCount++] = s;
						}
						
						//test for monomorphism
						if (alleleCounts[index.get(0)] > 20 && alleleCounts[index.get(1)] == 0 ) {
							monomorphicSites.add(new int[]{s, index.get(0)});
						}
					}
					s++;
				}
				
				int nSitesScanned = s - startS;
				
				//if indexCount is too small, create a window by adding the new sites on the end of the previous window
				if (indexCount < minWindow) {
					if (!isPreviousHapValid) {
						//at end of chromosome and no valid window has been found. Probably, the entire chromosome is homozygous.
						System.out.printf("No windows phased in chromosome %d%n", c + 1);
						continue; 
					}
					int combinedCount = previousSiteIndex.length + indexCount;
					int[] combinedIndex = new int[combinedCount];
					System.arraycopy(previousSiteIndex, 0, combinedIndex, 0, previousSiteIndex.length);
					System.arraycopy(siteIndex, 0, combinedIndex, previousSiteIndex.length, indexCount);
					siteIndex = combinedIndex;
					indexCount = siteIndex.length;
				}
				
				//create a list of these segments
				List<byte[]> seglist = new ArrayList<>();
				for (byte[] haps : phasedHaplotypes) {
					byte[] seg = new byte[indexCount];
					for (int i = 0; i < indexCount; i++) seg[i] = haps[siteIndex[i]];
					seglist.add(seg);
				}
				
				//calculate pairwise distances as mismatch proportion
				double[][] dist = mismatchDistanceMatrix(seglist);
				List<Taxon> dummyTaxa = IntStream.range(0,  dist.length).mapToObj(i -> new Taxon(Integer.toString(i))).collect(Collectors.toList());
				DistanceMatrix dm = new DistanceMatrix(dist, new TaxaListBuilder().addAll(dummyTaxa).build());
				
				//test for single haplotype
				ClassicMds mds = new ClassicMds(dm);
				System.out.printf("Eigenvalue ratio = %1.3f, Eigenvalues: %1.3f, %1.3f, %1.3f, %1.3f\n", mds.getEigenvalue(0)/mds.getEigenvalue(1), mds.getEigenvalue(0), mds.getEigenvalue(1), mds.getEigenvalue(2), mds.getEigenvalue(3));
				System.out.printf("%d sites scanned to generate this interval\n", nSitesScanned);
				if (mds.getEigenvalue(0)/mds.getEigenvalue(1) < minEigenRatio) continue; //skip this interval
				
				UPGMATree myTree = new UPGMATree(dm);
				TreeClusters myClusters = new TreeClusters(myTree);
				int[] myGroups = myClusters.getGroups(2);
				
				//make two groups
				List<Integer> hapIndices1 = new ArrayList<>();
				List<Integer> hapIndices2 = new ArrayList<>();
				int n = myGroups.length;
				for (int i = 0; i < n; i++) {
					String name = myTree.getExternalNode(i).getIdentifier().getName();
					if (myGroups[i] == 0) hapIndices1.add(Integer.parseInt(name));
					else hapIndices2.add(Integer.parseInt(name));
				}
				
				//keep track of which haplotypes are assigned to which windows
				boolean reverseHaps;
				if (prevHapIndices1 == null) reverseHaps = false;
				else {
					int[][] matches = new int[2][2];
					matches[0][0] = countSharedMembers(hapIndices1, prevHapIndices1);
					matches[0][1] = countSharedMembers(hapIndices1, prevHapIndices2);
					matches[1][0] = countSharedMembers(hapIndices2, prevHapIndices1);
					matches[1][1] = countSharedMembers(hapIndices2, prevHapIndices2);
					int mainDiagSum = matches[0][0] + matches[1][1];
					int offDiagSum = matches[0][1] + matches[1][0];
					if (mainDiagSum > 2 * offDiagSum) {
						reverseHaps = false;
						isPreviousHapValid = true;
					}
					else if (offDiagSum > 2 * mainDiagSum) {
						reverseHaps = true;
						isPreviousHapValid = true;
					}
					else {
						//if the previous haplotype is the first one, it has not been validated (by making sure that it can be matched to the adjacent haplotype)
						//if the previous haplotype has not been validated, invalidate it (do not use it)
						if (!isPreviousHapValid) {
							//set previous Hap sites to missing
							for (int prevSite : previousSiteIndex) {
								phasedParent[0][prevSite] = N;
								phasedParent[1][prevSite] = N;
							}
							reverseHaps = false;
						} else {
							continue; //do not process this interval but go on to next
						}
						
					}
					System.out.printf("haplotype matches at chr %d, site %d: %d, %d, %d, %d, reverse = %b\n", c + 1, s, matches[0][0],matches[0][1],matches[1][0],matches[1][1], reverseHaps);
				}

				System.out.printf("hap list 1 has %d members, 2 has %d members (chr %d, site %d)\n", hapIndices1.size(), hapIndices2.size(), c + 1, s);
				List<byte[]> hapList1, hapList2;
				if (reverseHaps) {
					hapList1 = hapIndices2.stream().map(I -> seglist.get(I)).collect(Collectors.toList());
					hapList2 = hapIndices1.stream().map(I -> seglist.get(I)).collect(Collectors.toList());
				} else {
					hapList1 = hapIndices1.stream().map(I -> seglist.get(I)).collect(Collectors.toList());
					hapList2 = hapIndices2.stream().map(I -> seglist.get(I)).collect(Collectors.toList());
				}
				
				System.out.println(haplotypeAsString(consensusHaplotype(hapList1)));
				System.out.println(haplotypeAsString(consensusHaplotype(hapList2)));

				//update phasedHaplotypes with sites that segregate with haplotypes
				//update phasedHaplotypes[0] with hapList1, phasedHaplotypes[1] with hapList2
				for (int i = 0; i < indexCount; i++) {
					int[] alleleCount1 = countAllelesAtSite(hapList1, i);
					int[] alleleCount2 = countAllelesAtSite(hapList2, i);
					long[][] counts = new long[2][2];
					int which = 0;
					for (int j = 0; j < 4; j++) {
						if (alleleCount1[j] > 0 || alleleCount2[j] > 0) {
							if (which == 2) {
								System.out.printf("Site %d has more than 2 alleles\n", siteIndex[i]);
								break;
							}
							counts[0][which] = alleleCount1[j];
							counts[1][which] = alleleCount2[j];
							which++;
						}
					}
					
					double[] ratio = new double[]{counts[0][0] / (double) counts[0][1], counts[1][0] / (double) counts[1][1]};
					if ((ratio[0] >= 2 && ratio[1] <= 0.5) || (ratio[1] >= 2 && ratio[0] <= 0.5)) {
						double testval = chisqTest.chiSquare(counts);
						if (testval >= chisqLimit) { //update haplotypes with this site
							phasedParent[0][siteIndex[i]] = maxAllele(alleleCount1);
							phasedParent[1][siteIndex[i]] = maxAllele(alleleCount2);
						}
					}
					
				}
				
				if (reverseHaps) {
					prevHapIndices1 = hapIndices2;
					prevHapIndices2 = hapIndices1;
				} else {
					prevHapIndices1 = hapIndices1;
					prevHapIndices2 = hapIndices2;
				}
				previousSiteIndex = siteIndex;
			}
		}
		
		//If there are not at least 1500 phased polymorphic sites, return null
		int siteCount = 0;
		int nsites = myGeno.numberOfSites();
		for (int s = 0; s < nsites; s++) {
			if (phasedParent[0][s] != N && phasedParent[1][s] != N) siteCount++;
		}
		System.out.printf("There were %d polymorphic sites phased for %s\n", siteCount, parent);
		if (siteCount < 1500) return null;
		
		//add in the monomorphic sites
		for (int[] site : monomorphicSites) {
			phasedParent[0][site[0]] = (byte) site[1];
			phasedParent[1][site[0]] = (byte) site[1];
		}
		return phasedParent;
	}
	
	private byte maxAllele(int[] alleleCounts) {
		int max = 0;
		for (int i = 1; i < alleleCounts.length; i++) {
			if (alleleCounts[i] > alleleCounts[max]) max = i;
		}
		return (byte) max;
	}
	
	private int countSharedMembers(List<Integer> list1, List<Integer> list2) {
		List<Integer> refList = new ArrayList(list1);
		Collections.sort(refList);
		int count = 0;
		for (Integer I : list2) 
			if (Collections.binarySearch(refList, I) > -1) count++;
		return count;
	}
	
	private byte[] consensusHaplotype(List<byte[]> haplotypes) {
		int nsites = haplotypes.get(0).length;
		byte[] consensus = new byte[nsites];
		for (int i = 0; i < nsites; i++) {
			int[] alleleCounts = new int[6];
			for (byte[] hap:haplotypes) {
				if (hap[i] < 6) alleleCounts[hap[i]]++;
			}
			int ndx = 0;
			for (int j = 1; j < 6; j++) {
				if (alleleCounts[j] > alleleCounts[ndx]) ndx = j;
			}
			consensus[i] = (byte) ndx;
		}
		return consensus;
	}
	
	private String haplotypeAsString(byte[] hap) {
		StringBuilder sb = new StringBuilder();
		for (byte b:hap) sb.append(NucleotideAlignmentConstants.getHaplotypeNucleotide(b));
		return sb.toString();
	}
	
	private int[] countAllelesAtSite(List<byte[]> haps, int site) {
		int[] alleleCounts = new int[6];
		for (byte[] hap:haps) {
			if (hap[site] < 6)  alleleCounts[hap[site]]++;
		}
		return alleleCounts;
	}
	
	public byte[][] phaseParentUsingOneProgeny(String parent, String otherParent, String progeny, GenotypeTable gt) {
		//phase only at sites with sufficient depth for parents and progeny (all hets are okay, homozygotes with depth >= minDepth)
		int minDepth = 7;
		int nsites = gt.numberOfSites();
		byte[][] phasedGenotype = new byte[2][nsites];
		Arrays.fill(phasedGenotype[0], N);
		Arrays.fill(phasedGenotype[1], N);
		
		int parentIndex = gt.taxa().indexOf(parent);
		int otherParentIndex = gt.taxa().indexOf(otherParent);
		int progenyIndex = gt.taxa().indexOf(progeny);
		
		for (int s = 0; s < nsites; s++) {
			byte parentGenotype = gt.genotype(parentIndex, s);
			byte otherParentGenotype = gt.genotype(otherParentIndex, s);
			byte progenyGenotype = gt.genotype(progenyIndex, s);
			
			boolean parentHomozygous = !GenotypeTableUtils.isHeterozygous(parentGenotype) && gt.depth().depth(parentIndex, s) >= minDepth;
			boolean otherParentHomozygous = !GenotypeTableUtils.isHeterozygous(otherParentGenotype) && gt.depth().depth(otherParentIndex, s) >= minDepth;
			boolean progenyHomozygous = !GenotypeTableUtils.isHeterozygous(progenyGenotype) && gt.depth().depth(progenyIndex, s) >= minDepth;
			
			//1. if parent is homozygous (minDepth) then haplotype = parent allele
			//2. if progeny is homozygous (minDepth) then haplotype = progeny allele
			//3. if other parent is homozygous (minDepth) then
			//3a. if progeny is not homozygous other parent then parent is opposite allele of other parent
			//otherwise haplotype == N
			
			if (parentHomozygous) {
				phasedGenotype[0][s] = phasedGenotype[1][s] = GenotypeTableUtils.getDiploidValues(parentGenotype)[0];
			} else if (progenyHomozygous) {
				phasedGenotype[0][s] = GenotypeTableUtils.getDiploidValues(progenyGenotype)[0];
				if (parentGenotype != NN) {
					byte progenyAllele = GenotypeTableUtils.getDiploidValues(progenyGenotype)[0];
					byte[] parentAlleles = GenotypeTableUtils.getDiploidValues(parentGenotype);
					if (parentAlleles[0] != progenyAllele) phasedGenotype[1][s] = parentAlleles[0];
					else if (parentAlleles[1] != progenyAllele) phasedGenotype[1][s] = parentAlleles[1];
				}
			} else if (otherParentHomozygous) {
				byte[] progenyAlleles = GenotypeTableUtils.getDiploidValues(progenyGenotype);
				byte otherAllele = GenotypeTableUtils.getDiploidValues(otherParentGenotype)[0];
				if (progenyAlleles[0] != otherAllele || progenyAlleles[1] != otherAllele) {
					if (progenyAlleles[0] != otherAllele) phasedGenotype[0][s] = progenyAlleles[0];
					else phasedGenotype[0][s] = progenyAlleles[1];
					byte[] parentAlleles = GenotypeTableUtils.getDiploidValues(parentGenotype);
					if (parentAlleles[0] != phasedGenotype[0][s]) phasedGenotype[1][s] = parentAlleles[0];
					else if (parentAlleles[1] != phasedGenotype[0][s]) phasedGenotype[1][s] = parentAlleles[1];
				}
			} 
		}
		
		return phasedGenotype;
	}
	
	private double[][] mismatchDistanceMatrix(List<byte[]> segments) {
		int n = segments.size();
		double[][] dist = new double[n][n];
		for (int i = 0; i < n - 1; i++) {
			byte[] seg1 = segments.get(i);
			for (int j = i + 1; j < n; j++) {
				byte[] seg2 = segments.get(j);
				int notMissingCount = 0;
				int notMatchCount = 0;
				for (int k = 0; k < seg1.length; k++) {
					if (seg1[k] != N && seg2[k] != N) {
						notMissingCount++;
						if ((seg1[k] != seg2[k])) notMatchCount++;
					}
				}
				if (notMissingCount > 0) {
					dist[i][j] = dist[j][i] = notMatchCount / (double) notMissingCount;
				}
			}
		}
		return dist;
	}

	private double averageDistanceToCluster(double[][] distance, List<Integer> cluster, int ndx) {
		int n = cluster.size();
		double total = 0;
		for (int  i = 0; i < n; i++) {
			total += distance[cluster.get(i)][ndx];
		}
		return total / (double) n;
	}
	
	private void saveSeglist(List<byte[]> seglist, String filename) {
		try (BufferedWriter bw = Files.newBufferedWriter(Paths.get(filename))) {
			for (byte[] seg : seglist) {
				for (byte b : seg) bw.write(NucleotideAlignmentConstants.getHaplotypeNucleotide(b));
			}
			bw.write("\n");
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	private void saveDistanceMatrix(DistanceMatrix dm, String filename) {
		int n = dm.numberOfTaxa();
		try (BufferedWriter bw = Files.newBufferedWriter(Paths.get(filename))) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					bw.write(String.format("%1.5f ", dm.getDistance(i, j)));
				}
			}
		} catch(IOException e) {
			e.printStackTrace();
		}

	}
	
	public void setOuthapsAllProgeny(String filename) {
		outhapsAllProgeny = Paths.get(filename);
	}
	
	public void setParentage(String filename) {
		parentagePath = Paths.get(filename);
	}
	
	public void setGenotypeTable(GenotypeTable genotype) {
		myGenotypeTable = genotype;
	}
}
