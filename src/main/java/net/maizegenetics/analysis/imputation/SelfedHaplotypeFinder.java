package net.maizegenetics.analysis.imputation;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.analysis.clustering.Haplotype;
import net.maizegenetics.analysis.clustering.HaplotypeCluster;
import net.maizegenetics.analysis.clustering.HaplotypeClusterer;
import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;

public class SelfedHaplotypeFinder {
	private static final Logger myLogger = LogManager.getLogger(SelfedHaplotypeFinder.class);
	private static byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
	private static byte N = GenotypeTable.UNKNOWN_ALLELE;
	private static byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
	private static byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
	private static byte GG = NucleotideAlignmentConstants.getNucleotideDiploidByte("GG");
	private static byte TT = NucleotideAlignmentConstants.getNucleotideDiploidByte("TT");

	private GenotypeTable myGenotype;
	private int window;
	private double minNotMissingProportion;
	
	public SelfedHaplotypeFinder(int window, double minProportionNotMissing) {
		this.window = window;
		minNotMissingProportion = minProportionNotMissing;
	}
	
	public static void main(String[] args) {
		SelfedHaplotypeFinder shf = new SelfedHaplotypeFinder(1,1);
		GenotypeTable genoTable = (GenotypeTable) new FileLoadPlugin(null, false).runPlugin("/Users/pbradbury/Documents/projects/landraces/data/Combined_LR13_LR14_parents2_AGPv4.h5");
		shf.setGenotype(genoTable);
		shf.correctSelfPhaseUsingCross("/Users/pbradbury/temp/test_phasedSelfParents_lr_agpv4_dec12.bin", 
				"/Users/pbradbury/Documents/projects/landraces/output/dec12/phaseHighCoverParents_lr_agpv4_dec12.bin", 
				"165_7_Tuxpeno_El_Aguacatito_Mex_:250358062", 9, "/Users/pbradbury/temp/test_adjusted_phasedSelfParents_lr_agpv4_dec12.bin");
	}
	
    public Map<String, byte[][]> phaseSelfedParents(Path parentpath, Path savepath) {
    	boolean singleParentChrom = false;
    	int minFamilySize = 10;
    	int nsites = myGenotype.numberOfSites();
    	long startTime = System.currentTimeMillis();
    	List<String[]> plotList = new ArrayList<>();
    	try (BufferedReader br = Files.newBufferedReader(parentpath)) {
    		br.readLine();
    		String input;
    		while ((input = br.readLine()) != null) {
    			String[] plot = input.split("\t");
    			if (plot[3].equals("self")) plotList.add(plot);
    		}
    	} catch (IOException e) {
    		throw new RuntimeException("Failed to read parentage.", e);
    	}
    	
		TreeSet<String> parentSet = plotList.stream().map(p -> p[1]).collect(Collectors.toCollection(TreeSet::new));
		
		if (singleParentChrom) {
			parentSet.clear();
			parentSet.add("165_7_Tuxpeno_El_Aguacatito_Mex_:250358062");
		}
		
		//process selfs
		Map<String, byte[][]> phasedHaplotypes = new HashMap<>();
		for (String parent : parentSet) {
			
			System.out.printf("Starting %s at %d ms.\n", parent, System.currentTimeMillis() - startTime);
			TaxaListBuilder taxaBuilder = new TaxaListBuilder();
			for (String[] plot : plotList) {
				if (plot[1].equals(parent)) taxaBuilder.add(new Taxon(plot[0])); 
			}
			TaxaList familyTaxa = taxaBuilder.build();
			
			//do not phase if the family size < minFamilySize
			if (familyTaxa.numberOfTaxa() < minFamilySize) {
				myLogger.info(String.format("%s not phase because family size is %d.", parent, familyTaxa.numberOfTaxa()));
				continue;
			}
			
			myLogger.info(String.format("Phasing %s, %d self progeny.", parent, familyTaxa.numberOfTaxa()));
			GenotypeTable familyGeno = FilterGenotypeTable.getInstance(myGenotype, familyTaxa, false);
			myLogger.info(String.format("%s, number of sites = %d, number of taxa = %d", parent, familyGeno.numberOfSites(), familyGeno.numberOfTaxa()));
			
			//phase polymorphic sites
			//filter out monomorhpic sites and low coverage sites
			int polyCount = 0;
			double minMaf = 0.1;
			int minCount = (int) (0.5 * familyGeno.numberOfTaxa());

			//create the phasedHaplotype array
			byte[][] myPhasedHap = new byte[2][nsites];
			Arrays.fill(myPhasedHap[0], N);
			Arrays.fill(myPhasedHap[1], N);
			
			//phase by chromosome
			int ntaxa = familyGeno.numberOfTaxa();
			Chromosome[] myChromosomes = myGenotype.chromosomes();
			int startIndex = 0;
			int limitIndex = myChromosomes.length;
			if (singleParentChrom) {
				startIndex = 8;
				limitIndex = 9;
			}
			
			for (int c = startIndex;c < limitIndex; c++) {
				Chromosome myChr = myChromosomes[c];
				
				GenotypeTable chrGeno = FilterGenotypeTable.getInstance(familyGeno, myChr);
				
				myLogger.debug(String.format("%s, chr %s has %d sites before filtering", parent, myChr.getName(), chrGeno.numberOfSites()));
				int[] sitesToKeep = new int[chrGeno.numberOfSites()];
				int nSitesKept = 0;
				for (int s = 0; s < chrGeno.numberOfSites(); s++) {
					double maf = chrGeno.minorAlleleFrequency(s);
					if (maf >= minMaf &&  chrGeno.totalNonMissingForSite(s) >= minCount) {
						sitesToKeep[nSitesKept++] = s;
					}
				}
				sitesToKeep = Arrays.copyOf(sitesToKeep, nSitesKept);
				myLogger.debug(String.format("%s, chr %s has %d sites retained", parent, myChr.getName(), nSitesKept));
				chrGeno = FilterGenotypeTable.getInstance(chrGeno, sitesToKeep);
				
				int minClusterSize = Math.max(ntaxa/10, 3);
				int ibdClusterSize = ntaxa/2;
				int chrSites = chrGeno.numberOfSites();
				
				if (singleParentChrom) {
					//debug
					System.out.printf("Chr %s has %d sites\n", myChr.getName(), chrSites);
					System.out.printf("first position %d, last position %d\n", chrGeno.chromosomalPosition(0), chrGeno.chromosomalPosition(chrSites - 1));
					int[] chrStart = myGenotype.chromosomesOffsets();
					System.out.printf("myGenotype start pos = %d, end pos = %d\n", myGenotype.chromosomalPosition(chrStart[c]), myGenotype.chromosomalPosition(chrStart[c + 1] - 1) );
					//end debug
				}
				
				int maxhet = 10;
				HaplotypeCluster chr0clusters = null;
				HaplotypeCluster chr1clusters = null;
				HaplotypeCluster previousChr0clusters = null;
				HaplotypeCluster previousChr1clusters = null;
				
				int startIncr = window;
				int diff = 6;
				
				for (int start = 0; start < chrSites; start += startIncr) {
					int windowSize = window;
					if (chrSites - start < 2 * window) { //at the end of the chromosome include final sites in the last full haplotype
						windowSize = chrSites - start;
						startIncr = windowSize;
						diff = diff * windowSize / window;
					}
					
					int minNotMissingAdjusted = (int) (windowSize * minNotMissingProportion);

					//find two distinct clusters that are different at every site (or missing at "bad" sites), if possible
					//if not find a single cluster that is homozygous or missing at every site

					HaplotypeClusterer myClusterMaker = clusterWindow(chrGeno, start, windowSize, diff, minNotMissingAdjusted); //1
					
					myClusterMaker.sortClusters(); //2
					myClusterMaker.moveAllHaplotypesToBiggestCluster(diff); //3
					
					myClusterMaker.removeHeterozygousClusters(maxhet);
					
					int nclus = 0, pos = 0;
					if (singleParentChrom) {
						//debug
						pos = chrGeno.chromosomalPosition(start);
						nclus = myClusterMaker.getNumberOfClusters();
						System.out.printf("at chr %s, %d, %d clusters after remove hets.\n", myChr.getName(), pos, nclus);
						if (pos > 146000000) {
							System.out.print("");
						}
						//end debug
					}

					if (myClusterMaker.getNumberOfClusters() == 0) continue; //too many bad sites in window
					
					//pick out good haplotypes
					List<HaplotypeCluster> clusterList = myClusterMaker.getClusterList();
					Iterator<HaplotypeCluster> hapit = clusterList.iterator();
					HaplotypeCluster hchead = hapit.next();
					byte[] head = hchead.getHaplotype();
					while (hapit.hasNext()) {
						HaplotypeCluster hcComp = hapit.next();
						if (Haplotype.getDistance(head, hcComp.getHaplotype()) < window) hapit.remove();
					}

					//if more than two haplotypes remain, use the first two (the two largest because the cluster list has been sorted)
					if (clusterList.size() > 2) {
						for (int i = clusterList.size() - 1; i > 1; i--) {
							clusterList.remove(i);
						}
					}

					//if there is only one haplotype left, is it good?
					boolean twoGoodClusters = false;
					boolean oneGoodCluster = false;
					boolean isIbd = false;
					boolean inChrOrder = true;
					
					if (singleParentChrom) {
						//debug
						nclus = clusterList.size();
						System.out.printf("at chr %s, %d, %d clusters: ", myChr.getName(), pos, nclus);
						for (int i = 0; i < nclus; i++) {
							System.out.printf(" %d", clusterList.get(i).getSize());
							System.out.print(" " + clusterTaxa(clusterList.get(i)));
						}
						System.out.println();
						//end debug
					}
					
					//test for total clusters
					int[] clusterSizes = new int[2];
					int totalSize = 0;
					int limit = Math.min(2, clusterList.size());
					for (int i = 0; i < limit; i++) {
						if (clusterList.get(i) != null) {
							clusterSizes[i] = clusterList.get(i).getSize();
							totalSize += clusterSizes[i];
						}
					}
					
					if (totalSize >= minClusterSize) {
						if (clusterList.size() == 2 && clusterSizes[0] > 1 && clusterSizes[1] > 1) {
							twoGoodClusters = true;
							if (previousChr0clusters == null || previousChr1clusters == null) {
								previousChr0clusters = chr0clusters = clusterList.get(0);
								previousChr1clusters = chr1clusters = clusterList.get(1);
							} else {
								int order = sameOrder(previousChr0clusters, previousChr1clusters, clusterList);
								if (order == 1 || order == 0) {
									previousChr0clusters = chr0clusters = clusterList.get(0);
									previousChr1clusters = chr1clusters = clusterList.get(1);
								} else {
									previousChr0clusters = chr0clusters = clusterList.get(1);
									previousChr1clusters = chr1clusters = clusterList.get(0);
									inChrOrder = false;
								}
							}
						} else if (clusterList.size() == 1 && clusterSizes[0] > ibdClusterSize) { //the section is ibd, do not update previous clusters
							oneGoodCluster = true;
							isIbd = true;
							chr0clusters = clusterList.get(0);
							chr1clusters = null;
						} else if ((clusterList.size() == 1 && clusterSizes[0] > minClusterSize) || (clusterList.size() == 2 && clusterSizes[1] == 1 && clusterList.get(0).getSize() > minClusterSize)) {
							oneGoodCluster = true;
							int order = sameOrder(previousChr0clusters, previousChr1clusters, clusterList);
							if (order > -1) {
								chr0clusters = clusterList.get(0);
								previousChr0clusters = chr0clusters;
								chr1clusters = null;
							} else {
								chr0clusters = null;
								chr1clusters = clusterList.get(0);
								previousChr1clusters = chr1clusters;
								inChrOrder = false;
							}
						}
					}
					
					if (singleParentChrom) {
						System.out.printf("twoGoodClusters, oneGoodCluster: %b, %b\n", twoGoodClusters, oneGoodCluster);
					}
					
					if (twoGoodClusters) {
						//for each site
						//if both haplotypes are in (A,C,G,T) then assign alleles to myPhasedHap
						//map chrGeno position to original geno site number
						byte[] hap0 = chr0clusters.getCensoredMajorityHaplotype(0.05, 1);
						byte[] hap1 = chr1clusters.getCensoredMajorityHaplotype(0.05, 1);
						
						for (int w = 0; w < windowSize; w++) {
							byte allele0 = homozygousDiploidToAllele(hap0[w]);
							byte allele1 = homozygousDiploidToAllele(hap1[w]);
							if (allele0 > -1 && allele1 > -1 && allele0 != allele1) {
								//get the site number from myGenotype
								int filterSite = start + w;
								int mygenoSite = myGenotype.siteOfPhysicalPosition(chrGeno.chromosomalPosition(filterSite), myChr);
								myPhasedHap[0][mygenoSite] = allele0;
								myPhasedHap[1][mygenoSite] = allele1;

								polyCount++;
							}
						}
					} else if (oneGoodCluster && !isIbd) {
						int ndx;
						byte[] hap;
						if (chr0clusters != null) {
							ndx = 0;
							hap = chr0clusters.getCensoredMajorityHaplotype(0.05, 1);
						}
						else {
							ndx = 1;
							hap = chr1clusters.getCensoredMajorityHaplotype(0.05, 1);
						}
						
						for (int w = 0; w < windowSize; w++) {
							byte allele = homozygousDiploidToAllele(hap[w]);
							if (allele > -1) {
								//get the site number from myGenotype
								int filterSite = start + w;
								int mygenoSite = myGenotype.siteOfPhysicalPosition(chrGeno.chromosomalPosition(filterSite), myChr);
								myPhasedHap[ndx][mygenoSite] = allele;
								polyCount++;
							}
						}
					} else if (oneGoodCluster && !isIbd) {
						byte[] hap = chr0clusters.getCensoredMajorityHaplotype(0.05, 1);
						for (int w = 0; w < windowSize; w++) {
							byte allele = homozygousDiploidToAllele(hap[w]);
							if (allele > -1) {
								//get the site number from myGenotype
								int filterSite = start + w;
								int mygenoSite = myGenotype.siteOfPhysicalPosition(chrGeno.chromosomalPosition(filterSite), myChr);
								myPhasedHap[0][mygenoSite] = allele;
								myPhasedHap[1][mygenoSite] = allele;
								polyCount++;
							}
						}
						
					}

				}
			}
			
			//add back the monomorphic sites for this parent for the whole genome
			nsites = familyGeno.numberOfSites();
			ntaxa = familyGeno.numberOfTaxa();
			int minCoverage = 10;
			int monoCount = 0;
			for (int s = 0; s < nsites; s++) {
				if (familyGeno.totalNonMissingForSite(s) >= minCoverage && familyGeno.majorAlleleFrequency(s) > 0.9999) {
					for (int i = 0; i < 2; i++) myPhasedHap[i][s] = familyGeno.majorAllele(s);
					monoCount++;
				}
			}
			
			phasedHaplotypes.put(parent, myPhasedHap);
			myLogger.info(String.format("parent %s: %d polymorphic sites and %d monomorphic sites add to haplotypes.", parent, polyCount, monoCount));
		}
		
		//save the phased haplotypes
		if (savepath != null) serializePhasedHaplotypes(phasedHaplotypes, savepath);
		return phasedHaplotypes;
    }

    public HaplotypeClusterer clusterWindow(GenotypeTable a, int start, int length, int maxdif, int minNotMissingPerHaplotype) {
    	int ntaxa = a.numberOfTaxa();
    	int end = start + length;
    	ArrayList<Haplotype> haps = new ArrayList<Haplotype>();

    	for (int t = 0; t < ntaxa; t++) {
    		Haplotype hap = new Haplotype(a.genotypeRange(t, start, end), t);
    		if (hap.notMissingCount >= minNotMissingPerHaplotype) haps.add(hap);
    	}


    	HaplotypeClusterer hc = new HaplotypeClusterer(haps);
    	hc.makeClusters();
    	if (maxdif > 0) hc.mergeClusters(maxdif);
    	return hc;
    }

    private int sameOrder(HaplotypeCluster chr0, HaplotypeCluster chr1, List<HaplotypeCluster> next) {
    	//1 : same order, -1 : opposite order, 0 : diff=same, inconclusive
    	int nhaps = next.size();
    	if (nhaps == 1) {
    		int count0 = taxaInCommon(chr0, next.get(0));
    		int count1  = taxaInCommon(chr1, next.get(0));
    		if (count1 > count0) return -1;
    		if (count0 > count1) return 1;
    		return 0;
    	}
    	
    	//if chr0 = null or chr1 = null, same = 0, diff = 0 and method will return true
    	int same = taxaInCommon(chr0, next.get(0)) + taxaInCommon(chr1, next.get(1));
    	int diff = taxaInCommon(chr0, next.get(1)) + taxaInCommon(chr1, next.get(0));

    	if (same > diff) return 1;
    	if (diff > same) return -1;
    	return 0;
    }

    private String clusterTaxa(HaplotypeCluster hc) {
    	return hc.getHaplotypeList().stream().map(h -> Integer.toString(h.taxonIndex)).collect(Collectors.joining(",", "(", ")"));
    }
    
    private int taxaInCommon(HaplotypeCluster hc0, HaplotypeCluster hc1) {
    	if (hc0 == null || hc1 == null) return 0;
    	int commonCount = 0;
    	for (Haplotype hap0 : hc0.getHaplotypeList()) {
    		int ndx = hap0.taxonIndex;
    		for (Haplotype hap1 : hc1.getHaplotypeList()) {
    			if (ndx == hap1.taxonIndex) {
    				commonCount++;
    				break;
    			}
    		}
    	}
    	return commonCount;
    }

    private byte homozygousDiploidToAllele(byte value) {
    	if (value == AA) return NucleotideAlignmentConstants.A_ALLELE;
    	if (value == CC) return NucleotideAlignmentConstants.C_ALLELE;
    	if (value == GG) return NucleotideAlignmentConstants.G_ALLELE;
    	if (value == TT) return NucleotideAlignmentConstants.T_ALLELE;
    	return -1;
    }

    public static void serializePhasedHaplotypes(Map<String, byte[][]> phasedHaps, Path savePath) {
		try {
			FileOutputStream fos = new FileOutputStream(savePath.toFile());
			ObjectOutputStream oos = new ObjectOutputStream(fos);
			oos.writeObject(phasedHaps);
			oos.close();
		} catch (IOException e) {
			throw new RuntimeException("Unable to save phased haplotypes.", e);
		}
    }

    public static Map<String, byte[][]> restorePhasedHaplotypes(Path restorePath) {
    	try {
    		FileInputStream fis = new FileInputStream(restorePath.toFile());
            ObjectInputStream ois = new ObjectInputStream(fis);
            Map<String, byte[][]> phasedHaps = (Map<String, byte[][]>) ois.readObject();
            ois.close();
            return phasedHaps;
		} catch (IOException | ClassNotFoundException e) {
			throw new RuntimeException("Unable to restore phased haplotypes.", e);
		}
    }
    
    public void setGenotype(GenotypeTable genotype) {
    	myGenotype = genotype;
    }
    
    public void correctSelfPhaseUsingCross(String selfPhaseFilename, String crossPhaseFilename, String parent, int chr, String outfileName) {
    	Map<String, byte[][]> selfMap = restorePhasedHaplotypes(Paths.get(selfPhaseFilename));
    	Map<String, byte[][]> crossMap = restorePhasedHaplotypes(Paths.get(crossPhaseFilename));
    	Map<String, byte[][]> outMap = new HashMap<>();
    	
    	int[] chrOffsets = myGenotype.chromosomesOffsets();
    	int start = chrOffsets[chr - 1];
    	int end;
    	
    	if (chr != 10) {
    		end = chrOffsets[chr];
    	} else {
    		end = myGenotype.numberOfSites();
    	}
    	int nchrSites = end - start;
    	
    	for (String taxonName : selfMap.keySet()) {
    		if (taxonName.equals(parent)) {
    			byte[][] selfhaps = selfMap.get(taxonName);
    			byte[][] crosshaps = crossMap.get(taxonName);
    			
    			// 1 = same phase, -1 = different phase, 0 = cannot tell, missing data
    			int[] same = new int[nchrSites];
    			for (int s = start; s < end; s++) {
    				if (selfhaps[0][s] != N && selfhaps[1][s] != N && crosshaps[0][s] != N && crosshaps[1][s] != N) {
    					if (selfhaps[0][s] != selfhaps[1][s] && crosshaps[0][s] != crosshaps[1][s]) {
        					int ndx = s - start;
        					if (selfhaps[0][s] == crosshaps[0][s]) {
        						same[ndx] = 1;
        					}
        					else {
//        						if (selfhaps[0][s] == crosshaps[1][s] && selfhaps[1][s] == crosshaps[0][s]) same[ndx] = -1;
        						same[ndx] = -1;
        					}
    					}
    				}
    			}
    			
    			//Any sites flanked by 1's can be left as is
    			//Sites flanked by -1, should have haplotypes reversed
    			//Sites flanked by 1 and -1 should be set to missing
    			//Treat ends of chromosomes as flanked by first (or last) 1 or -1
    			
    			//report counts
    			int[] sameCount = new int[3];
    			for (int s = 1; s < nchrSites; s++) {
    				sameCount[same[s] + 1]++;
    			}
    			System.out.printf("For %s chr %d same counts -1: %d, 0: %d, 1: %d\n", parent, chr, sameCount[0], sameCount[1], sameCount[2]);
    			
    			//fill gaps in same array
    			int currentStart = 0;
    			for (int s = 1; s < nchrSites; s++) {
    				if (same[s] != 0) {
    					if (same[currentStart] == 0) {
    						for (int i = currentStart; i < s; i++) same[i] = same[s];
    					} else if (same[currentStart] == same[s]) {
    						for (int i = currentStart + 1; i < s; i++) same[i] = same[s];
    					}
    					currentStart = s;
    				}
    			}
    			
    			//finish filling end of chromosome
    			for (int i = currentStart + 1; i < nchrSites; i++) same[i] = same[currentStart];
    			
    			//report counts
    			sameCount = new int[3];
    			for (int s = 1; s < nchrSites; s++) {
    				sameCount[same[s] + 1]++;
    			}
    			System.out.printf("For %s chr %d same counts: -1 -> %d, 0 -> %d, 1 -> %d\n", parent, chr, sameCount[0], sameCount[1], sameCount[2]);
    			
    			//flip haplotypes as necessary
    			for (int s = start; s < end; s++) {
    				int ndx = s - start;
    				if (same[ndx] == -1) { //flip haplotypes
    					byte tmp = selfhaps[0][s];
    					selfhaps[0][s] = selfhaps[1][s];
    					selfhaps[1][s] = tmp;
    				} else if (same[ndx] == 0) { //set to missing
    					selfhaps[0][s] = selfhaps[1][s] = N;
    				}
    			}
    			outMap.put(taxonName, selfhaps);
    			
    			
    		} else {
    			outMap.put(taxonName, selfMap.get(taxonName));
    		}
    	}
    	
    	
    	//save outMap
    	serializePhasedHaplotypes(outMap, Paths.get(outfileName));
    }
    
    
}
