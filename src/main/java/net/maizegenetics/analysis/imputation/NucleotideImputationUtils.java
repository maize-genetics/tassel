package net.maizegenetics.analysis.imputation;

import net.maizegenetics.analysis.clustering.Haplotype;
import net.maizegenetics.analysis.clustering.HaplotypeCluster;
import net.maizegenetics.analysis.clustering.HaplotypeClusterer;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import com.google.common.collect.Multiset.Entry;

import net.maizegenetics.stats.math.GammaFunction;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.TaxaListUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;

import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.*;

public class NucleotideImputationUtils {
	private static final Logger myLogger = LogManager.getLogger(NucleotideImputationUtils.class);
	static final byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
	static final byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
	static final byte GG = NucleotideAlignmentConstants.getNucleotideDiploidByte("GG");
	static final byte TT = NucleotideAlignmentConstants.getNucleotideDiploidByte("TT");
	static final byte AC = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
	static final byte AG = NucleotideAlignmentConstants.getNucleotideDiploidByte("AG");
	static final byte AT = NucleotideAlignmentConstants.getNucleotideDiploidByte("AT");
	static final byte CG = NucleotideAlignmentConstants.getNucleotideDiploidByte("CG");
	static final byte CT = NucleotideAlignmentConstants.getNucleotideDiploidByte("CT");
	static final byte GT = NucleotideAlignmentConstants.getNucleotideDiploidByte("GT");
	static final byte NN = NucleotideAlignmentConstants.getNucleotideDiploidByte("NN");
	static final byte CA = NucleotideAlignmentConstants.getNucleotideDiploidByte("CA");

	static final byte[] byteval = new byte[] {AA,CC,GG,TT,AC};
	final static HashMap<Byte, Integer> genotypeMap = new HashMap<Byte, Integer>();

	static{
		genotypeMap.put(AA, 0);
		genotypeMap.put(CC, 1);
		genotypeMap.put(GG, 2);
		genotypeMap.put(TT, 3);
		genotypeMap.put(AC, 4);
		genotypeMap.put(CA, 4);
	}

	static final byte[][] genoval = new byte[][]{{AA,AC,AG,AT},{AC,CC,CG,CT},{AG,CG,GG,GT},{AT,CT,GT,TT}};

	//prevents instantiation of this class
	private NucleotideImputationUtils() {}

	public static void callParentAllelesByWindow(PopulationData popdata, double maxMissing, double minMaf, int windowSize, double minR) {
		byte N = GenotypeTable.UNKNOWN_ALLELE;
		BitSet polybits;
		double segratio = popdata.contribution1;
		if (minMaf < 0 && (segratio == 0.5 || segratio == 0.25 || segratio == 0.75)) {
			polybits = whichSitesSegregateCorrectly(popdata.original, maxMissing, segratio);
		} else {
			if (minMaf < 0) minMaf = 0;
			polybits = whichSitesArePolymorphic(popdata.original, maxMissing, minMaf);
		}

		myLogger.info("polybits cardinality = " + polybits.cardinality());

		OpenBitSet filteredBits = whichSnpsAreFromSameTag(popdata.original, polybits);
		myLogger.info("filteredBits.cardinality = " + filteredBits.cardinality());

		BitSet ldFilteredBits;
		if (minR > 0) {
			int halfWindow = windowSize / 2;
			ldFilteredBits = ldfilter(popdata.original, halfWindow, minR, filteredBits);
		} else {
			ldFilteredBits = filteredBits;
		}
		myLogger.info("ldFilteredBits.cardinality = " + ldFilteredBits.cardinality());

		int nsites = popdata.original.numberOfSites();
		int ntaxa = popdata.original.numberOfTaxa();
		popdata.alleleA = new byte[nsites];
		popdata.alleleC = new byte[nsites];
		popdata.snpIndex = new OpenBitSet(nsites);
		for (int s = 0; s < nsites; s++) {
			popdata.alleleA[s] = N;
			popdata.alleleC[s] = N;
		}

		int[] parentIndex = new int[2];
		parentIndex[0] = popdata.original.taxa().indexOf(popdata.parent1);
		parentIndex[1] = popdata.original.taxa().indexOf(popdata.parent2);

		//iterate through windows
		GenotypeTable[] prevAlignment = null;
		int[][] snpIndices = getWindows(ldFilteredBits, windowSize);

		boolean append = false;

		int nWindows = snpIndices.length;
		for (int w = 0; w < nWindows; w++) {
			int[] snpIndex;
			if (append) {
				int n1 = snpIndices[w-1].length;
				int n2 = snpIndices[w].length;
				snpIndex = new int[n1 + n2];
				System.arraycopy(snpIndices[w-1], 0, snpIndex, 0, n1);
				System.arraycopy(snpIndices[w], 0, snpIndex, n1, n2);
				append = false;
			} else {
				snpIndex = snpIndices[w];
			}

			//mask the hets
			//			MutableNucleotideAlignment mna = MutableNucleotideAlignment.getInstance(FilterAlignment.getInstance(popdata.original, snpIndex));
			//			int n = snpIndex.length;
			//			byte missing = NucleotideAlignmentConstants.getNucleotideDiploidByte('N');
			//			for (int i = 0; i < n; i++) {
			//				for (int t = 0; t < ntaxa; t++) {
			//					if (hetMask[t].fastGet(snpIndex[i])) mna.setBase(t, i, missing);
			//				}
			//			}
			//			mna.clean();

			GenotypeTable windowAlignment = FilterGenotypeTable.getInstance(popdata.original, snpIndex);
			LinkedList<Integer> snpList = new LinkedList<Integer>(); //snpList is a list of snps (indices) in this window
			for (int s:snpIndex) snpList.add(s);

			GenotypeTable[] taxaAlignments = getTaxaGroupAlignments(windowAlignment, parentIndex, snpList);

			if (taxaAlignments == null) {
				append = true;
			} else {
				//are groups in this alignment correlated with groups in the previous alignment
				double r = 0;
				if (prevAlignment != null) {
					r = getIdCorrelation(new TaxaList[][] {{prevAlignment[0].taxa(), prevAlignment[1].taxa()},{taxaAlignments[0].taxa(), taxaAlignments[1].taxa()}});
					myLogger.info("For " + popdata.name + " the window starting at " + popdata.original.siteName(snpIndex[0]) + ", r = " + r + " , # of snps in alignment = " + snpList.size());
				} else {
					myLogger.info("For " + popdata.name + " the window starting at " + popdata.original.siteName(snpIndex[0]) + ", # of snps in alignment = " + snpList.size());
				}

				checkAlignmentOrderIgnoringParents(taxaAlignments, popdata, r); //if r is negative switch alignment order (which will result in a positive r) 

				//debug -check upgma tree
				//				int[] selectSnps = new int[snpList.size()];
				//				int cnt = 0;
				//				for (Integer s : snpList) selectSnps[cnt++] = s;
				//				SBitAlignment sba = SBitAlignment.getInstance(FilterAlignment.getInstance(popdata.original, selectSnps));
				//				IBSDistanceMatrix dm = new IBSDistanceMatrix(sba);
				//				estimateMissingDistances(dm);
				//				Tree myTree = new UPGMATree(dm);
				//				TreeDisplayPlugin tdp = new TreeDisplayPlugin(null, true);
				//				tdp.performFunction(new DataSet(new Datum("Snp Tree", myTree, "Snp Tree"), null));

				prevAlignment = taxaAlignments;
				callParentAllelesUsingTaxaGroups(popdata, taxaAlignments, snpList);
			}

		}

		myLogger.info("number of called snps = " + popdata.snpIndex.cardinality());

		//create the imputed array with A/C calls
		int nsnps = (int) popdata.snpIndex.cardinality();
		ntaxa = popdata.original.numberOfTaxa();
		nsites = popdata.original.numberOfSites();
		int[] snpIndex = new int[nsnps];
		int snpcount = 0;
		PositionListBuilder posBuilder = new PositionListBuilder();
		for (int s = 0; s < nsites; s++) {
			if (popdata.snpIndex.fastGet(s)) {
				snpIndex[snpcount++] = s;
				posBuilder.add(popdata.original.positions().get(s));
			}
		}

		snpIndex = Arrays.copyOf(snpIndex, snpcount);
		GenotypeTableBuilder builder = GenotypeTableBuilder.getTaxaIncremental(posBuilder.build());

		nsnps = snpIndex.length;
		for (int t = 0; t < ntaxa; t++) {
			byte[] taxonGeno = new byte[nsnps];
			for (int s = 0; s < nsnps; s++) {
				byte Aallele = popdata.alleleA[snpIndex[s]];
				byte Callele = popdata.alleleC[snpIndex[s]];
				byte[] val = GenotypeTableUtils.getDiploidValues(popdata.original.genotype(t, snpIndex[s]));
				if (val[0] == Aallele && val[1] == Aallele) taxonGeno[s] = AA;
				else if (val[0] == Callele && val[1] == Callele) taxonGeno[s] = CC;
				else if ((val[0] == Aallele && val[1] == Callele) || (val[0] == Callele && val[1] == Aallele)) taxonGeno[s] = AC;
				else taxonGeno[s] = NN;
			}
			builder.addTaxon(popdata.original.taxa().get(t), taxonGeno);
		}

		popdata.imputed = builder.build(); 
	}

		public static void callParentAllelesUsingClusters(PopulationData popdata, double maxMissing, double minMaf, int windowSize, boolean checkSubpops) {
			HaplotypeCluster.ReturnHaplotype = HaplotypeCluster.TYPE.majority;
			int ntaxa = popdata.original.numberOfTaxa();
			int nOriginalSites = popdata.original.numberOfSites();
			popdata.alleleA = new byte[nOriginalSites];
			popdata.alleleC = new byte[nOriginalSites];
			popdata.snpIndex = new OpenBitSet(nOriginalSites);

			//		int minCount = (int) Math.ceil(1 - (maxMissing * ntaxa));
			//		Alignment a = AlignmentUtils.removeSitesBasedOnFreqIgnoreMissing(popdata.original, minMaf, 1.0, minCount);
			//		a = deleteSnpsFromSameTag(a);
			double maxHet = 0.06;
			if (minMaf < 0) minMaf = 0.05;
			GenotypeTable a = filterSnpsByTag(popdata.original, minMaf, maxMissing, maxHet);
			int nFilteredSites = a.numberOfSites();
			int[] siteTranslationArray = new int[nFilteredSites];
			//FilterGenotypeTable fa = (FilterGenotypeTable) a;
			//for (int s = 0; s < nFilteredSites; s++) siteTranslationArray[s] = fa.translateSite(s);
                        int[] translateSites = a.siteTranslations();
                        for (int s = 0; s < nFilteredSites; s++) siteTranslationArray[s] = translateSites[s];

			//find a window with only two distinct haplotype clusters
			int maxdist = 100;
			ArrayList<HaplotypeCluster> clusters = null;
			int start = 0;
			while (maxdist > 10) {
				clusters = clusterWindow(a, start, windowSize, 4);
				int mindist2 = 0, mindist3 = 0;
				if (clusters.size() > 2) mindist2 = Math.min(HaplotypeClusterer.clusterDistanceMaxPairDiff(clusters.get(0), clusters.get(2)), 
						HaplotypeClusterer.clusterDistanceMaxPairDiff(clusters.get(1), clusters.get(2)));
				if (clusters.size() > 3) mindist3 = Math.min(HaplotypeClusterer.clusterDistanceMaxPairDiff(clusters.get(0), clusters.get(3)),
						HaplotypeClusterer.clusterDistanceMaxPairDiff(clusters.get(1), clusters.get(3)));
				maxdist = Math.max(mindist2, mindist3);

				start += windowSize; 
			}
			HaplotypeCluster.ReturnHaplotype = HaplotypeCluster.TYPE.majority;
			for (int i = 0; i < 4; i++) if (clusters.size() > i) System.out.println(clusters.get(i));

			int[] sites = new int[windowSize];
			int currentSite = start - windowSize;
			for (int i = 0; i < windowSize; i++) sites[i] = currentSite++;

			//extend to the right
			ArrayList<HaplotypeCluster> seedClusters = new ArrayList<HaplotypeCluster>(clusters);
			int[] siteNumber = new int[]{windowSize};

			int[] firstWindowSites = null;
			while (siteNumber[0] == windowSize && sites[windowSize - 1] < nFilteredSites - 1) {
				seedClusters = extendClusters(a, seedClusters, sites, siteNumber, true);
				if (firstWindowSites == null) firstWindowSites = Arrays.copyOf(sites, sites.length);  //start here to extend to the left
				if (siteNumber[0] > 0) {
					if (siteNumber[0] < windowSize) sites = Arrays.copyOf(sites, siteNumber[0]);
					int[] origSites = translateSitesBackToOriginal(sites, siteTranslationArray);
					updateAlleleCalls(seedClusters.get(0), origSites, popdata.alleleA);
					updateAlleleCalls(seedClusters.get(1), origSites, popdata.alleleC);
					for (int site : origSites) popdata.snpIndex.fastSet(site);
				}
			}

			//extend to the left
			seedClusters = new ArrayList<HaplotypeCluster>(clusters);
			siteNumber = new int[]{windowSize};
			sites = firstWindowSites;

			while (siteNumber[0] == windowSize && sites[0] > 0) {
				seedClusters = extendClusters(a, seedClusters, sites, siteNumber, false);
				if (siteNumber[0] > 0) {
					if (siteNumber[0] < windowSize) sites = Arrays.copyOfRange(sites, windowSize - siteNumber[0], windowSize);
					int[] origSites = translateSitesBackToOriginal(sites, siteTranslationArray);
					updateAlleleCalls(seedClusters.get(0), origSites, popdata.alleleA);
					updateAlleleCalls(seedClusters.get(1), origSites, popdata.alleleC);
					for (int site : origSites) popdata.snpIndex.fastSet(site);
				}
			}

			//check parents to see if alleles A and C should be exchanged
			checkParentage(popdata);

			//create the imputed array with A/C calls
			int nsnps = (int) popdata.snpIndex.cardinality();
			ntaxa = popdata.original.numberOfTaxa();
			int nsites = popdata.original.numberOfSites();
			int[] snpIndex = new int[nsnps];
			int snpcount = 0;
			for (int s = 0; s < nsites; s++) {
				if (popdata.snpIndex.fastGet(s)) snpIndex[snpcount++] = s;
			}

			snpIndex = Arrays.copyOf(snpIndex, snpcount);
			nsnps = snpIndex.length;
			GenotypeTable target = FilterGenotypeTable.getInstance(popdata.original, snpIndex);

			createSubPops(popdata, checkSubpops);
			checksubpops(popdata, 20);

			//filter on taxa ids
			int npops = popdata.subpopulationGroups.size();
			TaxaList[] taxaLists = new TaxaList[popdata.subpopulationGroups.size()];
			popdata.subpopulationGroups.toArray(taxaLists);
			TaxaList allTaxa = TaxaListUtils.getAllTaxa(taxaLists);
			target = FilterGenotypeTable.getInstance(target, allTaxa);
			GenotypeTableBuilder builder = GenotypeTableBuilder.getTaxaIncremental(target.positions());

			for (int p = 0; p < npops; p++) {
				//create an index back to target for the taxa
				TaxaList popTaxaList = popdata.subpopulationGroups.get(p);
				ntaxa = popTaxaList.numberOfTaxa();
				int[] taxaIndex = new int[ntaxa];
				for (int t = 0; t < ntaxa; t++) taxaIndex[t] = target.taxa().indexOf(popTaxaList.get(t));

				BitSet subpopSnpUse = popdata.subpoplationSiteIndex.get(p);
				nsnps = target.numberOfSites();
				for (int t = 0; t < ntaxa; t++) {
					byte[] taxonGeno = new byte[nsnps];
					String gstr = target.genotypeAsStringRow(t);
					for (int s = 0; s < nsnps; s++) {
						byte Aallele = popdata.alleleA[snpIndex[s]];
						byte Callele = popdata.alleleC[snpIndex[s]];
						byte val = target.genotype(taxaIndex[t], s);
						if (val == Aallele) taxonGeno[s] = AA;
						else if (val == Callele) taxonGeno[s] = CC;
						else if (GenotypeTableUtils.isHeterozygous(val)) taxonGeno[s] = AC;
						else taxonGeno[s] = NN;
					}
					builder.addTaxon(popTaxaList.get(t), taxonGeno);
				}
			}

			popdata.imputed = builder.build(); 
		}

		public static void callParentAllelesUsingClustersOnly(PopulationData popdata, double maxMissing, double minMaf, int windowSize, boolean checkSubpops ) {
			double mindiff = 0.7; //the minimum distance for which two clusters will be considered to be from different haplotypes
			int overlap = 10; //overlap between windows used to link adjacent windows
			HaplotypeCluster.ReturnHaplotype = HaplotypeCluster.TYPE.majority;
			int ntaxa = popdata.original.numberOfTaxa();
			int nOriginalSites = popdata.original.numberOfSites();
			popdata.alleleA = new byte[nOriginalSites];
			popdata.alleleC = new byte[nOriginalSites];
			popdata.snpIndex = new OpenBitSet(nOriginalSites);

			double maxHet = 0.06;
			GenotypeTable a = filterSnpsByTag(popdata.original, minMaf, maxMissing, maxHet);
			int nFilteredSites = a.numberOfSites();
			int[] siteTranslationArray = new int[nFilteredSites];
			//FilterGenotypeTable fa = (FilterGenotypeTable) a;
			//for (int s = 0; s < nFilteredSites; s++) siteTranslationArray[s] = fa.translateSite(s);
                        int[] translateSites = a.siteTranslations();
                        for (int s = 0; s < nFilteredSites; s++) siteTranslationArray[s] = translateSites[s];

			//iterate through windows
			int maxsite = nFilteredSites - windowSize - overlap;
			ArrayList<HaplotypeCluster> clusters = null;
			HaplotypeCluster.ReturnHaplotype = HaplotypeCluster.TYPE.unanimous;
			byte[] prevA = null;
			byte[] prevC = null;
			for (int start = 0; start <= maxsite; start += windowSize - overlap) {
				if (start + windowSize > maxsite) clusters  = clusterWindow(a, start, nFilteredSites - start, 2);
				else clusters = clusterWindow(a, start, windowSize, 2);

				/*
				 *  Strategy for separating clusters, populations with homozygous parents
				 *  Cluster taxa in overlapping windows (3 sites is enough), maxdiff=4
				 *  If cluster diff proportion between biggest cluster and next biggest is < 0.8 (mindiff), try third biggest. If diff < 0.8 raise an error and quit.
				 *  Get haplotype from biggest cluster and the next biggest differnt cluster. Match the overlap to tell which haplotype it goes with.
				 *  If the overlap does not match perfectly raise an error.
				 * 
				 * */			
				double diff1 = HaplotypeClusterer.clusterDistanceClusterDiffProportion(clusters.get(0), clusters.get(1));
				double diff2 = HaplotypeClusterer.clusterDistanceClusterDiffProportion(clusters.get(0), clusters.get(2));
				double diff3 = HaplotypeClusterer.clusterDistanceClusterDiffProportion(clusters.get(0), clusters.get(3));
				double diff4 = HaplotypeClusterer.clusterDistanceClusterDiffProportion(clusters.get(0), clusters.get(4));

				byte[] hap0 = clusters.get(0).getHaplotype();
				byte[] hap1;
				if (diff1 >= mindiff) {
					hap1 = clusters.get(1).getHaplotype();
				} else if (diff2 >= mindiff){
					hap1 = clusters.get(2).getHaplotype();
				} else if (diff3 >= mindiff){
					hap1 = clusters.get(3).getHaplotype();
				} else if (diff4 >= mindiff){
					hap1 = clusters.get(4).getHaplotype();
				} else {
					throw new IllegalArgumentException(String.format("No second haplotype at site %d, position %d.", start, a.chromosomalPosition(start)));
				}

				if (start == 0) {
					//update sites that are different in hap0 and 1
					for (int s = 0; s < windowSize; s++) {
						if (hap0[s] != NN && hap1[s] != NN && hap0[s] != hap1[s]) {
							int origSite = siteTranslationArray[s];
							popdata.snpIndex.fastSet(origSite);
							popdata.alleleA[origSite] = hap0[s];
							popdata.alleleC[origSite] = hap1[s];
						}
					}
					prevA = hap0;
					prevC = hap1;
				} else {
					//use overlap to match this cluster to existing haplotypes
					//if there is a .8 or better match, invalidate any non-matching sites
					int seqlength = hap0.length;
					if (haplotypeOverlapSimilarity(prevA, hap0, overlap, windowSize) > 0.8 && haplotypeOverlapSimilarity(prevC, hap1, overlap, windowSize) > 0.8) {
						for (int s = 0; s < seqlength; s++) {
							boolean overlapMatchingSite = true;
							if (s < overlap) {
								int prevSite = s + windowSize - overlap;
								if (prevA[prevSite] != NN && hap0[s] != NN && (prevA[prevSite] != hap0[s] || prevC[prevSite] != hap1[s])) overlapMatchingSite = false;
							}
							if (hap0[s] != NN && hap1[s] != NN && hap0[s] != hap1[s] && overlapMatchingSite) {
								int origSite = siteTranslationArray[start + s];
								popdata.snpIndex.fastSet(origSite);
								popdata.alleleA[origSite] = hap0[s];
								popdata.alleleC[origSite] = hap1[s];
							}
						}
						prevA = hap0;
						prevC = hap1;	
					} else if (haplotypeOverlapSimilarity(prevA, hap1, overlap, windowSize) > 0.8 && haplotypeOverlapSimilarity(prevC, hap0, overlap, windowSize) > 0.8) {
						for (int s = 0; s < seqlength; s++) {
							boolean overlapMatchingSite = true;
							if (s < overlap) {
								int prevSite = s + windowSize - overlap;
								if (prevA[prevSite] != NN && hap1[s] != NN && (prevA[prevSite] != hap1[s] || prevC[prevSite] != hap0[s])) overlapMatchingSite = false;
							}
							if (hap0[s] != NN && hap1[s] != NN && hap0[s] != hap1[s] && overlapMatchingSite) {
								int origSite = siteTranslationArray[start + s];
								popdata.snpIndex.fastSet(origSite);
								popdata.alleleA[origSite] = hap1[s];
								popdata.alleleC[origSite] = hap0[s];
							}
						}
						prevA = hap1;
						prevC = hap0;	
					} else {
						throw new IllegalArgumentException(String.format("Haplotypes fail to approximately match previous haplotypes at site %d, position %d.", start, a.chromosomalPosition(start)));
					}
				}
			}

			//test each called site for LD with neighbors
			createSubPops(popdata, checkSubpops);
			checksubpops(popdata, 20);

			int nsnps = (int) popdata.snpIndex.cardinality();
			//		ntaxa = popdata.original.getSequenceCount();
			int nsites = popdata.original.numberOfSites();
			int[] snpIndex = new int[nsnps];
			int snpcount = 0;
			for (int s = 0; s < nsites; s++) {
				if (popdata.snpIndex.fastGet(s)) snpIndex[snpcount++] = s;
			}

			snpIndex = Arrays.copyOf(snpIndex, snpcount);
			nsnps = snpIndex.length;
			GenotypeTable target = FilterGenotypeTable.getInstance(popdata.original, snpIndex);
			GenotypeTableBuilder builder = GenotypeTableBuilder.getTaxaIncremental(target.positions());

			//assign parent calls to imputed alignment
			int nsubpops = popdata.subpopulationGroups.size();
			for (int sp = 0; sp < nsubpops; sp++) {
				BitSet SI = popdata.subpoplationSiteIndex.get(sp);
				TaxaList subpopTaxa = popdata.subpopulationGroups.get(sp);
				ntaxa = subpopTaxa.numberOfTaxa();
				for (int t = 0; t < ntaxa; t++) {
					int taxaTargetIndex = target.taxa().indexOf(subpopTaxa.get(t));
					byte[] taxonGeno = new byte[nsnps];
					for (int s = 0; s < nsnps; s++) {
						if (SI.fastGet(snpIndex[s])) {
							byte Aallele = popdata.alleleA[snpIndex[s]];
							byte Callele = popdata.alleleC[snpIndex[s]];
							byte[] val = GenotypeTableUtils.getDiploidValues(target.genotype(taxaTargetIndex, s));
							if (val[0] == Aallele && val[1] == Aallele) taxonGeno[s] = AA;
							else if (val[0] == Callele && val[1] == Callele) taxonGeno[s] = CC;
							else if ((val[0] == Aallele && val[1] == Callele) || (val[0] == Callele && val[1] == Aallele)) taxonGeno[s] = AC;
							else taxonGeno[s] = NN;
						}
					}
					builder.addTaxon(subpopTaxa.get(t), taxonGeno);
				}
			}


			popdata.imputed = builder.build(); 

		}

		public static boolean compareHaplotypes(byte[] h0, byte[] h1, int overlap, int haplength) {
			boolean same = true;
			int ptr = 0;
			int start = haplength - overlap;
			while (same && ptr < overlap) {
				if (h0[start + ptr] != h1[ptr] && h0[start + ptr] != NN && h1[ptr] != NN) same = false;
				ptr++;
			}
			return same;
		}

		public static double haplotypeOverlapSimilarity(byte[] h0, byte[] h1, int overlap, int haplength) {
			int ptr = 0;
			int start = haplength - overlap;
			int matchCount = 0;
			int nonmissingCount = 0;
			while (ptr < overlap) {
				if (h0[start + ptr] != NN && h1[ptr] != NN) {
					nonmissingCount++;
					if(h0[start + ptr] == h1[ptr]) matchCount++;
				}
				ptr++;
			}
			return ((double) matchCount) / ((double) nonmissingCount);
		}

		public static void createSubPops(PopulationData popdata, boolean isNAM) {
			popdata.subpopulationGroups = new ArrayList<TaxaList>();

			if (isNAM) {
				//create the population subgoups
				TaxaListBuilder[] sublists = new TaxaListBuilder[3];
				for (int i = 0; i < 3; i++) sublists[i] = new TaxaListBuilder();
				for (Taxon taxon:popdata.original.taxa()) {
					int pop = SubpopulationFinder.getNamSubPopulation(taxon);
					if (pop >= 0) sublists[pop].add(taxon);
				}

				//add the subgroups to popdata
				for (int i = 0; i < 3; i++) {
					TaxaList poplist = sublists[i].build();
					if (poplist.size() > 0) popdata.subpopulationGroups.add(poplist);
				}

			} else {
				popdata.subpopulationGroups.add(popdata.original.taxa());
			}
		}

		public static void checksubpops(PopulationData popdata, int halfWindowSize) {
			double minMatchScore = 0.9;

			int nsubpops = popdata.subpopulationGroups.size();
			int totalsites = popdata.original.numberOfSites();
			popdata.subpoplationSiteIndex = new ArrayList<BitSet>();
			for (int i = 0; i < nsubpops; i++) {
				GenotypeTable a = FilterGenotypeTable.getInstance(popdata.original, popdata.subpopulationGroups.get(i) );
				OpenBitSet snpUseIndex = new OpenBitSet(totalsites);
				popdata.subpoplationSiteIndex.add(snpUseIndex);

				//goodsites are the sites to be used for this population, that is the sites that have been assigned parent allele calls
				//remove sites with maf = 0 from goodsites
				int nsites = (int) popdata.snpIndex.cardinality();
				int[] goodsites = new int[nsites];
				int ptr = 0;
				for (int ts = 0; ts < totalsites; ts++) if (popdata.snpIndex.fastGet(ts) && a.minorAlleleFrequency(ts) > 0) goodsites[ptr++] = ts;
				nsites = ptr;
				goodsites = Arrays.copyOf(goodsites, ptr);

				//calculate site scores for the goodsites
				int ntaxa = a.numberOfTaxa();
				double[] score = new double[nsites];
				for (int s = 0; s < nsites; s++) {

					int firstSite = s - halfWindowSize;
					int lastSite = s + halfWindowSize;
					if (firstSite < 0) {
						lastSite -= firstSite;
						firstSite = 0;
					} else if (lastSite > nsites - 1) {
						firstSite -= lastSite - nsites + 1;
						lastSite = nsites - 1;
					}

					byte Aallele = popdata.alleleA[goodsites[s]];
					byte Callele = popdata.alleleC[goodsites[s]];
					int matchCount = 0;
					int totalCount = 0;
					for (int t = 0; t < ntaxa; t++) {
						byte tgeno = a.genotype(t, goodsites[s]);
						byte[] matchAllele;
						if (tgeno == Aallele || tgeno == Callele) {
							if (tgeno == Aallele) matchAllele = popdata.alleleA;
							else matchAllele = popdata.alleleC;
							for (int s2 = firstSite; s2 <= lastSite; s2++) if (s2 != s) {
								byte tgeno2 = a.genotype(t, goodsites[s2]);
								if (tgeno2 != NN) {
									totalCount += 1;
									if (tgeno2 == matchAllele[goodsites[s2]]) matchCount++;
								}
							}
						}
					}
					score[s] = ((double) matchCount) / ((double) totalCount);

				}

				//iterative improvement in this block (experimental)
				int iter = 1;
				double[] cuts = new double[]{0.7, 0.8, 0.9};
				for (double cutoff:cuts) {
					int[] bestsites = new int[nsites];
					ptr = 0;
					for (int s = 0; s < nsites; s++) {
						if (score[s] > cutoff) bestsites[ptr++] = goodsites[s];
					}
					bestsites = Arrays.copyOf(bestsites, ptr);
					int nbest = ptr;

					for (int s = 0; s < nsites; s++) {
						int closeSite = Arrays.binarySearch(bestsites, goodsites[s]);
						if (closeSite < 0) closeSite = -closeSite - 1;

						int firstSite = closeSite - halfWindowSize;
						int lastSite = closeSite + halfWindowSize;
						if (firstSite < 0) {
							lastSite -= firstSite;
							firstSite = 0;
						} else if (lastSite > nbest - 1) {
							firstSite -= lastSite - nbest + 1;
							lastSite = nbest - 1;
						}

						byte Aallele = popdata.alleleA[goodsites[s]];
						byte Callele = popdata.alleleC[goodsites[s]];
						int matchCount = 0;
						int totalCount = 0;
						for (int t = 0; t < ntaxa; t++) {
							byte tgeno = a.genotype(t, goodsites[s]);
							byte[] matchAllele;
							if (tgeno == Aallele || tgeno == Callele) {
								if (tgeno == Aallele) matchAllele = popdata.alleleA;
								else matchAllele = popdata.alleleC;
								for (int s2 = firstSite; s2 <= lastSite; s2++) if (bestsites[s2] != goodsites[s]) {
									byte tgeno2 = a.genotype(t, bestsites[s2]);
									if (tgeno2 != NN) {
										totalCount += 1;
										if (tgeno2 == matchAllele[bestsites[s2]]) matchCount++;
									}
								}
							}
						}

						double matchScore = ((double) matchCount) / ((double) totalCount);
						score[s] = matchScore;
					}
				}

				//use score to set use index
				//note that scores are the goodsite scores, so the snp index is contained in the goodsite index
				for (int s = 0; s < nsites; s++) {
					if (score[s] >= minMatchScore) {
						snpUseIndex.fastSet(goodsites[s]);
					}
				}
			}

		}

		public static int[] translateSitesBackToOriginal(int[] sites, int[] translation) {
			int n = sites.length;
			int[] orig = new int[n];
			for (int i = 0; i < n; i++) {
				orig[i] = translation[sites[i]];
			}
			return orig;
		}

		public static ArrayList<HaplotypeCluster> extendClusters(GenotypeTable a, ArrayList<HaplotypeCluster> clusters, int[] sites, int[] windowSize, boolean forward) {
			int maxDiff = 4;
			int totalSites = a.numberOfSites();
			int targetSize = windowSize[0];

			int direction;
			int nextSite;
			int nSelected = 0;
			if (forward) {
				direction = 1;
				nextSite = sites[sites.length - 1] + 1;
			} else {
				direction = -1;
				nextSite = sites[0] - 1;
			}

			String[] hap0 = new String[targetSize];
			String[] hap1 = new String[targetSize];
			while (nSelected < targetSize && nextSite < totalSites && nextSite > -1) {
				//snpset0 provides a count of alleles for cluster0, snpset1 for cluster1
				Multiset<String> snpset0 = HashMultiset.create();
				Multiset<String> snpset1 = HashMultiset.create();

				HaplotypeCluster hc = clusters.get(0);
				for (Iterator<Haplotype> iterator = hc.getIterator(); iterator.hasNext();) {
					Haplotype hap = iterator.next();
					int taxon = hap.taxonIndex;
					byte allele = a.genotype(taxon, nextSite);
					if (allele != Haplotype.N) {
						snpset0.add(NucleotideAlignmentConstants.getNucleotideIUPAC(allele));
					}
				}

				hc = clusters.get(1);
				for (Iterator<Haplotype> iterator = hc.getIterator(); iterator.hasNext();) {
					Haplotype hap = iterator.next();
					int taxon = hap.taxonIndex;
					byte allele = a.genotype(taxon, nextSite);
					if (allele != Haplotype.N) {
						snpset1.add(NucleotideAlignmentConstants.getNucleotideIUPAC(allele));
					}
				}

				//if either cluster has two alleles with roughly equal classes, then do not use
				String mj0 = getMajorAlleleFromSnpset(snpset0);
				if (mj0 != null) {  /*only one big class in snpset0*/
					String mj1 = getMajorAlleleFromSnpset(snpset1);
					if (mj1 != null && !mj0.equals(mj1)) {  /*only one big class in snpset1*/ /*the alleles are not equal*/
						if (forward) {
							hap0[nSelected] = mj0;
							hap1[nSelected] = mj1;
							sites[nSelected] = nextSite;
						} else {
							int ptr = targetSize - 1 - nSelected;
							hap0[ptr] = mj0;
							hap1[ptr] = mj1;
							sites[ptr] = nextSite;
						}
						nSelected ++;
					}
				}

				nextSite += direction;
			}

			windowSize[0] = nSelected;

			//create 2 new clusters
			//add taxa to each cluster that have diff <= 4 compared to either hap0 or hap1 but not both
			int[] selectedSites;
			if (forward) {
				selectedSites = Arrays.copyOf(sites, nSelected);
			} else {
				selectedSites = Arrays.copyOfRange(sites, targetSize - nSelected, targetSize);
			}

			int numberOfSites = selectedSites.length;
			byte[] seq = new byte[numberOfSites];
			if (forward) {
				for (int s = 0; s < numberOfSites; s++) seq[s] = NucleotideAlignmentConstants.getNucleotideDiploidByte(hap0[s]);
			} else {
				for (int s = 0; s < numberOfSites; s++) seq[s] = NucleotideAlignmentConstants.getNucleotideDiploidByte(hap0[s + targetSize - numberOfSites]);
			}
			Haplotype haplotype0 = new Haplotype(seq);

			seq = new byte[numberOfSites];
			if (forward || numberOfSites == targetSize) {
				for (int s = 0; s < numberOfSites; s++) seq[s] = NucleotideAlignmentConstants.getNucleotideDiploidByte(hap1[s]);
			} else {
				int dif = targetSize - numberOfSites;
				for (int s = 0; s < numberOfSites; s++) seq[s] = NucleotideAlignmentConstants.getNucleotideDiploidByte(hap1[s + dif]);
			}
			Haplotype haplotype1 = new Haplotype(seq);

			ArrayList<Haplotype> cluster0 = new ArrayList<>();
			ArrayList<Haplotype> cluster1 = new ArrayList<>();

			int ntaxa = a.numberOfTaxa();
			GenotypeTable b = FilterGenotypeTable.getInstance(a, selectedSites);
			for (int t = 0; t < ntaxa; t++) {
				seq = b.genotypeAllSites(t);
				Haplotype hap = new Haplotype(seq, t);
				int dist0 = hap.distanceFrom(haplotype0);
				int dist1 = hap.distanceFrom(haplotype1);
				if (dist0 <= maxDiff && dist1 > maxDiff) cluster0.add(hap);
				else if (dist0 > maxDiff && dist1 <= maxDiff) cluster1.add(hap);
				else if (numberOfSites < maxDiff) {
					if (dist0 == 0 && dist1 > 1) cluster0.add(hap);
					else if (dist0 > 1 && dist1 == 0) cluster1.add(hap);
				}
			}

			//return the clusters

			ArrayList<HaplotypeCluster> newClusters = new ArrayList<>();
			newClusters.add(new HaplotypeCluster(cluster0));
			newClusters.add(new HaplotypeCluster(cluster1));
			return newClusters;
		}

		public static String getMajorAlleleFromSnpset(Multiset<String> snpset, double alpha) {
			ArrayList<Multiset.Entry<String>> snplist = new ArrayList<>(snpset.entrySet());
			if (snplist.size() == 0) return null;
			if (snplist.size() == 1) return snplist.get(0).getElement();
			Collections.sort(snplist, new Comparator<Multiset.Entry<String>>(){

				@Override
				public int compare(Entry<String> entry1, Entry<String> entry2) {
					if (entry1.getCount() > entry2.getCount()) return -1;
					else if (entry1.getCount() < entry2.getCount()) return 1;
					return 0;
				}

			});

			long[] observed = new long[]{snplist.get(0).getCount(), snplist.get(1).getCount()};
			double total = observed[0] + observed[1];
			double[] expected = new double[]{ total/2.0, total/2.0 };

			boolean different = false;
			try {
				different =  TestUtils.chiSquareTest(expected, observed, alpha);
			} catch (Exception e) {
				System.out.println("Exception calculating chi-square in NucleotideImputationUtils.testClassSize(): ");
				e.printStackTrace();
			}

			if (different) {
				return snplist.get(0).getElement();
			}

			return null;
		}

		public static String getMajorAlleleFromSnpset(Multiset<String> snpset) {
			return getMajorAlleleFromSnpset(snpset, 0.05);
		}

		public static void updateAlleleCalls(HaplotypeCluster cluster, int[] sites, byte[] alleleCalls) {
			byte[] haplotype = cluster.getHaplotype();
			int n = haplotype.length;
			for (int i = 0; i < n; i++) alleleCalls[sites[i]] = haplotype[i];
		}

		public static void callParentAllelesByWindowForBackcrosses(PopulationData popdata, double maxMissing, double minMaf, int windowSize, double minR) {

			BitSet polybits = whichSitesSegregateCorrectly(popdata.original, maxMissing, .25);
			myLogger.info("polybits cardinality = " + polybits.cardinality());

			OpenBitSet filteredBits = whichSnpsAreFromSameTag(popdata.original, 0.8);
			filteredBits.and(polybits);
			System.out.println("filteredBits.cardinality = " + filteredBits.cardinality());

			BitSet ldFilteredBits;
			if (minR > 0) {
				int halfWindow = windowSize / 2;
				ldFilteredBits = ldfilter(popdata.original, halfWindow, minR, filteredBits);
			} else {
				ldFilteredBits = filteredBits;
			}
			myLogger.info("ldFilteredBits.cardinality = " + ldFilteredBits.cardinality());

			int nsites = popdata.original.numberOfSites();
			int ntaxa = popdata.original.numberOfTaxa();
			popdata.alleleA = new byte[nsites];
			popdata.alleleC = new byte[nsites];
			popdata.snpIndex = new OpenBitSet(nsites);
			for (int s = 0; s < nsites; s++) {
				if(ldFilteredBits.fastGet(s)) {
					popdata.alleleA[s] = popdata.original.majorAllele(s);
					popdata.alleleC[s] = popdata.original.minorAllele(s);
					popdata.snpIndex.fastSet(s);
				} 
			}


			myLogger.info("number of called snps = " + popdata.snpIndex.cardinality());

			//create the imputed array with A/C calls
			int nsnps = (int) popdata.snpIndex.cardinality();
			ntaxa = popdata.original.numberOfTaxa();
			nsites = popdata.original.numberOfSites();
			int[] snpIndex = new int[nsnps];
			int snpcount = 0;
			for (int s = 0; s < nsites; s++) {
				if (popdata.snpIndex.fastGet(s)) snpIndex[snpcount++] = s;
			}

			snpIndex = Arrays.copyOf(snpIndex, snpcount);
			GenotypeTable target = FilterGenotypeTable.getInstance(popdata.original, snpIndex);
			myLogger.info("filtered on snps");

			//remove taxa with low coverage
			TaxaListBuilder taxaBuilder = new TaxaListBuilder();
			int minGametes = 200;
			for (int t = 0; t < ntaxa; t++) {
				if (target.totalGametesNonMissingForTaxon(t) > minGametes) taxaBuilder.add(target.taxa().get(t));
			}

			GenotypeTable target2 = FilterGenotypeTable.getInstance(target, taxaBuilder.build());
			myLogger.info("identified low coverage taxa");

			nsnps = snpIndex.length;
			ntaxa = target2.numberOfTaxa();
			GenotypeTableBuilder builder = GenotypeTableBuilder.getSiteIncremental(target2.taxa());

			for (int s = 0; s < nsnps; s++) {
				byte Aallele = popdata.alleleA[snpIndex[s]];
				byte Callele = popdata.alleleC[snpIndex[s]];
				byte[] siteGeno = new byte[ntaxa];
				for (int t = 0; t < ntaxa; t++) {
					byte[] val = GenotypeTableUtils.getDiploidValues(target2.genotype(t, s));
					if (val[0] == Aallele && val[1] == Aallele) siteGeno[t] = AA;
					else if (val[0] == Callele && val[1] == Callele) siteGeno[t] = CC;
					else if ((val[0] == Aallele && val[1] == Callele)||(val[0] == Callele && val[1] == Aallele)) siteGeno[t] = AC;
					else siteGeno[t] = NN;
				}
				builder.addSite(target2.positions().get(s), siteGeno);
			}

			myLogger.info("called alleles");
			popdata.imputed = builder.build(); 
		}

		public static void callParentAllelesByWindowForMultipleBC(PopulationData popdata, double maxMissing, int  minMinorAlleleCount, int windowSize) {

			BitSet polybits = whichSitesArePolymorphic(popdata.original, maxMissing, minMinorAlleleCount);
			myLogger.info("polybits cardinality = " + polybits.cardinality());

			OpenBitSet filteredBits = whichSnpsAreFromSameTag(popdata.original, 0.8);
			filteredBits.and(polybits);
			System.out.println("filteredBits.cardinality = " + filteredBits.cardinality());

			int nsites = popdata.original.numberOfSites();
			int ntaxa = popdata.original.numberOfTaxa();
			popdata.alleleA = new byte[nsites];
			popdata.alleleC = new byte[nsites];
			popdata.snpIndex = new OpenBitSet(nsites);
			for (int s = 0; s < nsites; s++) {
				if(filteredBits.fastGet(s)) {
					popdata.alleleA[s] = popdata.original.majorAllele(s);
					popdata.alleleC[s] = popdata.original.minorAllele(s);
					popdata.snpIndex.fastSet(s);
				} 
			}

			myLogger.info("number of called snps = " + popdata.snpIndex.cardinality());

			//create the imputed array with A/C calls
			int nsnps = (int) popdata.snpIndex.cardinality();
			ntaxa = popdata.original.numberOfTaxa();
			nsites = popdata.original.numberOfSites();
			int[] snpIndex = new int[nsnps];
			int snpcount = 0;
			for (int s = 0; s < nsites; s++) {
				if (popdata.snpIndex.fastGet(s)) snpIndex[snpcount++] = s;
			}

			snpIndex = Arrays.copyOf(snpIndex, snpcount);
			GenotypeTable target = FilterGenotypeTable.getInstance(popdata.original, snpIndex);
			//		MutableNucleotideAlignment mna = MutableNucleotideAlignment.getInstance(target);

			nsnps = snpIndex.length;
			GenotypeTableBuilder builder = GenotypeTableBuilder.getSiteIncremental(target.taxa());

			for (int s = 0; s < nsnps; s++) {
				byte Aallele = popdata.alleleA[snpIndex[s]];
				byte Callele = popdata.alleleC[snpIndex[s]];
				byte[] siteGeno = new byte[ntaxa];
				for (int t = 0; t < ntaxa; t++) {
					byte[] val = GenotypeTableUtils.getDiploidValues(target.genotype(t, s));
					if (val[0] == Aallele && val[1] == Aallele) siteGeno[t] = AA;
					else if (val[0] == Callele && val[1] == Callele) siteGeno[t] = CC;
					else if ((val[0] == Aallele && val[1] == Callele)||(val[0] == Callele && val[1] == Aallele)) siteGeno[t] = AC;
					else siteGeno[t] = NN;
				}
				builder.addSite(target.positions().get(s), siteGeno);
			}

			myLogger.info("called alleles");
			popdata.imputed = builder.build(); 
		}

		public static void checkAlignmentOrderIgnoringParents(GenotypeTable[] alignments, PopulationData family, double r) {
			boolean swapAlignments = false;
			double minR = -0.05;

			int p1group, p2group;

			//r is the correlation between taxa assignments in this window and the previous one
			//if r is negative then parental assignments should be swapped, which will make it positive
			//this can happen when the parents are not in the data set and the assignment of parents to group is arbitrary
			if (r < minR) swapAlignments = true;
			if (swapAlignments) {
				GenotypeTable temp = alignments[0];
				alignments[0] = alignments[1];
				alignments[1] = temp;
			}

		}

		public static int[][] getWindows(BitSet ispoly, int windowSize) {
			int npoly = (int) ispoly.cardinality();
			int nsnps = (int) ispoly.size();
			int nwindows = npoly/windowSize;
			int remainder = npoly % windowSize;
			if (remainder > windowSize/2) nwindows++; //round up
			int[][] windows = new int[nwindows][];
			int setsize = npoly/nwindows;

			int windowCount = 0;
			int snpCount = 0;
			int polyCount = 0;
			while (snpCount < nsnps && windowCount < nwindows) {
				int numberLeft = npoly - polyCount;
				if (numberLeft < setsize * 2) setsize = numberLeft;
				int[] set = new int[setsize];
				int setcount = 0;
				while (setcount < setsize && snpCount < nsnps) {
					if (ispoly.fastGet(snpCount)) {
						set[setcount++] = snpCount;
						polyCount++;
					}
					snpCount++;
				}
				windows[windowCount++] = set;
			}

			return windows;
		}

		/**
		 * @param family	a PopulationData object containing information for this family
		 * @param taxaGroups	an array of two alignments corresponding to two clusters of taxa
		 * @param snpList	the list of snps to be called
		 */
		public static void callParentAllelesUsingTaxaGroups(PopulationData family, GenotypeTable[] taxaGroups, LinkedList<Integer> snpList) {
			int nsnps = taxaGroups[0].numberOfSites();
			Iterator<Integer> snpit = snpList.iterator();
			for ( int s = 0; s < nsnps; s++) {
				byte[] major = new byte[2];
				major[0] = taxaGroups[0].majorAllele(s);
				major[1] = taxaGroups[1].majorAllele(s);
				Integer snpIndex = snpit.next();
				if(major[0] != GenotypeTable.UNKNOWN_ALLELE && major[1] != GenotypeTable.UNKNOWN_ALLELE && major[0] != major[1]) {
					family.alleleA[snpIndex] = major[0];
					family.alleleC[snpIndex] = major[1];
					family.snpIndex.fastSet(snpIndex);
				}
			}
		}

		public static double getIdCorrelation(TaxaList[][] id) {
			double[][] counts = new double[2][2];
			counts[0][0] = TaxaListUtils.getCommonTaxa(id[0][0], id[1][0]).numberOfTaxa();
			counts[0][1] = TaxaListUtils.getCommonTaxa(id[0][0], id[1][1]).numberOfTaxa();
			counts[1][0] = TaxaListUtils.getCommonTaxa(id[0][1], id[1][0]).numberOfTaxa();
			counts[1][1] = TaxaListUtils.getCommonTaxa(id[0][1], id[1][1]).numberOfTaxa();
			double num = counts[0][0] * counts[1][1] - counts[0][1] * counts[1][0];
			double p1 = counts[0][0] + counts[0][1];
			double q1 = counts[1][0] + counts[1][1];
			double p2 =  counts[0][0] + counts[1][0];
			double q2 =  counts[0][1] + counts[1][1];
			return num / Math.sqrt(p1 * q1 * p2 * q2);
		}

		public static BitSet whichSitesArePolymorphic(GenotypeTable a, double maxMissing, int minMinorAlleleCount) {
			int ntaxa = a.numberOfTaxa();
			double minMaf = ((double) minMinorAlleleCount) / ((double) (2 * ntaxa));
			return whichSitesArePolymorphic(a, maxMissing, minMaf);
		}

		public static BitSet whichSitesArePolymorphic(GenotypeTable a, double maxMissing, double minMaf) {
			//which sites are polymorphic? minor allele count > 2 and exceed the minimum allele count
			int nsites = a.numberOfSites();
			int ntaxa = a.numberOfTaxa();
			OpenBitSet polybits = new OpenBitSet(nsites);
			for (int s = 0; s < nsites; s++) {
				int notMissingCount = a.totalNonMissingForSite(s);
				int minorAlleleCount = a.minorAlleleCount(s);
				double maf = a.minorAlleleFrequency(s);
				double pMissing = ((double) (ntaxa - notMissingCount)) / ((double) ntaxa);
				if (pMissing <= maxMissing && maf >= minMaf) polybits.fastSet(s);
			}
			return polybits;
		}

		public static BitSet whichSitesSegregateCorrectly(GenotypeTable a, double maxMissing, double ratio) {
			int nsites = a.numberOfSites();
			int ntaxa = a.numberOfTaxa();
			OpenBitSet polybits = new OpenBitSet(nsites);
			for (int s = 0; s < nsites; s++) {
				int[][] freq = a.allelesSortedByFrequency(s);
				int nMissing = ntaxa - a.totalNonMissingForSite(s);
				double pMissing = ((double) nMissing) / ((double) ntaxa);
				if (freq[1].length > 1 && pMissing <= maxMissing) {
					int Mj = freq[1][0];
					int Mn = freq[1][1];
					double pmono = binomialProbability(Mj + Mn, Mn, 0.002);
					double pquarter = binomialProbability(Mj + Mn, Mn, 0.25);
					double phalf = binomialProbability(Mj + Mn, Mn, 0.5);
					if (ratio == 0.25 || ratio == 0.75) {
						//					if (pquarter / (pmono + phalf) > 2) polybits.fastSet(s);
						//					double poneseven = binomialProbability(Mj + Mn, Mn, .125);
						//					double pthreefive = binomialProbability(Mj + Mn, Mn, .375);
						//					if (pquarter > poneseven && pquarter > pthreefive) polybits.fastSet(s);
						if (pquarter > phalf && pquarter > pmono) polybits.fastSet(s);

					} else {
						if (phalf / (pmono + pquarter) > 2) polybits.fastSet(s);
					}
				}
			}
			return polybits;

		}

		private static double binomialProbability(int trials, int successes, double pSuccess) {
			double n = trials;
			double k = successes;

			double logprob = GammaFunction.lnGamma(n + 1.0) -
					GammaFunction.lnGamma(k + 1.0) - GammaFunction.lnGamma(n - k + 1.0) + k * Math.log(pSuccess) + (n - k) * Math.log(1 - pSuccess);

			return Math.exp(logprob);
		}

		public static GenotypeTable[] getTaxaGroupAlignments(GenotypeTable a, int[] parentIndex, LinkedList<Integer> snpIndices) {

			//cluster taxa for these snps to find parental haplotypes (cluster on taxa)
			GenotypeTable[] taxaClusters = ImputationUtils.getTwoClusters(a, 20);
			LinkedList<Integer> originalList = new LinkedList<Integer>(snpIndices);
			int nsites = a.numberOfSites();
			boolean[] include = new boolean[nsites];
			int[] includedSnps = new int[nsites];
			int snpcount = 0;
			for (int s = 0; s < nsites; s++) {
				Integer snpIndex = snpIndices.remove();
				if ( taxaClusters[0].majorAllele(s) != taxaClusters[1].majorAllele(s) ) {
					include[s] = true;
					includedSnps[snpcount++] = s;
					snpIndices.add(snpIndex);
				} else {
					include[s] = false;
				}
			}

			if (snpcount > 5) {
				includedSnps = Arrays.copyOf(includedSnps, snpcount);
				if (snpcount == nsites) return taxaClusters;
				else return ImputationUtils.getTwoClusters(FilterGenotypeTable.getInstance(a, includedSnps), 20);
			} else {
				snpIndices.clear();
				snpIndices.addAll(originalList);
				return taxaClusters;
			}
		}

		public static double computeRForAlleles(int site1, int site2, GenotypeTable a) {
			int s1Count = 0;
			int s2Count = 0;
			int prodCount = 0;
			int totalCount = 0;

			long[] m11 = a.allelePresenceForAllTaxa(site1, WHICH_ALLELE.Major).getBits();
			long[] m12 = a.allelePresenceForAllTaxa(site1, WHICH_ALLELE.Minor).getBits();
			long[] m21 = a.allelePresenceForAllTaxa(site2, WHICH_ALLELE.Major).getBits();
			long[] m22 = a.allelePresenceForAllTaxa(site2, WHICH_ALLELE.Minor).getBits();
			int n = m11.length;
			for (int i = 0; i < n; i++) {

				long valid = (m11[i] ^ m12[i]) & (m21[i] ^ m22[i]); //only count non-het & non-missing
				long s1major = m11[i] & valid;
				long s2major = m21[i] & valid;
				long s1s2major = s1major & s2major & valid;
				s1Count += BitUtil.pop(s1major);
				s2Count += BitUtil.pop(s2major);
				prodCount += BitUtil.pop(s1s2major);
				totalCount += BitUtil.pop(valid);
			}

			if (totalCount < 2) return Double.NaN;

			//Explanation of method:
			//  if major site one is x=1, minor site one is x = 0 and for site 2 y = 1 or 0
			//  r = [sum(xy) - sum(x)sum(y)/N] / sqrt[(sum(x) - sum(x)*sum(x)/N) * ((sum(y) - sum(y)*sum(y)/N)]
			//  and sum(x) - sum(x)*sum(x)/N = sum(x)(N - sum(x))/N
			//  because sum(x^2) = sum(x)
			double num = ((double) prodCount - ((double) s1Count * s2Count) / ((double) totalCount));
			double denom = ((double) (s1Count * (totalCount - s1Count))) / ((double) totalCount);
			denom *= ((double) (s2Count * (totalCount - s2Count))) / ((double) totalCount);
			if (denom == 0) return Double.NaN;
			return  num / Math.sqrt(denom);
		}

		public static double computeRForMissingness(int site1, int site2, GenotypeTable a) {
			//		int s1Count = 0;
			//		int s2Count = 0;
			//		int prodCount = 0;

			int totalCount = a.numberOfTaxa();

			BitSet m11 = a.allelePresenceForAllTaxa(site1, WHICH_ALLELE.Major);
			BitSet m12 = a.allelePresenceForAllTaxa(site1, WHICH_ALLELE.Minor);
			BitSet m21 = a.allelePresenceForAllTaxa(site2, WHICH_ALLELE.Major);
			BitSet m22 = a.allelePresenceForAllTaxa(site2, WHICH_ALLELE.Minor);

			OpenBitSet s1present = new OpenBitSet(m11.getBits(), m11.getNumWords());
			s1present.union(m12);

			OpenBitSet s2present = new OpenBitSet(m21.getBits(), m21.getNumWords());
			s2present.union(m22);
			long s1Count = s1present.cardinality();
			long s2Count = s2present.cardinality();
			long prodCount = OpenBitSet.intersectionCount(s1present, s2present);

			//Explanation of method:
			//  if major site one is x=1, minor site one is x = 0 and for site 2 y = 1 or 0
			//  r = [sum(xy) - sum(x)sum(y)/N] / sqrt[(sum(x) - sum(x)*sum(x)/N) * ((sum(y) - sum(y)*sum(y)/N)]
			//  and sum(x) - sum(x)*sum(x)/N = sum(x)(N - sum(x))/N
			//  because sum(x^2) = sum(x)
			double num = ((double) prodCount - ((double) s1Count * s2Count) / ((double) totalCount));
			double denom = ((double) (s1Count * (totalCount - s1Count))) / ((double) totalCount);
			denom *= ((double) (s2Count * (totalCount - s2Count))) / ((double) totalCount);
			if (denom == 0) return Double.NaN;
			return  num / Math.sqrt(denom);
		}

		public static boolean areSameTag(int site1, int site2, GenotypeTable a) {
			int pos1 = a.chromosomalPosition(site1);
			int pos2 = a.chromosomalPosition(site2);
			if (Math.abs(pos2 - pos1) > 64) return false;
			long[] matchCount = matchCount(site1, site2, a);
			long numberNotMatching = matchCount[1] - matchCount[0];
			if (numberNotMatching > 5) return false;
			return true;
		}

		public static long[] matchCount(int site1, int site2, GenotypeTable a) {

			int totalCount = a.numberOfTaxa();

			BitSet m11 = a.allelePresenceForAllTaxa(site1, WHICH_ALLELE.Major);
			BitSet m12 = a.allelePresenceForAllTaxa(site1, WHICH_ALLELE.Minor);
			BitSet m21 = a.allelePresenceForAllTaxa(site2, WHICH_ALLELE.Major);
			BitSet m22 = a.allelePresenceForAllTaxa(site2, WHICH_ALLELE.Minor);

			OpenBitSet majorSame = new OpenBitSet(m11);
			majorSame.notXor(m21);
			OpenBitSet minorSame = new OpenBitSet(m12);
			minorSame.notXor(m22);
			OpenBitSet equal1 = new OpenBitSet(majorSame);
			equal1.and(minorSame);

			OpenBitSet majorMinorSame = new OpenBitSet(m11.getBits());
			majorMinorSame.notXor(m22);
			OpenBitSet minorMajorSame = new OpenBitSet(m12.getBits());
			minorMajorSame.notXor(m21);
			OpenBitSet equal2 = new OpenBitSet(majorMinorSame.getBits());
			equal2.and(minorMajorSame);

			//Sites for which both the minor and minor alleles have the same value
			long majorCardinality = equal1.cardinality(totalCount - 1);
			long majorMinorCardinality = equal2.cardinality(totalCount - 1);

			if (majorCardinality >= majorMinorCardinality) return new long[]{majorCardinality, totalCount};
			return new long[]{majorMinorCardinality, totalCount}; 
		}

		public static double computeGenotypeR(int site1, int site2, GenotypeTable a) throws IllegalStateException {

			BitSet s1mj = a.allelePresenceForAllTaxa(site1, WHICH_ALLELE.Major);
			BitSet s1mn = a.allelePresenceForAllTaxa(site1, WHICH_ALLELE.Minor);
			BitSet s2mj = a.allelePresenceForAllTaxa(site2, WHICH_ALLELE.Major);
			BitSet s2mn = a.allelePresenceForAllTaxa(site2, WHICH_ALLELE.Minor);

			OpenBitSet bothpresent = new OpenBitSet(s1mj);
			bothpresent.union(s1mn);
			OpenBitSet s2present = new OpenBitSet(s2mj);
			s2present.union(s2mn);
			bothpresent.intersect(s2present);

			long nsites = a.numberOfTaxa();
			double sum1 = 0;
			double sum2 = 0;
			double sumsq1 = 0;
			double sumsq2 = 0;
			double sumprod = 0;
			double count = 0;

			for (long i = 0; i < nsites; i++) {
				if (bothpresent.fastGet(i)) {
					double val1 = 0; 
					double val2 = 0;

					if (s1mj.fastGet(i)) {
						if (s1mn.fastGet(i)) {
							val1++;
						} else {
							val1 += 2;
						}
					}

					if (s2mj.fastGet(i)) {
						if (s2mn.fastGet(i)) {
							val2++;
						} else {
							val2 += 2;
						}
					}

					sum1 += val1;
					sumsq1 += val1*val1;
					sum2 += val2;
					sumsq2 += val2*val2;
					sumprod += val1*val2;
					count++;
				}
			}

			//r = sum(xy)/sqrt[sum(x2)*sum(y2)], where x = X - Xbar and y = Y - Ybar
			//num.r = sum(XY) - 1/n * sum(X)sum(y)
			//denom.r = sqrt[ (sum(XX) - 1/n*sum(X)sum(X)) * (sum(YY) - 1/n*sum(Y)sum(Y)) ]
			double num = sumprod - sum1 / count * sum2;
			double denom = Math.sqrt( (sumsq1 - sum1 / count * sum1) * (sumsq2 - sum2 / count * sum2) );
			if (denom == 0) return Double.NaN;
			return num / denom;
		}

		public static GenotypeTable imputeUsingViterbiFiveState(GenotypeTable a, double probHeterozygous, String familyName) {
			return imputeUsingViterbiFiveState(a, probHeterozygous, familyName, false);
		}

		public static GenotypeTable imputeUsingViterbiFiveState(GenotypeTable a, double probHeterozygous, String familyName, boolean useVariableRecombitionRates) {
			//states are in {all A; 3A:1C; 1A:1C, 1A:3C; all C}
			//obs are in {A, C, M}, where M is heterozygote A/C
			int maxIterations = 50;
			HashMap<Byte, Byte> obsMap = new HashMap<Byte, Byte>();
			obsMap.put(AA, (byte) 0);
			obsMap.put(AC, (byte) 1);
			obsMap.put(CA, (byte) 1);
			obsMap.put(CC, (byte) 2);

			int ntaxa = a.numberOfTaxa();
			int nsites = a.numberOfSites();

			//initialize the transition matrix
			double[][] transition = new double[][] {
					{.999,.0001,.0003,.0001,.0005},
					{.0002,.999,.00005,.00005,.0002},
					{.0002,.00005,.999,.00005,.0002},
					{.0002,.00005,.00005,.999,.0002},
					{.0005,.0001,.0003,.0001,.999}
			};

			TransitionProbability tp;
			if (useVariableRecombitionRates) {
				tp = new TransitionProbabilityWithVariableRecombination(a.chromosomeName(0));
			} else {
				tp = new TransitionProbability();
			}

			tp.setTransitionProbability(transition);
			int chrlength = a.chromosomalPosition(nsites - 1) - a.chromosomalPosition(0);
			tp.setAverageSegmentLength( chrlength / nsites );

			//initialize the emission matrix, states (5) in rows, observations (3) in columns
			double[][] emission = new double[][] {
					{.998,.001,.001},
					{.6,.2,.2},
					{.4,.2,.4},
					{.2,.2,.6},
					{.001,.001,.998}
			};

			EmissionProbability ep = new EmissionProbability();
			ep.setEmissionProbability(emission);

			//set up indices to non-missing data
			ArrayList<BitSet> notMissingIndex = new ArrayList<BitSet>();
			int[] notMissingCount = new int[ntaxa];
			ArrayList<byte[]> nonMissingObs = new ArrayList<byte[]>();
			ArrayList<int[]> snpPositions = new ArrayList<int[]>();

			for (int t = 0; t < ntaxa; t++) {
				long[] bits = a.allelePresenceForAllSites(t, WHICH_ALLELE.Major).getBits();
				BitSet notMiss = new OpenBitSet(bits, bits.length);
				notMiss.or(a.allelePresenceForAllSites(t, WHICH_ALLELE.Minor));
				notMissingIndex.add(notMiss);
				notMissingCount[t] = (int) notMiss.cardinality();
			}

			for (int t = 0; t < ntaxa; t++) {
				byte[] obs = new byte[notMissingCount[t]];
				int[] pos = new int[notMissingCount[t]];
				nonMissingObs.add(obs);
				snpPositions.add(pos);
				BitSet isNotMissing = notMissingIndex.get(t);
				int nmcount = 0;
				for (int s = 0; s < nsites; s++) {
					byte base = a.genotype(t, s);
					if (isNotMissing.fastGet(s) && obsMap.get(a.genotype(t, s)) == null) {
						myLogger.info("null from " + Byte.toString(base));
					}
					if (isNotMissing.fastGet(s)) {
						obs[nmcount] = obsMap.get(a.genotype(t, s));
						pos[nmcount++] = a.chromosomalPosition(s);
					}

				}

			}

			double phom = (1 - probHeterozygous) / 2;
			double[] pTrue = new double[]{phom, .25*probHeterozygous ,.5 * probHeterozygous, .25*probHeterozygous, phom};

			//iterate
			ArrayList<byte[]> bestStates = new ArrayList<byte[]>();
			int[][] previousStateCount = new int[5][3];
			int iter = 0;
			boolean hasNotConverged = true;
			while (iter < maxIterations && hasNotConverged) {
				//apply Viterbi
				myLogger.info("Iteration " + iter++ + " for " + familyName);
				bestStates.clear();
				for (int t = 0; t < ntaxa; t++) {
					tp.setPositions(snpPositions.get(t));
					int nobs = notMissingCount[t];
					if (nobs >= 20) {
						ViterbiAlgorithm va = new ViterbiAlgorithm(nonMissingObs.get(t), tp, ep, pTrue);
						va.calculate();
						bestStates.add(va.getMostProbableStateSequence());
					} else { //do not impute if obs < 20
						myLogger.info("Fewer then 20 observations for " + a.taxa().taxaName(t));
						byte[] states = new byte[nobs];
						byte[] obs = nonMissingObs.get(t);
						for (int i = 0; i < nobs; i++) {
							if (obs[i] == AA) states[i] = 0;
							else if (obs[i] == CC) states[i] = 4;
							else states[i] = 2;
						}
						bestStates.add(states);
					}
				}

				//re-estimate transition probabilities
				int[][] transitionCounts = new int[5][5];
				double[][] transitionProb = new double[5][5];
				for (int t = 0; t < ntaxa; t++) {
					byte[] states = bestStates.get(t);
					for (int s = 1; s < notMissingCount[t]; s++) {
						transitionCounts[states[s-1]][states[s]]++;
					}
				}

				//transition is prob(state2 | state1) = count(cell)/count(row)
				for (int row = 0; row < 5; row++) {
					double rowsum = 0;
					for (int col = 0; col < 5; col++) rowsum += transitionCounts[row][col];
					for (int col = 0; col < 5; col++) transitionProb[row][col] = ((double) transitionCounts[row][col]) / rowsum;
				}
				tp.setTransitionCounts(transitionCounts, chrlength, ntaxa);

				//re-estimate emission probabilities
				int[][] emissionCounts = new int[5][3];
				double[][] emissionProb = new double[5][3];
				for (int t = 0; t < ntaxa; t++) {
					byte[] obs = nonMissingObs.get(t);
					byte[] states = bestStates.get(t);
					for (int s = 0; s < notMissingCount[t]; s++) {
						emissionCounts[states[s]][obs[s]]++;
					}
				}

				//debug - print observation/state counts
				//			StringBuilder strb = new StringBuilder("Imputation counts, rows=states, columns=observations:\n");
				//			for (int[] row:emissionCounts) {
				//				for (int cell:row) {
				//					strb.append(cell).append("\t");
				//				}
				//				strb.append("\n");
				//			}
				//			strb.append("\n");
				//			myLogger.info(strb.toString());

				//check to see if there is a change in the observation/state counts
				hasNotConverged = false;
				for (int r = 0; r < 5; r++) {
					for (int c = 0; c < 3; c++) {
						if (previousStateCount[r][c] != emissionCounts[r][c]) {
							hasNotConverged = true;
							previousStateCount[r][c] = emissionCounts[r][c];
						}
					}
				}

				//emission is prob(obs | state) = count(cell)/count(row)
				double[] rowSums = new double[5];
				double total = 0;
				for (int row = 0; row < 5; row++) {
					double rowsum = 0;
					for (int col = 0; col < 3; col++) rowsum += emissionCounts[row][col];
					for (int col = 0; col < 3; col++) emissionProb[row][col] = ((double) emissionCounts[row][col]) / rowsum;
					rowSums[row] = rowsum;
					total += rowsum;
				}
				ep.setEmissionProbability(emissionProb);

				//re-estimate pTrue
				for (int i = 0; i < 5; i++) {
					pTrue[i
					      ] = rowSums[i] / total;
				}

				//if the model has converged  or if the max iterations has been reached print tables
				if (!hasNotConverged || iter == maxIterations) {
					StringBuilder sb = new StringBuilder("Family ");
					sb.append(familyName).append(", chromosome ").append(a.chromosomeName(0));
					if (iter < maxIterations) {
						sb.append(": EM algorithm converged at iteration ").append(iter).append(".\n");
					} else {
						sb.append(": EM algorithm failed to converge after ").append(iter).append(" iterations.\n");
					}

					//print transition counts
					sb = new StringBuilder("Transition counts from row to column:\n");
					for (int[] row:transitionCounts) {
						for (int cell:row) {
							sb.append(cell).append("\t");
						}
						sb.append("\n");
					}
					sb.append("\n");
					myLogger.info(sb.toString());

					//print transition probabilities
					sb = new StringBuilder("Transition probabilities:\n");
					for (double[] row:transitionProb) {
						for (double cell:row) {
							sb.append(cell).append("\t");
						}
						sb.append("\n");
					}
					sb.append("\n");
					myLogger.info(sb.toString());

					//print observation/state counts
					sb = new StringBuilder("Imputation counts, rows=states, columns=observations:\n");
					for (int[] row:emissionCounts) {
						for (int cell:row) {
							sb.append(cell).append("\t");
						}
						sb.append("\n");
					}
					sb.append("\n");
					myLogger.info(sb.toString());

					//print emission probabilities
					sb = new StringBuilder("Emission probabilities:\n");
					for (double[] row:emissionProb) {
						for (double cell:row) {
							sb.append(cell).append("\t");
						}
						sb.append("\n");
					}
					sb.append("\n");
					myLogger.info(sb.toString());

				}


			}

			GenotypeTableBuilder builder = GenotypeTableBuilder.getTaxaIncremental(a.positions());
			nsites = a.numberOfSites();
			for (int t = 0; t < ntaxa; t++) {
				BitSet hasData = notMissingIndex.get(t);
				byte[] states = bestStates.get(t);
				int stateCount = 0;
				byte[] taxonGeno = new byte[nsites];
				Arrays.fill(taxonGeno, NN);
				for (int s = 0; s < nsites; s++) {
					if (hasData.fastGet(s)) {
						if (states[stateCount] == 0) taxonGeno[s] = AA;
						else if (states[stateCount] < 4) taxonGeno[s] = AC;
						else if (states[stateCount] == 4) taxonGeno[s] = CC;
						stateCount++;
					}
				}
				builder.addTaxon(a.taxa().get(t), taxonGeno);
			}

			return builder.build();
		}

		public static void fillGapsInAlignment(PopulationData popdata) {
			popdata.imputed = fillGapsInImputedAlignment(popdata);		
		}

		public static GenotypeTable fillGapsInImputedAlignment(PopulationData popdata) {
			int ntaxa = popdata.imputed.numberOfTaxa();
			int nsites = popdata.imputed.numberOfSites();
			GenotypeTableBuilder builder = GenotypeTableBuilder.getTaxaIncremental(popdata.imputed.positions());
			for (int t = 0; t < ntaxa; t++) {
				int prevsite = -1;
				byte prevValue = -1;
				byte[] taxonGeno = new byte[nsites];
				Arrays.fill(taxonGeno, NN);
				for (int s = 0; s < nsites; s++) {
					byte val = popdata.imputed.genotype(t, s);
					if (val != NN) {
						if (prevsite == -1) {
							prevsite = s;
							prevValue = val;
							taxonGeno[s] = val;
						} else if(val == prevValue) {
							for (int site = prevsite + 1; site <= s; site++) {
								taxonGeno[site] = prevValue;
							}
							prevsite = s;
						} else {
							prevsite = s;
							prevValue = val;
							taxonGeno[s] = val;
						}
					}
				}
				builder.addTaxon(popdata.imputed.taxa().get(t), taxonGeno);
			}

			return builder.build();		
		}

		public static GenotypeTable convertParentCallsToNucleotides(PopulationData popdata) {
			//set monomorphic sites to major (or only allele) (or not)
			//set polymorphic sites consistent with flanking markers if equal, unchanged otherwise
			//do not change sites that are not clearly monomorhpic or polymorphic

			GenotypeTableBuilder builder = GenotypeTableBuilder.getSiteIncremental(popdata.imputed.taxa());
			BitSet isPopSnp = popdata.snpIndex;
			int nsites = popdata.original.numberOfSites();
			int ntaxa = popdata.imputed.taxa().numberOfTaxa();
			int imputedSnpCount = 0;

			for (int s = 0; s < nsites; s++) {
				if (isPopSnp.fastGet(s)) {
					byte AAcall = popdata.alleleA[s];
					byte CCcall = popdata.alleleC[s];
					byte ACcall;
					if (GenotypeTableUtils.isHeterozygous(AAcall)) ACcall = AAcall;
					else if (GenotypeTableUtils.isHeterozygous(CCcall)) ACcall = CCcall; 
					else ACcall = GenotypeTableUtils.getDiploidValue(AAcall, CCcall);
					byte[] geno = new byte[ntaxa];
					
					for (int t = 0; t < ntaxa; t++) {
						byte parentCall = popdata.imputed.genotype(t, imputedSnpCount);
						if (parentCall == AA) {
							geno[t] = AAcall;
						} else if (parentCall == CC) {
							geno[t] = CCcall;
						} else if (parentCall == AC || parentCall == CA) {
							geno[t] = ACcall;
						} else {
							geno[t] = NN;
						}
					}
					builder.addSite(popdata.imputed.positions().get(imputedSnpCount), geno);
					imputedSnpCount++;
				}
			}

			myLogger.info("Original alignment updated for family " + popdata.name + " chromosome " + popdata.original.chromosomeName(0) + ".\n");
			return builder.build();
		}

		public static OpenBitSet whichSnpsAreFromSameTag(GenotypeTable geno, double minRsq) {
			//if (alignIn instanceof SBitAlignment) sba = (SBitAlignment) alignIn;
			//else sba = SBitAlignment.getInstance(alignIn);

			int nsites = geno.numberOfSites();
			int firstSite = 0;
			OpenBitSet isSelected = new OpenBitSet(nsites);
			isSelected.fastSet(0);
			Chromosome firstSnpChr = geno.chromosome(0);
			int firstSnpPos = geno.chromosomalPosition(0);
			while (firstSite < nsites - 1) {
				int nextSite = firstSite + 1;
				int nextSnpPos = geno.chromosomalPosition(nextSite);
				Chromosome nextSnpChr = geno.chromosome(nextSite);
				while (firstSnpChr.equals(nextSnpChr) && nextSnpPos - firstSnpPos < 64) {
					//calculate r^2 between snps
					BitSet rMj = geno.allelePresenceForAllTaxa(firstSite, WHICH_ALLELE.Major);
					BitSet rMn = geno.allelePresenceForAllTaxa(firstSite, WHICH_ALLELE.Minor);
					BitSet cMj = geno.allelePresenceForAllTaxa(nextSite, WHICH_ALLELE.Major);
					BitSet cMn = geno.allelePresenceForAllTaxa(nextSite, WHICH_ALLELE.Minor);
					int[][] contig = new int[2][2];
					contig[0][0] = (int) OpenBitSet.intersectionCount(rMj, cMj);
					contig[1][0] = (int) OpenBitSet.intersectionCount(rMn, cMj);
					contig[0][1] = (int) OpenBitSet.intersectionCount(rMj, cMn);
					contig[1][1] = (int) OpenBitSet.intersectionCount(rMn, cMn);

					double rsq = calculateRSqr(contig[0][0], contig[0][1], contig[1][0], contig[1][1], 2);
					//if rsq cannot be calculated or rsq is less than the minimum rsq for a snp to be considered highly correlated, select this snp
					if (Double.isNaN(rsq) || rsq < minRsq) {
						isSelected.fastSet(nextSite);
						break;
					}
					nextSite++;
					if (nextSite >= nsites) break;
					nextSnpPos = geno.chromosomalPosition(nextSite);
					nextSnpChr = geno.chromosome(nextSite);
				}
				firstSite = nextSite;
				if (firstSite < nsites) {
					firstSnpChr = nextSnpChr;
					firstSnpPos = nextSnpPos;
					isSelected.fastSet(firstSite);
				}
			}

			return isSelected;
		}

		public static OpenBitSet whichSnpsAreFromSameTag(GenotypeTable geno, BitSet polybits) {
			if (polybits.cardinality() == 0) {
				if (polybits instanceof OpenBitSet) return (OpenBitSet) polybits;
				else return new OpenBitSet(polybits.getBits(), polybits.getNumWords());
			}

			int nsites = geno.numberOfSites();
			int firstSite = 0;
			while (!polybits.fastGet(firstSite)) firstSite++;

			OpenBitSet isSelected = new OpenBitSet(nsites);
			isSelected.fastSet(firstSite);
			int firstSnpPos = geno.chromosomalPosition(firstSite);
			Chromosome firstChr = geno.chromosome(firstSite);
			while (firstSite < nsites - 1) {
				int nextSite = firstSite + 1;
				int nextSnpPos = geno.chromosomalPosition(nextSite);
				Chromosome nextChr = geno.chromosome(nextSite);

				boolean skip = !polybits.fastGet(nextSite) || ( firstChr.equals(nextChr) && nextSnpPos - firstSnpPos < 64 /*&& computeRForMissingness(firstSite, nextSite, geno)  > 0.8*/) ;

				while (skip)  {
					boolean validSite = nextSite < nsites;
					while (validSite && !polybits.fastGet(nextSite)) {
						nextSite++;
						validSite = nextSite < nsites;
					}
					if (validSite) {
						nextSnpPos = geno.chromosomalPosition(nextSite);
						nextChr = geno.chromosome(nextSite);
						int dist = nextSnpPos - firstSnpPos;
						boolean nearbySite = firstChr.equals(nextChr) && dist < 64;

						double r = computeRForMissingness(firstSite, nextSite, geno);

						boolean correlatedSite = Double.isNaN(r) || r  > 0.7;
						if (nearbySite && correlatedSite) {
							skip = true;
							nextSite++;
							validSite = nextSite < nsites;
						} else skip = false;

					} else skip = false;
				}

				firstSite = nextSite;
				if (firstSite < nsites) {
					firstChr = nextChr;
					firstSnpPos = nextSnpPos;
					isSelected.fastSet(firstSite);
				}
			}

			return isSelected;
		}

		public static GenotypeTable filterSnpsByTag(GenotypeTable a, double minMaf, double maxMissing, double maxHet) {
			int nOriginalSites = a.numberOfSites();
			int[] selectedSites = new int[nOriginalSites];
			int ntaxa = a.numberOfTaxa();

			//start at the first site with sufficient data
			int headSite = -1;
			double maf, pmissing, phet;
			int npresent, nhet;
			do {
				headSite++;
                                if (headSite >= nOriginalSites) {
                                    myLogger.error("Unable to locate a starting site with MAF > " + minMaf + " and fraction missing < " +
                                    maxMissing + " and fraction heterozygous < " + maxHet);
                                    break;
                                }
				maf = a.minorAlleleFrequency(headSite);
				npresent = a.totalNonMissingForSite(headSite);
				nhet = a.heterozygousCount(headSite);
				pmissing = ((double) (ntaxa - npresent)) /((double) ntaxa);
				if (npresent == 0) phet = 0; 
				else phet = ((double) nhet) / ((double) npresent);
			} while (maf < minMaf || pmissing > maxMissing || phet > maxHet);
				
			selectedSites[0] = headSite;
			int selectedCount = 1;
			for (int s = 1; s < nOriginalSites; s++) {
				int dist = a.chromosomalPosition(s) - a.chromosomalPosition(headSite);
				maf = a.minorAlleleFrequency(s);
				npresent = a.totalNonMissingForSite(s);
				nhet = a.heterozygousCount(s);
				pmissing = ((double) (ntaxa - npresent)) /((double) ntaxa);
				if (npresent == 0) phet = 0; 
				else phet = ((double) nhet) / ((double) npresent);
				if ((dist >= 64 || computeRForMissingness(headSite, s, a) < 0.7) && maf >= minMaf && pmissing <= maxMissing && phet <= maxHet) {
					selectedSites[selectedCount++] = s;
					headSite = s;
				} 
			}
			selectedSites = Arrays.copyOf(selectedSites, selectedCount);
			return FilterGenotypeTable.getInstance(a, selectedSites);

		}

		static double calculateRSqr(int countAB, int countAb, int countaB, int countab, int minTaxaForEstimate) {
			//this is the Hill & Robertson measure as used in Awadella Science 1999 286:2524
			double freqA, freqB, rsqr, nonmissingSampleSize;
			nonmissingSampleSize = countAB + countAb + countaB + countab;
			if (nonmissingSampleSize < minTaxaForEstimate) {
				return Double.NaN;
			}
			freqA = (double) (countAB + countAb) / nonmissingSampleSize;
			freqB = (double) (countAB + countaB) / nonmissingSampleSize;

			//Through missing data & incomplete datasets some alleles can be fixed this returns missing value
			if ((freqA == 0) || (freqB == 0) || (freqA == 1) || (freqB == 1)) {
				return Double.NaN;
			}

			rsqr = ((double) countAB / nonmissingSampleSize) * ((double) countab / nonmissingSampleSize);
			rsqr -= ((double) countaB / nonmissingSampleSize) * ((double) countAb / nonmissingSampleSize);
			rsqr *= rsqr;
			rsqr /= freqA * (1 - freqA) * freqB * (1 - freqB);
			return rsqr;
		}

		public static BitSet ldfilter(GenotypeTable a, int window, double minR, BitSet filterBits) {
			int nsites = a.numberOfSites();
			OpenBitSet ldbits = new OpenBitSet(nsites);
			for (int s = 0; s < nsites; s++) {
				if (filterBits.fastGet(s)) {
					double avgr = neighborLD(a, s, window, filterBits);
					if (avgr >= minR) {
						ldbits.fastSet(s);
					}
				}
			}
			return ldbits;
		}

		public static double neighborLD(GenotypeTable a, int snp, int window, BitSet filterBits) {
			double sumr = 0;
			int nsites = a.numberOfSites();

			//test next <window> sites or to the end of the chromosome
			int site = snp + 10;
			int siteCount = 0;
			int sitesTestedCount = 0;
			while (siteCount < window && site < nsites) {
				if (filterBits.fastGet(site)) {
					double r = computeGenotypeR(snp, site, a);
					if (!Double.isNaN(r)) {
						sumr += Math.abs(computeGenotypeR(snp, site, a));
						siteCount++;
						sitesTestedCount++;

					}
				}
				site++;
			}

			//test previous <window> sites or to the beginning of the chromosome
			site = snp - 10;
			siteCount = 0;
			while (siteCount < window && site >= 0) {
				if (filterBits.fastGet(site)) {
					sumr += Math.abs(computeGenotypeR(snp, site, a));
					siteCount++;
					sitesTestedCount++;
				}
				site--;
			}

			return sumr / sitesTestedCount;
		}

		public static ArrayList<HaplotypeCluster> clusterWindow(GenotypeTable a, int start, int length, int maxdif) {
			int ntaxa = a.numberOfTaxa();
			int end = start + length;
			ArrayList<Haplotype> haps = new ArrayList<Haplotype>();

			for (int t = 0; t < ntaxa; t++) {
				haps.add(new Haplotype(a.genotypeRange(t, start, end), t));
			}

			HaplotypeClusterer hc = new HaplotypeClusterer(haps);
			hc.makeClusters();
			ArrayList<HaplotypeCluster> mergedClusters = HaplotypeClusterer.getMergedClusters(hc.getClusterList(), maxdif);

			//get sequential haplotypes
			//		System.out.println("Sequential haplotypes:");
			//		hc.sortClusters();
			//
			//		for (int i = 0; i < 5; i++) {
			//			ArrayList<HaplotypeCluster> clusters = hc.getClusterList();
			//			if (clusters.size() > 0) {
			//				double score = clusters.get(0).getScore();
			//				String hapstring = clusters.get(0).getHaplotypeAsString();
			//				HaplotypeCluster hapclust = hc.removeFirstHaplotypes(2);
			//				System.out.printf("Count = %d, score = %1.2f: %s\n",hapclust.getSize(), score, hapstring);
			//			}
			//		}

			return mergedClusters;
		}

		/**
		 * Checks whether alleleA is parent1 and alleleC is parent2. If not, switches them.
		 * @param popdata
		 */
		public static void checkParentage(PopulationData popdata) {
			myLogger.info(String.format("Checking parents of %s", popdata.name));
			//		IdGroup ids = popdata.original.getIdGroup();
			TaxaList myTaxaList = popdata.original.taxa();
			LinkedList<Integer> p1list = new LinkedList<Integer>();
			LinkedList<Integer> p2list = new LinkedList<Integer>();
			int nsites = popdata.original.numberOfSites();
			int ntaxa = popdata.original.numberOfTaxa();

			for (int t = 0; t < ntaxa; t++) {
				String tname = myTaxaList.taxaName(t);
				if (tname.toUpperCase().contains(popdata.parent1.toUpperCase())) p1list.add(t);
				if (tname.toUpperCase().contains(popdata.parent2.toUpperCase())) p2list.add(t);
			}

			int consistent = 0;
			int notConsistent = 0;

			for (Integer tndx : p1list) {
				int total = 0;
				int amatches = 0;
				int cmatches = 0;
				int other = 0;
				for (int s = 0; s < nsites; s++) {
					if (popdata.snpIndex.fastGet(s)) {
						byte thisAllele = popdata.original.genotype(tndx, s);
						if (thisAllele != NN) {
							total++;
							if (thisAllele == popdata.alleleA[s]) amatches++;
							else if (thisAllele == popdata.alleleC[s]) cmatches++;
							else other++;
						}
					}
				}
				consistent += amatches;
				notConsistent += cmatches;
				myLogger.info(String.format("%s: total = %d, amatches = %d, cmatches = %d, other = %d", myTaxaList.taxaName(tndx), total, amatches, cmatches, other));
			}

			for (Integer tndx : p2list) {
				int total = 0;
				int amatches = 0;
				int cmatches = 0;
				int other = 0;
				for (int s = 0; s < nsites; s++) {
					if (popdata.snpIndex.fastGet(s)) {
						byte thisAllele = popdata.original.genotype(tndx, s);
						if (thisAllele != NN) {
							total++;
							if (thisAllele == popdata.alleleA[s]) amatches++;
							else if (thisAllele == popdata.alleleC[s]) cmatches++;
							else other++;
						}
					}
				}
				consistent += cmatches;
				notConsistent += amatches;
				myLogger.info(String.format("%s: total = %d, amatches = %d, cmatches = %d, other = %d", myTaxaList.taxaName(tndx), total, amatches, cmatches, other));
			}

			double probNotConsistent = 0;
			if (notConsistent > 0) {
				probNotConsistent = ((double) notConsistent) / ((double) (notConsistent + consistent));
				if (probNotConsistent > 0.75) {
					byte[] temp = popdata.alleleA;
					popdata.alleleA = popdata.alleleC;
					popdata.alleleC = temp; 
				}
			}
			myLogger.info(String.format("Finished checking parents. probNotConsistent = " + probNotConsistent));
		}
	}


