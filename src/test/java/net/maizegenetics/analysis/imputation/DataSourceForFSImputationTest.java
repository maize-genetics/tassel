package net.maizegenetics.analysis.imputation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import net.maizegenetics.analysis.filter.FilterTaxaBuilderPlugin;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

public class DataSourceForFSImputationTest {

	private GenotypeTable myGenotypeTable;
	private GenotypeTable myGenotypeTableWithError;
	private GenotypeTable myGenotypeTableWithErrorNucleotides;
	private ArrayList<byte[]> parentGenotypes;
	private ArrayList<PopulationData> popdata;
	private double percentMissing = 0.5;
	private double homhomRate = .001;
	private double homhetRate = .001;
	private double hethomRate = 0.8;
	private int ntaxa = 400;
	private int nsites = 5000;
	private int chrLength = 100000000;
	private int xoPerChr = 2;
	private int nSelfs = 5;
	private int numberOfFamilies = 2;
	private int[] taxonFamilyIndex;
	private static long randomSeed = 12345;
	private static Random ran;
	
	public static void testcode() {
		long start = System.currentTimeMillis();
		System.out.println("generating data.");
		DataSourceForFSImputationTest myData = new DataSourceForFSImputationTest();
		myData.simulateGenotypeData();
		System.out.printf("data generated in %d\n", System.currentTimeMillis() - start);
		
		GenotypeTable gt = myData.getGenotypes();
		GenotypeTable gterr = myData.getGenotypesWithError();
		System.out.printf("Original: number of taxa = %d, number of sites = %d\n", gt.numberOfTaxa(), gt.numberOfSites());
		System.out.printf("With Error: number of taxa = %d, number of sites = %d\n", gterr.numberOfTaxa(), gterr.numberOfSites());
		
		int notMissing = 0;
		int notMissingErr = 0;
		int hetCount = 0;
		int hetCountErr = 0;
		for (int t = 0; t < myData.ntaxa; t++) {
			notMissing += gt.totalNonMissingForTaxon(t);
			notMissingErr += gterr.totalNonMissingForTaxon(t);
			hetCount += gt.heterozygousCountForTaxon(t);
			hetCountErr += gterr.heterozygousCountForTaxon(t);
		}
		System.out.printf("percent not Missing = %1.3f\n", ((double) notMissingErr)/((double) notMissing));
		System.out.printf("count of het orig = %d, error = %d\n", hetCount, hetCountErr);
		ExportUtils.writeToHapmap(myData.myGenotypeTable, "/Volumes/Macintosh HD 2/temp/simGenotype");
		ExportUtils.writeToHapmap(myData.myGenotypeTableWithError, "/Volumes/Macintosh HD 2/temp/simGenotypeWithErr");
		ExportUtils.writeToHapmap(myData.myGenotypeTableWithErrorNucleotides, "/Volumes/Macintosh HD 2/temp/simGenotypeWithErrNuc");

		errorRate(myData.myGenotypeTable, myData.myGenotypeTableWithError);
	}
	
	public DataSourceForFSImputationTest() {
		ran = new Random(randomSeed);
	}
	
	public GenotypeTable selfUniformRecombination(GenotypeTable geno, double recombinationRate) {
		//recombinationRate is the average per base pair recombination rate
		
		//Haldane: c = [1 - exp(-2x)]/2
		GenotypeTableBuilder builder = GenotypeTableBuilder.getTaxaIncremental(geno.positions());
		int ntaxa = geno.numberOfTaxa();
		int nsites = geno.numberOfSites();
		PositionList posList = geno.positions();
		TaxaList myTaxa = geno.taxa();
		
		//for each site calculate the probability of observing a xo between that site and the previous site
		double[] probR = new double[nsites];
		int prevPos = posList.get(0).getPosition();
		for (int s = 1; s < nsites; s++) {
			int pos = posList.get(s).getPosition();
			double dist = pos - prevPos;
			probR[s] = (1 - Math.exp(-2 * dist * recombinationRate)) / 2;
			prevPos = pos;
		}
		
		for (int t = 0; t < ntaxa; t++) {
			byte[] oldTaxonGeno = geno.genotypeAllSites(t);
			byte[] newTaxonGeno = new byte[nsites];
			
			//generate two gametes
			int g0 = ran.nextInt(2);
			int g1 = ran.nextInt(2);
			for (int s = 0; s < nsites; s++) {
				if (ran.nextDouble() < probR[s]) {
					if (g0 == 0) g0 = 1;
					else g0 = 0;
				}
				if (ran.nextDouble() < probR[s]) {
					if (g1 == 0) g1 = 1;
					else g1 = 0;
				}
				byte[] val = GenotypeTableUtils.getDiploidValues(oldTaxonGeno[s]);
				newTaxonGeno[s] = GenotypeTableUtils.getDiploidValue(val[g0], val[g1]);
			}
			builder.addTaxon(myTaxa.get(t), newTaxonGeno);
		}
		
		return builder.build();
	}
	
	public GenotypeTable makeF1s(int ntaxa, int nsites, int chrlength) {
		PositionListBuilder posBuilder =  new PositionListBuilder();
		Chromosome chr = new Chromosome("George");
		double interval = ((double) chrlength) / ((double) nsites);
		for (int s = 0; s < nsites; s++) {
			int physPos = (int) ((s + 1) * interval);
			Position pos = new GeneralPosition.Builder(chr, physPos).build();
			posBuilder.add(pos);
		}
		GenotypeTableBuilder builder = GenotypeTableBuilder.getTaxaIncremental(posBuilder.build());
		
		byte[] universalGenotype = new byte[nsites];
		Arrays.fill(universalGenotype, GenotypeTableUtils.getDiploidValue(NucleotideAlignmentConstants.A_ALLELE, NucleotideAlignmentConstants.C_ALLELE));
		
		for (int t = 0; t < ntaxa; t++) {
			String name = "t" + t;
			Taxon taxon = new Taxon(name);
			builder.addTaxon(taxon, universalGenotype);
		}
		
		return builder.build();
	}
	
	public GenotypeTable introduceErrorAndMissing(GenotypeTable inputGenotype) {
		GenotypeTableBuilder builder = GenotypeTableBuilder.getTaxaIncremental(inputGenotype.positions());
		
		int ntaxa = inputGenotype.numberOfTaxa();
		int nsites = inputGenotype.numberOfSites();
		byte A = NucleotideAlignmentConstants.A_ALLELE;
		byte C = NucleotideAlignmentConstants.C_ALLELE;
		byte AA = GenotypeTableUtils.getDiploidValuePhased(A,A);
		byte CC = GenotypeTableUtils.getDiploidValuePhased(C,C);
		byte AC = GenotypeTableUtils.getDiploidValuePhased(A,C);
		byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
		
		for (int t = 0; t < ntaxa; t++) {
			byte[] taxonGeno = inputGenotype.genotypeAllSites(t);
			for (int s = 0; s < nsites; s++) {
				if (ran.nextDouble() < percentMissing) {
					taxonGeno[s] = NN;
				} else {
					boolean isHet = GenotypeTableUtils.isHeterozygous(taxonGeno[s]);
					if (isHet) {
						double randbl = ran.nextDouble();
						if (randbl < hethomRate/2) taxonGeno[s] = AA;
						else if (randbl < hethomRate) taxonGeno[s] = CC;
					} else {
						double randbl = ran.nextDouble();
						if (randbl < homhomRate) {
							if (taxonGeno[s] == AA) taxonGeno[s] = CC;
							else taxonGeno[s] = AA;
						} else if (randbl < homhomRate + homhetRate) taxonGeno[s] = AC;
					}
				}
			}
			builder.addTaxon(inputGenotype.taxa().get(t), taxonGeno);
		}
		return builder.build();
	}
	
	public void simulateGenotypeData(boolean withHets) {
		GenotypeTable gt = makeF1s(ntaxa, nsites, chrLength);
		double recombinationRate = ((double) xoPerChr)/chrLength;
		for (int i = 0; i < nSelfs; i++) gt = selfUniformRecombination(gt, recombinationRate);
		myGenotypeTable = gt;
		createPopulationData(gt);
		
		if (withHets) {
			createParentGenotypesWithHet(1000, 2000);
			myGenotypeTableWithError = introduceErrorAndMissing(gt);
			convertGenotypesWithErrorToNucleotidesResidualHets();
		} else {
			createParentGenotypes();
			myGenotypeTableWithError = introduceErrorAndMissing(gt);
			convertGenotypesWithErrorToNucleotides();
		}
	}

	public void simulateGenotypeData() {
		simulateGenotypeData(false);
	}
	
	public void createPopulationData(GenotypeTable gt) {
		popdata = new ArrayList<PopulationData>();
		int ntaxa = gt.numberOfTaxa();
		taxonFamilyIndex = new int[ntaxa];
		TaxaList myTaxa = gt.taxa();
		int numberTaxaPerFamily = ntaxa/numberOfFamilies;
		for (int i = 0; i < numberOfFamilies; i++) {
			PopulationData family = new PopulationData();
			family.name = "family" + i;
			family.members = new ArrayList<String>();
			family.parent1 = "parent1";
			family.parent2 = "parent2";
			family.contribution1 = 0.5;
			family.contribution2 = 0.5;
			family.inbredCoef = 1 - Math.pow(0.5, nSelfs);
			int firstTaxon = i * numberTaxaPerFamily;
			int lastTaxon;
			if (i == numberOfFamilies - 1)  lastTaxon = ntaxa;
			else lastTaxon = firstTaxon + numberTaxaPerFamily;
			
			TaxaListBuilder familyTaxaBuilder = new TaxaListBuilder();
			familyTaxaBuilder.addAll(myTaxa.subList(firstTaxon, lastTaxon));
			family.original = FilterGenotypeTable.getInstance(myGenotypeTable, familyTaxaBuilder.build());
			
			for (int t = firstTaxon; t < lastTaxon; t++) {
				family.members.add(myTaxa.taxaName(t));
				taxonFamilyIndex[t] = i;
			}
			
			popdata.add(family);
		}
	}
	
	public void createParentGenotypes() {
		parentGenotypes = new ArrayList<byte[]>();
		ArrayList<Byte> nukeList = new ArrayList<Byte>();
		nukeList.add(NucleotideAlignmentConstants.A_ALLELE);
		nukeList.add(NucleotideAlignmentConstants.C_ALLELE);
		nukeList.add(NucleotideAlignmentConstants.G_ALLELE);
		nukeList.add(NucleotideAlignmentConstants.T_ALLELE);

		for (int f = 0; f < numberOfFamilies; f++) {
			byte[] p1 = new byte[nsites];
			byte[] p2 = new byte[nsites];
			for (int s = 0; s < nsites; s++) {
				Collections.shuffle(nukeList);
				byte n0 = nukeList.get(0);
				byte n1 = nukeList.get(1);
				p1[s] = GenotypeTableUtils.getDiploidValuePhased(n0, n0);
				p2[s] = GenotypeTableUtils.getDiploidValuePhased(n1, n1);
			}
			parentGenotypes.add(p1);
			parentGenotypes.add(p2);
		}
	}
	
	public void createParentGenotypesWithHet(int startSite, int endSite) {
		parentGenotypes = new ArrayList<byte[]>();
		ArrayList<Byte> nukeList = new ArrayList<Byte>();
		nukeList.add(NucleotideAlignmentConstants.A_ALLELE);
		nukeList.add(NucleotideAlignmentConstants.C_ALLELE);
		nukeList.add(NucleotideAlignmentConstants.G_ALLELE);
		nukeList.add(NucleotideAlignmentConstants.T_ALLELE);

		for (int f = 0; f < numberOfFamilies; f++) {
			byte[] p1 = new byte[nsites];
			byte[] p2 = new byte[nsites];
			for (int s = 0; s < nsites; s++) {
				Collections.shuffle(nukeList);
				byte n0 = nukeList.get(0);
				byte n1 = nukeList.get(1);
				p1[s] = GenotypeTableUtils.getDiploidValuePhased(n0, n0);
				if (s >= startSite && s < endSite && ran.nextBoolean()) {
					if (ran.nextBoolean()) p2[s] = GenotypeTableUtils.getDiploidValuePhased(n0, n1);
					else p2[s] = GenotypeTableUtils.getDiploidValuePhased(n1, n0);
				} else p2[s] = GenotypeTableUtils.getDiploidValuePhased(n1, n1);
			}
			parentGenotypes.add(p1);
			parentGenotypes.add(p2);
		}
	}
	
	public void convertGenotypesWithErrorToNucleotides() {
		byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
		byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
		byte AC = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
		byte CA = NucleotideAlignmentConstants.getNucleotideDiploidByte("CA");
		GenotypeTableBuilder builder = GenotypeTableBuilder.getTaxaIncremental(myGenotypeTableWithError.positions());
		TaxaList myTaxa = myGenotypeTableWithError.taxa();
		for (int t = 0; t < ntaxa; t++) {
			int firstParent = 2 * taxonFamilyIndex[t];
			byte[] parentA = parentGenotypes.get(firstParent);
			byte[] parentC = parentGenotypes.get(firstParent + 1);
			byte[] taxonGeno = myGenotypeTableWithError.genotypeAllSites(t);
			for (int s = 0; s < nsites; s++) {
				if (taxonGeno[s] == AA) taxonGeno[s] = GenotypeTableUtils.getDiploidValue(parentA[s], parentA[s]);
				else if (taxonGeno[s] == CC) taxonGeno[s] = GenotypeTableUtils.getDiploidValue(parentC[s], parentC[s]);
				else if (taxonGeno[s] == AC) taxonGeno[s] = GenotypeTableUtils.getDiploidValue(parentA[s], parentC[s]);
				else if (taxonGeno[s] == CA) taxonGeno[s] = GenotypeTableUtils.getDiploidValue(parentC[s], parentA[s]);
			}
			builder.addTaxon(myTaxa.get(t), taxonGeno);
		}
		myGenotypeTableWithErrorNucleotides = builder.build();
	}
	
	public void convertGenotypesWithErrorToNucleotidesResidualHets() {
		byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
		byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
		byte AC = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
		byte CA = NucleotideAlignmentConstants.getNucleotideDiploidByte("CA");
		GenotypeTableBuilder builder = GenotypeTableBuilder.getTaxaIncremental(myGenotypeTableWithError.positions());
		TaxaList myTaxa = myGenotypeTableWithError.taxa();
		for (int t = 0; t < ntaxa; t++) {
			int whichGrandParent;
			if (ran.nextBoolean()) whichGrandParent = 0;
			else whichGrandParent = 1;
			int firstParent = 2 * taxonFamilyIndex[t];
			byte[] parentA = parentGenotypes.get(firstParent);
			byte[] parentC = parentGenotypes.get(firstParent + 1);
			byte[] taxonGeno = myGenotypeTableWithError.genotypeAllSites(t);
			for (int s = 0; s < nsites; s++) {
				if (taxonGeno[s] == AA) taxonGeno[s] = parentA[s];
				else {
					byte Cval = GenotypeTableUtils.getDiploidValues(parentC[s])[whichGrandParent];
					if (taxonGeno[s] == CC) taxonGeno[s] = GenotypeTableUtils.getDiploidValue(Cval, Cval);
					else if (taxonGeno[s] == AC) taxonGeno[s] = GenotypeTableUtils.getDiploidValue(parentA[s], Cval);
					else if (taxonGeno[s] == CA) taxonGeno[s] = GenotypeTableUtils.getDiploidValue(Cval, parentA[s]);
				}
				
			}
			builder.addTaxon(myTaxa.get(t), taxonGeno);
		}
		myGenotypeTableWithErrorNucleotides = builder.build();
	}
	
	public static int[] errorRate(GenotypeTable refTable, GenotypeTable compTable) {
		//assumes taxa and sites are the same for both tables
		//calculations only include sites not missing in both tables
		//calculate homHomErr = number of wrong hom/total hom in ref 
		//	homHetErr = number of hom in ref scored as het in comp/total hom in ref
		//	hetHomErr = number of het in ref scored as hom in comp/total het in ref
		
		int ntaxa = refTable.numberOfTaxa();
		int nHomHom = 0;
		int nHomHet = 0;
		int nHetHom = 0;
		int nHomRef = 0;
		int nHetRef = 0;
		
		//are the major alleles the same?
		int nsites = refTable.numberOfSites();
		OpenBitSet majorSame = new OpenBitSet(nsites);
		for (int s = 0; s < nsites; s++) {
			if (refTable.majorAllele(s) == compTable.majorAllele(s)) majorSame.fastSet(s);
		}
		
		for (int t = 0; t < ntaxa; t++) {
			BitSet mjref = refTable.allelePresenceForAllSites(t, WHICH_ALLELE.Major);
			BitSet mnref = refTable.allelePresenceForAllSites(t, WHICH_ALLELE.Minor);
			BitSet mjcomp = compTable.allelePresenceForAllSites(t, WHICH_ALLELE.Major);
			BitSet mncomp = compTable.allelePresenceForAllSites(t, WHICH_ALLELE.Minor);
			
			OpenBitSet presentBoth = new OpenBitSet(mjref);
			presentBoth.union(mnref);
			OpenBitSet presentComp = new OpenBitSet(mjcomp);
			presentComp.union(mncomp);
			presentBoth.intersect(presentComp);
			
			//homozygotes
			OpenBitSet homref = new OpenBitSet(mjref);
			homref.xor(mnref);
			OpenBitSet homcomp = new OpenBitSet(mjcomp);
			homcomp.xor(mncomp);
			nHomRef += (int) OpenBitSet.intersectionCount(homref, presentBoth);
			
			//heterozygotes
			OpenBitSet hetref = new OpenBitSet(mjref);
			hetref.and(mnref);
			OpenBitSet hetcomp = new OpenBitSet(mjcomp);
			hetcomp.and(mncomp);
			nHetRef += (int) OpenBitSet.intersectionCount(hetref, presentBoth);
			
			//wrong homozygotes
			OpenBitSet wrongHom = new OpenBitSet(mjref);
			wrongHom.and(mjcomp);
			wrongHom.xor(majorSame);
			wrongHom.and(homref);
			wrongHom.and(homcomp);
			nHomHom = (int) wrongHom.cardinality();
			
			//homozygous ref, heterozygous comp
			nHomHet += (int) OpenBitSet.intersectionCount(homref, hetcomp);
			
			//heterozygous ref, homozygous comp
			nHetHom += (int) OpenBitSet.intersectionCount(hetref, homcomp);
		}
		
		System.out.printf("number of hom-hom errors = %d, total hom = %d, ratio = %1.6f\n", nHomHom, nHomRef, ((double) nHomHom)/nHomRef);
		System.out.printf("number of hom-het errors = %d, total hom = %d, ratio = %1.6f\n", nHomHet, nHomRef, ((double) nHomHet)/nHomRef);
		System.out.printf("number of het-hom errors = %d, total het = %d, ratio = %1.6f\n", nHetHom, nHetRef, ((double) nHetHom)/nHetRef);
		System.out.printf("percent het in ref = %1.5f\n", ((double) nHetRef) / (nHetRef+nHomRef));
		return new int[]{nHomHom, nHomHet, nHetHom, nHomRef, nHetRef};
	}
	
	public void setPercentMissing(double percentMissing) {
		this.percentMissing = percentMissing;
	}

	public void setHomhomRate(double homhomRate) {
		this.homhomRate = homhomRate;
	}

	public void setHomhetRate(double homhetRate) {
		this.homhetRate = homhetRate;
	}

	public void setHethomRate(double hethomRate) {
		this.hethomRate = hethomRate;
	}

	public void setNtaxa(int ntaxa) {
		this.ntaxa = ntaxa;
	}

	public void setNsites(int nsites) {
		this.nsites = nsites;
	}

	public void setXoPerChr(int xoPerChr) {
		this.xoPerChr = xoPerChr;
	}

	public void setnSelfs(int nSelfs) {
		this.nSelfs = nSelfs;
	}

	public GenotypeTable getGenotypes() {
		return myGenotypeTable;
	}

	public GenotypeTable getGenotypes(int family) {
		TaxaListBuilder taxaInFamily = new TaxaListBuilder();
		TaxaList allTaxa = myGenotypeTable.taxa();
		for (int i = 0; i < ntaxa; i++) {
			if (taxonFamilyIndex[i] == family) taxaInFamily.add(allTaxa.get(i));
		}
                return new FilterTaxaBuilderPlugin().taxaList(taxaInFamily.build()).runPlugin(myGenotypeTable);
	}
	
	public GenotypeTable getGenotypesWithError() {
		return myGenotypeTableWithError;
	}
	
	public ArrayList<PopulationData> getPopulationData() {
		return popdata;
	}

	public void setNumberOfFamilies(int numberOfFamilies) {
		this.numberOfFamilies = numberOfFamilies;
	}

	public GenotypeTable getMyGenotypeTableWithErrorNucleotides() {
		return myGenotypeTableWithErrorNucleotides;
	}

	public int getNumberOfFamilies() {
		return numberOfFamilies;
	}

	/**
	 * @param family	zero based family index
	 * @param parent	zero based parent index. Valid values are 0 and 1.
	 * @return	the genotype for this parent and family used to generate myGenotypeTableWithErrorNucleotides
	 */
	public byte[] getParentGenotype(int family, int parent) {
		return parentGenotypes.get(family * 2 + parent);
	}
	
	public void setRandomSeed(long seed) {
		randomSeed = seed;
		ran = new Random(seed);
	}
}
