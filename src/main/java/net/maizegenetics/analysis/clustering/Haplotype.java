package net.maizegenetics.analysis.clustering;

import java.util.Arrays;

import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

/**
 * A Haplotype is a sequence which can be either a true haplotype or a genotype.
 * The natural ordering is by number of non-missing values, descending.
 * @author Peter Bradbury
 */
public class Haplotype implements Comparable<Haplotype>{
	//byte values for hets
	private static byte AC = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
	private static byte AG = NucleotideAlignmentConstants.getNucleotideDiploidByte("AG");
	private static byte AT = NucleotideAlignmentConstants.getNucleotideDiploidByte("AT");
	private static byte CG = NucleotideAlignmentConstants.getNucleotideDiploidByte("CG");
	private static byte CT = NucleotideAlignmentConstants.getNucleotideDiploidByte("CT");
	private static byte GT = NucleotideAlignmentConstants.getNucleotideDiploidByte("GT");
	private static byte CA = NucleotideAlignmentConstants.getNucleotideDiploidByte("CA");
	private static byte GA = NucleotideAlignmentConstants.getNucleotideDiploidByte("GA");
	private static byte TA = NucleotideAlignmentConstants.getNucleotideDiploidByte("TA");
	private static byte GC = NucleotideAlignmentConstants.getNucleotideDiploidByte("GC");
	private static byte TC = NucleotideAlignmentConstants.getNucleotideDiploidByte("TC");
	private static byte TG = NucleotideAlignmentConstants.getNucleotideDiploidByte("TG");
	
	/**
	 * The haplotype (or genotype) sequence
	 */
	public byte[] seq;
	
	/**
	 * The number of non missing values in the sequence
	 */
	public int notMissingCount;
	
	/**
	 * The length of the sequence
	 */
	public int seqlen;
	
	/**
	 * The byte value for missing
	 */
	public static final byte N = NucleotideAlignmentConstants.getNucleotideDiploidByte("N");
	
	/**
	 * The index of the taxon from the originating alignment from which this sequence was taken
	 */
	public int taxonIndex;

	/**
	 * @param hap	the haplotype or genotype sequence
	 * @param taxon	the taxon index
	 */
	public Haplotype(byte[] hap, int taxon) {
		seq = makeHetsConsistent(hap);
		seqlen = seq.length;
		
		//convert heterozygotes to consistent value
		countNotMissing();
		taxonIndex = taxon;
	}
	
	/**
	 * @param hap	the haplotype or genotype sequence
	 */
	public Haplotype(byte[] hap) {
		this(hap, -1);
	}
	
	public byte[] makeHetsConsistent(byte[] seqIn) {
		int n = seqIn.length;
		byte[] seqOut = new byte[n];
		for (int i = 0; i < n; i++) {
			byte val = seqIn[i];
			if (val == CA) val = AC;
			else if (val == GA) val = AG;
			else if (val == TA) val = AT;
			else if (val == GC) val = CG;
			else if (val == TC) val = CT;
			else if (val == TG) val = GT;
			seqOut[i] = val;
		}
		return seqOut;
	}
	
	/**
	 * 
	 */
	public void countNotMissing() {
		notMissingCount = 0;
		for (byte b:seq) if (b != N) notMissingCount++;
	}

	@Override
	public int compareTo(Haplotype arg0) {
		if (notMissingCount > arg0.notMissingCount) return -1;
		if (notMissingCount < arg0.notMissingCount) return 1;
		return 0;
	}
	
	public int numberOfHets() {
		int nhets = 0;
		for (byte b : seq) {
			if (GenotypeTableUtils.isHeterozygous(b)) nhets++;
		}
		return nhets;
	}
	
	/**
	 * Distance is defined as the sum of the distances between each pair of sites. If one or both of the sites are missing then
	 * the distance is zero. If the sites are equal the distance is zero. If both sites are homozygous but different, the distance is 2.
	 * If one site is heterozygous, the distance is 1.
	 * @param h0	another Haplotype. Both haplotypes must have the same sequence length.
	 * @return	the distance between the haplotypes 
	 */
	public int distanceFrom(Haplotype h0) {
		return getDistance(seq, h0.seq);
	}

	public Haplotype subHaplotype(int size, boolean fromStart) {
		byte[] newseq;
		if (fromStart) newseq = Arrays.copyOf(seq, size);
		else newseq = Arrays.copyOfRange(seq, seqlen - size, seqlen);
		return new Haplotype(newseq, taxonIndex);
	}
	

	/**
	 * @param hap0	the first Haplotype
	 * @param hap1	the second Haplotype
	 * @return distance as defined by distanceFrom(Haplotype h0)
	 */
	public static int getDistance(Haplotype hap0, Haplotype hap1) {
		return hap0.distanceFrom(hap1);
	}
	
	/**
	 * @param hap0	the first Haplotype
	 * @param hap1	the second Haplotype
	 * @return distance as defined by distanceFrom(Haplotype h0)
	 * @see #distanceFrom(Haplotype h0)
	 */
	public static int getDistance(byte[] hap0, byte[] hap1) {
		int d = 0;
		int n = hap0.length;
		for (int s = 0; s < n; s++) {
			byte b0 = hap0[s];
			byte b1 = hap1[s];
			if (b0 != b1 && b0 != N && b1 != N) {
				if (GenotypeTableUtils.isHeterozygous(b0)) {
					if (!GenotypeTableUtils.isHeterozygous(b1)) d++;
				} else if (GenotypeTableUtils.isHeterozygous(b1)) d++;
				else d += 2;
			}
		}
		return d;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (byte b : seq) sb.append(NucleotideAlignmentConstants.getNucleotideIUPAC(b));
		return sb.toString();
	}
	
}
