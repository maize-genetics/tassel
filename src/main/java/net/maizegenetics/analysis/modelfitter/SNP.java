package net.maizegenetics.analysis.modelfitter;

import java.util.ArrayList;

import net.maizegenetics.dna.map.Chromosome;

public class SNP {
	public String name;
	public Chromosome locus;
	public int position;
	public ArrayList<Object> alleles;
	public int index;
	
	public SNP(String name, Chromosome locus, int positionInLocus, int index) {
		this.name = name;
		this.locus = locus;
		position = positionInLocus;
		this.index = index;
	}
	
	public SNP() {
		
	}

	@Override
	public String toString() {
		return name;
	}
	
	
}
