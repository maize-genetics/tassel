package net.maizegenetics.phenotype;

import net.maizegenetics.analysis.filter.FilterTaxaBuilderPlugin;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListUtils;

/**
 * @author Peter Bradbury
 *
 */
public class GenotypePhenotypeBuilder {
	Phenotype myPhenotype = null;
	GenotypeTable myGenotype = null;
	boolean isUnion = false;
	String myName = null;
	
	/**
	 * @param thePhenotype	the Phenotype to be used by this builder
	 * @return	this builder
	 */
	public GenotypePhenotypeBuilder phenotype(Phenotype thePhenotype) {
		myPhenotype = thePhenotype;
		return this;
	}
	
	/**
	 * @param theGenotype	the GenotypeTable to be used by this builder
	 * @return	this builder
	 */
	public GenotypePhenotypeBuilder genotype(GenotypeTable theGenotype) {
		myGenotype = theGenotype;
		return this;
	}
	
	/**
	 * Indicates that a union join should be performed. The GenotypeTable and Phenotype will be used as is.
	 * If union is not specified, an intersect join will be performed by default
	 * @return	this builder
	 */
	public GenotypePhenotypeBuilder union() {
		isUnion = true;
		return this;
	}
	
	/**
	 * Indicates that an intersect join should be performed. If necessary the GenotypeTable and Phenotype will be filtered before building the GenotypePhenotype.
	 * The filtered tables will have only the taxa common to both.
	 * @return	this builder
	 */
	public GenotypePhenotypeBuilder intersect() {
		isUnion = false;
		return this;
	}
	
	/**
	 * @param name	the name that will be used for the resulting GenotypePhenotype
	 * @return	this builder
	 */
	public GenotypePhenotypeBuilder name(String name) {
		myName = name;
		return this;
	}
	
	/**
	 * @return	a GenotypePhenotype built using the supplied objects and parameters
	 */
	public GenotypePhenotype build() {
		if (myPhenotype == null) throw new IllegalArgumentException("Error: no phenotype data set was specified.");
		if (myGenotype == null) throw new IllegalArgumentException("Error: no genotype data set was specified.");
		if (myName == null) {
			myName = myPhenotype.name() + "_with_genotypes";
		}
		
		if (!isUnion) {
			TaxaList commonTaxa = TaxaListUtils.getCommonTaxa(myGenotype.taxa(), myPhenotype.taxa());
			if (myGenotype.taxa().numberOfTaxa() > commonTaxa.numberOfTaxa()) {
                                myGenotype = new FilterTaxaBuilderPlugin().taxaList(commonTaxa).runPlugin(myGenotype);
			}
			if (myPhenotype.taxa().numberOfTaxa() > commonTaxa.numberOfTaxa()) {
				myPhenotype = new PhenotypeBuilder().fromPhenotype(myPhenotype).keepTaxa(commonTaxa).build().get(0);
			}
		}
		
		return new GenotypePhenotype(myGenotype, myPhenotype, myName);
	}
}
