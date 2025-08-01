package net.maizegenetics.phenotype;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import net.maizegenetics.analysis.numericaltransform.NumericalGenotypePlugin;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.OpenBitSet;

public class GenotypePhenotypeBuilderTest {

    Phenotype pheno1, pheno2;
    GenotypeTable testGenotype;
    
    @Before
    public void setUp() throws Exception {
        List<Taxon> taxaList = new ArrayList<Taxon>();
        taxaList.add(new Taxon("t1"));
        taxaList.add(new Taxon("t2"));
        taxaList.add(new Taxon("t3"));
        taxaList.add(new Taxon("t4"));
        taxaList.add(new Taxon("t5"));
        taxaList.add(new Taxon("t6"));
        
        TaxaAttribute ta = new TaxaAttribute(taxaList);
        NumericAttribute trait = new NumericAttribute("trait1", new float[]{1,2,3,4,5,6}, new OpenBitSet(6));
        
        List<PhenotypeAttribute> attributes = new ArrayList<PhenotypeAttribute>();
        attributes.add(ta);
        attributes.add(trait);
        List<ATTRIBUTE_TYPE> types = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
        types.add(ATTRIBUTE_TYPE.taxa);
        types.add(ATTRIBUTE_TYPE.data);
        
        pheno1 = new CorePhenotype(attributes, types, "pheno1");
        
        List<Taxon> taxaList2 = new ArrayList<Taxon>();
        taxaList2.add(new Taxon("t1"));
        taxaList2.add(new Taxon("t3"));
        taxaList2.add(new Taxon("t5"));
        taxaList2.add(new Taxon("t2"));
        taxaList2.add(new Taxon("t4"));
        taxaList2.add(new Taxon("t6"));

        TaxaAttribute ta2 = new TaxaAttribute(taxaList2);
        NumericAttribute trait2 = new NumericAttribute("trait1", new float[]{1,2,3,4,5,6}, new OpenBitSet(6));
        
        List<PhenotypeAttribute> attributes2 = new ArrayList<PhenotypeAttribute>();
        attributes2.add(ta2);
        attributes2.add(trait2);
        List<ATTRIBUTE_TYPE> types2 = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
        types2.add(ATTRIBUTE_TYPE.taxa);
        types2.add(ATTRIBUTE_TYPE.data);
        
        pheno2 = new CorePhenotype(attributes2, types2, "pheno2");
        
        PositionListBuilder poslistBuilder = new PositionListBuilder();
        Chromosome chr = Chromosome.instance(1);
        poslistBuilder.add(new GeneralPosition.Builder(chr, 1).build());
        poslistBuilder.add(new GeneralPosition.Builder(chr, 2).build());
        poslistBuilder.add(new GeneralPosition.Builder(chr, 3).build());
        poslistBuilder.add(new GeneralPosition.Builder(chr, 4).build());
        poslistBuilder.add(new GeneralPosition.Builder(chr, 5).build());
        
        GenotypeTableBuilder genoBuilder = GenotypeTableBuilder.getTaxaIncremental(poslistBuilder.build());
        byte AA = GenotypeTableUtils.getDiploidValue(NucleotideAlignmentConstants.A_ALLELE, NucleotideAlignmentConstants.A_ALLELE);
        byte CC = GenotypeTableUtils.getDiploidValue(NucleotideAlignmentConstants.C_ALLELE, NucleotideAlignmentConstants.C_ALLELE);
        
        genoBuilder.addTaxon(taxaList.get(0), new byte[] {CC,AA,AA,AA,AA});
        genoBuilder.addTaxon(taxaList.get(1), new byte[] {AA,CC,AA,AA,AA});
        genoBuilder.addTaxon(taxaList.get(2), new byte[] {AA,AA,CC,AA,AA});
        genoBuilder.addTaxon(taxaList.get(3), new byte[] {AA,AA,AA,CC,AA});
        genoBuilder.addTaxon(taxaList.get(4), new byte[] {AA,AA,AA,AA,CC});
        genoBuilder.addTaxon(taxaList.get(5), new byte[] {CC,CC,AA,AA,AA});
        
        GenotypeTable gt = genoBuilder.build();
        NumericalGenotypePlugin ngt = new NumericalGenotypePlugin();
        testGenotype = (GenotypeTable) ngt.performFunction(DataSet.getDataSet(gt)).getData(0).getData();
    }

    @Test
    public void testSortedGenotype() {
        GenotypePhenotype sortedGenoPheno = new GenotypePhenotypeBuilder().genotype(testGenotype).phenotype(pheno1).build();
        testGenotype(sortedGenoPheno, new int[] {0,1,2,3,4,5});
    }
    
    @Test
    public void testUnsortedGenotype() {
        GenotypePhenotype unsortedGenoPheno = new GenotypePhenotypeBuilder().genotype(testGenotype).phenotype(pheno2).build();
        testGenotype(unsortedGenoPheno, new int[] {0,2,4,1,3,5});
    }
    
    private void testGenotype(GenotypePhenotype genoPheno, int[] taxaIndex) {
        int nobs = genoPheno.numberOfObservations();
        assertEquals("Number of observations", 6, nobs);
        int nsites = genoPheno.genotypeTable().numberOfSites();
        assertEquals("Number of sites", 5, nsites);
        for (int obs = 0; obs < nobs; obs++) {
            for (int s = 0; s < nsites; s++) {
                String label = String.format("Obs = %d, Site = %d", obs, s);
                assertEquals(label, genoPheno.genotype(obs, s), testGenotype.genotype(taxaIndex[obs], s));
            }
        }
    }
    
   
    @Test
    public void testSortedReferenceProbability() {
        GenotypePhenotype sortedGenoPheno = new GenotypePhenotypeBuilder().genotype(testGenotype).phenotype(pheno1).build();
        testReferenceProbability(sortedGenoPheno, new int[] {0,1,2,3,4,5});
    }
    
    @Test
    public void testUnsortedReferenceProbability() {
        GenotypePhenotype unsortedGenoPheno = new GenotypePhenotypeBuilder().genotype(testGenotype).phenotype(pheno2).build();
        testReferenceProbability(unsortedGenoPheno, new int[] {0,2,4,1,3,5});
    }

    private void testReferenceProbability(GenotypePhenotype genoPheno, int[] taxaIndex) {
        int nsites = genoPheno.genotypeTable().numberOfSites();
        int nobs = genoPheno.numberOfObservations();
        for (int s = 0; s < nsites; s++ ) {
            float[] refProbGenoPheno = genoPheno.referenceProb(s);
            for (int obs = 0; obs < nobs; obs++) {
                String label = String.format("Obs = %d, Site = %d", obs, s);
                assertEquals(label, refProbGenoPheno[obs], testGenotype.referenceProbability(taxaIndex[obs], s), 1e-5);
            }
        }
    }
}
