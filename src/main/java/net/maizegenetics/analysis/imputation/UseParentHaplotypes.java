package net.maizegenetics.analysis.imputation;

import java.util.OptionalDouble;
import java.util.stream.IntStream;

import net.maizegenetics.analysis.popgen.LDResult;
import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.stats.statistics.FisherExact;
import net.maizegenetics.util.BitSet;

public class UseParentHaplotypes {
    PopulationData myFamily;
    double minMaf = 0.1;
    double minCoverage = 0.5;
    double maxHet = 0.15;
    
    public UseParentHaplotypes(PopulationData family) {
        myFamily = family;
    }
    
    public void assignHaplotypes() {
        prefilterSites();
        setAllelesToParents();
        System.out.printf("imputed genotype table has %d sites and %d taxa", myFamily.imputed.numberOfSites(), myFamily.imputed.numberOfTaxa());
        System.out.println();
    }
    
    private void setAllelesToParents() {
        byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
        byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
        byte AC = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
        int parent1Index = myFamily.imputed.taxa().indexOf(myFamily.parent1);
        int parent2Index = myFamily.imputed.taxa().indexOf(myFamily.parent2);
        myFamily.alleleA = myFamily.imputed.genotypeAllSites(parent1Index);
        myFamily.alleleC = myFamily.imputed.genotypeAllSites(parent2Index);
        byte[] genoA = myFamily.alleleA;
        byte[] genoC = myFamily.alleleC;
        int nsites = myFamily.imputed.numberOfSites();
        int ntaxa = myFamily.imputed.numberOfTaxa();
        
        //if any sites are missing in A assign the opposite allele
        for (int s = 0; s < nsites; s++) {
            if (genoA[s] == NN && genoC[s] != NN) {
                byte major = myFamily.imputed.majorAllele(s);
                byte majorgeno = GenotypeTableUtils.getDiploidValue(major, major);
                byte minor = myFamily.imputed.minorAllele(s);
                byte minorgeno = GenotypeTableUtils.getDiploidValue(minor, minor);
                if (genoC[s] == majorgeno) {
                    genoA[s] = minorgeno;
                } else if (genoC[s] == minorgeno) {
                    genoA[s] = majorgeno;
                }
            } else if (genoA[s] != NN && genoC[s] == NN) {
                byte major = myFamily.imputed.majorAllele(s);
                byte majorgeno = GenotypeTableUtils.getDiploidValue(major, major);
                byte minor = myFamily.imputed.minorAllele(s);
                byte minorgeno = GenotypeTableUtils.getDiploidValue(minor, minor);
                if (genoA[s] == majorgeno) {
                    genoC[s] = minorgeno;
                } else if (genoA[s] == minorgeno) {
                    genoC[s] = majorgeno;
                }
            }
        }
        
        //code genotypes as the parent alleles
        GenotypeTableBuilder genoBuilder = GenotypeTableBuilder.getTaxaIncremental(myFamily.imputed.positions());
        for (int t = 0; t < ntaxa; t++) {
            byte[] taxonGeno = myFamily.imputed.genotypeAllSites(t);
            for (int s = 0; s < nsites; s++) {
                if (genoA[s] == genoC[s] || GenotypeTableUtils.isHeterozygous(genoA[s]) || GenotypeTableUtils.isHeterozygous(genoC[s])) taxonGeno[s] = NN;
                else {
                    if (taxonGeno[s] == genoA[s]) taxonGeno[s] = AA;
                    else if (taxonGeno[s] == genoC[s]) taxonGeno[s] = CC;
                    else if (GenotypeTableUtils.isHeterozygous(taxonGeno[s])) taxonGeno[s] = AC;
                    else taxonGeno[s] = NN;
                }
            }
            genoBuilder.addTaxon(myFamily.imputed.taxa().get(t), taxonGeno);
        }
        
        myFamily.imputed = genoBuilder.build();
    }
    
    private void validateParentGenotypes() {
        //parent1
        int p1 = myFamily.original.taxa().indexOf(myFamily.parent1);
        byte[] parent1Genotype = myFamily.original.genotypeAllSites(p1);
        
        //parent2
        int p2 = myFamily.original.taxa().indexOf(myFamily.parent1);
        byte[] parent2Genotype = myFamily.original.genotypeAllSites(p2);
        
        
    }
    
    private void prefilterSites() {
        //apply minMaf, minCoverage, and maxHet filters
        final GenotypeTable geno = myFamily.original;
        int nsites = geno.numberOfSites();
        final int minNonMissing = (int) minCoverage * nsites;
        
        int[] sites = IntStream.range(0, nsites).filter(s -> {
            if (geno.minorAlleleFrequency(s) < minMaf) return false;
            int numberNotMissing = geno.totalNonMissingForSite(s);
            if (numberNotMissing < minNonMissing) return false;
            double proportionHet = geno.heterozygousCount(s) / ((double) numberNotMissing);
            if (proportionHet > maxHet) return false;
            return true;
        }).toArray();
        
        GenotypeTable filteredGeno = FilterGenotypeTable.getInstance(geno, sites);
        GenotypeTable copy = GenotypeTableBuilder.getGenotypeCopyInstance(filteredGeno);
        
        //find sites that are not in LD
        
        int[] goodSites = IntStream.range(0,  copy.numberOfSites())
                .filter(s -> isSiteInLD(copy, s))
                .toArray();
       
        if (goodSites.length < sites.length) {
            myFamily.imputed = FilterGenotypeTable.getInstance(copy, goodSites);
        } else {
            myFamily.imputed = filteredGeno;
        }
        
    }
    
    private boolean isSiteInLD(GenotypeTable geno, int site) {
        int window = 50;
        double minr2 = 0.8;
        int nsites = geno.numberOfSites();
        int startSite = Math.max(0, site - window);
        int endSite = Math.min(nsites, site + window + 1);
        BitSet rMj = geno.allelePresenceForAllTaxa(site, WHICH_ALLELE.Major);
        BitSet rMn = geno.allelePresenceForAllTaxa(site, WHICH_ALLELE.Minor);
        FisherExact fisherExact = FisherExact.getInstance((2 * geno.numberOfTaxa()) + 10);
        
        double maxR2 = 0;
        for (int s = startSite; s < endSite; s++) {
            if (s != site) {
                BitSet cMj = geno.allelePresenceForAllTaxa(s, WHICH_ALLELE.Major);
                BitSet cMn = geno.allelePresenceForAllTaxa(s, WHICH_ALLELE.Minor);
                float r2 = LinkageDisequilibrium.getLDForSitePair(rMj, rMn, cMj, cMn, 2, 10, -1.0f, fisherExact, site, s).r2();
                maxR2 = Math.max(maxR2, r2);
            }
        }
        if (maxR2 < minr2) return false;
        return true;
    }
    
}
