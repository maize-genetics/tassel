package net.maizegenetics.analysis.modelfitter;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.List;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;

/**
 * GenotypeAdditiveSite is an AdditiveSite that takes a byte[] array of genotypes as input and converts them
 * to a double[] array. It codes the homozygous major as 2, homozygous minor as 0 and het as 1. It then subtracts the mean
 * value from all elements of the array and sets missing to zero (now the mean). For efficient storage, the values of 0,1,2,3
 * are stored in a byte array which indexes the double values.
 */
public class GenotypeAdditiveSite extends AbstractAdditiveSite {

    private static final long serialVersionUID = -7891486608129027827L;
    private int[] bitStore;
    private double[] byteConversion;
    protected int ntaxa;
    protected int[] taxaIndex = null;

    /**
     * @param site						the site index from the originating GenotypeTable
     * @param genotype					genotypeAllTaxa for this site
     * @param majorAllele				major allele
     * @param majorAlleleFrequency		major allele frequency
     */
    public GenotypeAdditiveSite(int site, String chr, int pos, String id,
            CRITERION selectionCriterion, byte[] genotype,
            byte majorAllele, double majorAlleleFrequency) {
        super(site, chr, pos, id, selectionCriterion);
        byte unknown = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        ntaxa = genotype.length;
        double mean = 2 * majorAlleleFrequency;
//        byteConversion = new double[] { -mean, 1 - mean, 2 - mean, 0 };
        byteConversion = new double[] { 0,1,2,mean };  //changed for AddPlusDom

        int numberOfInts = ntaxa / 16;
        int remainder = ntaxa % 16;
        if (remainder > 0)
            numberOfInts++;

        bitStore = new int[numberOfInts];
        int intCount = 0;
        int genoCount = 0;

        while (genoCount < ntaxa) {
            int intStore = 0;
            for (int i = 0; i < 32; i += 2) {
                //calculate genoval
                int genoIndex;
                byte geno = genotype[genoCount++];
                if (geno == unknown) {
                    genoIndex = 3;
                } else {
                    genoIndex = 0;
                    byte[] alleles = GenotypeTableUtils.getDiploidValues(geno);
                    if (alleles[0] == majorAllele)
                        genoIndex++;
                    if (alleles[1] == majorAllele)
                        genoIndex++;
                }

                intStore = intStore | (genoIndex << i);
                if (genoCount == ntaxa)
                    break;
            }

            bitStore[intCount++] = intStore;
        }
    }

    //private constructor for making a copy
    private GenotypeAdditiveSite(int site, String chr, int pos, String id,
                                 CRITERION selectionCriterion, int[] bitStore, int ntaxa, double[] byteConversion) {
        super(site, chr, pos, id, selectionCriterion);
        this.bitStore = bitStore;
        this.ntaxa = ntaxa;
        this.byteConversion = byteConversion;
    }

    @Override
    public void reindexTaxa(int[] taxaIndex, List<Integer> uniqueTaxa) {
        this.taxaIndex = taxaIndex;
        //count 0,1,2,3
        int[] counts = new int[4];
        for (Integer Ndx : uniqueTaxa) {
            counts[genotypeIndex(Ndx)]++;
        }

        double numerator = (double) (2 * counts[2] + counts[1]);
        double denominator = numerator + (counts[1] + 2 * counts[0]);
        double majorAlleleFreq = numerator / denominator;
        double mean = 2 * majorAlleleFreq;
//        byteConversion = new double[] { -mean, 1 - mean, 2 - mean, 0 };
        byteConversion = new double[] { 0,1,2,mean };  //changed for AddPlusDom
    }

    private int genotypeIndex(int n) {
        n = n * 2;
        int i = n / 32;
        int j = n % 32;
        int intStore = bitStore[i];
        return (intStore >> j) & 3;
    }

    @Override
    public double[] getCovariate() {
        if (taxaIndex == null)
            return getCovariateNoReindex();
        else
            return getCovariateWithReindex();
    }

    @Override
    public double[] getCovariate(int[] subset) {
        if (taxaIndex == null)
            return getCovariateNoReindex(subset);
        else
            return getCovariateWithReindex(subset);
    }

    public double[] getCovariateNoReindex() {
        double[] cov = new double[ntaxa];
        int intCount = 0;
        int genoCount = 0;
        while (genoCount < ntaxa) {
            int intStore = bitStore[intCount++];
            for (int i = 0; i < 32; i += 2) {
                int genoIndex = (intStore >> i) & 3;
                cov[genoCount++] = byteConversion[genoIndex];

                if (genoCount == ntaxa)
                    break;
            }
        }
        return cov;
    }

    public double[] getCovariateNoReindex(int[] subset) {
        int nobs = subset.length;
        double[] allCov = getCovariate();
        double[] cov = new double[nobs];
        for (int i = 0; i < nobs; i++)
            cov[i] = allCov[subset[i]];
        return cov;
    }

    public double[] getCovariateWithReindex() {
        int nval = taxaIndex.length;
        double[] cov = new double[nval];
        for (int i = 0; i < nval; i++) {
            cov[i] = byteConversion[genotypeIndex(taxaIndex[i])];
        }
        return cov;
    }

    public double[] getCovariateWithReindex(int[] subset) {
        int nval = subset.length;
        double[] cov = new double[nval];
        for (int i = 0; i < nval; i++) {
            cov[i] = byteConversion[genotypeIndex(taxaIndex[subset[i]])];
        }
        return cov;
    }
    
    public static void serializeAdditiveSites(GenotypeTable geno, String outFile) {
        long start = System.nanoTime();

        try {
            ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(outFile));
            int nsites = geno.numberOfSites();
            int ntaxa = geno.numberOfTaxa();
            
            out.writeObject(new Integer(ntaxa));
            
            //for now just serialize taxon name, since Taxon is not serializeable
            TaxaList myTaxa = geno.taxa();
            for (Taxon t : myTaxa) 
                out.writeObject(t.getName());
            
            out.writeObject(new Integer(nsites));
            
            for (int s = 0; s < nsites; s++) {
                GenotypeAdditiveSite mySite = new GenotypeAdditiveSite(s, geno.chromosomeName(s), geno.chromosomalPosition(s), geno.siteName(s), AdditiveSite.CRITERION.pval, geno.genotypeAllTaxa(s), geno.majorAllele(s), geno.majorAlleleFrequency(s));
                out.writeObject(mySite);
            }
            out.close();
        } catch (FileNotFoundException e) {
            throw new RuntimeException(String.format("Error writing additive sites to %s.", outFile), e);
        } catch (IOException e) {
            throw new RuntimeException(String.format("Error writing additive sites to %s.", outFile), e);
        }

        System.out.printf("%d sites written to %s at %d ms.\n", geno.numberOfSites(), outFile, (System.nanoTime() - start) / 1000000);
    }

    @Override
    public AdditiveSite copy() {
        return new GenotypeAdditiveSite(siteIndex, chrName, position, name,
                selectionCriterion, Arrays.copyOf(bitStore, bitStore.length), ntaxa,
                Arrays.copyOf(byteConversion, byteConversion.length));
    }
}
