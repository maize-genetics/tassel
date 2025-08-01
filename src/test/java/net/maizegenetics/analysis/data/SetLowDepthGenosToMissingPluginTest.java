/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.analysis.data;

import java.util.ArrayList;
import java.util.Random;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import static net.maizegenetics.dna.snp.ExportUtils.writeToVCF;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.genotypecall.BasicGenotypeMergeRule;
import net.maizegenetics.taxa.Taxon;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Christopher Bottoms
 * @author Jeff Glaubitz
 */
public class SetLowDepthGenosToMissingPluginTest {

    boolean DEBUG = false; //Always set this to false before committing the code
    
    // Random number generator (seeded to allow for comparing results between runs)
    Random randGen = new Random(123456789);

    // alleles can be 0=A, 1=C, 2=G, or 3=T
    int maxAlleleConstant = NucleotideAlignmentConstants.T_ALLELE;

    int minDepth;
    ArrayList<int[]> genosBelowMinDepth; // int[0] = taxonIndex; int[1] = siteIndex

    //CONSTANTS
    int lowestMinDepthTested = 2;
    int highestMinDepthTested = 12;
    int nSitesPerChr = 25;
    int nChrs = 2;
    int totalNumSites = nChrs*nSitesPerChr;
    int nTaxa = 40;
    int nFixedDepthSites = 4;
    
    public SetLowDepthGenosToMissingPluginTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testRunPlugin() {
        
        // Make a sample GenotypeTable to test the plugin with
        GenotypeTable testGenos = makeTestGenos();

        // test plugin on range of minimum depths
        for (minDepth = lowestMinDepthTested; minDepth <= highestMinDepthTested; minDepth++){
            testPluginAtMinDepth(testGenos); 
        }
    }

    private void testPluginAtMinDepth(GenotypeTable testGenos) {
        //initialize (or reinitialize) list of genotypes that should be missing
        genosBelowMinDepth = new ArrayList<>(); // int[0] = taxonIndex; int[1] = siteIndex
        
        // Run plugin for given minDepth
        GenotypeTable outputGenos =
                new SetLowDepthGenosToMissingPlugin()
                        .minDepth(minDepth)
                        .runPlugin(testGenos);
        
        //Output the resulting genotype table as a VCF file (only for debugging).
        if (DEBUG) { writeToVCF(outputGenos,"testGenosMinDepth" + minDepth + ".vcf", true); }
        
        //==================================================================
        //Test all of the "FixedDepth" sites
        
        // Determine how many of the fixed depth sites should have all missing genotypes
        int nSitesWithAllMissingGenos = Math.min(nFixedDepthSites, minDepth);
        boolean someFixedSiteNotMissing = nSitesWithAllMissingGenos < nFixedDepthSites;
        
        if (DEBUG){ System.out.println("Depths after minDepth of "+minDepth+" is applied:"); }

        for (int taxonIndex=0; taxonIndex < nTaxa; taxonIndex++){
            
            if (DEBUG) {
                byte[][] byteDepths = outputGenos.depth().valuesForTaxonByte(taxonIndex);
                byte[] genos = outputGenos.genotypeAllSites(taxonIndex);
                printGenosAndDepthsForTaxon(genos, byteDepths);
            }
            
            // Check that every single site in each taxon is missing for the sites that we know should be missing
            for (int missingSiteIndex = 0; missingSiteIndex < nSitesWithAllMissingGenos; missingSiteIndex++){
                assertEquals(
                        "All genotypes for siteIndex "+missingSiteIndex+" should be missing",
                        GenotypeTable.UNKNOWN_DIPLOID_ALLELE,
                        outputGenos.genotype(taxonIndex, missingSiteIndex)
                );
            }
            
            // For fixed sites that aren't missing (if any) check that their genotype hasn't changed
            if (someFixedSiteNotMissing) {
                for (int missingSiteIndex = nSitesWithAllMissingGenos; missingSiteIndex < nFixedDepthSites; missingSiteIndex++){
                    String expectedAllele = NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) missingSiteIndex);
                    String message = "All genotypes for siteIndex "+missingSiteIndex +" should be " + expectedAllele;
                    
                    assertEquals(
                            message,
                            expectedAllele,
                            outputGenos.genotypeAsString(taxonIndex, missingSiteIndex)
                    );
                }
            }
        }
        //------------------------------------------------------------------
        
        //==================================================================
        //Test that randomly generated genotypes with depths below minDepth are missing
        for (int[] taxonSiteIndices : genosBelowMinDepth){
            int taxonIndex = taxonSiteIndices[0];
            int siteIndex = taxonSiteIndices[1];
            
            assertEquals(
                    "Genotype at taxonIndex "+taxonIndex+" and siteIndex "+siteIndex+" should be missing",
                    GenotypeTable.UNKNOWN_DIPLOID_ALLELE,
                    outputGenos.genotype(taxonSiteIndices[0], taxonSiteIndices[1])
            );
        }
        //------------------------------------------------------------------
    }
    
    private GenotypeTable makeTestGenos() {

        //Start off GenotypeTable by adding a position list first
        PositionList posits = makePositionList();
        GenotypeTableBuilder gtb = GenotypeTableBuilder.getTaxaIncremental(posits);

        //Then incrementally add each taxon along with its data for each position
        addTaxaToGenotypeTableBuilder(gtb);

        //A GenotypeTable is immutable, which is why we use a builder to create it
        GenotypeTable gt = gtb.build();

        //Output the genotype table as a VCF file (only for debugging).
        if (DEBUG) { writeToVCF(gt,"testGenos.vcf", true);}

        return gt;
    }

    private PositionList makePositionList() {
        PositionListBuilder plb = new PositionListBuilder();
        for (int chr = 1; chr <= nChrs; chr++) {
            Chromosome myChr = new Chromosome("chr"+chr);
            for (int positIndex = 1; positIndex <= nSitesPerChr; positIndex++) {
                GeneralPosition posit =             //For variety, positions increment by 10's for chr 1, by 20's for chr 2
                        new GeneralPosition.Builder(myChr,(positIndex)*chr*10).build();
                plb.add(posit);
            }
        }
        return plb.build();
    }

    private void addTaxaToGenotypeTableBuilder(GenotypeTableBuilder gtb) {

        int[][] alleles = generateAlleleJaggedArray();

        for (int taxonIndex = 0; taxonIndex < nTaxa; taxonIndex++) {
            int[][] depths = makeDepthsForTaxon(alleles, taxonIndex);
            byte[] genos = callGenosForTaxonBasedOnDepths(depths);

            //Now that genotypes have been called, let's add some depths not related to genotype
            depths = addNongenotypicDepths(depths);

            if (DEBUG) { printDepthsForTaxon(depths); }

            gtb.addTaxon(new Taxon.Builder("Taxon"+taxonIndex).build(), depths, genos);
        }
    }

    // For debugging
    private void printDepthsForTaxon ( int[][] depths ){
        for (int allele=0; allele < depths.length; allele++){
            System.out.print(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)allele) + ": ");
            for (int innerDepth : depths[allele]){
                System.out.print(innerDepth + "\t");
            }
                System.out.println("");
        }
    }

    // For debugging
    private void printGenosAndDepthsForTaxon ( byte[] genos, byte[][] depths ){

        //Print the genotypes in one line
        System.out.print("   ");
        for (byte geno: genos){
            if( geno < 0){
                System.out.print(GenotypeTable.UNKNOWN_ALLELE_STR+"\t");
            } else{
                System.out.print(NucleotideAlignmentConstants.getNucleotideIUPAC(geno)+"\t");
            }
        }
        System.out.println("");

        //Print the depths for each allele in its own line
        for (int allele=0; allele < depths.length; allele++){
            System.out.print(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)allele) + ": ");
            for (int innerDepth : depths[allele]){
                System.out.print(innerDepth + "\t");
            }
                System.out.println("");
        }
        System.out.println("");
    }

    private int[][] generateAlleleJaggedArray () {
        int[][] returnAlleles = new int[totalNumSites][];
        int siteIndex;

        // First, second, third, and fourth sites are all A's, C's, G's, and T's, respectively
        for (siteIndex=0; siteIndex < nFixedDepthSites; siteIndex++) {
            returnAlleles[siteIndex] =  new int[1];
            returnAlleles[siteIndex][0] =  siteIndex; //0 is A, 1 is C, 2 is G, 3 is T
        }

        for (siteIndex=nFixedDepthSites; siteIndex < totalNumSites; siteIndex++){

            int numAlleles = numRandomAlleles();
            returnAlleles[siteIndex] =  new int[numAlleles];

            //Randomly choose each allele (as currently implemented, it can be same allele more than once)
            for (int alleleIndex = 0; alleleIndex < numAlleles; alleleIndex++) {
                returnAlleles[siteIndex][alleleIndex] = nextRandomAllele();
            }
        } 
        return returnAlleles; 
    }

    private int nextRandInt( int minInt, int maxInt ) {
                                             // Assume an input of 2, 10
        int span    = maxInt - minInt + 1;   // then span would be 9
        int randInt = randGen.nextInt(span); // randInt would be a number from 0 to 8
        return randInt + minInt;             // returned would be a number in the inclusive range of 2 to 10 
    }

    private int nextRandomAllele () {
       return nextRandInt(0, maxAlleleConstant);
    }

    private int numRandomAlleles () {
        int rollTen = nextRandInt(1,10);

        if (rollTen == 1){      //10% chance quadrallelic
            return 4;   
        }
        else if (rollTen ==2 ){ //10% chance triallelic
            return 3;   
        }
        else if (rollTen < 6){  //30% chance biallelic
            return 2;   
        }
        else {                  //50% chance monomorphic
            return 1;  
        }
    }

    // Depths range from none for any allele to some for four alleles
    private int[][] makeDepthsForTaxon(int[][] alleles, int taxonIndex){
        int depthLimit = highestMinDepthTested - 2;
        
        // First dimension is which allele (0=A, 1=C, 2=G, 3=T (and we're currently all but ignoring 4='+' and 5='-' ))
        // Second dimension is which site
        int[][] depthsForTaxon = 
            new int[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][totalNumSites];

        //Depth for first, second, third, and fourth sites are always the same as their siteIndex (0, 1, 2, and 3, respectively)
        //Also, the allele for the first four sites are determined by their index (A for the first, C for the second, etc.)
        for (int siteIndex = 0; siteIndex < nFixedDepthSites ; siteIndex++){

            // Set all sites to zero
            for (int allele = 0; allele <= NucleotideAlignmentConstants.T_ALLELE; allele++){
                depthsForTaxon[allele][siteIndex] = 0;
            }

            for (int alleleIndex = 0; alleleIndex < alleles[siteIndex].length; alleleIndex++){
                int currentAllele = alleles[siteIndex][alleleIndex];
                depthsForTaxon[currentAllele][siteIndex] = siteIndex;
            }
        }

        //Depths for remainder of table are random. They can be mono-, bi-, tri-, or even quad- allelic.
        for (int siteIndex = nFixedDepthSites; siteIndex < totalNumSites; siteIndex++){

            int numAlleles = alleles[siteIndex].length;

            // Decide how many of the alleles possible for this site are represented in this taxon
            boolean[] presentAllele = randAllelePresent(numAlleles);

            int maxDepthPerAllele   = (int) depthLimit/numAlleles;

            // Assign depths and keep track of total depth
            int sumOfDepthForTaxonSite = 0;
            for (int alleleIndex = 0; alleleIndex < numAlleles; alleleIndex++){
                if(presentAllele[alleleIndex]){
                    int currentAllele = alleles[siteIndex][alleleIndex];
                    depthsForTaxon[currentAllele][siteIndex] = nextRandInt(1,maxDepthPerAllele);
                    sumOfDepthForTaxonSite+=depthsForTaxon[currentAllele][siteIndex];
                }
            }

            // Keep track of each taxon/site combination that has a depth below the minimum
            // (later this will be used to check that the genotypes are missing in the final GenotypeTable)
            if (sumOfDepthForTaxonSite < minDepth){
                genosBelowMinDepth.add(new int[]{taxonIndex, siteIndex});
            }
        }
        return depthsForTaxon;
    }

    // Randomly determine which alleles are present (50% chance each)
    private boolean[] randAllelePresent (int numAlleles) {
        boolean[] allelePresentAtIndex = new boolean[numAlleles]; 

        for (int i=0; i < numAlleles; i++){
            boolean coinFlip = nextRandInt(0,1) > 0;
            allelePresentAtIndex[i] = coinFlip;
        }

        return allelePresentAtIndex; 
    }

    private byte[] callGenosForTaxonBasedOnDepths(int[][] depthsForTaxon) {  // [allele][site]
        BasicGenotypeMergeRule genoMergeRule = new BasicGenotypeMergeRule(0.01);
        int nAlleles = depthsForTaxon.length;
        int[] depthsAtSite = new int[nAlleles];
        int nSites = depthsForTaxon[0].length;
        byte[] genos = new byte[nSites];
        for (int site = 0; site < nSites; site++) {
            for (int allele = 0; allele < nAlleles; allele++) {
                depthsAtSite[allele] = depthsForTaxon[allele][site];
            }
            genos[site] = genoMergeRule.callBasedOnDepth(depthsAtSite);
        }
        return genos;
    }

    // The following two tests happen because of adding nongenotype-related depth in the first and third sites.
    // x Test something already missing but has depth
    // x Test biallelic with depth in "nonalleles" (depth for "T" when it's an "AC" snp)
    private int[][] addNongenotypicDepths(int[][] depths) {

        //Add depth for T allele to first and third sites (shouldn't affect
        // outcome, since T is not part of the genotype until the fourth site)
        depths[NucleotideAlignmentConstants.T_ALLELE][0] = 10; //These would already be missing
        depths[NucleotideAlignmentConstants.T_ALLELE][2] = 10; //These should have already been called 'G'
        return depths;
    }
}