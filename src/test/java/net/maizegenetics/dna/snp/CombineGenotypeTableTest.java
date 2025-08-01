/*
 *  CombineGenotypeTableTest
 * 
 *  Created on Feb 23, 2017
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.analysis.filter.FilterTaxaBuilderPlugin;
import net.maizegenetics.constants.TutorialConstants;
import org.junit.*;

/**
 *
 * @author Terry Casstevens
 */
public class CombineGenotypeTableTest {

    private GenotypeTable myExpectedAlignment;
    private GenotypeTable[] mySeparateChrs;

    public CombineGenotypeTableTest() {
    }

    @Before
    public void setUp() {

        myExpectedAlignment = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME);

        mySeparateChrs = new GenotypeTable[10];
        for (int c = 1; c < 11; c++) {
            mySeparateChrs[c - 1] = ImportUtils.readFromHapmap(TutorialConstants.TUTORIAL_DIR + "mdp_genotype_chrom" + c + ".hmp.txt");
        }

    }

    //@Test
    // CombineGenotypeTable isn't fully implemented and may not be needed.
    public void testCombiningGenotypes() {

        System.out.println("Testing Combine Genotypes...");

        GenotypeTable combineAlign = CombineGenotypeTable.getInstance(mySeparateChrs, false);

        AlignmentTestingUtils.alignmentsEqual(myExpectedAlignment, combineAlign);

    }

    //@Test
    // CombineGenotypeTable isn't fully implemented and may not be needed.
    public void testCombiningGenotypesFewerTaxa() {

        System.out.println("Testing Combine Genotypes with fewer taxa in Chromosome 1");

        FilterTaxaBuilderPlugin filter = new FilterTaxaBuilderPlugin()
                .includeTaxa(false)
                .taxaList("A188,A272,B79");
        GenotypeTable filterGenotype = filter.runPlugin(myExpectedAlignment);

        GenotypeTable[] individualChrs = new GenotypeTable[10];
        System.arraycopy(mySeparateChrs, 0, individualChrs, 0, 10);
        individualChrs[0] = ImportUtils.readFromHapmap(TutorialConstants.TUTORIAL_DIR + "mdp_genotype_chrom1_A188_A272_B79_Missing.hmp.txt");

        GenotypeTable combineAlign = CombineGenotypeTable.getInstance(individualChrs, false);

        AlignmentTestingUtils.alignmentsEqual(filterGenotype, combineAlign);

    }

    //@Test
    // CombineGenotypeTable isn't fully implemented and may not be needed.
    public void testCombiningGenotypesTaxaDiffOrder() {

        System.out.println("Testing Combine Genotypes with taxa different order in Chromosome 3");

        GenotypeTable[] individualChrs = new GenotypeTable[10];
        System.arraycopy(mySeparateChrs, 0, individualChrs, 0, 10);
        individualChrs[2] = ImportUtils.readFromHapmap(TutorialConstants.TUTORIAL_DIR + "mdp_genotype_chrom3_taxa_diff_order.hmp.txt");

        GenotypeTable combineAlign = CombineGenotypeTable.getInstance(individualChrs, false);

        AlignmentTestingUtils.alignmentsEqual(myExpectedAlignment, combineAlign);

    }

    //@Test
    // CombineGenotypeTable isn't fully implemented and may not be needed.
    public void testCombiningGenotypesUnion() {

        System.out.println("Testing Combine Genotypes Union...");

        GenotypeTable combineAlign = CombineGenotypeTable.getInstance(mySeparateChrs, true);

        AlignmentTestingUtils.alignmentsEqual(myExpectedAlignment, combineAlign);

    }

    //@Test
    // CombineGenotypeTable isn't fully implemented and may not be needed.
    public void testCombiningGenotypesTaxaDiffOrderUnion() {

        System.out.println("Testing Combine Genotypes Union with taxa different order in Chromosome 3");

        GenotypeTable[] individualChrs = new GenotypeTable[10];
        System.arraycopy(mySeparateChrs, 0, individualChrs, 0, 10);
        individualChrs[2] = ImportUtils.readFromHapmap(TutorialConstants.TUTORIAL_DIR + "mdp_genotype_chrom3_taxa_diff_order.hmp.txt");

        GenotypeTable combineAlign = CombineGenotypeTable.getInstance(individualChrs, true);

        AlignmentTestingUtils.alignmentsEqual(myExpectedAlignment, combineAlign);

    }

    //@Test
    // CombineGenotypeTable isn't fully implemented and may not be needed.
    public void testCombiningGenotypesFewerTaxaUnion() {

        System.out.println("Testing Combine Genotypes with fewer taxa in Chromosome 1");

        FilterTaxaBuilderPlugin filter = new FilterTaxaBuilderPlugin()
                .includeTaxa(false)
                .taxaList("WD,WF9,YU796NS");
        GenotypeTable filterGenotype = filter.runPlugin(mySeparateChrs[0]);

        GenotypeTable[] individualChrs = new GenotypeTable[2];
        individualChrs[0] = filterGenotype;
        individualChrs[1] = mySeparateChrs[1];

        GenotypeTable combineAlign = CombineGenotypeTable.getInstance(individualChrs, true);

        GenotypeTable expectedAlignment = ImportUtils.readFromHapmap(TutorialConstants.TUTORIAL_DIR + "mdp_genotype_chrom1_WD_WF9_YU79NS_Missing_chrom2.hmp.txt");

        AlignmentTestingUtils.alignmentsEqual(expectedAlignment, combineAlign);

    }

}
