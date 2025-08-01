package net.maizegenetics.analysis.imputation;

import net.maizegenetics.analysis.popgen.DonorHypoth;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.util.BitSet;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import static net.maizegenetics.analysis.imputation.FILLINImputationUtils.fillInc;
import static net.maizegenetics.analysis.imputation.FILLINImputationUtils.sumOf;
import static net.maizegenetics.dna.WHICH_ALLELE.Major;
import static net.maizegenetics.dna.WHICH_ALLELE.Minor;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class FILLINImputationUtilsTest {
    private GenotypeTable origAlign;


    @Before
    public void setUp() throws Exception {
        origAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
    }

    @Test
    public void testCalcAllelePresenceCountsBtwTargetAndDonors() throws Exception {
        int indexOfB73=origAlign.taxa().indexOf("B73");
        System.out.println(origAlign.heterozygousCountForTaxon(indexOfB73));
        BitSet[] bs=new BitSet[2];
        bs[0]=origAlign.allelePresenceForAllSites(indexOfB73, Major);
        bs[1]=origAlign.allelePresenceForAllSites(indexOfB73, Minor);
        byte[][][] t=FILLINImputationUtils.calcAllelePresenceCountsBtwTargetAndDonors(bs,origAlign);
//        for (int i=0; i<t.length; i++) {
//            System.out.println(origAlign.taxaName(i)+":"+sumOf(t[i][2])+":"+Arrays.toString(t[i][2]));
//        }

        int b73Present=sumOf(t[indexOfB73][0]);
        int b73Hets=sumOf(t[indexOfB73][3]);
        int b73DiffsWithSelf=sumOf(t[indexOfB73][3]);
        //Assert.assertEquals(origAlign.totalNonMissingForTaxon(indexOfB73),b73Present);
        Assert.assertEquals(b73Hets,origAlign.heterozygousCountForTaxon(indexOfB73));  //test differences
        Assert.assertEquals(b73DiffsWithSelf,b73Hets);  //identity on hets should be differences.
    }



    @Test
    public void testGetBestInbredDonors() throws Exception {
        int indexOfB73=origAlign.taxa().indexOf("B73");
        System.out.println(origAlign.heterozygousCountForTaxon(indexOfB73));
        BitSet[] bs=new BitSet[2];
        bs[0]=origAlign.allelePresenceForAllSites(indexOfB73, Major);
        bs[1]=origAlign.allelePresenceForAllSites(indexOfB73, Minor);
        byte[][][] t=FILLINImputationUtils.calcAllelePresenceCountsBtwTargetAndDonors(bs,origAlign);
        int[] donorIndices=fillInc(0,origAlign.numberOfTaxa()-1);
        DonorHypoth[] bestDonors=FILLINImputationUtils.findHomozygousDonorHypoth(indexOfB73, 0, 1, 2, donorIndices, t, 40, 10) ;
        for (DonorHypoth bestDonor : bestDonors) {
            System.out.printf("%s %s ",origAlign.taxaName(bestDonor.targetTaxon),origAlign.taxaName(bestDonor.donor1Taxon));
            System.out.println(bestDonor);
        }
        Assert.assertEquals(bestDonors[0].donor1Taxon,indexOfB73);
//        bestDonors=FILLINImputationUtils.findHomozygousDonorHypoth(indexOfB73,0,1,1,t,40,10) ;
//        for (DonorHypoth bestDonor : bestDonors) {
//            System.out.printf("Que: %s %s ",origAlign.taxaName(bestDonor.targetTaxon),origAlign.taxaName(bestDonor.donor1Taxon));
//            System.out.println(bestDonor);
//        }
    }
}
