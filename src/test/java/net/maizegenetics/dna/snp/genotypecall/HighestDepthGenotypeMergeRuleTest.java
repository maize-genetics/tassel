package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import org.junit.Assert;
import org.junit.Test;

public class HighestDepthGenotypeMergeRuleTest {
    @Test
    public void testMergeCalls() throws Exception {
        HighestDepthGenotypeMergeRule m=new HighestDepthGenotypeMergeRule();
        byte[] a=new byte[]{0,2,0,0,0,0};
        byte[] b=new byte[]{0,2,0,0,0,0};
        byte[] apb=m.mergeWithDepth(a,b);
        System.out.printf("%d %d %d %n",m.callBasedOnDepth(a),m.callBasedOnDepth(b),m.callBasedOnDepth(apb));
        a=new byte[]{4,0,0,4,0,0};
        b=new byte[]{4,0,0,4,0,0};
        apb=m.mergeWithDepth(a,b);
        System.out.printf("%d %d %d %n",m.callBasedOnDepth(a),m.callBasedOnDepth(b),m.callBasedOnDepth(apb));
        a=new byte[]{0, 4, 0, 4, 0, 0};
        b=new byte[]{0, 8, 0, 8, 0, 0};
        apb=m.mergeWithDepth(a,b);
        System.out.printf("%d %d %d %n",m.callBasedOnDepth(a),m.callBasedOnDepth(b),m.callBasedOnDepth(apb));
    }

    //This one should set to N
    @Test
    public void testNonDepthCalls() {
        HighestDepthGenotypeMergeRule mergeRule = new HighestDepthGenotypeMergeRule();

        byte homGenotype1 = (NucleotideAlignmentConstants.G_ALLELE << 4) | NucleotideAlignmentConstants.G_ALLELE;
        byte homGenotype2 = (NucleotideAlignmentConstants.G_ALLELE << 4) | NucleotideAlignmentConstants.G_ALLELE;
        byte homGenotype3 = (NucleotideAlignmentConstants.T_ALLELE << 4) | NucleotideAlignmentConstants.T_ALLELE;

        byte merged1_2 = mergeRule.mergeCalls(homGenotype1,homGenotype2);
        byte merged1_3 = mergeRule.mergeCalls(homGenotype1,homGenotype3);

        //merged1_2 should be GG call
        Assert.assertEquals("HighestDepthGenotypeMergeRuleTest testNonDepthCalls: Incorrect merged GG->GG",homGenotype1,merged1_2);
        //merged1_3 should be NN
        Assert.assertEquals("HighestDepthGenotypeMergeRuleTest testNonDepthCalls: Incorrect merged GT->NN", GenotypeTable.UNKNOWN_DIPLOID_ALLELE,merged1_3);
    }

    //Should resolve GT to GG
    @Test
    public void testDepthCalls() {
        HighestDepthGenotypeMergeRule mergeRule = new HighestDepthGenotypeMergeRule();
        byte[] homGenotypeDepths1 = new byte[]{0,0,10,0,0,0};
        byte homGenotype1 = (NucleotideAlignmentConstants.G_ALLELE << 4) | NucleotideAlignmentConstants.G_ALLELE;
        byte homGenotype2 = (NucleotideAlignmentConstants.T_ALLELE << 4) | NucleotideAlignmentConstants.T_ALLELE;

        byte[] homGenotypeDepths2 = new byte[]{0,0,4,3,0,0};
        byte[] homGenotypeDepths3 = new byte[]{0,1,4,5,0,0};

        byte merge1 = mergeRule.callBasedOnDepth(homGenotypeDepths1);
        byte merge2 = mergeRule.callBasedOnDepth(homGenotypeDepths2);
        byte merge3 = mergeRule.callBasedOnDepth(homGenotypeDepths3);

        Assert.assertEquals("HighestDepthGenotypeMergeRuleTest testNonDepthCalls: Incorrect merged GG->GG",homGenotype1,merge1);
        //merged1_2 should be GG
        Assert.assertEquals("HighestDepthGenotypeMergeRuleTest testNonDepthCalls: Incorrect merged GT->GG", homGenotype1,merge2);
        Assert.assertEquals("HighestDepthGenotypeMergeRuleTest testNonDepthCalls: Incorrect merged GT->TT", homGenotype2,merge3);

    }
}
