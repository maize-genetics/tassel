package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import org.junit.Assert;
import org.junit.Test;

public class SetToNGenotypeMergeRuleTest {
    @Test
    public void testMergeCalls() throws Exception {
        SetToNGenotypeMergeRule m=new SetToNGenotypeMergeRule();
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

    @Test
    public void testNonDepthCalls() {
        SetToNGenotypeMergeRule mergeRule = new SetToNGenotypeMergeRule();

        byte homGenotype1 = (NucleotideAlignmentConstants.G_ALLELE << 4) | NucleotideAlignmentConstants.G_ALLELE;
        byte homGenotype2 = (NucleotideAlignmentConstants.G_ALLELE << 4) | NucleotideAlignmentConstants.G_ALLELE;
        byte homGenotype3 = (NucleotideAlignmentConstants.T_ALLELE << 4) | NucleotideAlignmentConstants.T_ALLELE;

        byte merged1_2 = mergeRule.mergeCalls(homGenotype1,homGenotype2);
        byte merged1_3 = mergeRule.mergeCalls(homGenotype1,homGenotype3);

        //merged1_2 should be GG call
        Assert.assertEquals("SetToNGenotypeMergeRuleTest testNonDepthCalls: Incorrect merged GG->GG",homGenotype1,merged1_2);
        //merged1_3 should be NN
        Assert.assertEquals("SetToNGenotypeMergeRuleTest testNonDepthCalls: Incorrect merged GT->NN", GenotypeTable.UNKNOWN_DIPLOID_ALLELE,merged1_3);
    }

    @Test
    public void testDepthCalls() {
        SetToNGenotypeMergeRule mergeRule = new SetToNGenotypeMergeRule();
        byte[] homGenotypeDepths1 = new byte[]{0,0,10,0,0,0};
        byte homGenotype1 = (NucleotideAlignmentConstants.G_ALLELE << 4) | NucleotideAlignmentConstants.G_ALLELE;

        byte[] homGenotypeDepths2 = new byte[]{0,0,4,3,0,0};

        byte merge1 = mergeRule.callBasedOnDepth(homGenotypeDepths1);
        byte merge2 = mergeRule.callBasedOnDepth(homGenotypeDepths2);

        Assert.assertEquals("SetToNGenotypeMergeRuleTest testNonDepthCalls: Incorrect merged GG->GG",homGenotype1,merge1);
        //merged1_3 should be NN
        Assert.assertEquals("SetToNGenotypeMergeRuleTest testNonDepthCalls: Incorrect merged GT->NN", GenotypeTable.UNKNOWN_DIPLOID_ALLELE,merge2);

    }
}
