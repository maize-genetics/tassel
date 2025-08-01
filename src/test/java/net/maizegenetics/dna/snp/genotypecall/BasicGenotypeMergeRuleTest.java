package net.maizegenetics.dna.snp.genotypecall;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class BasicGenotypeMergeRuleTest {
    @Test
    public void testMergeCalls() throws Exception {
        BasicGenotypeMergeRule m=new BasicGenotypeMergeRule(0.01);
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

    /**
     * Test of resovleHetGeno method, of class BasicGenotypeMergeRule.
     */
    @Test
    public void testResolveHetGeno() {
        System.out.println("Begin testResolveHetGeno ...");
        
        // Testing these values from user issue, TAS-1076
        // GT calls are wrong for the provided depths - they
        // should all be hets (1/0, 0/1, 0/1, 0/1 respectively)
        // 
        // GT:AD:DP:GQ:PL 0/0:29,51:80:100:255,0,255  (
        // GT:AD:DP:GQ:PL 0/0:67,41:108:100:255,0,255  
        // GT:AD:DP:GQ:PL 1/1:203,80:283:100:255,0,255
        // GT:AD:DP:GQ:PL 1/1:90,82:172:99:255,0,255
        //
        
        double aveSeqErrorRate = 0.01; // default in ProductionSNPCallerPluginV2
        GenotypeMergeRule genoMergeRule = new BasicGenotypeMergeRule(aveSeqErrorRate);
        
        int[] depthsAtSite = new int[2];
        byte genos ;
        byte highByte;
        byte lowByte;
        
        // This array holds the depth of allele 1 at this site, and the depth
        // of allele 2 at this site.  Number of entries equals the number of alleles
        // Value of each entry is the number of that allele present across all the
        // taxa at this site.
        depthsAtSite[0] = 29;
        depthsAtSite[1] = 51;
        // How many times an allele shows up at a site across all the taxa tells
        // us whether this site is homoz or het.   Calls stored as HALF bytes in genoMergeRule
        
        // expecting het 1/0
        genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
        System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
                depthsAtSite[1], genos);
        highByte = (byte)((genos & 0xf0) >> 4);
        lowByte = (byte) (genos & 0x0f);
        System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
        assertEquals(highByte,(byte)1);
        assertEquals(lowByte,(byte)0);
                
        // expecting het 0/1
        depthsAtSite[0] = 67;
        depthsAtSite[1] = 41;
        genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
        System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
                depthsAtSite[1], genos);
        highByte = (byte)((genos & 0xf0) >> 4);
        lowByte = (byte) (genos & 0x0f);
        System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
        assertEquals(highByte,(byte)0);
        assertEquals(lowByte,(byte)1);

        
        // expecting het 0/1
        depthsAtSite[0] = 203;
        depthsAtSite[1] = 80;
        genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
        System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
                depthsAtSite[1], genos);
        highByte = (byte)((genos & 0xf0) >> 4);
        lowByte = (byte) (genos & 0x0f);
        System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
        assertEquals(highByte,(byte)0);
        assertEquals(lowByte,(byte)1);

        // expecting het 0/1
        depthsAtSite[0] = 90;
        depthsAtSite[1] = 82;
        genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
        System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
                depthsAtSite[1], genos);
        highByte = (byte)((genos & 0xf0) >> 4);
        lowByte = (byte) (genos & 0x0f);
        System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
        assertEquals(highByte,(byte)0);
        assertEquals(lowByte,(byte)1);
        
        
        // Expecting homozygous 1/1
        depthsAtSite[0] = 0;
        depthsAtSite[1] = 82;
        genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
        System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
                depthsAtSite[1], genos);
        highByte = (byte)((genos & 0xf0) >> 4);
        lowByte = (byte) (genos & 0x0f);
        System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
        assertEquals(highByte,(byte)1);
        assertEquals(lowByte,(byte)1);
        
        
        // Expecting homozygous 0/0
        depthsAtSite[0] = 82;
        depthsAtSite[1] = 0;
        genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
        System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
                depthsAtSite[1], genos);
        highByte = (byte)((genos & 0xf0) >> 4);
        lowByte = (byte) (genos & 0x0f);
        System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
        assertEquals(highByte,(byte)0);
        assertEquals(lowByte,(byte)0);
        
        // Expecting homozygous 0/0
        depthsAtSite[0] = 82;
        depthsAtSite[1] = 2;
        genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
        System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
                depthsAtSite[1], genos);
        highByte = (byte)((genos & 0xf0) >> 4);
        lowByte = (byte) (genos & 0x0f);
        System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
        assertEquals(highByte,(byte)0);
        assertEquals(lowByte,(byte)0);
        
        // Expecting homozygous 1/1
        depthsAtSite[0] = 3;
        depthsAtSite[1] = 82;
        genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
        System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
                depthsAtSite[1], genos);
        highByte = (byte)((genos & 0xf0) >> 4);
        lowByte = (byte) (genos & 0x0f);
        System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
        assertEquals(highByte,(byte)1);
        assertEquals(lowByte,(byte)1);
        
        // Try with 3 alleles: Expecting homozygous 1/1
        int[] threeAlleles = new int[3];
        threeAlleles[0] = 3;
        threeAlleles[1] = 82;
        threeAlleles[2] = 1;
        genos = genoMergeRule.callBasedOnDepth(threeAlleles);
        System.out.printf("Site with depths %d and %d and %d is called as %x\n", threeAlleles[0],
                threeAlleles[1], threeAlleles[2],genos);
        highByte = (byte)((genos & 0xf0) >> 4);
        lowByte = (byte) (genos & 0x0f);
        System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
        assertEquals(highByte,(byte)1);
        assertEquals(lowByte,(byte)1);
        
        // THe next series return homozygotes because the likelihood ratio isn't quite met       
        depthsAtSite[0] = 100;
        depthsAtSite[1] = 338;
        genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
        System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
                depthsAtSite[1], genos);
        highByte = (byte)((genos & 0xf0) >> 4);
        lowByte = (byte) (genos & 0x0f);
        System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
        assertEquals(highByte,(byte)1);
        assertEquals(lowByte,(byte)0);
 
        
        // not enough minor allele - requires 31
        // THis now is correctly called as het since changing to max 100 for binomial dist.
        depthsAtSite[0] = 134;
        depthsAtSite[1] = 29;
        genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
        System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
                depthsAtSite[1], genos);
        highByte = (byte)((genos & 0xf0) >> 4);
        lowByte = (byte) (genos & 0x0f);
        System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
        assertEquals(highByte,(byte)0);
        assertEquals(lowByte,(byte)1);
        
        // bumped up, minor allele count now results in het
        depthsAtSite[0] = 134;
        depthsAtSite[1] = 32;
        genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
        System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
                depthsAtSite[1], genos);
        highByte = (byte)((genos & 0xf0) >> 4);
        lowByte = (byte) (genos & 0x0f);
        System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
        assertEquals(highByte,(byte)0);
        assertEquals(lowByte,(byte)1);
        
      depthsAtSite[0] = 85;
      depthsAtSite[1] = 12; // needs 13
      genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
      System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
              depthsAtSite[1], genos);
      highByte = (byte)((genos & 0xf0) >> 4);
      lowByte = (byte) (genos & 0x0f);
      System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
      assertEquals(highByte,(byte)0);
      assertEquals(lowByte,(byte)0);
           
      // Testing border case: Likelihood ratio wants 18 here, not 17
      // expecting homo 0/0
      depthsAtSite[0] = 98;
      depthsAtSite[1] = 18;
      genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
      System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
              depthsAtSite[1], genos);
      highByte = (byte)((genos & 0xf0) >> 4);
      lowByte = (byte) (genos & 0x0f);
      System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
      assertEquals(highByte,(byte)0);
      assertEquals(lowByte,(byte)1);     
      
      // Test total counts > 500
      
      // expecting het 0/1 
      depthsAtSite[0] = 323;
      depthsAtSite[1] = 323;
      genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
      System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
              depthsAtSite[1], genos);
      highByte = (byte)((genos & 0xf0) >> 4);
      lowByte = (byte) (genos & 0x0f);
      System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
      assertEquals(highByte,(byte)0);
      assertEquals(lowByte,(byte)1);
      
      depthsAtSite[0] = 458;
      depthsAtSite[1] = 256;
      genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
      System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
              depthsAtSite[1], genos);
      highByte = (byte)((genos & 0xf0) >> 4);
      lowByte = (byte) (genos & 0x0f);
      System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
      assertEquals(highByte,(byte)0);
      assertEquals(lowByte,(byte)1);
      
      // expecting homo 0/0
      depthsAtSite[0] = 672;
      depthsAtSite[1] = 40;
      genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
      System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
              depthsAtSite[1], genos);
      highByte = (byte)((genos & 0xf0) >> 4);
      lowByte = (byte) (genos & 0x0f);
      System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
      assertEquals(highByte,(byte)0);
      assertEquals(lowByte,(byte)0);
      
      // expecting homo 1/1
      depthsAtSite[0] = 56;
      depthsAtSite[1] = 872;
      genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
      System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
              depthsAtSite[1], genos);
      highByte = (byte)((genos & 0xf0) >> 4);
      lowByte = (byte) (genos & 0x0f);
      System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
      assertEquals(highByte,(byte)1);
      assertEquals(lowByte,(byte)1);
      
      // These last 2 test cases show an inconsistency in how we look at the data.
      // We call hets based on a Binomial Distribution formula when our totals counts
      // are 0-499.  Once we hit 500, the software calls hets based on the formula
      //    minorAlleleCount/totalCount < 0.1 = homozygous
      // THis means we need 157 out 0f 497 to call a het, but only 80 out of 400 to call it het.
      
      //NOTE: this was resolved by only calculating Binomial Dist for totals up to 200.
      // It will be revisited when BasicGenotypeMergeRule is changed to match coming
      // changes to ExportUtils:writeToVCF().
      
      // Compare highest calculated range with .01 specified.
      // We calculate up to 499, the rest we merely ensure minorAlleleNumber/totalCount < .1 to be homo
      // NOTE: - change to 199 from 499, then .1 as base for het
      depthsAtSite[0] = 157;
      depthsAtSite[1] = 341;
      genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
      System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
              depthsAtSite[1], genos);
      highByte = (byte)((genos & 0xf0) >> 4);
      lowByte = (byte) (genos & 0x0f);
      System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
      assertEquals(highByte,(byte)1);
      assertEquals(lowByte,(byte)0);
      
      depthsAtSite[0] = 80;
      depthsAtSite[1] = 420;
      genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
      System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
              depthsAtSite[1], genos);
      highByte = (byte)((genos & 0xf0) >> 4);
      lowByte = (byte) (genos & 0x0f);
      System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
      assertEquals(highByte,(byte)1);
      assertEquals(lowByte,(byte)0);
      
      System.out.println("\nLCJ - look at this one ...");
      depthsAtSite[0] = 310;
      depthsAtSite[1] = 122;
      genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
      System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
              depthsAtSite[1], genos);
      highByte = (byte)((genos & 0xf0) >> 4);
      lowByte = (byte) (genos & 0x0f);
      System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
      assertEquals(highByte,(byte)0);
      assertEquals(lowByte,(byte)1);
        
    }
    
    @Test
    public void create500MatrixValues() {
        System.out.println("Begin create500MatrixValues ...");
        // THis is run when we have the debug line in setLIkelihoodThresh of BasicGenotypeMergeRule
        // That method is called when the constructor is executed, which prints out
        // the matrix of infor for us (because I put in code to do that).  Ed already
        // had code - I uncommented it and added additoinal info  To run this test and
        // get the table results, uncomment the print statement in setLikelihoodThresh()
        
        // Testing first 500 values
        
        double aveSeqErrorRate = 0.01; // default in ProductionSNPCallerPluginV2
        GenotypeMergeRule genoMergeRule = new BasicGenotypeMergeRule(aveSeqErrorRate);
        
        int[] depthsAtSite = new int[2];
        byte genos ;
        byte highByte;
        byte lowByte;
        
        // This array holds the depth of allele 1 at this site, and the depth
        // of allele 2 at this site.  Number of entries equals the number of alleles
        // Value of each entry is the number of that allele present across all the
        // taxa at this site.
        depthsAtSite[0] = 29;
        depthsAtSite[1] = 51;
        // How many times an allele shows up at a site across all the taxa tells
        // us whether this site is homoz or het.   Calls stored as HALF bytes in genoMergeRule
        
        // expecting het 1/0
        genos = genoMergeRule.callBasedOnDepth(depthsAtSite);
        System.out.printf("Site with depths %d and %d is called as %x\n", depthsAtSite[0],
                depthsAtSite[1], genos);
        highByte = (byte)((genos & 0xf0) >> 4);
        lowByte = (byte) (genos & 0x0f);
        System.out.println("   high byte: " + highByte + " lowByte " + lowByte);
    }
}
