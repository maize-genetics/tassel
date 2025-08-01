/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author edbuckler
 */
public class GenotypeTableUtilsTest {
    
    public GenotypeTableUtilsTest() {
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

 

    /**
     * Test of isHeterozygous method, of class GenotypeTableUtils.
     */
    @Test
    public void testIsHeterozygous() {
        System.out.println("isHeterozygous");        
        boolean result = GenotypeTableUtils.isHeterozygous(NucleotideAlignmentConstants.getNucleotideDiploidByte("R"));
        assertEquals(true, result);
        result = GenotypeTableUtils.isHeterozygous(NucleotideAlignmentConstants.getNucleotideDiploidByte("A"));
        assertEquals(false, result);
        result = GenotypeTableUtils.isHeterozygous(NucleotideAlignmentConstants.getNucleotideDiploidByte("M"));
        assertEquals(true, result);
        result = GenotypeTableUtils.isHeterozygous(NucleotideAlignmentConstants.getNucleotideDiploidByte("-"));
        assertEquals(false, result);
        result = GenotypeTableUtils.isHeterozygous(NucleotideAlignmentConstants.getNucleotideDiploidByte("0"));
        assertEquals(true, result);
    }


    /**
     * Test of getDiploidValue method, of class GenotypeTableUtils.
     */
    @Test
    public void testGetDiploidValuePhased() {
        System.out.println("getDiploidValuePhased");
        byte a = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        byte b = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        byte expResult = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        byte result = GenotypeTableUtils.getDiploidValuePhased(a, b);
        assertEquals(expResult, result);
        a = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("T");
        expResult = NucleotideAlignmentConstants.getNucleotideDiploidByte("AT");
        result = GenotypeTableUtils.getDiploidValuePhased(a, b);
        System.out.printf("A+T=AT:%d %d %d %d %n",a,b,expResult,result);
        assertEquals(expResult, result);
        a = NucleotideAlignmentConstants.getNucleotideDiploidByte("T");
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        expResult = NucleotideAlignmentConstants.getNucleotideDiploidByte("TA");
        result = GenotypeTableUtils.getDiploidValuePhased(a, b);
        a = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("T");
        expResult = NucleotideAlignmentConstants.getNucleotideDiploidByte("W");
        result = GenotypeTableUtils.getDiploidValuePhased(a, b);
        System.out.printf("T+A=TA: %d %d %d %d %n",a,b,expResult,result);
        assertEquals(expResult, result);
    }
    /**
     * Test of getUnphasedDiploidValueUnPhased method, of class GenotypeTableUtils.
     */
    @Test
    public void testGetUnphasedDiploidValue() {
        System.out.println("getUnphasedDiploidValue");
        byte a = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        byte b = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        byte expResult = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        byte result = GenotypeTableUtils.getUnphasedDiploidValue(a, b);
        assertEquals(expResult, result);
        a = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("T");
        expResult = NucleotideAlignmentConstants.getNucleotideDiploidByte("AT");
        result = GenotypeTableUtils.getUnphasedDiploidValue(a, b);
        System.out.printf("%d %d %d %d %n",a,b,expResult,result);
        assertEquals(expResult, result);
        a = NucleotideAlignmentConstants.getNucleotideDiploidByte("T");
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        expResult = NucleotideAlignmentConstants.getNucleotideDiploidByte("AT");
        result = GenotypeTableUtils.getUnphasedDiploidValue(a, b);
        System.out.printf("%d %d %d %d %n",a,b,expResult,result);
        assertEquals(expResult, result);
        a = NucleotideAlignmentConstants.getNucleotideDiploidByte("T");
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("C");
        expResult = NucleotideAlignmentConstants.getNucleotideDiploidByte("Y");
        result = GenotypeTableUtils.getUnphasedDiploidValue(a, b);
        System.out.printf("%d %d %d %d %n",a,b,expResult,result);
        assertEquals(expResult, result);
        
    }

    /**
     * Test of getDiploidValueForPotentialHets method, of class GenotypeTableUtils.
     */
    @Test
    public void testGetUnphasedDiploidValueNoHets() {
        System.out.println("getUnphasedDiploidValueNoHets");
        byte a = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        byte b = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        byte expResult = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        byte result = GenotypeTableUtils.getUnphasedDiploidValueNoHets(a, b);
        System.out.printf("%d %d %d %d %n",a,b,expResult,result);
        assertEquals(expResult, result);
        a = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("T");
        expResult = NucleotideAlignmentConstants.getNucleotideDiploidByte("W");
        result = GenotypeTableUtils.getUnphasedDiploidValueNoHets(a, b);
        System.out.printf("%d %d %d %d %n",a,b,expResult,result);
        assertEquals(expResult, result);
        a = NucleotideAlignmentConstants.getNucleotideDiploidByte("T");
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        expResult = NucleotideAlignmentConstants.getNucleotideDiploidByte("W");
        result = GenotypeTableUtils.getUnphasedDiploidValueNoHets(a, b);
        System.out.printf("%d %d %d %d %n",a,b,expResult,result);
        assertEquals(expResult, result);
        a = NucleotideAlignmentConstants.getNucleotideDiploidByte("W");
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        expResult = NucleotideAlignmentConstants.getNucleotideDiploidByte("N");
        result = GenotypeTableUtils.getUnphasedDiploidValueNoHets(a, b);
        System.out.printf("%d %d %d %d %n",a,b,expResult,result);
        assertEquals(expResult, result);
        a = NucleotideAlignmentConstants.getNucleotideDiploidByte("N");
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("A");
        expResult = NucleotideAlignmentConstants.getNucleotideDiploidByte("N");
        result = GenotypeTableUtils.getUnphasedDiploidValueNoHets(a, b);
        System.out.printf("%d %d %d %d %n",a,b,expResult,result);
        assertEquals(expResult, result);
    }
    
    @Test
    public void testIsPartiallyEqual() {
        System.out.println("isPartiallyEqual");
        byte a = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
        byte b = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
        assertFalse(GenotypeTableUtils.isPartiallyEqual(a, b));
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
        assertTrue(GenotypeTableUtils.isPartiallyEqual(a, b));
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("CA");
        assertTrue(GenotypeTableUtils.isPartiallyEqual(a, b));
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
        assertTrue(GenotypeTableUtils.isPartiallyEqual(a, b));
        b = NucleotideAlignmentConstants.getNucleotideDiploidByte("NN");
        assertFalse(GenotypeTableUtils.isPartiallyEqual(a, b));
    }
    
    @Test
    public void testCalcBitPresenceFromGenotype() {
        System.out.println("test calcBitPresenceFromGenotype");
        GenotypeTable inputAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
        byte[] majSeq=new byte[inputAlign.numberOfSites()];
        byte[] minSeq=new byte[inputAlign.numberOfSites()];
        for (int i = 0; i < minSeq.length; i++) {
            majSeq[i]=inputAlign.majorAllele(i);
            minSeq[i]=inputAlign.minorAllele(i);  
        }
        int taxon=6;  //6 is good
        byte[] genotype=inputAlign.genotypeAllSites(taxon);
        System.out.print("MAJ   :");
        for (int i = 0; i < 64; i++) {
            System.out.print(NucleotideAlignmentConstants.getNucleotideIUPAC(GenotypeTableUtils.getDiploidValuePhased(majSeq[i],majSeq[i])));
        }
        System.out.println();
        System.out.print("MIN   :");
        for (int i = 0; i < 64; i++) {
            System.out.print(NucleotideAlignmentConstants.getNucleotideIUPAC(GenotypeTableUtils.getDiploidValuePhased(minSeq[i],minSeq[i])));
        }
        System.out.println();
        System.out.print("Seq   :");
        for (int i = 0; i < 64; i++) {
            System.out.print(NucleotideAlignmentConstants.getNucleotideIUPAC(genotype[i]));
        }
        System.out.println();
        BitSet mjBS=inputAlign.allelePresenceForAllSites(taxon, WHICH_ALLELE.Major);
        BitSet mnBS=inputAlign.allelePresenceForAllSites(taxon, WHICH_ALLELE.Minor);
        System.out.println("mjBS  :"+ BitUtil.toPadStringLowSiteToHighSite(mjBS.getBits()[0]));
        System.out.println("mnBS  :"+ BitUtil.toPadStringLowSiteToHighSite(mnBS.getBits()[0]));
 //       BitSet[] gBS=GenotypeTableUtils.calcBitPresenceFromGenotype(genotype,majSeq,minSeq);
        BitSet[] gBS=GenotypeTableUtils.calcBitPresenceFromGenotype(genotype,majSeq,minSeq);
        assertArrayEquals(mjBS.getBits(), gBS[0].getBits());
        assertArrayEquals(mnBS.getBits(), gBS[1].getBits());
        gBS=GenotypeTableUtils.calcBitPresenceFromGenotype15(genotype,majSeq,minSeq);
        System.out.println("gBS[0]:"+ BitUtil.toPadStringLowSiteToHighSite(gBS[0].getBits()[0]));
        System.out.println("gBS[1]:"+ BitUtil.toPadStringLowSiteToHighSite(gBS[1].getBits()[0]));
        System.out.println(OpenBitSet.andNotCount(mjBS, gBS[0]));
        System.out.println(OpenBitSet.andNotCount(mnBS, gBS[1]));
        assertArrayEquals(mjBS.getBits(), gBS[0].getBits());
        assertArrayEquals(mnBS.getBits(), gBS[1].getBits());
//        gBS=GenotypeTableUtils.calcBitPresenceFromGenotype3(genotype,majSeq,minSeq);
//        assertArrayEquals(mjBS.getBits(), gBS[0].getBits());
//        assertArrayEquals(mnBS.getBits(), gBS[1].getBits());
 
//        for (int i = 40; i < 50; i++) {
//            System.out.println("A"+i+":"+Arrays.deepToString(inputAlign.allelesSortedByFrequency(i)));
//        }

//        
//        System.out.print("MjGet :");
//        for (int i = 0; i < 64; i++) {
//            if(mjBS.get(i)) {System.out.print(1);} else System.out.print(0);
//        }
//        System.out.println();
//        System.out.print("MjGet :");
//        for (int i = 0; i < 64; i++) {
//            if(mnBS.get(i)) {System.out.print(1);} else System.out.print(0);
//        }
//        System.out.println();
//
        System.out.println("mjBS  :"+ BitUtil.toPadStringLowSiteToHighSite(mjBS.getBits()[0]));
        System.out.println("gBS[0]:"+ BitUtil.toPadStringLowSiteToHighSite(gBS[0].getBits()[0]));
        System.out.println("mnBS  :"+ BitUtil.toPadStringLowSiteToHighSite(mnBS.getBits()[0]));
        System.out.println("gBS[1]:"+ BitUtil.toPadStringLowSiteToHighSite(gBS[1].getBits()[0]));
  
        long time=System.currentTimeMillis();
        long card=0;
        for (int i = 0; i < 5000; i++) {
            byte temp=genotype[i%3000];
            genotype[i%3000]=genotype[i%1000];
            genotype[i%1000]=temp;
            gBS=GenotypeTableUtils.calcBitPresenceFromGenotype(genotype,majSeq,minSeq);
            card+=gBS[1].cardinality();
            temp=genotype[i%1000];
            genotype[i%1000]=genotype[i%3000];
            genotype[i%3000]=temp;
        }
        System.out.printf("Reg Time: %d Card:%d %n",System.currentTimeMillis()-time, card);
        time=System.currentTimeMillis();
        for (int i = 0; i < 5000; i++) {
            byte temp=genotype[i%3000];
            genotype[i%3000]=genotype[i%1000];
            genotype[i%1000]=temp;
             gBS=GenotypeTableUtils.calcBitPresenceFromGenotype15(genotype,majSeq,minSeq);
            card+=gBS[1].cardinality();
            temp=genotype[i%1000];
            genotype[i%1000]=genotype[i%3000];
            genotype[i%3000]=temp;
        }
        System.out.printf("IndividBitSet Time: %d Card:%d %n",System.currentTimeMillis()-time, card);
        int scalar=1000000/3000;
        byte[] taxaL=new byte[inputAlign.numberOfSites()*scalar];
        byte[] majSeqL=new byte[inputAlign.numberOfSites()*scalar];
        byte[] minSeqL=new byte[inputAlign.numberOfSites()*scalar];
        for (int i = 0; i < minSeq.length*scalar; i++) {
            taxaL[i]=inputAlign.majorAllele(i%minSeq.length);
            majSeqL[i%minSeq.length]=inputAlign.majorAllele(i%minSeq.length);
            minSeqL[i%minSeq.length]=inputAlign.minorAllele(i%minSeq.length);  
        }
        time=System.currentTimeMillis();
        for (int i = 0; i < 40; i++) {
            byte temp=genotype[i%3000];
            genotype[i%3000]=genotype[i%1000];
            genotype[i%1000]=temp;
            gBS=GenotypeTableUtils.calcBitPresenceFromGenotype(taxaL,majSeqL,minSeqL);

            card+=gBS[1].cardinality();
            temp=genotype[i%1000];
            genotype[i%1000]=genotype[i%3000];
            genotype[i%3000]=temp;
        }
        System.out.printf("TestLong15 Time: %d Card:%d %n",System.currentTimeMillis()-time, card);
    }

 

  
}
