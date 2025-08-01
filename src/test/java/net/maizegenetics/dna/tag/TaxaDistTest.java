package net.maizegenetics.dna.tag;

import com.google.common.collect.Multiset;
import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.Random;

public class TaxaDistTest {
    private TaxaDistribution basicLrgTD;
    private int[] tagsScored ={0, 5, 10000, 32766, 32767, 32768, 32769, 70000, 70000};
    private int maxTaxa=75000;


    @Before
    public void setUp() throws Exception {
        basicLrgTD=TaxaDistBuilder.create(maxTaxa);
        for (int i : tagsScored) {
            basicLrgTD.increment(i);
        }

    }

    @Test
    public void testIncrement() throws Exception {
        TaxaDistribution copyTD=TaxaDistBuilder.create(basicLrgTD);
        System.out.println();
        System.out.println(copyTD.toString());

    }

    @Test
    public void testDepths() throws Exception {
        int[] depths=basicLrgTD.depths();
        for (int i : tagsScored) {
            Assert.assertTrue(depths[i]>0);
        }
        Assert.assertTrue(depths[70000]==2);
    }

    @Test
    public void testTaxaDepthArray() throws Exception {
        System.out.println(basicLrgTD.toString());
        Multiset<Integer> map=basicLrgTD.taxaDepthMap();
        Assert.assertEquals(9,map.size());
        Assert.assertEquals(8,map.elementSet().size());
    }

    @Test
    public void testRatesOfAccess() throws Exception {
        long addTime=0, mapTime=0, mapPrimTime=0, arrayTime=0;
        int sumDepth=0, sumEntry=0, sumEntry2=0;
        int taxaNum=80000;
        int trials=1000;
        Random r=new Random();
        int[] reads=new int[taxaNum/10];
        for (int i = 0; i < reads.length; i++) {
            reads[i]=r.nextInt(taxaNum);
        }
        for (int i = 0; i < trials; i++) {
            long time=System.nanoTime();
            TaxaDistribution td=TaxaDistBuilder.create(taxaNum);
            for (int read : reads) {
                td.increment(read);
            }
            addTime+=System.nanoTime()-time;
            time=System.nanoTime();
            int[] depthArray=td.depths();
            sumDepth+=depthArray[0];
            arrayTime+=System.nanoTime()-time;
            time=System.nanoTime();
            Multiset<Integer> map=td.taxaDepthMap();
            sumEntry+=map.size();
            mapTime+=System.nanoTime()-time;
            time=System.nanoTime();
            int[][] mapPrim=td.taxaWithDepths();
            sumEntry2+=mapPrim.length;
            mapPrimTime+=System.nanoTime()-time;
        }
        System.out.println("addTime = " + addTime/trials);
        System.out.println("mapTime = " + mapTime/trials + " sumEntry="+sumEntry);
        System.out.println("mapPrimTime = " + mapPrimTime/trials + " sumEntry="+sumEntry2);
        System.out.println("arrayTime = " + arrayTime/trials + " sumDepth="+sumDepth);
    }

    @Test
    public void testTotalDepth() throws Exception {
        Assert.assertEquals("Total Depth Error", tagsScored.length,basicLrgTD.totalDepth());
    }

    @Test
    public void testEncoding() throws Exception {
        byte[] encodeDepths=basicLrgTD.encodeTaxaDepth();
        TaxaDistribution newTD=TaxaDistBuilder.create(encodeDepths);
        int[] depths=basicLrgTD.depths();
        int[] newDepths=newTD.depths();
        Assert.assertEquals(basicLrgTD.totalDepth(),newTD.totalDepth());
        org.junit.Assert.assertArrayEquals(depths,newDepths);
    }

    @Test
    public void testHighDepthEncoding() throws Exception {
        int[] taxaWithTags={1,5,8,9,10};
        int[] depthTags={1,1,2,255,2};
        TaxaDistribution hdTD=TaxaDistBuilder.create(12,taxaWithTags,depthTags);
        byte[] encodeDepths=hdTD.encodeTaxaDepth();
        TaxaDistribution newTD=TaxaDistBuilder.create(encodeDepths);
        int[] depths=hdTD.depths();
        int[] newDepths=newTD.depths();
        Assert.assertEquals(hdTD.totalDepth(),newTD.totalDepth());
        org.junit.Assert.assertArrayEquals(depths,newDepths);
    }

    @Test
    public void testTaxaSpacingEncoding() throws Exception {
        int[] taxaWithTags={1,256,258,259,280};
        int[] depthTags={1,1,2,255,2};
        TaxaDistribution hdTD=TaxaDistBuilder.create(1025,taxaWithTags,depthTags);
        byte[] encodeDepths=hdTD.encodeTaxaDepth();
        TaxaDistribution newTD=TaxaDistBuilder.create(encodeDepths);
        int[] depths=hdTD.depths();
        int[] newDepths=newTD.depths();
        Assert.assertEquals(hdTD.totalDepth(),newTD.totalDepth());
        org.junit.Assert.assertArrayEquals(depths,newDepths);
    }

}