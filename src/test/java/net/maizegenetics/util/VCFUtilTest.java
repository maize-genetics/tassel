package net.maizegenetics.util;

import net.maizegenetics.dna.snp.io.VCFUtil;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

/**
 * Unit test for the VCF quality scores
 *
 * @author Ed Buckler
 */
public class VCFUtilTest {
    @Test
    public void testGetScore() throws Exception {
        //currently both of these seem wrong, but I don't know what the answer is supposed to be.
        int[] missingScore= VCFUtil.getScore(0,0);
        System.out.println(Arrays.toString(missingScore));
        int[] lowScore=VCFUtil.getScore(1,1);
        System.out.println(Arrays.toString(missingScore));

    }
    
    @Test
    public void testRefResolve() {
        int[] sortedAlleles = {2,1,4};
        int[] expectedAlleles = {0,2,1,4};
        int[] resolvedAlleles = VCFUtil.resolveRefSorted(sortedAlleles, (byte)0);
        
        Assert.assertArrayEquals( "Failure Unable to Resolve RefAlleles for VCF: Test 1",expectedAlleles, resolvedAlleles );
        
        int[] sortedAlleles2 = {2,1,4};
        int[] expectedAlleles2 = {4,2,1};
        int[] resolvedAlleles2 = VCFUtil.resolveRefSorted(sortedAlleles2, (byte)4);
        
        Assert.assertArrayEquals( "Failure Unable to Resolve RefAlleles for VCF: Test 2",expectedAlleles2, resolvedAlleles2 );
        
        int[] sortedAlleles3 = {2,1,4};
        int[] expectedAlleles3 = {1,2,4};
        int[] resolvedAlleles3 = VCFUtil.resolveRefSorted(sortedAlleles3, (byte)1);
        
        Assert.assertArrayEquals( "Failure Unable to Resolve RefAlleles for VCF: Test3",expectedAlleles3, resolvedAlleles3 );
        
        int[] sortedAlleles4 = {2,1,4};
        int[] expectedAlleles4 = {2,1,4};
        int[] resolvedAlleles4 = VCFUtil.resolveRefSorted(sortedAlleles4, (byte)2);
        
        Assert.assertArrayEquals( "Failure Unable to Resolve RefAlleles for VCF: Test4",expectedAlleles4, resolvedAlleles4 );
    
        //Now to test when using an unknown Reference
        
        
        //Test when Unknown is in the sorted alleles.  Should never be, but could cause some problems.
    }
    
    @Test
    public void testRefResolveUnknownRef() {
        //Testing when we have an Unknown Reference.  This can happen from H5 files if 0x15 is used in the REF file
        int[] sortedAlleles = {2,1,4};
        int[] expectedAlleles = {2,1,4};
        int[] resolvedAlleles = VCFUtil.resolveRefSorted(sortedAlleles, (byte)15);
        
        Assert.assertArrayEquals( "Failure Unable to Resolve Unknown RefAlleles for VCF: Test 1",expectedAlleles, resolvedAlleles );
        
        
    }
    
    @Test 
    public void testIndelInKnownVariant() {
        String[] knownVariants = {"AC","A","AG"};
        Assert.assertTrue("Indel not found in knownVariants: "+Arrays.toString(knownVariants),VCFUtil.indelInKnownVariant(knownVariants));
        
        String[] knownVariants2 = {"C","A","G"};
        Assert.assertFalse("Indel found in knownVariants: "+Arrays.toString(knownVariants),VCFUtil.indelInKnownVariant(knownVariants2));
    
        String[] knownVariants3 = {"A","AC","AG"};
        Assert.assertTrue("Indel not found in knownVariants: "+Arrays.toString(knownVariants3),VCFUtil.indelInKnownVariant(knownVariants3));
    }
    
    /*****
     * comment out speed tests from PL generations
    @Test
    public void testOldVsNewPL(){
        //Enumerate all possibilites for 2 Values
        int[][] depths = new int[16129][2];
        int depthCounter = 0;
        for(int i = 0; i < 127; i++) {
            for(int j = 0; j < 127; j++) {
                depths[depthCounter][0] = i;
                depths[depthCounter][1] = j;
                depthCounter++;
            }
        }
        
        for(int i = 0; i < depths.length; i++) {
            int[] oldVersion = VCFUtil.getScore(depths[i][0], depths[i][1]);
            int[] newVersion = VCFUtil2.getScore(depths[i]);
            //System.out.println(Arrays.toString(depths[i])+": "+Arrays.toString(oldVersion) + ": "+Arrays.toString(newVersion));
            Assert.assertArrayEquals("Arrays Do not match:"+i,newVersion,oldVersion);
        }
    }
    
    @Ignore
    @Test
    public void plSpeedTest() {
        //Make 10million depths and compute the depths to see which is faster
        int[][] depths = new int[10000000][2];
        Random rand = new Random();
        for(int i = 0; i<depths.length; i++) {
            depths[i][0] = rand.nextInt(256);
            depths[i][1] = rand.nextInt(256);
        }
        long startTimeOld = System.currentTimeMillis();
        for(int[] depthArray : depths) {
            VCFUtil.getScore(depthArray[0], depthArray[1]);
        }
        long endTimeOld   = System.currentTimeMillis();
        long totalTimeOld = endTimeOld - startTimeOld;
        
//        long startTimeNew = System.currentTimeMillis();
//        for(int[] depthArray : depths) {
//            VCFUtil2.getScore(depthArray);
//        }
//        long endTimeNew   = System.currentTimeMillis();
//        long totalTimeNew = endTimeNew - startTimeNew;
//        
        long startTimeCache = System.currentTimeMillis();
        for(int[] depthArray : depths) {
            VCFUtil2.getCachedScore(depthArray);
        }
        long endTimeCache   = System.currentTimeMillis();
        long totalTimeCache = endTimeCache - startTimeCache;
        
        
        System.out.println("Total Time for Old: "+totalTimeOld);
       // System.out.println("Total Time for New: "+totalTimeNew);
        System.out.println("Total Time for cache version: "+totalTimeCache);
        
        
    }
    */
    
}
