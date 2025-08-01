package net.maizegenetics.dna.map;

import net.maizegenetics.dna.WHICH_ALLELE;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.Arrays;
import java.util.Random;

/**
 *
 */
public class PositionArrayListTest {
    private static byte[] refSeq;
    private static int[] pos;
    private static byte[] chrIndex;
    private static int size=2000000;
    private static Chromosome[] chr;
    private static PositionList instance;

    @BeforeClass
    public static void setUpClass() throws Exception {
        long time=System.nanoTime(), totalTime;
        Random r=new Random(1);
        refSeq=new byte[size];
        chrIndex=new byte[size];
        pos=new int[size];
        chr=new Chromosome[10];
        for (int i = 0; i <chr.length; i++) {chr[i]=new Chromosome(""+i);}
        PositionListBuilder b= new PositionListBuilder();
        for (int i = 0; i <size; i++) {
            refSeq[i]=(byte)r.nextInt(4);
            chrIndex[i]=(byte)(i%10);
            pos[i]=r.nextInt(200000000);
            Position ap=new GeneralPosition.Builder(chr[chrIndex[i]],pos[i]).allele(WHICH_ALLELE.Reference, refSeq[i]).build();
            b.add(ap);
        }
        instance=b.build();
        totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)size;
        System.out.println("Object Creation & Sort Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
    }

    @Test
    public void testSpeedOfPositionLookup() throws Exception {
        long sum=0;
        long time=System.nanoTime(), totalTime;
       // Chromosome chr=new Chromosome("1");
        for (int i=0; i<size; i++) {
            int r=instance.siteOfPhysicalPosition(instance.chromosomalPosition(i),
                    instance.chromosome(i));
            //System.out.printf("pos[i]:%d chr[chrIndex[i]]%s %n",pos[i],chr[chrIndex[i]]);
            Assert.assertTrue(i-r<2);
            sum+=r;
        }
        totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)size;
        System.out.println("SpeedOfPositionLookup Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
    }

    @Test
    public void testGetReferenceAllele() throws Exception {
        for (int i=0; i<10; i++) {
            System.out.println(instance.allele(WHICH_ALLELE.Reference, i));
        }
    }

    @Test
    public void testGetSNPID() throws Exception {
        for (int i=0; i<10; i++) {
            System.out.println(instance.siteName(i));
        }
    }

    @Test
    public void testNullChromosomeGet() throws Exception {
        //if a null chromosome is passed the first chromosome should be used
        int r=instance.siteOfPhysicalPosition(instance.chromosomalPosition(5),
                instance.chromosome(0));
        int nullResult=instance.siteOfPhysicalPosition(instance.chromosomalPosition(5),
                null);
        Assert.assertEquals(r,nullResult);
    }

    @Test
    public void testChromosomeAspects() throws Exception {
        Assert.assertEquals(instance.numChromosomes(),10);
        Assert.assertEquals(instance.chromosome(0).getChromosomeNumber(),0);
        Assert.assertEquals(instance.chromosome(size-1).getChromosomeNumber(),9);
        Assert.assertEquals(instance.chromosomeSiteCount(new Chromosome("0")),size/10);
       // Assert.assertEquals(instance.(new Chromosome("0")),size/10);
        System.out.println("Offsets"+Arrays.toString(instance.chromosomesOffsets()));
    }
}
