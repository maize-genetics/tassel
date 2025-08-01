package net.maizegenetics.dna.map;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.Sizeof;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created with IntelliJ IDEA.
 * User: edbuckler
 * Date: 8/2/13
 * Time: 5:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class AnnotatedPositionTest {
    @Test
    public void testBuilder() throws Exception {
        System.out.println("Testing Builder");
        Chromosome aC=new Chromosome("1");
        Position ap= new GeneralPosition.Builder(aC,1232)
                .knownVariants(new String[]{"A","C", "ACGTACA"})
                .maf(0.05f)
                .allele(WHICH_ALLELE.Ancestral,NucleotideAlignmentConstants.C_ALLELE)
                .addAnno("DistToGene", 12323)
                .addAnno("NearestGene","GM1234")
                .addAnno("ImportantGene",true)
                .build();
        assertEquals(ap.getPosition(), 1232);
        assertEquals(ap.isNucleotide(), true);
        assertEquals(ap.getSNPID(), "S1_1232");
        assertEquals(ap.getAllele(WHICH_ALLELE.Ancestral), NucleotideAlignmentConstants.C_ALLELE);
        assertEquals(ap.getAllele(WHICH_ALLELE.Reference), GenotypeTable.UNKNOWN_ALLELE);
        assertEquals(ap.getGlobalMAF(), 0.05f, 0.0001);
        assertEquals(ap.getKnownVariants()[0], "A");
        assertEquals(ap.getKnownVariants()[1], "C");
        assertEquals(ap.getKnownVariants()[2], "ACGTACA");
        double[] tD= ap.getAnnotation().getQuantAnnotation("DistToGene");
        assertEquals(tD[0], 12323.0, 0.0001);
        assertEquals(ap.getAnnotation().getTextAnnotation("NearestGene")[0],"GM1234");
        assertEquals(ap.getAnnotation().getTextAnnotation("ImportantGene")[0],"true");
    }

    @Test
    public void testCompareTo() throws Exception {
        System.out.println("Testing compareTo");
        Chromosome aC=new Chromosome("1");
        Position cp1= new GeneralPosition.Builder(aC,1232).build();
        Position cp2= new GeneralPosition.Builder(aC,1232).build();
        assertEquals(0, cp1.compareTo(cp2));
        cp2= new GeneralPosition.Builder(aC,2000).build();
        assertEquals(-1, cp1.compareTo(cp2));
        cp2= new GeneralPosition.Builder(new Chromosome("-1"),1232).build();
        assertTrue(cp1.compareTo(cp2) > 0);
        //text chromosomes sort to the end
        cp2= new GeneralPosition.Builder(new Chromosome("X"),1232).build();
        assertTrue(cp1.compareTo(cp2) < 0);
    }

    @Ignore
    @Test
    public void testTime() throws Exception {
        int objNum=2_000_000;
        Random r=new Random();
        ArrayList<Position> ps=new ArrayList<Position>(objNum);
        long time=System.nanoTime();
        Chromosome aChr=new Chromosome("1");
        for (int i = 0; i < objNum; i++) {
            Position cp1= new GeneralPosition.Builder(aChr,r.nextInt(200000000))
                    .nucleotide(true)
                    .knownVariants(new String[]{"A", "C"})
                    .maf(0.05f)
                    .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.C_ALLELE)
                    .build();
            ps.add(cp1);
        }

        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)objNum;
        System.out.println("Object Creation Time:" + totalTime / 1e9 + "s  AvgPerObj:" + timePerObj + "ns");
        time=System.nanoTime();
        Collections.sort(ps);
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)objNum;
        System.out.println("Sort Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        for (int i = 0; i < 10; i++) {
            System.out.println(ps.get(i).toString());
        }
    }

    @Ignore
    @Test
    public void testMemory() throws Exception {
        int objNum=2000000;
        long initMem=Sizeof.getMemoryUse();
        Random r=new Random();
        ArrayList<Position> ps=new ArrayList<Position>(objNum);
        long time=System.nanoTime();
        Chromosome aChr=new Chromosome("1");
        for (int i = 0; i < objNum; i++) {
            Position cp1= new GeneralPosition.Builder(aChr,r.nextInt(200000000))
                    .nucleotide(true)
                    .knownVariants(new String[]{"A","C"})
                    .snpName("S1_"+r.nextInt(20000000))
                    .maf(0.05f)
                    .allele(WHICH_ALLELE.Reference,NucleotideAlignmentConstants.C_ALLELE)
                    .allele(WHICH_ALLELE.GlobalMajor,(byte)(i%4))
                    .build();
            ps.add(cp1);
        }

        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)objNum;
        System.out.println("Object Creation Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        time=System.nanoTime();
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)objNum;
        System.out.println("Sort Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        long finalMem=Sizeof.getMemoryUse();
        double memPerAP=(double)(finalMem-initMem)/(double)objNum;
        System.out.println("Memory Per object:"+memPerAP+" bytes");
        for (int i = 0; i < 10; i++) {
            System.out.println(ps.get(i).toString());
        }
    }


    @Test
    public void testGeneralPositionMemory() throws Exception {
        int objNum=2000000;
        long initMem=Sizeof.getMemoryUse();
        Random r=new Random();
        ArrayList<Position> ps=new ArrayList<>(objNum);
        long time=System.nanoTime();
        Chromosome aChr=new Chromosome("1");
        for (int i = 0; i < objNum; i++) {
            int position=r.nextInt(200000000);
            Position cp1= new GeneralPosition.Builder(aChr,position)
                    .nucleotide(true)
                    .knownVariants(new String[]{"A","C"})
 //                   .snpName("S1_"+r.nextInt(20000000)
                    .snpName("S1_"+position)
                    .maf(0.05f)
                    .allele(WHICH_ALLELE.Reference,NucleotideAlignmentConstants.C_ALLELE)
                    .allele(WHICH_ALLELE.GlobalMajor,(byte)(i%4))
//                    .addAnno("NearGene","G123"+r.nextInt(1000))
//                    .addAnno("PromoterType","Element"+r.nextInt(10))
                    .build();
            ps.add(cp1);
        }
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)objNum;
        System.out.println("GeneralPosition Creation Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        time=System.nanoTime();
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)objNum;
        System.out.println("GeneralPosition Sort Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        long finalMem=Sizeof.getMemoryUse();
        double memPerAP=(double)(finalMem-initMem)/(double)objNum;
        System.out.println("GeneralPosition Memory Per object:"+memPerAP+" bytes");
        for (int i = 0; i < 10; i++) {
            System.out.println(ps.get(i).toString()+" Variants:"+ Arrays.toString(ps.get(i).getKnownVariants()));
        }
    }
}
