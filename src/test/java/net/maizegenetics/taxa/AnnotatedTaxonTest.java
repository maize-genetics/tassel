package net.maizegenetics.taxa;


import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import net.maizegenetics.util.GeneralAnnotation;

/**
 *Test Taxon
 */
public class AnnotatedTaxonTest {
    @Test
    public void testBuilder() throws Exception {
        System.out.println("Testing Taxon Builder");
        Taxon at= new Taxon.Builder("Z001E0001:Line:mays:Zea")
                .inbreedF(0.99f)
                .parents("B73","B97")
                .pedigree("(B73xB97)S6I1")
                .addAnno("POOL","BSSS")
                .addAnno("SELFGEN",6)
                .build();
        GeneralAnnotation annotation = at.getAnnotation();
        Assert.assertEquals(annotation.getTextAnnotation(Taxon.MotherKey)[0], "B73");
        Assert.assertEquals(annotation.getTextAnnotation(Taxon.FatherKey)[0], "B97");
        Assert.assertEquals(annotation.getTextAnnotation(Taxon.PedigreeKey)[0], "(B73xB97)S6I1");;
        Assert.assertEquals(at.getName(), "Z001E0001:Line:mays:Zea");
        Assert.assertEquals(annotation.getTextAnnotation("POOL")[0], "BSSS");
        Assert.assertEquals(annotation.getQuantAnnotation("SELFGEN")[0], 6.0,0.001);
        Assert.assertEquals(annotation.getQuantAnnotation(Taxon.InbreedFKey)[0], 0.99f,0.0001);
    }


    @Test
    public void testTime() throws Exception {
        System.out.println("Testing Taxon Sorting and Creation Time");
        int objNum=50000;
        Random r=new Random();
        long time=System.nanoTime();
        ArrayList<Taxon> ps=new ArrayList<Taxon>(objNum);
        for (int i = 0; i < objNum; i++) {
            Taxon at= new Taxon.Builder("Z"+r.nextInt()+":Line:mays:Zea")
                    .inbreedF(0.99f)
                    .parents("B73","B97")
                    .pedigree("(B73xB97)S6I1")
                    .build();
            ps.add(at);
            //psA[i]=cp1;
        }
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)objNum;
        System.out.println("Object Creation Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        time=System.nanoTime();
        Collections.sort(ps);
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)objNum;
        System.out.println("Sort Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        for (int i = 0; i < 10; i++) {
            System.out.println(ps.get(i).toString());
        }
    }

}
