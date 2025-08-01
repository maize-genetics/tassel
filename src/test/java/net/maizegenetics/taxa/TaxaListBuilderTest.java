package net.maizegenetics.taxa;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.ImportUtils;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.Random;

/**
 * Tests
 */
public class TaxaListBuilderTest {
    static String testMutFile = GeneralConstants.TEMP_DIR + "test.hmp.h5";
    private static GenotypeTable inputAlign;

    @Before
    public void setUp() throws Exception {
        File tempDir = new File(GeneralConstants.TEMP_DIR);
        tempDir.mkdir();
        File f = new File(testMutFile);
        System.out.println(f.getAbsolutePath());
        if (f.getParentFile().listFiles() != null) {
            for (File file : f.getParentFile().listFiles()) {
                if (file.getName().endsWith("hmp.h5")) {
                    file.delete();
                    System.out.println("Deleting:" + file.toString());
                }
            }
        }
        System.out.println("Reading: " + TutorialConstants.HAPMAP_FILENAME);
        inputAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
        System.out.println("Writing to mutable HDF: " + testMutFile);
        ExportUtils.writeGenotypeHDF5(inputAlign, testMutFile);

    }

    @After
    public void tearDown() throws Exception {

    }

    @Test
    public void testBuilder() throws Exception {
        System.out.println("Testing TaxaListBuilder with de novo creation");
        int objNum=50000;
        Random r=new Random();
        long time=System.nanoTime();
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < objNum; i++) {
            Taxon at= new Taxon.Builder("Z"+(objNum-i)+":Line:mays:Zea")
                    .inbreedF(0.99f)
                    .parents("B73","B97")
                    .pedigree("(B73xB97)S6I1")
                    .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)objNum;
        System.out.println("Object Creation Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        for (int i = 0; i < 10; i++) {
            System.out.println(tl.get(i).toString());
        }
        time=System.nanoTime();
        tl=tlb.sortTaxaAlphabetically().build();
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)objNum;
        System.out.println("Sort Time:" + totalTime / 1e9 + "s  AvgPerObj:" + timePerObj + "ns");
        for (int i = 0; i < 10; i++) {
            System.out.println(tl.get(i).toString());
        }
    }

    @Test
    public void testBuilderMergeTaxonAnnotation() throws Exception {
        System.out.println("Testing TaxaListBuilder with de novo creation");
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 10; i++) {
            Taxon at= new Taxon.Builder("Z"+i+":Line:mays:Zea")
                    .inbreedF(0.99f)
                    .parents("B73","B97")
                    .pedigree("(B73xB97)S6I1")
                    .build();
            tlb.add(at);
        }
        for (int i = 0; i < 10; i++) {
            Taxon at= new Taxon.Builder("Z"+i+":Line:mays:Zea")
                    .inbreedF(0.99f)
                    .parents("B73","B97")
                    .pedigree("MessedUpPedigree")
                    .build();
            tlb.addOrMerge(at);
        }
        TaxaList tl=tlb.build();
        for (int i = 0; i < 10; i++) {
            Assert.assertEquals(tl.get(i).getAnnotation().getQuantAnnotation("INBREEDF").length,1);
            Assert.assertEquals(tl.get(i).getAnnotation().getTextAnnotation("PEDIGREE").length,2);
            //System.out.println(tl.get(i).toStringWithVCFAnnotation());
        }
    }


    @Test
    public void testTaxonListAgainstAlignment() throws Exception {
        System.out.println("Testing TaxaListBuilder with Existing Alignment");
        TaxaList tl=new TaxaListBuilder().addAll(inputAlign).build();
        Assert.assertEquals(inputAlign.numberOfTaxa(),tl.numberOfTaxa());
        Assert.assertEquals(inputAlign.numberOfTaxa(),tl.size());

        for (int i = 0; i < tl.size(); i++) {
            Assert.assertEquals(inputAlign.taxaName(i),tl.taxaName(i));
            Assert.assertEquals(inputAlign.taxaName(i),tl.taxaName(i));
            Assert.assertEquals((int)i,(int)tl.indexOf(inputAlign.taxaName(i)));
        }
    }

    @Test
    public void testTaxonListWithDuplicatedTaxa() throws Exception {
        System.out.println("Testing TaxaListBuilder with Existing Alignment");
        TaxaListBuilder tlb=new TaxaListBuilder().addAll(inputAlign);
        for (Taxon taxon : inputAlign.taxa()) {
            try{
                tlb.add(taxon);
                Assert.assertTrue("Duplicate taxon was added",false);
            } catch(IllegalStateException e) {
                //duplicate caught appropriate
                Assert.assertTrue("Duplicate taxon was rejected",true);
            }
        }
        TaxaList tl= tlb.build();

        Assert.assertEquals(inputAlign.numberOfTaxa(),tl.numberOfTaxa());
        Assert.assertEquals(inputAlign.numberOfTaxa(),tl.size());
        for (int i = 0; i < tl.size(); i++) {
            Assert.assertEquals(inputAlign.taxaName(i),tl.taxaName(i));
            Assert.assertEquals(inputAlign.taxaName(i),tl.taxaName(i));
            Assert.assertEquals((int)i,(int)tl.indexOf(inputAlign.taxaName(i)));
        }
    }

    @Test
    public void testTaxonListFromHDF5() throws Exception {
        System.out.println("Testing TaxaListBuilder with HDF5 Alignment");
        System.out.println("Reading: " + TutorialConstants.HAPMAP_FILENAME);
        TaxaList tl=new TaxaListBuilder().buildFromHDF5Genotypes(HDF5Factory.openForReading(testMutFile));
        TaxaList inputAlignTL=new TaxaListBuilder().addAll(inputAlign).build();      //cannot directly compare as the sort on the input align is not the same as java comparable
        Assert.assertEquals(inputAlignTL.numberOfTaxa(),tl.numberOfTaxa());
        Assert.assertEquals(inputAlignTL.numberOfTaxa(),tl.size());

        for (int i = 0; i < tl.size(); i++) {
            Assert.assertEquals(inputAlignTL.taxaName(i),tl.taxaName(i));
            Assert.assertEquals(inputAlignTL.taxaName(i),tl.taxaName(i));
            Assert.assertEquals(i,tl.indexOf(inputAlignTL.taxaName(i)));
        }
    }


}
