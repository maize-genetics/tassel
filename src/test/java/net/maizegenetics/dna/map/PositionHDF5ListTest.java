/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.dna.map;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.*;

import java.io.File;
import java.util.Random;

/**
 * Test file for PositionHDF5List
 *
 * @author edbuckler
 */
public class PositionHDF5ListTest {

    private static final Logger myLogger = LogManager.getLogger(PositionHDF5ListTest.class);
    static String testMutFile = GeneralConstants.TEMP_DIR + "testPosGDF5List.hmp.h5";
    static String inFile = TutorialConstants.HAPMAP_FILENAME;
 //   static String inFile = GeneralConstants.DATA_DIR + "/CandidateTests/Ames105v26.chr10.hmp.txt.gz";
     private static GenotypeTable inputAlign;
    private static PositionList instance;

    public PositionHDF5ListTest() {
    }

    @BeforeClass
    public static void setUpClass() {
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
        System.out.println("Reading: " + inFile);
        inputAlign = ImportUtils.readFromHapmap(inFile, null);
        //todo TAS-160 add reference and ancestral alleles
        PositionListBuilder plb=new PositionListBuilder();
        for (Position position : inputAlign.positions()) {
            Position np=(new GeneralPosition.Builder(position))
                    .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                    .allele(WHICH_ALLELE.Ancestral,NucleotideAlignmentConstants.C_ALLELE)
                    .build();
            plb.add(np);
        }
        inputAlign= GenotypeTableBuilder.getInstance(inputAlign.genotypeMatrix(),plb.build(),inputAlign.taxa());
        System.out.println("Writing to mutable HDF: " + testMutFile);
        ExportUtils.writeGenotypeHDF5(inputAlign, testMutFile);
        System.out.println("Opening Mutable HDF5: " + testMutFile);
        instance=PositionListBuilder.getInstance(testMutFile);
    }

    @AfterClass
    public static void tearDownClass() {
        //File f = new File(testMutFile);
        //System.out.println("Deleting:" + testMutFile);
        //  f.delete();
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    private double testPositionLookupSpeed(PositionList a) {
        Random r=new Random(0);
        long time=System.nanoTime();
        int tests=10_000_000;
        long cnt=0;
      //  Chromosome chr=a.chromosome(50);
        for (int i=0; i<tests; i++) {
            int rSite=r.nextInt(a.numberOfSites());
            Chromosome chr=a.chromosome(rSite);
            int pos=a.chromosomalPosition(rSite);
           // int pos=r.nextInt(100_000_000);
            int lSite=a.siteOfPhysicalPosition(pos,chr);
            cnt+=lSite;
          //  Assert.assertEquals(rSite,lSite);
        }
        double rate=(System.nanoTime()-time)/tests;
        System.out.println(cnt);
        return rate;
    }

    @Test
    public void testPositionLookupSpeed() {
        System.out.println("Testing rate of position lookup for:"+inFile);
        System.out.println("Rate:"+testPositionLookupSpeed(inputAlign.positions()));
        System.out.println("Testing rate of position lookup for:"+testMutFile);
        System.out.println("Rate:"+testPositionLookupSpeed(instance));
    }


    @Test
    public void testPositionReferenceAncestral() {
        long time=System.nanoTime();
        for (Position position : instance) {
            Assert.assertEquals(NucleotideAlignmentConstants.T_ALLELE,position.getAllele(WHICH_ALLELE.Reference));
            Assert.assertEquals(NucleotideAlignmentConstants.C_ALLELE,position.getAllele(WHICH_ALLELE.Ancestral));
        }
        System.out.println("Time to cycle through each:"+(System.nanoTime()-time));
        time=System.nanoTime();
        byte[] refAlleles=instance.alleleForAllSites(WHICH_ALLELE.Reference);
        byte[] ancAlleles=instance.alleleForAllSites(WHICH_ALLELE.Ancestral);
        for (int i = 0; i < instance.numberOfSites(); i++) {
            Assert.assertEquals(NucleotideAlignmentConstants.T_ALLELE,refAlleles[i]);
            Assert.assertEquals(NucleotideAlignmentConstants.C_ALLELE,ancAlleles[i]);
        }
        System.out.println("Time to read all REF and ANC alleles:"+(System.nanoTime()-time));
        time=System.nanoTime();
        byte[] mjAlleles=instance.alleleForAllSites(WHICH_ALLELE.GlobalMajor);
        byte[] mnAlleles=instance.alleleForAllSites(WHICH_ALLELE.GlobalMinor);
        for (int i = 0; i < instance.numberOfSites(); i++) {
            Assert.assertEquals(instance.get(i).getAllele(WHICH_ALLELE.GlobalMajor),mjAlleles[i]);
            Assert.assertEquals(instance.get(i).getAllele(WHICH_ALLELE.GlobalMinor),mnAlleles[i]);
        }
        System.out.println("Time to read all MAJOR and MINOR alleles:"+(System.nanoTime()-time));
        //Assert.assertTrue(false);
    }


    @Test
    public void testPositionHDF5ListWithAlignmentV4() {
        System.out.println("Testing HDF5 Alignment...");
        //System.out.println(Arrays.toString(instance.physicalPositions()));
        Assert.assertEquals(inputAlign.numberOfSites(),instance.numberOfSites());
        Assert.assertArrayEquals(inputAlign.physicalPositions(),instance.physicalPositions());
        long time=System.nanoTime();
        Assert.assertEquals(inputAlign.siteName(0),instance.siteName(0));
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)inputAlign.numberOfSites();
        System.out.println("Cache Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        time=System.nanoTime();
        for (int ii=0; ii<inputAlign.numberOfSites()*100; ii++) {
            int i=ii%inputAlign.numberOfSites();
            Assert.assertEquals(inputAlign.chromosomalPosition(i),instance.chromosomalPosition(i));  //209ns
            Assert.assertEquals(inputAlign.chromosomeName(i),instance.chromosomeName(i));      //158ns
            Assert.assertEquals(inputAlign.siteName(i),instance.siteName(i));                        //456ns
            //Assert.assertEquals(inputAlign.majorAllele(i),instance.referenceAllele(i));                        //456ns
        }
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)(inputAlign.numberOfSites()*100);
        System.out.println("Post Cache Comparison Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        time=System.nanoTime();
        for (int ii=0; ii<inputAlign.numberOfSites()*100; ii++) {
            int i=ii%inputAlign.numberOfSites();
            Position ap=instance.get(i);   //136
            Assert.assertEquals(inputAlign.chromosomalPosition(i),ap.getPosition());  //174ns
            Assert.assertEquals(inputAlign.chromosome(i).getName(),ap.getChromosome().getName());      //157ns
            Assert.assertEquals(inputAlign.siteName(i),ap.getSNPID());                        //174ns
            Assert.assertEquals(inputAlign.majorAllele(i),ap.getAllele(WHICH_ALLELE.GlobalMajor));                        //153ns
            //>50% of this time is actually from the old in memory alignment
        }
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)(inputAlign.numberOfSites()*100);
        System.out.println("Post Cache Grab AP Comparison Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
    }

}
