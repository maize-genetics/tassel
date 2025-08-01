package net.maizegenetics.dna.snp;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.io.BuilderFromHapMap;
import net.maizegenetics.dna.snp.io.ProjectionGenotypeIO;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.DonorHaplotypes;
import net.maizegenetics.taxa.Taxon;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.NavigableSet;
import java.util.Random;
import java.util.TreeSet;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class ProjectionGenotypeTest {
    private GenotypeTable projInst;
    private GenotypeTable homoOrig;
    static String inFile = TutorialConstants.HAPMAP_FILENAME;
//    static String inFile = GeneralConstants.DATA_DIR + "/CandidateTests/Ames105v26.chr10.hmp.txt.gz";
    static String testProjFile = GeneralConstants.TEMP_DIR + "testProj.pa.txt";

    @Before
    public void setUp() throws Exception {
        File tempDir = new File(GeneralConstants.TEMP_DIR);
        tempDir.mkdir();
        File f = new File(testProjFile);
        System.out.println(f.getAbsolutePath());
        if (f.getParentFile().listFiles() != null) {
            for (File file : f.getParentFile().listFiles()) {
                if (file.getName().endsWith("pa.txt")) {
                    file.delete();
                    System.out.println("Deleting:" + file.toString());
                }
            }
        }

        GenotypeTable a=BuilderFromHapMap.getBuilder(inFile).build();
        homoOrig=a=GenotypeTableBuilder.getHomozygousInstance(a);   //only homozygous sites will be identical

        ProjectionBuilder pb=new ProjectionBuilder(a);
        for (Taxon taxon : a.taxa()) {
            NavigableSet<DonorHaplotypes> bp=new TreeSet<>();
            int taxonIndex=a.taxa().indexOf(taxon);
            for (Chromosome chromosome : a.chromosomes()) {
                int[] se=a.firstLastSiteOfChromosome(chromosome);
                int endPosition=a.chromosomalPosition(se[1]);
                DonorHaplotypes dh=new DonorHaplotypes(chromosome,0,endPosition,taxonIndex,taxonIndex);
                bp.add(dh);
            }
            pb.addTaxon(taxon,bp);
        }
        projInst=pb.build();
        for (int s = 0; s < homoOrig.numberOfSites(); s++) {
            for (int t = 0; t < homoOrig.numberOfTaxa(); t++) {
                byte o=homoOrig.genotype(t,s);
                byte p=projInst.genotype(t,s);
                if(GenotypeTableUtils.isHeterozygous(p)) System.out.println("s:"+s+" t:"+t+" o:"+o+" p:"+p);
            }
        }

//        try{
//            Thread.sleep(100000);
//        } catch(Exception e) {}
    }

    private GenotypeTable createLargeProjAlignment(int size, double F) {
        GenotypeTable a=BuilderFromHapMap.getBuilder(inFile).build();
        a=GenotypeTableBuilder.getHomozygousInstance(a);   //only homozygous sites will be identical
       // a.genotypeMatrix().transposeData(false);
        //a=GenotypeTableBuilder.

        ProjectionBuilder pb=new ProjectionBuilder(a);
        Random r=new Random();
        for (int i = 0; i < size; i++) {
            Taxon taxon=new Taxon("T"+i);
            NavigableSet<DonorHaplotypes> bp=new TreeSet<>();
            for (Chromosome chromosome : a.chromosomes()) {
                int[] se=a.firstLastSiteOfChromosome(chromosome);
                int endPosition=a.chromosomalPosition(se[1]);
                int taxonIndex=r.nextInt(a.numberOfTaxa());
                int taxon2Index=r.nextInt(a.numberOfTaxa());
                if(r.nextDouble()<F)  taxon2Index=taxonIndex;
                DonorHaplotypes dh=new DonorHaplotypes(chromosome,0,endPosition,taxonIndex,taxon2Index);
                bp.add(dh);
            }
            pb.addTaxon(taxon,bp);
        }
        return pb.build();
    }


    @Test
    public void testIdentityOfAlignments() throws Exception {
        System.out.println("Alignment.UNKNOWN_DIPLOID_ALLELE:"+GenotypeTable.UNKNOWN_DIPLOID_ALLELE);
        System.out.println(NucleotideAlignmentConstants.getNucleotideIUPAC((byte)16));
        System.out.println("testIdentityOfAlignments");
        AlignmentTestingUtils.alignmentsEqual(homoOrig, projInst);

    }

    @Test
    public void testRates() throws Exception {
        System.out.println("Evaluate speed in general mode");
        System.out.println("...accessed in taxa mode");
        scoreRateOfAccessingBases(homoOrig, projInst);
        System.out.println("...accessed in site mode");
        scoreRateOfAccessingBasesSites(homoOrig, projInst);
        System.out.println();
//        System.out.println("Evaluate speed in taxa mode");
//        System.out.println("...accessed in taxa mode");
//        projInst.genotypeMatrix().transposeData(false);
//        scoreRateOfAccessingBases(homoOrig, projInst);
//        System.out.println("...accessed in site mode");
//        scoreRateOfAccessingBasesSites(homoOrig, projInst);
        System.out.println();
        System.out.println("Evaluate speed in site mode");
        System.out.println("...accessed in taxa mode");
        projInst.genotypeMatrix().transposeData(true);
        scoreRateOfAccessingBases(homoOrig, projInst);
        System.out.println("...accessed in site mode");
        scoreRateOfAccessingBasesSites(homoOrig, projInst);
        System.out.println();

        System.out.println("Evaluate speed of superbyte alignment comparison to itself");
        System.out.println("...accessed in taxa mode");
        homoOrig.genotypeMatrix().transposeData(false);
        scoreRateOfAccessingBases(homoOrig, homoOrig);
        System.out.println("...accessed in site mode");
        homoOrig.genotypeMatrix().transposeData(true);
        scoreRateOfAccessingBasesSites(homoOrig, projInst);

        System.out.println();
        System.out.println("Evaluate speed of large projection");
        GenotypeTable largeA=createLargeProjAlignment(5000,1.8);
        largeA.genotypeMatrix().transposeData(false);
        for (int i=0; i<2; i++) {
            System.out.println("...General Access Round "+i);
            scoreSpeedInGWASContext(largeA);
        }
        largeA.genotypeMatrix().transposeData(true);
        for (int i=0; i<2; i++) {
            System.out.println("...Site Access Round "+i);
            scoreSpeedInGWASContext(largeA);
        }

    }

    @Test
    public void testSavingAndReading() throws Exception {
        System.out.println("testSavingAndReading");
        ProjectionGenotypeIO.writeToFile(testProjFile, projInst);

        GenotypeTable fromFile=ProjectionGenotypeIO.getInstance(testProjFile, inFile);
        for (int s = 0; s < homoOrig.numberOfSites(); s++) {
            for (int t = 0; t < homoOrig.numberOfTaxa(); t++) {
                byte o=homoOrig.genotype(t,s);
                byte p=fromFile.genotype(t,s);
                if(GenotypeTableUtils.isHeterozygous(p)) System.out.println("s:"+s+" t:"+t+" o:"+o+" p:"+p);
            }
        }
        AlignmentTestingUtils.alignmentsEqual(homoOrig, fromFile);
        AlignmentTestingUtils.alignmentsEqual(projInst, fromFile);


    }

    private double scoreRateOfAccessingBases(GenotypeTable a, GenotypeTable b) {
        long same=0, diff=0;
        long time=System.nanoTime();
        for (int i=0; i<a.numberOfTaxa(); i++) {
            for (int j=0, n=a.numberOfSites(); j<n; j++) {
               if(a.genotype(i,j)==b.genotype(i,j)) {same++;} else {diff++;}
            }
        }
        long totalTime=System.nanoTime()-time;
        double rate=(double)totalTime/(double)(same+diff);
        System.out.printf("Same:%d Diff:%d Time:%d ms Rate:%g ns/bp %n",same, diff, totalTime/1_000_000, rate);
        return rate;
    }

    private double scoreRateOfAccessingBasesSites(GenotypeTable a, GenotypeTable b) {
        long same=0, diff=0;
        long time=System.nanoTime();
        for (int i=0; i<a.numberOfSites(); i++) {
            for (int j=0, n=a.numberOfTaxa(); j<n; j++) {
                 if(a.genotype(j,i)==b.genotype(j,i)) {same++;} else {diff++;}
            }
        }
        long totalTime=System.nanoTime()-time;
        double rate=(double)totalTime/(double)(same+diff);
        System.out.printf("Same:%d Diff:%d Time:%d ms Rate:%g ns/bp %n",same, diff, totalTime/1_000_000, rate);
        return rate;
    }

    private double scoreSpeedInGWASContext(GenotypeTable b) {
        long major=0, diff=0;
        long time=System.nanoTime();
        double[] traits=new double[b.numberOfTaxa()];
        Random random=new Random();
        for (int i=0; i<traits.length; i++) {
            traits[i]=random.nextGaussian();
        }
        double sSS=0, dSS=0;
        for (int i=0; i<b.numberOfSites(); i++) {
            byte first=b.genotype(0, i);
            for (int j=0, n=b.numberOfTaxa(); j<n; j++) {
                if(first==b.genotype(j,i))
                    {major++; sSS+=traits[j]*traits[j];} else {diff++; dSS+=traits[j]*traits[j];}
            }
        }
        long totalTime=System.nanoTime()-time;
        double rate=(double)totalTime/(double)(major+diff);
        System.out.printf("Major:%d NotMajor:%d Time:%d ms Rate:%g ns/bp %n",major, diff, totalTime/1_000_000, rate);
        System.out.printf("sSS:%g dSS:%g %n", sSS,dSS);
        return rate;
    }

    @Test
    public void testCreateProjectionAlignment() throws Exception {


    }
}
