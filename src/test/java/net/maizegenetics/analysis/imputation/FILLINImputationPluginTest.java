/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.analysis.imputation;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.taxa.Taxon;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.*;

import java.io.File;
import java.util.Random;

import static org.junit.Assert.assertTrue;

/**
 * Basic approach of this test is create a series of test genotypes based on the mdp_genotype.hmp.txt (maize 282
 * association panel). The test creates mapping populations, synthetics, and association panels.  It then tests
 * a range of parameter to optimize, and it evaluates when the error rate is sufficently low.
 * @author edbuckler
 */
public class FILLINImputationPluginTest {
    private static final Logger myLogger = LogManager.getLogger(FILLINImputationPluginTest.class);
    static String testMutFile=GeneralConstants.TEMP_DIR+"test.hmp.h5";
    static String baseFile=GeneralConstants.TEMP_DIR+"FMunmaskedbase.hmp.h5";
    static String maskedFile=GeneralConstants.TEMP_DIR+"FMmasked.hmp.h5";
    static String impFile=GeneralConstants.TEMP_DIR+"FMmasked.imp.hmp.h5";
    static String haploFiles=GeneralConstants.TEMP_DIR+"FMhaps.gX.hmp.txt.gz";
    public static GenotypeTable inputAlign;
 //   static MutableNucleotideAlignmentHDF5 instance;

    public FILLINImputationPluginTest() {

    }

    
    private static void clearTestFolder() {
        File f=new File(GeneralConstants.TEMP_DIR);
        System.out.println("Deleting files from "+GeneralConstants.TEMP_DIR);
        //System.out.println(f.getAbsolutePath());
        if(f.listFiles()!=null) {
            for (File file : f.listFiles()) {
             if(file.getName().contains("hmp.")) {
                 file.delete();
                 //System.out.println("Deleting:"+file.toString());
             }
            }
         }
    }
    
    public static void createTestDataForImputation(int numberOfOrig, int numberOfRILs, int numberOfRILparents, boolean createHetChromosomes) {
        //Create masked dataset
        Random r=new Random(1);
        System.out.println("Reading: "+TutorialConstants.HAPMAP_FILENAME);
        File f=new File(TutorialConstants.HAPMAP_FILENAME);
        System.out.println("Absolute Reading: "+f.getAbsolutePath());
        GenotypeTable origAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
        GenotypeTableBuilder ab=GenotypeTableBuilder.getTaxaIncremental(origAlign.positions());
        int origNumTaxa=origAlign.numberOfTaxa();
        int numberOfStates=(createHetChromosomes)?3:2;
        for (int i=0; i<numberOfOrig; i++) {
            ab.addTaxon(origAlign.taxa().get(i),origAlign.genotypeAllSites(i));
        }
        for (int nt = numberOfOrig /*set above zero to keep orig*/; nt < numberOfOrig+numberOfRILs; nt++) {
            int parent1=r.nextInt(numberOfRILparents);
            int parent2=r.nextInt(numberOfRILparents);
            Taxon tx=new Taxon(origAlign.taxaName(parent2)+origAlign.taxaName(parent1)+"_"+nt);
            int modeCnt;
            byte[] genotype=new byte[origAlign.numberOfSites()];
            for (Chromosome aL : origAlign.chromosomes()) {
//                System.out.printf("Chr:%s Mode:%d %n", aL.getName(), modeCnt);
                modeCnt=r.nextInt(numberOfStates);
                int[] startEnd=origAlign.positions().startAndEndOfChromosome(aL);
                for (int i = startEnd[0]; i <= startEnd[1]; i++) {
                    switch (modeCnt) {
                        case 0: genotype[i]=origAlign.genotype(parent2, i); break;
                        case 1: genotype[i]=origAlign.genotype(parent1, i); break;
                        case 2: genotype[i]= GenotypeTableUtils.getUnphasedDiploidValueNoHets(origAlign.genotype(parent2, i),
                                origAlign.genotype(parent1, i)); break;
                    }
                }
            }
            ab.addTaxon(tx,genotype);
        }
        GenotypeTable baseA=ab.build();
        ExportUtils.writeGenotypeHDF5(baseA, baseFile);
        ExportUtils.writeToHapmap(baseA, false, baseFile, '\t', null);
        GenotypeCallTableBuilder gBMask=GenotypeCallTableBuilder.getInstance(baseA.numberOfTaxa(),baseA.numberOfSites());
        for (int t = 0; t < baseA.numberOfTaxa(); t++) {
            for (int s = 0; s < baseA.numberOfSites(); s++) {
                if(r.nextDouble()<0.3) {gBMask.setBase(t, s, GenotypeTable.UNKNOWN_DIPLOID_ALLELE);}
                else {gBMask.setBase(t, s, baseA.genotype(t,s));}
            }
        }
        GenotypeTable maskA=GenotypeTableBuilder.getInstance(gBMask.build(),baseA.positions(), baseA.taxa());
        ExportUtils.writeToHapmap(maskA, false, maskedFile, '\t', null);
        ExportUtils.writeGenotypeHDF5(maskA, maskedFile);
        
    }

    @Test
    public void testRunFindMergeHaplotypesPlugin() {
        double[] bestError=new double[4];
      //  double[] expError={0.0015,0.055,0.01,0.0001};
        double[] expError={0.0016,0.1,0.011,0.0001};
        int testNumber=0;
        clearTestFolder();
        System.out.println("Starting Mixed diverse and 10 line synthetic");
        createTestDataForImputation(100,600,10,true);
        createHaplotypes();
        bestError[testNumber]=imputeMissing(25, 26, 0.01, 0.025);
        System.out.println("Mixed diverse and 10 line synthetic bestError:"+bestError);
        assertTrue("Mixed diverse and 10 line synthetic", bestError[testNumber]<expError[testNumber]);
        testNumber++;
        
        clearTestFolder();
        System.out.println("Starting Diverse Inbreds No Reps");
        createTestDataForImputation(281,0,10,true);
        createHaplotypes();
        bestError[testNumber]=imputeMissing(25, 26, 0.01, 0.025);
        System.out.println("Diverse Inbreds No Reps bestError:"+bestError);
        assertTrue("Diverse Inbreds No Reps", bestError[testNumber]<expError[testNumber]);
        testNumber++;
        
        clearTestFolder();
        System.out.println("Starting Diverse Inbreds 3-reps");
        createTestDataForImputation(281,281*3,281,true);
        createHaplotypes();
        bestError[testNumber]=imputeMissing(25, 26, 0.01, 0.025);
        System.out.println("Diverse Inbreds 3-reps bestError:"+bestError);
        assertTrue("Diverse Inbreds 3-reps", bestError[testNumber]<expError[testNumber]);
        testNumber++;
        
        clearTestFolder();
        System.out.println("Starting 10 lines synthetic");
        createTestDataForImputation(10,600,10,true);
        createHaplotypes();
        bestError[testNumber]=imputeMissing(25, 26, 0.01, 0.025);
        System.out.println("RILs bestError:"+bestError);
        assertTrue("10 lines synthetic", bestError[testNumber]<expError[testNumber]);
        testNumber++;

        System.out.println("Summary of tests");
        for (int i=0; i<bestError.length; i++) {
            System.out.printf("Test %d ExpectedError %.5g Observed %.5g \n",i,expError[i],bestError[i]);
        }


    }
    
    private void createHaplotypes() {
        String[] args2 = new String[]{
            "-hmp", maskedFile,
            "-o", haploFiles,
            "-mxDiv", "0.1",  //these values are much to high for GBS
            "-mxHet", "0.05",  //these values are much to high for GBS
            "-hapSize", "200",
            "-minPres", "50",
            "-maxOutMiss","0.4",
            "-maxHap", "100",
            "-nV"
        };
        FILLINFindHaplotypesPlugin plugin = new FILLINFindHaplotypesPlugin();
        plugin.setParameters(args2);
        plugin.performFunction(null);
    }
    
    private double imputeMissing(int minMinMinor, int maxMinMinor, double minmxInbErr, double maxmxInbErr) {
        double minError=1;
        int fileCnt=0;
        for (int mm = minMinMinor; mm < maxMinMinor; mm+=5) {
            for (double mie = minmxInbErr; mie < maxmxInbErr; mie+=0.01) {
            int minMinor=mm;//15;
            double mxInbErr=mie;//0.02; 
            String currImpFile=impFile.replace("FMmask", "FM"+fileCnt+"mask");
            String[] args3 = new String[]{
                "-hmp", maskedFile,
                "-d",haploFiles,
                "-o", currImpFile,
//                "-sC","1",
//                "-eC","1",
                "-mxHet","0.02",
                "-minMnCnt", ""+minMinor,
                "-mxInbErr",""+mxInbErr,
                "-mxHybErr","0.02",
                "-mnTestSite","20",
                "-mxDonH","10",
     //           "-projA",
                "-nV"
            };
            System.out.println("TasselPrefs:"+TasselPrefs.getAlignmentRetainRareAlleles());
            TasselPrefs.putAlignmentRetainRareAlleles(false);
            FILLINImputationPlugin impPlugin = new FILLINImputationPlugin();
            impPlugin.setParameters(args3);
            impPlugin.performFunction(null);
            
            int[] results=FILLINImputationAccuracy.compareAlignment(baseFile, maskedFile, currImpFile, false);
            System.out.println("Result\tminMinor\tmxInbErr\tGap\tUnimp\tUnimpHets\tCorrect\tErrors\tRate");
            double errorRate=(double)results[4]/(double)(results[3]+results[4]);
            System.out.printf("AccSummary\t%d\t%g\t%d\t%d\t%d\t%d\t%d\t%g%n",minMinor, mxInbErr, results[0],results[1],results[2],results[3],results[4],errorRate);
            fileCnt++;
            if(errorRate<minError) minError=errorRate;
            }
        }
        return minError;
    }


    
}
