package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.map.*;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.Arrays;
import java.util.Random;

/**
 * Provides testing for DynamicBitStorage conversion and storage.  There are routines to tests the speed of conversion, and
 * ensure agreement with the known results.
 * <p></p>
 * There are also several Ignored tests that compare the speed of transposing byte[taxa][site] to byte[site][taxa].
 * Additionally there are booleans flag that can be set to compare and evaluate speed in writing and reading to
 * HDF5 files.
 *
 * @author Ed Buckler
 */
public class GenotypeCalcSpeedTest {
    private static int taxa=500; //Timing tests for 500
    private static int sites=1_000_000;//Timing test for 1Million
    private static byte[][] myGenotype; //[taxa][sites]
    private static GenotypeTable myAlign;
    private static byte[] myMajor; //[taxa][sites]
    private static byte[] myMinor; //[taxa][sites]
    private static double baseLineCPUTime;


    @BeforeClass
    public static void setUp() throws Exception {
        long time=System.nanoTime();
        System.out.println("Creating GenotypeCalcSpeedTest test data");
        myGenotype=new byte[taxa][sites];
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < taxa; i++) {
            Taxon at= new Taxon.Builder("Z"+i)
                    .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();
        Chromosome chr=new Chromosome("1");
        PositionListBuilder b= new PositionListBuilder();
        GenotypeCallTableBuilder gb=GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(taxa,sites);
        myMajor=new byte[sites];
        myMinor=new byte[sites];
        Random r=new Random(0);
        for (int s=0; s<sites; s++) {
            myMajor[s]=(byte)r.nextInt(4);
            myMinor[s]=(byte)r.nextInt(4);
            if(myMajor[s]==myMinor[s]) {myMinor[s]=GenotypeTable.UNKNOWN_DIPLOID_ALLELE;}
            Position ap=new GeneralPosition.Builder(chr,s*10)
                    .allele(WHICH_ALLELE.GlobalMajor,myMajor[s])
                    .allele(WHICH_ALLELE.GlobalMinor,myMinor[s]).build();
            b.add(ap);
        }
        for (int t=0; t<taxa; t++) {
            for (int s=0; s<sites; s++) {
                if(r.nextInt(100)<10) {myGenotype[t][s]=GenotypeTableUtils.getUnphasedDiploidValue(myMinor[s], myMinor[s]);}
                else {myGenotype[t][s]=GenotypeTableUtils.getUnphasedDiploidValue(myMajor[s],myMajor[s]);}
                gb.setBase(t,s,myGenotype[t][s]);
            }
        }
        PositionList apl=b.build();
        GenotypeCallTable g=gb.build();
       // BitStorage bst2=new DynamicBitStorage(g, Alignment.ALLELE_SCOPE_TYPE.Frequency, myMajor, myMinor);
        myAlign=GenotypeTableBuilder.getInstance(g, apl, tl, null, null, null, null, null);
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("Completed GenotypeCalcSpeedTest test data in "+timePerObj+"ns/site");
        baseLineCPUTime=timePerObj;

    }

    @Test
    public void testBasicSpeed() throws Exception {
        int[][] af=calculateAlleleFreqThroughInterface(myAlign);
        for (int i=0; i<af.length; i++) {
            System.out.println(af[i][0]);
        }
        int[][] af2=calculateAlleleFreq(myAlign);
        for (int i=0; i<af2.length; i++) {
            System.out.println(af2[i][0]);
        }

    }

    private static int[][] calculateAlleleFreq(GenotypeTable ma) {
        long time=System.nanoTime();
        int sites=ma.numberOfSites();
        int taxa=ma.numberOfTaxa();
        int[][] af=new int[6][sites];
        byte[][] afOrder=new byte[6][sites];
        float[] coverage=new float[ma.numberOfTaxa()];
        float[] hets=new float[ma.numberOfTaxa()];
        for (int taxon = 0; taxon < taxa; taxon++) {
            byte[] genotype=ma.genotypeAllSites(taxon);
            int covSum=0;  //coverage of the taxon
            int hetSum=0;
            for (int s = 0; s < sites; s++) {
                byte[] b = GenotypeTableUtils.getDiploidValues(genotype[s]);
                if(b[0]<af.length) af[b[0]][s]++;
                if(b[1]<af.length) af[b[1]][s]++;
                if(GenotypeTableUtils.isHeterozygous(genotype[s])) hetSum++;
                if(genotype[s]!=GenotypeTable.UNKNOWN_DIPLOID_ALLELE) covSum++;
            }
            coverage[taxon]=(float)covSum/(float)sites;
            hets[taxon]=(float)hetSum/(float)covSum;
        }
        float[] maf=new float[sites];
        float[] paf=new float[sites];
        int baseMask=0xF;
        for (int s = 0; s < sites; s++) {
            int sum=0;
            int[] cntAndAllele=new int[6];
            for (byte i = 0; i < 6; i++) {
                cntAndAllele[i]=(af[i][s]<<4)|(5-i);  //size | allele (the 5-i is to get the sort right, so if case of ties A is first)
                sum+=af[i][s];
            }
            Arrays.sort(cntAndAllele);  //ascending quick sort, there are faster ways
            //http://stackoverflow.com/questions/2786899/fastest-sort-of-fixed-length-6-int-array
            for (byte i = 0; i < 6; i++) {
                afOrder[5-i][s]=(cntAndAllele[i]>0xF)?((byte)(5-(baseMask&cntAndAllele[i]))):GenotypeTable.UNKNOWN_ALLELE;
            }
            if(afOrder[1][s]!=GenotypeTable.UNKNOWN_ALLELE) maf[s]=(float)af[afOrder[1][s]][s]/(float)sum;
            paf[s]=(float)sum/(float)(2*taxa);
        }
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("calculateAlleleFreq Access Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        return af;
        //ma.setCalcAlleleFreq(af, afOrder, maf, paf, coverage, hets);
    }

    private static int[][] calculateAlleleFreqThroughInterface(GenotypeTable ma) {
        long time=System.nanoTime();
        int sites=ma.numberOfSites();
        int taxa=ma.numberOfTaxa();
        byte[][] afOrder=new byte[6][sites];
        int[][] af=new int[6][sites];
        for (int i=0; i<sites; i++) {
            int[][] laf=ma.allelesSortedByFrequency(i);
            for (int j=0; j<laf[0].length; j++) {
                afOrder[j][i]=(byte)laf[0][j];
                af[j][i]=laf[1][j];
            }
        }
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("calculateAlleleFreqThroughInterface Access Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        return af;
        //ma.setCalcAlleleFreq(af, afOrder, maf, paf, coverage, hets);
    }

}
