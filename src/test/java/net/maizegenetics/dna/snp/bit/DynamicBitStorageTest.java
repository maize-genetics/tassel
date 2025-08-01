package net.maizegenetics.dna.snp.bit;

import ch.systemsx.cisd.hdf5.*;
import junit.framework.Assert;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.HapMapHDF5Constants;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.map.*;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.OpenBitSet;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.io.File;
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
public class DynamicBitStorageTest {
    private static int taxa=500; //Timing tests for 500
    private static int sites=1_000_000;//Timing test for 1Million
    private static byte[][] myGenotype; //[taxa][sites]
    private static GenotypeTable myAlign;
    private static byte[] myMajor; //[taxa][sites]
    private static byte[] myMinor; //[taxa][sites]
    private static final int defLevel=2;
    private static final boolean testHDF5Byte=false;
    private static final boolean testHDF5Bit=false;
    private static double baseLineCPUTime;


    @BeforeClass
    public static void setUp() throws Exception {
        long time=System.nanoTime();
        System.out.println("Creating DynamicBitStorageTest test data");
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
        System.out.println("Completed DynamicBitStorageTest test data in "+timePerObj+"ns/site");
        baseLineCPUTime=timePerObj;
        String outFile=GeneralConstants.TEMP_DIR+"raw.h5";
        if(testHDF5Byte) writeBytesToTestHDF5(outFile);
        if(testHDF5Byte) readBytesToTestHDF5(outFile);
    }

    @Ignore
    @Test
    public void testGetAllelePresenceForAllSites() throws Exception {
       long time=System.nanoTime();
       long mjCard=0, mnCard=0;
       for (int t=0; t<taxa; t++) {
           mjCard+=myAlign.allelePresenceForAllSites(t, WHICH_ALLELE.Major).cardinality();
           mnCard+=myAlign.allelePresenceForAllSites(t, WHICH_ALLELE.Minor).cardinality();
        }
        System.out.printf("Cardinality1 Mj:%d Mn:%d %n", mjCard, mnCard);
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("TBit Creation & Access Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        double relativeRate=baseLineCPUTime/timePerObj;
        Assert.assertTrue("TBit Creation Time Longer Than Expected:"+timePerObj+"ns and relative rate:"+relativeRate, relativeRate>4);
        time=System.nanoTime();
        mjCard=0;
        mnCard=0;
        for (int t=0; t<taxa; t++) {
            mjCard+=myAlign.allelePresenceForAllSites(t, WHICH_ALLELE.Major).cardinality();
            mnCard+=myAlign.allelePresenceForAllSites(t, WHICH_ALLELE.Minor).cardinality();
        }
        System.out.printf("Cardinality2 Mj:%d Mn:%d %n", mjCard, mnCard);
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("TBit Access Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        relativeRate=baseLineCPUTime/timePerObj;
        Assert.assertTrue("SBit Access Time Longer Than Expected:"+timePerObj+"ns and relative rate:"+relativeRate,relativeRate>200);
        if(testHDF5Bit) writeBitsToTestHDF5(GeneralConstants.TEMP_DIR+"raw.h5",myAlign);
        if(testHDF5Bit) readBitsToTestHDF5(GeneralConstants.TEMP_DIR+"raw.h5");
    }

    @Test
    public void testTBitConversion() {
        int taxa=5;
        BitSet mj=myAlign.allelePresenceForAllSites(taxa, WHICH_ALLELE.Major);
        BitSet mn=myAlign.allelePresenceForAllSites(taxa, WHICH_ALLELE.Minor);
        BitSet mjhomo=new OpenBitSet(mj);
        System.out.println(mjhomo);
        mjhomo.andNot(mn);
        for (int i: mjhomo.getIndicesOfSetBits()) {
             Assert.assertEquals("TBit Error in Homozygous Major",myGenotype[taxa][i],GenotypeTableUtils.getUnphasedDiploidValue(myMajor[i],myMajor[i]));
        }
        BitSet mnhomo=new OpenBitSet(mn);
        mnhomo.andNot(mj);
        for (int i: mnhomo.getIndicesOfSetBits()) {
            Assert.assertEquals("TBit Error in Homozygous Minor",myGenotype[taxa][i],GenotypeTableUtils.getUnphasedDiploidValue(myMinor[i],myMinor[i]));
        }
        BitSet het=new OpenBitSet(mj);
        het.and(mn);
        for (int i: het.getIndicesOfSetBits()) {
            Assert.assertEquals("TBit Error in Hets",myGenotype[taxa][i],GenotypeTableUtils.getUnphasedDiploidValue(myMajor[i],myMinor[i]));
        }
    }

    @Test
    public void testSBitConversion() {
        int site=5;
        BitSet mj=myAlign.allelePresenceForAllTaxa(site, WHICH_ALLELE.Major);
        BitSet mn=myAlign.allelePresenceForAllTaxa(site, WHICH_ALLELE.Minor);
        BitSet mjhomo=new OpenBitSet(mj);
        System.out.println(mjhomo);
        mjhomo.andNot(mn);
        for (int i: mjhomo.getIndicesOfSetBits()) {
            Assert.assertEquals("SBit Error in Homozygous Major",myGenotype[i][site],GenotypeTableUtils.getUnphasedDiploidValue(myMajor[site],myMajor[site]));
        }
        BitSet mnhomo=new OpenBitSet(mn);
        mnhomo.andNot(mj);
        for (int i: mnhomo.getIndicesOfSetBits()) {
            Assert.assertEquals("SBit Error in Homozygous Minor",myGenotype[i][site],GenotypeTableUtils.getUnphasedDiploidValue(myMinor[site],myMinor[site]));
        }
        BitSet het=new OpenBitSet(mj);
        het.and(mn);
        for (int i: het.getIndicesOfSetBits()) {
            Assert.assertEquals("SBit Error in Hets",myGenotype[i][site],GenotypeTableUtils.getUnphasedDiploidValue(myMajor[site],myMinor[site]));
        }
    }

    public void testBitTranspose(BitStorage mjr, BitStorage min) {
        BitSet[][] myTBitData = new BitSet[2][taxa];
        for (int t=0; t<taxa; t++) {
            myTBitData[0][t]=mjr.allelePresenceForAllSites(t);
            myTBitData[1][t]=min.allelePresenceForAllSites(t);
        }
        long time=System.nanoTime();
        BitSet[][] mySBitData = BitUtil.transpose(myTBitData, 2, taxa, sites, null);
        long mjCard=0, mnCard=0;
        for (int s=0; s<sites; s++) {
           mjCard+=mySBitData[0][s].cardinality();
           mnCard+=mySBitData[1][s].cardinality();
        }
        System.out.printf("BitTranspose Cardinality Mj:%d Mn:%d %n", mjCard, mnCard);
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("BitTranspose Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
    }

    @Ignore
    @Test
    public void testGetAllelePresenceForAllTaxa() throws Exception {
        long time=System.nanoTime();
        long mjCard=0, mnCard=0;
        for (int s=0; s<sites; s++) {
            mjCard+=myAlign.allelePresenceForAllTaxa(s, WHICH_ALLELE.Major).cardinality();
            mnCard+=myAlign.allelePresenceForAllTaxa(s, WHICH_ALLELE.Minor).cardinality();
        }
        System.out.printf("CardinalityT1 Mj:%d Mn:%d %n", mjCard, mnCard);
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("Sbit Creation and Access Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        double relativeRate=baseLineCPUTime/timePerObj;
        System.out.println("GetAllelePresenceForAllTaxa relativeRate:"+relativeRate);
        Assert.assertTrue("SBit Creation Time Longer Than Expected:"+timePerObj+"ns and relative rate:"+relativeRate,relativeRate>2);
        time=System.nanoTime();
        mjCard=0;
        mnCard=0;
        for (int s=0; s<sites; s++) {
            mjCard+=myAlign.allelePresenceForAllTaxa(s, WHICH_ALLELE.Major).cardinality();
            mnCard+=myAlign.allelePresenceForAllTaxa(s, WHICH_ALLELE.Minor).cardinality();
        }
        System.out.printf("CardinalityT2 Mj:%d Mn:%d %n", mjCard, mnCard);
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("SBit Access Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        relativeRate=baseLineCPUTime/timePerObj;
        Assert.assertTrue("SBit Access Time Longer Than Expected:"+timePerObj+"ns and relative rate:"+relativeRate,relativeRate>20);
    }

   @Ignore
    @Test
    public void pivotMatrixTo1Dim() {
        long time=System.nanoTime();
        byte[] myGenotypePivot=new byte[sites*taxa];
        int cnt=0;
        for (int t=0; t<taxa; t++) {
            for (int s=0; s<sites; s++) {
                myGenotypePivot[(t*sites)+s]=myGenotype[t][s];
            }
        }
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("Copy to 1DByte Matrix Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        time=System.nanoTime();
        byte[] myGenotypePivot2=new byte[sites*taxa];
        for (int t=0; t<taxa; t++) {
            for (int s=0; s<sites; s++) {
                myGenotypePivot[(s*taxa)+t]=myGenotype[t][s];
            }
        }
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("Pivot 1DByte Matrix Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
    }

    @Ignore
    @Test
    public void pivotMatrixTo2Dim() {
        long time=System.nanoTime();
        byte[][] myGenotypePivot=new byte[sites][taxa];
        for (int t=0; t<taxa; t++) {
            for (int s=0; s<sites; s++) {
                myGenotypePivot[s][t]=myGenotype[t][s];
            }
        }
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("Pivot Byte 2D>2D Matrix Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        time=System.nanoTime();
        for (int bigS=0; bigS<sites; bigS+=64) {
            int length=(sites-bigS<64)?sites-bigS:64;
            byte[][] myGenotypeTBlock=new byte[length][taxa];
            for (int t=0; t<taxa; t++) {
                for (int s=0; s<myGenotypeTBlock.length; s++) {
                    myGenotypeTBlock[s][t]=myGenotype[t][bigS+s];
                }
            }
        }
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("Transpose 64sites Matrix Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        //return myGenotypePivot;
    }


    @Ignore
    @Test
    public void pivotMatrixTo2Dim64() {
        long time=System.nanoTime();
        int n=10000;
        byte[] dst=new byte[n*n];
        byte[] src=new byte[n*n];
        Random r=new Random(0);
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                src[i+j*n]=(byte)r.nextInt(4);
            }
        }
        int blocksize=8;
        for (int i = 0; i < n; i += blocksize) {
            int iend=(i + blocksize>n)?n:i + blocksize;
            for (int j = 0; j < n; j += blocksize) {
                // transpose the block beginning at [i,j]
                int jend=(j + blocksize>n)?n:j + blocksize;
                for (int k = i; k < iend; ++k) {
                    for (int l = j; l < jend; ++l) {
                        dst[k + l*n] = src[l + k*n];
                    }
                }
            }
        }
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(n*n);
        System.out.println("Pivot Byte 2D>2D Matrix Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        //return myGenotypePivot;
    }

    private static void writeBytesToTestHDF5(String newHDF5file) {
        System.out.println("DynamicBitStorageTest: Writing to HDF5 Bit");
        long time=System.nanoTime();
        File hdf5File = new File(newHDF5file);
        IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
        config.overwrite();
        config.dontUseExtendableDataTypes();
        IHDF5Writer h5w = config.writer();
        HDF5IntStorageFeatures genoFeatures = HDF5IntStorageFeatures.createDeflation(defLevel);
        h5w.object().createGroup(HapMapHDF5Constants.GENOTYPES);
        for (int t = 0; t < taxa; t++) {
            String basesPath = HapMapHDF5Constants.GENOTYPES + "/T" + t;
            h5w.int8().createArray(basesPath, sites,1<<16, genoFeatures);
            HDF5Utils.writeHDF5EntireArray(basesPath, h5w, sites, 1<<16, myGenotype[t]);
        }
        h5w.close();
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("HDF5 Writing Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
    }

    private static void readBytesToTestHDF5(String newHDF5file) {
        System.out.println("DynamicBitStorageTest: Reading from HDF5 Byte");
        long time=System.nanoTime();
        IHDF5Reader reader = HDF5Factory.openForReading(newHDF5file);
        myGenotype=new byte[taxa][sites];
        for (int t = 0; t < taxa; t++) {
            String basesPath = HapMapHDF5Constants.GENOTYPES + "/T" + t;
            myGenotype[t]=reader.readAsByteArray(basesPath);
        }
        reader.close();
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("HDF5 Byte Reading Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
    }

    private static void writeBitsToTestHDF5(String newHDF5file, GenotypeTable bst) {
        System.out.println("DynamicBitStorageTest: Writing to HDF5 Bit");
        long time=System.nanoTime();
        IHDF5Writer h5w = HDF5Factory.open(newHDF5file);
        HDF5IntStorageFeatures genoFeatures = HDF5IntStorageFeatures.createDeflation(defLevel);
        h5w.object().createGroup(HapMapHDF5Constants.TBIT);
        for (int t = 0; t < taxa; t++) {
            String basesPath = HapMapHDF5Constants.TBIT + "/T" + t+"_0";
            h5w.int8().createArray(basesPath, sites,genoFeatures);
            h5w.int64().writeArray(basesPath,bst.allelePresenceForAllSites(t, WHICH_ALLELE.Major).getBits(), genoFeatures);
            basesPath = HapMapHDF5Constants.TBIT + "/T" + t+"_1";
            h5w.int8().createArray(basesPath, sites,genoFeatures);
            h5w.int64().writeArray(basesPath,bst.allelePresenceForAllSites(t, WHICH_ALLELE.Minor).getBits(), genoFeatures);
        }
        h5w.close();
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("HDF5 Bit Writing Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
    }

    private static void readBitsToTestHDF5(String newHDF5file) {
        System.out.println("DynamicBitStorageTest: Reading from HDF5 Bit");
        long time=System.nanoTime();
        IHDF5Reader reader = HDF5Factory.openForReading(newHDF5file);
        long mjCard=0, mnCard=0;
        for (int t = 0; t < taxa; t++) {
            String basesPath = HapMapHDF5Constants.TBIT + "/T" + t+"_0";
            BitSet t0=new OpenBitSet(reader.readLongArray(basesPath));
            basesPath = HapMapHDF5Constants.TBIT + "/T" + t+"_1";
            BitSet t1=new OpenBitSet(reader.readLongArray(basesPath));
            mjCard+=t0.cardinality();
            mnCard+=t1.cardinality();
        }
        System.out.printf("HDF5Bit Reading Cardinality Mj:%d Mn:%d %n", mjCard, mnCard);
        reader.close();
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(taxa*sites);
        System.out.println("HDF5Bit Reading Time:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
    }


}
