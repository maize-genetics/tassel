package net.maizegenetics.dna.snp;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.io.BuilderFromHapMap;
import net.maizegenetics.dna.snp.io.BuilderFromVCF;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Sizeof;
import net.maizegenetics.util.Utils;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Random;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertEquals;

/**
* Tests for Alignment Import builders.  May need to rename test, once we settle on package structure.
 */
public class ImportUtilsTest {
    public static final String outName= GeneralConstants.TEMP_DIR+"/Z001v27.hmp.t.txt.gz";
    public static final String infileNameNoZip="/Users/edbuckler/SolexaAnal/GBS/GBS27/test/Z001v27.hmp.txt";
   public static final String infileName= TutorialConstants.HAPMAP_FILENAME;
    public static final String complexVCF= GeneralConstants.DATA_DIR+"/CandidateTests/VCF4_2.vcf";
    private static final String largeHapMap=GeneralConstants.TEMP_DIR+"manySiteHapMap.hmp.txt.gz";
    private static final String largeVCF=GeneralConstants.TEMP_DIR+"manySiteVCF.vcf.gz";
    private static final String swapOrderHapMap=GeneralConstants.DATA_DIR+"CandidateTests/sitesOutOfOrder.hmp.txt.gz";


    private void createManySiteGenotypeTable() {
        if(!Files.exists(Paths.get(largeHapMap))) {
            GenotypeTable inGT=ImportUtils.readFromHapmap(infileName);
            int inflationSize=100;
            GenotypeTableBuilder outGTB=GenotypeTableBuilder.getSiteIncremental(inGT.taxa());
            for (int i=0; i<inflationSize; i++) {
                int chrBase=100*i;
                for (int s=0; s<inGT.numberOfSites(); s++) {
                    Position p=new GeneralPosition.Builder(new Chromosome(""+(chrBase+inGT.chromosome(s).getChromosomeNumber())),inGT.chromosomalPosition(s)).build();
                    outGTB.addSite(p,inGT.genotypeAllTaxa(s));
                }
            }
            ExportUtils.writeToHapmap(outGTB.build(),largeHapMap);
            ExportUtils.writeToVCF(outGTB.build(),largeVCF,false);
        }
    }

    @Test
    public void testBuilderFromHapMap() throws Exception {
        long time=System.nanoTime();
        GenotypeTable a=BuilderFromHapMap.getBuilder(infileName).build();
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(a.numberOfTaxa()*a.numberOfSites());
        System.out.printf("BuilderFromHapMap in Alignment 4.5 in %gs with site rate of %gns/bp %n", totalTime / 1e9, timePerObj);
        time=System.nanoTime();
        GenotypeTable a3=BuilderFromHapMap.getBuilder(infileName).build();
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)(a3.numberOfTaxa()*a3.numberOfSites());
        System.out.printf("BuilderFromHapMap1Again in Alignment 4.5 in %gs with site rate of %gns/bp %n", totalTime / 1e9, timePerObj);

    }


    @Test
    public void testBuilderFromHapMapV2() throws Exception {
        long time=System.nanoTime();
        GenotypeTable[] a=new GenotypeTable[10];
        for (int i = 1; i <11 ; i++) {
   //         System.out.println(Sizeof.getMemoryUse());
   //         if(i>1) Thread.sleep(1000_000);
             String fileName=infileName.replace("XXXX",""+i);

            System.out.println("Reading:"+fileName);
            a[i-1]=BuilderFromHapMap.getBuilder(fileName).build();
            System.out.println("Done Reading:"+fileName);
            System.out.println("Sites:"+a[i-1].numberOfSites());

        }
        int sites=0;
        for (GenotypeTable an: a) sites+=an.numberOfSites();
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)((long)a[0].numberOfTaxa()*(long)sites);
        System.out.println(Sizeof.getMemoryUse());
        System.out.printf("testBuilderFromHapMapV2 in Alignment 4.5 in %gs wth site rate of %gns/bp %n", totalTime / 1e9, timePerObj);

    }

    @Test
    public void testReadHapMapWithSwappedSite() throws Exception {
        System.out.println("Testing whether paired site swapping fixes minor order errors");
        System.out.println(infileName);
        GenotypeTable a=ImportUtils.readFromHapmap(infileName);
        System.out.println(swapOrderHapMap);
        GenotypeTable swapGT=ImportUtils.readFromHapmap(swapOrderHapMap);
        AlignmentTestingUtils.alignmentsEqual(a,swapGT,true);
    }


    @Ignore
    @Test
    public void testReadFromHapmap() throws Exception {
        long time=System.nanoTime();
        GenotypeTable a=ImportUtils.readFromHapmap(infileName, null);
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(a.numberOfTaxa()*a.numberOfSites());
        System.out.printf("ReadFromHapmap in Alignment 4.0 in %gs wiht site rate of %gns/bp %n", totalTime / 1e9, timePerObj);

    }

   @Ignore
    @Test
    public void testBuilderFromVCF() throws Exception {
        long time=System.nanoTime();
        String vcfIn="/Users/edbuckler/SolexaAnal/HapMap3/genotyping/HapMap2_taxa_bqc20_q5_bwamem_rlxd/thr.all_c10_bqc20_q5_bwamem_rlxd.inGBS.gz";
        GenotypeTable a=BuilderFromVCF.getBuilder(vcfIn).build();
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(a.numberOfTaxa()*a.numberOfSites());
        System.out.printf("BuilderFromVCF in Alignment 5 in %gs with site rate of %gns/bp %n", totalTime / 1e9, timePerObj);
        time=System.nanoTime();
        GenotypeTable a3=BuilderFromVCF.getBuilder(vcfIn).build();
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)(a3.numberOfTaxa()*a3.numberOfSites());
        System.out.printf("BuilderFromVCF in Alignment 5 in %gs with site rate of %gns/bp %n", totalTime / 1e9, timePerObj);

    }

    @Test
    public void testVCFWithStrangeIndels() throws Exception {
        long time=System.nanoTime();
        GenotypeTable a=BuilderFromVCF.getBuilder(complexVCF).build();
        for (Position position : a.positions()) {
            System.out.println(position.toString());
        }
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(a.numberOfTaxa()*a.numberOfSites());
        System.out.printf("BuilderFromVCF in Alignment 5 in %gs with site rate of %gns/bp %n", totalTime / 1e9, timePerObj);
//        time=System.nanoTime();
//        GenotypeTable a3=BuilderFromVCF.getBuilder(vcfIn).build();
//        totalTime=System.nanoTime()-time;
//        timePerObj=(double)totalTime/(double)(a3.numberOfTaxa()*a3.numberOfSites());
//        System.out.printf("BuilderFromVCF in Alignment 5 in %gs with site rate of %gns/bp %n", totalTime / 1e9, timePerObj);

    }

    @Test
//    @Ignore("HDF5 exception")
    public void testDirectVCFToHDF5() throws Exception {
        createManySiteGenotypeTable();
        Random r=new Random();
        String exportH5FileName=GeneralConstants.TEMP_DIR+r.nextInt(1000)+"manySiteVCFtoH5.h5";
        BuilderFromVCF.getBuilder(largeVCF).convertToHDF5(exportH5FileName).build();
        GenotypeTable largeHMP=ImportUtils.readGuessFormat(largeHapMap);
        GenotypeTable largeH5=ImportUtils.readGuessFormat(exportH5FileName);
        for (int i=0; i<largeHMP.numberOfTaxa(); i++) {
            System.out.println(i+":"+largeHMP.taxaName(i)+" ?= "+largeH5.taxaName(i));
        }
        AlignmentTestingUtils.alignmentsEqual(largeHMP,largeH5,true);

//        long time=System.nanoTime();
//        GenotypeTable a=BuilderFromVCF.getBuilder(complexVCF).build();
//        for (Position position : a.positions()) {
//            System.out.println(position.toString());
//        }
//        long totalTime=System.nanoTime()-time;
//        double timePerObj=(double)totalTime/(double)(a.numberOfTaxa()*a.numberOfSites());
//        System.out.printf("BuilderFromVCF in Alignment 5 in %gs with site rate of %gns/bp %n", totalTime / 1e9, timePerObj);
//        time=System.nanoTime();
//        GenotypeTable a3=BuilderFromVCF.getBuilder(vcfIn).build();
//        totalTime=System.nanoTime()-time;
//        timePerObj=(double)totalTime/(double)(a3.numberOfTaxa()*a3.numberOfSites());
//        System.out.printf("BuilderFromVCF in Alignment 5 in %gs with site rate of %gns/bp %n", totalTime / 1e9, timePerObj);

    }

    @Test
    public void testRoundTripVCF() throws Exception {
        long time=System.nanoTime();
        GenotypeTable gt=ImportUtils.readFromHapmap(infileName);
        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(gt.numberOfTaxa()*gt.numberOfSites());
        System.out.printf("testRoundTripVCF in Reading HapMap in %gs with site rate of %gns/bp %n", totalTime / 1e9, timePerObj);
        //annotate taxa
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (Taxon taxon : gt.taxa()) {
            Taxon.Builder tb=new Taxon.Builder(taxon.getName());
            tb.parents(taxon.getName(),taxon.getName());
            tb.inbreedF(0.97f);
            tlb.add(tb.build());
        }
        GenotypeTable nextA=GenotypeTableBuilder.getInstance(gt.genotypeMatrix(),gt.positions(),tlb.build(),AlignmentTestingUtils.createRandomDepth(gt, 1));
        String outName= GeneralConstants.TEMP_DIR+"/tempWithDepth.vcf.gz";
        time=System.nanoTime();
        ExportUtils.writeToVCF(nextA,outName,true);
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)(nextA.numberOfTaxa()*nextA.numberOfSites());
        System.out.printf("testRoundTripVCF in Writing To VCF with depth in %gs with site rate of %gns/bp %n", totalTime / 1e9, timePerObj);
        time=System.nanoTime();
        GenotypeTable nextARead=ImportUtils.readFromVCF(outName, null);
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)(nextARead.numberOfTaxa()*nextARead.numberOfSites());
        System.out.printf("testRoundTripVCF in Reading VCF with depth in %gs with site rate of %gns/bp %n", totalTime / 1e9, timePerObj);
        for (int i=0; i<nextA.numberOfTaxa(); i++) {
            assertEquals("Taxa annotation different",nextA.taxa().get(i).toStringWithVCFAnnotation(),
                    nextARead.taxa().get(i).toStringWithVCFAnnotation());
        }
//        for (Taxon taxon : nextARead.taxa()) {
//            System.out.println(taxon.toStringWithVCFAnnotation());
//        }
        AlignmentTestingUtils.alignmentsEqual(nextA,nextARead);
        AlignmentTestingUtils.alignmentsEqual(gt,nextARead);
        for (int t=0; t<gt.numberOfTaxa(); t++) {
            for (int s=0; s<gt.numberOfSites(); s++) {
                int[] dOrig=nextA.depthForAlleles(t,s);
                int[] dCopy=nextARead.depthForAlleles(t,s);
                Assert.assertArrayEquals(dOrig,dCopy);
            }
        }
    }



//    @Ignore
//    @Test
//    public void testCompareRead() throws Exception {
//        Alignment a=ImportUtils.readFromHapmap(infileName, null);
//        Alignment aN=BuilderFromHapMap.getBuilder(infileName).build();
//        for (int i = 0; i <a.numberOfTaxa(); i++) {
//            Assert.assertEquals("Taxa Names not Equal",a.taxaName(i),aN.taxaName(i));
//        }
//        for (int i = 0; i <a.numberOfSites(); i++) {
//            Assert.assertEquals("SNPID Names not Equal",a.getSNPID(i),aN.getSNPID(i));
//            Assert.assertEquals("Chr Names not Equal",a.getLocusName(i),aN.chromosomeName(i));
//            Assert.assertEquals("Position Names not Equal",a.getPositionInLocus(i),aN.chromosomalPosition(i));
//        }
//        for (int i = 0; i <a.numberOfTaxa(); i++) {
//            for (int s = 0; s <a.numberOfSites(); s++) {
//                //All the assert error are coming from hets being processed differently, I am doing
//                if(a.genotype(i, s)!=aN.genotype(i, s)) System.out.printf("t%d s%d %d %d %s %s %n", i, s,
//                        a.genotype(i, s), aN.genotype(i, s), a.genotypeAsString(i,s), aN.genotypeAsString(i,s));
//                //System.out.printf("t%d s%d %d %d %n", i, s, a.genotype(i, s), aN.genotype(i, s));
//     //           if(AlignmentUtils.isEqualOrUnknown(a.genotype(i, s), aN.genotype(i, s))) continue;
//                Assert.assertTrue("Base not Equal at taxa" + i + " site" + s, AlignmentUtils.isEqual(a.genotype(i, s), aN.genotype(i, s)));
//            }
//
//        }
//
//    }

    @Ignore
    @Test
    public void speedFileTests() {
        long time=System.nanoTime();
        try {
            final byte nlb=(byte) '\n';
            final char nlc='\n';
            BufferedReader r=Utils.getBufferedReader(infileName, -1);
            //  BufferedReader r=Utils.getBufferedReader(infileNameNoZip);
            String s;
            ArrayList<String> stuff=new ArrayList<>();
            int lines=0;
            while((s=r.readLine())!=null) {
                stuff.add(s);
                lines++;
            }
            r.close();
            if(false) {
                OutputStream os=new GZIPOutputStream(new FileOutputStream(new File(outName))) {
                    {
                        def.setLevel(2);
                    }
                };
                BufferedWriter w=new BufferedWriter(new OutputStreamWriter(os));
                for (String sx: stuff) {
                    w.write(sx);
                }
                w.close(); }
            long totalTime=System.nanoTime()-time;
            System.out.printf("ImportUtil ReadText data timing %gs %d %n", totalTime/1e9, lines);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Ignore
    @Test
    public void speedFileTests2() {
        long time=System.nanoTime();
        try {
            InputStream gzipOS=new GZIPInputStream(new FileInputStream(infileName));

            //gzipOS=new FileInputStream(infileNameNoZip);
            ReadableByteChannel in = Channels.newChannel(gzipOS);
            ByteBuffer buf = ByteBuffer.allocate(100_000_000);
            int bytesRead = in.read(buf);
            int lines=0;
            //time=System.nanoTime();
            while (bytesRead != -1) {
                System.out.println("Read "+bytesRead);
                buf.flip();
                byte[] b=new byte[100_000_001];
                int i=0;
                buf.get(b,0,buf.limit());
                buf.clear();
                bytesRead = in.read(buf);
                lines++;
            }

            long totalTime=System.nanoTime()-time;
            System.out.printf("ImportUtil ReadText data timing %gs %d %n", totalTime/1e9, lines);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
