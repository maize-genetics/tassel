package net.maizegenetics.dna.snp.io;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.zip.Deflater;
import java.util.zip.Inflater;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class BuilderFromGenotypeHDF5Test {

    @Test@Ignore
    public void testMergeToMutableHDF5() throws Exception {
        String root="/Users/edbuckler/SolexaAnal/HapMapV2/";
     //   String infileRoot="/Users/edbuckler/SolexaAnal/GBS/GBS27/hmp/maizeHapMapV2_B73RefGenV2_201203028_ALLpXXXX.hmp.h5";
     //   String outfile="/Users/edbuckler/SolexaAnal/GBS/GBS27/hmp/maizeHapMapV2_B73RefGenV2_201203028_ALLt2.hmp.h5";
        int r=new Random().nextInt(1000);

        String infileRoot=root+"test/maizeHapMapV2_B73RefGenV2_201203028_chrXXXX.hmp.h5";
        String outfile=root+"test/merge"+r+".hmp.h5";;
        int startChr=8;
        int endChr=10;
        System.out.println("combineHapMapV2");
        String[] infiles=new String[endChr-startChr+1];
        for (int i = startChr; i <= endChr; i++) {
            infiles[i-startChr] = infileRoot.replace("XXXX", ""+i);
            System.out.println("Including "+infiles[i-startChr]);
        }
        System.out.println("Starting Export");
        BuilderFromGenotypeHDF5.mergeToMutableHDF5(infiles, outfile);


    }

    @Test@Ignore
    public void testDeflationBySite() throws Exception {
        String root="/Users/edbuckler/SolexaAnal/HapMapV2/";

        String infileRoot=root+"test/maizeHapMapV2_B73RefGenV2_201203028_chrXXXX.hmp.h5";
        int startChr=8;
        GenotypeTable a=GenotypeTableBuilder.getInstance(infileRoot.replace("XXXX", ""+startChr));
        int bS=1<<16;
        for (int block=0; block<10; block++) {
            long time=System.nanoTime();
            byte[] input=new byte[(1<<16)*a.numberOfTaxa()];
            int index=0;
            byte lastBase=a.genotype(0,0);
            for (int s=(bS*block); s<(bS*block)+bS; s++) {
                for (int t=0; t<a.numberOfTaxa(); t++) {
                    byte b=a.genotype(t,s);
                    if(b==GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                        input[index++]=lastBase;
                    } else {
                        input[index++]=lastBase=b;
                    }
                }
            }

            long getTime=(System.nanoTime()-time)/1000_000; time=System.nanoTime();
            // Compress the bytes
            byte[] output = new byte[input.length];
            Deflater compresser = new Deflater(1);
            compresser.setInput(input);
            compresser.finish();
            int compressedDataLength = compresser.deflate(output);
            compresser.end();
            long compTime=(System.nanoTime()-time)/1000_000; time=System.nanoTime();
//
            // Decompress the bytes
            Inflater decompresser = new Inflater();
            decompresser.setInput(output, 0, compressedDataLength);
            byte[] result = new byte[input.length];
            int resultLength = decompresser.inflate(result);
            decompresser.end();
            long deflTime=(System.nanoTime()-time)/1000_000; time=System.nanoTime();
            System.out.printf("Block:%s Orig:%d Compressed:%d Result:%d getTime:%d compTime:%d inflTime:%d %n",
                    block,output.length,compressedDataLength,resultLength,getTime,compTime,deflTime);
            Assert.assertArrayEquals(input,result);
//            for (int i=0; i<208; i++) {System.out.print(result[i]+",");}
//            System.out.println();
//            // Decode the bytes into a String
//            String outputString = new String(result, 0, resultLength, "UTF-8");
        }
//        compresser.end();

    }

    @Test@Ignore
    public void testDeflation() throws Exception {
        String root="/Users/edbuckler/SolexaAnal/HapMapV2/";
        //   String infileRoot="/Users/edbuckler/SolexaAnal/GBS/GBS27/hmp/maizeHapMapV2_B73RefGenV2_201203028_ALLpXXXX.hmp.h5";
        //   String outfile="/Users/edbuckler/SolexaAnal/GBS/GBS27/hmp/maizeHapMapV2_B73RefGenV2_201203028_ALLt2.hmp.h5";
        int r=new Random().nextInt(1000);

        String infileRoot=root+"test/maizeHapMapV2_B73RefGenV2_201203028_chrXXXX.hmp.h5";
        int startChr=8;
        GenotypeTable a=GenotypeTableBuilder.getInstance(infileRoot.replace("XXXX", ""+startChr));
  //      Deflater compresser = new Deflater(1);
        byte[] base=Arrays.copyOf(a.genotypeAllSites(0),1000);

        for (int i=0; i<a.numberOfTaxa(); i++) {
            long time=System.nanoTime();
            byte[] input = a.genotypeAllSites(i);
            long getTime=(System.nanoTime()-time)/1000_000; time=System.nanoTime();
            // Compress the bytes
            byte[] output = new byte[input.length];
            Deflater compresser = new Deflater(1);
    //        compresser.setDictionary(base);
            compresser.setInput(input);
            compresser.finish();
            int compressedDataLength = compresser.deflate(output);
            compresser.end();
            long compTime=(System.nanoTime()-time)/1000_000; time=System.nanoTime();
//
            // Decompress the bytes
            Inflater decompresser = new Inflater();
        //    decompresser.setDictionary(base);
            decompresser.setInput(output, 0, compressedDataLength);
            byte[] result = new byte[input.length];
            int resultLength = decompresser.inflate(result);
            decompresser.end();
            long deflTime=(System.nanoTime()-time)/1000_000; time=System.nanoTime();
            System.out.printf("Taxa:%s Orig:%d Compressed:%d Result:%d getTime:%d compTime:%d inflTime:%d %n",
                    a.taxaName(i),output.length,compressedDataLength,resultLength,getTime,compTime,deflTime);
//            // Decode the bytes into a String
//            String outputString = new String(result, 0, resultLength, "UTF-8");
        }
//        compresser.end();

    }


    @Test@Ignore
    public void testHMPtoHDF5() throws Exception {
       // String root="/Volumes/LaCie/HapMapV2/";
        String root="/Users/edbuckler/SolexaAnal/HapMapV2/";
        String hmpFiles=root+"AGP2/maizeHapMapV2_B73RefGenV2_201203028_chrXXXX.hmp.txt.gz";
        String h5Files=root+"test/maizeHapMapV2_B73RefGenV2_201203028_chrXXXX.hmp.h5";
        int startChr=9;
        int endChr=10;
        System.out.println("combineHapMapV2");

        for (int i = startChr; i <= endChr; i++) {
            String infile = hmpFiles.replace("XXXX", ""+i);
            String outfile = h5Files.replace("XXXX", ""+i);
            System.out.println("Reading "+infile);
            GenotypeTable a=BuilderFromHapMap.getBuilder(infile).build();
            TaxaListBuilder tlb=new TaxaListBuilder();
            Map<Taxon, Integer> taxaWithDups= new HashMap();
            for (int t = 0; t < a.taxa().size(); t++) {
                Taxon theT=a.taxa().get(t);
                if(taxaWithDups.containsKey(theT)) {
                    int cnt=taxaWithDups.get(theT);
                    cnt++;
                    Taxon newT=new Taxon.Builder(theT).name(theT.getName()+"C"+cnt).build();
                    tlb.add(newT);
                    taxaWithDups.put(theT,cnt);
                } else {
                    tlb.add(theT);
                    taxaWithDups.put(theT,1);
                }
            }
            System.out.println("Writing "+outfile);
            GenotypeTableBuilder.getInstance(a.genotypeMatrix(), a.positions(), tlb.build(), outfile);
        }

    }

    public static void main(String[] args) {
        ///Users/edbuckler/SolexaAnal/HapMapV2/test/maizeHapMapV2_B73RefGenV2_201203028_chrXXXX.hmp.h5
        //Users/edbuckler/SolexaAnal/HapMapV2/test/out1.hmp.h5
        String infileRoot=args[0];
        String outfile=args[1];
        int startChr=1;
        int endChr=10;
        System.out.println("combineHapMapV2");
        String[] infiles=new String[endChr-startChr+1];
        for (int i = startChr; i <= endChr; i++) {
            infiles[i-startChr] = infileRoot.replace("XXXX", ""+i);
            System.out.println("Including "+infiles[i-startChr]);
        }
        System.out.println("Starting Export");
        BuilderFromGenotypeHDF5.mergeToMutableHDF5(infiles, outfile);

    }
}


