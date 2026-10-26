package net.maizegenetics.dna.snp.io;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.zip.Deflater;
import java.util.zip.Inflater;

/**
 * @author Ed Buckler
 */
public class BuilderFromGenotypeHDF5Test {

    //Not in tassel_test_data; see docs/developer/missing-test-data.md
    private static final String HAPMAP_V2_CHR8_H5 = GeneralConstants.DATA_DIR + "CandidateTests/maizeHapMapV2_B73RefGenV2_201203028_chr8.hmp.h5";

    @Test
    public void testMergeToMutableHDF5() throws Exception {
        //chrom3 lists the same taxa as the others in a different order, so the merge has to remap taxa
        String[] hapmaps = {"mdp_genotype_chrom3_taxa_diff_order.hmp.txt", "mdp_genotype_chrom8.hmp.txt", "mdp_genotype_chrom10.hmp.txt"};
        GenotypeTable[] sources = new GenotypeTable[hapmaps.length];
        String[] infiles = new String[hapmaps.length];
        int expectedSites = 0;
        for (int i = 0; i < hapmaps.length; i++) {
            sources[i] = BuilderFromHapMap.getBuilder(TutorialConstants.TUTORIAL_DIR + hapmaps[i]).build();
            infiles[i] = GeneralConstants.TEMP_DIR + "mergeToMutableSource" + i + ".h5";
            Files.deleteIfExists(Paths.get(infiles[i]));
            GenotypeTableBuilder.getInstance(sources[i], infiles[i]);
            expectedSites += sources[i].numberOfSites();
        }
        String outfile = GeneralConstants.TEMP_DIR + "mergeToMutable.h5";
        Files.deleteIfExists(Paths.get(outfile));

        BuilderFromGenotypeHDF5.mergeToMutableHDF5(infiles, outfile);

        GenotypeTable merged = GenotypeTableBuilder.getInstance(outfile);
        Assert.assertEquals("Merged site count", expectedSites, merged.numberOfSites());
        Assert.assertEquals("Merged taxa count", sources[0].numberOfTaxa(), merged.numberOfTaxa());
        for (GenotypeTable source : sources) {
            int[] mergedTaxon = new int[source.numberOfTaxa()];
            for (int t = 0; t < source.numberOfTaxa(); t++) {
                mergedTaxon[t] = merged.taxa().indexOf(source.taxa().get(t));
                Assert.assertTrue("Taxon missing from merge: " + source.taxaName(t), mergedTaxon[t] >= 0);
            }
            for (int s = 0; s < source.numberOfSites(); s++) {
                int mergedSite = merged.positions().indexOf(source.positions().get(s));
                Assert.assertTrue("Site missing from merge: " + source.positions().get(s), mergedSite >= 0);
                for (int t = 0; t < source.numberOfTaxa(); t++) {
                    Assert.assertEquals("Genotype differs for " + source.taxaName(t) + " at " + source.positions().get(s),
                            source.genotype(t, s), merged.genotype(mergedTaxon[t], mergedSite));
                }
            }
        }
    }

    @Ignore("benchmark: needs the maize HapMap v2 chr8 HDF5 file, which is not in tassel_test_data (#36)")
    @Test
    public void testDeflationBySite() throws Exception {
        GenotypeTable a=GenotypeTableBuilder.getInstance(HAPMAP_V2_CHR8_H5);
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

    @Ignore("benchmark: needs the maize HapMap v2 chr8 HDF5 file, which is not in tassel_test_data (#36)")
    @Test
    public void testDeflation() throws Exception {
        GenotypeTable a=GenotypeTableBuilder.getInstance(HAPMAP_V2_CHR8_H5);
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
}
