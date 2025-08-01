package net.maizegenetics.dna.tag;

import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.constants.GeneralConstants;
import org.junit.Ignore;
import org.junit.Test;
import org.xerial.snappy.Snappy;
import org.xerial.snappy.SnappyOutputStream;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.zip.GZIPOutputStream;

public class TagsByTaxaStorageTest {
    private static final String gbs4TBTTag= GBSConstants.GBS_EXPECTED_MODIFY_TBT_HDF5_PLUGIN_PIVOTED_FILE;
    private static final String outName= GeneralConstants.TEMP_DIR+"RANDTaxaTest.h5";
    private static String outFile;

    static {
        Random r=new Random();
        outFile=outName.replace("RAND",""+r.nextInt(1000));
    }


    private static HashMap<Tag,TaxaDistribution> createBigTBT(int taxaExpansion, int tagExp) throws Exception {
        System.out.println("gbs4TBTTaxa = " + gbs4TBTTag);
        TagsByTaxaByteHDF5TagGroups aTBT4=new TagsByTaxaByteHDF5TagGroups(gbs4TBTTag);
        System.out.println("aTBT4.getTaxaCount() = " + aTBT4.getTaxaCount());
        System.out.println("aTBT4.getTagCount() = " + aTBT4.getTagCount());
        System.out.printf("Original File size:%d%n", Files.size(Paths.get(gbs4TBTTag)));
        HashMap<Tag,TaxaDistribution> bh=new HashMap<>(aTBT4.tagCount*tagExp);
        for (int tagExp_i = 0; tagExp_i < tagExp; tagExp_i++) {
            for (int tgi = 0; tgi < aTBT4.getTagCount(); tgi++) {
                long[] seqs = aTBT4.getTag(tgi);
                seqs[0]+=tagExp_i;
                Tag tag = TagBuilder.instance(seqs, (short) aTBT4.getTagLength(tgi)).build();
                TaxaDistribution td = TaxaDistBuilder.create(aTBT4.getTaxaCount() * taxaExpansion);
               // if(tgi<10) System.out.printf("Initial tag %d: %s%n",tgi,Arrays.toString(aTBT4.getTaxaReadCountsForTag(tgi)));
                for (int i = 0; i < taxaExpansion; i++) {
                    for (int txi = 0; txi < aTBT4.getTaxaCount(); txi++) {
                        int numTagsForTaxon = aTBT4.getReadCountForTagTaxon(tgi, txi);
                        for (int k = 0; k < numTagsForTaxon; k++) {
                            td.increment(txi + (i*aTBT4.getTaxaCount()));
                        }
                    }
                }
                //if(tgi<10) System.out.printf("Initial tag %d: %s%n",tgi,Arrays.toString(td.encodeTaxaDepth()));
                bh.put(tag, td);
            }
            System.out.println("tagExp_i = " + tagExp_i);
        }
       return bh;
    }

    private long[] getTaxaInfoCounts(HashMap<Tag,TaxaDistribution> map) {
        int infoCnt=0, totalReads=0;
        for (Map.Entry<Tag, TaxaDistribution> entry : map.entrySet()) {
            int[][] taxaDist=entry.getValue().taxaWithDepths();
            infoCnt+=taxaDist[0].length;
            for (int i : taxaDist[1]) {
                totalReads+=i;
            }

        }
        int taxa=map.entrySet().iterator().next().getValue().maxTaxa();
        return new long[]{map.size(),taxa,infoCnt,totalReads};
    }


    @Test
    public void testSerialization() throws Exception {
        System.out.println("TagsByTaxaHDF5Test.testSerialization");

        int taxaExpansion=5;
        int tagExp=1;
        System.out.printf("Expand number of taxa by %d, number of tags by %d%n", taxaExpansion, tagExp);
        //create map
        HashMap<Tag,TaxaDistribution> bh=createBigTBT(taxaExpansion,tagExp);
        long[] mapCounts=getTaxaInfoCounts(bh);


        long time = System.nanoTime();

        outFile=outFile.replace("h5","xser.gz");
        try (
                OutputStream file = new GZIPOutputStream(new FileOutputStream(outFile));
                //OutputStream file = new SnappyOutputStream(new FileOutputStream(outFile));
                OutputStream buffer = new BufferedOutputStream(file);
                ObjectOutput output = new ObjectOutputStream(buffer);
        ){
            output.writeObject(bh);
        }
        catch(IOException ex){
            ex.printStackTrace();
        }

        System.out.println("time: " + (System.nanoTime() - time)/1e9);
        System.out.printf("New File size:%d%n", Files.size(Paths.get(outFile)));
        System.out.printf("New Total bytes per Taxa-Tag combination:%g%n", (double)Files.size(Paths.get(outFile))/(double)mapCounts[2]);
        System.out.println("Total Tags:"+mapCounts[0]);
        System.out.println("Total Taxa:"+mapCounts[1]);
        System.out.println("Total Taxa-Tag combinations:"+mapCounts[2]);
        System.out.println("Total Reads:"+mapCounts[3]);
        System.out.println("Total Reads/Taxon:"+mapCounts[3]/(mapCounts[1]));


        time=0;

    }

    @Test
    public void testSerializationEncode() throws Exception {
        System.out.println("TagsByTaxaHDF5Test.testSerializationEncode");

        int taxaExpansion=5;
        int tagExp=1;
        System.out.printf("Expand number of taxa by %d, number of tags by %d%n", taxaExpansion, tagExp);
        //create map
        HashMap<Tag,TaxaDistribution> bh=createBigTBT(taxaExpansion,tagExp);
        long[] mapCounts=getTaxaInfoCounts(bh);


        long time = System.nanoTime();
        Map<Tag,byte[]> map = new HashMap<>();

        for (Map.Entry<Tag, TaxaDistribution> entry : bh.entrySet()) {
            map.put(entry.getKey(), Snappy.compress(entry.getValue().encodeTaxaDepth()));
        }

        outFile=outFile.replace("h5","xser.gz");
        try (

                OutputStream file = new GZIPOutputStream(new FileOutputStream(outFile));
                //OutputStream file = new SnappyOutputStream(new FileOutputStream(outFile));
                OutputStream buffer = new BufferedOutputStream(file);
                ObjectOutput output = new ObjectOutputStream(buffer);
        ){
            output.writeObject(map);
        }
        catch(IOException ex){
            ex.printStackTrace();
        }

        System.out.println("time: " + (System.nanoTime() - time)/1e9);
        System.out.printf("New File size:%d%n", Files.size(Paths.get(outFile)));
        System.out.printf("New Total bytes per Taxa-Tag combination:%g%n", (double)Files.size(Paths.get(outFile))/(double)mapCounts[2]);
        System.out.println("Total Tags:"+mapCounts[0]);
        System.out.println("Total Taxa:"+mapCounts[1]);
        System.out.println("Total Taxa-Tag combinations:"+mapCounts[2]);
        System.out.println("Total Reads:"+mapCounts[3]);
        System.out.println("Total Reads/Taxon:"+mapCounts[3]/(mapCounts[1]));


        time=0;

    }


}