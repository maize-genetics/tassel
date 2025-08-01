package net.maizegenetics.dna.tag;


import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import junit.framework.Assert;
import net.maizegenetics.dna.map.*;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.dna.snp.SimpleAllele;
import net.maizegenetics.dna.tag.*;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Tuple;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;
import org.junit.runners.Suite;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

public class TagDataWriterTest {
    static TagDataWriter tdw;
//    private static final String dbName= GeneralConstants.TEMP_DIR+"RANDSQLList.db";
//    private static final String gbs4TBTTag= GBSConstants.GBS_EXPECTED_MODIFY_TBT_HDF5_PLUGIN_PIVOTED_FILE;
    private static final String dbName= "tempDir/GBS/"+"RANDSQLList.db";

    public static final String RAW_SEQ_CHR_9_10_20000000 = "Chr9_10-20000000/";
    public static final String RAW_SEQ_CURRENT_TEST = RAW_SEQ_CHR_9_10_20000000;
    public static final String GBS_TEMP_DIR = "tempDir/GBS/" + RAW_SEQ_CURRENT_TEST;
    public static final String GBS_EXPECTED_DIR = "dataFiles/ExpectedResults/GBS/" + RAW_SEQ_CURRENT_TEST;

    private static final String gbs4TBTTag= GBS_EXPECTED_DIR + "ModifyTBTHDF5Plugin/TBT_Pivoted.h5";
    private static final String gbs4TOPM= GBS_EXPECTED_DIR + "ModifyTBTHDF5Plugin/TBT_Pivoted.h5";

    public static final String GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_DIR = GBS_EXPECTED_DIR + "DiscoverySNPCallerPlugin/";
    public static final String GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE = GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_DIR + "TOPM_with_Variants.topm";


    private static Map<Tag,TaxaDistribution> inputTDMap;
    private static Multimap<Tag,Position> tagCutMap;
    private static Multimap<Tag,Allele> tagAlleleMap;
    private static TaxaList inputTaxaList;


    @BeforeClass
    public static void  setUpBeforeClass() throws Exception {
        Files.deleteIfExists(Paths.get(dbName));
        tdw=new TagDataSQLite(dbName);
        int taxaExpansion=1;
        int tagExpansion=1;
        inputTaxaList=createTaxaList(taxaExpansion);
        inputTDMap=createBigTBT(taxaExpansion,tagExpansion);
        tagCutMap=createBigTOPM(tagExpansion);
        tagAlleleMap=createBigTOPMWithAllele(tagExpansion);
    }

    @Ignore
    @Test
    public void testPutTag() throws Exception {

        //Map<Tag,TaxaDistribution> inputTDMap=createBigTBT(1,1);
        long time=System.nanoTime();
        System.out.println("Map size:"+inputTDMap.size());
        tdw.putAllTag(inputTDMap.keySet());
        double timePerTag=(System.nanoTime()-time)/inputTDMap.size();
        System.out.println("insert tags:  timePerTag = " + timePerTag +" ns/tag");
        time=System.nanoTime();
        int count=0;
        for (Tag t : tdw.getTags()) {
            count++;
        }
        timePerTag=(System.nanoTime()-time)/inputTDMap.size();
        System.out.println("reads tags:  timePerTag = " + timePerTag +" ns/tag");
        time=System.nanoTime();
        System.out.println("Map size:"+count);
        tdw.putTaxaList(inputTaxaList);
        tdw.putTaxaDistribution(inputTDMap);
        timePerTag=(System.nanoTime()-time)/inputTDMap.size();
        System.out.println("add taxa distributions:  timePerTag = " + timePerTag +" ns/tag");

        time=System.nanoTime();
        int count2=0;
        for (Map.Entry<Tag, TaxaDistribution> entry : inputTDMap.entrySet()) {
            TaxaDistribution dbTD=null;
            try{
                dbTD=tdw.getTaxaDistribution(entry.getKey());
            } catch (Exception e) {
                TaxaDistribution tdOrig=entry.getValue();
                System.out.println(tdOrig.toString());
                byte[] enc=tdOrig.encodeTaxaDepth();
                System.out.println(Arrays.toString(enc));
                TaxaDistribution tdNew=TaxaDistBuilder.create(enc);
                System.out.println(tdNew.toString());
            }
            count2++;
            if(count2%10000==0) System.out.println(count2);
            Assert.assertTrue("Taxa Distribution not equal",entry.getValue().equals(dbTD));
        }
        timePerTag=(System.nanoTime()-time)/inputTDMap.size();
        System.out.println("access taxa distributions:  timePerTag = " + timePerTag +" ns/tag");
        time=System.nanoTime();

//        for (Tag tag : inputTDMap.keySet()) {
//            tdw.putTag(tag);
//        }
//        tdw.getTags();
    }

    @Ignore
    @Test
    public void testCutPositions() throws Exception {
        System.out.println("Adding tags from TOPM");
        tdw.putAllTag(tagCutMap.keySet());
        tdw.putTagAlignments(tagCutMap);
        HashSet<Position> origPos=tagCutMap.values().stream().collect(Collectors.toCollection(HashSet::new));
        PositionList dbPositionsAll=tdw.getTagCutPositions(true);
        for (Position origPosition : dbPositionsAll) {
            Assert.assertTrue(origPos.contains(origPosition));
        }
        PositionList dbPositionsSub=tdw.getTagCutPositions(new Chromosome("9"),500_000,750_000,true);
        dbPositionsSub.stream().forEach(p -> Assert.assertTrue(
                p.getChromosome().getName().equals("9")  && p.getPosition()>500_000 && p.getPosition()<750_000));
        Map<Position, Map<Tag,Tuple<Boolean,TaxaDistribution>>> posMap=tdw.getCutPositionTagTaxaMap(new Chromosome("9"),500_000,750_000);
        //posMap.entrySet().stream().forEach(e -> System.out.println(e.getKey().toString()+"="+e.getValue().size()));
    }

    @Ignore
    @Test
    public void testAlleles() throws Exception {

        System.out.println("Adding alleles from TOPM");
        System.out.println("tagAlleleMap.keySet().size() = " + tagAlleleMap.keySet().size());
        System.out.println("tagAlleleMap.size = " + tagAlleleMap.size());
        tdw.putTagAlleles(tagAlleleMap);

        Multimap<Tag,Allele> tagAlleleMapDB=tdw.getAlleleMap();
        Assert.assertEquals(tagAlleleMap.size(),tagAlleleMapDB.size());
        for (Map.Entry<Tag, Allele> tagAlleleEntry : tagAlleleMap.entries()) {
            Assert.assertTrue(tagAlleleMapDB.containsEntry(tagAlleleEntry.getKey(),tagAlleleEntry.getValue()));
        }

    }

    private static Multimap<Tag,Position> createBigTOPM(int tagExp) {
        TOPMInterface topm= new TagsOnPhysicalMap(GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE, true);
        System.out.println("Initial TOPM tags numbers ="+topm.getTagCount());
        System.out.println("Tag Expansion = " + tagExp);
        Random r=new Random(0);
        Multimap<Tag,Position> bh=HashMultimap.create(topm.getTagCount()*tagExp,1);
        for (int tagExp_i = 0; tagExp_i < tagExp; tagExp_i++) {
            for (int i = 0; i < topm.getTagCount(); i++) {
                long[] seqs = topm.getTag(i);
                if(tagExp_i>0) seqs[0] += r.nextLong();
                Tag tag = TagBuilder.instance(seqs, (short) topm.getTagLength(i)).build();
                Chromosome chr=topm.getLocus(i);
                if(chr==null) chr=Chromosome.UNKNOWN;
                Position position = new GeneralPosition
                        .Builder(chr, topm.getStartPosition(i))
                        .addAnno("mappingapproach", "BWA")
                        .addAnno("cigar", "64M")
                        .addAnno("supportvalue", topm.getMapP(i))
                        //.addAnno("best",topm)
                        .build();
                bh.put(tag,position);
            }
        }
        return bh;
    }

    private static Multimap<Tag,Allele> createBigTOPMWithAllele(int tagExp) {
        TOPMInterface topm= new TagsOnPhysicalMap(GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE, true);
        System.out.println("Initial TOPM tags numbers ="+topm.getTagCount());
        System.out.println("Tag Expansion = " + tagExp);
        Random r=new Random(0);
        Multimap<Tag,Allele> bh= HashMultimap.create(topm.getTagCount() * tagExp,3);
        int multiVariantsTag=0;
        for (int tagExp_i = 0; tagExp_i < tagExp; tagExp_i++) {
            for (int i = 0; i < topm.getTagCount(); i++) {
                long[] seqs = topm.getTag(i);
                if(tagExp_i>0) seqs[0] += r.nextLong();
                Tag tag = TagBuilder.instance(seqs, (short) topm.getTagLength(i)).build();
                Chromosome chr=topm.getLocus(i);
                if(chr==null) chr=Chromosome.UNKNOWN;
                byte[] varDef=topm.getVariantDefArray(i);
                byte[] varOffset=topm.getVariantDefArray(i);
                if(varDef==null) continue;
                if(varDef.length>1) multiVariantsTag++;
                for (int alleleNum = 0; alleleNum < varDef.length; alleleNum++) {
                    Position position = new GeneralPosition
                            .Builder(chr, topm.getStartPosition(i)+varOffset[alleleNum])
                            .addAnno("mappingapproach", "BWA")
                            .addAnno("cigar", "64M")
                            .addAnno("supportvalue", topm.getMapP(i))
                                    //.addAnno("best",topm)
                            .build();
                    Allele allele=new SimpleAllele(varDef[alleleNum],position);
                    bh.put(tag,allele);
                }

            }
            System.out.println("multiVariantsTag = " + multiVariantsTag);
        }
        return bh;
    }



    private static Map<Tag,TaxaDistribution> createBigTBT(int taxaExpansion, int tagExp) throws Exception {
        System.out.println("gbs4TBTTaxa = " + gbs4TBTTag);
        TagsByTaxaByteHDF5TagGroups aTBT4=new TagsByTaxaByteHDF5TagGroups(gbs4TBTTag);
        System.out.println("aTBT4.getTaxaCount() = " + aTBT4.getTaxaCount());
        System.out.println("aTBT4.getTagCount() = " + aTBT4.getTagCount());
        System.out.printf("Original File size:%d%n", Files.size(Paths.get(gbs4TBTTag)));
        Map<Tag,TaxaDistribution> bh=new HashMap<>(aTBT4.getTagCount()*tagExp);
        Random r=new Random(0);
        long origReadCnt=0,readCnt=0;
        for (int tagExp_i = 0; tagExp_i < tagExp; tagExp_i++) {
            for (int tgi = 0; tgi < aTBT4.getTagCount(); tgi++) {
                long[] seqs = aTBT4.getTag(tgi);
                if(tagExp_i>0) seqs[0]+=r.nextLong();

                Tag tag = TagBuilder.instance(seqs, (short) aTBT4.getTagLength(tgi)).build();
                TaxaDistribution td = TaxaDistBuilder.create(aTBT4.getTaxaCount() * taxaExpansion);
                // if(tgi<10) System.out.printf("Initial tag %d: %s%n",tgi,Arrays.toString(aTBT4.getTaxaReadCountsForTag(tgi)));
                for (int i = 0; i < taxaExpansion; i++) {
                    for (int txi = 0; txi < aTBT4.getTaxaCount(); txi++) {
                        int numTagsForTaxon = aTBT4.getReadCountForTagTaxon(tgi, txi);
                        if(tagExp_i==0 && i==0) origReadCnt+=numTagsForTaxon;
                        int taxon=txi + (i*aTBT4.getTaxaCount());
                        if(i>0) taxon-=r.nextInt(50);
                        if(numTagsForTaxon>1) {
                            numTagsForTaxon+=(r.nextInt(3)-1);
                        }
                        for (int k = 0; k < numTagsForTaxon; k++) {
                            td.increment(taxon);
                            readCnt++;
                        }
                    }
                }
                //if(tgi<10) System.out.printf("Initial tag %d: %s%n",tgi,Arrays.toString(td.encodeTaxaDepth()));
                bh.put(tag, td);
            }
            System.out.println("taxaExpansion = [" + taxaExpansion + "], tagExp = [" + tagExp + "]");
            System.out.println("origReadCnt = " + origReadCnt);
            System.out.println("readCnt = " + readCnt);
        }
        return bh;
    }

    private static TaxaList createTaxaList(int taxaExpansion) throws Exception {
        System.out.println("gbs4TBTTaxa = " + gbs4TBTTag);
        TagsByTaxaByteHDF5TagGroups aTBT4 = new TagsByTaxaByteHDF5TagGroups(gbs4TBTTag);
        final AtomicInteger counter=new AtomicInteger();
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int txExp = 0; txExp < taxaExpansion; txExp++) {
            aTBT4.getTaxaList().stream().forEachOrdered(t -> tlb.add(new Taxon(t.getName()+counter.incrementAndGet())));
        }
        return tlb.build();
    }

}