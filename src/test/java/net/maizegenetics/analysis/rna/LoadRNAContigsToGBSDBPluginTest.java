package net.maizegenetics.analysis.rna;

import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.tag.*;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportUtils;
import net.maizegenetics.util.Tuple;

import org.junit.Test;

import java.io.File;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.junit.Assert.*;

/**
 * Created by edbuckler on 10/19/15.
 */
public class LoadRNAContigsToGBSDBPluginTest {
    private static String myDBName="/Users/edbuckler/Development/temp/rnaDB.db";

    @Test
    public void testPlugin() throws Exception {
        LoggingUtils.setupDebugLogging();
        String myDBName= new LoadRNAContigsToGBSDBPlugin()
                .deleteOldData(true)
                .contigFile("/Users/edbuckler/Development/PrivMaizeGen/src/net/maizegenetics/karl/B73_GShoot_contigscopy.fasta.gz")
                .outputDB("/Users/edbuckler/Development/temp/rnaDB.db")
                .runPlugin(null);
        System.out.println(myDBName);

        TagDataWriter tdw =new TagDataSQLite(myDBName);
        Set<Tag> tags=tdw.getTags();
        assertEquals(36923, tags.size());


        tdw.getTagsNameMap().entrySet().stream().forEach(entry -> {
            System.out.println(entry.getKey().toString() + "=" + entry.getValue());
        });

        ((TagDataSQLite)tdw).close();
    }

    @Test
    public void testKmerMatch() throws Exception {
        //10 billion reads
        TagDataWriter tdw =new TagDataSQLite(myDBName);
        Set<Tag> allContigs=tdw.getTags();
        FindMatchByWordHash matcher= FindMatchByWordHash.getBuilder(allContigs)
                .matchType(FindMatchByWordHash.MatchType.MODAL_KMER)
                .maxWordCopies(100)
                .wordLength(16)
                .build();

        LongAdder cntHits=new LongAdder();
        LongAdder cntMisses=new LongAdder();
        long time=System.nanoTime();
        for (int i = 0; i < 3; i++) {
            tdw.getTags().stream().forEach(t -> {
                        if (matcher.match(t.sequence()).isPresent()) {
                            cntHits.increment();
                        } else {
                            cntMisses.increment();
                          //  if (t.seqLength() > 30) System.out.println(t.toString());
                        }
                    }
            );
        }
        ((TagDataSQLite)tdw).close();
        long totalTime=System.nanoTime()-time;
        double timePerRead=(double)totalTime/(double)(cntMisses.longValue()+cntHits.longValue());
        System.out.println("cntMisses = " + cntMisses.longValue());
        System.out.println("cntHits = " + cntHits.longValue());
        System.out.println("timePerRead = " + timePerRead);
    }

    @Test
    public void testKmerMatchAccuracy() throws Exception {
        //10 billion reads
        TagDataWriter tdw = new TagDataSQLite(myDBName);
        Set<Tag> allContigs = tdw.getTags();
        FindMatchByWordHash matcher= FindMatchByWordHash.getBuilder(allContigs)
                .matchType(FindMatchByWordHash.MatchType.MODAL_KMER)
                .wordLength(16)
                .maxWordCopies(100)
                .build();

        LongAdder cntHits = new LongAdder();
        LongAdder cntIncorrectHits = new LongAdder();
        LongAdder cntMisses = new LongAdder();
        final Random r=new Random(0);

        long time = System.nanoTime();
        tdw.getTags().stream().forEach(t -> {
                    String querySequence = t.sequence();
                    String mutatedSequence = mutateString(querySequence, r);
                    if (t.seqLength() < 30) return;
                    if (maxBaseProportion(querySequence) > 0.50) return;
                    FindMatchByWordHash.Match matchTag = matcher.match(querySequence);
//                    FindMatchByWordHash.Match matchTag = matcher.match(mutatedSequence);
                    if (matchTag.isEmpty()) {
                        cntMisses.increment();
                        if (t.seqLength() > 30) System.out.println(t.toString());
                    } else {
                        if (matchTag.sequence().equals(querySequence)) {
                            cntHits.increment();
                        } else {
                            System.out.println("querySequence = " + querySequence + " !=" + matchTag.sequence());
                            cntIncorrectHits.increment();
                        }
                    }
                }
        );
        ((TagDataSQLite)tdw).close();
        long totalTime = System.nanoTime() - time;
        double timePerRead = (double) totalTime / (double) allContigs.size();
        System.out.println("cntMisses = " + cntMisses.longValue());
        System.out.println("cntIncorrectHits = " + cntIncorrectHits.longValue());
        System.out.println("cntHits = " + cntHits.longValue());
        System.out.println("timePerRead = " + timePerRead);
    }

    @Test
    public void testRNADeMultiplexProductionPlugin() throws Exception {

        LoggingUtils.setupDebugLogging();

        new RNADeMultiplexProductionPlugin()
                .batchSize(1)
                .keyFile("/Users/edbuckler/Development/temp/RNAGBS_tempdir/Pipeline_Testing_Simplified_Key_RNA.txt")
                .inputDir("/Users/edbuckler/Development/temp/RNAGBS_tempdir/")
                .inputDB("/Users/edbuckler/Development/temp/rnaDB.db")
                .batchSize(4)
                .wordSizeMatching(16)
                .maxWordRepeats(100)
                .matchingType(FindMatchByWordHash.MatchType.MODAL_KMER)
                .runPlugin(null);
        System.out.println(myDBName);

        TagData tdw =new TagDataSQLite(myDBName);
        TableReport p0 = tdw.getAllCountsForTaxonTissue(new Taxon("B73ind1"), "Gshoot");
        if (p0 != null)
           Stream.of(p0.toStringTabDelim().split("\n")).limit(10).forEach(System.out::println);
        
        TaxaList myTaxonList = tdw.getTaxaList();
        Set<String> myTissues = tdw.getAllTissue();
        int taxonNum = 0;
        // Print values to tab delimited file.  Shu-Yun has just 5 taxon, 2 tissues,
        // so not too many files.  Others may have many more taxa/tissue pairs.
        // Code is commented out by default to avoid a mess of
        // files being created.  It is  here as an example of how to get a tab-delimited
        // file from the TableReport, and could be run selectively for specific taxon/tissue
        // pairs.
//        for (Taxon taxon : myTaxonList) {
//            for (String tissue: myTissues) {
//                TableReport report=tdw.getAllCountsForTaxonTissue(taxon, tissue);
//                System.out.println("Table name: " + report.getTableTitle() + " rowCount: " + report.getRowCount());
//                if (report == null || report.getRowCount() == 0) {
//                    System.out.println("No data for Tissue " + tissue + ", taxon " + taxonNum);
//                }
//                
//                // Choose your own output file directory
//                String tissueTaxonFile = GBSConstants.GBS_TEMP_DIR + tissue + "_" + taxonNum + ".txt";
//                try {
//                    File theFile = new File(tissueTaxonFile);
//                    TableReportUtils.saveDelimitedTableReport(report, "\t", theFile);
//                } catch (Exception e) {
//                    ((TagDataSQLite)tdw).close();
//                    e.printStackTrace();
//                    throw new IllegalStateException("ExportPlugin: performFunctionForTableReport: Problem writing file: " + tissueTaxonFile);
//                }
//            }
//        }

        ((TagDataSQLite)tdw).close();
        System.out.println();
        Phenotype p = tdw.getAllCountsForTagTissue(TagBuilder.instance("CAGAAGATTGAAATTTTCTCACTCAAAAAAAAAAAAAAAAAA").build(),"Gshoot");

        if (p != null)
           System.out.println(p.toStringTabDelim());

    }

    @Test
    public void testExpressionData() throws Exception {
        LoggingUtils.setupDebugLogging();

        TagData tdw =new TagDataSQLite(myDBName);
        TableReport p0 = tdw.getAllCountsForTaxonTissue(new Taxon("B73ind1"), "Gshoot");
        Stream.of(p0.toStringTabDelim().split("\n")).limit(10).forEach(System.out::println);


        System.out.println();
        Phenotype p=tdw.getAllCountsForTagTissue(TagBuilder.instance("CAGAAGATTGAAATTTTCTCACTCAAAAAAAAAAAAAAAAAA").build(),"Gshoot");

        System.out.println(p.toStringTabDelim());

        Phenotype p3=tdw.getAllCountsForTissue("Gshoot");
        System.out.println(p3.toStringTabDelim());
        ((TagDataSQLite)tdw).close();

    }

    private String mutateString(String inputString, Random r) {
        final char[] bases={'A','C','G','T'};
        StringBuilder sb = new StringBuilder(inputString);
        sb.setCharAt(r.nextInt(sb.length()),bases[r.nextInt(4)]);
        return sb.toString();
    }

    private double maxBaseProportion(String sequence) {
        long baseA=sequence.chars().filter(c -> c=='A').count();
        long baseT=sequence.chars().filter(c -> c=='T').count();
        return (double)Math.max(baseA,baseT)/(double)sequence.length();
    }

}