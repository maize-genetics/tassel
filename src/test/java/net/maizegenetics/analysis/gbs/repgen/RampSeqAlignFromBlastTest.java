/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.Ignore;
import org.junit.Test;

import com.google.common.collect.Multimap;

import net.maizegenetics.analysis.gbs.v2.GBSUtils;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.util.DirectoryCrawler;

/**
 * Builds reference tags from BLAST hits of the {@link RepGenSimData} tags. The BLAST file is the
 * filtered three-column form the plugin expects (chromosome, subject start, subject end), with one hit
 * per tag at its amplicon and one hit reported on the reverse strand. Each reference tag should run
 * from the forward primer to the end of the reverse primer, i.e. be exactly the amplicon.
 *
 * @author lcj34
 *
 */
public class RampSeqAlignFromBlastTest {

    private static final int PERFECT_SCORE = 2 * RepGenSimData.KMER_LENGTH;
    private static final int ONE_MISMATCH_SCORE = 2 * (RepGenSimData.KMER_LENGTH - 1) - 1;

    @Test
    public void testRampSeqAlignFromBlastPlugin() throws Exception {
        RepGenSimData data = RepGenSimData.createUnder("rampSeqBlast");
        String db = RepGenLoadSeqToDBPluginTest.loadDB(data, 5);

        List<String> blastHits = new ArrayList<>();
        for (RepGenSimData.ExpectedTag tag : data.expectedTags.values()) {
            int start = tag.locus.start;
            blastHits.add(tag.locus.chromosome + "\t" + start + "\t" + (start + RepGenSimData.KMER_LENGTH - 1));
        }
        RepGenSimData.Locus first = data.loci.get(0);
        blastHits.add(first.chromosome + "\t" + (first.start + RepGenSimData.KMER_LENGTH - 1) + "\t" + first.start);
        Path blastFile = Files.write(data.dir.resolve("blast_98per_3cols.txt"), blastHits);

        new RampSeqAlignFromBlastTags()
                .dBFile(db)
                .blastFile(blastFile.toString())
                .refGenome(data.referenceFile.toString())
                .forwardp(RepGenSimData.FORWARD_PRIMER)
                .reversep(RepGenSimData.REVERSE_PRIMER)
                .performFunction(null);

        RepGenSQLite repGenData = new RepGenSQLite(db);
        try {
            Set<RefTagData> refTags = repGenData.getRefTags();
            assertEquals("One reference tag per amplicon, however many hits it has", data.loci.size(), refTags.size());
            Map<RepGenSimData.Locus, RefTagData> refTagByLocus = data.loci.stream().collect(Collectors.toMap(locus -> locus,
                    locus -> refTags.stream()
                            .filter(rtd -> rtd.chromosome().equals(locus.chromosome) && rtd.position() == locus.start)
                            .findFirst()
                            .orElseThrow(() -> new AssertionError("No reference tag at " + locus.chromosome + ":" + locus.start))));
            for (Map.Entry<RepGenSimData.Locus, RefTagData> entry : refTagByLocus.entrySet()) {
                assertEquals("Reference tag runs from the forward primer through the reverse primer",
                        entry.getKey().refAmplicon, entry.getValue().tag().sequence());
            }

            List<Tag> tags = new ArrayList<>(repGenData.getTags());
            int numTags = tags.size();
            assertEquals(data.expectedTags.size(), numTags);

            Multimap<Tag, AlignmentInfo> tagAlignments = repGenData.getAllNonRefTagAlignments(0);
            assertEquals("Each pair of tags is aligned once", numTags * (numTags - 1) / 2, tagAlignments.size());
            List<String> highScoringPairs = tagAlignments.entries().stream()
                    .filter(e -> e.getValue().score() >= ONE_MISMATCH_SCORE)
                    .map(e -> pairKey(e.getKey().sequence(), e.getValue().tag2().sequence()))
                    .sorted().collect(Collectors.toList());
            List<String> allelePairs = data.loci.stream()
                    .map(locus -> pairKey(locus.tag(false), locus.tag(true)))
                    .sorted().collect(Collectors.toList());
            assertEquals("Only the two alleles of an amplicon align with at most one mismatch", allelePairs, highScoringPairs);

            Multimap<Tag, AlignmentInfo> refAlignments = repGenData.getRefAlignmentsForTags(tags, 0);
            for (Tag tag : tags) {
                Collection<AlignmentInfo> alignments = refAlignments.get(tag);
                assertEquals("Each tag is aligned to both strands of every reference tag", 2 * refTags.size(), alignments.size());
                RepGenSimData.ExpectedTag expected = data.expectedTags.get(tag.sequence());
                AlignmentInfo best = alignments.stream().max(Comparator.comparingInt(AlignmentInfo::score)).get();
                assertEquals(refTagByLocus.get(expected.locus).tag(), best.tag2());
                assertEquals("Forward strand", 1, best.ref_strand());
                assertEquals("Alignment starts at the amplicon start", expected.locus.start, best.alignmentPos());
                assertEquals(expected.alt ? ONE_MISMATCH_SCORE : PERFECT_SCORE, best.score());
            }
        } finally {
            repGenData.close();
        }
    }

    @Test
    public void testRetrievingAlignments() throws Exception {
        RepGenSimData data = RepGenSimData.createUnder("rampSeqRetrieve");
        String db = RepGenLoadSeqToDBPluginTest.loadDB(data, 5);
        RepGenSimData.Locus locus = data.loci.get(4);
        Path blastFile = Files.write(data.dir.resolve("blast_98per_3cols.txt"),
                List.of(locus.chromosome + "\t" + locus.start + "\t" + (locus.start + RepGenSimData.KMER_LENGTH - 1)));

        new RampSeqAlignFromBlastTags()
                .dBFile(db)
                .blastFile(blastFile.toString())
                .refGenome(data.referenceFile.toString())
                .forwardp(RepGenSimData.FORWARD_PRIMER)
                .reversep(RepGenSimData.REVERSE_PRIMER)
                .performFunction(null);

        RepGenSQLite repGenData = new RepGenSQLite(db);
        try {
            List<RefTagData> refTags = new ArrayList<>(repGenData.getRefTags());
            assertEquals(1, refTags.size());
            Multimap<Tag, AlignmentInfo> aligned = repGenData.getTagAlignmentsForRefTag(refTags.get(0), ONE_MISMATCH_SCORE);
            Set<String> alignedSequences = aligned.keySet().stream().map(Tag::sequence).collect(Collectors.toSet());
            assertEquals("Both alleles of the amplicon, and nothing else, align well to its reference tag",
                    Set.of(locus.tag(false), locus.tag(true)), alignedSequences);
            for (AlignmentInfo ai : aligned.values()) {
                assertEquals(locus.start, ai.alignmentPos());
                assertEquals(1, ai.ref_strand());
            }
        } finally {
            repGenData.close();
        }
    }

    private static String pairKey(String sequence1, String sequence2) {
        return sequence1.compareTo(sequence2) < 0 ? sequence1 + "|" + sequence2 : sequence2 + "|" + sequence1;
    }

    @Ignore("#43: runs a PEAR binary from a developer's machine on the ANP68 IJE R1/R2 FASTQ files, which are not "
            + "in tassel_test_data; it only merges read pairs and asserts nothing")
    @Test
    public void testGetPEARAssemblies () {
        // This test executes the PEAR (Paired-End reAd mergeR).
        // It is assumed PEAR lives on the user's system.
        // user should adjust the forwardDir, reverseDir and outputDir's to
        // match their system.

        String forwardDir = GeneralConstants.DATA_DIR + "RepGen/anp68_ije_R1_fastq";
        String reverseDir = GeneralConstants.DATA_DIR + "RepGen/anp68_ije_R2_fastq";
        String outputDir = GeneralConstants.TEMP_DIR + "RepGen/pear_output/";

        // Get list of 5' (forward strand) files
        List<Path> forwardFiles= DirectoryCrawler.listPaths(GBSUtils.inputFileGlob, Paths.get(forwardDir).toAbsolutePath());
        if(forwardFiles.isEmpty()) {
            System.out.println("No forward files matching:"+ GBSUtils.inputFileGlob);
            return ;
        }

        if (forwardFiles.size() == 0) {
            System.out.println("testGetPEARAssemblies - found NO forward files to process.");
            System.out.println("Please verify your directory contains fastq files.");
            return ; // no files in this directory to process
        }
        Collections.sort(forwardFiles);
        System.out.println("Found " + forwardFiles.size() + " forward files to process");


        // Get list of 3' (reverse strand) fastq files
        List<Path> reverseFiles= DirectoryCrawler.listPaths(GBSUtils.inputFileGlob, Paths.get(reverseDir).toAbsolutePath());
        if(reverseFiles.isEmpty()) {
            System.out.println("No reverse files matching:"+ GBSUtils.inputFileGlob);
            return ;
        }

        if (reverseFiles.size() == 0) {
            System.out.println("testGetPEARAssemblies - found NO reverse files to process.");
            System.out.println("Please verify your reverse directory contains fastq files.");
            return ; // no files in this directory to process
        }
        System.out.println("Found " + reverseFiles.size() + "  reverse files to process");        
        Collections.sort(reverseFiles);

        // Need to be matching forward/reverse files. The assmption is the name is identical
        // on each except for R1 vs R2.  Sorting the list should result in the lists matching
        // such that traversing the lists together gives complementary files to send to PEAR.
        if (forwardFiles.size() != reverseFiles.size()) {
            System.out.println("ERROR - number of forward files " + forwardFiles.size() 
            + " does not match number of reverse files " + reverseFiles.size());          
            return;
        }

        long time=System.nanoTime();

        int fileCount = 0;
        int totalFilesProcessed = 0;
        for (int idx = 0; idx < forwardFiles.size(); idx++) {
            fileCount++;
            totalFilesProcessed++;
            Path forwardPath = forwardFiles.get(idx);
            Path reversePath = reverseFiles.get(idx);

            // Create name for output files.  the "pear_output" variable is the
            // prefix PEAR adds to each file.  For example:  If the prefix were "combined",
            // PEAR would create files called 
            //     combined.assembled.fastq
            //     combined.discarded.fastq
            //     combined.unassembled.forward.fastq
            //     combined.unassembled.reverse.fastq
            // It is the *.assembled.fastq files that will be loaded into RepGenLoadSeqToDBPlugin

            // THis code takes everything up to the last underscore.  Depending on your file,
            // you could take up to the index of the first period, or any other filtering you like.
            // This means the files look like:
            //   6970_1465_34282_ANP68_10357781_CI_7_Locus9_ACTGAGCG_CCTAGAGT.assembled.fastq
            // An alternative would be to do lastindexof ("_")
            String forwardName = forwardPath.getFileName().toString();
            String sampleName = forwardName.substring(0,forwardName.lastIndexOf("_")); // for writing the assembly

            String pear_output = outputDir + sampleName;

            String cmd="/usr/local/bin/pear  -j 8 -f " + forwardPath + " -r " + reversePath 
                    + " -o " + pear_output ;
            try {
                Process start = Runtime.getRuntime().exec(cmd);
                BufferedReader br = new BufferedReader(
                        new InputStreamReader(start.getErrorStream()));
                String line = null;
                while ((line = br.readLine()) != null)
                {
                    System.out.println(line);
                }
            } catch (Exception exc) {
                exc.printStackTrace();
            }
            if (fileCount > 19) {
                System.out.println("PEAR finished processing " + totalFilesProcessed);
                fileCount = 0;
            }
        }
        System.out.println("Finished PEAR processing " + forwardFiles.size() + " files to directory " + outputDir);

        System.out.printf("\nTotalTime for testGetPEARAssemblies: %g sec%n", (double) (System.nanoTime() - time) / 1e9);

    }
}
