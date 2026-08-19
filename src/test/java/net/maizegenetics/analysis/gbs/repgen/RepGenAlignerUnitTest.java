/**
 *
 */
package net.maizegenetics.analysis.gbs.repgen;

import static org.junit.Assert.*;

import java.lang.reflect.Method;
import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;

/**
 * Pure-unit coverage for {@link RepGenAlignerPlugin} extracted from {@code RepGenAlignerPluginTest}
 * (which is excluded because its other tests require hardcoded dev-machine paths and external data).
 * This class contains only the self-contained kmer-seed generation test and therefore runs in CI
 * without any external data.
 *
 * @author lcj34 (original), extracted for CI unit coverage
 */
public class RepGenAlignerUnitTest {

    @Test
    public void testCreateKmerSeedsFromDBTags() {
        // Create a list of tags, run them through the plugin's private kmer-seed builder, and verify
        // we get exactly the kmer set (forward + reverse complement) the reference implementation
        // below produces.
        String seq1 = "ACGTCGATCTAGGGGGTCTCGACGAAGGCAACCATCTTGGTATGCCGGAGGACGGTGATCCTCCTAGGCCTGCGCCTCGCGTTGACATCCTTCGGGAGCTAGCTGTGGTCCCAGTCCCTGCGAGGGGTCAGGACGCATAGCTCGAGCAAA";
        String seq2 = "ACGTCGATCTAGGGGGTCTCGACGAAGGAAACCATTTTGGTATGCCGGAGGACGGTGATCCCCCTAGGCCCGCGCCTCGCGTTGACATCCTTTGGGAGCTAGCTATGGTCCCAGTCCCTGTGGGGGGTCAGGACGACGCACAGCTCGAGC";
        String seq3 = "GCACAAGTTGTCCTGCTTCCTCGTCGAGCCTGGCGTGCATCTCGCGGATTTTCTCGAGCTGTGCGTCCTGACCCCCCGCAGGGACTGGGACCACAGCTAGCTCCCGAAGGATGTCAACGCGAGGCGCAGGCCTAGAGGGATCGCCGTCCT";
        String seq4 = "GCACAAGTTGTCCTGCTTCCTCGTCGAGCCTGGCCTGTACCTCGCGGATTTGCTCGAGCTGTGCGTCCTGACCTTCCGCAGGGACTGGGACCACAGCTAGCTCCCGAAGGATGTCAACGCGAGGGGCAGACCTAGGGGGATCACCGTCCT";
        String seq5 = "GCACAAGTTGTCCTGCTTCCTCGTCGAGCTTGGTCTGCATCTCACGGATTTGCTCGAGCTGTGTGTCCTGACCCCCCGCAGGGACTGGGACCACAGCTAGCTCCCGAAGGATGTCAACACGAGGTGCAGGCCTAGGGGGATCGTCGTCCT";

        Set<String> sequences = new HashSet<String>();
        sequences.add(seq1);
        sequences.add(seq2);
        sequences.add(seq3);
        sequences.add(seq4);
        sequences.add(seq5);

        Multimap<String, String> kmerSequenceMap = HashMultimap.create();
        int window = 20;
        int seedlen = 17;
        for (String seq : sequences) {
            for (int idx = 0; idx < seq.length() - window;) {
                String kmer = seq.substring(idx, idx + seedlen);
                kmerSequenceMap.put(kmer, seq);
                byte[] kmerRC = NucleotideAlignmentConstants.reverseComplementAlleleByteArray(kmer.getBytes());
                String kmerRCString = new String(kmerRC);
                if (kmerRCString.contains("N")) {
                    idx += window;
                    continue;
                }
                kmerSequenceMap.put(kmerRCString, seq);
                idx += window;
            }
        }

        Set<Tag> tags = new HashSet<Tag>();
        tags.add(TagBuilder.instance(seq1).build());
        tags.add(TagBuilder.instance(seq2).build());
        tags.add(TagBuilder.instance(seq3).build());
        tags.add(TagBuilder.instance(seq4).build());
        tags.add(TagBuilder.instance(seq5).build());

        Multimap<String, Tag> kmerTagMap = HashMultimap.create();
        RepGenAlignerPlugin repGenAlignerPlugin = new RepGenAlignerPlugin();
        Class theClass = repGenAlignerPlugin.getClass();
        try {
            Method theMethod = theClass.getDeclaredMethod("createKmerSeedsFromDBTags",
                    new Class[]{Set.class, Multimap.class, int.class});
            theMethod.setAccessible(true);
            theMethod.invoke(repGenAlignerPlugin, tags, kmerTagMap, 20);
        } catch (Exception exc) {
            exc.printStackTrace();
            fail("createKmerSeedsFromDBTags threw: " + exc.getMessage());
        }

        assertEquals(kmerSequenceMap.keySet().size(), kmerTagMap.keySet().size());

        for (String kmer : kmerSequenceMap.keySet()) {
            if (!(kmerTagMap.keySet().contains(kmer))) {
                fail("createKmerSeedsFromDBTags does not contain kmer from junit: " + kmer);
            }
        }
        for (String kmer : kmerTagMap.keySet()) {
            if (!(kmerSequenceMap.keySet().contains(kmer))) {
                fail("createKmerSeedsFromDBTags contains kmer not in junit array: " + kmer);
            }
        }
    }
}
