package net.maizegenetics.ed.t5rna;

import net.maizegenetics.analysis.rna.FindMatchByWordHash;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.util.Tuple;
import org.junit.Before;
import org.junit.Test;

import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.*;

/**
 * Created by edbuckler on 10/21/15.
 */
public class FindMatchByWordHashTest {
    private FindMatchByWordHash findMatchByWordHash;

    @Before
    public void setUp() throws Exception {
        Set<Tag> tagSet=new HashSet<>();
        tagSet.add(TagBuilder.instance("ACGTACGT").build());
        tagSet.add(TagBuilder.instance("AACCGGTT").build());
        tagSet.add(TagBuilder.instance("AAACCCGGG").build());
        tagSet.add(TagBuilder.instance("TTTTTTTTT").build());
        findMatchByWordHash = FindMatchByWordHash.getBuilder(tagSet)
                .matchType(FindMatchByWordHash.MatchType.MODAL_KMER)
                .wordLength(4)
                .searchBiDirectional(false)
                .build();
        System.out.println(findMatchByWordHash.toString());
    }

    @Test
    public void testFindFirstMatch() throws Exception {
        assertEquals("ACGTACGT", findMatchByWordHash.match("ACGTACGT").sequence());
        assertEquals("AACCGGTT", findMatchByWordHash.match("AACCGGTT").sequence());
        assertEquals("AAACCCGGG", findMatchByWordHash.match("AAACCGGG").sequence());
        assertEquals("TTTTTTTTT", findMatchByWordHash.match("TTTTTTTTT").sequence());

        //test subtle mismatches
        assertEquals("ACGTACGT", findMatchByWordHash.match("ACGTTCGT").sequence());
        assertEquals("AACCGGTT", findMatchByWordHash.match("AGCCGGTT").sequence());
        assertEquals("AAACCCGGG", findMatchByWordHash.match("AAACCGCG").sequence());
        assertEquals("TTTTTTTTT", findMatchByWordHash.match("TGTTTTTTT").sequence());

        assertEquals("AACCGGTT", findMatchByWordHash.match("CGGAGCCGGTTTGT").sequence());

    }

    @Test
    public void testSingleMismatchAliger() {
        String refSeq=   "ACGTAGCTACGATGCATCTAGCTGACTACT";
        String queryPerfSeq=   "CTACGATGCATCTAGCTGACT";
        String queryMisSeq=    "CTACGATTCATCTAGCTGACT";
        String queryIndelSeq=  "CTACGAGCATCTAGCTGACT";

        Tuple<Integer,Integer> matchLenAndIdentityCnt=FindMatchByWordHash.calcIdentity(queryPerfSeq,refSeq);
        System.out.println(matchLenAndIdentityCnt.toString());
        System.out.println(FindMatchByWordHash.LevenshteinDistance(refSeq.getBytes(), queryPerfSeq.getBytes()));
        matchLenAndIdentityCnt=FindMatchByWordHash.calcIdentity(queryMisSeq,refSeq);
        System.out.println(matchLenAndIdentityCnt.toString());
        System.out.println(FindMatchByWordHash.LevenshteinDistance(refSeq.getBytes(), queryMisSeq.getBytes()));
        matchLenAndIdentityCnt=FindMatchByWordHash.calcIdentity(queryIndelSeq,refSeq);
        System.out.println(matchLenAndIdentityCnt.toString());
    }

}