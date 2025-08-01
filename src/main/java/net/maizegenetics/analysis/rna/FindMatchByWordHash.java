package net.maizegenetics.analysis.rna;

import com.google.common.collect.*;
import com.google.common.primitives.Ints;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.util.Tuple;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

import java.util.*;
import java.util.stream.IntStream;

/**
 * Find Match by Words using a simple hashmap of fixed length word (kmers) to determine most similar sequence.
 *
 * @author Ed Buckler
 * @author Karl Kremling
 *
 */
public class FindMatchByWordHash {
    public enum MatchType{COMPLETE_SEQ, FIRST_KMER, MODAL_KMER};
    private final List<Tag> tags;
    private final Map<String,int[]> kmerMap;
    private final Map<String,Integer> completeSeqMap;
    private final int wordLength;
    private final int maxWordCopies;
    private final MatchType matchType;
    private final int[] emptyInts=new int[0];
    private final boolean searchBiDirectional;

    /**
     * Create kmer
     * @param tagSet unique set sequences to be indexed with kmers
     * @param wordLength length of all the kmers
     * @param maxWordCopies retain only those kmers appearing less then this value
     */
    private FindMatchByWordHash(Set<Tag> tagSet, MatchType matchType, int wordLength, int maxWordCopies,
                                boolean searchBiDirectional) { //find a contig which is the consensus source of the kmers from a read
        this.wordLength=wordLength; // length of kmers for matching
        this.maxWordCopies=maxWordCopies; // kmer appearance frequency cutoff. Kmers are later discarded if they appear in the contigs more than maxKmerCopies times
        this.matchType=matchType;
        this.searchBiDirectional=searchBiDirectional;
        tags=new ArrayList<>(tagSet); //holds all the sequences or contigs in an array
        long bpLength=0;
        Multimap<String,Integer> completeKmerMap = HashMultimap.create(4_000_000,2); //? hashmap pointing from each kmer to its Contig of origin
        for (int ti = 0; ti < tags.size(); ti++) { //iterates through the ArrayList of tags
            String seq=tags.get(ti).sequence(); //pulls the sequence of contig/tag ti
            bpLength+=seq.length(); //increase the bpLength counter by ti's length
            for (int i = 0; i <= seq.length() - wordLength; i++) { //slides through the contig to pull out each kmer
                String sub = seq.substring(i, i + wordLength);
                completeKmerMap.put(sub, (ti+1)); //appends a key value pair to the KmerMap in which the k=kmer and the v=contig of origin
                //reverse complement reads are given negative values
                if(searchBiDirectional==true) completeKmerMap.put(BaseEncoder.getReverseComplement(sub), -(ti+1));
            }
        }
        System.out.printf("Total bp %,d Kmer Map Size: %d %n", bpLength,completeKmerMap.size());
        //purge high copy kmers from the multimap
        //consider making it immutable
        kmerMap=new HashMap<>(completeKmerMap.size()*3/2); //? why make it 1.5x larger than the complete kmer multimap if we know we will have kmers <= CompleteKmerMap.size?
        for (String kmer : completeKmerMap.keySet()) {
            Collection<Integer> tagIndices=completeKmerMap.get(kmer);
            if(tagIndices.size()<= this.maxWordCopies) { //only record key-value pairs for kmers appearing <= maxKmerCopies times
                kmerMap.put(kmer, Ints.toArray(tagIndices));
            }
        }
        //Allow direct searching for the complete sequence
        completeSeqMap = new HashMap<>(tags.size()*3/2);
        for (int ti = 0; ti < tags.size(); ti++) {
            completeSeqMap.put(tags.get(ti).sequence(),(ti+1));
            if(searchBiDirectional==true) completeKmerMap.put(BaseEncoder.getReverseComplement(tags.get(ti).sequence()), -(ti+1));
        }
        System.out.println(kmerMap.size());

    }

    /**
     * Total number of tags (reference sequences)
     */
    public int totalNumberOfTags() {
        return tags.size();
    }

    /**
     * Returns a unmodifiable list of the tags (reference sequences)
     */
    public List<Tag> getTags() {
        return Collections.unmodifiableList(tags);
    }


    /**
     * Return the first unique kmer match found starting at the 5' end of the sequence.
     * @param seq query sequence - usually the reads
     * @return matching sequence
     */
    public Match match(String seq) {
        switch(matchType) {
            case COMPLETE_SEQ:return lookForCompleteMatch(seq);
            case FIRST_KMER:return getFirstUniqueMatchIndex(seq);
            case MODAL_KMER:return getMostCommonMatchIndex(seq);
            default: return new Match();
        }
    }

    /**
     * Return the first unique kmer match found starting at the 5' end of the sequence.
     * @param seq query sequence - usually the reads
     * @return matching sequence
     */
    private Match getFirstUniqueMatchIndex(String seq) {
        OptionalInt index= IntStream.range(0,seq.length() - wordLength)
                .mapToObj(i -> kmerMap.getOrDefault(seq.substring(i, i + wordLength), emptyInts)) //search the kmer map
                .filter(hits -> hits.length == 1)  //only return kmer hits with a single match
                .mapToInt(hits -> hits[0])
                .findFirst();
        if(index.isPresent()) return new Match(index.getAsInt(),Double.NaN);
        return new Match();
    }

    /**
     * Returns the sequence match with the most kmer hits across the entire length of the query sequence.
     * @param seq query sequence - usually the reads
     * @return matching sequence
     */
    private Match getMostCommonMatchIndex(String seq) {
        int[] allHits=IntStream.range(0, seq.length() - wordLength).filter(i -> i%3==0) //stream formatting of a loop to only consider kmers starting at evey 3rd position
                .flatMap(i -> IntStream.of(kmerMap.getOrDefault(seq.substring(i, i + wordLength), emptyInts))) // load kmer from pos i through wordLength
                .toArray();
        if(allHits.length==0) return new Match(); //if there are no hits for any of the kmers in a read, return an emtpy hit entry

        return new Match(mode(allHits),Double.NaN);
    }

    private Match lookForCompleteMatch(String seq) {
        Integer tagIndex=completeSeqMap.get(seq);
        if (tagIndex!=null) {return new Match(tagIndex,Double.NaN);}
        else {return new Match();}
    }


    private static GapPenalty penalty = new SimpleGapPenalty();
    private static SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
    static {
        penalty.setOpenPenalty((short) 10);
        penalty.setExtensionPenalty((short) 2);
    }

    /**
     * Tools for evaluating the quality of alignments.  Slow but useful for evaluating quality.
     */
    private double identity(String query, String target) {
        try {
            DNASequence querySeq = new DNASequence(query);
            DNASequence refSeq = new DNASequence(target);

            SequencePair<DNASequence, NucleotideCompound> pair = Alignments.getPairwiseAlignment(refSeq, querySeq,
                    Alignments.PairwiseSequenceAlignerType.LOCAL, penalty, matrix);
            long gaps = IntStream.rangeClosed(1, pair.getLength()).filter(pair::hasGap).count();
            double identity = (double) pair.getNumIdenticals() / (double) pair.getLength();
            double snpIdentity = (double) pair.getNumIdenticals() / (double) (pair.getLength() - gaps);
            boolean debug=false;
            if (debug) System.out.println("Ref :" + refSeq.getSequenceAsString());
            if (debug) System.out.println("Query:" + querySeq.getSequenceAsString());
            if (debug) System.out.println("Alignment\n" + pair.toString());
            if (debug) System.out.println("pair.getNumIdenticals() = " + pair.getNumIdenticals());
            if (debug) System.out.println("pair.getLength() = " + pair.getLength());
            if (debug) System.out.printf("identity:%.3g snpidentity: %.3g %n", identity, snpIdentity);
            if (debug) System.out.println("-----------------------");
            return snpIdentity;
        }catch (Exception e) {
            e.printStackTrace();
            return Double.NaN;
        }
    }

    @Override
    public String toString() {
        return "FindMatchByKmers{" +
                "searchBiDirectional=" + searchBiDirectional +
                ", matchType=" + matchType +
                ", maxWordCopies=" + maxWordCopies +
                ", wordLength=" + wordLength +
                '}';
    }

    public static Tuple<Integer,Integer> calcIdentity(String querySeq, String refSeq) {
        int seedLength=6;
        //forward
        int start=refSeq.indexOf(querySeq.substring(0,seedLength));
        int matchLength=0;
        int identical=0;
        for (int i = 0; i <querySeq.length(); i++) {
            if(querySeq.charAt(i)==refSeq.charAt(i+start)) identical++;
            matchLength++;
        }
        //reverse
        return new Tuple<>(matchLength,identical);
    }

    public static int LevenshteinDistance(byte[] s, byte[] t) {
        // degenerate cases
        if (s == t) return 0;
        if (s.length == 0) return t.length;
        if (t.length == 0) return s.length;

        // create two work vectors of integer distances
        int[] v0 = new int[t.length + 1];
        int[] v1 = new int[s.length + 1];

        // initialize v0 (the previous row of distances)
        // this row is A[0][i]: edit distance for an empty s
        // the distance is just the number of characters to delete from t
        FillConsecutive(v0);

        for (int i = 0; i < s.length; i++) {
            //System.out.println(Arrays.toString(v0));
            // calculate v1 (current row distances) from the previous row v0

            // first element of v1 is A[i+1][0]
            //   edit distance is delete (i+1) chars from s to match empty t
            v1[0] = i + 1;

            // use formula to fill in the rest of the row
            for (int j = 0; j < t.length; j++)
                v1[j + 1] = Math.min(v1[j] + 1, Math.min(v0[j + 1] + 1, v0[j] + mismatchCost(s, i, t, j)));

            // copy v1 (current row) to v0 (previous row) for next iteration
            //CopyArray(v0, v1);
            System.arraycopy(v1, 0, v0, 0, v0.length);

        }

        return v1[t.length];
    }

    private static void FillConsecutive(int[] v) {
        for (int i = 0; i < v.length; i++)
            v[i] = i;
    }

    private static int mismatchCost(byte[] s, int i, byte[] t, int j) {
        return (s[i] != t[j])?1:0;
    }


    private static int mode(int[] array) {
        HashMap<Integer, Integer> cntMap = new HashMap<>();
        int maxCnt = 1, mode = array[0];
        for (int i = 0; i < array.length; i++) {
            int count = cntMap.compute(array[i], (k, cnt) -> cnt == null ? 1 : cnt + 1);
            if (count > maxCnt) {
                maxCnt = count;
                mode = array[i];
            }
        }
        return mode;
    }

    public static Builder getBuilder(Set<Tag> tagSet) {
        return new Builder(tagSet);
    }

    public static class Builder {
        private Set<Tag> tagSet;
        private FindMatchByWordHash.MatchType matchType = MatchType.MODAL_KMER;
        private int wordLength = 16;
        private int maxWordCopies = 10;
        private boolean searchBiDirectional = true;

        public Builder(Set<Tag> tagSet) {
            this.tagSet = tagSet;
        }

        public Builder matchType(FindMatchByWordHash.MatchType matchType) {
            this.matchType = matchType;
            return this;
        }

        public Builder wordLength(int wordLength) {
            this.wordLength = wordLength;
            return this;
        }

        public Builder maxWordCopies(int maxWordCopies) {
            this.maxWordCopies = maxWordCopies;
            return this;
        }

        public Builder searchBiDirectional(boolean searchBiDirectional) {
            this.searchBiDirectional = searchBiDirectional;
            return this;
        }

        public FindMatchByWordHash build() {
            return new FindMatchByWordHash(tagSet, matchType, wordLength, maxWordCopies, searchBiDirectional);
        }
    }


    public class Match {
        private final Tag tag;
        private final int tagIndex;
        private final boolean direction;
        private final double quality;

        public Match(int mapIndex, double quality) {
            this.tagIndex = Math.abs(mapIndex)-1;
            this.tag=tags.get(tagIndex);
            this.direction = mapIndex>0;
            this.quality = quality;
        }

        public Match() {
            this.tagIndex = Integer.MIN_VALUE;
            this.tag= null;
            this.direction = true;
            this.quality = Double.NaN;
        }

        public boolean isEmpty() {
            return tag==null;
        }

        public boolean isPresent() {
            return tag!=null;
        }

        public Tag tag() {
            return tag;
        }

        public String sequence() {
            return isEmpty()?"":tag.sequence();
        }

        public int tagIndex() {
            return tagIndex;
        }

        public boolean direction() {
            return direction;
        }

        public double quality() {
            return quality;
        }
    }
}
