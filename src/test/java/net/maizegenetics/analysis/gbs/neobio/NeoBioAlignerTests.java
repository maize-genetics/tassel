/**
 * 
 */
package net.maizegenetics.analysis.gbs.neobio;

import static org.junit.Assert.*;

import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;

import org.junit.Test;

/**
 * @author lcj34
 *
 */
public class NeoBioAlignerTests {
    
    @Test
    public void neoBioIdenticalSequenceTest() {
        // I added row and column offset as something
        // that can be returned from neobio

        String seq1 = "CGGGTGTGACAGTCGTGCAGTCGACCGTTGGG"; // 32 length
        String seq2 = "CGGGTGTGACAGTCGTGCAGTCGACCGTTGGG";  
        Reader reader1 = new StringReader(seq1);
        Reader reader2 = new StringReader(seq2);

        PairwiseAlignmentAlgorithm  algorithm;
        ScoringScheme  scoring;
        PairwiseAlignment alignment;

        algorithm = new SmithWaterman();

        scoring = new BasicScoringScheme(2,-2,-1); // match_reward, mismatch_penalty, gap_cost
        algorithm.setScoringScheme(scoring);
        try {
            algorithm.loadSequences(reader1, reader2);

            // get score, perfect score is 64 (2 points for each match, 32 bp length
            int score = algorithm.getScore();
            System.out.println("SCORE: " + score);
            //Compute alignment
            alignment = algorithm.getPairwiseAlignment();
            System.out.println(seq1);
            System.out.println(seq2);
            System.out.println("\n\n");
            System.out.println("ALignment: \n" + alignment);

 
            int seq1Offset = alignment.getRowStart();
            int seq2Offset = alignment.getColStart();

            assertEquals(seq1Offset,0); // sequences were equal, alignment should be from beginning
            assertEquals(seq2Offset,0);
            assertEquals(score,64);
            
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (InvalidSequenceException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IncompatibleScoringSchemeException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
    
    @Test
    public void neoBioSeq2LongerTest() {
        // Sequences are identical, but sequence 2 is longer than sequence 1
        String seq1 = "CGGGTGTGACAGTCGTGCAGTCGACCGTTGGG";
        String seq2 = "XXXXCGGGTGTGACAGTCGTGCAGTCGACCGTTGGGXXXX";  
        Reader reader1 = new StringReader(seq1);
        Reader reader2 = new StringReader(seq2);

        PairwiseAlignmentAlgorithm  algorithm;
        ScoringScheme  scoring;
        PairwiseAlignment alignment;

        algorithm = new SmithWaterman();

        scoring = new BasicScoringScheme(2,-2,-1);// match_reward, mismatch_penalty, gap_cost
        algorithm.setScoringScheme(scoring);
        try {
            algorithm.loadSequences(reader1, reader2);

            // get score
            int score = algorithm.getScore();
            System.out.println("SCORE: " + score);
            //Compute alignment
            alignment = algorithm.getPairwiseAlignment();
            System.out.println(seq1);
            System.out.println(seq2);
            System.out.println("\n\n");
            System.out.println("ALignment: \n" + alignment);

            // first sequence given to loadSequences is stored in row, second is in columns
            int seq1Offset = alignment.getRowStart();
            int seq2Offset = alignment.getColStart();
 
            System.out.println("seq1offset - row: " + seq1Offset + ", seq2offset = col: " + seq2Offset);
            assertEquals(seq1Offset,0); // Sequence is aligned from the beginning
            assertEquals(seq2Offset,4); // sequence aligns past the initial XXXX
            assertEquals(score,64); // it matched perfectly from where it started
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (InvalidSequenceException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IncompatibleScoringSchemeException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
    
    @Test
    public void neoBioSeq2LongerBothMisMatchTest() {
        // Sequences are identical, but sequence 2 is longer than sequence 1
        String seq1 = "AACTGGTGTGACAGTCGTGCAGTCGACCGTTGGGTA";
        String seq2 = "XXXXCGGGTGTGACAGTCGTGCAGTCGACCGTTGGGXXXX";  
        Reader reader1 = new StringReader(seq1);
        Reader reader2 = new StringReader(seq2);

        PairwiseAlignmentAlgorithm  algorithm;
        ScoringScheme  scoring;
        PairwiseAlignment alignment;

        algorithm = new SmithWaterman();

        scoring = new BasicScoringScheme(2,-1,-1);// match_reward, mismatch_penalty, gap_cost
        algorithm.setScoringScheme(scoring);
        try {
            algorithm.loadSequences(reader1, reader2);

            // get score
            int score = algorithm.getScore();
            System.out.println("SCORE: " + score);
            //Compute alignment
            alignment = algorithm.getPairwiseAlignment();
            System.out.println(seq1);
            System.out.println(seq2);
            System.out.println("\n\n");
            System.out.println("ALignment: \n" + alignment);

            // first sequence given to loadSequences is stored in row, second is in columns
            int seq1Offset = alignment.getRowStart();
            int seq2Offset = alignment.getColStart();
 
            System.out.println("seq1offset - row: " + seq1Offset + ", seq2offset = col: " + seq2Offset);
            assertEquals(seq1Offset,2); // Sequence is aligned from 3rd char
            assertEquals(seq2Offset,4); // sequence aligns past the initial XXXX
            assertEquals(score,61); // one mismatch
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (InvalidSequenceException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IncompatibleScoringSchemeException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
 
    @Test
    public void neoBioSeqHasMiddleGapTest() {
        // Sequence1 has gap in middle, ends of each are different
        String seq1 = "AACTGGTGTGACAGTCGTGCAGTCGACCGTTGGGTA";
        String seq2 = "XXXXCTGGTGTGATATATATATATACAGTCGTGCAGTCGACCGTTGGGXXXX";  
        Reader reader1 = new StringReader(seq1);
        Reader reader2 = new StringReader(seq2);

        PairwiseAlignmentAlgorithm  algorithm;
        ScoringScheme  scoring;
        PairwiseAlignment alignment;

        algorithm = new SmithWaterman();

        // This does best with either mismatch higher, or the 2 penalties equal
        scoring = new BasicScoringScheme(2,-1,-1);// match_reward, mismatch_penalty, gap_cost
        algorithm.setScoringScheme(scoring);
        try {
            algorithm.loadSequences(reader1, reader2);

            // get score
            int score = algorithm.getScore();
            System.out.println("SCORE: " + score);
            //Compute alignment
            alignment = algorithm.getPairwiseAlignment();
            System.out.println(seq1);
            System.out.println(seq2);
            System.out.println("\n\n");
            System.out.println("ALignment: \n" + alignment);

            // first sequence given to loadSequences is stored in row, second is in columns
            int seq1Offset = alignment.getRowStart();
            int seq2Offset = alignment.getColStart();
 
            System.out.println("seq1offset - row: " + seq1Offset + ", seq2offset = col: " + seq2Offset);
            assertEquals(seq1Offset,2); // Sequence is aligned from 3rd char
            assertEquals(seq2Offset,4); // sequence aligns past the initial XXXX
            //assertEquals(score,61); // one mismatch
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (InvalidSequenceException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IncompatibleScoringSchemeException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
    @Test
    public void neoBioSequence1LongerTest() {

        String seq1 = "XXXXCGGGTGTGACAGTCGTGCAGTCGACCGTTGGGXXXX";
        String seq2 = "CGGGTGTGACAGTCGTGCAGTCGACCGTTGGG";
        
        Reader reader1 = new StringReader(seq1);
        Reader reader2 = new StringReader(seq2);

        PairwiseAlignmentAlgorithm  algorithm;
        ScoringScheme  scoring;
        PairwiseAlignment alignment;

        algorithm = new SmithWaterman();

        scoring = new BasicScoringScheme(2,-2,-1);
        algorithm.setScoringScheme(scoring);
        try {
            algorithm.loadSequences(reader1, reader2);

            // get score
            int score = algorithm.getScore();
            System.out.println("SCORE: " + score);
            //Compute alignment
            alignment = algorithm.getPairwiseAlignment();
            System.out.println(seq1);
            System.out.println(seq2);
            System.out.println("\n\n");
            System.out.println("ALignment: \n" + alignment);

            // first sequence given to loadSequences is stored in row, second is in columns
            int seq1Offset = alignment.getRowStart();
            int seq2Offset = alignment.getColStart();
 
            assertEquals(seq1Offset,4); // Sequence is aligned past the initial XXXX
            assertEquals(seq2Offset,0); // sequence aligns from the beginning
            assertEquals(score,64); // it matched perfectly from where it started (32 matches = score of 64)

        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (InvalidSequenceException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IncompatibleScoringSchemeException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
    
    @Test
    public void neoBio2BeginningMisMatchesTest() {
        // run example as per NeoBio's textui package,
        // file Neobio.java
        // THIs is nice, but I need to know the starting position
        // of where seq1 is against seq 2.  Can we change the code
        // to print this out?

        // Looks like the start position of both is just r and c,
        // which I printed in line in PariwiseAlignmentAlgorigtm.
        // Need to have these 2 values returned, perhaps stored in
        // PariwriseAlignment?

        String seq1 = "CGGGTGTGACAGTCGTGCAGTCGACCGTTGGG";
        String seq2 = "ACGGTGTGACAGTCGTGCAGTCGACCGTTGGG";
        
        Reader reader1 = new StringReader(seq1);
        Reader reader2 = new StringReader(seq2);

        PairwiseAlignmentAlgorithm  algorithm;
        ScoringScheme  scoring;
        PairwiseAlignment alignment;

        algorithm = new SmithWaterman();

        scoring = new BasicScoringScheme(2,-2,-1);
        algorithm.setScoringScheme(scoring);
        try {
            algorithm.loadSequences(reader1, reader2);

            // get score
            int score = algorithm.getScore();
            System.out.println("SCORE: " + score);
            //Compute alignment
            alignment = algorithm.getPairwiseAlignment();
            System.out.println(seq1);
            System.out.println(seq2);
            System.out.println("\n\n");
            System.out.println("ALignment: \n" + alignment);

            // first sequence given to loadSequences is stored in row, second is in columns
            int seq1Offset = alignment.getRowStart();
            int seq2Offset = alignment.getColStart();
 
            assertEquals(seq1Offset,0); // Sequence is aligned from the beginning
            assertEquals(seq2Offset,1); // sequence aligns past the initial A
            // 31 matches = score of 62, 1 gap = 62=1 = 61
            assertEquals(score,61); // it matched perfectly from where it started

        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (InvalidSequenceException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IncompatibleScoringSchemeException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
    
    @Test
    public void neoBio4MisMatchesInMiddleTest() {

        String seq1 = "CGGGTGTGGCAGTCGTGCAGTCGACCGTTGGG";
        String seq2 = "CCGGTGTGACAGCCGTGCAGTCGAACGTTGGG";
        
        Reader reader1 = new StringReader(seq1);
        Reader reader2 = new StringReader(seq2);

        PairwiseAlignmentAlgorithm  algorithm;
        ScoringScheme  scoring;
        PairwiseAlignment alignment;

        algorithm = new SmithWaterman();

        scoring = new BasicScoringScheme(2,-2,-1);
        algorithm.setScoringScheme(scoring);
        try {
            algorithm.loadSequences(reader1, reader2);

            // get score
            int score = algorithm.getScore();
            System.out.println("SCORE: " + score);
            //Compute alignment
            alignment = algorithm.getPairwiseAlignment();
            System.out.println(seq1);
            System.out.println(seq2);
            System.out.println("\n\n");
            System.out.println("ALignment: \n" + alignment);

            // first sequence given to loadSequences is stored in row, second is in columns
            int seq1Offset = alignment.getRowStart();
            int seq2Offset = alignment.getColStart();
 
            // THis gives gaps vs mismatches, starting at column 1 in sequence 2
            // ALignment: 
            //    CGGGTGTGG-CAGTC-GTGCAGTCGAC-CGTTGGG
            //    | ||||| | ||| | ||||||||||  |||||||
            //    C-GGTGT-GACAG-CCGTGCAGTCGA-ACGTTGGG
            //    Score: 49
            assertEquals(seq1Offset,0); // Sequence is aligned from the beginning
            assertEquals(seq2Offset,1); // sequence aligns past the initial A
            // 28 matches = score of 56, 7 gaps = 56-7 = 49
            assertEquals(score,49); // it matched perfectly from where it started

        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (InvalidSequenceException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IncompatibleScoringSchemeException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
    
    @Test
    public void neoBio3MisMatchesInMiddleTest() {

        // This test is similar to below, but gaps are given a higher
        // penalty than mismatches.
        String seq1 = "CGGGTGTGGCAGTCGTGCAGTCGACCGTTGGG";
        String seq2 = "CGGGTGTGCCAGCCGTGCAGTCGAACGTTGGG";
        
        Reader reader1 = new StringReader(seq1);
        Reader reader2 = new StringReader(seq2);

        PairwiseAlignmentAlgorithm  algorithm;
        ScoringScheme  scoring;
        PairwiseAlignment alignment;

        algorithm = new SmithWaterman();

        scoring = new BasicScoringScheme(2,-1,-2);
        algorithm.setScoringScheme(scoring);
        try {
            algorithm.loadSequences(reader1, reader2);

            // get score
            int score = algorithm.getScore();
            System.out.println("SCORE: " + score);
            //Compute alignment
            alignment = algorithm.getPairwiseAlignment();
            System.out.println(seq1);
            System.out.println(seq2);
            System.out.println("\n\n");
            System.out.println("ALignment: \n" + alignment);

            // first sequence given to loadSequences is stored in row, second is in columns
            int seq1Offset = alignment.getRowStart();
            int seq2Offset = alignment.getColStart();
 
            assertEquals(seq1Offset,0); // Sequence is aligned from the beginning
            assertEquals(seq2Offset,0); // sequence aligns past the initial A
            // 29 matches = score of 58, 3 mismatches,  = 58-3 = 55
            assertEquals(score,55); // it matched perfectly from where it started

        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (InvalidSequenceException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IncompatibleScoringSchemeException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
    
    @Test
    public void neoBio3GapsInMiddleTest() {

        // This test is similar to above, but mismatches are given a higher
        // penalty than mismatches.
        String seq1 = "CGGGTGTGGCAGTCGTGCAGTCGACCGTTGGG";
        String seq2 = "CGGGTGTGCCAGCCGTGCAGTCGAACGTTGGG";
        
        Reader reader1 = new StringReader(seq1);
        Reader reader2 = new StringReader(seq2);

        PairwiseAlignmentAlgorithm  algorithm;
        ScoringScheme  scoring;
        PairwiseAlignment alignment;

        algorithm = new SmithWaterman();

        scoring = new BasicScoringScheme(2,-2,-1);
        algorithm.setScoringScheme(scoring);
        try {
            algorithm.loadSequences(reader1, reader2);

            // get score
            int score = algorithm.getScore();
            System.out.println("SCORE: " + score);
            //Compute alignment
            alignment = algorithm.getPairwiseAlignment();
            System.out.println(seq1);
            System.out.println(seq2);
            System.out.println("\n\n");
            System.out.println("ALignment: \n" + alignment);

            // first sequence given to loadSequences is stored in row, second is in columns
            int seq1Offset = alignment.getRowStart();
            int seq2Offset = alignment.getColStart();
 
            assertEquals(seq1Offset,0); // Sequence is aligned from the beginning
            assertEquals(seq2Offset,0); // sequence aligns past the initial A
            // 29 matches = score of 58, 6 gaps,  = 58-6 = 52
            assertEquals(score,52); // it matched perfectly from where it started

        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (InvalidSequenceException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IncompatibleScoringSchemeException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
}
