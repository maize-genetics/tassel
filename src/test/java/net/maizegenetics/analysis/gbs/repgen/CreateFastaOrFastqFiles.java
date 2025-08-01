/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.GenomeSequence;
import net.maizegenetics.dna.map.GenomeSequenceBuilder;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;

/**
 * @author lcj34
 *
 */
public class CreateFastaOrFastqFiles {
    
    public static final String nextGenDir="/Users/lcj34/notes_files/repGen/";
    
    public static final String repGenDB_FILE=nextGenDir + "repGen.db";
    
    // public static final String referenceGenome = "/Volumes/Samsung_T1/machineLearningFiles/refGenomeFiles/Zea_mays.AGPv3.29.dna.genome.fa.gz";
     public static final String referenceGenome = "/Volumes/Samsung_T1/machineLearningFiles/refGenomeFiles/Zea_mays.AGPv4.dna.genome.fa";
    
     @Test
     public void verifyRefTag() {
         System.out.println("Begin verifyRefTag");
         GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(referenceGenome);
    
         String chromStartEndFile = "/Users/lcj34/notes_files/repGen/LDmetric_graphs/blast_files/justChromStartEndQuery_sorted2.txt";
         String indexNums = "/Users/lcj34/notes_files/repGen/LDmetric_graphs/blast_files/blast_nodups_fastaLineNums.txt";
         String tagSeqs = "/Users/lcj34/notes_files/repGen/LDmetric_graphs/blast_files/nfLynnnButAlignedByBlast.txt";
         BufferedReader rdIndex = Utils.getBufferedReader(indexNums);
         BufferedReader rdTags = Utils.getBufferedReader(tagSeqs);
         BufferedReader rdChroms = Utils.getBufferedReader(chromStartEndFile);
       //chromStarEndMap contains the index that matches indexSequeneMap, then a Tuple
         // which consists of the chrom number, then a string of <start:end>
         Map<Integer,Tuple<Integer,String>> chromStartEndMap = new HashMap<Integer,Tuple<Integer,String>>();
         // indexSequenceMap contains id number and tag sequence from fasta
         Map<Integer,String> indexSequenceMap = new HashMap<Integer,String>();
         try {
             String iLine;
             String tLine;
             while ( (iLine = rdIndex.readLine()) != null  && (tLine = rdTags.readLine()) != null) {
                 indexSequenceMap.put(Integer.parseInt(iLine), tLine);                 
             }
         } catch (Exception exc) {
             exc.printStackTrace();
             return;
         }
         
         // indexSequenceMap populated, populated the chromStartEndMap
         try {
             String chromLine;
             while ((chromLine = rdChroms.readLine()) != null) {
                 String tokens[] = chromLine.split("\t");
                 int index = Integer.parseInt(tokens[3]);
                 int chrom = Integer.parseInt(tokens[0]);
                 String startEnd = tokens[1] + ":" + tokens[2];
                 Tuple<Integer,String> chromSE = new Tuple<Integer,String>(chrom,startEnd);
                 chromStartEndMap.put(index, chromSE);
             }
         } catch (Exception exc) {
             exc.printStackTrace();
         }
         
         int notMatched = 0;
         int matched = 0;
         int numReverse = 0;
         System.out.println("LCJ - finished creating maps, on to aliging, size of indexSequenceMap keyset: " + 
            indexSequenceMap.keySet().size());
         // maps are created.  Now grab reference for each tag and see if it is a match as Blast claims
         for (int index : indexSequenceMap.keySet()) {
             String tagSeq = indexSequenceMap.get(index);
             Tuple<Integer,String> alignData = chromStartEndMap.get(index);
             Chromosome chrom = new Chromosome(Integer.toString(alignData.x));
             int colonIdx = alignData.y.indexOf(":");
             // are these 1-based values ????
             int tempstart = Integer.parseInt(alignData.y.substring(0,colonIdx));
             int tempend = Integer.parseInt(alignData.y.substring(colonIdx+1));
             int start;
             int end;
             if (tempstart > tempend) {
                 numReverse++;
                 continue; // skip reverse strand for now.
             }
//             if (tempstart < tempend) {
//                 start  = tempstart;
//                 end = tempend;
//             } else {
//                 start = tempend;
//                 end = tempstart;
//             }
             byte[] chromSeedBytes = myRefSequence.chromosomeSequence(chrom,tempstart, tempend);
             String chromKmerString = NucleotideAlignmentConstants.nucleotideBytetoString(chromSeedBytes);
             if (!tagSeq.equals(chromKmerString)) {
                 
                 notMatched++;
             } else {
                 System.out.println("\nIndex:" + index + " chrom: " + chrom.getName() + " start: "  + tempstart + " end: " + tempend + " TagSeq: " + tagSeq);                 
                 matched++;
             }
         }
         System.out.println("End of test: number tags matching ref: " + matched 
                 + ", number NOT matching: " + notMatched + ", number skipped as reverse: " + numReverse);
     }
     

    @Test
    public void createFastaFromDBTags() {

        LoggingUtils.setupDebugLogging();

        String repGenDB = nextGenDir + "repGen.db";
        //String repGenDB = "/Volumes/Samsung_T1/repgen/rampSeq_dbs/anp68ijeR2_LoadMinTaxa10.db"; // LCJ - correct this - put real db !!
        System.out.println("Begin createFastFromDBTags");
        
        RepGenDataWriter repGenData=new RepGenSQLite(repGenDB);
        Set<Tag> tagsToAlign = repGenData.getTags();
        List<Tag> tagList = new ArrayList<Tag>(tagsToAlign);
        
        String outputFile = "/Users/lcj34/notes_files/repgen/LDmetric_graphs/single_fastq.fasta";
       // String outputFile = "/Users/lcj34/notes_files/repgen/LDmetric_graphs/output_files/anp68_ijeR2Tags.fasta"; // LCJ - get a real output file!
        BufferedWriter bw = Utils.getBufferedWriter(outputFile);
        
        try {
            int count = 0;
            for (Tag tag: tagsToAlign){
                count++;
                String tagSeq = ">" + count + "\n" + tag.sequence() + "\n";
                bw.write(tagSeq);
            } 
            bw.close();
            System.out.println("FInished !!");
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }        
    }
    
    @Test
    public void createFastaFromFile() {
        LoggingUtils.setupDebugLogging();
        System.out.println("Begin createFastaFromFile");
 

        String inputFile = "/Users/lcj34/notes_files/repgen/LDmetric_graphs/output_files/Locus9R1Tags.fasta";
        String outputFile = "/Users/lcj34/notes_files/repgen/LDmetric_graphs/output_files/allTagsNotFoundR1.fasta"; // LCJ - get a real output file!
        BufferedWriter bw = Utils.getBufferedWriter(outputFile);
        BufferedReader rw = Utils.getBufferedReader(inputFile);
        
        try {
            String line;
            int count = 0;
            int linesNot150 = 0;
            while ((line = rw.readLine()) != null) {
                count++;
                if (line.length() != 150) linesNot150++;
                String tagSeq = ">" + count + "\n" + line + "\n";
                bw.write(tagSeq);
            } 
            bw.close();
            rw.close();
            System.out.println("LCJ - number of lines not equal to 150: " + linesNot150);
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }  

    }
    @Test
    public void countForwardPrimersFromTagFile() {
        LoggingUtils.setupDebugLogging();
        String repGenDB = "/Volumes/Samsung_T1/repgen/rampSeq_dbs/repGen_Locus9R1_justLoadMinTaxa2.db"; // LCJ - correct this - put real db !!
        System.out.println("Begin createFastFromDBTags");
        
        RepGenDataWriter repGenData=new RepGenSQLite(repGenDB);
        Set<Tag> tagsToAlign = repGenData.getTags();
        List<Tag> tagList = new ArrayList<Tag>(tagsToAlign);
        
        String outputFile = "/Users/lcj34/notes_files/repgen/LDmetric_graphs/output_files/tagsWithoutForwardPrimer.fasta"; // LCJ - get a real output file!
        BufferedWriter bw = Utils.getBufferedWriter(outputFile);
        
        try {
            int notForward = 0;
            for (Tag tag: tagsToAlign){
                
                String tagSeq = tag.sequence();
                if (!(tagSeq.startsWith("GCACAAGTTGTCCTGCTTCC"))) {
                    notForward++;
                    bw.write(tagSeq);
                }
            } 
            bw.close();
            System.out.println("LCJ: Nubmer of tags not starting with forward primer: " + notForward);
        } catch (IOException ioe) {
            ioe.printStackTrace();
        } 
       
    }
    
    @Test
    public void createTagFileFromFasta() {
        // I want to create a list of tag sequences, one per line
        // in a file.  These are sequences that blast was able to
        // align that we were not.  I looked at blast results after
        // filtering with awk to get down to just those fully aligned
        // (all 150 bp) with 0 mismatch and no dups.  I am pulling
        // those sequences by ID from the allTagsNotFoundR1.fasta file
        // created earlier
        LoggingUtils.setupDebugLogging();
        System.out.println("Begin createTagFileFromFasta");
 

        String tagFastaFile = "/Users/lcj34/notes_files/repgen/LDmetric_graphs/output_files/allTagsNotFoundR1.fasta";
        String indexFile = "/Users/lcj34/notes_files/repgen/LDmetric_graphs/blast_files/blast_nodups_fastaLineNums.txt";
        
        String outputFile = "/Users/lcj34/notes_files/repgen/LDmetric_graphs/blast_files/nfLynnnButAlignedByBlast.txt"; 
        BufferedWriter bw = Utils.getBufferedWriter(outputFile);
        BufferedReader rw = Utils.getBufferedReader(indexFile);
        
        try {
            String line;
            List<String> fastaList = new ArrayList<String>();
            while ((line = rw.readLine()) != null){
                String headerline = ">" + line;
                fastaList.add(headerline);
            }
            rw.close();
            rw = Utils.getBufferedReader(tagFastaFile);
            while ((line = rw.readLine()) != null) {
                String line2 = rw.readLine();
                if (fastaList.contains(line)) {
                    bw.write(line2);
                    bw.write("\n");
                }
            } 
            bw.close();
            rw.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }  

    }
}
