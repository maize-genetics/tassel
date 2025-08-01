/**
 * 
 */
package net.maizegenetics.dna.map;

import java.util.*;

import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.Tuple;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * @author Lynn Johnson
 *
 */
public class GenomeSequenceBuilderTest {

    public static final String REFERENCE_GENOME_EXPECTED_DIR = "dataFiles/ReferenceGenomeSequenceTests/";
    public static final String REFERENCE_GENOME_FILE = REFERENCE_GENOME_EXPECTED_DIR + "Zea_maysParsed.AGPv3.23.dna.genome.fa";
    public static final String REFERENCE_GENOME_KEYFILE = REFERENCE_GENOME_EXPECTED_DIR + "chromosomes.txt";

    public static final String REFERENCE_CHROM910_FILE = GBSConstants.GBS_DATA_DIR + "ZmB73_RefGen_v2_chr9_10_1st20MB.fasta";
    public static final String TEST_SHORT_CHROM3_FILE = REFERENCE_GENOME_EXPECTED_DIR + "testChrom3_short.fa";
    public static final String TEST_SHORT_CHROM3456_FILE = REFERENCE_GENOME_EXPECTED_DIR + "testChrom3456_short.fa";
    public static final String TEST_UNKNOWNA_FILE = REFERENCE_GENOME_EXPECTED_DIR + "testUnknownA.fa";

    String chromeAndPosition = "";

    /**
     * @throws java.lang.Exception
     */
    @BeforeClass
    public static void  setUpBeforeClass() throws Exception {

    }

    /**
     * Test method for {@link net.maizegenetics.dna.map.HalfByteGenomeSequence#chromosomes()}.
     */
    @Test
    public void testChromosomes() {
        // Tests that the requested chromosomes were put into the Chromosomes array	
        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(REFERENCE_CHROM910_FILE);

        Chromosome nine = new Chromosome("9");
        Chromosome ten = new Chromosome("10");

        byte[] refGenomeBytes = myRefSequence.chromosomeSequence(nine);
        assert(refGenomeBytes.length > 0);       
        refGenomeBytes = myRefSequence.chromosomeSequence(ten);
        assert(refGenomeBytes.length > 0);
        
        Set<Chromosome> chromosomeSet  = myRefSequence.chromosomes();
        assertTrue(chromosomeSet.contains(nine));
        assertTrue(chromosomeSet.contains(ten));       
    }

    /**
     * Test method for {@link net.maizegenetics.dna.map.HalfByteGenomeSequence#chromosomeSequence(net.maizegenetics.dna.map.Chromosome)}.
     */
    @Test
    public void testChromosomeSequenceChromosome() {
        // Verify we can read a chromosome and that when we pull the chromosome bytes
        // they match the bytes obtained from reading the file.
        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(TEST_SHORT_CHROM3_FILE);        
        Chromosome three = new Chromosome("3");
        String chrom3 = "AAGCTTGTGAAGGTTCTTCATCCCCACATGTGCTAAGCGGCGATGCCACAGCCAGCCCAT";
        byte[] chrom3bytes = myRefSequence.chromosomeSequence(three);

        assertEquals(chrom3.length(), chrom3bytes.length);    
    }

    /**
     * Test method for {@link net.maizegenetics.dna.map.HalfByteGenomeSequence#chromosomeSequence(net.maizegenetics.dna.map.Chromosome)}.
     */
    @Test
    public void testInstanceWithStringSequence() {
        // Verify we can read a chromosome and that when we pull the chromosome bytes
        // they match the bytes obtained from reading the file.
        String chromosomeSeq="AAAACCCCGGGGTTTT";
        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(new Chromosome("10"), chromosomeSeq);
        assertEquals("Length problem",chromosomeSeq.length(),myRefSequence.genomeSize());

        assertEquals("Sequence problem",chromosomeSeq,myRefSequence.genomeSequenceAsString(0,myRefSequence.genomeSize()-1));
        assertEquals("Correct result ",chromosomeSeq,NucleotideAlignmentConstants.nucleotideBytetoString(
                myRefSequence.chromosomeSequence(new Chromosome("10"),1,(int)myRefSequence.genomeSize())));
        // chromoseomSequence() expects 1-based parameters.  If start is < 1, it will return null.
        try {
            assertNull(myRefSequence.chromosomeSequence(new Chromosome("10"),0,(int)myRefSequence.genomeSize()));
            fail("myRefSequence.chromosomeSequence passed, but should have failed for 0 as starting parameter");
        } catch (IllegalArgumentException iae) {
            System.out.println("Test of 0-based chromosomeSequence passes.");
        }
        
//        assertEquals("Error in Zero based response",chromosomeSeq,NucleotideAlignmentConstants.nucleotideBytetoString(
//                myRefSequence.chromosomeSequence(new Chromosome("10"),0,(int)myRefSequence.genomeSize())));
//        System.out.println(Arrays.toString(myRefSequence.chromosomeSequence(new Chromosome("10"),0,(int)myRefSequence.genomeSize())));
    }

    @Test
    public void testChromosomeSequenceWithUnknownA() {

        // The String variables used for comparison in this junit are identical
        // to the chromosome sequences stored in the test file referenced.
        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(TEST_UNKNOWNA_FILE);

        Chromosome one = new Chromosome("1");
        String oneString = "AAGCTTGTGNNNNTTCTTCATC";
        byte[] chromBytes = myRefSequence.chromosomeSequence(one);
        //String byteToStr = nucleotideBytetoString(chromBytes);
        String byteToStr = NucleotideAlignmentConstants.nucleotideBytetoString(chromBytes);
        assertEquals(oneString.length(),chromBytes.length);
        assertTrue(oneString.equals(byteToStr));
    }

    /**
     * Test method for {@link net.maizegenetics.dna.map.HalfByteGenomeSequence#chromosomeSequence(net.maizegenetics.dna.map.Chromosome)}.
     */

    @Test
    public void testChromosomeSequenceChromosomeMultipleChromosome() {

        // The String variables used for comparison in this junit are identical
        // to the chromosome sequences stored in the test file referenced.
        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(TEST_SHORT_CHROM3456_FILE);

        Chromosome three = new Chromosome("3");
        Chromosome four = new Chromosome("4");
        Chromosome five = new Chromosome("5");
        Chromosome six = new Chromosome("6");

        String threeString = "AAGCTTGTGAAGGTTCTTCATCCCCACATGTGCTAAGCGGCGATGCCACAGCCAGCCCAT";
        byte[] chromBytes = myRefSequence.chromosomeSequence(three);
        String byte3ToStr = NucleotideAlignmentConstants.nucleotideBytetoString(chromBytes);
        assertEquals(threeString.length(),chromBytes.length);
        assertTrue(threeString.equals(byte3ToStr));

        String fourString = "AAGCTTCATGCAACTGCTGCCAGAAATGAGAGGTGAACTGCGTTCCTCTATCTGACACTA" +
                "TCTTCTTTGGCACACCATGAAGACAAACGATCCGAGACATATACAATTCTGCCAACACTG" +
                "CACTGCTATAGTTGGTCTTGACAGGTATGAAGTGGGCTGACTTGGTCAAGCGGTCCACTA" +
                "CTACCCAAATGGAATCGTAGCCGGCTCGAGTGCGAGGCAATCCGACTATGAAATCCATA";
        chromBytes = myRefSequence.chromosomeSequence(four);
        String byte4ToStr = NucleotideAlignmentConstants.nucleotideBytetoString(chromBytes);

        assertEquals(fourString.length(),chromBytes.length);  //should this return the real sequence?
        assertTrue(fourString.equals(byte4ToStr.toString()));

        String fiveString = "CAACATTTCTCCTATTCATATCTCTATCATGCTTATGGGAGGAACTAGAGGCAAACATGG" +
                "CATGAGAATCATAAACATGTGAATCAATATCATTATAAGCATTTCTAGCATTAGCATTTC" +
                "TAGCATGTCTCCTATCATTATACATAAAAGCATGGTTCTTTTTAGCACTACTAGCCATAG" +
                "GGGCCTTCCCCTTCCTCCTTGGTGGAGATGGGAGCCTTATGGCTTGTCAAGTTCAGGGAT" +
                "TCCCTCTTGAAGCCAAGACCATCCTTAATAGAGGGGTGTCTACCACTAGTGTAGGCATCC" +
                "CTTGCAAATTTTAGTTTATCAAAATCATTTTTGCTAGTCTTAAGTTGAGCATTAAGACTA" +
                "GCCAATTCATCATTCAATTTGGAAATTGAAACTAAGTGATCACTACAAGCATCAACATCG" +
                "AAATCTTTACATCTCTTACAGATAGTAACATGCTCTACGCAAGAGTTAGATTTACTAGCT" +
                "ACTTCTAGTTTAGCATTTAAATCATCATTAACACTTTTTAAAGTAGCAATGGTTTCATGA" +
                "CAAGTAGATAGTTCATAAGAAAGCATTTCATTTCTTTTAACTTCTAAAGCATGAGATTTT" +
                "TGTGCTTCTACAAATTTATCATGTTCTTCATATAACAGATCCTCTTGCTTTTCTAGTAAC" +
                "CTATTTTTGTCATTCAAGGCATCAATTAACTCGTTAATTTTGCCAATTTTAGTTCTATCT" +
                "AATCCTTTAAATAAACTAGAATAGTCTATTTCATCGTCGTCACTAGAATCATTATCACTA" +
                "GAAGAATCATAAGTAGTACTGTCTCGAGTACATACCTTCTTCTCCTTTGCCATGAGGCAT" +
                "GTGTGACGCTCGTTGGGGAAGAGGGAAAACTTGTTGAAGGCGGTGGCGGCGAGTCCTTCA" +
                "TTTTCGGAGTCGGAGGAGGAGCAATCCGAATCCCACTCCTTCCCGATATGTGCTTCGCCC" +
                "TTAGCCTTCTTGTAGTTCTTTTTCTTCTCCCTCTTGTTTCCTTTTTCATGGTCACTTTCA" +
                "TTATCGGGACAGTTAGCAATAAAATGACCAAGCTTACCACATTTGAAGCATGAGCGCTTC" +
                "CCCTTTGTCTTGGTCTTGCTTGGCTGTCCCTTGCGACCCTTAAGCGCCGTCTTGAAGCGC" +
                "TTGATGATGAGGGCCATCTCTTCATCATTGAGCCCTGCCGCCTCAACTTGCGCCACCTTG" +
                "CTAGGTAACGCCTCCTTGCTCCTTGTTGCCTTGAGAGCAATGGGTTGAGGCTCATGGATT" +
                "AGACCGTTCAATGCGTCGTCCACGTACCTCGCCTCCTTGATCATCATTCGCCCGCTTACG" +
                "AATTTTCCAAGAACTTCTTCGGGTGACATTTTGGTGTACCTGGGATTCTCACGAATATTA" +
                "TTCACCAAATGTGGATCAAGAACGGTAAATGACCTTAGCATTAGGTGGACGACGTCGTGG" +
                "TCCGTCCATCGCGTGCTTCCGTAGCTCCTTATTTTGTTGATAAGGGTCTTGAGCCGGTTG" +
                "TATGTTTGGGTTGGCTCCTCGCCCCTTATCATCGCGAATCTTCCAAGCTCTCCTTCTATC";
        chromBytes = myRefSequence.chromosomeSequence(five);
        String byte5ToStr = NucleotideAlignmentConstants.nucleotideBytetoString(chromBytes);
        assertEquals(fiveString.length(),chromBytes.length);
        assertTrue(fiveString.equals(byte5ToStr));

        String sixString = "TCAACCCTTTCCCTCTCTCAAACGGTCACCTAGACCGAGTGAGCTTCTCTTCTCAATCAA" +
                "ACGGAACACAAAGTTCCCGCAAGACCACCACACAATTGGTGTCTCTTGCCTTGGTTACAC" +
                "GAATCAAGCTGATTCATCAGAACCCTCCCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                "NNNNNNNNNAACCACCAAAATATGGCCGTTGGAAGACTGCTGTCGCATGGCGCACTGGAC" +
                "ACTGTCCCGTGCGCCAGCCACGTCAGCAGACCGTTGAGGTTCGACCGTTGGAGCTCTGAC";
        chromBytes = myRefSequence.chromosomeSequence(six);
        String byte6ToStr = NucleotideAlignmentConstants.nucleotideBytetoString(chromBytes);
        assertEquals(sixString.length(),chromBytes.length);
        assertTrue(sixString.equals(byte6ToStr));


        String genome=threeString+fourString+fiveString+sixString;
        assertEquals(genome.length(),myRefSequence.genomeSize());
        assertEquals(4,myRefSequence.numberOfChromosomes());
        assertEquals("AAGCTTGTGAA", NucleotideAlignmentConstants.nucleotideBytetoString(myRefSequence.genomeSequence(0, 10)));  //tests of zero base start of genome
        assertEquals("GCTTGTGAA", NucleotideAlignmentConstants.nucleotideBytetoString(myRefSequence.genomeSequence(2, 10)));
        int chr3Size=myRefSequence.chromosomeSize(new Chromosome("3"));
        assertEquals("CATAAG", NucleotideAlignmentConstants.nucleotideBytetoString(myRefSequence.genomeSequence(chr3Size-3, chr3Size+2)));  //boundary between chrom 3 & 4
        assertEquals(NucleotideAlignmentConstants.nucleotideBytetoString(myRefSequence.genomeSequence(0, myRefSequence.genomeSize()-1)),genome);
    }


    /**
     * Test method for {@link net.maizegenetics.dna.map.HalfByteGenomeSequence#chromosomeSequence(net.maizegenetics.dna.map.Chromosome, int, int)}.
     */
    @Test
    public void testChromosomeSequenceChromosomeIntInt() {
        // myShortFasta line is the sequence line from TEST_SHORT_CHROM3_FILE
        // that we are grabbing portions of for this test.  Displayed here for
        // ease in reading the start/end sequences that are tested.
        //
        // String myShortFasta = "AAGCTTGTGAAGGTTCTTCATCCCCACATGTGCTAAGCGGCGATGCCACAGCCAGCCCAT";
        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(TEST_SHORT_CHROM3_FILE);

        Chromosome three = new Chromosome("3");

        String fourToSixteen = "CTTGTGAAGGTTC"; 
        byte[] partial1 = myRefSequence.chromosomeSequence(three, 4,16);

        String partial1Str = NucleotideAlignmentConstants.nucleotideBytetoString(partial1);
        System.out.println("TestCSCIntInt: fourTOsixteen: "  + fourToSixteen + " partial1Str: " + partial1Str);
        assertTrue(fourToSixteen.equals(partial1Str));

        String fourToFifteen = "CTTGTGAAGGTT";
        byte[] partial2 = myRefSequence.chromosomeSequence(three, 4,15);
        String partial2Str = NucleotideAlignmentConstants.nucleotideBytetoString(partial2);
        System.out.println("TestCSCIntInt: fourTOFifteen: "  + fourToFifteen + " partial2Str: " + partial2Str);
        assertTrue(fourToFifteen.equals(partial2Str));

        String threeToFourteen = "GCTTGTGAAGGT";       
        byte[] partial3 = myRefSequence.chromosomeSequence(three, 3,14);


        String partial3Str = NucleotideAlignmentConstants.nucleotideBytetoString(partial3);

        System.out.println("TestCSCIntInt: threeToFourteen: "  + threeToFourteen + " partial3Str: " + partial3Str);
        assertTrue(threeToFourteen.equals(partial3Str));

        // String myShortFasta = "AAGCTTGTGAAGGTTCTTCATCCCCACATGTGCTAAGCGGCGATGCCACAGCCAGCCCAT";
        String threeToThirteen = "GCTTGTGAAGG";     
        byte[] partial4 = myRefSequence.chromosomeSequence(three, 3,13);

        String partial4Str = NucleotideAlignmentConstants.nucleotideBytetoString(partial4);
        System.out.println("TestCSCIntInt: threeToThirteen: "  + threeToThirteen + " partial4Str: " + partial4Str);
        assertTrue(threeToThirteen.equals(partial4Str));
        
        System.out.println("LCJ - before test of 1 allele byte array");
        // LCJ - added test 1 position, 12/7/15
        // "chromosomeSequence() shifts the request to 0 based.
        byte[] oneAllele = myRefSequence.chromosomeSequence(three, 3,13);
        byte[] myOneAllele = new byte[1];
        myOneAllele[0] = oneAllele[0]; // this is 0 based, but we grabbed chromsomeSequence based on 1
        String oneAlleleString = NucleotideAlignmentConstants.nucleotideBytetoString(myOneAllele);
        System.out.println("LCJ - oneALlele value:" + oneAlleleString);
    }

    /**
     * Test method for {@link net.maizegenetics.dna.map.HalfByteGenomeSequence#readReferenceGenomeChr(java.lang.String, int)}.
     */
    @Test
    public void testReadReferenceGenomeChr() {
        // The "myShortFast" string  matches the very short chromosome string from the TEST_SHORT_CHROM3_FILE
        // This test case is to ensure basic encoding/decoding is working
        String myShortFasta = "AAGCTTGTGAAGGTTCTTCATCCCCACATGTGCTAAGCGGCGATGCCACAGCCAGCCCAT";

        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(TEST_SHORT_CHROM3_FILE);

        Chromosome three = new Chromosome("3");
        byte[] chrom3bytes = myRefSequence.chromosomeSequence(three);
        if (chrom3bytes == null) {
            fail("Nothing read from REFERENCE_GENOME_FILE for chromosome 3");
        } 

        // Translate the returned bytes to the nucleotide string
        String chrom3bytesStr = NucleotideAlignmentConstants.nucleotideBytetoString(chrom3bytes);
        boolean seqEquals = myShortFasta.equals(chrom3bytesStr);
        System.out.println("Test CSC: seqEquals: " + seqEquals + " chrom3bytesStr: " + chrom3bytesStr + "  myShortFasta: "
                + myShortFasta);
        assertTrue(myShortFasta.equals(chrom3bytesStr));
    }

    @Test
    public void testFullRefCoordinateToChromCoordinate() {
        System.out.println("Begin junit testFulLRefCoordinateToChromCoordinate");
        Map<Long, Tuple<Chromosome, Integer>> mappedCoords = new HashMap<Long, Tuple<Chromosome, Integer>>();
        
        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(TEST_SHORT_CHROM3456_FILE);
        System.out.println("TestFullRefCtoCM: myRef genomeSize = " + myRefSequence.genomeSize());
        Chromosome three = new Chromosome("3");
        Chromosome four = new Chromosome("4");
        Chromosome five = new Chromosome("5");
        Chromosome six = new Chromosome("6");
        
        
        ArrayList<Long> csToMap = new ArrayList<Long>();
        int size3 = myRefSequence.chromosomeSize(three);
        int size4 = myRefSequence.chromosomeSize(four);
        int size5 = myRefSequence.chromosomeSize(five);
        int size6 = myRefSequence.chromosomeSize(six);
        System.out.println("Size of chrom sequences: chrom3: " + size3 + " chrom4: " 
          + size4 + " chrom5: " + size5 + " chrom6: " + size6);
        
        long size3Long = size3;
        long item1 = size3Long -1;
        csToMap.add(item1); // want the last value of chrom3, need to send it as 0 based
        mappedCoords = myRefSequence.fullRefCoordinateToChromCoordinate(csToMap);
        Tuple<Chromosome, Integer> mappedChromPos = mappedCoords.get(item1);
        assertTrue(mappedChromPos.x.equals(three));
        assertTrue(mappedChromPos.y.equals(size3-1));
        
        // this should get me chrom 4, index 5, because the size is 1 greater than the array
        // position, ie if size is 60, array position of last one is 59.  If you add 5,
        // you are indexing into 60+5, and the values 60 to X are in chrom 4, so 65
        // gets you 6th coordinate in a 0 based system. array[5] is the 6th coordinate.
        
        csToMap.clear();
        long item2 = size3Long + 5;
        csToMap.add(item2); // should map to Chromosome 4, position 6 (array[5])
        mappedCoords = myRefSequence.fullRefCoordinateToChromCoordinate(csToMap);
        mappedChromPos = mappedCoords.get(item2);
        assertTrue(mappedChromPos.x.equals(four));
        assertTrue(mappedChromPos.y.equals(5));
        
        csToMap.clear();
        long item3 = size3 + size4 + size5 + 10;
        csToMap.add(item3); // should map to chromosome 6, position 10               
        mappedCoords = myRefSequence.fullRefCoordinateToChromCoordinate(csToMap);
        mappedChromPos = mappedCoords.get(item3);
        assertTrue(mappedChromPos.x.equals(six));
        assertTrue(mappedChromPos.y.equals(10));
        
        csToMap.clear();
        
        // Tests above done individually to make the asserts easier.  Now
        // verify we can send in an array and get the correct number of responses.
        csToMap.add(item1);
        csToMap.add(item2);
        csToMap.add(item3);
        mappedCoords = myRefSequence.fullRefCoordinateToChromCoordinate(csToMap);
        assertEquals(mappedCoords.size(), 3);
        Tuple<Chromosome,Integer> chromPos1 = new Tuple<>(three,59);
        Tuple<Chromosome,Integer> chromPos2 = new Tuple<>(four,5);
        Tuple<Chromosome,Integer> chromPos3 = new Tuple<>(six,10);
        
        assertTrue(mappedCoords.get(item1).equals(chromPos1));
        assertTrue(mappedCoords.get(item2).equals(chromPos2));
        assertTrue(mappedCoords.get(item3).equals(chromPos3));
        
        System.out.println("\ntestFullRefCoordinateToChromCoordinate passed !!");
    }
    
    @Test
    public void testGenotypeByPosition() {

        Chromosome three = new Chromosome("3");
        String threeString = "AAGCTTGTGAAGGTTCTTCATCCCCACATGTGCTAAGCGGCGATGCCACAGCCAGCCCAT";        
        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(three, threeString);
        
        byte allele5byte = (byte)3;// tassel encoding of "T", which is allele at physical position 5
        byte genotypeByte5 = myRefSequence.genotype(three,5);
        assertEquals("Byte from genotype doesn't match t-3",allele5byte,genotypeByte5);
        
        byte allele1byte = (byte)0;// tassel encoding of "A", which is allele at physical position 1
        byte genotypeByte1 = myRefSequence.genotype(three, 1);
        assertEquals("First byte chrom3 failed to match",allele1byte,genotypeByte1);
        
        byte alleleLastByte = (byte)3;
        byte genotypeByteLast = myRefSequence.genotype(three, threeString.length());
        assertEquals("Last byte chrom3 failed to match",alleleLastByte,genotypeByteLast);
        
    }
    
    @Test
    public void testGenotypeByPositionObject() {

        Chromosome three = new Chromosome("3");
        String threeString = "AAGCTTGTGAAGGTTCTTCATCCCCACATGTGCTAAGCGGCGATGCCACAGCCAGCCCAT";        
        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(three, threeString);
        
        byte allele5byte = (byte)3;// tassel encoding of "T", which is allele at physical position 5
        Position allele5bytePos = new GeneralPosition.Builder(three, 5).build();       
        byte genotypeByte5 = myRefSequence.genotype(three,allele5bytePos);       
        assertEquals("Byte from genotype by Position doesn't match t-3",genotypeByte5,allele5byte);
        
        byte allele1byte = (byte)0;// tassel encoding of "A", which is allele at physical position 1
        Position allele1bytePos = new GeneralPosition.Builder(three,1).build();
        byte genotypeByte1 = myRefSequence.genotype(three, allele1bytePos);
        assertEquals("First byte chrom3 failed to match",allele1byte,genotypeByte1);
        
        byte alleleLastByte = (byte)3;
        Position alleleLastBytePos = new GeneralPosition.Builder(three,threeString.length()).build();
        byte genotypeByteLast = myRefSequence.genotype(three, alleleLastBytePos);
        assertEquals("Last byte chrom3 failed to match",alleleLastByte,genotypeByteLast);
 
    }
    
    @Test
    public void testGenotypeAsStringByPosition() {

        Chromosome three = new Chromosome("3");
        String threeString = "AAGCTTGTGAAGGTTCTTCATCCCCACATGTGCTAAGCGGCGATGCCACAGCCAGCCCAT";        
        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(three, threeString);
        
        String allele5 = "T";
        String genotypePos5 = myRefSequence.genotypeAsString(three,5);
        assertEquals("Byte from genotype doesn't match t-3",allele5,genotypePos5);
        
        String allele1byte = "A";
        String genotypePos1 = myRefSequence.genotypeAsString(three, 1);
        assertEquals("First byte chrom3 failed to match",allele1byte,genotypePos1);
        
        String alleleLast = "T";
        String genotypePosLast = myRefSequence.genotypeAsString(three, threeString.length());
        assertEquals("Last byte chrom3 failed to match",alleleLast,genotypePosLast);
        
    }
    
    @Test
    public void testGenotypeAsStringByPositionObject() {

        Chromosome three = new Chromosome("3");
        String threeString = "AAGCTTGTGAAGGTTCTTCATCCCCACATGTGCTAAGCGGCGATGCCACAGCCAGCCCAT";        
        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(three, threeString);
        
        String allele5byte = "T";// tassel encoding of "T", which is allele at physical position 5
        Position allele5bytePos = new GeneralPosition.Builder(three, 5).build();       
        String genotype5 = myRefSequence.genotypeAsString(three,allele5bytePos);       
        assertEquals("Byte from genotype by Position doesn't match t-3",genotype5,allele5byte);
        
        String allele1byte = "A";// tassel encoding of "A", which is allele at physical position 1
        Position allele1bytePos = new GeneralPosition.Builder(three,1).build();
        String genotype1 = myRefSequence.genotypeAsString(three, allele1bytePos);
        assertEquals("First byte chrom3 failed to match",allele1byte,genotype1);
        
        String alleleLastByte = "T";
        Position alleleLastPos = new GeneralPosition.Builder(three,threeString.length()).build();
        String genotypeLast = myRefSequence.genotypeAsString(three, alleleLastPos);
        assertEquals("Last byte chrom3 failed to match",alleleLastByte,genotypeLast);
 
    }
    
    @Test
    public void testGenotypeStartEndAsString() {

        Chromosome three = new Chromosome("3");
        String threeString = "AAGCTTGTGAAGGTTCTTCATCCCCACATGTGCTAAGCGGCGATGCCACAGCCAGCCCAT";        
        GenomeSequence myRefSequence = GenomeSequenceBuilder.instance(three, threeString);
        
        String alleleSevenToTen = "GTGA";
        String genotypeAsString = myRefSequence.genotypeAsString(three,7,10);
        assertEquals("String from genotypAsString doesn't match chrom seq:",alleleSevenToTen,genotypeAsString);
        
        String lastThree = "CAT";
        int start = threeString.length()-2;
        int end = threeString.length();
        genotypeAsString = myRefSequence.genotypeAsString(three, start,end);
        assertEquals("First 6 chrom sequence doesnt match genotypeAsString with start/end coordinates:",lastThree,genotypeAsString);
        
        String first5 = "AAGCT";
        genotypeAsString = myRefSequence.genotypeAsString(three, 1,5);
        assertEquals("Last byte chrom3 failed to match",first5,genotypeAsString);
 
    }
}
