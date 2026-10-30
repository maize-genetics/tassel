package net.maizegenetics.dna.snp.io;

import junit.framework.Assert;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.map.*;
import net.maizegenetics.dna.snp.FilterAndMaskGVCFGenomeSequence;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import static org.junit.Assert.assertNotNull;

/**
 * Created by zrm22 on 3/24/17.
 */
public class BuilderFromGVCFTest {

    //Class to Unit test GVCF import
    //Eventually the tests from BuilderFromVCFTest will reside here as well as GVCF is also valid VCF.  The import should be flexible to handle this

    String directoryName = GeneralConstants.DATA_DIR + "CandidateTests/VCFFiles/";
    String agpv3FileName = GeneralConstants.DATA_DIR + "CandidateTests/CHR8_Zea_mays.AGPv3.dna.toplevel.fa.gz";
    String agpv4FileName = GeneralConstants.DATA_DIR + "CandidateTests/CHR8_Zea_mays.AGPv4.dna.toplevel.fa.gz";
    String fullGVCFFileName = GeneralConstants.DATA_DIR + "CandidateTests/VCFFiles/B73_haplotype_caller_output.g.vcf.gz";

    //Made-up two-chromosome genome for the multi-chromosome tests. On chromosome 1 the sample is
    //homozygous T at site 11 and has no GVCF records for sites 21-24; chromosome 2 matches the reference.
    private static final String SIM_CHR1 = "ACGTACGTACGTACGTACGTACGTACGTAC";
    private static final String SIM_CHR2 = "TTGGCCAATTGGCCAATTGGCCAATTGGCC";
    private static final String SIM_SAMPLE_CHR1 = SIM_CHR1.substring(0, 10) + "T" + SIM_CHR1.substring(11);

    //Loading these takes about 40 seconds, so the tests share one copy. Tests must not modify them.
    private static GVCFGenomeSequence b73GVCFSequence;
    private static GenomeSequence agpv4Sequence;

    private GVCFGenomeSequence b73GVCFSequence() throws Exception {
        if (b73GVCFSequence == null) {
            b73GVCFSequence = (GVCFGenomeSequence) GVCFGenomeSequenceBuilder.instance(agpv4FileName, fullGVCFFileName);
        }
        return b73GVCFSequence;
    }

    private GenomeSequence agpv4Sequence() {
        if (agpv4Sequence == null) {
            agpv4Sequence = GenomeSequenceBuilder.instance(agpv4FileName);
        }
        return agpv4Sequence;
    }

    private GVCFGenomeSequence simulatedGVCFSequence() throws Exception {
        Path dir = Paths.get(GeneralConstants.TEMP_DIR, "BuilderFromGVCFTest");
        Files.createDirectories(dir);
        Path reference = dir.resolve("simulatedReference.fa");
        Files.write(reference, Arrays.asList(">1", SIM_CHR1, ">2", SIM_CHR2));

        String refBlockFormat = "GT:DP:GQ:MIN_DP:PL\t0/0:40:99:38:0,89,1385";
        Path gvcf = dir.resolve("simulatedSample.g.vcf");
        Files.write(gvcf, Arrays.asList(
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIM",
                "1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=10\t" + refBlockFormat,
                "1\t11\t.\tG\tT,<NON_REF>\t600\t.\tDP=40\tGT:AD:DP:GQ:PL\t1/1:0,40,0:40:99:641,0,456,691,519,1210",
                "1\t12\t.\tT\t<NON_REF>\t.\t.\tEND=20\t" + refBlockFormat,
                "1\t25\t.\tA\t<NON_REF>\t.\t.\tEND=30\t" + refBlockFormat,
                "2\t1\t.\tT\t<NON_REF>\t.\t.\tEND=30\t" + refBlockFormat));

        return (GVCFGenomeSequence) GVCFGenomeSequenceBuilder.instance(reference.toString(), gvcf.toString());
    }

    private static Chromosome chromosomeNamed(GenomeSequence sequence, String name) {
        return sequence.chromosomes().stream()
                .filter(chr -> chr.getName().equals(name))
                .findFirst()
                .orElseThrow(() -> new AssertionError("No chromosome " + name));
    }

    @Test
    public void testSmallGVCF() {
        String fileName = directoryName + "refBlockDepthAboveFive.g.vcf";
        String refFileName = agpv3FileName;
        try {
            GVCFGenomeSequence sequence =  (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, fileName);
            Set<Chromosome> chromosomes =  sequence.chromosomes();
            String exportedString = "";
            //Loop through each chromosome
            for(Chromosome chr : chromosomes) {
                if(chr.getChromosomeNumber()==8) {
                    byte[] currentSequence = sequence.chromosomeSequence(chr, 10000000,10000211);
                    exportedString = NucleotideAlignmentConstants.nucleotideBytetoString(currentSequence);
                }
            }


            String expectedString = "TCTAATAATAGTAATTTAGGCATATATCAATTAAGCTAATTTGCTTTTATGCAAAATATATTTGTATACTATTATTAACAACA" +
                    "TGTCGGAGATATTTATGTGCTACATTTTTACTATTGAGGAGTGGAATGAAGAGTGTCATGTAAGTTACAGAGTAGAAACAAATTCTAGTAATGTATAAA" +
                    "ATCATTTCCTATCCTCCACCCTATGAATTC";

            Assert.assertEquals("GVCF DepthAbove5 Round Trip test fails:",expectedString,exportedString);

        }
        catch(Exception e) {
            Assert.fail("GVCF DepthAbove5 Round Trip test fails: Exception thrown.");
            e.printStackTrace();
        }
    }

    @Test
    public void testSmallGVCFFirstSiteStartsInBlock() {
        String fileName = directoryName + "refBlockDepthAboveFiveBigFirstRecord.g.vcf";
        String refFileName = agpv3FileName;
        try {
            GVCFGenomeSequence sequence =  (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, fileName);
            Set<Chromosome> chromosomes =  sequence.chromosomes();
            String exportedString = "";
            //Loop through each chromosome
            for(Chromosome chr : chromosomes) {
                if(chr.getChromosomeNumber()==8) {
                    byte[] currentSequence = sequence.chromosomeSequence(chr, 10000000,10000211);
                    exportedString = NucleotideAlignmentConstants.nucleotideBytetoString(currentSequence);
                }
            }


            String expectedString = "TCTAATAATAGTAATTTAGGCATATATCAATTAAGCTAATTTGCTTTTATGCAAAATATATTTGTATACTATTATTAACAACA" +
                    "TGTCGGAGATATTTATGTGCTACATTTTTACTATTGAGGAGTGGAATGAAGAGTGTCATGTAAGTTACAGAGTAGAAACAAATTCTAGTAATGTATAAA" +
                    "ATCATTTCCTATCCTCCACCCTATGAATTC";

            Assert.assertEquals("GVCF DepthAbove5 Round Trip test fails:",expectedString,exportedString);

        }
        catch(Exception e) {
            Assert.fail("GVCF DepthAbove5 Round Trip test fails: Exception thrown.");

            e.printStackTrace();
        }
    }

    @Test
    public void testSmallGVCFNsAtStart() {
        String fileName = directoryName + "refBlockDepthAboveFive.g.vcf";
        String refFileName = agpv3FileName;
        try {
            GVCFGenomeSequence sequence =  (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, fileName);
            Set<Chromosome> chromosomes =  sequence.chromosomes();
            String exportedString = "";
            //Loop through each chromosome
            for(Chromosome chr : chromosomes) {
                if(chr.getChromosomeNumber()==8) {
                    byte[] currentSequence = sequence.chromosomeSequence(chr, 9999990,10000211);
                    exportedString = NucleotideAlignmentConstants.nucleotideBytetoString(currentSequence);
                }
            }


            String expectedString = "NNNNNNNNNNTCTAATAATAGTAATTTAGGCATATATCAATTAAGCTAATTTGCTTTTATGCAAAATATATTTGTATACTATTATTAACAACA" +
                    "TGTCGGAGATATTTATGTGCTACATTTTTACTATTGAGGAGTGGAATGAAGAGTGTCATGTAAGTTACAGAGTAGAAACAAATTCTAGTAATGTATAAA" +
                    "ATCATTTCCTATCCTCCACCCTATGAATTC";

            Assert.assertEquals("GVCF DepthAbove5 Round Trip test fails:",expectedString,exportedString);

        }
        catch(Exception e) {
            Assert.fail("GVCF DepthAbove5 Round Trip test fails: Exception thrown");
            e.printStackTrace();
        }
    }

    @Test
    public void testSmallGVCFInsert() {
        String fileName = directoryName + "insertionGVCFFile.g.vcf";
        String refFileName = agpv3FileName;
        try {
            GVCFGenomeSequence sequence =  (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, fileName);
            Set<Chromosome> chromosomes =  sequence.chromosomes();
            //Loop through each chromosome
            String exportedString = "";
            for(Chromosome chr : chromosomes) {
                if(chr.getChromosomeNumber()==8) {
                    byte[] currentSequence = sequence.chromosomeSequence(chr, 10000000,10000211);
                    exportedString = NucleotideAlignmentConstants.nucleotideBytetoString(currentSequence);
                }
            }

            String expectedResult = "TCTAATAATAGTAATTTAGGCATATATCAATTAAGCTAATTTGCTTTTATGCAAAATATATTTGTATACTATTATTAACAA" +
                    "CATGTCGGAGATATTTATGTGCTACATTTTTACTATACCTTGAGGAGTGGAATGAAGAGTGTCATGTAAGTTACAGAGTAGAAACAAATTCTAGTAA" +
                    "TGTATAAAATCATTTCCTATCCTCCACCCTATGAATTC";

            Assert.assertEquals("GVCF Insertion Test Fails:",expectedResult,exportedString);

        }
        catch(Exception e) {
            Assert.fail("GVCF Insertion Test Fails: Exception thrown");
            e.printStackTrace();
        }
    }

    @Test
    public void testSmallGVCFDelete() {
        String fileName = directoryName + "deletionGVCFFile.g.vcf";
        String refFileName = agpv3FileName;
        try {
            GVCFGenomeSequence sequence =  (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, fileName);
            Set<Chromosome> chromosomes =  sequence.chromosomes();
            //Loop through each chromosome
            String exportedString = "";
            for(Chromosome chr : chromosomes) {
                if(chr.getChromosomeNumber()==8) {
                    byte[] currentSequence = sequence.chromosomeSequence(chr, 10000000,10000211);
                    exportedString = NucleotideAlignmentConstants.nucleotideBytetoString(currentSequence);
                }
            }

            String expectedResult = "TCTAATAATAGTAATTTAGGCATATATCAATTAAGCTAATTTGCTTTTATGCAAAATATATTTGTATACTATTATTAACAA" +
                    "CATGTCGGAGATATTTATGTGCTACATTTTTACTATGAGGAGTGGAATGAAGAGTGTCATGTAAGTTACAGAGTAGAAACAAATTCTAGTAATGTAT" +
                    "AAAATCATTTCCTATCCTCCACCCTATGAATTC";

            Assert.assertEquals("GVCF Deletion test fails:",expectedResult,exportedString);

        }
        catch(Exception e) {
            Assert.fail("GVCF Deletion test fails: Exception thrown.");
            e.printStackTrace();
        }
    }


    @Test
    public void testConsecSequence() {
        String fileName = directoryName + "consecGVCFTest.g.vcf";
        String refFileName = agpv3FileName;
        try {
            GVCFGenomeSequence sequence = (GVCFGenomeSequence) GVCFGenomeSequenceBuilder.instance(refFileName, fileName);
            Set<Chromosome> chromosomes = sequence.chromosomes();

            //Expected results
            ArrayList<ArrayList<Integer>> expectedList = new ArrayList<>();
            expectedList.add(new ArrayList<Integer>(Arrays.asList(8,10000000,10000211)));
            expectedList.add(new ArrayList<Integer>(Arrays.asList(8,20000000,20000012)));
            expectedList.add(new ArrayList<Integer>(Arrays.asList(8,30000000,30000005)));
            expectedList.add(new ArrayList<Integer>(Arrays.asList(8,30000007,30000011)));
            expectedList.add(new ArrayList<Integer>(Arrays.asList(8,30000013,30000020)));

            HashMap<Chromosome, ArrayList<ArrayList<Integer>>> consecRegions = sequence.getConsecutiveRegions();
            for(Chromosome chr : chromosomes) {
                int expectedCounter = 0;
                for(ArrayList<Integer> currentRange: consecRegions.get(chr)) {

                    Assert.assertEquals("Consecutive Sequence Error Chromosome is incorrect:"+chr.getChromosomeNumber()+":"+currentRange.toString(),
                            (int)expectedList.get(expectedCounter).get(0),chr.getChromosomeNumber());
                    Assert.assertEquals("Consecutive Sequence Error StartPosition is incorrect:"+chr.getChromosomeNumber()+":"+currentRange.toString(),
                            (int) expectedList.get(expectedCounter).get(1),(int)currentRange.get(0));
                    Assert.assertEquals("Consecutive Sequence Error EndPosition is incorrect:"+chr.getChromosomeNumber()+":"+currentRange.toString(),
                            (int)expectedList.get(expectedCounter).get(2),(int)currentRange.get(1));
                    expectedCounter++;
                }
            }
            sequence.writeFASTA(GeneralConstants.TEMP_DIR+"sampleConsecGVCFExport.fasta");
            //TODO compare this file to an expected one
        }
        catch(Exception e) {
            Assert.fail("Consecutive Sequence Test Fail: Exception thrown");
            e.printStackTrace();
        }
    }


    @Test
    public void filterAboveFiveGVCF() {
        String fileName = directoryName + "refBlockDepthOne.g.vcf";
        String refFileName = agpv3FileName;
        try {
            GVCFGenomeSequence sequence =  (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, fileName);
            GVCFGenomeSequence filteredSequence = FilterAndMaskGVCFGenomeSequence.getInstance(sequence,"DP",3,5,false);

            Set<Chromosome> chromosomes =  filteredSequence.chromosomes();
            //Loop through each chromosome
            String exportedString = "";
            for(Chromosome chr : chromosomes) {
                if(chr.getChromosomeNumber()==8) {
                    byte[] currentSequence = filteredSequence.chromosomeSequence(chr, 10000000, 10000211);
                    exportedString = NucleotideAlignmentConstants.nucleotideBytetoString(currentSequence);
                }
            }

            String expectedResult = "TCTAATAATAGTAATTTAGGCATATATCAATTAAGCTAATTTGCTTTTATGCAAAATATATTTGTATACTATTATTAACA" +
                    "ACATGTCGGAGATATTTATGTGCTACATTTTTACTATTC";

            Assert.assertEquals("GVCF Filter Depth above 5 test fails:",expectedResult,exportedString);
        }
        catch(Exception e) {
            Assert.fail("GVCF FilterDepth above 5 test fails: Exception thrown");
            e.printStackTrace();
        }
    }

    @Test
    public void maskAboveFiveGVCF() {
        String fileName = directoryName + "refBlockDepthOne.g.vcf";
        String refFileName = agpv3FileName;
        try {
            GVCFGenomeSequence sequence =  (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, fileName);

            GVCFGenomeSequence filteredSequence = FilterAndMaskGVCFGenomeSequence.getInstance(sequence,"DP",3,5,true);

            Set<Chromosome> chromosomes =  filteredSequence.chromosomes();
            //Loop through each chromosome
            String exportedString = "";
            for(Chromosome chr : chromosomes) {
                if(chr.getChromosomeNumber()==8) {
                    byte[] currentSequence = filteredSequence.chromosomeSequence(chr, 10000000,10000211);
                    exportedString = NucleotideAlignmentConstants.nucleotideBytetoString(currentSequence);
                }
            }
            String expectedResult = "TCTAATAATAGTAATTTAGGCATATATCAATTAAGCTAATTTGCTTTTATGCAAAATATATTTGTATACTATTATTAACA" +
                    "ACATGTCGGAGATATTTATGTGCTACATTTTTACTATTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNC";

            Assert.assertEquals("GVCF Mask Depth>5 test fails:",expectedResult,exportedString);
        }
        catch(Exception e) {
            Assert.fail("GVCF Mask Depth>5 test fails: Exception thrown");
            e.printStackTrace();
        }
    }
    @Test
    public void testHetInsertions() {
        //Test region is here
        //8:105839363:105844126
        try {
            GVCFGenomeSequence sequence = b73GVCFSequence();
            GenomeSequence refSequence = agpv4Sequence();
            HashMap<String, Chromosome> chrMap = new HashMap<>();


            Object[] chrs = sequence.chromosomes().toArray();
            for(int i = 0; i < chrs.length; i++) {
                chrMap.put(""+((Chromosome)chrs[i]).getChromosomeNumber(),(Chromosome)chrs[i]);
            }
            String pulledSequence = ""+NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chrMap.get("8"),
                    105839363,105844126));
            String refString = ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    105839363,105839510));
            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    105839578,105840351));
            //Here is the insert
            refString += "N";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    105840353,105841913));

            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    105841968,105843107));

            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    105843455,105843595));

            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                    "NNNNNNNNNNNNNNNNNNNNNN";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    105843718,105844126));

            Assert.assertEquals("Het SNP Insertion test fails:",refString,pulledSequence);
        }
        catch(Exception e) {
            Assert.fail("Het SNP Insertion test fails: Exception thrown");
            e.printStackTrace();
        }
    }
    @Test
    public void testHetDeletions() {
        //8:83443248:83445925
        try {
            GVCFGenomeSequence sequence = b73GVCFSequence();
            GenomeSequence refSequence = agpv4Sequence();
            HashMap<String, Chromosome> chrMap = new HashMap<>();


            Object[] chrs = sequence.chromosomes().toArray();
            for(int i = 0; i < chrs.length; i++) {
                chrMap.put(""+((Chromosome)chrs[i]).getChromosomeNumber(),(Chromosome)chrs[i]);
            }
            String pulledSequence = ""+NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chrMap.get("8"),
                    83443248,83445925));
            String refString = ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    83443248,83443722));
            refString += "N";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    83443724,83444607));
            //here is a het SNP
            refString += "N";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    83444609,83445062));
            //84 Ns
            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    83445149,83445194));

            refString += "N";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    83445196,83445232));

            refString += "N";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    83445234,83445259));

            refString += "N";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    83445261,83445264));
            refString +="N";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    83445266,83445853));

            //Het deletion
            refString += "NN";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    83445856,83445925));

            Assert.assertEquals("Het SNP Deletion test fails:",refString,pulledSequence);
        }
        catch(Exception e) {
            Assert.fail("Het SNP Deletion test fails: Exception thrown");
            e.printStackTrace();
        }
    }
    @Test
    public void testHetSNPs() {
        //8:15917882:15920151
        try {
            GVCFGenomeSequence sequence = b73GVCFSequence();
            GenomeSequence refSequence = agpv4Sequence();
            HashMap<String, Chromosome> chrMap = new HashMap<>();


            Object[] chrs = sequence.chromosomes().toArray();
            for(int i = 0; i < chrs.length; i++) {
                chrMap.put(""+((Chromosome)chrs[i]).getChromosomeNumber(),(Chromosome)chrs[i]);
            }
            String pulledSequence = ""+NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chrMap.get("8"),
                    15917882,15920151));
            String refString = ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    15917882,15918042));
            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
            refString+=NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    15918109,15919048));
            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
            refString += NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    15919107,15919714));
            //ZRM22 added in an N for the het call
            refString += "N";
            refString += NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    15919716,15919821));

            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

            refString += NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    15919859,15919863));

            refString += "NN";

            refString += NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    15919866,15919869));

            refString += "N";

            refString += NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    15919871,15919997));

            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

            refString += NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    15920099,15920151));

            Assert.assertEquals("Het SNP test fails:",refString,pulledSequence);
        }
        catch(Exception e) {
            Assert.fail("Het SNP test fails: Exception thrown");
            e.printStackTrace();
        }
    }
    @Test
    public void test0Coverage() {
        //8:8937202-8939694 is inside the zero-depth block 8:8931596-8946845, so every site is exported as N
        try {
            GVCFGenomeSequence sequence = b73GVCFSequence();
            HashMap<String, Chromosome> chrMap = new HashMap<>();


            Object[] chrs = sequence.chromosomes().toArray();
            for(int i = 0; i < chrs.length; i++) {
                chrMap.put(""+((Chromosome)chrs[i]).getChromosomeNumber(),(Chromosome)chrs[i]);
            }
            String pulledSequence = ""+NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chrMap.get("8"),
                    8937202,8939694));

            Assert.assertEquals("Zero coverage test fails:","N".repeat(8939694 - 8937202 + 1),pulledSequence);

        }
        catch(Exception e) {
            Assert.fail("Zero coverage test fails: Exception thrown");
            e.printStackTrace();
        }
    }
    @Test
    public void testRefHomozygous() {
        //should match reference
        //8:15484032:15486280
        try {
            GVCFGenomeSequence sequence = b73GVCFSequence();
            GenomeSequence refSequence = agpv4Sequence();

            HashMap<String, Chromosome> chrMap = new HashMap<>();


            Object[] chrs = sequence.chromosomes().toArray();
            for(int i = 0; i < chrs.length; i++) {
                chrMap.put(""+((Chromosome)chrs[i]).getChromosomeNumber(),(Chromosome)chrs[i]);
            }
            String pulledSequence = ""+NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chrMap.get("8"),
                    15484032,15486280));

            String pulledRefSequence = ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    15484032,15486280));

            Assert.assertEquals("Ref Homozygous test fails:",pulledRefSequence,pulledSequence);

        }
        catch(Exception e) {
            Assert.fail("Ref Homozygous test fails: Excpetion thrown");
            e.printStackTrace();
        }
    }
    @Test
    public void testAltHomozygousDP1or2() {
        //put Ns in the where this is true
        //Good test sequence is here:
        //8:62915759:62918979
        try {
            GVCFGenomeSequence sequence = b73GVCFSequence();
            GenomeSequence refSequence = agpv4Sequence();
            HashMap<String, Chromosome> chrMap = new HashMap<>();


            Object[] chrs = sequence.chromosomes().toArray();
            for(int i = 0; i < chrs.length; i++) {
                chrMap.put(""+((Chromosome)chrs[i]).getChromosomeNumber(),(Chromosome)chrs[i]);
            }
            String pulledSequence = ""+NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chrMap.get("8"),
                    62915759,62918979));
            String refString = ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    62915759,62916511));


            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNN";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                   62916540,62917929));


            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    62918050,62918746));

            //Should be a TG but there is low depth for the alt
            refString += "NN";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    62918748,62918795));

            refString += "N";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    62918797,62918979));


            Assert.assertEquals("Low Alt Homozygous test fails:",refString,pulledSequence);
        }
        catch(Exception e) {
            Assert.fail("Low Alt Homozygous test fails: Exception thrown.");
            e.printStackTrace();
        }
    }
    @Test
    public void testAltHomozygousDPAbove3() {
        //Put alt alleles in
        //good test sequence is here:
        //8:103920155:103926901
        try {
            GVCFGenomeSequence sequence = b73GVCFSequence();
            GenomeSequence refSequence = agpv4Sequence();
            HashMap<String, Chromosome> chrMap = new HashMap<>();


            Object[] chrs = sequence.chromosomes().toArray();
            for(int i = 0; i < chrs.length; i++) {
                chrMap.put(""+((Chromosome)chrs[i]).getChromosomeNumber(),(Chromosome)chrs[i]);
            }
            String pulledSequence = ""+NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chrMap.get("8"),
                    103920155,103926901));
            //This region starts with Ns
            String refString = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNN";


            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    103920465,103920602));

            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    103920745,103924851));

            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    103924903,103926649));

            refString += "TC";

            refString += ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("8"),
                    103926651,103926868));

            refString += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

            Assert.assertEquals("Low Alt Homozygous test fails:",refString,pulledSequence);
        }
        catch(Exception e) {
            Assert.fail("Low Alt Homozygous test fails: Excpetion thrown");
            e.printStackTrace();
        }
    }

    //Pulls regions from both chromosomes of one loaded sequence, as the old anchor-region harnesses did
    //SuccessCondition: covered sites come from the sample, the uncovered sites 21-24 of chromosome 1 are N
    @Test
    public void pullRegionsAcrossChromosomes() throws Exception {
        GVCFGenomeSequence sequence = simulatedGVCFSequence();
        Chromosome chr1 = chromosomeNamed(sequence, "1");
        Chromosome chr2 = chromosomeNamed(sequence, "2");

        Assert.assertEquals("Region around the SNP", SIM_SAMPLE_CHR1.substring(4, 15),
                NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chr1, 5, 15)));
        Assert.assertEquals("Region across the uncovered sites", SIM_CHR1.substring(17, 20) + "NNNN" + SIM_CHR1.substring(24, 27),
                NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chr1, 18, 27)));
        Assert.assertEquals("Region on chromosome 2", SIM_CHR2.substring(2, 8),
                NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chr2, 3, 8)));
        Assert.assertEquals("Region on chromosome 1 pulled after chromosome 2", SIM_SAMPLE_CHR1.substring(0, 12),
                NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chr1, 1, 12)));
    }

    //Writes the sample's covered regions of every chromosome to FASTA
    //SuccessCondition: one record per run of consecutive covered sites, with the sample's sequence
    @Test
    public void writeFASTAOfAllChromosomes() throws Exception {
        GVCFGenomeSequence sequence = simulatedGVCFSequence();
        Path outFile = Paths.get(GeneralConstants.TEMP_DIR, "BuilderFromGVCFTest", "simulatedSample.fa");
        sequence.writeFASTA(outFile.toString());

        Map<String, String> expected = new HashMap<>();
        expected.put(">Chr_1_StartSite_1_EndSite_20", SIM_SAMPLE_CHR1.substring(0, 20));
        expected.put(">Chr_1_StartSite_25_EndSite_30", SIM_CHR1.substring(24, 30));
        expected.put(">Chr_2_StartSite_1_EndSite_30", SIM_CHR2);

        List<String> lines = Files.readAllLines(outFile);
        Map<String, String> written = new HashMap<>();
        for (int i = 0; i + 1 < lines.size(); i += 2) {
            written.put(lines.get(i), lines.get(i + 1));
        }
        Assert.assertEquals("FASTA should have one header and one sequence line per region", expected.size() * 2, lines.size());
        Assert.assertEquals(expected, written);
    }

//TODO clean up these old tests

//    @Test
//    public void attemptSmallVCFFileImport() {
//        //We should be able to import a vcf file using the GVCF import
//        //If not we have something wrong
//        //TODO replace this with BuilderFromGVCF
//        String fileName = directoryName+"correctVCF1.vcf";
//        BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
//        GenotypeTable table = builder.buildAndSortInMemory();
//        assertNotNull("Standard GVCF Test Table should not be null",table);
//        //TODO also check to make sure the calls are correct for each position and taxa
//    }
//
//    @Test
//    public void smallGVCFJustHeaderChange() {
//        //Test to make sure if we just add in the header, but the rest of it is normal vcf
//        //TODO replace this with BuilderFromGVCF
//        String fileName = directoryName+"correctGVCFOnlyHeader.g.vcf";
//        BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
//        GenotypeTable table = builder.buildAndSortInMemory();
//        assertNotNull("Standard GVCF Test Table should not be null",table);
//        //TODO also check to make sure the calls are correct for each position and taxa
//    }
//
//    @Test
//    public void smallGVCFNoBlockingImport() {
//        //Load in a small non-blocked GVCF file and make sure the GenotypeTable is correct(this might actually be the easiest test to hit)
//        //TODO replace this with BuilderFromGVCF
//        String fileName = directoryName+"correctGVCFNonBlocking.g.vcf";
//        BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
//        GenotypeTable table = builder.buildAndSortInMemory();
//        assertNotNull("Standard GVCF Test Table should not be null",table);
//        //TODO also check to make sure the calls are correct for each position and taxa
//    }
//
//    @Test
//    public void smallGVCFFileImport() {
//        //Load in a small blocked GVCF file and make sure the corresponding GenotypeTable has the correct calls for the correct positions
//        //TODO replace this with BuilderFromGVCF
//        String fileName = directoryName+"correctGVCF.g.vcf";
//        BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName);
//        GenotypeTable table = builder.buildAndSortInMemory();
//        assertNotNull("Standard GVCF Test Table should not be null",table);
//        //TODO also check to make sure the calls are correct for each position and taxa
//
//    }

}
