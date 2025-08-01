package net.maizegenetics.dna.snp.io;

import junit.framework.Assert;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.map.*;
import net.maizegenetics.dna.snp.FilterAndMaskGVCFGenomeSequence;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import org.junit.Ignore;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
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
        //String gvcfFileName = "/Users/zrm22/PanGenome/FASTAGeneration/B73_haplotype_caller_output.g.vcf";
        String gvcfFileName = fullGVCFFileName;
        String refFileName = agpv4FileName;
        try {
            System.out.println("Loading in GCVF");
            GVCFGenomeSequence sequence = (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, gvcfFileName);
            System.out.println("Loading in ref");
            GenomeSequence refSequence = (GenomeSequence)GenomeSequenceBuilder.instance(refFileName);
            System.out.println("Done loadin in files");
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
        //String gvcfFileName = "/Users/zrm22/PanGenome/FASTAGeneration/B73_haplotype_caller_output.g.vcf";
        String gvcfFileName = fullGVCFFileName;
        String refFileName = agpv4FileName;
        try {
            GVCFGenomeSequence sequence = (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, gvcfFileName);
            GenomeSequence refSequence = (GenomeSequence)GenomeSequenceBuilder.instance(refFileName);
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
//        String gvcfFileName = "/Users/zrm22/PanGenome/FASTAGeneration/B73_haplotype_caller_output.g.vcf";
        String gvcfFileName = fullGVCFFileName;
        String refFileName = agpv4FileName;
        try {
            GVCFGenomeSequence sequence = (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, gvcfFileName);
            GenomeSequence refSequence = (GenomeSequence)GenomeSequenceBuilder.instance(refFileName);
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
        //Should be 1 N exported..maybe we change this to a number of Ns
        //A good region is here:
        //8:8937202:8939694
//        String gvcfFileName = "/Users/zrm22/PanGenome/FASTAGeneration/B73_haplotype_caller_output.g.vcf";
        String gvcfFileName = fullGVCFFileName;
        String refFileName = agpv4FileName;
        try {
            GVCFGenomeSequence sequence = (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, gvcfFileName);
            HashMap<String, Chromosome> chrMap = new HashMap<>();


            Object[] chrs = sequence.chromosomes().toArray();
            for(int i = 0; i < chrs.length; i++) {
                chrMap.put(""+((Chromosome)chrs[i]).getChromosomeNumber(),(Chromosome)chrs[i]);
            }
            String pulledSequence = ""+NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chrMap.get("8"),
                    8937202,8939694));

            Assert.assertEquals("Zero coverage test fails:","N",pulledSequence);

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
//        String gvcfFileName = "/Users/zrm22/PanGenome/FASTAGeneration/B73_haplotype_caller_output.g.vcf";
        String gvcfFileName = fullGVCFFileName;
        String refFileName = agpv4FileName;
        try {
            GVCFGenomeSequence sequence = (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, gvcfFileName);
            GenomeSequence refSequence = (GenomeSequence)GenomeSequenceBuilder.instance(refFileName);

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
//        String gvcfFileName = "/Users/zrm22/PanGenome/FASTAGeneration/B73_haplotype_caller_output.g.vcf";
        String gvcfFileName = fullGVCFFileName;
        String refFileName = agpv4FileName;
        try {
            GVCFGenomeSequence sequence = (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, gvcfFileName);
            GenomeSequence refSequence = (GenomeSequence)GenomeSequenceBuilder.instance(refFileName);
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
//        String gvcfFileName = "/Users/zrm22/PanGenome/FASTAGeneration/B73_haplotype_caller_output.g.vcf";
        String gvcfFileName = fullGVCFFileName;
        String refFileName = agpv4FileName;
        try {
            GVCFGenomeSequence sequence = (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, gvcfFileName);
            System.out.println("Loading Reference sequence");
            GenomeSequence refSequence = (GenomeSequence)GenomeSequenceBuilder.instance(refFileName);
            System.out.println("Done loading Ref");
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

    //TODO fix for git upload
    @Ignore
    @Test
    public void pullRealV4Sequence() {
        String fileName = "/Users/zrm22/PanGenome/W22SentieonTest/W22_haplotype_caller_output.g.vcf";
        String refFileName = "/Users/zrm22/PanGenome/W22SentieonTest/Zea_mays.AGPv4.dna.toplevel.fa.gz";
        String outputFileName = "/Users/zrm22/PanGenome/W22SentieonTest/outputTest2.fa";
        try {
            GVCFGenomeSequence sequence =  (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, fileName);
            sequence.writeFASTA(outputFileName);
        }catch(Exception e) {
            e.printStackTrace();
        }
    }

    //TODO fix for git repo
    @Ignore
    @Test
    public void pullW22FASTASequence() {
        String refFileName = "/Users/zrm22/PanGenome/W22SentieonTest/Zea_mays.AGPv4.dna.toplevel.fa.gz";
        String outputFileName = "/Users/zrm22/PanGenome/FASTAGeneration/W22Chr8NoFilteringFullFileTry1.fa";
        String gvcfFileName = "/Users/zrm22/PanGenome/FASTAGeneration/W22_haplotype_caller_output.g.vcf";
        String intervalFileName = "/Users/zrm22/PanGenome/FASTAGeneration/anchorsFile_chr8_genesExactCoordinates.csv";

        try{
            System.out.println("Creating GVCFSequence object");
            GVCFGenomeSequence sequence = (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, gvcfFileName);
            HashMap<String, Chromosome> chrMap = new HashMap<>();


            Object[] chrs = sequence.chromosomes().toArray();
            for(int i = 0; i < chrs.length; i++) {
                chrMap.put(""+((Chromosome)chrs[i]).getChromosomeNumber(),(Chromosome)chrs[i]);
            }

            BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName));
            BufferedReader reader = new BufferedReader(new FileReader(intervalFileName));

            String currentLine = reader.readLine();
            int lineCounter = 0;

            while((currentLine = reader.readLine())!=null) {
                if(lineCounter%100==0) {
                    System.out.println("Counter: "+lineCounter);
                }
                lineCounter++;
                String[] currentLineSplit = currentLine.split(",");
                writer.write(">"+currentLineSplit[0]+":"+currentLineSplit[1]+":"+currentLineSplit[2]);
                writer.newLine();
                writer.write(""+NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chrMap.get(currentLineSplit[0]),
                        Integer.parseInt(currentLineSplit[1]),Integer.parseInt(currentLineSplit[2]))));
                writer.newLine();
            }
            writer.close();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }
    //TODO fix for git upload
    @Ignore
    @Test
    public void pullB73FASTASequence() {
        String refFileName = "/Users/zrm22/PanGenome/W22SentieonTest/Zea_mays.AGPv4.dna.toplevel.fa.gz";
        String outputFileName = "/Users/zrm22/PanGenome/FASTAGeneration/B73Chr8NoFilteringFullFileADTest6.fa";
        String outputFileName2 = "/Users/zrm22/PanGenome/FASTAGeneration/B73Chr8NoFilteringFullFileADTest6Ref.fa";

        String gvcfFileName = "/Users/zrm22/PanGenome/FASTAGeneration/B73_haplotype_caller_output.g.vcf";
        String intervalFileName = "/Users/zrm22/PanGenome/FASTAGeneration/anchorsFile_chr8_genesExactCoordinates.csv";

        try{
            System.out.println("Creating GVCFSequence object");
            GVCFGenomeSequence sequence = (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, gvcfFileName);

            GenomeSequence refSequence = (GenomeSequence)GenomeSequenceBuilder.instance(refFileName);
            HashMap<String, Chromosome> chrMap = new HashMap<>();

            Object[] chrs = sequence.chromosomes().toArray();
            for(int i = 0; i < chrs.length; i++) {
                chrMap.put(""+((Chromosome)chrs[i]).getChromosomeNumber(),(Chromosome)chrs[i]);
            }

            BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName));
            BufferedWriter writer2 = new BufferedWriter(new FileWriter(outputFileName2));
            BufferedReader reader = new BufferedReader(new FileReader(intervalFileName));

            String currentLine = reader.readLine();
            int lineCounter = 0;

            while((currentLine = reader.readLine())!=null) {
                if(lineCounter%10==0) {
                    System.out.println("Counter: "+lineCounter);
                }
                lineCounter++;
                String[] currentLineSplit = currentLine.split(",");
                String pulledSequence = ""+NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chrMap.get(currentLineSplit[0]),
                        Integer.parseInt(currentLineSplit[1]),Integer.parseInt(currentLineSplit[2])));

                writer.write(">" + currentLineSplit[0] + ":" + currentLineSplit[1] + ":" + currentLineSplit[2]);
                writer.newLine();
                writer.write(pulledSequence);
                String pulledSequenceRef = ""+NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get(currentLineSplit[0]),
                        Integer.parseInt(currentLineSplit[1]),Integer.parseInt(currentLineSplit[2])));
                writer.newLine();
                writer.write(pulledSequenceRef);
                writer.newLine();


                writer2.write(">" + currentLineSplit[0] + ":" + currentLineSplit[1] + ":" + currentLineSplit[2]);
                writer2.newLine();
                writer2.write(pulledSequenceRef);
                writer2.newLine();
                if(lineCounter==20) {
                    break;
                }
            }
            writer.close();
            writer2.close();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }

    //TODO fix for git upload
    @Ignore
    @Test
    public void pullB73FASTASequenceFullGVCF() {
        String refFileName = "/Users/zrm22/PanGenome/W22SentieonTest/Zea_mays.AGPv4.dna.toplevel.fa.gz";
        String outputFileName = "/Users/zrm22/PanGenome/FASTAGeneration/B73Chr8NoFilteringFullFileADMergedAnchors.fa";

        String gvcfFileName = "/Users/zrm22/PanGenome/FASTAGeneration/B73_haplotype_caller_output.g.vcf";
        String intervalFileName = "/Users/zrm22/PanGenome/FASTAGeneration/anchorsFile_chr8_MergedPlus1000orGapDiffcoordinates.csv";
        String statFileName = "/Users/zrm22/PanGenome/FASTAGeneration/statsFileMergedAnchors.csv";
        try{
            System.out.println("Creating GVCFSequence object");
            GVCFGenomeSequence sequence = (GVCFGenomeSequence)GVCFGenomeSequenceBuilder.instance(refFileName, gvcfFileName);
            System.out.println("Created Genome Sequence");
            HashMap<String, Chromosome> chrMap = new HashMap<>();

            Object[] chrs = sequence.chromosomes().toArray();
            for(int i = 0; i < chrs.length; i++) {
                chrMap.put(""+((Chromosome)chrs[i]).getChromosomeNumber(),(Chromosome)chrs[i]);
            }

            BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName));
            BufferedWriter statsWriter = new BufferedWriter(new FileWriter(statFileName));
            BufferedReader reader = new BufferedReader(new FileReader(intervalFileName));

            String currentLine = reader.readLine();
            int lineCounter = 0;
            statsWriter.write("ID,RequestedSize, ExportedSize,HetCount,AltCount,Depth,GQ,Min_Depth,ZeroCoverageCount,HomoRefCount,HomoAltLowDepthCount,HomoAltHighDepthCount");
            statsWriter.newLine();

            while((currentLine = reader.readLine())!=null) {
                if(lineCounter%10==0) {
                    System.out.println("Counter: "+lineCounter);
                }
                lineCounter++;
                String[] currentLineSplit = currentLine.split(",");
                String pulledSequence = ""+NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chrMap.get(currentLineSplit[0]),
                        Integer.parseInt(currentLineSplit[1]),Integer.parseInt(currentLineSplit[2])));

                writer.write(">" + currentLineSplit[0] + ":" + currentLineSplit[1] + ":" + currentLineSplit[2]);
                writer.newLine();
                writer.write(pulledSequence);
                writer.newLine();

                statsWriter.write(currentLineSplit[0] + ":" + currentLineSplit[1] + ":" + currentLineSplit[2]+",");
                HashMap<String, Integer> stats = sequence.getPreviousRegionStats();
                statsWriter.write(stats.get("RefSize")+",");
                statsWriter.write(stats.get("Size")+",");
                statsWriter.write(stats.get("HetCount")+",");
                statsWriter.write(stats.get("AltCount")+",");
                statsWriter.write(stats.get("Depth")+",");
                statsWriter.write(stats.get("GQ")+",");
                statsWriter.write(stats.get("Min_Depth")+",");
                statsWriter.write(stats.get("ZeroCoverageCount")+",");
                statsWriter.write(stats.get("HomoRefCount")+",");
                statsWriter.write(stats.get("HomoAltLowDepthCount")+",");
                statsWriter.write(""+stats.get("HomoAltHighDepthCount"));
                statsWriter.newLine();
                sequence.resetCounters();

            }
            writer.close();
            statsWriter.close();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }

    //TODO fix for git upload
    @Ignore
    @Test
    public void getChr8FromB73V4() {
        GenomeSequence sequence = GenomeSequenceBuilder.instance("/Users/zrm22/PanGenome/W22SentieonTest/Zea_mays.AGPv4.dna.toplevel.fa.gz");
        Object[] chrs = sequence.chromosomes().toArray();
        for(int i = 0; i < chrs.length; i++) {
            if(((Chromosome)chrs[i]).getChromosomeNumber()==8) {
                String sequenceString = ""+NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence((Chromosome)chrs[i],1,181122637));
                try {
                    BufferedWriter writer = new BufferedWriter(new FileWriter("/Users/zrm22/PanGenome/W22SentieonTest/CHR8_Zea_mays.AGPv4.dna.toplevel.fa"));

                    writer.write(">8 dna:chromosome chromosome:AGPv4:8:1:181122637:1 REF");
                    writer.newLine();
                    writer.write(sequenceString);
                    writer.newLine();

                    writer.close();
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    //TODO fix for git upload
    @Ignore
    @Test
    public void getChr8FromB73V3() {
        GenomeSequence sequence = GenomeSequenceBuilder.instance("/Users/zrm22/PanGenome/LynnsFiles/Zea_mays.AGPv3.29.dna.genome.fa.gz");
        Object[] chrs = sequence.chromosomes().toArray();
        for(int i = 0; i < chrs.length; i++) {
            if(((Chromosome)chrs[i]).getChromosomeNumber()==8) {
                String sequenceString = ""+NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence((Chromosome)chrs[i],1,175377492));
                try {
                    BufferedWriter writer = new BufferedWriter(new FileWriter("/Users/zrm22/PanGenome/W22SentieonTest/CHR8_Zea_mays.AGPv3.dna.toplevel.fa"));

                    writer.write(">8 dna:chromosome chromosome:AGPv3:8:1:175377492:1");
                    writer.newLine();
                    writer.write(sequenceString);
                    writer.newLine();

                    writer.close();
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    //TODO fix for git upload
    @Ignore
    @Test
    public void compareB73AnchorRegions() {
        String referenceHaplotypeFile = "/Users/zrm22/PanGenome/FASTAGeneration/ref_chr8_anchorsExactGenes.fa";
        String gatkHaplotypeFile = "/Users/zrm22/PanGenome/FASTAGeneration/B73Chr8NoFilteringFullFile.fa";
        try {
            BufferedReader referenceReader = new BufferedReader(new FileReader(referenceHaplotypeFile));
            BufferedReader gatkReader = new BufferedReader(new FileReader(gatkHaplotypeFile));


            HashMap<String, String> refTagToSequenceMap = new HashMap<>();
            HashMap<String, String> gatkTagToSequenceMap = new HashMap<>();

            String refString = "";
            while((refString = referenceReader.readLine())!=null) {
                //1-8:76622:81307
                String[] refTagSplitDash = refString.split("-");
                //String[] refTagSplitDashColon = refTagSplitDash[1].split(":");

                String tag = refTagSplitDash[1];

                String sequence = referenceReader.readLine();
                refTagToSequenceMap.put(tag,sequence);
            }

            String gatkString = "";
            while((gatkString=gatkReader.readLine())!=null) {
//                String[]
            }


        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }
    @Ignore
    @Test
    public void pullChr8GVCF() {
        try {
            BufferedReader reader = new BufferedReader(new FileReader("/Users/zrm22/PanGenome/FASTAGeneration/B73_haplotype_caller_output.g.vcf"));
            BufferedWriter writer = new BufferedWriter(new FileWriter("/Users/zrm22/PanGenome/FASTAGeneration/Chr8_B73_haplotype_caller_output.g.vcf"));

            String currLine = "";
            while((currLine = reader.readLine())!=null) {
                if(currLine.startsWith("#CHROM")) {
                    writer.write(currLine);
                    writer.newLine();
                    break;
                }
                else {
                    writer.write(currLine);
                    writer.newLine();
                }
            }
//            while((currLine = reader.readLine())!=null) {
//                if(currLine.startsWith("8"))
//            }
        }
        catch(Exception e) {
            e.printStackTrace();
        }
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

    @Test
    public void pullAnchor2kFromW22() {
        String gvcfFileName = "/Volumes/ZackBackup/Temp/Pangenome/DBTests/W22_haplotype_caller_output.g.vcf";
        String refFileName = "/Users/zrm22/ReferencePlus/Datasets/B73Ref/Zea_mays.AGPv4.dna.toplevel.fa";
        try {
            GVCFGenomeSequence sequence = (GVCFGenomeSequence) GVCFGenomeSequenceBuilder.instance(refFileName, gvcfFileName);
            System.out.println("Loading Reference sequence");
            GenomeSequence refSequence = (GenomeSequence) GenomeSequenceBuilder.instance(refFileName);
            System.out.println("Done loading Ref");
            HashMap<String, Chromosome> chrMap = new HashMap<>();


            Object[] chrs = sequence.chromosomes().toArray();
            for (int i = 0; i < chrs.length; i++) {
                chrMap.put("" + ((Chromosome) chrs[i]).getChromosomeNumber(), (Chromosome) chrs[i]);
            }
            String pulledSequence = "" + NucleotideAlignmentConstants.nucleotideBytetoString(sequence.chromosomeSequence(chrMap.get("1"),
                    37964752, 37965359));
//            37,964,752-37,965,359
            System.out.println(pulledSequence);

            String pulledSequenceRef = "" + NucleotideAlignmentConstants.nucleotideBytetoString(refSequence.chromosomeSequence(chrMap.get("1"),
                    37964752, 37965359));

            System.out.println(">B73Seq");
            System.out.println(pulledSequenceRef);
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }



}
