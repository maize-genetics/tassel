/*
 * ExportUtilsTest
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.score.AlleleDepth;
import net.maizegenetics.dna.snp.score.AlleleDepthBuilder;
import net.maizegenetics.dna.snp.score.AlleleDepthUtil;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.snp.io.BuilderFromHapMap;
import net.maizegenetics.dna.snp.io.BuilderFromVCF;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.dna.snp.io.VCFUtil;

import org.apache.commons.io.FileUtils;
import org.junit.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

/**
 *
 * @author terryc
 */
public class ExportUtilsTest {

    public static char DELIMIT_CHAR = '\t';
    public static String EXPORT_TEMP_DIR = GeneralConstants.TEMP_DIR + "ExportUtilsTest/";
    public static String HAPMAP_TEMP_FILENAME = "temp_hapmap.hmp.txt";

    public ExportUtilsTest() {
    }

    private static void clearTestFolder() {
        File f=new File(EXPORT_TEMP_DIR);
        System.out.println(f.getAbsolutePath());
        if(f.listFiles()!=null) {
            for (File file : f.listFiles()) {
                if(file.getName().contains("hmp.") || file.getName().contains(".h5") || file.getName().contains(".vcf")) {
                    file.delete();
                    System.out.println("Deleting:"+file.toString());
                }
                
            }
        }
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        File tempDir = new File(EXPORT_TEMP_DIR);
        if (!tempDir.exists()) {
            tempDir.mkdirs();
        }
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of writeToHapmap method, of class ExportUtils.
     */
    @Ignore
    @Test
    public void testWriteToHapmap() throws IOException {
        System.out.println("Testing Export to Hapmap...");
        String outputFile = EXPORT_TEMP_DIR + HAPMAP_TEMP_FILENAME;
        File input = new File(TutorialConstants.HAPMAP_FILENAME);
        System.out.println("   Input File: " + input.getCanonicalPath());
        if (!input.exists()) {
            fail("Input File: " + TutorialConstants.HAPMAP_FILENAME + " doesn't exist.");
        }
        GenotypeTable expResult = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
        ExportUtils.writeToHapmap(expResult, false, outputFile, DELIMIT_CHAR, null);
        GenotypeTable result = ImportUtils.readFromHapmap(outputFile, null);

        AlignmentTestingUtils.alignmentsEqual(expResult, result);
    }
    @Ignore
    @Test
    public void testHDF5WithDepth() throws Exception {
        clearTestFolder();
        GenotypeTable a=BuilderFromHapMap.getBuilder(TutorialConstants.HAPMAP_FILENAME).sortTaxa().build(); //need a sort taxa method to make the two approaches equal
        System.out.println("Testing roundtrip of HDF5 with Depth");
        GenotypeTable aWithDepth=AlignmentTestingUtils.createRandomDepthForGenotypeTable(a, 1);
        ExportUtils.writeGenotypeHDF5(aWithDepth,EXPORT_TEMP_DIR+"rtDepth.t5.h5",true);
        GenotypeTable rtWithDepth=ImportUtils.readGuessFormat(EXPORT_TEMP_DIR+"rtDepth.t5.h5");
        AlignmentTestingUtils.alignmentsEqual(a, rtWithDepth);

        for (int t=0; t<aWithDepth.numberOfTaxa(); t++) {
            for (int s=0; s<aWithDepth.numberOfSites(); s++) {
                int[] oD=aWithDepth.depthForAlleles(t,s);
                int[] x2D=rtWithDepth.depthForAlleles(t,s);
                boolean same=GenotypeTableUtils.isEqual(aWithDepth.genotype(t,s),rtWithDepth.genotype(t,s));
                Assert.assertTrue("Error"+aWithDepth.genotypeAsString(t,s)+" values"+Arrays.toString(oD),same);
                Assert.assertArrayEquals("Depth not  correct",oD,x2D);
            }
        }
    }
    
    @Test
    public void testVCFAllRef() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","A:A","T:T"};
        calls[1] = new String[]{"G:G","A:A","T:T"};
        calls[2] = new String[]{"G:G","A:A","T:T"};
        calls[3] = new String[]{"G:G","A:A","T:T"};
        calls[4] = new String[]{"G:G","A:A","T:T"};
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF file to the temp directory then load it and compare it
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFAllRef.vcf";
        loadAndCompareVCF(exportFileName,table);
    }
    
    @Test
    public void testVCFRefOneAlt() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A","T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T","G"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","A:A","G:G"};
        calls[1] = new String[]{"G:G","T:T","T:T"};
        calls[2] = new String[]{"G:G","A:A","G:G"};
        calls[3] = new String[]{"G:G","T:T","T:T"};
        calls[4] = new String[]{"G:G","A:A","T:T"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF file to the temp directory and load and compare
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAlt.vcf";
        loadAndCompareVCF(exportFileName,table);
    }
    @Test
    public void testVCFAllRefMissingTaxa() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","N:N","T:T"};
        calls[1] = new String[]{"G:G","A:A","T:T"};
        calls[2] = new String[]{"G:G","N:N","T:T"};
        calls[3] = new String[]{"G:G","A:A","N:N"};
        calls[4] = new String[]{"G:G","N:N","T:T"};
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF file to the temp directory
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFAllRefMissingTaxa.vcf";
        loadAndCompareVCF(exportFileName,table);
    }
    
    @Test
    public void testVCFRefOneAltMissingTaxa() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A","T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T","G"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 6; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[6][3];
        calls[0] = new String[]{"G:G","A:A","G:G"};
        calls[1] = new String[]{"G:G","T:T","T:T"};
        calls[2] = new String[]{"G:G","N:N","G:G"};
        calls[3] = new String[]{"G:G","T:T","T:T"};
        calls[4] = new String[]{"G:G","A:A","N:N"};
        calls[5] = new String[]{"G:G","A:A","T:T"};

        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF file to the temp directory and load and compare
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltMissingTaxa.vcf";
        loadAndCompareVCF(exportFileName,table);
    }
    
    @Test
    public void testVCFRefOneAltNoRefCalls() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A","T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T","G"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","T:T","G:G"};
        calls[1] = new String[]{"G:G","T:T","G:G"};
        calls[2] = new String[]{"G:G","T:T","G:G"};
        calls[3] = new String[]{"G:G","T:T","G:G"};
        calls[4] = new String[]{"G:G","T:T","G:G"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF file to the temp directory and load and compare
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCalls.vcf";
        loadAndCompareVCF(exportFileName,table);
    }
    
    @Test
    public void testVCFRefOneAltNoRefCallsMissingTaxa() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A","T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T","G"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","T:T","G:G"};
        calls[1] = new String[]{"G:G","T:T","N:N"};
        calls[2] = new String[]{"G:G","N:N","G:G"};
        calls[3] = new String[]{"G:G","T:T","N:N"};
        calls[4] = new String[]{"G:G","T:T","G:G"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF file to the temp directory and load and compare
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCallsMissingTaxa.vcf";
        loadAndCompareVCF(exportFileName,table);
    }
    
    
    
    //Indel Tests
    @Test
    public void testVCFRefOneAltIndel() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"NA","NT","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"NT","NG","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","A:A","G:G"};
        calls[1] = new String[]{"G:G","T:T","T:T"};
        calls[2] = new String[]{"G:G","A:A","-:-"};
        calls[3] = new String[]{"G:G","-:-","T:T"};
        calls[4] = new String[]{"G:G","A:A","T:T"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF file to the temp directory and load and compare
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltIndel.vcf";
        loadAndCompareVCF(exportFileName,table);
    }
    
    
    @Test
    public void testVCFRefOneAltMissingTaxaIndel() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"NA","NT","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"NT","NG","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","A:A","G:G"};
        calls[1] = new String[]{"G:G","T:T","T:T"};
        calls[2] = new String[]{"G:G","N:N","-:-"};
        calls[3] = new String[]{"G:G","-:-","T:T"};
        calls[4] = new String[]{"G:G","A:A","N:N"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF file to the temp directory and load and compare
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltMissingTaxaIndel.vcf";
        loadAndCompareVCF(exportFileName,table);
    }
    
    @Test
    public void testVCFRefOneAltNoRefCallsIndel() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            //.knownVariants(new String[]{"NA","NT","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            //.knownVariants(new String[]{"NT","NG","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","T:T","G:G"};
        calls[1] = new String[]{"G:G","T:T","-:-"};
        calls[2] = new String[]{"G:G","-:-","G:G"};
        calls[3] = new String[]{"G:G","T:T","-:-"};
        calls[4] = new String[]{"G:G","T:T","G:G"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF file to the temp directory and load and compare
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCallsIndel.vcf";
        loadAndCompareVCF(exportFileName,table);
    }
    
    @Test
    public void testVCFRefOneAltNoRefCallsMissingTaxaIndel() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"NA","NT","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"NT","NG"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","T:T","G:G"};
        calls[1] = new String[]{"G:G","-:-","N:N"};
        calls[2] = new String[]{"G:G","N:N","-:-"};
        calls[3] = new String[]{"G:G","-:-","N:N"};
        calls[4] = new String[]{"G:G","T:T","G:G"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF file to the temp directory and load and compare
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCallsMissingTaxaIndel.vcf";
        loadAndCompareVCF(exportFileName,table);
    }
    
    
    
    
    @Test
    public void testVCFRefOneAltQiTest() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            //.knownVariants(new String[]{"A","T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            //.knownVariants(new String[]{"T","G"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","A:A","G:G"};
        calls[1] = new String[]{"G:G","T:T","T:T"};
        calls[2] = new String[]{"G:G","A:A","G:G"};
        calls[3] = new String[]{"G:G","T:-","T:T"};
        calls[4] = new String[]{"G:G","A:A","T:T"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF file to the temp directory and load and compare
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAlt.vcf";
        loadAndCompareVCF(exportFileName,table);
    }
    
    //Create tests which export both a hapmap and vcf file
    //Load them in and compare Genotype Tables
    @Test
    public void testVCFvsHMPAllRef() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","A:A","T:T"};
        calls[1] = new String[]{"G:G","A:A","T:T"};
        calls[2] = new String[]{"G:G","A:A","T:T"};
        calls[3] = new String[]{"G:G","A:A","T:T"};
        calls[4] = new String[]{"G:G","A:A","T:T"};
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF and HMP files to the temp directory then load it and compare it
        String exportVCFFileName = EXPORT_TEMP_DIR+"/testVCFAllRef.vcf";
        String exportHMPFileName = EXPORT_TEMP_DIR+"/testHMPAllRef.hmp.txt";
        loadAndCompareVCFandHMP(exportVCFFileName,exportHMPFileName,table);
        
        //Test HMP->VCF
        String exportVCFFileName2 = EXPORT_TEMP_DIR+"/testVCFAllRef2.vcf";
        String exportHMPFileName2 = EXPORT_TEMP_DIR+"/testHMPAllRef2.hmp.txt";
        loadAndCompareHMPtoVCF(exportVCFFileName2,exportHMPFileName2,table);
    }
    
    @Test
    public void testVCFvsHMPRefOneAlt() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A","T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T","G"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","A:A","G:G"};
        calls[1] = new String[]{"G:G","T:T","T:T"};
        calls[2] = new String[]{"G:G","A:A","G:G"};
        calls[3] = new String[]{"G:G","T:T","T:T"};
        calls[4] = new String[]{"G:G","A:A","T:T"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        
        //Write out the VCF and HMP files to the temp directory then load it and compare it
        String exportVCFFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAlt.vcf";
        String exportHMPFileName = EXPORT_TEMP_DIR+"/testHMPRefOneAlt.hmp.txt";
        loadAndCompareVCFandHMP(exportVCFFileName,exportHMPFileName,table);
        
        //Test HMP->VCF
        String exportVCFFileName2 = EXPORT_TEMP_DIR+"/testVCFRefOneAlt2.vcf";
        String exportHMPFileName2 = EXPORT_TEMP_DIR+"/testHMPRefOneAlt2.hmp.txt";
        loadAndCompareHMPtoVCF(exportVCFFileName2,exportHMPFileName2,table);
    }
    
    @Test
    public void testVCFvsHMPAllRefMissingTaxa() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","N:N","T:T"};
        calls[1] = new String[]{"G:G","A:A","T:T"};
        calls[2] = new String[]{"G:G","N:N","T:T"};
        calls[3] = new String[]{"G:G","A:A","N:N"};
        calls[4] = new String[]{"G:G","N:N","T:T"};
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        
        //Write out the VCF and HMP files to the temp directory then load it and compare it
        String exportVCFFileName = EXPORT_TEMP_DIR+"/testVCFRefMissingTaxa.vcf";
        String exportHMPFileName = EXPORT_TEMP_DIR+"/testHMPRefMissingTaxa.hmp.txt";
        loadAndCompareVCFandHMP(exportVCFFileName,exportHMPFileName,table);
        
        //Test HMP->VCF
        String exportVCFFileName2 = EXPORT_TEMP_DIR+"/testVCFRefMissingTaxa2.vcf";
        String exportHMPFileName2 = EXPORT_TEMP_DIR+"/testHMPRefMissingTaxa2.hmp.txt";
        loadAndCompareHMPtoVCF(exportVCFFileName2,exportHMPFileName2,table);
        
    }
    
    @Test
    public void testVCFvsHMPRefOneAltMissingTaxa() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A","T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T","G"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 6; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[6][3];
        calls[0] = new String[]{"G:G","A:A","G:G"};
        calls[1] = new String[]{"G:G","T:T","T:T"};
        calls[2] = new String[]{"G:G","N:N","G:G"};
        calls[3] = new String[]{"G:G","T:T","T:T"};
        calls[4] = new String[]{"G:G","A:A","N:N"};
        calls[5] = new String[]{"G:G","A:A","T:T"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        
        //Write out the VCF and HMP files to the temp directory then load it and compare it
        String exportVCFFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltMissingTaxa.vcf";
        String exportHMPFileName = EXPORT_TEMP_DIR+"/testHMPRefOneAltMissingTaxa.hmp.txt";
        loadAndCompareVCFandHMP(exportVCFFileName,exportHMPFileName,table);
        
        //Test HMP->VCF
        String exportVCFFileName2 = EXPORT_TEMP_DIR+"/testVCFRefOneAltMissingTaxa2.vcf";
        String exportHMPFileName2 = EXPORT_TEMP_DIR+"/testHMPRefOneAltMissingTaxa2.hmp.txt";
        loadAndCompareHMPtoVCF(exportVCFFileName2,exportHMPFileName2,table);
    }
    
    @Test
    public void testVCFvsHMPRefOneAltNoRefCalls() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A","T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T","G"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","T:T","G:G"};
        calls[1] = new String[]{"G:G","T:T","G:G"};
        calls[2] = new String[]{"G:G","T:T","G:G"};
        calls[3] = new String[]{"G:G","T:T","G:G"};
        calls[4] = new String[]{"G:G","T:T","G:G"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF and HMP files to the temp directory then load it and compare it
        String exportVCFFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCalls.vcf";
        String exportHMPFileName = EXPORT_TEMP_DIR+"/testHMPRefOneAltNoRefCalls.hmp.txt";
        loadAndCompareVCFandHMP(exportVCFFileName,exportHMPFileName,table);
        
        //Test HMP->VCF
        String exportVCFFileName2 = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCalls2.vcf";
        String exportHMPFileName2 = EXPORT_TEMP_DIR+"/testHMPRefOneAltNoRefCalls2.hmp.txt";
        loadAndCompareHMPtoVCF(exportVCFFileName2,exportHMPFileName2,table);
    }
    
    @Test
    public void testVCFvsHMPRefOneAltNoRefCallsMissingTaxa() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A","T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T","G"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","T:T","G:G"};
        calls[1] = new String[]{"G:G","T:T","N:N"};
        calls[2] = new String[]{"G:G","N:N","G:G"};
        calls[3] = new String[]{"G:G","T:T","N:N"};
        calls[4] = new String[]{"G:G","T:T","G:G"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        
        //Write out the VCF and HMP files to the temp directory then load it and compare it
        String exportVCFFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCallsMissingTaxa.vcf";
        String exportHMPFileName = EXPORT_TEMP_DIR+"/testHMPRefOneAltNoRefCallsMissingTaxa.hmp.txt";
        loadAndCompareVCFandHMP(exportVCFFileName,exportHMPFileName,table);
        
        //Test HMP->VCF
        String exportVCFFileName2 = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCallsMissingTaxa2.vcf";
        String exportHMPFileName2 = EXPORT_TEMP_DIR+"/testHMPRefOneAltNoRefCallsMissingTaxa2.hmp.txt";
        loadAndCompareHMPtoVCF(exportVCFFileName2,exportHMPFileName2,table);
    }
    
    //Indel Tests
    @Test
    public void testVCFvsHMPRefOneAltIndel() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"NA","NT","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"NT","NG","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","A:A","G:G"};
        calls[1] = new String[]{"G:G","T:T","T:T"};
        calls[2] = new String[]{"G:G","A:A","-:-"};
        calls[3] = new String[]{"G:G","-:-","T:T"};
        calls[4] = new String[]{"G:G","A:A","T:T"};
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF and HMP files to the temp directory then load it and compare it
        String exportVCFFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltIndel.vcf";
        String exportHMPFileName = EXPORT_TEMP_DIR+"/testHMPRefOneAltIndel.hmp.txt";
        loadAndCompareVCFandHMP(exportVCFFileName,exportHMPFileName,table);
        
        //Test HMP->VCF
        String exportVCFFileName2 = EXPORT_TEMP_DIR+"/testVCFRefOneAltIndel2.vcf";
        String exportHMPFileName2 = EXPORT_TEMP_DIR+"/testHMPRefOneAltIndel2.hmp.txt";
        loadAndCompareHMPtoVCF(exportVCFFileName2,exportHMPFileName2,table);
    }
    
    @Test
    public void testVCFvsHMPRefOneAltMissingTaxaIndel() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"NA","NT","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"NT","NG","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","A:A","G:G"};
        calls[1] = new String[]{"G:G","T:T","T:T"};
        calls[2] = new String[]{"G:G","N:N","-:-"};
        calls[3] = new String[]{"G:G","-:-","T:T"};
        calls[4] = new String[]{"G:G","A:A","N:N"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF and HMP files to the temp directory then load it and compare it
        String exportVCFFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltMissingTaxaIndel.vcf";
        String exportHMPFileName = EXPORT_TEMP_DIR+"/testHMPRefOneAltMissingTaxaIndel.hmp.txt";
        loadAndCompareVCFandHMP(exportVCFFileName,exportHMPFileName,table);
        
        //Test HMP->VCF 
        String exportVCFFileName2 = EXPORT_TEMP_DIR+"/testVCFRefOneAltMissingTaxaIndel2.vcf";
        String exportHMPFileName2 = EXPORT_TEMP_DIR+"/testHMPRefOneAltMissingTaxaIndel2.hmp.txt";
        loadAndCompareHMPtoVCF(exportVCFFileName2,exportHMPFileName2,table);
    }
    
    @Test
    public void testVCFvsHMPRefOneAltNoRefCallsIndel() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"NA","NT","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"NT","NG","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","T:T","G:G"};
        calls[1] = new String[]{"G:G","T:T","-:-"};
        calls[2] = new String[]{"G:G","-:-","G:G"};
        calls[3] = new String[]{"G:G","T:T","-:-"};
        calls[4] = new String[]{"G:G","T:T","G:G"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF and HMP files to the temp directory then load it and compare it
        String exportVCFFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCallsIndel.vcf";
        String exportHMPFileName = EXPORT_TEMP_DIR+"/testHMPRefOneAltNoRefCallsIndel.hmp.txt";
        loadAndCompareVCFandHMP(exportVCFFileName,exportHMPFileName,table);
        
        //Test the HMP->VCF interaction
        String exportVCFFileName2 = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCallsIndel2.vcf";
        String exportHMPFileName2 = EXPORT_TEMP_DIR+"/testHMPRefOneAltNoRefCallsIndel2.hmp.txt";
        loadAndCompareHMPtoVCF(exportVCFFileName2,exportHMPFileName2,table);
    }
    
    @Test
    public void testVCFvsHMPRefOneAltNoRefCallsMissingTaxaIndel() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"NA","NT","N"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"NT","NG"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 5; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[5][3];
        calls[0] = new String[]{"G:G","T:T","G:G"};
        calls[1] = new String[]{"G:G","-:-","N:N"};
        calls[2] = new String[]{"G:G","N:N","-:-"};
        calls[3] = new String[]{"G:G","-:-","N:N"};
        calls[4] = new String[]{"G:G","T:T","G:G"};
        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl);
        
        //Write out the VCF and HMP files to the temp directory then load it and compare it
        String exportVCFFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCallsMissingTaxaIndel.vcf";
        String exportHMPFileName = EXPORT_TEMP_DIR+"/testHMPRefOneAltNoRefCallsMissingTaxaIndel.hmp.txt";
        loadAndCompareVCFandHMP(exportVCFFileName,exportHMPFileName,table);
        
        //Test the HMP->VCF interaction
        String exportVCFFileName2 = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCallsMissingTaxaIndel2.vcf";
        String exportHMPFileName2 = EXPORT_TEMP_DIR+"/testHMPRefOneAltNoRefCallsMissingTaxaIndel2.hmp.txt";
        loadAndCompareHMPtoVCF(exportVCFFileName2,exportHMPFileName2,table);
    }
    
    @Test
    public void testPLWithSingleAlleleSites() {
        clearTestFolder();
        
        String fileName = GeneralConstants.DATA_DIR+"CandidateTests/VCFFiles/testPLWithSingleAlleleSites.vcf";
        BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName).keepDepth();
        
        //Build GenotypeTable
        GenotypeTable table = builder.buildAndSortInMemory();
        
        //Export GenotypeTable to VCF
        String outputFileName = GeneralConstants.TEMP_DIR + "testPLWithSingleAlleleSitesExport.vcf";
        try {
            ExportUtils.writeToVCF(table, outputFileName, true);   
            
        }
        catch(Exception e) {
            Assert.fail("Error exporting due to PL calcuation with 1 allele(no ALTS)");
        }
       
        //Read in Exported and compare to a Input VCF file
        try {
            boolean filesMatch = FileUtils.contentEquals(new File(GeneralConstants.DATA_DIR+"ExpectedResults/VCFFiles/testPLWithSingleAlleleSitesExport.vcf"), new File(outputFileName));              
            assertTrue("Exported File does not exist.  Incorrect coding of PL likely an ArrayOutOfBounds",filesMatch);
        }
        catch(Exception e) {
            System.out.println(e);
        }
    }
    
    @Test
    public void testDPImported() {
        clearTestFolder();
        
        String fileName = GeneralConstants.DATA_DIR+"CandidateTests/VCFFiles/DPInfoImport.vcf";
        BuilderFromVCF builder = BuilderFromVCF.getBuilder(fileName).keepDepth();
        
        //Build GenotypeTable
        GenotypeTable table = builder.buildAndSortInMemory();
        
        //Export GenotypeTable to VCF
        String outputFileName = GeneralConstants.TEMP_DIR + "DPInfoImportExport.vcf";
        try {
            ExportUtils.writeToVCF(table, outputFileName, true);   
            
        }
        catch(Exception e) {
            Assert.fail("Error exporting due to PL calcuation with 1 allele(no ALTS)");
        }
       
        //Read in Exported and compare to a Input VCF file
        try {
            boolean filesMatch = FileUtils.contentEquals(new File(GeneralConstants.DATA_DIR+"CandidateTests/VCFFiles/DPInfoImport.vcf"), new File(outputFileName));              
            assertTrue("Exported File does not exist.  Incorrect coding of PL likely an ArrayOutOfBounds",filesMatch);
        }
        catch(Exception e) {
            System.out.println(e);
        }
    }
    
    //Here we are trying to cause the VCF export to reorder the 
    //alleles(by changing the order of the knownVariants or ignoring them) and make sure the depth is consistent
    @Test
    public void testVCFRefOneAltNoRefCallsDepthIssue() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             //.knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"T","A"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T","G"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 3; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[3][3];
        calls[0] = new String[]{"A:A","T:T","G:G"};
        calls[1] = new String[]{"G:A","A:A","G:G"};
        calls[2] = new String[]{"A:A","T:T","T:T"};
   
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        AlleleDepthBuilder adBuilder = AlleleDepthBuilder.getInstance(3, 3,tl);
        int[][][] values = new int[3][6][3];
   
        values[0][NucleotideAlignmentConstants.A_ALLELE][0] = 51;
        values[0][NucleotideAlignmentConstants.G_ALLELE][0] = 1;
        
        values[0][NucleotideAlignmentConstants.A_ALLELE][1] = 0;
        values[0][NucleotideAlignmentConstants.T_ALLELE][1] = 51;
        
        values[0][NucleotideAlignmentConstants.G_ALLELE][2] = 30;
        values[0][NucleotideAlignmentConstants.T_ALLELE][2] = 0;
        
        values[1][NucleotideAlignmentConstants.A_ALLELE][0] = 52;
        values[1][NucleotideAlignmentConstants.G_ALLELE][0] = 30;
        
        values[1][NucleotideAlignmentConstants.A_ALLELE][1] = 53;
        values[1][NucleotideAlignmentConstants.T_ALLELE][1] = 0;
        
        values[1][NucleotideAlignmentConstants.G_ALLELE][2] = 30;
        values[1][NucleotideAlignmentConstants.T_ALLELE][2] = 1;
        
        values[2][NucleotideAlignmentConstants.A_ALLELE][0] = 53;
        values[2][NucleotideAlignmentConstants.G_ALLELE][0] = 0;
        
        values[2][NucleotideAlignmentConstants.A_ALLELE][1] = 0;
        values[2][NucleotideAlignmentConstants.T_ALLELE][1] = 53;
        
        values[2][NucleotideAlignmentConstants.G_ALLELE][2] = 1;
        values[2][NucleotideAlignmentConstants.T_ALLELE][2] = 53;
    
        for(int taxonCounter = 0; taxonCounter < values.length; taxonCounter++) {
            adBuilder.addTaxon(taxonCounter, AlleleDepthUtil.depthIntToByte(values[taxonCounter]));
        }
       
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl,adBuilder.build());
        
        //Write out the VCF file to the temp directory and load and compare
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCallsDepth.vcf";
        loadAndCompareVCFDepth(exportFileName,table);
    }
    
    //This test is designed for making sure we are not converting the ints to bytes 
    //then getting the int depth from those bytes in the VCF export.
    //Before the depth 141 would be casted to a byte and then the depth would be set to 4190
    //Also Before depth 285 would be cast to a byte and the AlleleDepth would return 29
    @Test
    public void testVCFRefOneAltNoRefCallsByteDepthConversionIssue() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A","T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T","G"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 3; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[3][3];
        calls[0] = new String[]{"A:G","T:A","G:G"};
        calls[1] = new String[]{"G:A","A:A","G:T"};
        calls[2] = new String[]{"A:G","T:T","T:T"};
//        calls[3] = new String[]{"G:G","T:T","G:G"};
//        calls[4] = new String[]{"G:G","T:T","G:G"};
//        
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        AlleleDepthBuilder adBuilder = AlleleDepthBuilder.getInstance(3, 3,tl);
        int[][][] values = new int[3][6][3];
   
        values[0][NucleotideAlignmentConstants.A_ALLELE][0] = 28;
        values[0][NucleotideAlignmentConstants.G_ALLELE][0] = 141;
        
        values[0][NucleotideAlignmentConstants.A_ALLELE][1] = 285;
        values[0][NucleotideAlignmentConstants.T_ALLELE][1] = 51;
        
        values[0][NucleotideAlignmentConstants.G_ALLELE][2] = 30;
        values[0][NucleotideAlignmentConstants.T_ALLELE][2] = 1;
        
        values[1][NucleotideAlignmentConstants.A_ALLELE][0] = 52;
        values[1][NucleotideAlignmentConstants.G_ALLELE][0] = 30;
        
        values[1][NucleotideAlignmentConstants.A_ALLELE][1] = 30;
        values[1][NucleotideAlignmentConstants.T_ALLELE][1] = 0;
        
        values[1][NucleotideAlignmentConstants.G_ALLELE][2] = 30;
        values[1][NucleotideAlignmentConstants.T_ALLELE][2] = 52;
        
        values[2][NucleotideAlignmentConstants.A_ALLELE][0] = 53;
        values[2][NucleotideAlignmentConstants.G_ALLELE][0] = 31;
        
        values[2][NucleotideAlignmentConstants.A_ALLELE][1] = 1;
        values[2][NucleotideAlignmentConstants.T_ALLELE][1] = 53;
        
        values[2][NucleotideAlignmentConstants.G_ALLELE][2] = 2;
        values[2][NucleotideAlignmentConstants.T_ALLELE][2] = 53;
    
        for(int taxonCounter = 0; taxonCounter < values.length; taxonCounter++) {
            adBuilder.addTaxon(taxonCounter, AlleleDepthUtil.depthIntToByte(values[taxonCounter]));
        }
       
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl,adBuilder.build());
        
        //Write out the VCF file to the temp directory and load and compare
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCallsDepthConversion.vcf";
        loadAndCompareVCFDepth(exportFileName,table);
    }
    
    //This test is designed for making sure we are not converting the ints to bytes 
    //then getting the int depth from those bytes in the VCF export like the previous test.
    //Before the depth 141 would be casted to a byte and then the depth would be set to 4190
    //Also Before depth 285 would be cast to a byte and the AlleleDepth would return 29
    
    @Test
    public void testVCFRefOneAltNoRefCallsByteDepthConversionExtraDepth() throws Exception {
        clearTestFolder();
        //Create PositionList
        Chromosome chr=new Chromosome("1");
        Position pos1= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A","T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T","G"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions = new ArrayList<Position>();
        positions.add(pos1);
        positions.add(pos2);
        positions.add(pos3);
        PositionList pl = PositionListBuilder.getInstance(positions);
        
        //Create Taxa List
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < 3; i++) {
           Taxon at= new Taxon.Builder("Taxon_"+i+"")
                .build();
            tlb.add(at);
        }
        TaxaList tl=tlb.build();

        
        //Create call tables
        String[][] calls = new String[3][3];
        calls[0] = new String[]{"A:G","T:A","G:G"};
        calls[1] = new String[]{"G:A","A:A","G:T"};
        calls[2] = new String[]{"A:G","T:T","T:T"};
      
        
        GenotypeCallTable callTable = generateSmallCallTable(calls);
        
        AlleleDepthBuilder adBuilder = AlleleDepthBuilder.getInstance(3, 3,tl);
        int[][][] values = new int[3][6][3];
   
        values[0][NucleotideAlignmentConstants.A_ALLELE][0] = 28;
        values[0][NucleotideAlignmentConstants.G_ALLELE][0] = 141;
        
        values[0][NucleotideAlignmentConstants.A_ALLELE][1] = 285;
        values[0][NucleotideAlignmentConstants.T_ALLELE][1] = 51;
        
        values[0][NucleotideAlignmentConstants.G_ALLELE][2] = 30;
        values[0][NucleotideAlignmentConstants.T_ALLELE][2] = 1;
        
        values[1][NucleotideAlignmentConstants.A_ALLELE][0] = 52;
        values[1][NucleotideAlignmentConstants.G_ALLELE][0] = 30;
        
        values[1][NucleotideAlignmentConstants.A_ALLELE][1] = 30;
        values[1][NucleotideAlignmentConstants.T_ALLELE][1] = 0;
        
        values[1][NucleotideAlignmentConstants.G_ALLELE][2] = 30;
        values[1][NucleotideAlignmentConstants.T_ALLELE][2] = 52;
        
        values[2][NucleotideAlignmentConstants.A_ALLELE][0] = 53;
        values[2][NucleotideAlignmentConstants.G_ALLELE][0] = 31;
        
        values[2][NucleotideAlignmentConstants.A_ALLELE][1] = 1;
        values[2][NucleotideAlignmentConstants.T_ALLELE][1] = 53;
        
        values[2][NucleotideAlignmentConstants.G_ALLELE][2] = 2;
        values[2][NucleotideAlignmentConstants.T_ALLELE][2] = 53;
    
        for(int taxonCounter = 0; taxonCounter < values.length; taxonCounter++) {
            adBuilder.addTaxon(taxonCounter, AlleleDepthUtil.depthIntToByte(values[taxonCounter]));
        }
        AlleleDepth ad = adBuilder.build();
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList 
        GenotypeTable table = GenotypeTableBuilder.getInstance(callTable, pl, tl,ad);
        
        //Write out the VCF file to the temp directory and load and compare
        String exportFileName = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCallsDepthConversionExtraDepth.vcf";
        
//        loadAndCompareVCFSameTableDifferentDepth(exportFileName,table);
        loadAndCompareVCFDepth(exportFileName,table);
        
        //Now if we make a new PositionList with the correct Variants it should pass the full CompareVCF assertions
        //Create PositionList
        Position pos1_2= new GeneralPosition.Builder(chr,5)
             .knownVariants(new String[]{"G","C"}) //two SNPs and an indel
             .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.G_ALLELE)
             .build();
        
        Position pos2_2 = new GeneralPosition.Builder(chr,10)
                            .knownVariants(new String[]{"A","T"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.A_ALLELE)
                            .build();
        Position pos3_2 = new GeneralPosition.Builder(chr,15)
                            .knownVariants(new String[]{"T","G"})
                            .allele(WHICH_ALLELE.Reference, NucleotideAlignmentConstants.T_ALLELE)
                            .build();
        ArrayList<Position> positions_2 = new ArrayList<Position>();
        positions_2.add(pos1_2);
        positions_2.add(pos2_2);
        positions_2.add(pos3_2);
        PositionList pl_2 = PositionListBuilder.getInstance(positions_2);
        //Create a basic Genotype table with the signature (GenotypeCallTable, PositionList, TaxaList, AlleleDepth)
        GenotypeTable table_2 = GenotypeTableBuilder.getInstance(callTable, pl_2, tl,ad);
        String exportFileName_2 = EXPORT_TEMP_DIR+"/testVCFRefOneAltNoRefCallsDepthConversionExtraDepth_2.vcf";
        
        loadAndCompareVCFDepth(exportFileName_2,table_2);
    }
    
    @Test
    public void testDepths() {
        for(int i = -1; i < 1024; i++) {
            if(AlleleDepthUtil.depthByteToInt((byte)i)==128 ) {
                System.out.println(""+i);
            }
        }
    }
    
    @Ignore
    @Test
    public void find4190Depth() {
        for(int i = 0; i<1000000;i++) {
            if(AlleleDepthUtil.depthByteToInt((byte)i)==4190) {
                System.out.println(i);
            }
        }
        
    }
    @Ignore
    @Test
    public void find29Depth() {
        for(int i = 0; i<500;i++) {
            if(AlleleDepthUtil.depthByteToInt((byte)i)==29) {
                System.out.println(i);
            }
        }
        
    }
    
    @Ignore
    @Test
    public void calcScores() {
        clearTestFolder();
        
        try {
            BufferedWriter homo0HeatMapWriter = new BufferedWriter(new FileWriter(EXPORT_TEMP_DIR+"homo0HM.csv"));
            BufferedWriter hetHeatMapWriter = new BufferedWriter(new FileWriter(EXPORT_TEMP_DIR+"het01HM.csv"));
            BufferedWriter homo1HeatMapWriter = new BufferedWriter(new FileWriter(EXPORT_TEMP_DIR+"homo1HM.csv"));
            
            homo0HeatMapWriter.write(", ");
            hetHeatMapWriter.write(", ");
            homo1HeatMapWriter.write(", ");
            
            for(int i = 0; i<256; i++) {
                homo0HeatMapWriter.write(i+", ");
                hetHeatMapWriter.write(i+", ");
                homo1HeatMapWriter.write(i+", ");
            }
            
            homo0HeatMapWriter.newLine();
            hetHeatMapWriter.newLine();
            homo1HeatMapWriter.newLine();
            
            for(int i = 0; i < 256; i++) {
                homo0HeatMapWriter.write(i+", ");
                hetHeatMapWriter.write(i+", ");
                homo1HeatMapWriter.write(i+", ");
                for(int j = 0; j < 256; j++) {
                    int[] scores = VCFUtil.getScore(i,j);
                    homo0HeatMapWriter.write(scores[0]+", ");
                    hetHeatMapWriter.write(scores[1]+", ");
                    homo1HeatMapWriter.write(scores[2]+", ");
                }
                homo0HeatMapWriter.newLine();
                hetHeatMapWriter.newLine();
                homo1HeatMapWriter.newLine();
              
            }
            homo0HeatMapWriter.close();
            hetHeatMapWriter.close();
            homo1HeatMapWriter.close();
        }catch(IOException e) {
            System.out.println(e);
        }
        
        System.out.println((byte)4190);
        System.out.println((byte)-115);
        System.out.println(AlleleDepthUtil.depthIntToByte(4190));
        
        System.out.println(AlleleDepthUtil.depthByteToInt((byte)-115));
        
        System.out.println(Arrays.toString(VCFUtil.getScore(29, 51)));
        System.out.println(Arrays.toString(VCFUtil.getScore(500, 50)));
        System.out.println(Arrays.toString(VCFUtil.getScore(28,AlleleDepthUtil.depthByteToInt((byte)4190))));
        System.out.println(Arrays.toString(VCFUtil.getScore(28,AlleleDepthUtil.depthByteToInt((byte)338))));
        
        System.out.println(AlleleDepthUtil.depthIntToByte(4190));
    }
    
    public GenotypeCallTable generateSmallCallTable(String[][] calls) {
        GenotypeCallTableBuilder callTableBuilder = GenotypeCallTableBuilder.getInstance(calls.length,calls[0].length);
        callTableBuilder.setBases(calls);
        return callTableBuilder.build();
    }
    
    public void loadAndCompareVCF(String exportFileName, GenotypeTable table) {
        ExportUtils.writeToVCF(table, exportFileName , false);
        //Load up the VCF and compare it to the manually created table 
        GenotypeTable loadedGenotypeTable=ImportUtils.readGuessFormat(exportFileName);
        AlignmentTestingUtils.alignmentsEqual(table, loadedGenotypeTable);
    }
    public void loadAndCompareVCFDepth(String exportFileName, GenotypeTable table) {
        ExportUtils.writeToVCF(table, exportFileName , true);
        //Load up the VCF and compare it to the manually created table 
        GenotypeTable loadedGenotypeTable=ImportUtils.readGuessFormat(exportFileName);
        AlignmentTestingUtils.alignmentsEqual(table, loadedGenotypeTable);
        TaxaList loadedGenotypeTaxaList = loadedGenotypeTable.taxa();
        TaxaList generatedGenotypeTaxaList = table.taxa();
        
        System.out.println("Testing depths");
        PositionList loadedGenotypePosList = loadedGenotypeTable.positions();
        PositionList generatedGenotypePosList = table.positions();
        for(int taxon = 0; taxon< loadedGenotypeTaxaList.size(); taxon++) {
            for(int pos = 0; pos < loadedGenotypePosList.size(); pos++) {
                //Compare Depths
                int[] loadedDepths = loadedGenotypeTable.depthForAlleles(taxon, pos);
                int[] generatedDepths = table.depthForAlleles(taxon, pos);
                Assert.assertArrayEquals("Depths do not match for Taxon:"+taxon+" Pos:"+pos+".",loadedDepths,generatedDepths);
                
            }
        }
        
    }
    //Auxilary method to check that the GenotypeTables are the same but the depths are different
    public void loadAndCompareVCFSameTableDifferentDepth(String exportFileName, GenotypeTable table) {
        ExportUtils.writeToVCF(table, exportFileName , true);
        //Load up the VCF and compare it to the manually created table 
        GenotypeTable loadedGenotypeTable=ImportUtils.readGuessFormat(exportFileName);
        AlignmentTestingUtils.alignmentsEqual(table, loadedGenotypeTable);
        TaxaList loadedGenotypeTaxaList = loadedGenotypeTable.taxa();
        TaxaList generatedGenotypeTaxaList = table.taxa();
        
        System.out.println("Testing depths");
        PositionList loadedGenotypePosList = loadedGenotypeTable.positions();
        PositionList generatedGenotypePosList = table.positions();
        boolean differingDepths = false;
        ArrayList<Integer> positionHolder = new ArrayList<Integer>();
        for(int taxon = 0; taxon< loadedGenotypeTaxaList.size(); taxon++) {
            for(int pos = 0; pos < loadedGenotypePosList.size(); pos++) {
                //Compare Depths
                int[] loadedDepths = loadedGenotypeTable.depthForAlleles(taxon, pos);
                int[] generatedDepths = table.depthForAlleles(taxon, pos);
                for(int depthIndex = 0; depthIndex < loadedDepths.length; depthIndex++) {
                    if(loadedDepths[depthIndex] != generatedDepths[depthIndex]) {
                        differingDepths = true;
                        positionHolder.add(pos);
                    }
                }
            }
        }
        Assert.assertTrue("Depths are incorrectly matching at positions:"+positionHolder.toString(), differingDepths);    
    }
    
    
    public void loadAndCompareVCFandHMP(String exportVCFFileName, String exportHMPFileName,GenotypeTable table) {
        ExportUtils.writeToVCF(table,exportVCFFileName, false);
        ExportUtils.writeToHapmap(table, exportHMPFileName);
        
        GenotypeTable loadedVCF = ImportUtils.readGuessFormat(exportVCFFileName);
        GenotypeTable loadedHMP = ImportUtils.readGuessFormat(exportHMPFileName);
        
        AlignmentTestingUtils.alignmentsEqual(table, loadedVCF);
        AlignmentTestingUtils.alignmentsEqual(table, loadedHMP);
        AlignmentTestingUtils.alignmentsEqual(loadedHMP, loadedVCF);
        
    }
    
    public void loadAndCompareHMPtoVCF(String exportVCFFileName, String exportHMPFileName,GenotypeTable table) {
        //Export as HMP
        ExportUtils.writeToHapmap(table, exportHMPFileName);
        
        //Load in the HMP
        GenotypeTable loadedHMP = ImportUtils.readGuessFormat(exportHMPFileName);
        
        //Export the loadedHMP to VCF
        ExportUtils.writeToVCF(table,exportVCFFileName, false);
        
        //Load in the VCF
        GenotypeTable loadedVCF = ImportUtils.readGuessFormat(exportVCFFileName);
        
        AlignmentTestingUtils.alignmentsEqual(table, loadedVCF);
        AlignmentTestingUtils.alignmentsEqual(table, loadedHMP);
        AlignmentTestingUtils.alignmentsEqual(loadedHMP, loadedVCF);
    }
    
    
}