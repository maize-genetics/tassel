package net.maizegenetics.analysis.imputation;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.TableReport;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Random;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;

import static net.maizegenetics.dna.snp.NucleotideAlignmentConstants.getNucleotideDiploidByte;
import static net.maizegenetics.dna.snp.NucleotideAlignmentConstants.getNucleotideIUPAC;
import static org.junit.Assert.assertTrue;

public class LDKNNiImputationPluginTest {

    private static final Logger myLogger = LogManager.getLogger(FILLINImputationPluginTest.class);
    static String baseFile=GeneralConstants.TEMP_DIR+"LDKNNbase.hmp.txt.gz";
    static String maskedFile=GeneralConstants.TEMP_DIR+"LDKNNmasked.hmp.txt.gz";
    static String outputFile=GeneralConstants.TEMP_DIR+"LDjunk.hmp.txt.gz";
    private static String infile=TutorialConstants.HAPMAP_FILENAME;

    public static GenotypeTable inputAlign;

    public LDKNNiImputationPluginTest() {

    }

    @Test
    public void testBasicPlugin() throws Exception {
        clearTestFolder();
        System.out.println("Starting Mixed diverse and 10 line synthetic");
        createTestDataForImputation(100,600,10,true);
        imputeMissing(30, 10);

        System.out.println(testAccuracyPlugin().toStringTabDelim());

        TableReport tr=testAccuracyPlugin();
        assertTrue((Double) testAccuracyPlugin().getValueAt(7, 1) < 0.035);


        System.out.println("Starting Mixed diverse and 10 line RILS");
        createTestDataForImputation(10,600,10,false);
        imputeMissing(30, 10);

        System.out.println(testAccuracyPlugin().toStringTabDelim());
        assertTrue((Double) testAccuracyPlugin().getValueAt(7, 1) < 5e-4);

        System.out.println("281 diverse inbreds lines");
        createTestDataForImputation(281,0,0,false);
        imputeMissing(30, 10);

        System.out.println(testAccuracyPlugin().toStringTabDelim());
        assertTrue((Double) testAccuracyPlugin().getValueAt(7, 1) < 0.17);

    }


    private TableReport testAccuracyPlugin() throws Exception {
        DataSet maskedDataSet = new DataSet(Arrays.asList(
                new Datum("OriginalFile", ImportUtils.readGuessFormat(baseFile), ""),
                new Datum("MaskedFile", ImportUtils.readGuessFormat(maskedFile), ""),
                new Datum("ImputedFile", ImportUtils.readGuessFormat(outputFile), "")
        ),null);

        return new ImputationAccuracyPlugin()
                .runPlugin(maskedDataSet);
    }

    private static void clearTestFolder() {
        File f=new File(GeneralConstants.TEMP_DIR);
        System.out.println("Deleting files from "+GeneralConstants.TEMP_DIR);
        //System.out.println(f.getAbsolutePath());
        if(f.listFiles()!=null) {
            for (File file : f.listFiles()) {
                if(file.getName().contains("hmp.")) {
                    file.delete();
                    //System.out.println("Deleting:"+file.toString());
                }
            }
        }
    }

    public static void createTestDataForImputation(int numberOfOrig, int numberOfRILs, int numberOfRILparents,
                                                   boolean createHetChromosomes) {
        //Create masked dataset
        Random r=new Random(1);
        System.out.println("Reading: "+ infile);
        File f=new File(infile);
        System.out.println("Absolute Reading: "+f.getAbsolutePath());
        GenotypeTable origAlign = ImportUtils.readFromHapmap(infile, null);
        GenotypeTableBuilder ab=GenotypeTableBuilder.getTaxaIncremental(origAlign.positions());
        int origNumTaxa=origAlign.numberOfTaxa();
        int numberOfStates=(createHetChromosomes)?3:2;
        for (int i=0; i<numberOfOrig; i++) {
            ab.addTaxon(origAlign.taxa().get(i),origAlign.genotypeAllSites(i));
        }
        for (int nt = numberOfOrig /*set above zero to keep orig*/; nt < numberOfOrig+numberOfRILs; nt++) {
            int parent1=r.nextInt(numberOfRILparents);
            int parent2=r.nextInt(numberOfRILparents);
            Taxon tx=new Taxon(origAlign.taxaName(parent2)+origAlign.taxaName(parent1)+"_"+nt);
            int modeCnt;
            byte[] genotype=new byte[origAlign.numberOfSites()];
            for (Chromosome aL : origAlign.chromosomes()) {
//                System.out.printf("Chr:%s Mode:%d %n", aL.getName(), modeCnt);
                modeCnt=r.nextInt(numberOfStates);
                int[] startEnd=origAlign.positions().startAndEndOfChromosome(aL);
                for (int i = startEnd[0]; i <= startEnd[1]; i++) {
                   // if(r.nextInt(100)<1) modeCnt=r.nextInt(numberOfStates);
                    switch (modeCnt) {
                        case 0: genotype[i]=origAlign.genotype(parent2, i); break;
                        case 1: genotype[i]=origAlign.genotype(parent1, i); break;
                        case 2: genotype[i]= GenotypeTableUtils.getUnphasedDiploidValueNoHets(origAlign.genotype(parent2, i),
                                origAlign.genotype(parent1, i)); break;
                    }
                }
            }
            ab.addTaxon(tx,genotype);
        }
        GenotypeTable baseA=ab.build();
        //baseA=FilterGenotypeTableBuilder.getInstance(baseA).minorAlleleFreqForSite(0.01,.99).build();
        baseA=new FilterSiteBuilderPlugin().siteMinAlleleFreq(0.01).siteMaxAlleleFreq(0.99).runPlugin(baseA);
        System.out.println(baseFile);

        ExportUtils.writeToHapmap(baseA, false, baseFile, '\t', null);
        GenotypeCallTableBuilder gBMask=GenotypeCallTableBuilder.getInstance(baseA.numberOfTaxa(),baseA.numberOfSites());
        for (int t = 0; t < baseA.numberOfTaxa(); t++) {
            for (int s = 0; s < baseA.numberOfSites(); s++) {
                if(r.nextDouble()<0.1) {gBMask.setBase(t, s, GenotypeTable.UNKNOWN_DIPLOID_ALLELE);}
                else {gBMask.setBase(t, s, baseA.genotype(t,s));}
            }
        }
        GenotypeTable maskA=GenotypeTableBuilder.getInstance(gBMask.build(),baseA.positions(), baseA.taxa());
        ExportUtils.writeToHapmap(maskA, false, maskedFile, '\t', null);

    }

    private GenotypeTable imputeMissing(int maxHighLDSites, int maxKNNTaxa) {
        GenotypeTable hetGenotypeTable = ImportUtils.readGuessFormat(maskedFile);
        DataSet maskedDataSet = new DataSet(new Datum("MaskedFile", hetGenotypeTable,""),null);
       // DataSet maskedDataSet=new FileLoadPlugin(null, true).guessAtUnknowns(maskedFile);

        GenotypeTable impGenotype= new LDKNNiImputationPlugin()
            .highLDSSites(maxHighLDSites)
            .knnTaxa(maxKNNTaxa)
            //.maxDistance(300_000_000)
            .runPlugin(maskedDataSet)
        ;
        ExportUtils.writeToHapmap(impGenotype,outputFile);
        return impGenotype;
    }


    @Test
    public void testImpute() throws Exception {
        Multimap<Double,Byte> homoMap= ArrayListMultimap.create();
        homoMap.put(0.0, getNucleotideDiploidByte("A"));
        homoMap.put(0.0, getNucleotideDiploidByte("T"));
        homoMap.put(0.0, getNucleotideDiploidByte("T"));

        System.out.println(getNucleotideIUPAC(LDKNNiImputationPlugin.impute(homoMap, 20)));
       // assertEquals(getNucleotideDiploidByte("T"),LDKNNiImputationPlugin.impute(homoMap, 20));

        homoMap= ArrayListMultimap.create();
        homoMap.put(0.0, getNucleotideDiploidByte("A"));
        homoMap.put(0.0, getNucleotideDiploidByte("K"));
        homoMap.put(0.0, getNucleotideDiploidByte("K"));

        System.out.println(getNucleotideIUPAC(LDKNNiImputationPlugin.impute(homoMap, 20)));
     //   assertEquals(getNucleotideDiploidByte("K"),LDKNNiImputationPlugin.impute(homoMap, 20));

        homoMap= ArrayListMultimap.create();
        homoMap.put(0.0, getNucleotideDiploidByte("A"));
        homoMap.put(0.25, getNucleotideDiploidByte("T"));
        homoMap.put(0.5, getNucleotideDiploidByte("T"));
        homoMap.put(0.5, getNucleotideDiploidByte("T"));

        System.out.println(getNucleotideIUPAC(LDKNNiImputationPlugin.impute(homoMap, 20)));
    //    assertEquals(getNucleotideDiploidByte("A"),LDKNNiImputationPlugin.impute(homoMap, 20));

        homoMap= ArrayListMultimap.create();
        homoMap.put(0.0, getNucleotideDiploidByte("A"));
        homoMap.put(0.0, getNucleotideDiploidByte("T"));

        System.out.println(getNucleotideIUPAC(LDKNNiImputationPlugin.impute(homoMap, 20)));
  //      assertEquals(getNucleotideDiploidByte("W"),LDKNNiImputationPlugin.impute(homoMap, 20));

        homoMap= ArrayListMultimap.create();
        homoMap.put(0.0, getNucleotideDiploidByte("A"));
        homoMap.put(0.0, getNucleotideDiploidByte("T"));
        homoMap.put(0.0, getNucleotideDiploidByte("W"));

        System.out.println(getNucleotideIUPAC(LDKNNiImputationPlugin.impute(homoMap, 20)));
    //    assertEquals(getNucleotideDiploidByte("W"),LDKNNiImputationPlugin.impute(homoMap, 20));

    }
}