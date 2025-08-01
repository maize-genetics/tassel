package net.maizegenetics.taxa;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Multimap;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Set;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class TaxaListIOUtilsTest {
    static String taxaAnnotationFile=GeneralConstants.DATA_DIR+"CandidateTests/AnnotationForGoodman282.txt";
    static String gbsTestingKeyFile=GBSConstants.GBS_TESTING_KEY_FILE;
    static String testMutFile = GeneralConstants.TEMP_DIR + "test.hmp.h5";
    private static GenotypeTable inputAlign;
    private static TaxaList inputTaxaList;
    private static TaxaList gbsTaxaList;

    @Before
    public void setUp() throws Exception {
        inputTaxaList=TaxaListIOUtils.readTaxaAnnotationFile(taxaAnnotationFile,"<NAME>");
        gbsTaxaList=TaxaListIOUtils.readTaxaAnnotationFile(gbsTestingKeyFile,"FullSampleName");
    }

    @Test
    public void testReadTaxaAnnotationFile() throws Exception {
        System.out.println("Testing ReadTaxaAnnotationFile");
        for (Taxon taxon : inputTaxaList) {
            System.out.println(taxon.toStringWithVCFAnnotation());
        }
        Assert.assertEquals("Number of taxa read not expected number",inputTaxaList.numberOfTaxa(),281);
        int b73Index=inputTaxaList.indexOf("B73");
        Assert.assertTrue("B73 not found", b73Index>=0);
        Taxon theB73=inputTaxaList.get(inputTaxaList.indexOf("B73"));
        Assert.assertEquals("B73 not found",theB73.getAnnotation().getQuantAnnotation(Taxon.InbreedFKey)[0],0.98,0.001);
    }

    @Test
    public void testMapOfTaxonByAnnotation() throws Exception {
        System.out.println("Testing MapOfTaxonByAnnotation");
        Multimap<String,Taxon> setMap=TaxaListIOUtils.getMapOfTaxonByAnnotation(inputTaxaList,"SET");
        Assert.assertEquals("Number of sets is incorrect", 6, setMap.keySet().size(), 6);
        Assert.assertEquals("Number of sets is incorrect",39, setMap.get("CIMMYT").size());
//        for (String s : setMap.keySet()) {
//            System.out.print(s+":size="+setMap.get(s).size()+":");
//            for (Taxon taxon : setMap.get(s)) {
//                System.out.print(taxon.toString()+"  ");
//            }
//            System.out.println();
//        }
    }

    @Test
    public void testFilteringOnReading() throws Exception {
        System.out.println("Testing Filtering of KeyFile");
        TaxaList tl=TaxaListIOUtils.readTaxaAnnotationFile(gbsTestingKeyFile,"FullSampleName",ImmutableMap.of("Barcode","CTCC"),true);
        Assert.assertEquals("Number filtered taxa is wrong", 2, tl.size());
        tl=TaxaListIOUtils.readTaxaAnnotationFile(gbsTestingKeyFile,"FullSampleName",ImmutableMap.of("Col","2","Row","A"),true);
        Assert.assertEquals("Number filtered taxa is wrong", 2, tl.size());
    }

    @Test
    public void testSetSelection() throws Exception {
        System.out.println("Testing SetSelection");
        Set<String> setMap=TaxaListIOUtils.allAnnotationKeys(inputTaxaList);
        Assert.assertEquals("Number of keys", 4, setMap.size());
        for (String s : setMap) {
            System.out.println(s);
        }
        TaxaList inbredFTL=TaxaListIOUtils.retainSpecificAnnotations(inputTaxaList,new String[]{"INBREEDF"});
        setMap=TaxaListIOUtils.allAnnotationKeys(inbredFTL);
        Assert.assertEquals("Number of keys", 1, setMap.size());
        for (int i=0; i<inputTaxaList.numberOfTaxa(); i++) {
            Assert.assertEquals("INBREED keys equal", inputTaxaList.get(i).getAnnotation().getQuantAnnotation("INBREEDF")[0],
                    inbredFTL.get(i).getAnnotation().getQuantAnnotation("INBREEDF")[0],0.0001);
        }
        for (String s : setMap) {
            System.out.println(s);
        }
        TaxaList synFTL=TaxaListIOUtils.removeSpecificAnnotations(inputTaxaList,new String[]{"INBREEDF","SET","GERMTYPE"});
        setMap=TaxaListIOUtils.allAnnotationKeys(inbredFTL);
        Assert.assertEquals("Number of keys", 1, setMap.size());
        for (int i=0; i<inputTaxaList.numberOfTaxa(); i++) {
            Assert.assertEquals("GERMTYPE not deleted", 0, synFTL.get(i).getAnnotation().getTextAnnotation("GERMTYPE").length);
            for (String synonym : inputTaxaList.get(i).getAnnotation().getTextAnnotation("SYNONYM")) {
                Assert.assertEquals("SYNONYM keys not equal", synonym, synFTL.get(i).getAnnotation().getTextAnnotation("SYNONYM")[0]);
            }
        }
    }


    @Test
    public void testHapMapVCFTaxaAnnotation() throws Exception {
        GenotypeTable gtOrig=ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME);
        //should check for same order
        GenotypeTable gtAnno=GenotypeTableBuilder.getInstance(gtOrig.genotypeMatrix(),gtOrig.positions(),inputTaxaList);
        ExportUtils.writeToVCF(gtAnno,GeneralConstants.TEMP_DIR + "annoGt.vcf", false);
        GenotypeTable rtVCF=ImportUtils.readFromVCF(GeneralConstants.TEMP_DIR + "annoGt.vcf",null,true);
        Assert.assertTrue("Annotations changed on VCF roundtrip",compareAnnotations(gtAnno.taxa(),rtVCF.taxa()));
        ExportUtils.writeToHapmap(gtAnno, GeneralConstants.TEMP_DIR+"annoGt.hmp.txt");
        GenotypeTable rtHMP=ImportUtils.readFromHapmap(GeneralConstants.TEMP_DIR + "annoGt.hmp.txt");
        Assert.assertTrue("Annotations changed on HMP roundtrip",compareAnnotations(gtAnno.taxa(),rtHMP.taxa()));
    }

    @Test
//    @Ignore("HDF5 exception")
    public void testHDF5TaxaAnnotation() throws Exception {
        GenotypeTable gtOrig=ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME);
        //should check for same order
        GenotypeTable gtAnno=GenotypeTableBuilder.getInstance(gtOrig.genotypeMatrix(),gtOrig.positions(),inputTaxaList);
        if(Files.exists(Paths.get(GeneralConstants.TEMP_DIR+"annoGt.hmp.hd5").toAbsolutePath())) Files.delete(Paths.get(GeneralConstants.TEMP_DIR+"annoGt.hmp.hd5").toAbsolutePath());
        ExportUtils.writeGenotypeHDF5(gtAnno, GeneralConstants.TEMP_DIR+"annoGt.hmp.hd5");
        IHDF5Reader reader=HDF5Factory.openForReading(GeneralConstants.TEMP_DIR + "annoGt.hmp.hd5");
        TaxaList tl=new TaxaListBuilder().buildFromHDF5(reader);
        GenotypeTable rtHDF5=ImportUtils.readGuessFormat(GeneralConstants.TEMP_DIR + "annoGt.hmp.hd5");
//        GenotypeTable rtVCF=ImportUtils.readFromVCF(GeneralConstants.TEMP_DIR + "annoGt.vcf",null,true);
//        Assert.assertTrue("Annotations changed on VCF roundtrip",compareAnnotations(gtAnno.taxa(),rtVCF.taxa()));
//        ExportUtils.writeToHapmap(gtAnno, GeneralConstants.TEMP_DIR+"annoGt.hmp.txt");
//        GenotypeTable rtHMP=ImportUtils.readFromHapmap(GeneralConstants.TEMP_DIR + "annoGt.hmp.txt");
//        Assert.assertTrue("Annotations changed on HMP roundtrip",compareAnnotations(gtAnno.taxa(),rtHMP.taxa()));
    }

    private static boolean compareAnnotations(TaxaList tl1, TaxaList tl2) {
        boolean flagSame=true;
        for (int i=0; i<tl1.numberOfTaxa(); i++) {
            if(!tl1.get(i).toStringWithVCFAnnotation().equals(tl2.get(i).toStringWithVCFAnnotation())) {
                System.out.println("Annotation changed:"+tl1.get(i).toStringWithVCFAnnotation()+" versus"+
                        tl2.get(i).toStringWithVCFAnnotation());
                flagSame=false;
            }
        }
        return flagSame;
    }


}
