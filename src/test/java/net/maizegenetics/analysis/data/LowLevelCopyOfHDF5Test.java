package net.maizegenetics.analysis.data;

import com.google.common.collect.ImmutableMap;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.AlignmentTestingUtils;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.io.BuilderFromHapMap;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class LowLevelCopyOfHDF5Test {
    private static final String GBS_DEPTH_HDF5=GeneralConstants.DATA_DIR+"GenotypeTableTests/maize_chr9_10.h5";
    public static String EXPORT_TEMP_DIR = GeneralConstants.TEMP_DIR + "ExportUtilsTest/";
    public static String HAPMAP_TEMP_FILENAME = "temp_hapmap.hmp.txt";

    private static void clearTestFolder() {
        File f=new File(EXPORT_TEMP_DIR);
        f.mkdirs();
        System.out.println(f.getAbsolutePath());
        if(f.listFiles()!=null) {
            for (File file : f.listFiles()) {
                if(file.getName().contains("hmp.") || file.getName().contains(".h5")) {
                    file.delete();
                    System.out.println("Deleting:"+file.toString());
                }
            }
        }
    }



    @Test
    public void testSmallGenotypesToNewFile() throws Exception {
        clearTestFolder();
        int reps=2;
        GenotypeTable a=BuilderFromHapMap.getBuilder(TutorialConstants.HAPMAP_FILENAME).sortTaxa().build(); //need a sort taxa method to make the two approaches equal
     //   GenotypeTable a=BuilderFromHapMap.getBuilder("/Users/edbuckler/Development/Tassel4TestingGit/dataFiles/CandidateTests/Ames105v26.chr10.hmp.txt.gz").sortTaxa().build();
        System.out.println("Testing low level copy of HDF5 with Depth");
        GenotypeTable aWithDepth=AlignmentTestingUtils.createRandomDepthForGenotypeTable(a, 1);
        ExportUtils.writeGenotypeHDF5(aWithDepth, EXPORT_TEMP_DIR+"smDepth.t5.h5", true);

        //GenotypeTable gt=ImportUtils.readGuessFormat(EXPORT_TEMP_DIR+"smDepth.t5.h5");
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i=0; i<a.numberOfTaxa(); i++) {
            if(i%100!=0) tlb.add(a.taxa().get(i));
        }
        TaxaList tl=tlb.build();
        long time=System.nanoTime();
        for (int i=0; i<reps; i++) {
            LowLevelCopyOfHDF5.subsetGenotypesToNewFile(EXPORT_TEMP_DIR+"smDepth.t5.h5", EXPORT_TEMP_DIR+"smSubsetDepth"+i+".t5.h5",tl);
        }
        double bpPerNS=(double)(reps*tl.numberOfTaxa()*a.numberOfSites())/(double)(System.nanoTime()-time);
        System.out.printf("Rate of low level copy %g in bp/ns%n", bpPerNS);

        // The low-level subset copy must retain the requested taxa and all sites.
        GenotypeTable subset=ImportUtils.readGuessFormat(EXPORT_TEMP_DIR+"smSubsetDepth0.t5.h5");
        assertEquals("Subset file has an unexpected number of taxa", tl.numberOfTaxa(), subset.numberOfTaxa());
        assertEquals("Subset file has an unexpected number of sites", a.numberOfSites(), subset.numberOfSites());

        time=System.nanoTime();
        for (int i=0; i<reps; i++) {
            ExportUtils.writeGenotypeHDF5(aWithDepth,EXPORT_TEMP_DIR+"smExportDepth"+i+".t5.h5",true);
        }
        double exportbpPerNS=(double)(reps*aWithDepth.numberOfTaxa()*aWithDepth.numberOfSites())/(double)(System.nanoTime()-time);
        System.out.printf("Rate of ExportUtils copy %g in bp/ns%n", exportbpPerNS);
        time=System.nanoTime();
        for (int i=0; i<reps; i++) {
            GenotypeTable rt=ImportUtils.readGuessFormat(EXPORT_TEMP_DIR+"smDepth.t5.h5");
            ExportUtils.writeGenotypeHDF5(rt,EXPORT_TEMP_DIR+"smImpExpDepth"+i+".t5.h5",true);
        }
        double impExpbpPerNS=(double)(reps*aWithDepth.numberOfTaxa()*aWithDepth.numberOfSites())/(double)(System.nanoTime()-time);
        System.out.printf("Rate of ExportUtils copy %g in bp/ns%n", impExpbpPerNS);
        //  System.out.println(gt.numberOfTaxa());
        System.out.printf("Relative rates Low/Exp:%g  Low/ImpExp:%g %n",bpPerNS/exportbpPerNS,bpPerNS/impExpbpPerNS);

//        time=System.nanoTime();
//        for (int i=0; i<reps; i++) {
//            LowLevelCopyOfHDF5.subsetGenotypesToNewFile(EXPORT_TEMP_DIR+"smDepth.t5.h5", EXPORT_TEMP_DIR+"smSubsetDepth2"+i+".t5.h5",tl);
//        }
//        double bpPerNS2=(double)(reps*tl.numberOfTaxa()*a.numberOfSites())/(double)(System.nanoTime()-time);
//        System.out.printf("Rate of low level copy %g in bp/ns%n", bpPerNS2);
    }

    @Test
    public void testSubsetGenotypesToNewFile() throws Exception {
        clearTestFolder();
        String subsetFile=EXPORT_TEMP_DIR+"gbsSubsetDepth.t5.h5";
        // One flowcell's samples from the GBS key file; the key lists one sample that the HDF5 file lacks.
        TaxaList tl=TaxaListIOUtils.readTaxaAnnotationFile(GBSConstants.GBS_TESTING_KEY_FILE, "FullSampleName",
                ImmutableMap.of("Flowcell", "C08L7ACXX"), true);
        LowLevelCopyOfHDF5.subsetGenotypesToNewFile(GBS_DEPTH_HDF5, subsetFile, tl);

        GenotypeTable orig=ImportUtils.readGuessFormat(GBS_DEPTH_HDF5);
        int expectedTaxa=0;
        for (Taxon taxon : tl) {
            if (orig.taxa().indexOf(taxon.getName())>=0) expectedTaxa++;
        }
        assertTrue("Expected the key file to select a strict subset of the HDF5 taxa",
                expectedTaxa>0 && expectedTaxa<orig.numberOfTaxa() && expectedTaxa<tl.numberOfTaxa());

        GenotypeTable subset=ImportUtils.readGuessFormat(subsetFile);
        assertEquals("Subset file has an unexpected number of taxa", expectedTaxa, subset.numberOfTaxa());
        assertTrue("Subset file lost its depth", subset.hasDepth());
        AlignmentTestingUtils.alignmentsEqual(FilterGenotypeTable.getInstance(orig, subset.taxa()), subset);
    }
}
