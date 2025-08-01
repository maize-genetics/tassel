package net.maizegenetics.analysis.data;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.AlignmentTestingUtils;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.io.BuilderFromHapMap;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.TaxaListIOUtils;
import org.junit.Ignore;
import org.junit.Test;

import java.io.File;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class LowLevelCopyOfHDF5Test {
    private String bigTest="/Volumes/BucklerFastBig5/GBS27play/hmp/AllZeaGBSv27wDepth_tassel5.t5.h5";
    private String bigTestOut="/Volumes/BucklerFastBig5/GBS27play/hmp/MARS27.t5.h5";
    private String bigSubTaxaList="/Volumes/BucklerFastBig5/GBS27play/hmp/jasonTaxa.txt";
    public static String EXPORT_TEMP_DIR = GeneralConstants.TEMP_DIR + "ExportUtilsTest/";
    public static String HAPMAP_TEMP_FILENAME = "temp_hapmap.hmp.txt";

    private static void clearTestFolder() {
        File f=new File(EXPORT_TEMP_DIR);
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

    @Ignore
    @Test
    public void testSubsetGenotypesToNewFile() throws Exception {
        //  GenotypeTable gt=ImportUtils.readGuessFormat(bigTest);
        TaxaList tl=TaxaListIOUtils.readTaxaAnnotationFile(bigSubTaxaList, "<Name>");
        System.out.println(tl.numberOfTaxa());
        LowLevelCopyOfHDF5.subsetGenotypesToNewFile(bigTest, bigTestOut,tl);
        //  System.out.println(gt.numberOfTaxa());

    }
}
