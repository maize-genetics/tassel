package net.maizegenetics.dna.snp;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.genotypecall.BasicGenotypeMergeRule;
import net.maizegenetics.dna.snp.io.BuilderFromHapMap;
import net.maizegenetics.taxa.Taxon;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class AlignmentBuilderTest {

    private static void clearTestFolder() {
        File f=new File(GeneralConstants.TEMP_DIR);
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
    public void testSorting() throws Exception {
        GenotypeTable a=BuilderFromHapMap.getBuilder(TutorialConstants.HAPMAP_FILENAME).build();
        GenotypeTableBuilder gtb=GenotypeTableBuilder.getTaxaIncremental(a.positions());
        for (int i=0; i<a.numberOfTaxa(); i++) {
            gtb.addTaxon(a.taxa().get(i),a.genotypeAllSites(i));
        }
        GenotypeTable b=gtb.build();
        AlignmentTestingUtils.alignmentsEqual(a,b);
        GenotypeTable bs=gtb.sortTaxa().build();
        int aIndex=0;
        for (Taxon taxon : a.taxa()) {
            int bIndex=bs.taxa().indexOf(taxon);
            Assert.assertArrayEquals(a.genotypeAllSites(aIndex),bs.genotypeAllSites(bIndex));
            aIndex++;
        }

        GenotypeTable as=BuilderFromHapMap.getBuilder(TutorialConstants.HAPMAP_FILENAME).sortTaxa().build();
        AlignmentTestingUtils.alignmentsEqual(as,bs);
    }

    @Test
    public void testGetInstanceFromHDF5() throws Exception {
        clearTestFolder();
        GenotypeTable a=BuilderFromHapMap.getBuilder(TutorialConstants.HAPMAP_FILENAME).build(); //need a sort taxa method to make the two approaches equal
        GenotypeTable hdfExport=GenotypeTableBuilder.getInstance(a,GeneralConstants.TEMP_DIR+"regHmpToH5.t5.h5");
        GenotypeTable a2=GenotypeTableBuilder.getInstance(GeneralConstants.TEMP_DIR+"regHmpToH5.t5.h5");
        ExportUtils.writeToHapmap(a2,false,GeneralConstants.TEMP_DIR+"regHmpToH5.hmp.txt.gz",'\t',null);
        GenotypeTable a3=BuilderFromHapMap.getBuilder(GeneralConstants.TEMP_DIR+"regHmpToH5.hmp.txt.gz").build(); //need a sort taxa method to make the two approaches equal

//        for (int i=0; i<a2.getSiteCount(); i++) {
//            System.out.println(i+":"+Arrays.deepToString(a2.getAllelesSortedByFrequency(i)));
//            System.out.println(i+":"+Arrays.deepToString(a.getAllelesSortedByFrequency(i)));
//        }
//        System.out.println(a.getBaseAsStringRow(0));
//        System.out.println(a2.getBaseAsStringRow(0));
//        System.out.println(a.getBaseAsStringRow(52));
//        System.out.println(a2.getBaseAsStringRow(52));
//        for (int i=0; i<a.getTaxaCount(); i++) {
//            System.out.println(i+":"+a.getFullTaxaName(i) +">>"+a2.getFullTaxaName(i));
//        }
        AlignmentTestingUtils.alignmentsEqual(a2,a3);
    }

    @Test
//    @Ignore("HDF5 exception")
    public void testDuplicateTaxaAddition() throws Exception {
        clearTestFolder();
        GenotypeTable a=BuilderFromHapMap.getBuilder(TutorialConstants.HAPMAP_FILENAME).build(); //need a sort taxa method to make the two approaches equal
        System.out.println("Basic replication of GenotypeTable using GenotypeTableBuilder");
        GenotypeTableBuilder gtb=GenotypeTableBuilder.getTaxaIncremental(a.positions());
        for (int i=0; i<a.numberOfTaxa(); i++) {
            gtb.addTaxon(a.taxa().get(i),a.genotypeAllSites(i));
        }
        AlignmentTestingUtils.alignmentsEqual(a,gtb.build());
        System.out.println("Building of GenotypeTable in two step with append using GenotypeTableBuilder");
        GenotypeTableBuilder gtb2=GenotypeTableBuilder.getTaxaIncremental(a.positions());
        for (int i=0; i<100; i++) {
            gtb2.addTaxon(a.taxa().get(i),a.genotypeAllSites(i));
        }
        GenotypeTable gt2=gtb2.build();
        GenotypeTableBuilder gtbAdd=GenotypeTableBuilder.getTaxaIncremental(gt2, new BasicGenotypeMergeRule(0.01));
        for (int i=100; i<a.numberOfTaxa(); i++) {
            gtbAdd.addTaxon(a.taxa().get(i),a.genotypeAllSites(i));
        }
        AlignmentTestingUtils.alignmentsEqual(a,gtbAdd.build());
        System.out.println("Building of GenotypeTable in two step with duplicate identical taxa using GenotypeTableBuilder");
        GenotypeTableBuilder gtbAddMerge=GenotypeTableBuilder.getTaxaIncremental(gt2, new BasicGenotypeMergeRule(0.01));
        for (int i=0; i<a.numberOfTaxa(); i++) {
          //  System.out.println(i+":"+a.taxa().get(i).toString());
            gtbAddMerge.addTaxon(a.taxa().get(i), a.genotypeAllSites(i));
        }
        AlignmentTestingUtils.alignmentsEqual(a,gtbAddMerge.build());

        //AlleleDepth randomAd=AlignmentTestingUtils.createRandomDepth(a);
        GenotypeTable aWithDepth=AlignmentTestingUtils.createRandomDepthForGenotypeTable(a,1);
        GenotypeTableBuilder gtbAddWithDepth=GenotypeTableBuilder.getTaxaIncremental(aWithDepth, new BasicGenotypeMergeRule(0.01));
        for (int i=0; i<a.numberOfTaxa(); i++) {
            gtbAddWithDepth.addTaxon(aWithDepth.taxa().get(i), aWithDepth.genotypeAllSites(i), aWithDepth.depth().valuesForTaxonByte(i));
        }
        GenotypeTable aWith2XDepth=gtbAddWithDepth.build();
        for (int t=0; t<aWithDepth.numberOfTaxa(); t++) {
            for (int s=0; s<aWithDepth.numberOfSites(); s++) {
                int[] oD=aWithDepth.depthForAlleles(t,s);
                for (int i=0; i<oD.length; i++) {oD[i]*=2;}
                int[] x2D=aWith2XDepth.depthForAlleles(t,s);
//                System.out.println(Arrays.toString(oD)+":"+Arrays.toString(x2D));
                boolean same=GenotypeTableUtils.isEqual(aWithDepth.genotype(t,s),aWith2XDepth.genotype(t,s));
                Assert.assertTrue("Error"+aWithDepth.genotypeAsString(t,s)+" values"+Arrays.toString(oD),same);
                Assert.assertArrayEquals("Depth not doubled correctly",oD,x2D);
            }
        }
    }

    @Test
//    @Ignore("HDF5 exception")
    public void testDuplicateTaxaAdditionHDF5() throws Exception {
        clearTestFolder();
        GenotypeTable a=BuilderFromHapMap.getBuilder(TutorialConstants.HAPMAP_FILENAME).sortTaxa().build(); //need a sort taxa method to make the two approaches equal
        System.out.println("Basic replication of GenotypeTable using GenotypeTableBuilder");
//        AlleleDepth randomAd=AlignmentTestingUtils.createRandomDepth(a);
        GenotypeTable aWithDepth=AlignmentTestingUtils.createRandomDepthForGenotypeTable(a,1);
        GenotypeTableBuilder gtb=GenotypeTableBuilder.getTaxaIncremental(a.positions(),GeneralConstants.TEMP_DIR+"depHmpToH5.t5.h5");
        for (int i=0; i<aWithDepth.numberOfTaxa(); i++) {
            gtb.addTaxon(aWithDepth.taxa().get(i),aWithDepth.genotypeAllSites(i),aWithDepth.depth().valuesForTaxonByte(i));
        }
        gtb.closeUnfinished();
        //TODO following line throws error as TaxaBuilderList has an open read only.  Better strategy?
        GenotypeTableBuilder gtb2=GenotypeTableBuilder.mergeTaxaIncremental(GeneralConstants.TEMP_DIR+"depHmpToH5.t5.h5",
                new BasicGenotypeMergeRule(0.001));
        for (int i=0; i<aWithDepth.numberOfTaxa(); i++) {
            gtb2.addTaxon(aWithDepth.taxa().get(i),aWithDepth.genotypeAllSites(i),aWithDepth.depth().valuesForTaxonByte(i));
        }
        GenotypeTable aWith2XDepth=gtb2.build();
        for (int t=0; t<aWithDepth.numberOfTaxa(); t++) {
 //           System.out.println(aWithDepth.taxaName(t)+"->"+aWith2XDepth.taxaName(t));
            for (int s=0; s<aWithDepth.numberOfSites(); s++) {
                int[] oD=aWithDepth.depthForAlleles(t,s);
                for (int i=0; i<oD.length; i++) {oD[i]*=2;}
                int[] x2D=aWith2XDepth.depthForAlleles(t,s);
 //               System.out.println(s+":"+Arrays.toString(oD)+":"+Arrays.toString(x2D));
                boolean same=GenotypeTableUtils.isEqual(aWithDepth.genotype(t,s),aWith2XDepth.genotype(t,s));
                Assert.assertTrue("Error"+aWithDepth.genotypeAsString(t,s)+" values"+Arrays.toString(oD),same);
                Assert.assertArrayEquals("Depth not doubled correctly",oD,x2D);
            }
        }

    }



    @Test
    public void testSingleSiteChr() throws Exception {
        GenotypeTable a=BuilderFromHapMap.getBuilder(TutorialConstants.HAPMAP_FILENAME).build();
        GenotypeTableBuilder siteBuilder=GenotypeTableBuilder.getSiteIncremental(a.taxa());
        int countOfCurrChr=0;
        int chrMax=10;
        for (int i=1; i<a.numberOfSites(); i++) {
            if(a.chromosome(i-1)!=a.chromosome(i)) {
                countOfCurrChr=0;
                chrMax=(a.chromosome(i).getChromosomeNumber()%2==0)?1:10;
                System.out.println(a.chromosome(i).toString()+":"+chrMax);
            }
            countOfCurrChr++;
            if(countOfCurrChr<=chrMax) {
                siteBuilder.addSite(a.positions().get(i),a.genotypeAllTaxa(i));
            }
        }
        GenotypeTable gt=siteBuilder.build();

    }
//
//    @Test
//    public void testGetInstanceOnlyMajorMinor() throws Exception {
//
//    }
//
//    @Test
//    public void testGetHomozygousInstance() throws Exception {
//
//    }
//
//    @Test
//    public void testGetGenotypeCopyInstance() throws Exception {
//
//    }
//
//    @Test
//    public void testGetGenotypeCopyInstance() throws Exception {
//
//    }
}
