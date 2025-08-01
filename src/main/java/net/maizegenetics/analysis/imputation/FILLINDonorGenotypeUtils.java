package net.maizegenetics.analysis.imputation;

import com.google.common.primitives.Ints;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import static net.maizegenetics.dna.WHICH_ALLELE.Major;
import static net.maizegenetics.dna.WHICH_ALLELE.Minor;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

/**
 * Methods for loading the donor haplotype files and for arranging the bit states (Major versus Minor) if they differ between
 * the donor genotypes and the target genotypes
 *
 * @author Ed Buckler
 */
public class FILLINDonorGenotypeUtils {

    public static GenotypeTable[] loadDonors(String donorFile, GenotypeTable unimpAlign, int minTestSites,
                                             boolean verboseOutput, int appoxSitesPerDonorGenotypeTable) {
        try {
            File d= new File(donorFile);
            boolean containsDonors= false;
            if (d.isDirectory()) for (File file:new File(donorFile).listFiles()) {if (file.getName().contains(".gc")) containsDonors= true;}
            if (containsDonors) {
                return loadDonors(donorFile, unimpAlign, minTestSites, verboseOutput);}
            else { return loadDonorAndChunk(donorFile, unimpAlign, appoxSitesPerDonorGenotypeTable, verboseOutput);}
        }
        catch (IllegalArgumentException e){
            throw new IllegalArgumentException("Incorrect donor file directory supplied. Must contain files with '.gc' within the file name,\n" +
                    "and match a set of files output from FILLINFindHaplotypesPlugin()");
        }

    }

    public static GenotypeTable[] loadDonors(String donorFileRoot, GenotypeTable unimpAlign, int minTestSites, boolean verboseOutput){
        File theDF=new File(donorFileRoot);
//        String prefilter=theDF.getName().split(".gX.")[0]+".gc"; //grabs the left side of the file
//        String prefilterOld=theDF.getName().split("s\\+")[0]+"s"; //grabs the left side of the file
        ArrayList<File> d=new ArrayList<File>();
        for (File file : theDF.listFiles()) {
            if(file.getName().contains(".gc")) {d.add(file);}
//            if(file.getName().startsWith(prefilter)) {d.add(file);}
//            if(file.getName().startsWith(prefilterOld)) {d.add(file);}
        }
        PositionList targetPositions= unimpAlign.positions();
        ArrayList<GenotypeTable> donorList=new ArrayList<>();
        for (int i = 0; i < d.size(); i++) {
            if(verboseOutput) System.out.println("Starting Read");
            GenotypeTable donorAlign=ImportUtils.readFromHapmap(d.get(i).getPath());
            ArrayList<Integer> subSites= new ArrayList<>();
            PositionList donorPositions= donorAlign.positions();
            for (int j = 0; j < donorAlign.numberOfSites(); j++) {if (targetPositions.siteOfPhysicalPosition(donorPositions.physicalPositions()[j],
                    donorPositions.chromosome(j)) > -1) subSites.add(j);} //if unimputed contains donorAlign position keep in donor align
            if (subSites.size()==donorAlign.numberOfSites()) {
                donorList.add(donorAlign);
                if (verboseOutput)
                    System.out.printf("Donor file shares all sites with target:%s taxa:%d sites:%d %n", d.get(i).getPath(), donorAlign.numberOfTaxa(), donorAlign.numberOfSites());
                continue;
            }
            if (subSites.size()<2) {
                if(verboseOutput) System.out.printf("Donor file contains <2 matching sites and will not be used:%s",d.get(i).getPath());
                continue;
            } else {
                donorAlign= GenotypeTableBuilder.getGenotypeCopyInstance(FilterGenotypeTable.getInstance(donorAlign,Ints.toArray(subSites)));
                donorList.add(donorAlign);
                if(verboseOutput) System.out.printf("Donor file sites filtered to match target:%s taxa:%d sites:%d %n",
                        d.get(i).getPath(), donorAlign.numberOfTaxa(),donorAlign.numberOfSites());
                if (subSites.size() < minTestSites*2 && verboseOutput) System.out.println("This donor alignment contains " +
                        "marginally sufficient matching snp positions. Region unlikely to impute well.");
            }
        }
        return donorList.toArray(new GenotypeTable[0]);
    }

    public static GenotypeTable[] loadDonorAndChunk(String donorFile, GenotypeTable unimpAlign, int appoxSitesPerHaplotype, boolean verboseOutput){
        GenotypeTable donorMasterGT=ImportUtils.read(donorFile);
        donorMasterGT=GenotypeTableBuilder.getHomozygousInstance(donorMasterGT);
        int[][] donorFirstLastSites=FILLINFindHaplotypesPlugin.divideChromosome(donorMasterGT,appoxSitesPerHaplotype,verboseOutput);
        GenotypeTable[] donorAlign=new GenotypeTable[donorFirstLastSites.length];
        for (int i = 0; i < donorAlign.length; i++) {
            if(verboseOutput) System.out.println("Starting Read");
            donorAlign[i]=GenotypeTableBuilder.getGenotypeCopyInstance(FilterGenotypeTable.getInstance(donorMasterGT, donorFirstLastSites[i][0], donorFirstLastSites[i][1]));
        }
        return donorAlign;
    }
    
    /**
     * Remove sites in Donor and Unimp that don't match minMaj. Assumes same 
     * coordinate system. Designed to match sites between GBS and hapmap. Removes
     * sites where minMaj is indel. Keeps sites with same minMaj even if third allele differs.
     * @param donor
     * @param unimp
     * @return 
     */
    public static GenotypeTable RemoveSitesThatDoNotMatchMinMaj(String donorFile, GenotypeTable unimp, boolean verboseOutput) {
        File newDonor= new File(donorFile.subSequence(0, donorFile.lastIndexOf("."))+".matchMinMaj.hmp.txt.gz");
        if (newDonor.exists()) {
            System.out.println("Donor already filtered for minMaj");
            return unimp;
        }
        if (verboseOutput) System.out.println("Starting filter of donor and unimputed file for sites that match for minMaj at same physical position");
        GenotypeTable donor= ImportUtils.readGuessFormat(donorFile);
        ArrayList<String> sitesToKeepDonor= new ArrayList<>(); ArrayList<String> sitesToKeepUnimp= new ArrayList<>();
        int cntBadMinMaj= 0; int cntBadMinMajToAlleleNoIndel= 0; int matchPos= 0;
        byte indel= GenotypeTableUtils.getDiploidValue(NucleotideAlignmentConstants.GAP_ALLELE,NucleotideAlignmentConstants.INSERT_ALLELE);
        for (Chromosome chr:donor.chromosomes()) {
            if (Arrays.binarySearch(unimp.chromosomes(),chr)<0) continue;
            int[] pos= Arrays.copyOfRange(donor.physicalPositions(), 
                    donor.firstLastSiteOfChromosome(chr)[0],
                    donor.firstLastSiteOfChromosome(chr)[1]+1);
            for (int unimpSite = unimp.firstLastSiteOfChromosome(chr)[0]; unimpSite < unimp.firstLastSiteOfChromosome(chr)[1]+1; unimpSite++) {
//                if (genosKeep.alleles(keepSite).length>2) continue;
                int donorSiteOnChr= Arrays.binarySearch(pos, unimp.physicalPositions()[unimpSite]);
                if (donorSiteOnChr<0) continue;
                matchPos++;
                int donorSite= donor.firstLastSiteOfChromosome(chr)[0]+donorSiteOnChr;//get the actual site for Donor
//                if (genosMod.alleles(modSite).length>2) continue;
                byte minMajDonor= GenotypeTableUtils.getDiploidValue(donor.majorAllele(donorSite), donor.minorAllele(donorSite));
                byte minMajUnimp= GenotypeTableUtils.getDiploidValue(unimp.majorAllele(unimpSite), unimp.minorAllele(unimpSite));
                if (GenotypeTableUtils.isEqual(indel, minMajDonor) || GenotypeTableUtils.isEqual(indel, minMajUnimp)) continue;
                if (!GenotypeTableUtils.isEqual(minMajDonor, minMajUnimp)) {
                    cntBadMinMaj++;
                    boolean biallelicNoIndel= true;
                    if (donor.alleles(donorSite).length>2 || unimp.alleles(unimpSite).length>2 || 
                            GenotypeTableUtils.isPartiallyEqual(minMajDonor, indel) || GenotypeTableUtils.isPartiallyEqual(minMajUnimp, indel)) 
                    {biallelicNoIndel= false; cntBadMinMajToAlleleNoIndel++;}
//                  continue;
                }
                else {
                    sitesToKeepDonor.add(donor.siteName(donorSite));
                    sitesToKeepUnimp.add(unimp.siteName(unimpSite));
                }
            }
        }
        if (sitesToKeepDonor.size()<100 || sitesToKeepUnimp.size()<100) {
            if (verboseOutput) System.out.println("Not enough sites to impute on that match minMaj at identical positions");
            return null;
        }
        GenotypeTable filterDonor= FilterGenotypeTable.getInstance(donor, sitesToKeepDonor);
        GenotypeTable filterUnimp= FilterGenotypeTable.getInstance(unimp, sitesToKeepUnimp);
        ExportUtils.writeToHapmap(filterDonor, true, newDonor.getAbsolutePath(), '\t', null);
        if (verboseOutput) System.out.println("Wrote donor file to "+newDonor.getAbsolutePath());
        if (verboseOutput) System.out.println("Sites that match positionally: "+matchPos+
                "\nSites removed because of inconsistent min/maj: "+cntBadMinMaj+
                "\nSites removed because of inconsistent min/maj at sites without indels or third alleles: "+cntBadMinMajToAlleleNoIndel+
                "\nSites retained: "+filterDonor.numberOfSites()+
                "\nAccuracy: "+(double)filterDonor.numberOfSites()/((double)filterDonor.numberOfSites()+(double)cntBadMinMaj+(double)cntBadMinMajToAlleleNoIndel));
        return filterUnimp;
    }
    
    /**
     * Create mask for all sites where major & minor are swapped between GenotypeTables
     * returns in this order [goodMask, swapMjMnMask, errorMask, invariantMask]
     */
    public static OpenBitSet[][] createMaskForAlignmentConflicts(GenotypeTable unimpAlign, GenotypeTable[] donorAlign,
                                                                 boolean print) {
        OpenBitSet[][] result=new OpenBitSet[donorAlign.length][4];
        for (int da = 0; da < result.length; da++) {
//            int donorOffset=unimpAlign.siteOfPhysicalPosition(donorAlign[da].chromosomalPosition(0), donorAlign[da].chromosome(0), donorAlign[da].siteName(0));
            int donorOffset=unimpAlign.positions().siteOfPhysicalPosition(donorAlign[da].positions().physicalPositions()[0], donorAlign[da].positions().chromosome(0));
            OpenBitSet goodMask=new OpenBitSet(donorAlign[da].numberOfSites());
            OpenBitSet swapMjMnMask=new OpenBitSet(donorAlign[da].numberOfSites());
            OpenBitSet errorMask=new OpenBitSet(donorAlign[da].numberOfSites());
            OpenBitSet invariantMask=new OpenBitSet(donorAlign[da].numberOfSites());
            int siteConflicts=0, swaps=0, invariant=0, good=0;
            for (int i = 0; i < donorAlign[da].numberOfSites(); i++) {
                /*we have three classes of data:  invariant in one alignment, conflicts about minor and minor,
                *swaps of major and minor.  Adding the invariant reduces imputation accuracy.
                *the major/minor swaps should be flipped in the comparisons
                */
                byte tMj=unimpAlign.majorAllele(i + donorOffset);
                byte tMn=unimpAlign.minorAllele(i + donorOffset);
                byte daMj=donorAlign[da].majorAllele(i);
                byte daMn=donorAlign[da].minorAllele(i);
                if(daMn==GenotypeTable.UNKNOWN_ALLELE) {
                    invariant++;
                    invariantMask.set(i);
                    goodMask.set(i);
                } else
                if((daMj==tMn)&&(daMn==tMj)) {
                    swaps++;
                    swapMjMnMask.set(i);
                    goodMask.set(i);
                } else
                if((daMj!=tMj)) {
                    siteConflicts++;
                    errorMask.set(i);
                    goodMask.set(i);
                }

            }
            goodMask.not();
            if(print) System.out.println("Donor:"+da+"invariant in donor:"+invariant+" swapConflicts:"+swaps+" errors:"+siteConflicts);
            result[da]=new OpenBitSet[] {goodMask, swapMjMnMask, errorMask, invariantMask};
        }
        return result;
    }


    /**
     * Major and minor allele get swapped between GenotypeTables.  This method flips these, and sets the bad sites to missing
     */
    public static BitSet[] arrangeMajorMinorBtwAlignments(GenotypeTable unimpAlign, int bt, int donorOffset, int donorLength,
                                                    OpenBitSet goodMask, OpenBitSet swapMjMnMask, boolean isSwapMajorMinor) {
        int unimpAlignStartBlock=donorOffset/64;
        int shift=(donorOffset-(unimpAlignStartBlock*64));
        int unimpAlignEndBlock=unimpAlignStartBlock+((donorLength+shift-1)/64);
        OpenBitSet mjUnImp=new OpenBitSet(unimpAlign.allelePresenceForSitesBlock(bt, Major, unimpAlignStartBlock, unimpAlignEndBlock + 1));
        OpenBitSet mnUnImp=new OpenBitSet(unimpAlign.allelePresenceForSitesBlock(bt, Minor, unimpAlignStartBlock, unimpAlignEndBlock + 1));
        OpenBitSet mjTbs=new OpenBitSet(donorLength);
        OpenBitSet mnTbs=new OpenBitSet(donorLength);
        for (int i = 0; i < donorLength; i++) {
            if(mjUnImp.fastGet(i+shift)) mjTbs.set(i);
            if(mnUnImp.fastGet(i+shift)) mnTbs.set(i);
        }
        OpenBitSet newmj=new OpenBitSet(mjTbs);
        OpenBitSet newmn=new OpenBitSet(mnTbs);
        mjTbs.and(goodMask);
        mnTbs.and(goodMask);
        //       System.out.printf("mjTbs:%d Goodmask:%d swapMjMnMask:%d",mjTbs.getNumWords(),goodMask.getNumWords(), swapMjMnMask.getNumWords());
        if(isSwapMajorMinor) {
            newmj.and(swapMjMnMask);
            newmn.and(swapMjMnMask);
            mjTbs.or(newmn);
            mnTbs.or(newmj);
        }
//        System.out.printf("Arrange shift:%d finalEnd:%d totalBits:%d DesiredLength:%d %n", shift, finalEnd, totalBits, donorLength);
        BitSet[] result={mjTbs,mnTbs};
        return result;
    }


}
