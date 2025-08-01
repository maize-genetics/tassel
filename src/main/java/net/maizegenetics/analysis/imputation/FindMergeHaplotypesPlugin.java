/*
 * FindMergeHaplotypesPlugin
 */
package net.maizegenetics.analysis.imputation;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.*;
import net.maizegenetics.util.BitSet;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedWriter;
import java.util.*;
import java.util.List;
import java.util.Map.Entry;

/**
 * Creates haplotypes by finding large IBS regions within GBS data.  Starts with the 
 * highest coverage taxa and looks within windows of near perfect matches.  Combines
 * all matches together into one haplotype.  The haplotype is named for the highest coverage
 * sample.  
 * 
 * TODO:
 * 1.  plus add short inbred segments not present full ones
 * 2.  Cluster and choose haplotypes by cluster and number of new minor alleles (or information)
 * 3.  Set max het frequency as a setting
 * 
 * @author edbuckler
 */
public class FindMergeHaplotypesPlugin extends AbstractPlugin {
    private int startChr, endChr;
    private int startDiv=-1, endDiv=-1;
    private String hmpFile;
    private String outFileBase;
    private String errFile=null;
    private double minJointGapProb=0.01;
    private boolean callGaps=false;
    private double maxDistFromFounder=0.01;
    private int appoxSitesPerHaplotype=8192;
    private int minSitesPresentPerHap=500;

    private double maximumMissing=0.4;
    private int maxHaplotypes=100;
    private int minSitesForSectionComp=50;
    private double maxHetFreq=0.01;
    private double maxErrorInCreatingConsensus=0.05;
    private int minTaxaInGroup=2;
    
    
    private double[] propMissing;
    private int[] siteErrors, siteCallCnt;
    private BitSet badMask=null;
    
    private static ArgsEngine engine = new ArgsEngine();
    private static final Logger myLogger = LogManager.getLogger(FindMergeHaplotypesPlugin.class);
    
    public FindMergeHaplotypesPlugin() {
        super(null, false);
    }

    public FindMergeHaplotypesPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }
    
    /**
     * Produces donor haplotypes that can be used for imputation.
     * @param inFile Input DNA sequence alignment file name in either HapMap or HDF5
     * @param exportFile
     * @param errorExportFile
     * @param maxDistance maximum distance between founder haplotype and other haplotype
     * @param minSites
     * @param appoxSitesPerHaplotype 
     */
    public void runFindMergeHaplotypes(String inFile, String exportFile,
            String errorExportFile, double maxDistance, int minSites, int appoxSitesPerHaplotype) {
        System.out.println("Reading: "+inFile);
        GenotypeTable baseAlign=ImportUtils.readGuessFormat(inFile);

        int[][] divisions=divideChromosome(baseAlign, appoxSitesPerHaplotype);  
        System.out.printf("In taxa:%d sites:%d %n",baseAlign.numberOfTaxa(),baseAlign.numberOfSites());
        siteErrors=new int[baseAlign.numberOfSites()];
        siteCallCnt=new int[baseAlign.numberOfSites()];
        if(startDiv==-1) startDiv=0;
        if(endDiv==-1) endDiv=divisions.length-1;
        for (int i = startDiv; i <=endDiv; i++) {
            GenotypeTable mna=createHaplotypeAlignment(divisions[i][0], divisions[i][1], baseAlign,
             minSites,  maxDistance);
            String newExport=exportFile.replace("sX.hmp", "s"+i+".hmp");
            newExport=newExport.replace("gX", "gc"+mna.chromosomeName(0)+"s"+i);
            ExportUtils.writeToHapmap(mna, false, newExport, '\t', null);
            if(errorExportFile!=null) exportBadSites(baseAlign, errorExportFile, 0.01);  
            mna=null;
            System.gc();
        }
        
    }
    
    private GenotypeTable createHaplotypeAlignment(int startSite, int endSite, GenotypeTable baseAlign,
            int minSites, double maxDistance) {
        GenotypeTable fa=FilterGenotypeTable.getInstance(baseAlign, startSite, endSite);
        GenotypeTable inAlign=GenotypeTableBuilder.getGenotypeCopyInstance(fa);
        int sites=inAlign.numberOfSites();
        System.out.printf("SubInAlign Locus:%s StartPos:%d taxa:%d sites:%d %n",inAlign.chromosome(0),
                inAlign.chromosomalPosition(0),inAlign.numberOfTaxa(),inAlign.numberOfSites());

        propMissing=new double[inAlign.numberOfTaxa()];
        int startBlock=0;
        int endBlock=inAlign.allelePresenceForAllSites(0, WHICH_ALLELE.Major).getNumWords();
        TreeMap<Integer,Integer> presentRanking=createPresentRankingForWindow(inAlign, startBlock, endBlock, minSites, maxHetFreq);
        System.out.printf("Block %d Inbred and modest coverage:%d %n",startBlock,presentRanking.size());
        System.out.printf("Current Site %d Current block %d EndBlock: %d %n",startSite, startBlock, endBlock);
        TreeMap<Integer,byte[][]> results=mergeWithinWindow(inAlign, presentRanking, startBlock, endBlock, maxDistance, startSite);
        TaxaListBuilder tLB=new TaxaListBuilder();
        GenotypeCallTableBuilder gB=GenotypeCallTableBuilder.getInstance(results.size(),inAlign.numberOfSites());
        int index=0;
        for (byte[][] calls : results.values()) {
            tLB.add(new Taxon("h"+index+(new String(calls[1]))));
            gB.setBaseRangeForTaxon(index,0,calls[0]);
            index++;
        }
        return GenotypeTableBuilder.getInstance(gB.build(),inAlign.positions(),tLB.build());
    }
    
    public static int[][] divideChromosome(GenotypeTable a, int appoxSitesPerHaplotype) {
        Chromosome[] theL=a.chromosomes();
        ArrayList<int[]> allDivisions=new ArrayList<int[]>();
        for (Chromosome aL: theL) {
            System.out.println("");
            int[] startEnd=a.positions().startAndEndOfChromosome(aL);
            //todo chromosome offsets will be need to replace this
            int locusSites=startEnd[1]-startEnd[0]+1;
            int subAlignCnt=(int)Math.round((double)locusSites/(double)appoxSitesPerHaplotype);
            if(subAlignCnt==0) subAlignCnt++;
            int prefBlocks=(locusSites/(subAlignCnt*64));
            System.out.printf("Chr:%s Alignment Sites:%d subAlignCnt:%d RealSites:%d %n",
                    aL.getName(),locusSites, subAlignCnt, prefBlocks*64);
            for (int i = 0; i < subAlignCnt; i++) {
                int[] divs=new int[2];
                divs[0]=(i*prefBlocks*64)+startEnd[0];
                divs[1]=divs[0]+(prefBlocks*64)-1;
                if(i==subAlignCnt-1) divs[1]=startEnd[1];
                allDivisions.add(divs);
            }
        }
        int[][] result=new int[allDivisions.size()][2];
        for (int i = 0; i < result.length; i++) {
            result[i]=allDivisions.get(i);
           // 
            System.out.printf("Chromosome Divisions: %s start:%d end:%d %n", a.chromosome(result[i][0]).getName(),
                    result[i][0], result[i][1]);
        }
        return result;
    }
    
    private TreeMap<Integer,Integer> createPresentRankingForWindow(GenotypeTable inAlign, int startBlock, int endBlock, 
            int minSites, double maxHetFreq) {
        int sites=64*(endBlock-startBlock+1);
        TreeMap<Integer,Integer> presentRanking=new TreeMap<Integer,Integer>(Collections.reverseOrder());
        for (int i = 0; i < inAlign.numberOfTaxa(); i++) {
            long[] mj=inAlign.allelePresenceForSitesBlock(i, WHICH_ALLELE.Major, startBlock, endBlock);
            long[] mn=inAlign.allelePresenceForSitesBlock(i, WHICH_ALLELE.Minor, startBlock, endBlock);
            int totalSitesNotMissing = 0;
            int hetCnt=0;
            for (int j = 0; j < mj.length; j++) {
                totalSitesNotMissing+=BitUtil.pop(mj[j]|mn[j]);
                hetCnt+=BitUtil.pop(mj[j]&mn[j]);
            }
            double hetFreq=(double)hetCnt/(double)totalSitesNotMissing;
            propMissing[i]=(double)(1+sites-totalSitesNotMissing)/(double)sites; //1 prevents error in joint missing calculations
            double propPresent=1.0-propMissing[i];
            if((hetFreq>maxHetFreq)||(totalSitesNotMissing<minSites)) continue;
            int index=(1000000*((int)(propPresent*100)))+i;
//            System.out.println(index);
            presentRanking.put(index, i);
        }
        return presentRanking;
    }
    
    private void exportBadSites(GenotypeTable baseAlign, String exportMap, double errorThreshold) {
        BufferedWriter bw = null;
        try {
            String fullFileName = Utils.addSuffixIfNeeded(exportMap, ".txt", new String[]{".txt"});
            bw = Utils.getBufferedWriter(fullFileName);
            bw.write("<Map>\n");
            for (int i = 0; i < baseAlign.numberOfSites(); i++) {
                double errorsRate=(double)siteErrors[i]/(double)siteCallCnt[i];
                if(errorsRate<errorThreshold) continue;
                bw.write(baseAlign.siteName(i)+"\t");
                bw.write(baseAlign.chromosomeName(i) +"\t");
                bw.write(i+"\t"); //dummy for genetic position
                bw.write(baseAlign.chromosomalPosition(i) +"\n"); //dummy for genetic position
            } 
            bw.close();
            
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error writing GeneticMap file: " + exportMap + ": " + ExceptionUtils.getExceptionCauses(e));
        }
    }
    
    private TreeMap<Integer,byte[][]> mergeWithinWindow(GenotypeTable inAlign, TreeMap<Integer,Integer> presentRanking,
            int startBlock, int endBlock, double maxDistance, int siteOffsetForError){
        int startSite=startBlock*64;
        int endSite=63+(endBlock*64);
        if(endSite>=inAlign.numberOfSites()) endSite=inAlign.numberOfSites()-1;
        TreeMap<Integer,ArrayList> mergeSets=new TreeMap<Integer,ArrayList>();
        TreeMap<Integer,byte[][]> results=new TreeMap<Integer,byte[][]>(Collections.reverseOrder());
        TreeSet<Integer> unmatched=new TreeSet<Integer>(presentRanking.values());
        TaxaList inIDG=inAlign.taxa();
        for (Entry<Integer,Integer> e : presentRanking.entrySet()) {
            int taxon1=e.getValue();
            if(unmatched.contains(taxon1)==false) continue;//already included in another group
            ArrayList<Integer> hits=new ArrayList<Integer>();
            unmatched.remove(taxon1);
            for(int taxon2 : unmatched) {
               double[] dist=IBSDistanceMatrix.computeHetBitDistances(inAlign, 
                       taxon1, taxon2, minSitesForSectionComp, startBlock, endBlock, badMask);
               if((!Double.isNaN(dist[0]))&&(dist[0]<maxDistance)) {
                   hits.add(taxon2);
               }
            }
            byte[] calls=null;
            if(((hits.size()+1)<this.minTaxaInGroup)) continue;
            if(hits.size()>0) {
                ArrayList<String> mergeNames=new ArrayList<String>();
                mergeNames.add(inIDG.taxaName(taxon1));
                mergeSets.put(taxon1, hits);         
               // System.out.print(inAlign.taxaName(taxon1)+"=");
                for (Integer taxon2 : hits) {
                    unmatched.remove(taxon2);
                   // System.out.print(inAlign.taxaName(taxon2)+"=");
                    mergeNames.add(inIDG.taxaName(taxon2));
                }
              //  System.out.println("");              
                calls=consensusGameteCalls(inAlign, mergeNames, startSite, endSite, maxErrorInCreatingConsensus, siteOffsetForError);
            } else {
                calls=inAlign.genotypeRange(taxon1, startSite, endSite+1);
            }
            int[] unkCnt=countUnknown(calls);
            double missingFreq=(double)unkCnt[0]/(double)inAlign.numberOfSites();
            double hetFreq=(double)unkCnt[1]/(double)(inAlign.numberOfSites()-unkCnt[0]);
            if(((missingFreq<maximumMissing)&&(hetFreq<maxHetFreq))) {
                int index=(hits.size()*200000)+taxon1;
                System.out.printf("Output %s plus %d missingF:%g hetF:%g index: %d %n",inIDG.taxaName(taxon1),
                        hits.size(), missingFreq, hetFreq, index);
                byte[][] callPlusNames=new byte[2][];
                callPlusNames[0]=calls;
                String newName=inIDG.get(taxon1).getName()+":d"+(hits.size()+1);
                callPlusNames[1]=newName.getBytes();
                results.put(index, callPlusNames);
            }
            if(results.size()>=maxHaplotypes) break;
        }
        return results;
    }
    
    private byte[] consensusGameteCalls(GenotypeTable a, List<String> taxa, int startSite, 
            int endSite, double maxError, int siteOffsetForError) {
        int[] taxaIndex = new int[taxa.size()];
        for (int t = 0; t < taxaIndex.length; t++) {  //why are we working with names rather than numbers
            taxaIndex[t] = a.taxa().indexOf(taxa.get(t));
        }
        byte[] calls = new byte[endSite-startSite+1];
        Arrays.fill(calls, GenotypeTable.UNKNOWN_DIPLOID_ALLELE);
        for (int s = startSite; s <= endSite; s++) {
            byte mjAllele = a.majorAllele(s);
            byte mnAllele = a.minorAllele(s);
            byte mj=GenotypeTableUtils.getUnphasedDiploidValue(mjAllele,mjAllele);
            byte mn=GenotypeTableUtils.getUnphasedDiploidValue(mnAllele,mnAllele);
            byte het = GenotypeTableUtils.getUnphasedDiploidValue(mjAllele, mnAllele);
            int mjCnt=0, mnCnt=0;
            for (int t = 0; t < taxaIndex.length; t++) {
                byte ob = a.genotype(taxaIndex[t], s);
                if (ob == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    continue;
                }
                if (ob == mj) {
                    mjCnt++;
                } else if (ob == mn) {
                    mnCnt++;
                } else if (GenotypeTableUtils.isEqual(ob, het)) {
                    mjCnt++;
                    mnCnt++;
                }
            }
            int totalCnt = mjCnt + mnCnt;
            
            if (totalCnt == 0) {
                double missingProp=1.0;
                for (int t : taxaIndex) {missingProp*=propMissing[t];}
                if(callGaps&(missingProp<minJointGapProb)) calls[s-startSite]=NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE;
                continue;
            }
            double errorRate;
            if(totalCnt>1) siteCallCnt[s+siteOffsetForError]+=totalCnt;
            if(mjCnt < mnCnt) {
                errorRate=(double)mjCnt/(double)totalCnt;
                if(errorRate<maxError) {calls[s-startSite] = mn;}
                else {siteErrors[s+siteOffsetForError]+=mjCnt;}
            } else {
                errorRate=(double)mnCnt/(double)totalCnt;
                if(errorRate<maxError) {calls[s-startSite] = mj;}
                else {siteErrors[s+siteOffsetForError]+=mnCnt;}
            }
        }
        return calls;
    }
    
    public static ArrayList<Integer> maxMajorAllelesTaxa(GenotypeTable a, int numMaxTaxa, WHICH_ALLELE alleleNumber) {
        ArrayList<Integer> maxTaxa=new ArrayList<Integer>();
        OpenBitSet curMj=new OpenBitSet(a.numberOfSites());
        long maxMjCnt=curMj.cardinality();
        for (int i = 0; i < numMaxTaxa; i++) {
            long bestCnt=0;
            int bestAddTaxa=-1;
            for (int t = 0; t < a.numberOfTaxa(); t++) {
                BitSet testMj=new OpenBitSet(a.allelePresenceForAllSites(t, alleleNumber));
                testMj.union(curMj);
                long cnt=testMj.cardinality();
                if(cnt>bestCnt) {
                    bestCnt=cnt;
                    bestAddTaxa=t;
                }
            }
            if(maxMjCnt==bestCnt) continue;
            curMj.union(a.allelePresenceForAllSites(bestAddTaxa, alleleNumber));
            maxMjCnt=curMj.cardinality();
            maxTaxa.add(bestAddTaxa);
            System.out.printf("Allele:%d Taxa: %s %d %n",alleleNumber,a.taxaName(bestAddTaxa),maxMjCnt);
        }
        return maxTaxa;
    }
    
    private int[] countUnknown(byte[] b) {
        int cnt=0, cntHet=0;
        for (int i = 0; i < b.length; i++) {
            if(b[i]==GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {cnt++;}
            else if(GenotypeTableUtils.isHeterozygous(b[i])) {cntHet++;}
        }
        int[] result={cnt,cntHet};
        return result;
    }
    
        @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        engine.add("-hmp", "-hmpFile", true);
        engine.add("-o", "--outFile", true);
        engine.add("-oE", "--outErrorFile", true);
        engine.add("-sC", "--startChrom", true);
        engine.add("-eC", "--endChrom", true);
        engine.add("-mxDiv", "--mxDiv", true);
        engine.add("-mxHet", "--mxHet", true);
        engine.add("-hapSize", "--hapSize", true);
        engine.add("-minPres", "--minPres", true);
        engine.add("-maxHap", "--maxHap", true);
        engine.add("-maxOutMiss", "--maxOutMiss", true);
        engine.add("-sD", "--startDivision", true);
        engine.add("-eD", "--endDivision", true);
        engine.parse(args);
        if (engine.getBoolean("-sC")) {
            startChr = Integer.parseInt(engine.getString("-sC"));
        }
        if (engine.getBoolean("-eC")) {
            endChr = Integer.parseInt(engine.getString("-eC"));
        }
        if (engine.getBoolean("-sD")) {
            startDiv = Integer.parseInt(engine.getString("-sD"));
        }
        if (engine.getBoolean("-eD")) {
            endDiv = Integer.parseInt(engine.getString("-eD"));
        }
        hmpFile = engine.getString("-hmp");
        outFileBase = engine.getString("-o");
        errFile = engine.getString("-oE");
        if (engine.getBoolean("-mxDiv")) {
            maxDistFromFounder = Double.parseDouble(engine.getString("-mxDiv"));
        }
        if (engine.getBoolean("-mxHet")) {
            maxHetFreq = Double.parseDouble(engine.getString("-mxHet"));
        }
        if (engine.getBoolean("-maxOutMiss")) {
            maximumMissing = Double.parseDouble(engine.getString("-maxOutMiss"));
        }
        if (engine.getBoolean("-hapSize")) {
            appoxSitesPerHaplotype = Integer.parseInt(engine.getString("-hapSize"));
        }
        if (engine.getBoolean("-minPres")) {
            minSitesPresentPerHap = Integer.parseInt(engine.getString("-minPres"));
        }
        if (engine.getBoolean("-maxHap")) {
            maxHaplotypes = Integer.parseInt(engine.getString("-maxHap"));
        }
    }



    private void printUsage() {
        myLogger.info(
                "\n\n\nAvailable options for the FindMergeHaplotypesPlugin are as follows:\n"
                + "-hmp   Input HapMap file (either hmp.txt.gz or hmp.h5)\n"
                + "-o     Output file(s) must include 's+.' plus will be replace by segment (0..(~sites/hapSize)\n"
                + "-oE  Optional file to record site by sites errors as the haplotypes are developed\n"
                + "-sC    Start chromosome\n"
                + "-eC    End chromosome\n"
                + "-mxDiv    Maximum divergence from founder haplotype\n"
                + "-mxHet    Maximum heterozygosity of haplotype to even scanned\n"
                + "-hapSize    Preferred haplotype block size in sites (minimum 64); will use the closest multiple of 64 at or below the supplied value\n"
                + "-minPres    Minimum number of present sites within input sequence to do the search\n"
                + "-maxHap    Maximum number of haplotypes per segment\n"
                + "-maxOutMiss  Maximum frequency of missing data in the output haplotype"
                );
    }

   @Override
    public DataSet performFunction(DataSet input) {
       if(outFileBase.contains(".gX.")) {
           runFindMergeHaplotypes(hmpFile, outFileBase, errFile, maxDistFromFounder, minSitesPresentPerHap, appoxSitesPerHaplotype);
           return null;
       }
       for (int chr = startChr; chr <=endChr; chr++) {
           String chrHmpFile=hmpFile.replace("chrX", "chr"+chr);
           chrHmpFile=chrHmpFile.replace("cX", "c"+chr);
           String chrOutfiles=outFileBase.replace("chrX", "chr"+chr);
           chrOutfiles=chrOutfiles.replace("cX", "c"+chr);
           runFindMergeHaplotypes(chrHmpFile, chrOutfiles, errFile, maxDistFromFounder, minSitesPresentPerHap, appoxSitesPerHaplotype);
       }     
       return null;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "ExtractInbredHaplotypes";
    }

    @Override
    public String getToolTipText() {
        return "Creates haplotype alignments based on long IBD regions of inbred lines";
    }
   
    
    /**
     *
     * @param args
     */
    public static void main(String[] args) {
       // String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/impOrig/";
//        String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/IMP26/orig/";
//        String rootOut="/Users/edbuckler/SolexaAnal/GBS/build20120701/IMP26/haplos/";
        String root="/Volumes/LaCie/build20120701/IMP26/orig/";
        String rootOut="/Volumes/LaCie/build20120701/IMP26/haplos/";

        
        String infileH5=root+"AllZeaGBS_v2.6.chrX.hmp.h5";
        String secFile26v=rootOut+"Tall26_8k.cXsX.hmp.txt.gz";

        String errorFile=rootOut+"mcErrorXMerge20130425.txt";

        String[] args2 = new String[]{
            "-hmp", infileH5,
            "-o", secFile26v,
            "-sC","8",
            "-eC","8",
            "-mxDiv", "0.01",
            "-mxHet", "0.01",
            "-hapSize", "8000",
            "-minPres", "500",
            "-maxOutMiss","0.4",
            "-maxHap", "2000",
        };

        FindMergeHaplotypesPlugin plugin = new FindMergeHaplotypesPlugin();
        plugin.setParameters(args2);
        plugin.performFunction(null);

    }
    
}
