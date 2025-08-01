/*
 * FILLINFindHaplotypesPlugin
 */
package net.maizegenetics.analysis.imputation;

import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.Utils;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Creates haplotypes by finding large IBS regions within GBS data.  Starts with the 
 * highest coverage taxa and looks within windows of near perfect matches.  Combines
 * all matches together into one haplotype.  The haplotype is named for the highest coverage
 * sample.  
 *
 * @author Ed Buckler
 * @author Kelly Swarts
 */
public class FILLINFindHaplotypesPlugin extends AbstractPlugin {
    //Plugin parameters
    private PluginParameter<String> hmpFile= new PluginParameter.Builder<>("hmp",null,String.class).guiName("Target file").inFile().required(true)
            .description("Input genotypes to generate haplotypes from. Usually best to use all available samples from a species. Accepts all file types supported by TASSEL5.").build();
    private PluginParameter<String> outFileBase= new PluginParameter.Builder<>("o",null,String.class).guiName("Donor dir/file basename").outDir().required(true)
            .description("Output file directory name, or new directory path; "
                    + "Directory will be created, if doesn't exist. Outfiles will be placed in the directory and given the same name and appended with the substring '.gc#s#.hmp.txt' to denote chromosome and section").build();

    private PluginParameter<Double> maxDistFromFounder= new PluginParameter.Builder<>("mxDiv",0.01,Double.class).guiName("Max divergence from founder")
            .description("Maximum genetic divergence from founder haplotype to cluster sequences").build();
    private PluginParameter<Double> maxHetFreq= new PluginParameter.Builder<>("mxHet",0.01,Double.class).guiName("Max heterozygosity of output haplotypes")
            .description("Maximum heterozygosity of output haplotype. Heterozygosity results from clustering sequences that either have residual heterozygosity "
                    + "or clustering sequences that do not share all minor alleles.").build();
    private PluginParameter<Integer> minSitesForSectionComp= new PluginParameter.Builder<>("minSites",50,Integer.class).guiName("Min sites to cluster")
            .description("The minimum number of sites present in two taxa to compare genetic distance to evaluate similarity for clustering").build();
    private PluginParameter<Double> maxErrorInCreatingConsensus= new PluginParameter.Builder<>("mxErr",0.05,Double.class).guiName("Max combined error to impute two donors")
            .description("The maximum genetic divergence allowable to cluster taxa").build();
    private PluginParameter<Integer> appoxSitesPerHaplotype= new PluginParameter.Builder<>("hapSize",8192,Integer.class).guiName("Preferred haplotype size")
            .description("Preferred haplotype block size in sites (minimum 64); will use the closest multiple of 64 at or below the supplied value").build();
    private PluginParameter<Integer> minSitesPresentPerHap= new PluginParameter.Builder<>("minPres",500,Integer.class).guiName("Min sites to test match")
            .description("Minimum number of present sites within input sequence to do the search").build();
    private PluginParameter<Integer> maxHaplotypes= new PluginParameter.Builder<>("maxHap",3000,Integer.class).guiName("Max haplotypes per segment")
            .description("Maximum number of haplotypes per segment").build();
    private PluginParameter<Integer> minTaxaInGroup= new PluginParameter.Builder<>("minTaxa",2,Integer.class).guiName("Min taxa to generate a haplotype")
            .description("Minimum number of taxa to generate a haplotype").build();
    private PluginParameter<Double> maximumMissing= new PluginParameter.Builder<>("maxOutMiss",0.4,Double.class).guiName("Max frequency missing per haplotype")
            .description("Maximum frequency of missing data in the output haplotype").build();
    private PluginParameter<Boolean> nonverboseOutput= new PluginParameter.Builder<>("nV",false,Boolean.class).guiName("Supress system out")
            .description("Supress system out").build();
    private PluginParameter<Boolean> extendedOutput= new PluginParameter.Builder<>("extOut",false,Boolean.class).guiName("Detailed system out on haplotypes")
            .description("Details of taxa included in each haplotype written to file").build();
    
    //other parameters
    private boolean verboseOutput= true;
    private int startDiv=-1, endDiv=-1;

    private double minJointGapProb=0.01;
    private boolean callGaps=false;
    private boolean anonymous= false;

    private double[] propMissing;
    private int[] siteErrors, siteCallCnt;
    private BitSet badMask=null;
    
    //these are for use with extended out. holds the new haplotypes and the taxon that go into them
    private TreeMap<String, ArrayList<String>> currHits= null; 
    private HashMap<String,ArrayList<String>> hapsForResults= null;
    private boolean uniqueHaps= true; //secret flag to make sure that all of the haplotypes are truely different from each other. For pop gen
    
    public FILLINFindHaplotypesPlugin() {
        super(null, false);
    }

    public  FILLINFindHaplotypesPlugin(Frame  parentFrame,  boolean  isInteractive)  {  
              super(parentFrame,  isInteractive);
    } 
    
    @Override
    protected void postProcessParameters() {
        System.out.println("Calling FILLINFindHaplotypesPlugin.postProcessParameters()...");
        if (nonverboseOutput.value()) verboseOutput= false;
        try {
//            if (new File(outFileBase.value()).isFile()) {
//                if (outFileBase.value().contains(".gX")) outputFilename(outputFilename().substring(0, outputFilename().indexOf(".gX")));
//                else throw new IOException();
//                new File(outFileBase.value()).mkdir();
//            }
            if (new File(outFileBase.value()).exists()) {
                if (!new File(outFileBase.value()).isDirectory()) throw new IOException();
            }
            else {
                new File(outFileBase.value()).mkdir();
                if (new File(outFileBase.value()).isFile()) {
                    new File(outFileBase.value()).delete();
                    throw new IOException();
                }
            }
        }
        catch (Exception e) {
            System.out.println("output directory must be an existing directory or new directory path");
        }
//        try {
//            if (outFileBase.value().contains("gX")==false) throw new IOException();
//                    }
//        catch (Exception e) {
//            System.out.println("output file name must contain gX, eg outfile.gX.hmp.txt");
//        }
//        if (!outFileBase.value().contains(".hmp")) outputFilename(outFileBase.value()+".hmp.txt");
    }
    
    @Override
    public String getCitation() {
        return "Swarts K, Li H, Romero Navarro JA, An D, Romay MC, Hearne S, Acharya C, Glaubitz JC, Mitchell S, Elshire RJ, Buckler ES, Bradbury PJ "
                + "(2014) "
                + "Novel methods to optimize genotypic imputation for low-coverage, next-generation sequence data in crop plants. "
                + "Plant Genome 7(3). doi:10.3835/plantgenome2014.05.0023";
    }
    
    @Override
    public DataSet processData(DataSet input) {
        System.out.println("Reading: "+hmpFile.value());
        GenotypeTable baseAlign=ImportUtils.read(hmpFile.value());
        int[][] divisions=divideChromosome(baseAlign, appoxSitesPerHaplotype.value(), verboseOutput);
        System.out.printf("In taxa:%d sites:%d %n",baseAlign.numberOfTaxa(),baseAlign.numberOfSites());
        siteErrors=new int[baseAlign.numberOfSites()];
        siteCallCnt=new int[baseAlign.numberOfSites()];
        if(startDiv==-1) startDiv=0;
        if(endDiv==-1) endDiv=divisions.length-1;
        if (extendedOutput.value()) new File(outFileBase.value()+"/extOut/").mkdir();
        for (int i = startDiv; i <=endDiv; i++) {
            String newExport= outFileBase.value()+"/"+new File(outFileBase.value()).getName()+".gc"+baseAlign.chromosomeName(divisions[i][0])+"s"+i+".hmp.txt";
            if (new File(newExport).exists()) continue;
            GenotypeTable mna=createHaplotypeAlignment(divisions[i][0], divisions[i][1], baseAlign,
             minSitesPresentPerHap.value(),  maxDistFromFounder.value());
            if (mna.taxa().isEmpty()) continue;
//            String newExport=outFileBase.value().replace("sX.hmp", "s"+i+".hmp");
//            newExport=newExport.replace("gX", "gc"+mna.chromosomeName(0)+"s"+i);
            if (extendedOutput.value()) writeHaplotypeMembersToFile(outFileBase.value()+"/extOut/"+new File(outFileBase.value()).getName()+".gc"+mna.chromosomeName(0)+"s"+i+".haplotypes.txt");
//            String newExport= outFileBase.value()+"/"+new File(outFileBase.value()).getName()+".gc"+mna.chromosomeName(0)+"s"+i+".hmp.txt";
            ExportUtils.writeToHapmap(mna, false, newExport, '\t', null);
            if(outFileBase.value()!=null) exportBadSites(baseAlign, outFileBase.value(), 0.01);  
            mna=null;
            System.gc();
        }
        return null;
     }
    
    private GenotypeTable createHaplotypeAlignment(int startSite, int endSite, GenotypeTable baseAlign,
            int minSites, double maxDistance) {
        GenotypeTable fa=FilterGenotypeTable.getInstance(baseAlign, startSite, endSite);
        GenotypeTable inAlign=GenotypeTableBuilder.getGenotypeCopyInstance(fa);
        int sites=inAlign.numberOfSites();
        if(verboseOutput) System.out.printf("SubInAlign Locus:%s StartPos:%d taxa:%d sites:%d %n",inAlign.chromosome(0),
                inAlign.chromosomalPosition(0),inAlign.numberOfTaxa(),inAlign.numberOfSites());

        propMissing=new double[inAlign.numberOfTaxa()];
        int startBlock=0;
        int lastBlock=inAlign.allelePresenceForAllSites(0, WHICH_ALLELE.Major).getNumWords()-1;
        TreeMap<Integer,Integer> presentRanking=createPresentRankingForWindow(inAlign, startBlock, lastBlock, minSites, maxHetFreq.value());
        if(verboseOutput) System.out.printf("\tBlock %d Inbred and modest coverage:%d %n",startBlock,presentRanking.size());
        if(verboseOutput) System.out.printf("\tCurrent Site %d Current block %d EndBlock: %d %n",startSite, startBlock, lastBlock);
        if (extendedOutput.value()) currHits= new TreeMap<>();
        TreeMap<Integer,byte[][]> results=mergeWithinWindow(inAlign, presentRanking, startBlock, lastBlock, maxDistance, startSite);
        TaxaListBuilder tLB=new TaxaListBuilder();
        GenotypeCallTableBuilder gB=GenotypeCallTableBuilder.getInstance(results.size(),inAlign.numberOfSites());
        int index=0;
        for (byte[][] calls : results.values()) {
            if (anonymous) tLB.add(new Taxon("h"+index));
            else tLB.add(new Taxon("h"+index+(new String(calls[1]))));
            gB.setBaseRangeForTaxon(index,0,calls[0]);
            if (extendedOutput.value()) currHits.put("h"+index+(new String(calls[1])), hapsForResults.get(new String(calls[1])));
            index++;
        }
        return GenotypeTableBuilder.getInstance(gB.build(),inAlign.positions(),tLB.build());
    }


    /**
     * Divides a genome up into roughly equally sized blocks of sites
     * @param a
     * @param appoxSitesPerHaplotype
     * @param verboseOutput
     * @return array of start and last positions
     */
    public static int[][] divideChromosome(GenotypeTable a, int appoxSitesPerHaplotype, boolean verboseOutput) {
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
            if(verboseOutput) System.out.printf("Chr:%s Alignment Sites:%d subAlignCnt:%d RealSites:%d %n",
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
            if(verboseOutput) System.out.printf("Chromosome Divisions: %s start:%d end:%d %n", a.chromosome(result[i][0]).getName(),
                    result[i][0], result[i][1]);
        }
        return result;
    }
    
    private TreeMap<Integer,Integer> createPresentRankingForWindow(GenotypeTable inAlign, int startBlock, int endBlock,
            int minSites, double maxHetFreq) {
        int sites=64*(endBlock-startBlock+1);
        TreeMap<Integer,Integer> presentRanking=new TreeMap<Integer,Integer>(Collections.reverseOrder());
        for (int i = 0; i < inAlign.numberOfTaxa(); i++) {
            long[] mj=inAlign.allelePresenceForSitesBlock(i, WHICH_ALLELE.Major, startBlock, endBlock+1);
            long[] mn=inAlign.allelePresenceForSitesBlock(i, WHICH_ALLELE.Minor, startBlock, endBlock+1);
            long[] mn2=inAlign.allelePresenceForSitesBlock(i, WHICH_ALLELE.Minor2, startBlock, endBlock+1);
            int totalSitesNotMissing = 0;
            int hetCnt=0;
            for (int j = 0; j < mj.length; j++) {
                totalSitesNotMissing+=BitUtil.pop(mj[j]|mn[j]|mn2[j]);
                hetCnt+=BitUtil.pop((mj[j]&mn[j]) | (mj[j]&mn2[j]) | (mn[j]&mn2[j]));
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


//    private BitSet maskBadSites(GeneticMap gm, GenotypeTable a) {
//        OpenBitSet obs=new OpenBitSet(a.numberOfSites());
//        int count=0;
//        for (int i = 0; i < gm.getNumberOfMarkers(); i++) {
//            int site=a.siteOfPhysicalPosition(gm.getPhysicalPosition(i), null);
//            if(site>0) {obs.set(site);}
//
//        }
//        System.out.println("Bad Sites matched:"+obs.cardinality());
//        obs.not();  //change all bad sites to 0, good to 1
//
//        return obs;
//    }
    
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
            int firstBlock, int lastBlock, double maxDistance, int siteOffsetForError){
        int startSite=firstBlock*64;
        int endSite=63+(lastBlock*64);
        if(endSite>=inAlign.numberOfSites()) endSite=inAlign.numberOfSites()-1;
        TreeMap<Integer,ArrayList> mergeSets=new TreeMap<Integer,ArrayList>();
        TreeMap<Integer,byte[][]> results=new TreeMap<Integer,byte[][]>(Collections.reverseOrder());
        if (extendedOutput.value()) hapsForResults= new HashMap<String,ArrayList<String>>();
        TreeSet<Integer> unmatched=new TreeSet<Integer>(presentRanking.values());
        TaxaList inIDG=inAlign.taxa();
        for (Entry<Integer,Integer> e : presentRanking.entrySet()) {
            int taxon1=e.getValue();
            if(unmatched.contains(taxon1)==false) continue;//already included in another group
            ArrayList<Integer> hits=new ArrayList<Integer>();
            unmatched.remove(taxon1);
            for(int taxon2 : unmatched) {
               double[] dist=IBSDistanceMatrix.computeHetBitDistances(inAlign, 
                       taxon1, taxon2, minSitesForSectionComp.value(), firstBlock, lastBlock, badMask);
               if((!Double.isNaN(dist[0]))&&(dist[0]<maxDistance)) {
                   hits.add(taxon2);
               }
            }
            byte[] calls=inAlign.genotypeRange(taxon1, startSite, endSite+1);//kls
            int[] unkCnt=countUnknown(calls);//kls
            double missingFreq=(double)unkCnt[0]/(double)inAlign.numberOfSites();//kls
            if(((hits.size()+1)<this.minTaxaInGroup.value())&&missingFreq>maximumMissing.value()) continue;//KLS changed this to not skip over taxa with no neighbors, but high coverage
            if(hits.size()>0) {
                ArrayList<String> mergeNames=new ArrayList<String>();
                mergeNames.add(inIDG.taxaName(taxon1));
                mergeSets.put(taxon1, hits);         
               // System.out.print(inAlign.getFullTaxaName(taxon1)+"=");
                for (Integer taxon2 : hits) {
                    unmatched.remove(taxon2);
                   // System.out.print(inAlign.getFullTaxaName(taxon2)+"=");
                    mergeNames.add(inIDG.taxaName(taxon2));
                }
              //  System.out.println("");              
                calls=consensusGameteCalls(inAlign, mergeNames, startSite, endSite, maxErrorInCreatingConsensus.value(), siteOffsetForError);
            } else {
                calls=inAlign.genotypeRange(taxon1, startSite, endSite+1);
            }
            unkCnt=countUnknown(calls);//kls
            missingFreq=(double)unkCnt[0]/(double)inAlign.numberOfSites();//kls
            double hetFreq=(double)unkCnt[1]/(double)(inAlign.numberOfSites()-unkCnt[0]);
            if(((missingFreq<maximumMissing.value())&&(hetFreq<maxHetFreq.value())&&(hits.size()>=(this.minTaxaInGroup.value()-1)))) {
                int index=(hits.size()*200000)+taxon1;
                if(verboseOutput) System.out.printf("\t\tOutput %s plus %d missingF:%g hetF:%g index: %d %n",inIDG.taxaName(taxon1),
                        hits.size(), missingFreq, hetFreq, index);
                if(extendedOutput.value()) {
                    ArrayList<String> hitNames= new ArrayList<>(); hitNames.add(inIDG.taxaName(taxon1));
                    for (Integer taxon:hits) {
                        hitNames.add(inIDG.taxaName(taxon));
                    }
                    hapsForResults.put(inIDG.taxaName(taxon1)+":d"+(hits.size()+1), hitNames);
                }
                byte[][] callPlusNames=new byte[2][];
                callPlusNames[0]=calls;
                String newName=inIDG.taxaName(taxon1)+":d"+(hits.size()+1);
                callPlusNames[1]=newName.getBytes();
                results.put(index, callPlusNames);
            }
            if(results.size()>=maxHaplotypes.value()) break;
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
            boolean third= a.minorAlleles(s).length>1?true:false;
            byte mjAllele = a.majorAllele(s);
            byte mnAllele = a.minorAllele(s);
            byte mj=GenotypeTableUtils.getUnphasedDiploidValue(mjAllele,mjAllele);
            byte mn=GenotypeTableUtils.getUnphasedDiploidValue(mnAllele,mnAllele);
            byte het1 = GenotypeTableUtils.getUnphasedDiploidValue(mjAllele, mnAllele);
            byte mn2= GenotypeTable.UNKNOWN_DIPLOID_ALLELE; 
            byte het2= GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
            byte het3= GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
            if (third) {
                byte mn2Allele = a.minorAlleles(s)[1];
                mn2=GenotypeTableUtils.getUnphasedDiploidValue(mn2Allele,mn2Allele);
                het2 = GenotypeTableUtils.getUnphasedDiploidValue(mjAllele, mn2Allele);
                het3 = GenotypeTableUtils.getUnphasedDiploidValue(mnAllele, mn2Allele);
            }
            int mjCnt=0, mnCnt=0, mn2Cnt=0;
            for (int t = 0; t < taxaIndex.length; t++) {
                byte ob = a.genotype(taxaIndex[t], s);
                if (ob == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    continue;
                }
                if (ob == mj) {
                    mjCnt+=2;
                } else if (ob == mn) {
                    mnCnt+=2;
                } else if (GenotypeTableUtils.isEqual(ob, het1)) {
                    mjCnt++;
                    mnCnt++;
                } else if (third && ob == mn2) {
                    mn2Cnt+=2;
                } else if (third && GenotypeTableUtils.isEqual(ob, het2)) {
                    mjCnt++;
                    mn2Cnt++;
                } else if (third && GenotypeTableUtils.isEqual(ob, het3)) {
                    mnCnt++;
                    mn2Cnt++;
                }
            }
            int totalCnt = mjCnt + mnCnt + mn2Cnt;
            
            if (totalCnt == 0) {
                double missingProp=1.0;
                for (int t : taxaIndex) {missingProp*=propMissing[t];}
                if(callGaps&(missingProp<minJointGapProb)) calls[s-startSite]=NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE;
                continue;
            }
            double errorRate;
            if(totalCnt>1) siteCallCnt[s+siteOffsetForError]+=totalCnt;
            if((mjCnt+mn2Cnt) < mnCnt) {
                errorRate=(double)(mjCnt+mn2Cnt)/(double)totalCnt;
                if(errorRate<maxError) {calls[s-startSite] = mn;}
                else {siteErrors[s+siteOffsetForError]+=(mjCnt+mn2Cnt);}
            } else if ((mnCnt+mn2Cnt) < mjCnt) {
                errorRate=(double)(mnCnt+mn2Cnt)/(double)totalCnt;
                if(errorRate<maxError) {calls[s-startSite] = mj;}
                else {siteErrors[s+siteOffsetForError]+=(mnCnt+mn2Cnt);}
            } else if ((mjCnt+mnCnt) < mn2Cnt) {
                errorRate=(double)(mjCnt+mnCnt)/(double)totalCnt;
                if(errorRate<maxError) {calls[s-startSite] = mn2;}
                else {siteErrors[s+siteOffsetForError]+=(mjCnt+mnCnt);}
            } else if (uniqueHaps) {
                if (mjCnt==mnCnt && mjCnt!=0) calls[s-startSite]= het1;
                else if (mjCnt==mn2Cnt && mjCnt!=0) calls[s-startSite]= het2;
                else if (mnCnt==mn2Cnt && mnCnt!=0) calls[s-startSite]= het3;
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

    
    private boolean writeHaplotypeMembersToFile(String outFileName) {
        if (currHits==null) {System.out.println("No haplotypes to write to "+outFileName); return false;}
        DataOutputStream outStream = Utils.getDataOutputStream(outFileName, 1040);
        try {
            for (String hap:currHits.keySet()) {
                outStream.writeBytes(hap);
                for (String taxon: currHits.get(hap)) {
                    outStream.writeBytes("\t"+taxon);
                }
                outStream.writeBytes("\n");
            }
            outStream.close();
        } catch(Exception e) {
            System.out.println("Problem writing haplotype information file: "+e);
        }
        currHits= null; hapsForResults= null;
        return true;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Extract Inbred Haplotypes by FILLIN";
    }

    @Override
    public String getToolTipText() {
        return "Creates haplotype alignments based on long IBD regions of inbred lines";
    }
    
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(FILLINFindHaplotypesPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    public Boolean runPlugin(DataSet input) {
        performFunction(input);
        return false;
//        return (Boolean) performFunction(input).getData(0).getData();
    }

    /**
     * Input HapMap file of target genotypes to impute or
     * all genotypes available for target taxon. Accepts all
     * file types supported by TASSEL5.
     *
     * @return Target file
     */
    public String targetFile() {
        return hmpFile.value();
    }

    /**
     * Set Target file. Input HapMap file of target genotypes
     * to impute or all genotypes available for target taxon.
     * Accepts all file types supported by TASSEL5.
     *
     * @param value Target file
     *
     * @return this plugin
     */
    public FILLINFindHaplotypesPlugin targetFile(String value) {
        hmpFile = new PluginParameter<>(hmpFile, value);
        return this;
    }

    /**
     * Output file; Use .gX in the output filename to denote
     * the substring .gc#s# found in donor files, ie 'out.gX.hmp.txt'.
     *
     * @return Output filename
     */
    public String outputFilename() {
        return outFileBase.value();
    }

    /**
     * Set Output filename. Output file; Use .gX in the output
     * filename to denote the substring .gc#s# found in donor
     * files, ie 'out.gX.hmp.txt'.
     *
     * @param value Output filename
     *
     * @return this plugin
     */
    public FILLINFindHaplotypesPlugin outputFilename(String value) {
        outFileBase = new PluginParameter<>(outFileBase, value);
        return this;
    }

    /**
     * Maximum genetic divergence from founder haplotype to
     * cluster sequences
     *
     * @return Max divergence from founder
     */
    public Double maxDivergenceFromFounder() {
        return maxDistFromFounder.value();
    }

    /**
     * Set Max divergence from founder. Maximum genetic divergence
     * from founder haplotype to cluster sequences
     *
     * @param value Max divergence from founder
     *
     * @return this plugin
     */
    public FILLINFindHaplotypesPlugin maxDivergenceFromFounder(Double value) {
        maxDistFromFounder = new PluginParameter<>(maxDistFromFounder, value);
        return this;
    }

    /**
     * Maximum heterozygosity of output haplotype. Heterozygosity
     * results from clustering sequences that either have
     * residual heterozygosity or clustering sequences that
     * do not share all minor alleles.
     *
     * @return Max heterozygosity of output haplotypes
     */
    public Double maxHeterozygosityOfOutputHaplotypes() {
        return maxHetFreq.value();
    }

    /**
     * Set Max heterozygosity of output haplotypes. Maximum
     * heterozygosity of output haplotype. Heterozygosity
     * results from clustering sequences that either have
     * residual heterozygosity or clustering sequences that
     * do not share all minor alleles.
     *
     * @param value Max heterozygosity of output haplotypes
     *
     * @return this plugin
     */
    public FILLINFindHaplotypesPlugin maxHeterozygosityOfOutputHaplotypes(Double value) {
        maxHetFreq = new PluginParameter<>(maxHetFreq, value);
        return this;
    }

    /**
     * The minimum number of sites present in two taxa to
     * compare genetic distance to evaluate similarity for
     * clustering
     *
     * @return Min sites to cluster
     */
    public Integer minSitesToCluster() {
        return minSitesForSectionComp.value();
    }

    /**
     * Set Min sites to cluster. The minimum number of sites
     * present in two taxa to compare genetic distance to
     * evaluate similarity for clustering
     *
     * @param value Min sites to cluster
     *
     * @return this plugin
     */
    public FILLINFindHaplotypesPlugin minSitesToCluster(Integer value) {
        minSitesForSectionComp = new PluginParameter<>(minSitesForSectionComp, value);
        return this;
    }

    /**
     * The maximum genetic divergence allowable to cluster
     * taxa
     *
     * @return Max combined error to impute two donors
     */
    public Double maxCombinedErrorToImputeTwoDonors() {
        return maxErrorInCreatingConsensus.value();
    }

    /**
     * Set Max combined error to impute two donors. The maximum
     * genetic divergence allowable to cluster taxa
     *
     * @param value Max combined error to impute two donors
     *
     * @return this plugin
     */
    public FILLINFindHaplotypesPlugin maxCombinedErrorToImputeTwoDonors(Double value) {
        maxErrorInCreatingConsensus = new PluginParameter<>(maxErrorInCreatingConsensus, value);
        return this;
    }

    /**
     * Preferred haplotype block size in sites (minimum 64);
     * will use the closest multiple of 64 at or below the
     * supplied value
     *
     * @return Preferred haplotype size
     */
    public Integer preferredHaplotypeSize() {
        return appoxSitesPerHaplotype.value();
    }

    /**
     * Set Preferred haplotype size. Preferred haplotype block
     * size in sites (minimum 64); will use the closest multiple
     * of 64 at or below the supplied value
     *
     * @param value Preferred haplotype size
     *
     * @return this plugin
     */
    public FILLINFindHaplotypesPlugin preferredHaplotypeSize(Integer value) {
        appoxSitesPerHaplotype = new PluginParameter<>(appoxSitesPerHaplotype, value);
        return this;
    }

    /**
     * Minimum number of present sites within input sequence
     * to do the search
     *
     * @return Min sites to test match
     */
    public Integer minSitesToTestMatch() {
        return minSitesPresentPerHap.value();
    }

    /**
     * Set Min sites to test match. Minimum number of present
     * sites within input sequence to do the search
     *
     * @param value Min sites to test match
     *
     * @return this plugin
     */
    public FILLINFindHaplotypesPlugin minSitesToTestMatch(Integer value) {
        minSitesPresentPerHap = new PluginParameter<>(minSitesPresentPerHap, value);
        return this;
    }

    /**
     * Maximum number of haplotypes per segment
     *
     * @return Max haplotypes per segment
     */
    public Integer maxHaplotypesPerSegment() {
        return maxHaplotypes.value();
    }

    /**
     * Set Max haplotypes per segment. Maximum number of haplotypes
     * per segment
     *
     * @param value Max haplotypes per segment
     *
     * @return this plugin
     */
    public FILLINFindHaplotypesPlugin maxHaplotypesPerSegment(Integer value) {
        maxHaplotypes = new PluginParameter<>(maxHaplotypes, value);
        return this;
    }

    /**
     * Minimum number of taxa to generate a haplotype
     *
     * @return Min taxa to generate a haplotype
     */
    public Integer minTaxaToGenerateAHaplotype() {
        return minTaxaInGroup.value();
    }

    /**
     * Set Min taxa to generate a haplotype. Minimum number
     * of taxa to generate a haplotype
     *
     * @param value Min taxa to generate a haplotype
     *
     * @return this plugin
     */
    public FILLINFindHaplotypesPlugin minTaxaToGenerateAHaplotype(Integer value) {
        minTaxaInGroup = new PluginParameter<>(minTaxaInGroup, value);
        return this;
    }

    /**
     * Maximum frequency of missing data in the output haplotype
     *
     * @return Max frequency missing per haplotype
     */
    public Double maxFrequencyMissingPerHaplotype() {
        return maximumMissing.value();
    }

    /**
     * Set Max frequency missing per haplotype. Maximum frequency
     * of missing data in the output haplotype
     *
     * @param value Max frequency missing per haplotype
     *
     * @return this plugin
     */
    public FILLINFindHaplotypesPlugin maxFrequencyMissingPerHaplotype(Double value) {
        maximumMissing = new PluginParameter<>(maximumMissing, value);
        return this;
    }

    /**
     * Supress system out
     *
     * @return Express system out
     */
    public Boolean supressSystemOut() {
        return nonverboseOutput.value();
    }

    /**
     * Set Supress system out. Express system out
     *
     * @param value Express system out
     *
     * @return this plugin
     */
    public FILLINFindHaplotypesPlugin supressSystemOut(Boolean value) {
        nonverboseOutput = new PluginParameter<>(nonverboseOutput, value);
        return this;
    } 
    
    /**
     * Details of taxa included in each haplotype to system
     * out
     *
     * @return Detailed system out on haplotypes
     */
    public Boolean detailedSystemOutOnHaplotypes() {
        return extendedOutput.value();
    }

    /**
     * Set Detailed system out on haplotypes. Details of taxa
     * included in each haplotype to system out
     *
     * @param value Detailed system out on haplotypes
     *
     * @return this plugin
     */
    public FILLINFindHaplotypesPlugin detailedSystemOutOnHaplotypes(Boolean value) {
        extendedOutput = new PluginParameter<>(extendedOutput, value);
        return this;
    }
}
