package net.maizegenetics.analysis.imputation;

import com.google.common.primitives.Bytes;
import com.google.common.primitives.Ints;
import net.maizegenetics.analysis.popgen.DonorHypoth;
import net.maizegenetics.dna.map.DonorHaplotypes;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.dna.snp.io.ProjectionGenotypeIO;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import static net.maizegenetics.dna.snp.GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
import static net.maizegenetics.dna.WHICH_ALLELE.Major;
import static net.maizegenetics.dna.WHICH_ALLELE.Minor;
import static net.maizegenetics.dna.snp.GenotypeTableUtils.isHeterozygous;
import static net.maizegenetics.dna.snp.NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;



/**
 * FILLIN imputation relies on a libary of haplotypes and uses nearest neighbor searches
 * followed by HMM Viterbi resolution or block-based resolution.
 * It is the best approach for substantially unrelated taxa in TASSEL.  BEAGLE4 is a better
 * approach currently for landraces, while FILLIN outperforms if there is a good reference
 * set of haplotypes.
 * <p></p>
 * The algorithm relies on having a series of donor haplotypes.  If phased haplotypes are
 * already known, they can be used directly.  If they need to be reconstructed, the
 * {@link FindMergeHaplotypesPlugin} can be used to create haplotypes for windows
 * across the genome.
 * <p></p>
 * Imputation is done one taxon at a time using the donor haplotypes.  The strategy is as
 *  follows:
 * <li></li>Every 64 sites is considered a block (or a word in the bitset terminology - length of a long).  There is
 * a separate initial search for nearest neighbors centered around each block. The flanks
 * of the scan window vary but always contain the number of minor alleles specified (a key parameter).
 * Minor alleles approximate the information content of the search.
 * <li></li> Calculate distance between the target taxon and the donor haplotypes for every window, and rank
 * the top 10 donor haplotypes for each 64 site focus block.
 * <li></li> Evaluate by Viterbi whether 1 or 2 haplotypes will explain all the sites across the
 * donor haplotype window.  If successful, move to the next region for donor haplotypes.
 * <li></li> Resolve each focus block by nearest neighbors.  If inbred NN are not close enough,
 * then do a hybrid search seeded with one parent being the initial 10 haplotypes.
 * Set bases based on what is resolved first inbred or hybrid.
 *
 *<p></p>
 * Error rates are bounded away from zero, but adding 0.5 error to all error
 * rates that that were observed to be zero.
 * <p></p>

 * @author Edward Buckler
 * @author Kelly Swarts
 * @cite
 */
//@Citation("Two papers: Viterbi from Bradbury, et al (in prep) Recombination patterns in maize\n"+
//        "NearestNeighborSearch Swarts,...,Buckler (in prep) Imputation with large genotyping by sequencing data\n")
public class FILLINImputationPlugin extends AbstractPlugin {
    private int minMajorRatioToMinorCnt=10;  //refinement of minMinorCnt to account for regions with all major
    //Plugin parameters
    private PluginParameter<String> hmpFile= new PluginParameter.Builder<>("hmp",null,String.class).guiName("Target file").inFile().required(true)
            .description("Input HapMap file of target genotypes to impute. Accepts all file types supported by TASSEL5. This file should include _ALL_ available sites, not " +
                    "just those segregating in your material. (ie: don't filter the input)").build();
    private PluginParameter<String> donorFile= new PluginParameter.Builder<>("d",null,String.class).guiName("Donor").required(true)
            .description("Directory containing donor haplotype files from output of FILLINFindHaplotypesPlugin. All files with '.gc' in the filename will be read in, "
                    + "only those with matching sites are used. Alternately, a single file to use as a donor, will be cut into sub genos in size specified (eg, high density"
                    + "SNP file for projection").build();
    private PluginParameter<String> outFileBase= new PluginParameter.Builder<>("o",null,String.class).guiName("Output filename").outFile().required(true)
            .description("Output file; hmp.txt.gz and .hmp.h5 accepted.").build();
    
    private PluginParameter<Integer> appoxSitesPerDonorGenotypeTable= new PluginParameter.Builder<>("hapSize",8000,Integer.class).guiName("Preferred haplotype size")
            .description("Preferred haplotype block size in sites (use same as in FILLINFindHaplotypesPlugin)").build();
    private PluginParameter<Double> hetThresh= new PluginParameter.Builder<>("mxHet",0.02,Double.class).guiName("Heterozygosity threshold")
            .description("Threshold per taxon heterozygosity for treating taxon as heterozygous (no Viterbi, het thresholds).").build();
    private PluginParameter<Double> maximumInbredError= new PluginParameter.Builder<>("mxInbErr",0.01,Double.class).guiName("Max error to impute one donor")
            .description("Maximum error rate for applying one haplotype to entire site window").build();
    private PluginParameter<Double> maxHybridErrorRate= new PluginParameter.Builder<>("mxHybErr",0.003,Double.class).guiName("Max combined error to impute two donors")
            .description("Maximum error rate for applying Viterbi with two haplotypes to entire site window").build();
    private PluginParameter<Integer> minTestSites= new PluginParameter.Builder<>("mnTestSite",20,Integer.class).guiName("Min sites to test match")
            .description("Minimum number of sites to test for IBS between haplotype and target in focus block").build();
    private PluginParameter<Integer> minMinorCnt= new PluginParameter.Builder<>("minMnCnt",20,Integer.class).guiName("Min num of minor alleles to compare")
            .description("Minimum number of informative minor alleles in the search window (or "+minMajorRatioToMinorCnt+"X major)").build();
    private PluginParameter<Integer> maxDonorHypotheses= new PluginParameter.Builder<>("mxDonH",20,Integer.class).guiName("Max donor hypotheses")
            .description("Maximum number of donor hypotheses to be explored. (Note: For heterozygous samples you may want to set this number higher, but the computation load goes up quickly.)").build();
    private PluginParameter<Boolean> imputeAllHets= new PluginParameter.Builder<>("impAllHets",false,Boolean.class).guiName("Impute all het calls")
            .description("Write all imputed heterozygous calls as such, even if the original file has a homozygous call. (Not recommended for inbred lines.)").build();


        //if true, uses combination mode in focus block, else set don't impute (default is true)
    private PluginParameter<Boolean> hybridNN= new PluginParameter.Builder<>("hybNN",true,Boolean.class).guiName("Combine two haplotypes as heterozygote")
            .description("If true, uses combination mode in focus block, else does not impute").build();
    private PluginParameter<Boolean> isOutputProjection= new PluginParameter.Builder<>("ProjA",false,Boolean.class).guiName("Output projection alignment")
            .description("Create a projection alignment for high density markers").build();
    private PluginParameter<Boolean> imputeDonorFile= new PluginParameter.Builder<>("impDonor",false,Boolean.class).guiName("Impute donor file")
            .description("Impute the donor file itself").build();
    private PluginParameter<Boolean> nonverboseOutput= new PluginParameter.Builder<>("nV",false,Boolean.class).guiName("Supress system out")
            .description("Supress system out").build();
    
            //for calculating accuracy
    private PluginParameter<Boolean> accuracy= new PluginParameter.Builder<>("accuracy",false,Boolean.class).guiName("Calculate accuracy")
            .description("Masks input file before imputation and calculates accuracy based on masked genotypes").build();
    private PluginParameter<Double> propSitesMask= new PluginParameter.Builder<>("propSitesMask",0.01,Double.class).guiName("Proportion of genotypes to mask if no depth")
            .description("Proportion of genotypes to mask for accuracy calculation if depth not available").dependentOnParameter(accuracy).build();
    private PluginParameter<Integer> depthToMask= new PluginParameter.Builder<>("depthMask",9,Integer.class).guiName("Depth of genotypes to mask")
            .description("Depth of genotypes to mask for accuracy calculation if depth information available").dependentOnParameter(accuracy).build();
    private PluginParameter<Double> propDepthSitesMask= new PluginParameter.Builder<>("propDepthSitesMask",0.2,Double.class).guiName("Proportion of depth genotypes to mask")
            .description("Proportion of genotypes of given depth to mask for accuracy calculation if depth available").dependentOnParameter(accuracy).build();
    private PluginParameter<String> maskKey= new PluginParameter.Builder<>("maskKey",null,String.class).inFile().required(false).guiName("Optional key to calculate accuracy")
            .description("Key to calculate accuracy. Genotypes missing (masked) in target file should be present in key, with all other sites set to missing. Overrides other"
                    + "accuracy options if present and all sites and taxa present in target file present").dependentOnParameter(accuracy).build();
    private PluginParameter<Boolean> byMAF= new PluginParameter.Builder<>("byMAF",false,Boolean.class).guiName("Calculate accuracy within MAF categories")
            .description("Calculate R2 accuracy within MAF categories based on donor file").dependentOnParameter(accuracy).build();
    
    //Additional variables
    private boolean verboseOutput= true;
    private boolean twoWayViterbi= true;//if true, the viterbi runs in both directions (the longest path length wins, if inconsistencies)
    private double minimumDonorDistance=maximumInbredError.value()*5; //used to prevent Viterbi errors with too similar sequences
    private double maxNonMedelian=maximumInbredError.value()*5; //if used avoid Viterbi if too many sites are non-Mendelian

        //options for focus blocks
    private double maxHybridErrFocusHomo= .3333*maxHybridErrorRate.value();////max error rate for discrepacy between two haplotypes for the focus block. it's default is higher because calculating for fewer sites
    private double maxInbredErrFocusHomo= .3*maximumInbredError.value();//.003;
    private double maxSmashErrFocusHomo= maximumInbredError.value();//.01;
    private double maxInbredErrFocusHet= .1*maximumInbredError.value();//.001;//the error rate for imputing one haplotype in focus block for a het taxon
    private double maxSmashErrFocusHet= maximumInbredError.value();//.01;

    public static GenotypeTable unimpAlign;  //the unimputed alignment to be imputed, unphased
    public static GenotypeTable maskKeyAlign= null;  //the key to the unimputed alignment mask
    public static double[] MAFClass= new double[]{0,.02,.05,.10,.20,.3,.4,.5,1};//must retain 0 and 1
    private int testing=0;  //level of reporting to stdout
        //major and minor alleles can be differ between the donor and unimp alignment
    private boolean isSwapMajorMinor=true;  //if swapped try to fix it

    private boolean resolveHetIfUndercalled=false;//if true, sets genotypes called to a het to a het, even if already a homozygote

    //initialize the transition matrix (T1)
    double[][] transition = new double[][] {
            {.999,.0001,.0003,.0001,.0005},
            {.0002,.999,.00005,.00005,.0002},
            {.0002,.00005,.999,.00005,.0002},
            {.0002,.00005,.00005,.999,.0002},
            {.0005,.0001,.0003,.0001,.999}
    };

    //initialize the emission matrix, states (5) in rows, observations (3) in columns
    double[][] emission = new double[][] {
                    {.998,.001,.001},
                    {.6,.2,.2},
                    {.4,.2,.4},
                    {.2,.2,.6},
                    {.001,.001,.998}
    };
    FILLINImputationAccuracy acc= null; //holds the accuracy information if accuracy flagged


    private static final Logger myLogger = LogManager.getLogger(FILLINImputationPlugin.class);

    @Override
    protected void postProcessParameters() {
        System.out.println("Calling FILLINImputationPlugin.postProcessParameters()...");
        File parent = new File(outputFilename()).getParentFile();
        if ((parent != null) && (!parent.exists())) {
            parent.mkdirs();
        }
        minimumDonorDistance=maximumInbredError.value()*5;
        maxNonMedelian=maximumInbredError.value()*5; 
        maxInbredErrFocusHomo= .3*maximumInbredError.value();//.003;
        maxSmashErrFocusHomo= maximumInbredError.value();//.01;
        maxInbredErrFocusHet= .1*maximumInbredError.value();//.001;//the error rate for imputing one haplotype in focus block for a het taxon
        maxSmashErrFocusHet= maximumInbredError.value();//.01;
        maxHybridErrFocusHomo= .3333*maxHybridErrorRate.value();
        resolveHetIfUndercalled = imputeAllHets.value();
        if (nonverboseOutput.value()) verboseOutput= false;
        if (byMAF.value()==false) MAFClass= null;
        if (!new File(donorFile.value()).exists()) System.out.println("Donor is not a file or folder: "+donorFile.value());
    }
    
    public FILLINImputationPlugin() {
        super(null, false);
    }

    public FILLINImputationPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
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
        long time=System.currentTimeMillis();
        String donor= donorFile.value();
        unimpAlign=ImportUtils.read(hmpFile.value());
        if (isOutputProjection.value()) {
            unimpAlign= FILLINDonorGenotypeUtils.RemoveSitesThatDoNotMatchMinMaj(donorFile.value(), unimpAlign,verboseOutput);
            donor= donor.subSequence(0, donor.lastIndexOf("."))+".matchMinMaj.hmp.txt.gz";
        }
        GenotypeTable[] donorAlign=FILLINDonorGenotypeUtils.loadDonors(donor, unimpAlign, minTestSites.value(),
                verboseOutput,appoxSitesPerDonorGenotypeTable.value());
        if (accuracy.value()) {
            time= System.currentTimeMillis()-time; //holds the time so far
            if (maskKey.value()!=null) maskKeyAlign= ImportUtils.readGuessFormat(maskKey.value());
            acc= new FILLINImputationAccuracy(unimpAlign,maskKeyAlign,donorAlign,propSitesMask.value(),depthToMask.value(), 
                propDepthSitesMask.value(),outFileBase.value(),MAFClass,verboseOutput);
            unimpAlign= acc.initiateAccuracy();
            time= System.currentTimeMillis()-time;//restarts the time, including the time to load donors and target but not including the time to set up accuracy
        }
        OpenBitSet[][] conflictMasks=FILLINDonorGenotypeUtils.createMaskForAlignmentConflicts(unimpAlign, donorAlign,
                verboseOutput);

        System.out.printf("Unimputed taxa:%d sites:%d %n",unimpAlign.numberOfTaxa(),unimpAlign.numberOfSites());
        System.out.println("Creating Export GenotypeTable:"+outFileBase.value());
        Object mna;
        if(isOutputProjection.value()) {
            mna=new ProjectionBuilder(ImportUtils.readGuessFormat(donorFile.value()));
        } else {
            if(outFileBase.value().contains(".h5")) {
                mna= GenotypeTableBuilder.getTaxaIncremental(this.unimpAlign.positions(),outFileBase.value());
            }else {
                mna= GenotypeTableBuilder.getTaxaIncremental(this.unimpAlign.positions());
            }
        }
        int numThreads = Runtime.getRuntime().availableProcessors();
        System.out.println("Time to read in files and generate masks: "+((System.currentTimeMillis()-time)/1000)+" sec");
        ExecutorService pool = Executors.newFixedThreadPool(numThreads);
        
        for (int taxon = 0; taxon < unimpAlign.numberOfTaxa(); taxon+=1) {
            int[] trackBlockNN= new int[5];//global variable to track number of focus blocks solved in NN search for system out; index 0 is inbred, 1 is viterbi, 2 is smash, 3 is not solved, 4 is total for all modes
            ImputeOneTaxon theTaxon= (((double)unimpAlign.heterozygousCountForTaxon(taxon)/(double)unimpAlign.totalNonMissingForTaxon(taxon))<hetThresh.value())?
                new ImputeOneTaxon(taxon, donorAlign, minTestSites.value(), conflictMasks,imputeDonorFile.value(), mna, trackBlockNN, maxInbredErrFocusHomo, maxHybridErrFocusHomo, maxSmashErrFocusHomo, true):
                    new ImputeOneTaxon(taxon, donorAlign, minTestSites.value(), conflictMasks,imputeDonorFile.value(), mna, trackBlockNN, maxInbredErrFocusHet, 0, maxSmashErrFocusHet, false);
            //        theTaxon.run(); //retained to provide a quick way to debug.  uncomment to help in debugging.
            pool.execute(theTaxon);
        }
        pool.shutdown();
        try{
            if (!pool.awaitTermination(48, TimeUnit.HOURS)) {
                System.out.println("processing threads timed out.");
            }
        }catch(Exception e) {
            System.out.println("Error processing threads");
        }
        System.out.println("");
        StringBuilder s=new StringBuilder();
        s.append(String.format("%s %s MinMinor:%d ", donorFile.value(), hmpFile.value(), minMinorCnt.value()));
        System.out.println(s.toString());
        
        double runtime= (double)(System.currentTimeMillis()-time)/(double)1000;
        System.out.printf("%d %g %d %n",minMinorCnt.value(), maximumInbredError.value(), maxDonorHypotheses.value());
        System.out.println("Runtime: "+runtime+" seconds");
        GenotypeTable out= null;
        if(isOutputProjection.value()) {
            out= ((ProjectionBuilder) mna).build();
            ProjectionGenotypeIO.writeToFile(outFileBase.value(), out);
        } else {
            GenotypeTableBuilder ab=(GenotypeTableBuilder)mna;
            ab.sortTaxa();
            out= ab.build();
            if(ab.isHDF5()) {
            } else {
                ExportUtils.writeToHapmap(out, false, outFileBase.value(), '\t', null);
            }
            if (accuracy.value()) {
                acc.calcAccuracy(ImportUtils.readGuessFormat(outFileBase.value()), runtime);
            }

        }
        return new DataSet(new Datum("outFile",out,null),null);
    }

    private class ImputeOneTaxon implements Runnable{
        int taxon;
        GenotypeTable[] donorAlign;
        int minSitesPresent;
        OpenBitSet[][] conflictMasks;
        boolean imputeDonorFile;
        GenotypeTableBuilder alignBuilder=null;
        ProjectionBuilder projBuilder=null;

        int[] trackBlockNN;//global variable to track number of focus blocks solved in NN search for system out; index 0 is inbred, 1 is viterbi, 2 is smash, 3 is not solved, 4 is total for all modes
        double focusInbredErr; //threshold for one haplotype imputation in focus block mode
        double focusHybridErr; //threshold for Viterbi imputation in focus block mode
        double focusSmashErr; //threshold for haplotype combination in focus block mode
        boolean hetsMiss; //for inbred lines in two haplotype combination, set hets to missing because likely error. for heterozygous, impute estimated hets in focus block mode
        
        public ImputeOneTaxon(int taxon, GenotypeTable[] donorAlign, int minSitesPresent, OpenBitSet[][] conflictMasks,
            boolean imputeDonorFile, Object mna, int[] trackBlockNN, double focusInbErr, double focusHybridErr, double focusSmashErr, boolean hetsToMissing) {
            this.taxon=taxon;
            this.donorAlign=donorAlign;
            this.minSitesPresent=minSitesPresent;
            this.conflictMasks=conflictMasks;
            this.imputeDonorFile=imputeDonorFile;
            this.trackBlockNN=trackBlockNN;
            this.focusInbredErr=focusInbErr;
            this.focusHybridErr=focusHybridErr;
            this.focusSmashErr=focusSmashErr;
            this.hetsMiss=hetsToMissing;
            if(mna instanceof GenotypeTableBuilder) {alignBuilder=(GenotypeTableBuilder)mna;}
            else if(mna instanceof ProjectionBuilder) {projBuilder=(ProjectionBuilder)mna;}
            else {throw new IllegalArgumentException("Only Aligmnent or Projection Builders may be used.");}
            
        }



        @Override
        public void run() {
            StringBuilder sb=new StringBuilder();
            String name=unimpAlign.taxaName(taxon);
            ImputedTaxon impTaxon=new ImputedTaxon(taxon, unimpAlign.genotypeAllSites(taxon),isOutputProjection.value());
            boolean het= (focusHybridErr==0)?true:false;
            int[] unkHets=FILLINImputationUtils.countUnknownAndHeterozygotes(impTaxon.getOrigGeno());
            sb.append(String.format("Imputed %d:%s AsHet:%b Mj:%d, Mn:%d Unk:%d Hets:%d... ", taxon,name,het,
                    unimpAlign.allelePresenceForAllSites(taxon, Major).cardinality(),
                    unimpAlign.allelePresenceForAllSites(taxon, Minor).cardinality(), unkHets[0], unkHets[1]));
            boolean enoughData=(unimpAlign.totalNonMissingForTaxon(taxon)>minSitesPresent);
            int countFullLength=0;
            int countByFocus= 0;
            for (int da = 0; (da < donorAlign.length)&&enoughData ; da++) {
                int donorOffset=unimpAlign.siteOfPhysicalPosition(donorAlign[da].chromosomalPosition(0), donorAlign[da].chromosome(0));
                int blocks=donorAlign[da].allelePresenceForAllSites(0, Major).getNumWords();
                BitSet[] maskedTargetBits=FILLINDonorGenotypeUtils.arrangeMajorMinorBtwAlignments(unimpAlign, taxon, donorOffset, donorAlign[da].numberOfSites(), conflictMasks[da][0], conflictMasks[da][1], isSwapMajorMinor);

                //if imputing the donor file, these donor indices prevent self imputation
                int[] donorIndices;
                if(imputeDonorFile){
                    donorIndices=new int[donorAlign[da].numberOfTaxa()-1];
                    for (int i = 0; i < donorIndices.length; i++) {donorIndices[i]=i; if(i>=taxon) donorIndices[i]++;}
                } else {
                    donorIndices=new int[donorAlign[da].numberOfTaxa()];
                    for (int i = 0; i < donorIndices.length; i++) {donorIndices[i]=i;}
                }

                //Finds the best haplotype donors for each focus block within a donorGenotypeTable
                DonorHypoth[][] regionHypthInbred=new DonorHypoth[blocks][maxDonorHypotheses.value()];
                byte[][][] targetToDonorDistances=FILLINImputationUtils.calcAllelePresenceCountsBtwTargetAndDonors(maskedTargetBits,
                        donorAlign[da]);
                for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
                    int[] resultRange=FILLINImputationUtils.getBlockWithMinMinorCount(maskedTargetBits[0].getBits(), maskedTargetBits[1].getBits(), focusBlock, minMinorCnt.value(), minMinorCnt.value()*minMajorRatioToMinorCnt);
                    if(resultRange==null) continue; //no data in the focus Block
                    //search for the best inbred donors for a segment
                    regionHypthInbred[focusBlock]=FILLINImputationUtils.findHomozygousDonorHypoth(taxon, resultRange[0], resultRange[2],
                            focusBlock, donorIndices, targetToDonorDistances, minTestSites.value(), maxDonorHypotheses.value());
                }
                impTaxon.setSegmentSolved(false);

                //tries to solve the entire donorAlign region by Virterbi or Inbred
                impTaxon=solveEntireDonorRegion(taxon, donorAlign[da], donorOffset, regionHypthInbred, impTaxon, maskedTargetBits, maxHybridErrorRate.value(), targetToDonorDistances);
                if(impTaxon.isSegmentSolved()) {countFullLength++; continue;}

                //resorts to solving block by block, first by inbred, then by viterbi, and then by hybrid
                impTaxon=solveByBlockNearestNeighbor(impTaxon, taxon, donorAlign[da], donorOffset, regionHypthInbred, hybridNN.value(), maskedTargetBits, minMinorCnt.value(), focusInbredErr, focusHybridErr, focusSmashErr, donorIndices, trackBlockNN, hetsMiss);
                if(impTaxon.isSegmentSolved()) {countByFocus++;}
            }
            double totalFocus= (double)trackBlockNN[3]+(double)trackBlockNN[4];
            sb.append(String.format("InbredOrViterbi:%d FocusBlock:%d PropFocusInbred:%f PropFocusViterbi:%f PropFocusSmash:%f PropFocusMissing:%f BlocksSolved:%d ",
                    countFullLength, countByFocus, (double)trackBlockNN[0]/totalFocus, (double)trackBlockNN[1]/totalFocus, (double)trackBlockNN[2]/totalFocus, (double)trackBlockNN[3]/totalFocus, impTaxon.getBlocksSolved()));
            //sb.append(String.format("InbredOrViterbi:%d FocusBlock:%d BlocksSolved:%d ", countFullLength, countByFocus, impTaxon.getBlocksSolved()));
            int[] unk=FILLINImputationUtils.countUnknownAndHeterozygotes(impTaxon.resolveGeno);
            sb.append(String.format("Unk:%d PropMissing:%g ", unk[0], (double) unk[0] / (double) impTaxon.getOrigGeno().length));
            sb.append(String.format("Het:%d PropHet:%g ", unk[1], (double)unk[1]/(double)impTaxon.getOrigGeno().length));
            sb.append(" BreakPoints:"+impTaxon.getBreakPoints().size());
            if(!isOutputProjection.value()) {
                alignBuilder.addTaxon(unimpAlign.taxa().get(taxon), impTaxon.resolveGeno);
            } else {
                projBuilder.addTaxon(unimpAlign.taxa().get(taxon),impTaxon.getBreakPoints());
            }
            if(verboseOutput) System.out.println(sb.toString());
        }
    }



    /**
     * Solve the entire donor alignment window with either one or two haplotypes.  Generally, it is solved with the
     * HMM Viterbi algorithm.  The Viterbi algorithm is relatively slow the big choice is how to choose the potential
     * donors.
     *
     *
     * @param taxon
     * @param donorAlign
     * @param donorOffset
     * @param regionHypoth
     * @param impT
     * @param maskedTargetBits
     * @param maxHybridErrorRate
     * @return
     */
    private ImputedTaxon solveEntireDonorRegion(int taxon, GenotypeTable donorAlign, int donorOffset,
                DonorHypoth[][] regionHypoth, ImputedTaxon impT, BitSet[] maskedTargetBits, double maxHybridErrorRate, byte[][][] targetToDonorDistances) {

        int blocks=maskedTargetBits[0].getNumWords();
        if(testing==1) System.out.println("Starting complete hybrid search");
        //create a list of the best donors based on showing up frequently high in many focus blocks
        //in the test data set the best results achieved with on the best hypothesis recovered.
        //todo consider whether to use this approach
//        int[] d=FILLINImputationUtils.mostFrequentDonorsAcrossFocusBlocks(regionHypoth, maxDonorHypotheses);
        //Alternative test is find best donors
       int[] d=FILLINImputationUtils.bestDonorsAcrossEntireRegion(targetToDonorDistances, minTestSites.value(),maxDonorHypotheses.value());
        int[] testList=FILLINImputationUtils.fillInc(0,donorAlign.numberOfTaxa()-1);
        int[] bestDonorList=Arrays.copyOfRange(d,0,Math.min(d.length,5));
        DonorHypoth[] bestDBasedOnBest=FILLINImputationUtils.findHeterozygousDonorHypoth(taxon, maskedTargetBits[0].getBits(),
                maskedTargetBits[1].getBits(), 0, blocks-1, blocks/2, donorAlign, bestDonorList, testList, maxDonorHypotheses.value(), minTestSites.value());

        //make all combinations of best donor and find the the pairs that minimize errors
        //with the true switch also will make inbreds
        DonorHypoth[] best2Dsearchdonors=FILLINImputationUtils.findHeterozygousDonorHypoth(taxon, maskedTargetBits[0].getBits(),
                maskedTargetBits[1].getBits(), 0, blocks-1, blocks/2, donorAlign, d, d, maxDonorHypotheses.value(), minTestSites.value());
        DonorHypoth[] best2donors=FILLINImputationUtils.combineDonorHypothArrays(maxDonorHypotheses.value(),bestDBasedOnBest,best2Dsearchdonors);
        if(testing==1) System.out.println(Arrays.toString(best2donors));
        ArrayList<DonorHypoth> goodDH=new ArrayList<DonorHypoth>();
        for (DonorHypoth dh : best2donors) {
            if(dh==null) continue;
            if(dh.isInbred() && (dh.getErrorRate()<maximumInbredError.value())) {
                goodDH.add(dh);
            } else if(dh.getErrorRate()<maxHybridErrorRate) {
                dh=getStateBasedOnViterbi(dh, donorOffset, donorAlign, twoWayViterbi, transition);
                if(dh!=null) goodDH.add(dh);
            }
        }
        if(goodDH.isEmpty()) {
            return impT;
        }
        DonorHypoth[] vdh=new DonorHypoth[goodDH.size()];
        for (int i = 0; i < vdh.length; i++) {vdh[i]=goodDH.get(i);}
        impT.setSegmentSolved(true);
        return setAlignmentWithDonors(donorAlign,vdh, donorOffset,false,impT, false, false);
    }






    /**
     * If the entire donor target regions has too many Mendelian errors that it looks for overlapping regional
     * solutions that are better.
     * @param impT
     * @param targetTaxon
     * @param regionHypth
     */
    private ImputedTaxon solveByBlockNearestNeighbor(ImputedTaxon impT, int targetTaxon, GenotypeTable donorAlign,
               int donorOffset, DonorHypoth[][] regionHypth, boolean hybridMode, BitSet[] maskedTargetBits, int minMinorCnt, double focusInbredErr, double focusHybridErr, double focusSmashErr, int[] donorIndices, int[] blockNN, boolean hetsToMiss) {
        int[] currBlocksSolved= new int[5];//track number of focus blocks solved in NN search for system out; index 0 is inbred, 1 is viterbi, 2 is smash, 3 is not solved, 4 is total for all modes
        int blocks=maskedTargetBits[0].getNumWords();
        for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
            boolean vit= false;//just if viterbi works for this block
            int[] d= getUniqueDonorsForBlock(regionHypth,focusBlock);//get the donor indices for that block from regionHypoth
            int[] resultRange=FILLINImputationUtils.getBlockWithMinMinorCount(maskedTargetBits[0].getBits(), maskedTargetBits[1].getBits(), focusBlock, minMinorCnt, minMinorCnt*minMajorRatioToMinorCnt);
            if(resultRange==null) {currBlocksSolved[3]++; continue; } //no data in the focus Block
            if(d==null) {currBlocksSolved[3]++; continue; } //no donors for the focus block
            //search for the best hybrid donors for a segment
//            DonorHypoth[] best2donors=getBestHybridDonors(targetTaxon, maskedTargetBits[0].getBits(resultRange[0], resultRange[2]),
//                        maskedTargetBits[1].getBits(resultRange[0], resultRange[2]), resultRange[0], resultRange[2], focusBlock, donorAlign, d, d, true);
            DonorHypoth[] best2donors=FILLINImputationUtils.findHeterozygousDonorHypoth(targetTaxon, maskedTargetBits[0].getBits(resultRange[0], resultRange[2]),
                    maskedTargetBits[1].getBits(resultRange[0], resultRange[2]), resultRange[0], resultRange[2], focusBlock, donorAlign, d, d, (int)maxDonorHypotheses.value(), (int)minTestSites.value());


            if(best2donors[0]==null) {currBlocksSolved[3]++; continue; } //no good hybrid donors for the focus block

            //check for Viterbi and inbred
            ArrayList<DonorHypoth> goodDH= new ArrayList<DonorHypoth>(); //KLS0201
            if (best2donors[0].getErrorRate()<focusInbredErr) {
                boolean top= true;
                for (DonorHypoth dh : best2donors) {
                    if(dh==null) continue;
                    if(dh.isInbred() && (dh.getErrorRate()<focusInbredErr)) {
                        goodDH.add(dh);
                    } else if(dh.getErrorRate()<focusHybridErr) {
                        dh=getStateBasedOnViterbi(dh, donorOffset, donorAlign, twoWayViterbi, transition);
                        if(dh!=null) goodDH.add(dh);
                        if (top) vit= true;
                    }
                    top= false;
                }
                if (goodDH.size()!=0) {
                    DonorHypoth[] vdh=new DonorHypoth[goodDH.size()];
                    for (int i = 0; i < vdh.length; i++) {vdh[i]=goodDH.get(i);}
                    regionHypth[focusBlock]= vdh;
                    impT= setAlignmentWithDonors(donorAlign, regionHypth[focusBlock], donorOffset, true, impT, false, false);
                        impT.incBlocksSolved();
                    if(vit==true) {currBlocksSolved[1]++;} else {currBlocksSolved[0]++;}
                    currBlocksSolved[4]++; continue;
                    }
                }
            //if fails, try smash mode
            else if (hybridMode&&(best2donors[0].getErrorRate()<focusSmashErr)) {//smash mode goes into hybrid block NN
                for (DonorHypoth dh:best2donors) {
                    if(dh!=null&&dh.getErrorRate()<focusSmashErr) {
                        goodDH.add(dh);
                    }
                }
                if (goodDH.size()!=0) {
                    DonorHypoth[] vdh=new DonorHypoth[goodDH.size()];
                    for (int i = 0; i < vdh.length; i++) {vdh[i]=goodDH.get(i);}
                    regionHypth[focusBlock]= vdh;
                    impT= setAlignmentWithDonors(donorAlign, regionHypth[focusBlock], donorOffset, true,impT, true, hetsToMiss);//only set donors for focus block //KLS0201
                    impT.incBlocksSolved(); currBlocksSolved[2]++; currBlocksSolved[4]++; continue;
                }
            //if fails, do not impute this focus block    
            } else currBlocksSolved[3]++;
        }
        
        if (currBlocksSolved[4]!=0) impT.setSegmentSolved(true);
        else impT.setSegmentSolved(false);
        int leftNullCnt=currBlocksSolved[3];
        if(testing==1)
            System.out.printf("targetTaxon:%d hybridError:%g block:%d proportionBlocksImputed:%d null:%d inbredDone:%d viterbiDone:%d hybridDone:%d noData:%d %n",
                targetTaxon, focusHybridErr, blocks,currBlocksSolved[4]/blocks, leftNullCnt, currBlocksSolved[0], currBlocksSolved[1], currBlocksSolved[2], currBlocksSolved[3]);
        for (int i = 0; i < currBlocksSolved.length; i++) {blockNN[i]+= currBlocksSolved[i];}
        return impT;
    }

    private DonorHypoth getStateBasedOnViterbi(DonorHypoth dh, int donorOffset, GenotypeTable donorAlign, boolean forwardReverse,
                                               double[][] trans) {
        //Get the coordinates
        int endSite=(dh.endSite>=donorAlign.numberOfSites())?donorAlign.numberOfSites()-1: dh.endSite;
        int sites=endSite-dh.startSite+1;
        //Find the informative sites
        StatePositionChain informative=createInformativeStateChainForViterbi(dh.startSite,
                maxNonMedelian, minimumDonorDistance,
                unimpAlign.genotypeRange(dh.targetTaxon, dh.startSite+donorOffset, endSite+1+donorOffset),
                donorAlign.genotypeRange(dh.donor1Taxon, dh.startSite, endSite+1),
                donorAlign.genotypeRange(dh.donor2Taxon, dh.startSite, endSite+1));
        if(informative==null) return null;
        //Find the most likely states by Virterbi
        int chrlength = donorAlign.chromosomalPosition(endSite) - donorAlign.chromosomalPosition(dh.startSite);
        byte[] callsF=callsFromViterbi(trans, chrlength/sites, informative);
        DonorHypoth dh2=new DonorHypoth(dh.targetTaxon,dh.donor1Taxon, dh.donor2Taxon, dh.startBlock, dh.focusBlock, dh.endBlock);
        dh2.phasedResults= callsF;
        if (forwardReverse) {
            byte[] callsR=callsFromViterbi(trans, chrlength/sites, StatePositionChain.reverseInstance(informative));
            ArrayUtils.reverse(callsR);  //These are now in same direction as callsF
            byte[] callsC=new byte[callsF.length];
            System.arraycopy(callsR,0,callsC,0,callsC.length/2);
            System.arraycopy(callsF,callsC.length/2,callsC,callsC.length/2,callsC.length-callsC.length/2);
            dh2.phasedResults= callsC;
        }
        return dh2;
    }

    private byte[] callsFromViterbi(double[][] trans, int avgChrLength, StatePositionChain informative) {
        //There is no EM optimization of Viterbi, clearly something that could be changed
        TransitionProbability tpF = new TransitionProbability();
        EmissionProbability ep = new EmissionProbability();
        tpF.setTransitionProbability(trans);
        ep.setEmissionProbability(emission);
        final double probHeterozygous=0.5;
        final double phom = (1 - probHeterozygous) / 2;
        final double[] pTrue = new double[]{phom, .25*probHeterozygous ,.5 * probHeterozygous, .25*probHeterozygous, phom};
        tpF.setAverageSegmentLength(avgChrLength);
        tpF.setPositions(informative.informSites);
        ViterbiAlgorithm vaF = new ViterbiAlgorithm(informative.informStates, tpF, ep, pTrue);
        vaF.calculate();
        byte[] resultStatesF=vaF.getMostProbableStateSequence();
        int currPos=0;
        //converts the informative states back to all states
        byte[] callsF=new byte[informative.totalSiteCnt];
        for(int cs=0; cs<informative.totalSiteCnt; cs++) {
            callsF[cs]=(resultStatesF[currPos]==1)?(byte)1:(byte)(resultStatesF[currPos]/2); //converts the scale back to 0,1,2 from 0..4
            if((informative.informSites[currPos]<cs+informative.startSite)&&(currPos<resultStatesF.length-1)) currPos++;
        }
        return callsF;
    }

    /**
     * Creates the state chain for three genotypes, after excluding the uninformative sites
     */
    private StatePositionChain createInformativeStateChainForViterbi(int startSite, double maxNonMendelian,
                                            double minDonorDistance ,byte[] targetGenotype, byte[] donor1Genotype, byte[] donor2Genotype) {
        int nonMendel=0, donorDifferences=0;
        ArrayList<Byte> nonMissingObs = new ArrayList<>();
        ArrayList<Integer> informSites = new ArrayList<>();
        //Selects only the informative sites for the Viterbi algorithm (known & gap free & homozygous in donor)
        for(int cs=0; cs<targetGenotype.length; cs++) {
            if(targetGenotype[cs]==UNKNOWN_DIPLOID_ALLELE) continue;
            if(donor1Genotype[cs]==UNKNOWN_DIPLOID_ALLELE) continue;
            if(donor2Genotype[cs]==UNKNOWN_DIPLOID_ALLELE) continue;
            if(targetGenotype[cs]==GAP_DIPLOID_ALLELE) continue;
            if(donor1Genotype[cs]==GAP_DIPLOID_ALLELE) continue;
            if(donor2Genotype[cs]==GAP_DIPLOID_ALLELE) continue;
            if(isHeterozygous(donor1Genotype[cs]) || isHeterozygous(donor2Genotype[cs])) continue;
            if(donor1Genotype[cs]==donor2Genotype[cs]) {
                if(targetGenotype[cs]!=donor1Genotype[cs]) nonMendel++;
                continue;
            } else {
                donorDifferences++;
            }
            byte state=1;
            if(targetGenotype[cs]==donor1Genotype[cs]) {state=0;}
            else if(targetGenotype[cs]==donor2Genotype[cs]) {state=2;}
            nonMissingObs.add(state);
            informSites.add(cs+startSite);
        }
        if(informSites.size()<10) return null;
        //if the target has too many unexplained sites then return false
        if((double)nonMendel/(double)informSites.size()>maxNonMendelian) return null;
        //if the donors are too similar Viterbi performs poorly.  Only accepted different donors
        if((double)donorDifferences/(double)informSites.size()<minDonorDistance) return null;
        return new StatePositionChain(startSite,targetGenotype.length,Bytes.toArray(nonMissingObs), Ints.toArray(informSites));
    }

    private static class StatePositionChain {
        final int startSite;
        final int totalSiteCnt;
        final byte[] informStates;
        final int[] informSites;

        private StatePositionChain(int startSite, int totalSiteCnt, byte[] informStates, int[] informSites) {
            this.startSite=startSite;
            this.totalSiteCnt=totalSiteCnt;
            this.informStates=informStates;
            this.informSites=informSites;
        }

        static protected StatePositionChain reverseInstance(StatePositionChain forwardSPC) {
            byte[] informStatesR=Arrays.copyOf(forwardSPC.informStates,forwardSPC.informStates.length);
            ArrayUtils.reverse(informStatesR);
            //this retains the distance between sites in the reverse
            int[] informSitesReverse= new int[forwardSPC.informSites.length];
            for (int site = 0; site < informSitesReverse.length; site++) {
                informSitesReverse[site]= forwardSPC.totalSiteCnt-1-forwardSPC.informSites[forwardSPC.informSites.length-site-1];
            }
            return new StatePositionChain(informSitesReverse[0], forwardSPC.totalSiteCnt, informStatesR, informSitesReverse);
        }
    }

    private int[] getUniqueDonorsForBlock(DonorHypoth[][] regionHypth, int block) {//change so only adding those donors that are not identical for target blocks
        Set<Integer> donors= new HashSet<>();
        for (int h = 0; h < regionHypth[block].length; h++) {
            if (regionHypth[block][h]!=null) {
                donors.add(regionHypth[block][h].donor1Taxon);
                donors.add(regionHypth[block][h].donor2Taxon);
            }
        }
        if (donors.isEmpty()) return null;
        return Ints.toArray(donors);
    }

//    /**
//     * Simple algorithm that tests every possible two donor combination to minimize
//     * the number of unmatched informative alleles.  Currently, there is little tie
//     * breaking, longer matches are favored.
//     * @param targetTaxon
//     * @param startBlock
//     * @param endBlock
//     * @param focusBlock
//     * @return int[] array of {donor1, donor2, testSites}  sorted by  testPropUnmatched
//     */
//    private DonorHypoth[] getBestHybridDonors(int targetTaxon, long[] mjT, long[] mnT,
//                                              int startBlock, int endBlock, int focusBlock, GenotypeTable donorAlign,
//                                              int[] donor1Indices, int[] donor2Indices,
//                                              boolean viterbiSearch) {
//        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
//        Set<Long> testedDonorPairs=new HashSet<>();
//        double lastKeytestPropUnmatched=1.0;
//        double inc=1e-9;
//        double[] donorDist;
//        for (int d1: donor1Indices) {
//            long[] mj1=donorAlign.allelePresenceForSitesBlock(d1, Major, startBlock, endBlock + 1);
//            long[] mn1=donorAlign.allelePresenceForSitesBlock(d1, Minor, startBlock, endBlock + 1);
//            for (int d2 : donor2Indices) {
//                if((!viterbiSearch)&&(d1==d2)) continue;
//                if(testedDonorPairs.contains(((long)d1<<32)+(long)d2)) continue;
//                if(testedDonorPairs.contains(((long)d2<<32)+(long)d1)) continue;
//                testedDonorPairs.add(((long)d1<<32)+(long)d2);
//                long[] mj2=donorAlign.allelePresenceForSitesBlock(d2, Major, startBlock, endBlock + 1);
//                long[] mn2=donorAlign.allelePresenceForSitesBlock(d2, Minor, startBlock, endBlock + 1);
//                if(viterbiSearch) {
//                    donorDist=IBSDistanceMatrix.computeHetBitDistances(mj1, mn1, mj2, mn2, minTestSites);
//                    if((d1!=d2)&&(donorDist[0]<this.minimumDonorDistance)) continue;
//                }
//                int[] mendErr=mendelErrorComparison(mjT, mnT, mj1, mn1, mj2, mn2);
//                if(mendErr[1]<minTestSites) continue;
//                double testPropUnmatched=(double)(mendErr[0])/(double)mendErr[1];
//                inc+=1e-9;
//                testPropUnmatched+=inc;
//                if(testPropUnmatched<lastKeytestPropUnmatched) {
//                    DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d2, startBlock,
//                            focusBlock, endBlock, mendErr[1], mendErr[0]);
//                    DonorHypoth prev=bestDonors.put(new Double(testPropUnmatched), theDH);
//                    if(prev!=null) {
//                        System.out.println("Hybrid TreeMap index crash:"+testPropUnmatched);
//                    }
//                    if(bestDonors.size()>maxDonorHypotheses) {
//                        bestDonors.remove(bestDonors.lastKey());
//                        lastKeytestPropUnmatched=bestDonors.lastKey();
//                    }
//                }
//            }
//        }
//        DonorHypoth[] result=new DonorHypoth[maxDonorHypotheses];
//        int count=0;
//        for (DonorHypoth dh : bestDonors.values()) {
//            result[count]=dh;
//            count++;
//        }
//        return result;
//    }
//
//    private int[] mendelErrorComparison(long[] mjT, long[] mnT, long[] mj1, long[] mn1,
//                                        long[] mj2, long[] mn2) {
//        int mjUnmatched=0;
//        int mnUnmatched=0;
//        int testSites=0;
//        for (int i = 0; i < mjT.length; i++) {
//            long siteMask=(mjT[i]|mnT[i])&(mj1[i]|mn1[i])&(mj2[i]|mn2[i]);
//            mjUnmatched+=Long.bitCount(siteMask&mjT[i]&(mjT[i]^mj1[i])&(mjT[i]^mj2[i]));
//            mnUnmatched+=Long.bitCount(siteMask&mnT[i]&(mnT[i]^mn1[i])&(mnT[i]^mn2[i]));
//            testSites+=Long.bitCount(siteMask);
//        }
//        int totalMendelianErrors=mjUnmatched+mnUnmatched;
//        // double testPropUnmatched=(double)(totalMendelianErrors)/(double)testSites;
//        return (new int[] {totalMendelianErrors, testSites});
//    }


    private ImputedTaxon setAlignmentWithDonors(GenotypeTable donorAlign, DonorHypoth[] theDH, int donorOffset,
                                                boolean setJustFocus, ImputedTaxon impT, boolean smashOn, boolean hetsMiss) {
        if(theDH[0].targetTaxon<0) return impT;
        boolean print=false;
        int startSite=(setJustFocus)?theDH[0].getFocusStartSite():theDH[0].startSite;
        int endSite=(setJustFocus)?theDH[0].getFocusEndSite():theDH[0].endSite;
        if(endSite>=donorAlign.numberOfSites()) endSite=donorAlign.numberOfSites()-1;
        //prevDonors is used to track recombination breakpoint start and stop
        int[] prevDonors=new int[]{-1, -1};
        if(theDH[0].getPhaseForSite(startSite)==0) {prevDonors=new int[]{theDH[0].donor1Taxon, theDH[0].donor1Taxon};}
        else if(theDH[0].getPhaseForSite(startSite)==2) {prevDonors=new int[]{theDH[0].donor2Taxon, theDH[0].donor2Taxon};}
        else if(theDH[0].getPhaseForSite(startSite)==1) {prevDonors=new int[]{theDH[0].donor1Taxon, theDH[0].donor2Taxon};}
        int prevDonorStart=startSite;

        int[] currDonors=prevDonors;

        for(int cs=startSite; cs<=endSite; cs++) {
            byte donorEst=UNKNOWN_DIPLOID_ALLELE;
            byte neighbor=0;
            //site imputation will try to use all donor hypotheses; breakpoints only use the best (e.g. i==0 tests)
            for (int i = 0; (i < theDH.length) && (donorEst==UNKNOWN_DIPLOID_ALLELE); i++) {
                neighbor++;
                if((theDH[i]==null)||(theDH[i].donor1Taxon<0)) continue;
                if(theDH[i].getErrorRate()>this.maximumInbredError.value()) continue;
                byte bD1=donorAlign.genotype(theDH[i].donor1Taxon, cs);
                if(theDH[i].getPhaseForSite(cs)==0) {
                    donorEst=bD1;
                    if(i==0) currDonors=new int[]{theDH[0].donor1Taxon, theDH[0].donor1Taxon};
                }
                else {
                    byte bD2=donorAlign.genotype(theDH[i].donor2Taxon, cs);
                    if(theDH[i].getPhaseForSite(cs)==2) {
                        donorEst=bD2;
                        if(i==0) currDonors=new int[]{theDH[0].donor2Taxon, theDH[0].donor2Taxon};
                    }
                    else {
                        donorEst=GenotypeTableUtils.getUnphasedDiploidValueNoHets(bD1, bD2);
                        if(i==0) currDonors=new int[]{theDH[0].donor1Taxon, theDH[0].donor2Taxon};
                    }
                }
            }
            byte knownBase=impT.getOrigGeno(cs+donorOffset);
            //record recombination breakpoint if it changes
            if(!Arrays.equals(prevDonors, currDonors)) {
                DonorHaplotypes dhaps=new DonorHaplotypes(donorAlign.chromosome(prevDonorStart), donorAlign.chromosomalPosition(prevDonorStart),
                        donorAlign.chromosomalPosition(cs),prevDonors[0],prevDonors[1]);
                impT.addBreakPoint(dhaps);
                prevDonors=currDonors;
                prevDonorStart=cs;
            }
            //chgHis records the section of the imputation pipeline that made the change, Viterbi negative, blockNN positive
            if(theDH[0].phasedResults==null) {impT.chgHis[cs+donorOffset]=(byte)-neighbor;}
            else {impT.chgHis[cs+donorOffset]=(byte)neighbor;}

            impT.impGeno[cs+donorOffset]= donorEst;  //predicted based on neighbor
            //if genotype is unknown or het undercall then resolves
            //todo there is currently not an option if the predicted disagrees with a homozygous call
            if(knownBase==UNKNOWN_DIPLOID_ALLELE) {
                if (isHeterozygous(donorEst)) {
                    if (smashOn && hetsMiss) {//if imputing a heterozygote, just set to missing
                        impT.resolveGeno[cs+donorOffset]= knownBase;
                    }
                    else impT.resolveGeno[cs+donorOffset]= donorEst; //if not in hybrid, set to het
                }
                else {//if not imputed to a het
                    impT.resolveGeno[cs+donorOffset]= donorEst;}
            } else if (isHeterozygous(donorEst)){
                if(resolveHetIfUndercalled&&GenotypeTableUtils.isPartiallyEqual(knownBase, donorEst)&&smashOn==false){//if smash off, set homozygotes imputed to het to het
                    //System.out.println("ResolveHet:"+theDH[0].targetTaxon+":"+cs+donorOffset+":"+knownBaseString+":"+donorEstString);
                    impT.resolveGeno[cs+donorOffset]= donorEst;
                }
            }
        } //end of cs loop
        //enter a stop of the DH at the beginning of the next block
        DonorHaplotypes dhaps=new DonorHaplotypes(donorAlign.chromosome(prevDonorStart), donorAlign.chromosomalPosition(prevDonorStart),
                donorAlign.chromosomalPosition(endSite),prevDonors[0],prevDonors[1]);
        impT.addBreakPoint(dhaps);
        return impT;

    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Impute By FILLIN";
    }

    @Override
    public String getToolTipText() {
        return "Imputation that relies on a combination of HMM and Nearest Neighbor using donors derived from FILLINFindHaplotypesPlugin";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
     public static void main(String[] args) {
         GeneratePluginCode.generate(FILLINImputationPlugin.class);
     }

    /**
     * Convenience method to run plugin with one return object.
     */
    /*// TODO: Replace <Type> with specific type.
    public <Type> runPlugin(DataSet input) {
        return (<Type>) performFunction(input).getData(0).getData();
    }*/

    /**
     * Input HapMap file of target genotypes to impute. Accepts
     * all file types supported by TASSEL5
     *
     * @return Target file
     */
    public String targetFile() {
        return hmpFile.value();
    }

    /**
     * Set Target file. Input HapMap file of target genotypes
     * to impute. Accepts all file types supported by TASSEL5
     *
     * @param value Target file
     *
     * @return this plugin
     */
    public FILLINImputationPlugin targetFile(String value) {
        hmpFile = new PluginParameter<>(hmpFile, value);
        return this;
    }

    /**
     * Directory containing donor haplotype files from output
     * of FILLINFindHaplotypesPlugin. All files with '.gc'
     * in the filename will be read in, only those with matching
     * sites are used
     *
     * @return Donor Dir
     */
    public String donorDir() {
        return donorFile.value();
    }

    /**
     * Set Donor Dir. Directory containing donor haplotype
     * files from output of FILLINFindHaplotypesPlugin. All
     * files with '.gc' in the filename will be read in, only
     * those with matching sites are used
     *
     * @param value Donor Dir
     *
     * @return this plugin
     */
    public FILLINImputationPlugin donorDir(String value) {
        donorFile = new PluginParameter<>(donorFile, value);
        return this;
    }

    /**
     * Output file; hmp.txt.gz and .hmp.h5 accepted.
     *
     * @return Output filename
     */
    public String outputFilename() {
        return outFileBase.value();
    }

    /**
     * Set Output filename. Output file; hmp.txt.gz and .hmp.h5
     * accepted.
     *
     * @param value Output filename
     *
     * @return this plugin
     */
    public FILLINImputationPlugin outputFilename(String value) {
        outFileBase = new PluginParameter<>(outFileBase, value);
        return this;
    }

    /**
     * Preferred haplotype block size in sites (use same as
     * in FILLINFindHaplotypesPlugin)
     *
     * @return Preferred haplotype size
     */
    public Integer preferredHaplotypeSize() {
        return appoxSitesPerDonorGenotypeTable.value();
    }

    /**
     * Set Preferred haplotype size. Preferred haplotype block
     * size in sites (use same as in FILLINFindHaplotypesPlugin)
     *
     * @param value Preferred haplotype size
     *
     * @return this plugin
     */
    public FILLINImputationPlugin preferredHaplotypeSize(Integer value) {
        appoxSitesPerDonorGenotypeTable = new PluginParameter<>(appoxSitesPerDonorGenotypeTable, value);
        return this;
    }

    /**
     * Threshold per taxon heterozygosity for treating taxon
     * as heterozygous (no Viterbi, het thresholds).
     *
     * @return Heterozygosity threshold
     */
    public Double heterozygosityThreshold() {
        return hetThresh.value();
    }

    /**
     * Set Heterozygosity threshold. Threshold per taxon heterozygosity
     * for treating taxon as heterozygous (no Viterbi, het
     * thresholds).
     *
     * @param value Heterozygosity threshold
     *
     * @return this plugin
     */
    public FILLINImputationPlugin heterozygosityThreshold(Double value) {
        hetThresh = new PluginParameter<>(hetThresh, value);
        return this;
    }

    /**
     * Maximum error rate for applying one haplotype to entire
     * site window
     *
     * @return Max error to impute one donor
     */
    public Double maxErrorToImputeOneDonor() {
        return maximumInbredError.value();
    }

    /**
     * Set Max error to impute one donor. Maximum error rate
     * for applying one haplotype to entire site window
     *
     * @param value Max error to impute one donor
     *
     * @return this plugin
     */
    public FILLINImputationPlugin maxErrorToImputeOneDonor(Double value) {
        maximumInbredError = new PluginParameter<>(maximumInbredError, value);
        return this;
    }

    /**
     * Maximum error rate for applying Viterbi with to haplotypes
     * to entire site window
     *
     * @return Max combined error to impute two donors
     */
    public Double maxCombinedErrorToImputeTwoDonors() {
        return maxHybridErrorRate.value();
    }

    /**
     * Set Max combined error to impute two donors. Maximum
     * error rate for applying Viterbi with to haplotypes
     * to entire site window
     *
     * @param value Max combined error to impute two donors
     *
     * @return this plugin
     */
    public FILLINImputationPlugin maxCombinedErrorToImputeTwoDonors(Double value) {
        maxHybridErrorRate = new PluginParameter<>(maxHybridErrorRate, value);
        return this;
    }

    /**
     * Minimum number of sites to test for IBS between haplotype
     * and target in focus block
     *
     * @return Min sites to test match
     */
    public Integer minSitesToTestMatch() {
        return minTestSites.value();
    }

    /**
     * Set Min sites to test match. Minimum number of sites
     * to test for IBS between haplotype and target in focus
     * block
     *
     * @param value Min sites to test match
     *
     * @return this plugin
     */
    public FILLINImputationPlugin minSitesToTestMatch(Integer value) {
        minTestSites = new PluginParameter<>(minTestSites, value);
        return this;
    }

    /**
     * Minimum number of informative minor alleles in the
     * search window (or 10X major)
     *
     * @return Min num of minor alleles to compare
     */
    public Integer minNumOfMinorAllelesToCompare() {
        return minMinorCnt.value();
    }

    /**
     * Set Min num of minor alleles to compare. Minimum number
     * of informative minor alleles in the search window (or
     * 10X major)
     *
     * @param value Min num of minor alleles to compare
     *
     * @return this plugin
     */
    public FILLINImputationPlugin minNumOfMinorAllelesToCompare(Integer value) {
        minMinorCnt = new PluginParameter<>(minMinorCnt, value);
        return this;
    }

    /**
     * Maximum number of donor hypotheses to be explored
     *
     * @return Max donor hypotheses
     */
    public Integer maxDonorHypotheses() {
        return maxDonorHypotheses.value();
    }

    /**
     * Set Max donor hypotheses. Maximum number of donor hypotheses
     * to be explored
     *
     * @param value Max donor hypotheses
     *
     * @return this plugin
     */
    public FILLINImputationPlugin maxDonorHypotheses(Integer value) {
        maxDonorHypotheses = new PluginParameter<>(maxDonorHypotheses, value);
        return this;
    }

    /**
     * Write all imputed heterozygous calls as such, even
     * if the original file has a homozygous call. (Not recommended
     * for inbred lines.)
     *
     * @return Impute all het calls
     */
    public Boolean imputeAllHetCalls() {
        return imputeAllHets.value();
    }

    /**
     * Set Impute all het calls. Write all imputed heterozygous
     * calls as such, even if the original file has a homozygous
     * call. (Not recommended for inbred lines.)
     *
     * @param value Impute all het calls
     *
     * @return this plugin
     */
    public FILLINImputationPlugin imputeAllHetCalls(Boolean value) {
        imputeAllHets = new PluginParameter<>(imputeAllHets, value);
        return this;
    }

    /**
     * If true, uses combination mode in focus block, else
     * does not impute
     *
     * @return Combine two haplotypes as heterozygote
     */
    public Boolean combineTwoHaplotypesAsHeterozygote() {
        return hybridNN.value();
    }

    /**
     * Set Combine two haplotypes as heterozygote. If true,
     * uses combination mode in focus block, else does not
     * impute
     *
     * @param value Combine two haplotypes as heterozygote
     *
     * @return this plugin
     */
    public FILLINImputationPlugin combineTwoHaplotypesAsHeterozygote(Boolean value) {
        hybridNN = new PluginParameter<>(hybridNN, value);
        return this;
    }

    /**
     * Create a projection alignment for high density markers
     *
     * @return Output projection alignment
     */
    public Boolean outputProjectionAlignment() {
        return isOutputProjection.value();
    }

    /**
     * Set Output projection alignment. Create a projection
     * alignment for high density markers
     *
     * @param value Output projection alignment
     *
     * @return this plugin
     */
    public FILLINImputationPlugin outputProjectionAlignment(Boolean value) {
        isOutputProjection = new PluginParameter<>(isOutputProjection, value);
        return this;
    }

    /**
     * Impute the donor file itself
     *
     * @return Impute donor file
     */
    public Boolean imputeDonorFile() {
        return imputeDonorFile.value();
    }

    /**
     * Set Impute donor file. Impute the donor file itself
     *
     * @param value Impute donor file
     *
     * @return this plugin
     */
    public FILLINImputationPlugin imputeDonorFile(Boolean value) {
        imputeDonorFile = new PluginParameter<>(imputeDonorFile, value);
        return this;
    }

    /**
     * Supress system out
     *
     * @return Supress system out
     */
    public Boolean supressSystemOut() {
        return nonverboseOutput.value();
    }

    /**
     * Set Supress system out. Supress system out
     *
     * @param value Supress system out
     *
     * @return this plugin
     */
    public FILLINImputationPlugin supressSystemOut(Boolean value) {
        nonverboseOutput = new PluginParameter<>(nonverboseOutput, value);
        return this;
    }

    /**
     * Masks input file before imputation and calculates accuracy
     * based on masked genotypes
     *
     * @return Calculate accuracy
     */
    public Boolean calculateAccuracy() {
        return accuracy.value();
    }

    /**
     * Set Calculate accuracy. Masks input file before imputation
     * and calculates accuracy based on masked genotypes
     *
     * @param value Calculate accuracy
     *
     * @return this plugin
     */
    public FILLINImputationPlugin calculateAccuracy(Boolean value) {
        accuracy = new PluginParameter<>(accuracy, value);
        return this;
    }

    /**
     * Proportion of genotypes to mask for accuracy calculation
     * if depth not available
     *
     * @return Proportion of genotypes to mask if no depth
     */
    public Double proportionOfGenotypesToMaskIfNoDepth() {
        return propSitesMask.value();
    }

    /**
     * Set Proportion of genotypes to mask if no depth. Proportion
     * of genotypes to mask for accuracy calculation if depth
     * not available
     *
     * @param value Proportion of genotypes to mask if no depth
     *
     * @return this plugin
     */
    public FILLINImputationPlugin proportionOfGenotypesToMaskIfNoDepth(Double value) {
        propSitesMask = new PluginParameter<>(propSitesMask, value);
        return this;
    }

    /**
     * Depth of genotypes to mask for accuracy calculation
     * if depth information available
     *
     * @return Depth of genotypes to mask
     */
    public Integer depthOfGenotypesToMask() {
        return depthToMask.value();
    }

    /**
     * Set Depth of genotypes to mask. Depth of genotypes
     * to mask for accuracy calculation if depth information
     * available
     *
     * @param value Depth of genotypes to mask
     *
     * @return this plugin
     */
    public FILLINImputationPlugin depthOfGenotypesToMask(Integer value) {
        depthToMask = new PluginParameter<>(depthToMask, value);
        return this;
    }

    /**
     * Proportion of genotypes of given depth to mask for
     * accuracy calculation if depth available
     *
     * @return Proportion of depth genotypes to mask
     */
    public Double proportionOfDepthGenotypesToMask() {
        return propDepthSitesMask.value();
    }

    /**
     * Set Proportion of depth genotypes to mask. Proportion
     * of genotypes of given depth to mask for accuracy calculation
     * if depth available
     *
     * @param value Proportion of depth genotypes to mask
     *
     * @return this plugin
     */
    public FILLINImputationPlugin proportionOfDepthGenotypesToMask(Double value) {
        propDepthSitesMask = new PluginParameter<>(propDepthSitesMask, value);
        return this;
    }

    /**
     * Key to calculate accuracy. Genotypes missing (masked)
     * in target file should be present in key, with all other
     * sites set to missing. Overrides otheraccuracy options
     * if present and all sites and taxa present in target
     * file present
     *
     * @return Optional key to calculate accuracy
     */
    public String optionalKeyToCalculateAccuracy() {
        return maskKey.value();
    }

    /**
     * Set Optional key to calculate accuracy. Key to calculate
     * accuracy. Genotypes missing (masked) in target file
     * should be present in key, with all other sites set
     * to missing. Overrides otheraccuracy options if present
     * and all sites and taxa present in target file present
     *
     * @param value Optional key to calculate accuracy
     *
     * @return this plugin
     */
    public FILLINImputationPlugin optionalKeyToCalculateAccuracy(String value) {
        maskKey = new PluginParameter<>(maskKey, value);
        return this;
    }

    /**
     * Calculate R2 accuracy within MAF categories based on
     * donor file
     *
     * @return Calculate accuracy within MAF categories
     */
    public Boolean calculateAccuracyWithinMAFCategories() {
        return byMAF.value();
    }

    /**
     * Set Calculate accuracy within MAF categories. Calculate
     * R2 accuracy within MAF categories based on donor file
     *
     * @param value Calculate accuracy within MAF categories
     *
     * @return this plugin
     */
    public FILLINImputationPlugin calculateAccuracyWithinMAFCategories(Boolean value) {
        byMAF = new PluginParameter<>(byMAF, value);
        return this;
    }
}
