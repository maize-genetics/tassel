/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.analysis.gbs.pana;


import java.io.File;
import java.util.Arrays;
import net.maizegenetics.analysis.gbs.MergeMultipleTagCountPlugin;

/**
 *
 * @author Fei Lu
 */
public class PanAUsageExample {
    
    public PanAUsageExample () {
        this.h5ToAnchorPlugin();
        //this.splitTBTPlugin();
        //this.buildTagBlockPositionPlugin();
        //this.splitTagBlockPositionPlugin();
        //this.GWASMappingPlugin();
        //this.mergeMappingResultPlugin();
        //this.mappingResultToTagGWASMapPlugin();
        //this.tagMapToFastaPlugin();
        //this.alignmentWithBowtie2();
        //this.samToMultiPositionTOPMPlugin();
        //this.addPosToTagMapPlugin();
        //this.buildTrainingSetPlugin();
        //this.modelTrainingPlugin();
        //this.predictionPlugin();
        //this.filterTagMapPlugin();
        
        
        this.readDigestPlugin();
        //this.mergeMultipleTagCountPlugin();
        //this.buildPivotTBTPlugin();
    }
    
    public void readDigestPlugin () {
        String rawSeqDirS = "M:\\pipelineTest\\PanA\\Illumina\\fastq\\";
        String keyFileS = "M:\\pipelineTest\\PanA\\key\\keyFastq.txt";
        String recSeq = "GCTG";
        int customTagLength = 96;
        String outputDirS = "M:\\pipelineTest\\PanA\\tagCount\\";
        String arguments = "-i " + rawSeqDirS + " -k " + keyFileS + " -s " + recSeq + " -l " + String.valueOf(customTagLength)+ " -f 0 -o " + outputDirS;
        String[] args = arguments.split(" ");
        PanAReadDigestPlugin rdp = new PanAReadDigestPlugin();
        rdp.setParameters(args);
        rdp.performFunction(null);
    }
    
    public void mergeMultipleTagCountPlugin () {
        String tagCountDirectory = "M:\\pipelineTest\\PanA\\tagCount\\";
        String masterTagCount = "M:\\pipelineTest\\PanA\\masterTagCount\\master.cnt";
        int minCount = 1;
        String arguments = "-i " + tagCountDirectory + " -o " + masterTagCount + " -c " + String.valueOf(minCount);
        String[] args = arguments.split(" ");
        MergeMultipleTagCountPlugin p = new MergeMultipleTagCountPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void buildPivotTBTPlugin () {
        String masterTagCountFileS = "M:\\pipelineTest\\PanA\\masterTagCount\\master.cnt";
        String tagCountDirS = "M:\\pipelineTest\\PanA\\tagCount";
        String tbtFileS = "M:\\pipelineTest\\PanA\\tbt\\TBT_shotgun_pivot.h5";
        String arguments = "-m " + masterTagCountFileS + " -d " + tagCountDirS + " -o " + tbtFileS;
        String[] args = arguments.split(" ");
        PanABuildPivotTBTPlugin p = new PanABuildPivotTBTPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void h5ToAnchorPlugin () {
        String h5GentoypeFileS = "M:\\pipelineTest\\PanA\\genotype\\GBS27.1024sites.T5.imp.hmp.h5";
        String sBitGenotypeFileS = "M:\\pipelineTest\\PanA\\genotype\\GBS27.1024sites.sBit.h5";
        String arguments = "-i " + h5GentoypeFileS +  " -o " + sBitGenotypeFileS;
        String[] args = arguments.split(" ");
        PanAH5ToAnchorPlugin hta = new PanAH5ToAnchorPlugin();
        hta.setParameters(args);
        hta.performFunction(null);
    }
    
    public void splitTBTPlugin () {
        String inputTBTS = "M:\\pipelineTest\\PanA\\tbt\\TBTHDF5_4096_tags_mergedtaxa_pivot_20120921.h5";
        String outputDirS = "M:\\pipelineTest\\PanA\\tbt\\subTBT\\";
        String chunkSize = "1000";
        String arguments = "-i " + inputTBTS + " -s " + chunkSize + " -o " + outputDirS;
        String[] args = arguments.split(" ");
        PanASplitTBTPlugin sbp = new PanASplitTBTPlugin();
        sbp.setParameters(args);
        sbp.performFunction(null);
    }
    
    public void buildTagBlockPositionPlugin () {
        String tbtHDF5 = "M:\\pipelineTest\\PanA\\tbt\\TBTHDF5_4096_tags_mergedtaxa_pivot_20120921.h5";
        String topmFileS = "M:\\production\\geneticMapping\\tagBlock\\AllZeaGBSv2.6ProdTOPM_20130605.topm.h5";
        String blockFileS = "M:\\pipelineTest\\PanA\\tbp\\tagBlock.tbp";
        int TOPMVersionValue = 1;
        String arguments = "-t " + tbtHDF5 + " -p " + topmFileS + " -v " + String.valueOf(TOPMVersionValue) + " -o " + blockFileS;
        String[] args = arguments.split(" ");
        PanABuildTagBlockPosPlugin tbp = new PanABuildTagBlockPosPlugin();
        tbp.setParameters(args);
        tbp.performFunction(null);
    }
    
    public void splitTagBlockPositionPlugin () {
        String blockFileS = "M:\\pipelineTest\\PanA\\tbp\\tagBlock.tbp";
        String outputDirS = "M:\\pipelineTest\\PanA\\tbp\\subTBP";
        String chunkSize = "1000";
        String arguments = "-i " + blockFileS + " -s " + chunkSize + " -o " + outputDirS;
        String[] args = arguments.split(" ");
        PanASplitTagBlockPosPlugin stb = new PanASplitTagBlockPosPlugin();
        stb.setParameters(args);
        stb.performFunction(null);
    }
    
    public void GWASMappingPlugin () {
        String sBitGenotypeFileS = "M:\\pipelineTest\\PanA\\genotype\\GBS27.1024sites.sBit.h5";
        String outDirS = "M:\\pipelineTest\\PanA\\gwasResult\\sub\\";
        String tbtDirS = "M:\\pipelineTest\\PanA\\tbt\\subTBT\\";
        String tbpDirS = "M:\\pipelineTest\\PanA\\tbp\\subTBP\\";
        File[] tbts = new File (tbtDirS).listFiles();
        File[] tbps = new File (tbpDirS).listFiles();
        Arrays.sort(tbts);
        Arrays.sort(tbps);
        for (int i = 0; i < tbts.length; i++) {
            String arguments = "-g " + sBitGenotypeFileS + " -t " + tbts[i].getAbsolutePath() + " -b " + tbps[i].getAbsolutePath() + " -o " + outDirS + " -c max -s 1000 -cs 0 -ce 1";
            
            String[] args = arguments.split(" ");
            PanATagGWASMappingPlugin tgm = new  PanATagGWASMappingPlugin();
            tgm.setParameters(args);
            tgm.performFunction(null);
        }
    }
    
    public void mergeMappingResultPlugin () {
        String subResultDirS = "M:\\pipelineTest\\PanA\\gwasResult\\sub\\";
        String mergedResultFileS = "M:\\pipelineTest\\PanA\\gwasResult\\pivotTBT.gwas.txt";
        String arguments = "-i " + subResultDirS + " -o " + mergedResultFileS;
        String[] args = arguments.split(" ");
        PanAMergeMappingResultPlugin mmr = new PanAMergeMappingResultPlugin();
        mmr.setParameters(args);
        mmr.performFunction(null);
    }
    
    
    public void mappingResultToTagGWASMapPlugin () {
        String mappingResultFileS = "M:\\pipelineTest\\PanA\\gwasResult\\pivotTBT.gwas.txt";
        String tagCountFileS = "M:/production/v3gbs/tagCount/AllZeaMasterTags_c10_20120606.cnt";
        String tagGWASMapFileS = "M:\\pipelineTest\\PanA\\tagMap\\tagGWASMap.h5";
        String arguments = "-i " + mappingResultFileS + " -t " + tagCountFileS + " -o " + tagGWASMapFileS;
        String[] args = arguments.split(" ");
        PanABuildTagGWASMapPlugin mrtg = new PanABuildTagGWASMapPlugin();
        mrtg.setParameters(args);
        mrtg.performFunction(null);
    }
    
    public void tagMapToFastaPlugin () {
        String tagGWASMapFileS = "M:\\pipelineTest\\PanA\\tagMap\\tagGWASMap.h5";
        String fastaFileS = "M:\\pipelineTest\\PanA\\alignment\\tagGWASMap.fa";
        String arguments = "-i " + tagGWASMapFileS +  " -o " + fastaFileS;
        String[] args = arguments.split(" ");
        PanATagMapToFastaPlugin tmtf = new PanATagMapToFastaPlugin();
        tmtf.setParameters(args);
        tmtf.performFunction(null);
    }
    
    public void alignmentWithBowtie2 () {
        String command = "bowtie2 -x ZmB73_RefGen_v2.fa -f tagGWASMap.fa -k 2 --very-sensitive-local -p 8 -S tagGWASMap.sam";
    }
    
    public void samToMultiPositionTOPMPlugin () {
        String samFileS = "M:\\pipelineTest\\PanA\\alignment\\tagGWASMap.sam";
        String tagGWASMapFileS = "M:\\pipelineTest\\PanA\\tagMap\\tagGWASMap.h5";
        String topmV3FileS = "M:\\pipelineTest\\PanA\\alignment\\tagGWASMap.v3.topm.h5";
        String arguments = "-i " + samFileS + " -t " + tagGWASMapFileS +  " -o " + topmV3FileS;
        String[] args = arguments.split(" ");
        PanASamToMultiPositionTOPMPlugin stt = new PanASamToMultiPositionTOPMPlugin();
        stt.setParameters(args);
        stt.performFunction(null);
    }
    
    public void addPosToTagMapPlugin () {
        String tagGWASMapFileS = "M:\\pipelineTest\\PanA\\tagMap\\tagGWASMap.h5";
        String topmV3FileS = "M:\\pipelineTest\\PanA\\alignment\\tagGWASMap.v3.topm.h5";
        String arguments = "-i " + tagGWASMapFileS +  " -t " + topmV3FileS;
        String[] args = arguments.split(" ");
        PanAAddPosToTagMapPlugin  p = new PanAAddPosToTagMapPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void buildTrainingSetPlugin () {
        String tagGWASMapFileS = "M:\\pipelineTest\\PanA\\tagMap\\tagGWASMap.h5";
        String trainingSetFileS = "M:\\pipelineTest\\PanA\\training\\uniqueRefTrain.arff";
        String rScriptPath = "C:\\Users\\fl262\\Documents\\R\\R-3.0.2\\bin\\Rscript.exe";
        String boxcoxParemeterFileS = "M:\\pipelineTest\\PanA\\training\\boxcoxParemeter.txt";
        String arguments = "-m " + tagGWASMapFileS +  " -t " + trainingSetFileS + " -r " + rScriptPath +  " -b " + boxcoxParemeterFileS;
        String[] args = arguments.split(" ");
        PanABuildTrainingSetPlugin p = new PanABuildTrainingSetPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void modelTrainingPlugin () {
        String trainingSetFileS = "M:\\pipelineTest\\PanA\\training\\uniqueRefTrain.arff";
        String wekaPath = "E:\\Database\\Weka-3-6\\weka.jar";
        String modelFileS = "M:\\pipelineTest\\PanA\\training\\m5.mod";
        String trainingReportDirS = "M:\\pipelineTest\\PanA\\training\\report\\";
        String arguments = "-t " + trainingSetFileS +  " -w " + wekaPath +  " -m " + modelFileS + " -r " + trainingReportDirS;
        String[] args = arguments.split(" ");
        PanAModelTrainingPlugin p = new PanAModelTrainingPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void predictionPlugin () {
        String wekaPath = "E:\\Database\\Weka-3-6\\weka.jar";
        String tagGWASMapFileS = "M:\\pipelineTest\\PanA\\tagMap\\tagGWASMap.h5";
        String modelFileS = "M:\\pipelineTest\\PanA\\training\\m5.mod";
        String boxcoxParemeterFileS = "M:\\pipelineTest\\PanA\\training\\boxcoxParemeter.txt";
        String arguments = "-t " + tagGWASMapFileS +  " -w " + wekaPath +  " -m " + modelFileS + " -b " + boxcoxParemeterFileS;
        String[] args = arguments.split(" ");
        PanAPredictionPlugin p = new PanAPredictionPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void filterTagMapPlugin () {
        String tagGWASMapFileS = "M:\\pipelineTest\\PanA\\tagMap\\tagGWASMap.h5";
        int distanceCutoff = Integer.MAX_VALUE;
        String anchorFileS = "M:\\pipelineTest\\PanA\\togm\\anchor.txt" ;
        String arguments = "-t " + tagGWASMapFileS +  " -a " + anchorFileS +  " -c " + String.valueOf(distanceCutoff);
        String[] args = arguments.split(" ");
        PanAFilteringTagMapPlugin p = new PanAFilteringTagMapPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public static void main (String[] args) {
        new PanAUsageExample();
    }
}
