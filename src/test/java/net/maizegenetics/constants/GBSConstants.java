/*
 * GBSConstants
 */
package net.maizegenetics.constants;

/**
 *
 * @author terry
 */
public class GBSConstants {

    // Input Data
    public static final String GBS_DATA_DIR = "dataFiles/GBS/";
    public static final String RAW_SEQ_CHR_9_10_200000 = "Chr9_10-200000/";
    public static final String RAW_SEQ_CHR_9_10_20000000 = "Chr9_10-20000000/";
    public static final String RAW_SEQ_CURRENT_TEST = RAW_SEQ_CHR_9_10_20000000;
    public static final String GBS_INPUT_DIR = GBS_DATA_DIR + RAW_SEQ_CURRENT_TEST;
    public static final String GBS_INPUT_FASTQ_FILE = GBS_INPUT_DIR + "C05F2ACXX_5_fastq.gz";
    public static final String GBS_REFERENCE_GENOME = GBS_DATA_DIR + "ZmB73_RefGen_v2_chr9_10_1st20MB.fasta";
    public static final String GBS_INPUT_TOPM = "dataFiles/ExpectedResults/GBS/" + RAW_SEQ_CHR_9_10_20000000 + "DiscoverySNPCallerPlugin/chr910.topm.h5";
    // Temp Data
    public static final String GBS_TEMP_DIR = "tempDir/GBS/" + RAW_SEQ_CURRENT_TEST;
    public static final String GBS_TEMP_FASTQ_TO_TAG_COUNT_PLUGIN_DIR = GBS_TEMP_DIR + "FastqToTagCountPlugin/";
    public static final String GBS_TEMP_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_DIR = GBS_TEMP_DIR + "MergeMultipleTagCountPlugin/";
    public static final String GBS_TEMP_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_FILE = GBS_TEMP_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_DIR + "Merged_Tag_Counts.cnt";
    public static final String GBS_TEMP_TAG_COUNT_TO_FASTQ_PLUGIN_DIR = GBS_TEMP_DIR + "TagCountToFastqPlugin/";
    public static final String GBS_TEMP_TAG_COUNT_TO_FASTQ_PLUGIN_FILE = GBS_TEMP_TAG_COUNT_TO_FASTQ_PLUGIN_DIR + "Fastq_For_Aligner.fq.gz";
    public static final String GBS_TEMP_SAM_CONVERTER_PLUGIN_DIR = GBSConstants.GBS_TEMP_DIR + "SAMConverterPlugin/";
    public static final String GBS_TEMP_SAM_CONVERTER_PLUGIN_FILE = GBSConstants.GBS_TEMP_SAM_CONVERTER_PLUGIN_DIR + "TOPM_from_SAM.topm";
    public static final String GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_DIR = GBSConstants.GBS_TEMP_DIR + "SeqToTBTHDF5Plugin/";
    public static final String GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_FILE = GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_DIR + "TBT_from_Raw_Seq.h5";
    public static final String GBS_TEMP_MODIFY_TBT_HDF5_PLUGIN_DIR = GBSConstants.GBS_TEMP_DIR + "ModifyTBTHDF5Plugin/";
    public static final String GBS_TEMP_MODIFY_TBT_HDF5_PLUGIN_PIVOTED_FILE = GBSConstants.GBS_TEMP_MODIFY_TBT_HDF5_PLUGIN_DIR + "TBT_Pivoted.h5";
    public static final String GBS_DIRECT_TBT_FILE = GBSConstants.GBS_TEMP_MODIFY_TBT_HDF5_PLUGIN_DIR + "TBT_Direct.h5";
    public static final String GBS_TEMP_DISCOVERY_SNP_CALLER_PLUGIN_DIR = GBSConstants.GBS_TEMP_DIR + "DiscoverySNPCallerPlugin/";
    public static final String GBS_TEMP_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE = GBSConstants.GBS_TEMP_DISCOVERY_SNP_CALLER_PLUGIN_DIR + "TOPM_with_Variants.topm";
    public static final String GBS_TEMP_DISCOVERY_SNP_CALLER_PLUGIN_HAPMAP_OUT_FILE = GBSConstants.GBS_TEMP_DISCOVERY_SNP_CALLER_PLUGIN_DIR + "GBSGenos.chr+.hmp.txt";
    public static final String GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_DIR = GBSConstants.GBS_TEMP_DIR + "ProductionSNPCallerPlugin/";
    public static final String GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_OUT_FILE = GBSConstants.GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_DIR + "PipelineTestingGenos.h5";
    public static final String GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_UNFINISHED_FILE = GBSConstants.GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_DIR + "PipelineTestingGenosKO.h5";
    public static final String GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_FINISHED_FILE = GBSConstants.GBS_TEMP_PRODUCTION_SNP_CALLER_PLUGIN_DIR + "PipelineTestingGenos__DATE__.h5";
    public static final String GBS_TEMP_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_DIR = GBSConstants.GBS_TEMP_DIR + "AddReferenceAlleleToHDF5Plugin/";
    public static final String GBS_TEMP_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_OUT_FILE = GBSConstants.GBS_TEMP_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_DIR + "PipelineTestingGenos_withRef.h5";
    public static final String GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_DIR = GBSConstants.GBS_TEMP_DIR + "ReimputeUpdatedTaxaByFillinPlugin/";
    public static final String GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_IMP_FILE = GBSConstants.GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_DIR + "PipelineTestingGenosWithRefImp.h5";
    public static final String GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_IMP_FILE2 = GBSConstants.GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_DIR + "PipelineTestingGenosWithRefImpCopy.h5";
    public static final String GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_RAW_FILE = GBSConstants.GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_DIR + "PipelineTestingGenosWithRef.h5";
    public static final String GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_MOD_FILE = GBSConstants.GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_DIR + "PipelineTestingGenosWithRefMod.h5";
    // Expected Data
    public static final String GBS_EXPECTED_DIR = "dataFiles/ExpectedResults/GBS/" + RAW_SEQ_CURRENT_TEST;
    public static final String GBS_EXPECTED_FASTQ_TO_TAG_COUNT_PLUGIN_DIR = GBS_EXPECTED_DIR + "FastqToTagCountPlugin/";
    public static final String GBS_EXPECTED_MERGE_MULTIPLE_TAG_COUNT_PLUGIN_FILE = GBS_EXPECTED_DIR + "MergeMultipleTagCountPlugin/" + "Merged_Tag_Counts.cnt";
    public static final String GBS_EXPECTED_TAG_COUNT_TO_FASTQ_PLUGIN_FILE = GBS_EXPECTED_DIR + "TagCountToFastqPlugin/" + "Fastq_For_Aligner.fq.gz";
    public static final String GBS_EXPECTED_BOWTIE_SAM_FILE = GBS_EXPECTED_DIR + "bowtie2.1/Aligned_B73.sam";
    public static final String GBS_EXPECTED_SAM_CONVERTER_PLUGIN_FILE = GBSConstants.GBS_EXPECTED_DIR + "SAMConverterPlugin/TOPM_from_SAM.topm";
    public static final String GBS_EXPECTED_SEQ_TO_TBT_HDF5_PLUGIN_FILE = GBSConstants.GBS_EXPECTED_DIR + "SeqToTBTHDF5Plugin/TBT_from_Raw_Seq.h5";
    public static final String GBS_EXPECTED_SEQ_TO_TBT_HDF5_PLUGIN__MIRROR_FILE = GBSConstants.GBS_EXPECTED_DIR + "SeqToTBTHDF5Plugin/TBT_from_Raw_Seq_Mirror.h5";
    public static final String GBS_EXPECTED_MODIFY_TBT_HDF5_PLUGIN_MERGE_FILE = GBSConstants.GBS_EXPECTED_DIR + "ModifyTBTHDF5Plugin/TBT_after_Merge.h5";
    public static final String GBS_EXPECTED_MODIFY_TBT_HDF5_PLUGIN_MERGE_TAXA_FILE = GBSConstants.GBS_EXPECTED_DIR + "ModifyTBTHDF5Plugin/TBT_after_Merge_Taxa.h5";
    public static final String GBS_EXPECTED_MODIFY_TBT_HDF5_PLUGIN_PIVOTED_FILE = GBSConstants.GBS_EXPECTED_DIR + "ModifyTBTHDF5Plugin/TBT_Pivoted.h5";
    public static final String GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_DIR = GBSConstants.GBS_EXPECTED_DIR + "DiscoverySNPCallerPlugin/";
    public static final String GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_TOPM_OUT_FILE = GBSConstants.GBS_EXPECTED_DISCOVERY_SNP_CALLER_PLUGIN_DIR + "TOPM_with_Variants.topm";
    public static final String GBS_EXPECTED_PRODUCTION_SNP_CALLER_PLUGIN_DIR = GBSConstants.GBS_EXPECTED_DIR + "ProductionSNPCallerPlugin/";
    public static final String GBS_EXPECTED_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_OUT_FILE = GBSConstants.GBS_EXPECTED_PRODUCTION_SNP_CALLER_PLUGIN_DIR + "PipelineTestingGenos.h5";
    public static final String GBS_EXPECTED_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_OUT_FILE = GBSConstants.GBS_EXPECTED_PRODUCTION_SNP_CALLER_PLUGIN_DIR + "PipelineTestingGenos_withRef.h5";
    public static final String GBS_FILLIN_DONOR_DIR = GBSConstants.GBS_EXPECTED_DIR + "maizeTestDataDonorHaplos/";
    // Key File
    public static final String GBS_TESTING_KEY_FILE = GBS_DATA_DIR + "Pipeline_Testing_key.txt";

    //GBS2 Files
    public static final String GBS_GBS2DB_FILE = GBSConstants.GBS_TEMP_DIR + "GBSv2.db";
    public static final String SAM_FLAG_EXPECTED_RESULTS_FILE = GBS_DATA_DIR + "SAMFlagExpectedResults.txt";

    private GBSConstants() {
        // utility
    }
}
