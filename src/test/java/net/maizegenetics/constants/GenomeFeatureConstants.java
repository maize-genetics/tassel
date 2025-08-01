package net.maizegenetics.constants;

/**
 * Created by jgw87 on 7/11/14.
 */
public class GenomeFeatureConstants {

    public static final String GENOME_FEATURE_DIR = GeneralConstants.DATA_DIR + "GenomeFeatures/";
    public static final String GTF_FILENAME = GENOME_FEATURE_DIR + "features_input1.gtf";       //GTF input format
    public static final String GFF_FILENAME = GENOME_FEATURE_DIR + "features_input2.gff";       //GFF input format
    public static final String JSON_FILENAME = GENOME_FEATURE_DIR + "features_input3.json";     //JSON input format
    public static final String FLATFILE_FILENAME = GENOME_FEATURE_DIR + "features_input4.txt";  //Tab-delimited input format

    private GenomeFeatureConstants() {
        // not needed
    }
}
