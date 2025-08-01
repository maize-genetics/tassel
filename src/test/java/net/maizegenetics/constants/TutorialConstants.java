/*
 * TutorialConstants
 */
package net.maizegenetics.constants;

/**
 *
 * @author terry
 */
public class TutorialConstants {

    public static final String TUTORIAL_DIR = GeneralConstants.DATA_DIR + "Tutorial/";
    public static final String HAPMAP_FILENAME = TUTORIAL_DIR + "mdp_genotype.hmp.txt";
    public static final String HAPMAP_CHR_9_10_FILENAME = TUTORIAL_DIR + "mdp_genotype_chr9_10.hmp.txt";
    public static final String PHYLIP_FILENAME = TUTORIAL_DIR + "d8_sequence.phy";
    public static final String POP_STRUCTURE_FILENAME = TUTORIAL_DIR + "mdp_population_structure.txt";
    public static final String TRAITS_FILENAME = TUTORIAL_DIR + "mdp_traits.txt";
    public static final String HAPMAP_BGZIP_FILENAME = TUTORIAL_DIR + "mdp_genotype.hmp.txt.gz";
    public static final String HAPMAP_LIX_FILENAME = TUTORIAL_DIR + "mdp_genotype.hmp.txt.gz.lix";

    private TutorialConstants() {
        // utility
    }
}
