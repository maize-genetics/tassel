package net.maizegenetics.dna.snp.io;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import org.junit.Assert;
import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * @author Ed Buckler
 */
public class BuilderFromGenotypeHDF5Test {

    @Test
    public void testMergeToMutableHDF5() throws Exception {
        //chrom3 lists the same taxa as the others in a different order, so the merge has to remap taxa
        String[] hapmaps = {"mdp_genotype_chrom3_taxa_diff_order.hmp.txt", "mdp_genotype_chrom8.hmp.txt", "mdp_genotype_chrom10.hmp.txt"};
        GenotypeTable[] sources = new GenotypeTable[hapmaps.length];
        String[] infiles = new String[hapmaps.length];
        int expectedSites = 0;
        for (int i = 0; i < hapmaps.length; i++) {
            sources[i] = BuilderFromHapMap.getBuilder(TutorialConstants.TUTORIAL_DIR + hapmaps[i]).build();
            infiles[i] = GeneralConstants.TEMP_DIR + "mergeToMutableSource" + i + ".h5";
            Files.deleteIfExists(Paths.get(infiles[i]));
            GenotypeTableBuilder.getInstance(sources[i], infiles[i]);
            expectedSites += sources[i].numberOfSites();
        }
        String outfile = GeneralConstants.TEMP_DIR + "mergeToMutable.h5";
        Files.deleteIfExists(Paths.get(outfile));

        BuilderFromGenotypeHDF5.mergeToMutableHDF5(infiles, outfile);

        GenotypeTable merged = GenotypeTableBuilder.getInstance(outfile);
        Assert.assertEquals("Merged site count", expectedSites, merged.numberOfSites());
        Assert.assertEquals("Merged taxa count", sources[0].numberOfTaxa(), merged.numberOfTaxa());
        for (GenotypeTable source : sources) {
            int[] mergedTaxon = new int[source.numberOfTaxa()];
            for (int t = 0; t < source.numberOfTaxa(); t++) {
                mergedTaxon[t] = merged.taxa().indexOf(source.taxa().get(t));
                Assert.assertTrue("Taxon missing from merge: " + source.taxaName(t), mergedTaxon[t] >= 0);
            }
            for (int s = 0; s < source.numberOfSites(); s++) {
                int mergedSite = merged.positions().indexOf(source.positions().get(s));
                Assert.assertTrue("Site missing from merge: " + source.positions().get(s), mergedSite >= 0);
                for (int t = 0; t < source.numberOfTaxa(); t++) {
                    Assert.assertEquals("Genotype differs for " + source.taxaName(t) + " at " + source.positions().get(s),
                            source.genotype(t, s), merged.genotype(mergedTaxon[t], mergedSite));
                }
            }
        }
    }
}
