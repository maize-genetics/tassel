/*
 * SetLowDepthGenosToMissingPlugin
 */
package net.maizegenetics.analysis.data;

import com.google.common.collect.Range;
import java.awt.Frame;
import java.net.URL;
import java.util.List;
import javax.swing.ImageIcon;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.MaskMatrixBuilder;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.Taxon;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Set genotypes below a minimum depth to missing. Set each genotype in the
 * input genotypes to missing if the underlying allelic depth is below a user-
 * specified minimum. Input: GenotypeTable stored as a Datum within a DataSet
 * Output: GenotypeTable stored as a Datum within a DataSet
 *
 * @author Christopher Bottoms
 * @author Jeff Glaubitz
 */
public class SetLowDepthGenosToMissingPlugin extends net.maizegenetics.plugindef.AbstractPlugin {

    @Override
    public String pluginDescription() {
        return "Sets each genotype in the input genotypes to missing if "
                + "the underlying allelic depth is below a specified minimum.";
    }

    private PluginParameter<Integer> minDepth
            = new PluginParameter.Builder<>("minDepth", null, Integer.class)
            .guiName("minimum genotype depth")
            .required(true)
            .description("Minimum depth, below which genotypes are set to missing. Must be between 2 and 127, inclusive.")
            .range(Range.closed(2, 127))
            .build();

    private static final Logger myLogger = LogManager.getLogger(SetLowDepthGenosToMissingPlugin.class);
    private GenotypeTable inputGenotypes = null;
    private String inputGenosName = null;

    public SetLowDepthGenosToMissingPlugin() {
        super(null, false);
    }

    public SetLowDepthGenosToMissingPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        if (input == null) {
            throw new IllegalArgumentException(
                    "SetLowDepathGenotypesToMissingPlugin: preProcessParameters: Please select one Genotype Table.");
        }
        List<Datum> genotypeTables = input.getDataOfType(GenotypeTable.class);

        inputGenosName = genotypeTables.get(0).getName();

        myLogger.info("\n" + pluginDescription() + "\n");

        myLogger.info("Input genotype name: " + inputGenosName);

        if (genotypeTables.size() == 1) {
            inputGenotypes
                    = (GenotypeTable) genotypeTables.get(0).getData();
            if (!inputGenotypes.hasDepth()) {
                throw new IllegalArgumentException("SetLowDepthGenotypesToMissingPlugin: preProcessParameters: Please select a Genotype Table with allele depth information.");
            }
        } else {
            throw new IllegalArgumentException(
                    "SetLowDepthGenotypesToMissingPlugin: preProcessParameters: Please select one Genotype Table.");
        }
    }

    /**
     * Main method for this plugin.
     *
     * @param input DataSet object containing the input GenotypeTable.
     * @return DataSet object containing the output GenotypeTable where
     * genotypes with depth below the specified minimum are set to missing.
     */
    @Override
    public DataSet processData(DataSet input) {
        int numberOfSites = inputGenotypes.numberOfSites();
        MaskMatrixBuilder mgtb = MaskMatrixBuilder.getInstance(inputGenotypes.numberOfTaxa(), inputGenotypes.numberOfSites(), true);
        int nGenosSetToMissing = 0;
        int nOrigNonMissing = 0;
        for (int siteIndex = 0; siteIndex < numberOfSites; siteIndex++) {
            nOrigNonMissing += inputGenotypes.totalNonMissingForSite(siteIndex);
            int[][] allelesByFreq = inputGenotypes.allelesSortedByFrequency(siteIndex);
            int numAlleles = allelesByFreq[0].length;
            int[] allelesAtSite = new int[numAlleles];

            for (int alleleIndex = 0; alleleIndex < numAlleles; alleleIndex++) {
                allelesAtSite[alleleIndex] = allelesByFreq[0][alleleIndex];
            }

            for (Taxon inTaxon : inputGenotypes.taxa()) {
                int taxonIndex = inputGenotypes.taxa().indexOf(inTaxon);
                int[] alleleDepths = inputGenotypes.depthForAlleles(taxonIndex, siteIndex);
                int depthSum = 0;
                for (int allele : allelesAtSite) {
                    depthSum += alleleDepths[allele];
                }
                if (depthSum < minDepth()) {
                    mgtb.set(taxonIndex, siteIndex);
                    nGenosSetToMissing++;
                }
            }
        }

        String outGenosName = inputGenosName + "MinDepth" + minDepth();
        GenotypeTable outGenos = GenotypeTableBuilder.getInstance(inputGenotypes, mgtb.build());
        double percentOfTotalSetToMissing = (double) 100 * nGenosSetToMissing / (numberOfSites * inputGenotypes.numberOfTaxa());
        double percentOfNonMissingSetToMissing = (double) 100 * nGenosSetToMissing / nOrigNonMissing;
        myLogger.info("\n" + nGenosSetToMissing + " genotypes with a depth less than "
                + minDepth() + " were set to missing \n"
                + "   = " + String.format("%,.2f", percentOfTotalSetToMissing) + "% of the original genotypes"
                + " (nSites x nSamples), or "
                + String.format("%,.2f", percentOfNonMissingSetToMissing) + "% of the original "
                + "nonMissing genotypes.\n");
        return new DataSet(new Datum(outGenosName, outGenos, null), null);
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = SetLowDepthGenosToMissingPlugin.class
                .getResource("/net/maizegenetics/analysis/images/lowDepthToMissing.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "SetLowDepthGenosToMissing";
    }

    @Override
    public String getToolTipText() {
        return "Set genotypes below a minimum depth to missing";
    }

    @Override
    public String getCitation() {
        return "Christopher Bottoms, Jeff Glaubitz (2015) TASSEL Hackathon Oct 2015";
    }

    /**
     * Convenience method to run plugin with input and output GenotypeTable
     * objects (rather than DataSets)
     *
     * @param inputGenos Input GenotypeTable.
     * @return GenotypeTable where genotypes with depth below the specified
     * minimum are set to missing.
     */
    public GenotypeTable runPlugin(GenotypeTable inputGenos) {
        DataSet input = new DataSet(new Datum("inputGenotypes", inputGenos, null), null);
        return runPlugin(input);
    }

    /**
     * Convenience method to run plugin and output a GenotypeTable object
     * (rather than a DataSet)
     *
     * @param input DataSet object containing the input GenotypeTable.
     * @return GenotypeTable where genotypes with depth below the specified
     * minimum are set to missing.
     */
    private GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Min Depth
     *
     * @return Min Depth
     */
    public Integer minDepth() {
        return minDepth.value();
    }

    /**
     * Set Min Depth. Min Depth
     *
     * @param value Min Depth
     *
     * @return this plugin
     */
    public SetLowDepthGenosToMissingPlugin minDepth(Integer value) {
        minDepth = new PluginParameter<>(minDepth, value);
        return this;
    }

}
