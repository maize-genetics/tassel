/*
 *  MaskGenotypePlugin
 * 
 *  Created on May 8, 2015
 */
package net.maizegenetics.analysis.data;

import com.google.common.collect.Range;
import java.awt.Frame;
import java.net.URL;
import javax.swing.*;
import java.util.List;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.MaskMatrixBuilder;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.score.AlleleDepth;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * The purpose of this class is to mask (make UNKNOWN) a portion of the
 * genotypes of the given GenotypeTable. The resulting GenotypeTable will be the
 * same except some genotypes will have been changed to Unknown.
 *
 * @author Terry Casstevens
 */
public class MaskGenotypePlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(MaskGenotypePlugin.class);

    private PluginParameter<Double> myPercentageMasked
            = new PluginParameter.Builder<>("percentageMasked", 0.01, Double.class)
            .range(Range.openClosed(0.0, 1.0))
            .description("Percentage of genotypes (not already unknown) to mask.")
            .build();

    private PluginParameter<Integer> myMinDepth
            = new PluginParameter.Builder<>("minDepth", 0, Integer.class)
            .range(Range.atLeast(0))
            .description("Minimum depth required before masking.")
            .build();

    public MaskGenotypePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> data = input.getDataOfType(GenotypeTable.class);
        if (data.size() != 1) {
            throw new IllegalArgumentException("MaskGenotypePlugin: preProcessParameters: must input 1 GenotypeTable.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        Datum inputDatum = input.getDataOfType(GenotypeTable.class).get(0);
        GenotypeTable original = (GenotypeTable) inputDatum.getData();

        if ((minDepth() > 0) && !original.hasDepth()) {
            throw new IllegalArgumentException("MaskGenotypePlugin: processData: input doesn't have depth information and you set minimum depth to " + minDepth());
        }

        MaskMatrixBuilder builder = MaskMatrixBuilder.getInstance(original.numberOfTaxa(), original.numberOfSites(), true);

        long numGenotypesEligibleToMask = markEligibleGenotypes(builder, original);

        if (numGenotypesEligibleToMask == 0) {
            DialogUtils.showWarning("No Genotypes match your criteria to be masked.", getParentFrame());
        }

        long numberActuallyMasked = builder.reduceMaskTo(percentageMasked());

        myLogger.info("Number of Genotypes Masked: " + numberActuallyMasked);

        GenotypeTable result = GenotypeTableBuilder.getInstance(original, builder.build());

        Datum genotype = new Datum(inputDatum.getName() + "_Masked", result, null);

        return new DataSet(genotype, this);

    }

    private long markEligibleGenotypes(MaskMatrixBuilder builder, GenotypeTable orig) {
        if (orig.genotypeMatrix().isSiteOptimized()) {
            return siteMarkEligibleGenotypes(builder, orig);
        } else {
            return taxaMarkEligibleGenotypes(builder, orig);
        }
    }

    private long siteMarkEligibleGenotypes(MaskMatrixBuilder builder, GenotypeTable orig) {

        GenotypeCallTable origCalls = orig.genotypeMatrix();

        int numTaxa = origCalls.numberOfTaxa();
        int numSites = origCalls.numberOfSites();

        myLogger.info("Number of taxa: " + numTaxa);
        myLogger.info("Number of sites: " + numSites);
        myLogger.info("Number of genotypes: " + ((long) numTaxa * (long) numSites));

        long numGenotypesEligibleToMask = 0;

        if (minDepth() > 0) {
            AlleleDepth depth = orig.depth();
            for (int s = 0; s < numSites; s++) {
                byte[] genotypes = origCalls.genotypeForAllTaxa(s);
                for (int t = 0; t < numTaxa; t++) {
                    if (genotypes[t] != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                        byte[] alleles = GenotypeTableUtils.getDiploidValues(genotypes[t]);
                        int currentDepth = 0;
                        if (alleles[0] < 6) {
                            currentDepth = depth.value(t, s, AlleleDepth.ALLELE_DEPTH_TYPES[alleles[0]]);
                        }
                        if (alleles[1] < 6) {
                            currentDepth += depth.value(t, s, AlleleDepth.ALLELE_DEPTH_TYPES[alleles[1]]);
                        }
                        if (currentDepth >= minDepth()) {
                            builder.set(t, s);
                            numGenotypesEligibleToMask++;
                        }
                    }
                }
            }
        } else {
            for (int s = 0; s < numSites; s++) {
                byte[] genotypes = origCalls.genotypeForAllTaxa(s);
                for (int t = 0; t < numTaxa; t++) {
                    if (genotypes[t] != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                        builder.set(t, s);
                        numGenotypesEligibleToMask++;
                    }
                }
            }
        }

        myLogger.info("Number of genotypes eligible to mask: " + numGenotypesEligibleToMask);

        return numGenotypesEligibleToMask;

    }

    private long taxaMarkEligibleGenotypes(MaskMatrixBuilder builder, GenotypeTable orig) {

        GenotypeCallTable origCalls = orig.genotypeMatrix();

        int numTaxa = origCalls.numberOfTaxa();
        int numSites = origCalls.numberOfSites();

        myLogger.info("Number of taxa: " + numTaxa);
        myLogger.info("Number of sites: " + numSites);
        myLogger.info("Number of genotypes: " + ((long) numTaxa * (long) numSites));

        long numGenotypesEligibleToMask = 0;

        if (minDepth() > 0) {
            AlleleDepth depth = orig.depth();
            for (int t = 0; t < numTaxa; t++) {
                byte[] genotypes = origCalls.genotypeForAllSites(t);
                for (int s = 0; s < numSites; s++) {
                    if (genotypes[s] != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                        byte[] alleles = GenotypeTableUtils.getDiploidValues(genotypes[s]);
                        int currentDepth = 0;
                        if (alleles[0] < 6) {
                            currentDepth = depth.value(t, s, AlleleDepth.ALLELE_DEPTH_TYPES[alleles[0]]);
                        }
                        if (alleles[1] < 6) {
                            currentDepth += depth.value(t, s, AlleleDepth.ALLELE_DEPTH_TYPES[alleles[1]]);
                        }
                        if (currentDepth >= minDepth()) {
                            builder.set(t, s);
                            numGenotypesEligibleToMask++;
                        }
                    }
                }
            }
        } else {
            for (int t = 0; t < numTaxa; t++) {
                byte[] genotypes = origCalls.genotypeForAllSites(t);
                for (int s = 0; s < numSites; s++) {
                    if (genotypes[s] != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                        builder.set(t, s);
                        numGenotypesEligibleToMask++;
                    }
                }
            }
        }

        myLogger.info("Number of genotypes eligible to mask: " + numGenotypesEligibleToMask);

        return numGenotypesEligibleToMask;

    }

    /**
     * Convenience method to run plugin with one return object.
     */
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Percentage of genotypes to mask.
     *
     * @return Percentage Masked
     */
    public Double percentageMasked() {
        return myPercentageMasked.value();
    }

    /**
     * Set Percentage Masked. Percentage of genotypes to mask.
     *
     * @param value Percentage Masked
     *
     * @return this plugin
     */
    public MaskGenotypePlugin percentageMasked(Double value) {
        myPercentageMasked = new PluginParameter<>(myPercentageMasked, value);
        return this;
    }

    /**
     * Minimum depth required before masking.
     *
     * @return Min Depth
     */
    public Integer minDepth() {
        return myMinDepth.value();
    }

    /**
     * Set Min Depth. Minimum depth required before masking.
     *
     * @param value Min Depth
     *
     * @return this plugin
     */
    public MaskGenotypePlugin minDepth(Integer value) {
        myMinDepth = new PluginParameter<>(myMinDepth, value);
        return this;
    }

    @Override
    public String getToolTipText() {
        return "Mask Genotype";
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = MaskGenotypePlugin.class.getResource("/net/maizegenetics/analysis/images/mask.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Mask Genotype";
    }

}
