package net.maizegenetics.analysis.imputation;

import net.maizegenetics.dna.snp.*;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportBuilder;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * This evaluates the accuracy of an imputed genotype. First input genotype
 * should be the original. Second should be the masked genotype. Third should be
 * the imputed genotype.
 *
 * @author Ed Buckler
 */
public class ImputationAccuracyPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(ImputationAccuracyPlugin.class);

	private PluginParameter<Boolean> isSubset = new PluginParameter.Builder<>("subset", false, Boolean.class)
            .guiName("Imputed is subset of original")
            .description("The imputed data does not contain all of the sites or taxa in the original data.").build();

    public ImputationAccuracyPlugin() {
        super(null, false);
    }

    public ImputationAccuracyPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
        if (alignInList.size() != 3) {
            throw new IllegalArgumentException("ImputationAccuracyPlugin: preProcessParameters: Please select three Genotype Table.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        if (isSubset.value()) return processDataForSubset(input);

        // Load in the genotype table
        GenotypeTable origGenoTable = (GenotypeTable) input.getDataOfType(GenotypeTable.class).get(0).getData();
        myLogger.info("Original Genotype: " + input.getDataOfType(GenotypeTable.class).get(0).getName());
        GenotypeTable maskGenoTable = (GenotypeTable) input.getDataOfType(GenotypeTable.class).get(1).getData();
        myLogger.info("Masked Genotype: " + input.getDataOfType(GenotypeTable.class).get(1).getName());
        GenotypeTable impGenoTable = (GenotypeTable) input.getDataOfType(GenotypeTable.class).get(2).getData();
        myLogger.info("Imputed Genotype: " + input.getDataOfType(GenotypeTable.class).get(2).getName());

        int[][] cnts = new int[5][5];
        for (int site = 0; site < origGenoTable.numberOfSites(); site++) {
            byte majorAllele = origGenoTable.majorAllele(site);
            byte minorAllele = origGenoTable.minorAllele(site);
            Map<Byte, Integer> genotypeToIndexMap = new HashMap<>();
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(majorAllele, majorAllele), 0);
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(majorAllele, minorAllele), 1);
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(minorAllele, majorAllele), 1);
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(minorAllele, minorAllele), 2);
            genotypeToIndexMap.put(GenotypeTable.UNKNOWN_DIPLOID_ALLELE, 3);
            for (int taxonIdx = 0; taxonIdx < origGenoTable.numberOfTaxa(); taxonIdx++) {
                byte origGeno = origGenoTable.genotype(taxonIdx, site);
                if (origGeno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    continue;  //skip if unknown in the original data
                }
                if (maskGenoTable.genotype(taxonIdx, site) != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    continue;  //skip if not missing in the masked data
                }
                int originalIndex = genotypeToIndexMap.getOrDefault(origGeno, 4);
                int impIndex = genotypeToIndexMap.getOrDefault(impGenoTable.genotype(taxonIdx, site), 4);
                if (impIndex == 4) {
                    System.out.println(site + ":" + impGenoTable.taxa().get(taxonIdx).toString() + ":" + impGenoTable.genotype(taxonIdx, site) + NucleotideAlignmentConstants.getNucleotideIUPAC(impGenoTable.genotype(taxonIdx, site)));
                }
                cnts[originalIndex][impIndex]++;
            }
        }

        DataSet dataSet = new DataSet(new Datum("AccuracyReport", makeTableReport(cnts), ""), this);
        return dataSet;
    }

    private DataSet processDataForSubset(DataSet input) {
        // Load in the genotype table
        GenotypeTable origGenoTable = (GenotypeTable) input.getDataOfType(GenotypeTable.class).get(0).getData();
        myLogger.info("Original Genotype: " + input.getDataOfType(GenotypeTable.class).get(0).getName());
        GenotypeTable maskGenoTable = (GenotypeTable) input.getDataOfType(GenotypeTable.class).get(1).getData();
        myLogger.info("Masked Genotype: " + input.getDataOfType(GenotypeTable.class).get(1).getName());
        GenotypeTable impGenoTable = (GenotypeTable) input.getDataOfType(GenotypeTable.class).get(2).getData();
        myLogger.info("Imputed Genotype: " + input.getDataOfType(GenotypeTable.class).get(2).getName());

        int[][] cnts = new int[5][5];
        for (int site = 0; site < origGenoTable.numberOfSites(); site++) {
        		int imputedSite = impGenoTable.positions().siteOfPhysicalPosition(origGenoTable.chromosomalPosition(site), origGenoTable.chromosome(site));
        		if (imputedSite < 0) continue;
            byte majorAllele = origGenoTable.majorAllele(site);
            byte minorAllele = origGenoTable.minorAllele(site);
            Map<Byte, Integer> genotypeToIndexMap = new HashMap<>();
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(majorAllele, majorAllele), 0);
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(majorAllele, minorAllele), 1);
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(minorAllele, majorAllele), 1);
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(minorAllele, minorAllele), 2);
            genotypeToIndexMap.put(GenotypeTable.UNKNOWN_DIPLOID_ALLELE, 3);
            for (int taxonIdx = 0; taxonIdx < origGenoTable.numberOfTaxa(); taxonIdx++) {
            		int imputedTaxonIdx = impGenoTable.taxa().indexOf(origGenoTable.taxa().get(taxonIdx));
            		if (imputedTaxonIdx < 0) continue;
                byte origGeno = origGenoTable.genotype(taxonIdx, site);
                if (origGeno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    continue;  //skip if unknown in the original data
                }
                if (maskGenoTable.genotype(taxonIdx, site) != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    continue;  //skip if not missing in the masked data
                }
                int originalIndex = genotypeToIndexMap.getOrDefault(origGeno, 4);
                
                try {
                		byte impGeno = impGenoTable.genotype(imputedTaxonIdx, imputedSite);
                } catch (Exception e) {
                		System.out.printf("Error at orig site = %d, imputed site = %d%n", site, imputedSite);
                }
                
                int impIndex = genotypeToIndexMap.getOrDefault(impGenoTable.genotype(imputedTaxonIdx, imputedSite), 4);
                if (impIndex == 4) {
                    System.out.println(site + ":" + impGenoTable.taxa().get(taxonIdx).toString() + ":" + impGenoTable.genotype(imputedTaxonIdx, imputedSite) + NucleotideAlignmentConstants.getNucleotideIUPAC(impGenoTable.genotype(imputedTaxonIdx, imputedSite)));
                }
                cnts[originalIndex][impIndex]++;
            }
        }

        DataSet dataSet = new DataSet(new Datum("AccuracyReport", makeTableReport(cnts), ""), this);
        return dataSet;
    }
    
    private TableReport makeTableReport(int[][] cnts) {
        String[] headers = {"Original/Imputed", "AA", "Aa", "aa", "N", "Other"};
        TableReportBuilder reportBuilder = TableReportBuilder.getInstance("ImputationAccuracy", headers);
        int errors = 0, correct = 0;
        for (int i = 0; i < cnts.length; i++) {
            reportBuilder.addElements(headers[i + 1], ArrayUtils.toObject(cnts[i]));
            correct += cnts[i][i];
            errors += cnts[i][0] + cnts[i][1] + cnts[i][2] - cnts[i][i];
        }
        double errorRate = (double) errors / (double) (correct + errors);
        reportBuilder.addElements("Correct", correct);
        reportBuilder.addElements("Errors", errors);
        reportBuilder.addElements("ErrorRate", errorRate);
        return reportBuilder.build();
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Evaluate Imputation Accuracy";
    }

    @Override
    public String getToolTipText() {
        return "Evaluate Imputation Accuracy";
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    public TableReport runPlugin(DataSet input) {
        return (TableReport) performFunction(input).getData(0).getData();
    }

}
