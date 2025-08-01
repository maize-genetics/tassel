/*
 * FilterSiteNamePlugin.java
 *
 * Created on November 2, 2011
 *
 */
package net.maizegenetics.analysis.filter;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.util.OpenBitSet;

import javax.swing.*;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import net.maizegenetics.gui.SelectFromAvailableSitesDialog;
import net.maizegenetics.gui.SiteNamesAvailableListModel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * @deprecated 
 * Please use FilterSiteBuilderPlugin instead
 * 
 * @author Terry Casstevens
 */
public class FilterSiteNamePlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(FilterSiteNamePlugin.class);
    private int[] mySitesToKeep = null;
    private int[] mySitesToCovariates = null;
    private int[] mySitesToFactors = null;
    private String[] mySiteNamesToKeep = null;
    private String[] mySiteNamesToRemove = null;
    private String[] mySiteNamesToCovariates = null;
    private String[] mySiteNamesToFactors = null;

    /**
     * Creates a new instance of FilterSiteNamePlugin
     */
    public FilterSiteNamePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {

            List inputData = input.getDataOfType(GenotypeTable.class);
            if (inputData.size() != 1) {
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection. Please select a single alignment.");
                } else {
                    myLogger.error("performFunction: Please input a single alignment.");
                }
                return null;
            }

            List<Datum> td = processDatum((Datum) inputData.get(0), isInteractive());
            if (td == null) {
                return null;
            }

            DataSet output = new DataSet(td, this);

            fireDataSetReturned(new PluginEvent(output, FilterSiteNamePlugin.class));

            return output;

        } finally {
            fireProgress(100);
        }
    }

    private List<Datum> processDatum(Datum inDatum, boolean isInteractive) {

        final GenotypeTable alignment = (GenotypeTable) inDatum.getData();
        List<Datum> resultList = new ArrayList<Datum>();

        if (isInteractive) {
            SiteNamesAvailableListModel listModel = new SiteNamesAvailableListModel(alignment.positions());
            SelectFromAvailableSitesDialog dialog = new SelectFromAvailableSitesDialog(getParentFrame(), "Site Name Filter", listModel);
            dialog.setLocationRelativeTo(getParentFrame());
            dialog.setVisible(true);
            if (dialog.isCanceled()) {
                return null;
            }
            
            for (int[] indices : dialog.listOfSelectedIndices()) {
            	mySitesToKeep = indices;
            	addResults(resultList, alignment, inDatum.getName());
            	mySitesToKeep = null;
            }
            
            for (int[] indices : dialog.listOfCovariateIndices()) {
            	mySitesToCovariates = indices;
            	addResults(resultList, alignment, inDatum.getName());
            	mySitesToCovariates = null;
            }
            
            for (int[] indices : dialog.listOfFactorIndices()) {
            	mySitesToFactors = indices;
            	addResults(resultList, alignment, inDatum.getName());
            	mySitesToFactors = null;
            }
            
            dialog.dispose();
        } else addResults(resultList, alignment, inDatum.getName());

        return resultList;

    }

    private void addResults(List<Datum> resultList, GenotypeTable alignment, String datasetName) {
        GenotypeTable result = null;

        if (((mySitesToKeep != null) && (mySitesToKeep.length != 0))) {
            result = FilterGenotypeTable.getInstance(alignment, mySitesToKeep);
        } else if (((mySiteNamesToKeep != null) && (mySiteNamesToKeep.length != 0))) {
            result = FilterGenotypeTable.getInstance(alignment, mySiteNamesToKeep);
        } else if (((mySiteNamesToRemove != null) && (mySiteNamesToRemove.length != 0))) {
            result = FilterGenotypeTable.getInstanceRemoveSiteNames(alignment, mySiteNamesToRemove);
        } 

        if (result != null) {
            String theName, theComment;
            theName = datasetName + "_" + result.numberOfSites() + "_Sites";
            theComment = "Subset of " + result.numberOfSites() + " from " + datasetName;
            resultList.add(new Datum(theName, result, theComment));
        }
        
        if (mySiteNamesToCovariates != null) {
        	mySitesToCovariates = siteNumbersFromNames(mySiteNamesToCovariates, alignment);
        }
        if (mySitesToCovariates != null) {
        	List<ATTRIBUTE_TYPE> typeList = new ArrayList<>();
        	typeList.add(ATTRIBUTE_TYPE.taxa);
        	for (int site : mySitesToCovariates) typeList.add(ATTRIBUTE_TYPE.covariate);
        	
        	List<PhenotypeAttribute> attrList = new ArrayList<>();
        	attrList.add(new TaxaAttribute(alignment.taxa()));
        	for (int site : mySitesToCovariates) attrList.add(convertSiteToCovariate(site, alignment));
        	
        	Phenotype myPhenotype = new PhenotypeBuilder().fromAttributeList(attrList, typeList).build().get(0);
        	String name = "Site_Covariates_" + datasetName;
        	StringBuilder commentBuilder = new StringBuilder("Sites as covariates\n");
        	commentBuilder.append("from ").append(datasetName).append("\n");
        	for (int siteNumber : mySitesToCovariates) commentBuilder.append(alignment.siteName(siteNumber)).append("\n");
        	resultList.add(new Datum(name, myPhenotype, commentBuilder.toString() ));
        } 

        if (mySiteNamesToFactors != null) {
        	mySitesToFactors = siteNumbersFromNames(mySiteNamesToFactors, alignment);
        }
        if (mySitesToFactors != null) {
        	List<ATTRIBUTE_TYPE> typeList = new ArrayList<>();
        	typeList.add(ATTRIBUTE_TYPE.taxa);
        	for (int site : mySitesToFactors) typeList.add(ATTRIBUTE_TYPE.factor);
        	
        	List<PhenotypeAttribute> attrList = new ArrayList<>();
        	attrList.add(new TaxaAttribute(alignment.taxa()));
        	for (int site : mySitesToFactors) attrList.add(convertSiteToFactor(site, alignment));
        	
        	Phenotype myPhenotype = new PhenotypeBuilder().fromAttributeList(attrList, typeList).build().get(0);
        	String name = "Site_Factors_" + datasetName;
        	StringBuilder commentBuilder = new StringBuilder("Sites as factors\n");
        	commentBuilder.append("from ").append(datasetName).append("\n");
        	for (int siteNumber : mySitesToFactors) commentBuilder.append(alignment.siteName(siteNumber)).append("\n");
        	resultList.add(new Datum(name, myPhenotype, commentBuilder.toString() ));
        } 
    }
    
    private PhenotypeAttribute convertSiteToCovariate(int site, GenotypeTable geno) {
    	final byte major = geno.majorAllele(site);
    	final float sitemean = (float) (geno.majorAlleleFrequency(site) * 2);
    	final int ntaxa = geno.numberOfTaxa();
    	
    	float[] siteScore = new float[ntaxa];
    	for (int t = 0; t < ntaxa; t++) {
    		byte taxonGeno = geno.genotype(t, site);
    		if (taxonGeno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) siteScore[t] = sitemean;
    		else {
    			siteScore[t] = 0;
    			byte[] tgArray = GenotypeTableUtils.getDiploidValues(taxonGeno);
    			if (tgArray[0] == major) siteScore[t]++;
    			if (tgArray[1] == major) siteScore[t]++;
    		}
    	}
    	
    	return new NumericAttribute(geno.siteName(site), siteScore, new OpenBitSet(ntaxa));
    }
    
    private PhenotypeAttribute convertSiteToFactor(int site, GenotypeTable geno) {
    	final int ntaxa = geno.numberOfTaxa();
    	final Byte unknownByte = new Byte(GenotypeTable.UNKNOWN_DIPLOID_ALLELE);
    	byte[] byteGeno = geno.genotypeAllTaxa(site);
    	Map<Byte, Long> genoMap = IntStream.range(0, ntaxa).mapToObj(i -> new Byte(byteGeno[i]))
    			.filter(g -> !g.equals(unknownByte))
    			.collect(Collectors.groupingBy(g -> g, Collectors.counting()));
    	double totalGenoCount = genoMap.values().stream().mapToInt(v -> v.intValue()).sum();
    	double[] freq = genoMap.values().stream().mapToDouble(v -> v.doubleValue() / totalGenoCount).toArray();
    	double[] cumfreq = new double[freq.length];
    	cumfreq[0] = freq[0];
    	for (int i = 1; i < freq.length; i++) cumfreq[i] = cumfreq[i - 1] + freq[i];
    	
    	List<Byte> byteList = new ArrayList<Byte>(genoMap.keySet());
    	
    	Random ran = new Random();
    	String[] strGeno = IntStream.range(0, ntaxa)
    			.mapToObj(t -> {
    				if (byteGeno[t] == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
    					double p = ran.nextDouble();
    					int bin = 0;
    					while (p > cumfreq[bin]) bin++;
    					return NucleotideAlignmentConstants.getNucleotideIUPAC(byteList.get(bin));
    				}
    				return NucleotideAlignmentConstants.getNucleotideIUPAC(byteGeno[t]);
    			})
    			.toArray(String[]::new);
    	
    	return new CategoricalAttribute(geno.siteName(site), strGeno);
    }
    
    private int[] siteNumbersFromNames(String[] siteNames, GenotypeTable a) {
        Arrays.sort(siteNames);
        int[] temp = new int[siteNames.length];
        int count = 0;
        for (int i = 0, n = a.numberOfSites(); i < n; i++) {
            if (Arrays.binarySearch(siteNames, a.siteName(i)) >= 0) {
                temp[count++] = i;
                if (count == siteNames.length) {
                    break;
                }
            }
        }

        int[] result = null;
        if (count == siteNames.length) {
            result = temp;
        } else {
            result = new int[count];
            System.arraycopy(temp, 0, result, 0, count);
        }
        
        return result;
    }
    
    public int[] getSitesToKeep() {
        return mySitesToKeep;
    }

    public void setSitesToKeep(int[] sitesToKeep) {
        mySitesToKeep = sitesToKeep;
        validItemsSet();
    }

    public String[] getSiteNamesToKeep() {
        return mySiteNamesToKeep;
    }

    public void setSiteNamesToKeep(String[] sitesToKeep) {
        mySiteNamesToKeep = sitesToKeep;
        validItemsSet();
    }

    public String[] getSiteNamesToRemove() {
        return mySiteNamesToRemove;
    }

    public void setSiteNamesToRemove(String[] sitesToRemove) {
        mySiteNamesToRemove = sitesToRemove;
        validItemsSet();
    }

    public String[] getSiteNamesToCovariates() {
		return mySiteNamesToCovariates;
	}

	public void setSiteNamesToCovariates(String[] mySitesToCovariates) {
		mySiteNamesToCovariates = mySitesToCovariates;
	}

	public int[] getSitesToCovariatesIndex() {
		return mySitesToCovariates;
	}
	
	public void setSitesToCovariatesIndex(int[] index) {
		mySitesToCovariates = index;
	}
	
	public String[] getSiteNamesToFactors() {
		return mySiteNamesToFactors;
	}

	public void setSiteNamesToFactors(String[] mySitesToFactors) {
		this.mySiteNamesToFactors = mySitesToFactors;
	}

	public int[] getSitesToFactorsIndex() {
		return mySitesToFactors;
	}
	
	public void setSitesToFactorsIndex(int[] index) {
		mySitesToFactors = index;
	}
	
	private void validItemsSet() {

        int count = 0;
        if ((mySitesToKeep != null) && (mySitesToKeep.length != 0)) {
            count++;
        }
        if ((mySiteNamesToKeep != null) && (mySiteNamesToKeep.length != 0)) {
            count++;
        }
        if ((mySiteNamesToRemove != null) && (mySiteNamesToRemove.length != 0)) {
            count++;
        }

        if (count > 1) {
            throw new IllegalStateException("FilterSiteNamePlugin: validItemsSet: Can only set one of the following: sites to keep, site names to keep, or site names to remove.");
        }

    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = FilterSiteNamePlugin.class.getResource("/net/maizegenetics/analysis/images/Filter.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Site Names";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Select Site Names Within Dataset";
    }
}
