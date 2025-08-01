package net.maizegenetics.analysis.association;

import java.awt.Frame;
import java.io.File;
import java.net.URL;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable.GENOTYPE_TABLE_COMPONENT;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.TableReport;

import net.maizegenetics.util.TableReportUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.google.common.collect.Range;


public class FixedEffectLMPlugin extends AbstractPlugin {
	
    private static final Logger myLogger = LogManager.getLogger(FixedEffectLMPlugin.class);
    String baseOutFileName = "";
    
	enum GENOTYPE_DATA_TYPE { genotype, probability, allele_probabilities, none };
	private GenotypeTable.GENOTYPE_TABLE_COMPONENT[] GENOTYPE_COMP = new GenotypeTable.GENOTYPE_TABLE_COMPONENT[]{
	        GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability, GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability};
	
    //parameters
	private PluginParameter<Boolean> phenotypeOnly = new PluginParameter.Builder<>("phenoOnly", false, Boolean.class)
			.description("Should the phenotype be analyzed with no markers and BLUEs generated? (BLUE = best linear unbiased estimate)")
			.guiName("Analyze Phenotype Only")
			.build();
    private PluginParameter<Boolean> saveAsFile = new PluginParameter.Builder<>("saveToFile", false, Boolean.class)
    		.description("Should the results be saved to a file rather than stored in memory? It true, the results will be written to a file as each SNP is analyzed in order to reduce memory requirements"
    				+ "and the results will NOT be saved to the data tree. Default = false.")
    		.guiName("Save to file")
    		.build();
    private PluginParameter<String> siteReportFilename = new PluginParameter.Builder<>("siteFile", null, String.class)
    		.outFile()
    		.dependentOnParameter(saveAsFile)
    		.description("The name of the file to which these results will be saved.")
    		.guiName("Statistics File")
    		.build();
    private PluginParameter<String> alleleReportFilename = new PluginParameter.Builder<>("alleleFile", null, String.class)
    		.outFile()
    		.dependentOnParameter(saveAsFile)
    		.description("The name of the file to which these results will be saved.")
    		.guiName("Genotype Effect File")
    		.build();
    private PluginParameter<String> bluesReportFilename = new PluginParameter.Builder<>("bluesFile", null, String.class)
            .description("Name of file to which BLUEs values will be saved")
            .outFile()
            .dependentOnParameter(phenotypeOnly)
            .guiName("BLUEs File")
            .build();
    private PluginParameter<String> anovaReportFilename = new PluginParameter.Builder<>("anovaFile", null, String.class)
            .description("Name of file to which ANOVA report will be saved")
            .outFile()
            .dependentOnParameter(phenotypeOnly)
            .guiName("ANOVA File")
            .build();
    private PluginParameter<Double> maxPvalue = new PluginParameter.Builder<>("maxP", 1.0, Double.class)
    		.description("Only results with p <= maxPvalue will be reported. Default = 1.0.")
    		.dependentOnParameter(phenotypeOnly, false)
    		.range(Range.closed(0.0, 1.0))
    		.guiName("max P value")
    		.build();
    private PluginParameter<Boolean> permute = new PluginParameter.Builder<>("permute", false, Boolean.class)
    		.description("Should a permutation analysis be run? The permutation analysis controls the experiment-wise error rate for individual phenotypes.")
    		.dependentOnParameter(phenotypeOnly, false)
    		.guiName("Run Permutations")
    		.build();
    private PluginParameter<Integer> numberOfPermutations = new PluginParameter.Builder<>("nperm", 0, Integer.class)
    		.description("The number of permutations to be run for the permutation analysis.")
    		.dependentOnParameter(permute)
    		.guiName("Number of Permutations")
    		.build();
	private PluginParameter<GenotypeTable.GENOTYPE_TABLE_COMPONENT> myGenotypeTable = new PluginParameter.Builder<>("genotypeComponent", GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT.class)
			.genotypeTable()
	        .range(GENOTYPE_COMP)
	        .description("If the genotype table contains more than one type of genotype data, choose the type to use for the analysis.")
	        .build();
	
	private PluginParameter<Integer> minClassSize = new PluginParameter.Builder<>("minClassSize", 0, Integer.class)
			.description("The minimum acceptable genotype class size. Genotypes in a class with a smaller size will be set to missing.")
			.guiName("Minimum Class Size")
			.build();
	private PluginParameter<Boolean> biallelicOnly = new PluginParameter.Builder<>("biallelicOnly", false, Boolean.class)
			.description("Only test sites that are bi-allelic. The alternative is to test sites with two or more alleles.")
			.guiName("Bi-Allelic Sites Only")
			.build();
    private PluginParameter<Boolean> siteStatsOutput = new PluginParameter.Builder<>("siteStatsOut", false, Boolean.class)
    		.description("Generate an output dataset with only p-val, F statistic, and number of obs per site for all sites.")
    		.guiName("Output Site Stats")
    		.build();
    private PluginParameter<String> siteStatFilename = new PluginParameter.Builder<>("siteStatFile", null, String.class)
    		.description("")
    		.guiName("")
    		.dependentOnParameter(siteStatsOutput)
    		.outFile()
    		.build();
    private PluginParameter<Boolean> appendAddDom = new PluginParameter.Builder<>("appendAddDom", false, Boolean.class)
    		.description("If true, additive and dominance effect estimates will be added to the stats report for bi-allelic sites only. The effect will only be estimated when the data source is genotype (not a probability). The additive effect will always be non-negative.")
    		.guiName("Append Effect Estimates to Stats")
//    		.dependentOnParameter(myGenotypeTable, GENOTYPE_COMP[0])
    		.build();
	
    public FixedEffectLMPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public String getButtonName() {
        return "GLM";
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = FixedEffectLMPlugin.class.getResource("/net/maizegenetics/analysis/images/LinearAssociation.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getToolTipText() {
        return "Use fixed effect model to test associations";
    }

    protected void preProcessParameters(DataSet data) {
    	List<Datum> genoPhenoList = data.getDataOfType(GenotypePhenotype.class);
    	
    	if (genoPhenoList.size() == 0){
    		List<Datum> phenoList = data.getDataOfType(Phenotype.class);
    		if (phenoList.size() == 0) throw new IllegalArgumentException("A dataset that can be analyzed by GLM has not been selected.");
    		else if (phenoList.size() == 1) {
        		phenotypeOnly = new PluginParameter.Builder<>("phenoOnly", true, Boolean.class)
            			.description("Should the phenotype be analyzed with no markers and BLUEs generated? (BLUE = best linear unbiased estimate)")
            			.guiName("Analyze Phenotype Only")
            			.build();
    		} else  throw new IllegalArgumentException("GLM can only process one data set at a time.");
    	} 
    	else if (genoPhenoList.size() > 1)  throw new IllegalArgumentException("GLM can only process one data set at a time.");
    	else {
        	//code to handle Tassel 4 pipeline style commands
    		
    	}
    	
    }
    
    public DataSet processData(DataSet data) {
    	if (phenotypeOnly.value()) {
    	    List<Datum> datumList = data.getDataOfType(Phenotype.class);
    	    if (datumList.size() == 1) {
    	        Datum myDatum = datumList.get(0);
    	        PhenotypeLM plm = new PhenotypeLM(myDatum);
    	        if (saveAsFile()) {
                    TableReportUtils.saveDelimitedTableReport(plm.myBlues, new File(bluesReportFilename()));
                    TableReportUtils.saveDelimitedTableReport(plm.report(), new File(anovaReportFilename()));
                    return null;
                }
    	        return new DataSet(plm.datumList(), this);
    	    } else if (datumList.size() == 0) {
    	        throw new IllegalArgumentException("The phenotype only option was selected, but no phenotype data set was provided as input.");
    	    } else {
    	        throw new IllegalArgumentException("Multiple phenotype data sets were provided. Only one is allowed.");
    	    }

    	} else {
        	FixedEffectLM myLM; 
    		Datum myDatum = data.getDataOfType(GenotypePhenotype.class).get(0);
        	if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype) {
        		myLM = new DiscreteSitesFELM(myDatum, this);
        	} else if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability) {
        		myLM = new ReferenceProbabilityFELM(myDatum, this);
        	} else if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability) {
        		myLM = new AlleleProbabilityFELM(myDatum, this);
        	} else return null;
        	if (permute.value()) myLM.permutationTest(true, numberOfPermutations.value());
        	if (saveAsFile.value()) {
        		myLM.siteReportFilepath(siteReportFilename.value());
        		myLM.alleleReportFilepath(alleleReportFilename.value());
        	}
        	myLM.maxP(maxPvalue.value());
        	myLM.biallelicOnly(biallelicOnly.value());
        	myLM.minimumClassSize(minClassSize.value());
        	myLM.appendAddDom(appendAddDom.value());
        	myLM.solve();
        	if (saveAsFile.value()) return null;
        	else return new DataSet(myLM.datumList(), this);
    	} 
    	
    }
    
    public void updateProgress(int percent) {
    	fireProgress(percent);
    }
    
    //setters needed for compatability with Tassel 4.0 pipeline commands
    public void setOutputFile(String name) {
    	baseOutFileName = name;
    }
    
    public void setMaxP(double maxp) {
    	maxPvalue(maxp);
    }
    
    public void setPermute(boolean permute) {
    	permute(permute);
    }
    
    public void setNumberOfPermutations(int nperm) {
    	numberOfPermutations(nperm);
    }
    
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(FixedEffectLMPlugin.class);
    // }

    /**
     * Should the phenotype be analyzed with no markers and
     * BLUEs generated? (BLUE = best linear unbiased estimate)
     *
     * @return Analyze Phenotype Only
     */
    public Boolean phenotypeOnly() {
        return phenotypeOnly.value();
    }

    /**
     * Set Analyze Phenotype Only. Should the phenotype be
     * analyzed with no markers and BLUEs generated? (BLUE
     * = best linear unbiased estimate)
     *
     * @param value Analyze Phenotype Only
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin phenotypeOnly(Boolean value) {
        phenotypeOnly = new PluginParameter<>(phenotypeOnly, value);
        return this;
    }

    /**
     * Should the results be saved to a file rather than stored
     * in memory? It true, the results will be written to
     * a file as each SNP is analyzed in order to reduce memory
     * requirementsand the results will NOT be saved to the
     * data tree. Default = false.
     *
     * @return Save to file
     */
    public Boolean saveAsFile() {
        return saveAsFile.value();
    }

    /**
     * Set Save to file. Should the results be saved to a
     * file rather than stored in memory? It true, the results
     * will be written to a file as each SNP is analyzed in
     * order to reduce memory requirementsand the results
     * will NOT be saved to the data tree. Default = false.
     *
     * @param value Save to file
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin saveAsFile(Boolean value) {
        saveAsFile = new PluginParameter<>(saveAsFile, value);
        return this;
    }

    /**
     * The name of the file to which these results will be
     * saved.
     *
     * @return Statistics File
     */
    public String siteReportFilename() {
        return siteReportFilename.value();
    }

    /**
     * Set Statistics File. The name of the file to which
     * these results will be saved.
     *
     * @param value Statistics File
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin siteReportFilename(String value) {
        siteReportFilename = new PluginParameter<>(siteReportFilename, value);
        return this;
    }

    /**
     * The name of the file to which these results will be
     * saved.
     *
     * @return Genotype Effect File
     */
    public String alleleReportFilename() {
        return alleleReportFilename.value();
    }

    /**
     * Set Genotype Effect File. The name of the file to which
     * these results will be saved.
     *
     * @param value Genotype Effect File
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin alleleReportFilename(String value) {
        alleleReportFilename = new PluginParameter<>(alleleReportFilename, value);
        return this;
    }

    /**
     * Only results with p <= maxPvalue will be reported.
     * Default = 1.0.
     *
     * @return max P value
     */
    public Double maxPvalue() {
        return maxPvalue.value();
    }

    /**
     * Set max P value. Only results with p <= maxPvalue will
     * be reported. Default = 1.0.
     *
     * @param value max P value
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin maxPvalue(Double value) {
        maxPvalue = new PluginParameter<>(maxPvalue, value);
        return this;
    }

    /**
     * Should a permutation analysis be run? The permutation
     * analysis controls the experiment-wise error rate for
     * individual phenotypes.
     *
     * @return Run Permutations
     */
    public Boolean permute() {
        return permute.value();
    }

    /**
     * Set Run Permutations. Should a permutation analysis
     * be run? The permutation analysis controls the experiment-wise
     * error rate for individual phenotypes.
     *
     * @param value Run Permutations
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin permute(Boolean value) {
        permute = new PluginParameter<>(permute, value);
        return this;
    }

    /**
     * The number of permutations to be run for the permutation
     * analysis.
     *
     * @return Number of Permutations
     */
    public Integer numberOfPermutations() {
        return numberOfPermutations.value();
    }

    /**
     * Set Number of Permutations. The number of permutations
     * to be run for the permutation analysis.
     *
     * @param value Number of Permutations
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin numberOfPermutations(Integer value) {
        numberOfPermutations = new PluginParameter<>(numberOfPermutations, value);
        return this;
    }

    /**
     * If the genotype table contains more than one type of
     * genotype data, choose the type to use for the analysis.
     *
     * @return Genotype Component
     */
    public GENOTYPE_TABLE_COMPONENT genotypeTable() {
        return myGenotypeTable.value();
    }

    /**
     * Set Genotype Component. If the genotype table contains
     * more than one type of genotype data, choose the type
     * to use for the analysis.
     *
     * @param value Genotype Component
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin genotypeTable(GENOTYPE_TABLE_COMPONENT value) {
        myGenotypeTable = new PluginParameter<>(myGenotypeTable, value);
        return this;
    }

    /**
     * The minimum acceptable genotype class size. Genotypes
     * in a class with a smaller size will be set to missing.
     *
     * @return Minimum Class Size
     */
    public Integer minClassSize() {
        return minClassSize.value();
    }

    /**
     * Set Minimum Class Size. The minimum acceptable genotype
     * class size. Genotypes in a class with a smaller size
     * will be set to missing.
     *
     * @param value Minimum Class Size
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin minClassSize(Integer value) {
        minClassSize = new PluginParameter<>(minClassSize, value);
        return this;
    }

    /**
     * Only test sites that are bi-allelic. The alternative
     * is to test sites with two or more alleles.
     *
     * @return Bi-Allelic Sites Only
     */
    public Boolean biallelicOnly() {
        return biallelicOnly.value();
    }

    /**
     * Set Bi-Allelic Sites Only. Only test sites that are
     * bi-allelic. The alternative is to test sites with two
     * or more alleles.
     *
     * @param value Bi-Allelic Sites Only
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin biallelicOnly(Boolean value) {
        biallelicOnly = new PluginParameter<>(biallelicOnly, value);
        return this;
    }

    /**
     * Generate an output dataset with only p-val, F statistic,
     * and number of obs per site for all sites.
     *
     * @return Output Site Stats
     */
    public Boolean siteStatsOutput() {
        return siteStatsOutput.value();
    }

    /**
     * Set Output Site Stats. Generate an output dataset with
     * only p-val, F statistic, and number of obs per site
     * for all sites.
     *
     * @param value Output Site Stats
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin siteStatsOutput(Boolean value) {
        siteStatsOutput = new PluginParameter<>(siteStatsOutput, value);
        return this;
    }

    /**
     * Site Stat File
     *
     * @return Site Stat File
     */
    public String siteStatFilename() {
        return siteStatFilename.value();
    }

    /**
     * Set Site Stat File. Site Stat File
     *
     * @param value Site Stat File
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin siteStatFilename(String value) {
        siteStatFilename = new PluginParameter<>(siteStatFilename, value);
        return this;
    }

    /**
     * If true, additive and dominance effect estimates will
     * be added to the stats report for bi-allelic sites only.
     * The effect will only be estimated when the data source
     * is genotype (not a probability). The additive effect
     * will always be non-negative.
     *
     * @return Append Effect Estimates to Stats
     */
    public Boolean appendAddDom() {
        return appendAddDom.value();
    }

    /**
     * Set Append Effect Estimates to Stats. If true, additive
     * and dominance effect estimates will be added to the
     * stats report for bi-allelic sites only. The effect
     * will only be estimated when the data source is genotype
     * (not a probability). The additive effect will always
     * be non-negative.
     *
     * @param value Append Effect Estimates to Stats
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin appendAddDom(Boolean value) {
        appendAddDom = new PluginParameter<>(appendAddDom, value);
        return this;
    }

    /**
     * Name of file to which BLUEs values will be saved
     *
     * @return BLUEs File
     */
    public String bluesReportFilename() {
        return bluesReportFilename.value();
    }

    /**
     * Set BLUEs File. Name of file to which BLUEs values
     * will be saved
     *
     * @param value BLUEs File
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin bluesReportFilename(String value) {
        bluesReportFilename = new PluginParameter<>(bluesReportFilename, value);
        return this;
    }

    /**
     * Name of file to which ANOVA report will be saved
     *
     * @return ANOVA File
     */
    public String anovaReportFilename() {
        return anovaReportFilename.value();
    }

    /**
     * Set ANOVA File. Name of file to which ANOVA report
     * will be saved
     *
     * @param value ANOVA File
     *
     * @return this plugin
     */
    public FixedEffectLMPlugin anovaReportFilename(String value) {
        anovaReportFilename = new PluginParameter<>(anovaReportFilename, value);
        return this;
    }
}


