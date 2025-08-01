package net.maizegenetics.analysis.imputation;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ImageIcon;

import com.google.common.collect.Range;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.TableReport;

public class FSFHapImputationPlugin extends AbstractPlugin {
	
	//parameters for CallParentAllelesPlugin
	private PluginParameter<String> pedigreeFilename = new PluginParameter.Builder<>("pedigrees", null, String.class)
			.description("the pedigree file name")
                        .guiName("Pedigree File")
			.inFile().required(true).build();
	private PluginParameter<String> logFilename = new PluginParameter.Builder<>("logfile", null, String.class)
			.description("the name of a log file for runtime messages")
                        .guiName("Log File")
			.outFile().build();
	private PluginParameter<Boolean> useClusterAlgorithm = new PluginParameter.Builder<>("cluster", false, Boolean.class)
                        .guiName("Use Cluster Algorithm")
			.description("use the cluster algorithm").build();
	private PluginParameter<Boolean> useWindowLD = new PluginParameter.Builder<>("windowLD", false, Boolean.class)
                        .guiName("Use Window LD Algorithm")
			.description("use the windowLD algorithm").build();
	private PluginParameter<Boolean> useBCFilter = new PluginParameter.Builder<>("bc", true, Boolean.class)
                        .guiName("Use Single Backcross Algorithm")
			.description("use the single backcross algorithm").build();
	private PluginParameter<Boolean> useMultipleBCFilter = new PluginParameter.Builder<>("multbc", false, Boolean.class)
                        .guiName("Use Multiple Backcross Algorithm")
			.description("use the multiple backcross algorithm").build();
	private PluginParameter<Double> minMinorAlleleFreq = new PluginParameter.Builder<>("minMaf", 0.1, Double.class)
                        .guiName("Min Minor Allele Frequency")
			.range(Range.closed(0.0, 1.0)).description("filter out sites with less than minimumMinorAlleleFrequency").build();
	private PluginParameter<Integer> windowSize = new PluginParameter.Builder<>("window", 50, Integer.class)
                        .guiName("Window Size")
			.description("Window Size").build();
	private PluginParameter<Double> minRforSnps = new PluginParameter.Builder<>("minR", 0.2, Double.class)
			.range(Range.closed(0.0, 1.0)).description("filter out sites not correlated with neighboring sites").build();
	private PluginParameter<Double> maxMissing = new PluginParameter.Builder<>("maxMissing", 0.8, Double.class)
			.range(Range.closed(0.0, 1.0)).description("filter out sites with proportion missing > maxMissing").build();
	private PluginParameter<Boolean> noHets = new PluginParameter.Builder<>("nohets", false, Boolean.class)
                        .guiName("Don't Use Heterozygous Calls")
			.description("delete heterozygous calls before imputing").build();
	private PluginParameter<Integer> maxDifference = new PluginParameter.Builder<>("maxDiff", 0, Integer.class)
			.description("use to decide if two haplotypes are equivalent").build();
	private PluginParameter<Integer> minHaplotypeCluster = new PluginParameter.Builder<>("minHap", 5, Integer.class)
			.description("haplotype must be observed at least this often").build();
	private PluginParameter<Integer> overlap = new PluginParameter.Builder<>("overlap", 25, Integer.class)
                        .guiName("Window Overlap")
			.description("overlap between adjacent windows").build();
	
	//parameters for ViterbiAlgorithmPlugin
	private PluginParameter<Boolean> fillgaps = new PluginParameter.Builder<>("fillgaps", false, Boolean.class)
                        .guiName("Fill Gaps")
			.description("replace missing values with flanking values if equal").build();
	private PluginParameter<Double> probHeterozygous = new PluginParameter.Builder<>("phet", 0.07, Double.class)
                        .guiName("Proportion Heterozygous")
			.range(Range.closed(0.0, 1.0)).description("proportion of sites that are heterozygous").build();
	
	//parameters for WritePopulationAlignmentPlugin
	private PluginParameter<Boolean> mergeAlignments = new PluginParameter.Builder<>("merge", false, Boolean.class)
			.description("merge families and chromosomes").build();
	private PluginParameter<Boolean> outParentCalls = new PluginParameter.Builder<>("outParents", true, Boolean.class)
			.description("replace missing values with flanking values if equal").build();
	private PluginParameter<Boolean> outNucleotides = new PluginParameter.Builder<>("outNuc", true, Boolean.class)
                        .guiName("Out Nucleotides")
			.description("replace missing values with flanking values if equal").build();
	private PluginParameter<Boolean> outIUPAC = new PluginParameter.Builder<>("outIUPAC", true, Boolean.class)
                        .guiName("Output IUPAC Codes")
			.description("use IUPAC ambiguity codes for output").build();
	
	
	public FSFHapImputationPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
	}
	
	@Override
	public String getCitation() {
        return "Swarts K, Li H, Romero Navarro JA, Romay-Alvarez MC, Hearne S, Acharya C, "
                + "Glaubitz JC, Mitchell S, Elshire RJ, Buckler ES, Bradbury PJ (2014) "
                + "FSFHap (Full-Sib Family Haplotype Imputation) and FILLIN "
                + "(Fast, Inbred Line Library ImputatioN) optimize genotypic imputation "
                + "for low-coverage, next-generation sequence data in crop plants. "
                + "Plant Genome (in review)";
	}
	
	@Override
	public String pluginDescription() {
		return "The FSFHapImputation Plugin infers parental haplotypes for a full sib family then uses those haplotypes in an HMM to impute variants. "
				+ "It is effective at correctly imputing heterzygotes in GBS data. To use from the command line, use TASSEL's default syntax that " +
                "passes data from one plugin to another (Note that this creates 2 files, one of just parental calls (A/C) and one of imputed genotypes):\n\n" +
                "\trun_pipeline.pl -h input.hmp.txt -FSFHapImputationPlugin [options] -endPLugin -export output.hmp.txt";
	}

	public DataSet processData(DataSet input) {
		
		try {
			CallParentAllelesPlugin cpa = new CallParentAllelesPlugin(null);
			cpa.setPedfileName(pedigreeFilename.value());
			cpa.setLogFile(logFilename.value());
			cpa.setUseClusterAlgorithm(useClusterAlgorithm.value());
			cpa.setUseWindowLD(useWindowLD.value());
			cpa.setUseBCFilter(useBCFilter.value());
			cpa.setUseMultipleBCFilter(useMultipleBCFilter.value());
			cpa.setMinMinorAlleleFrequency(minMinorAlleleFreq.value());
			cpa.setWindowSize(windowSize.value());
			cpa.setMinRforSnps(minRforSnps.value());
			cpa.setMaxMissing(maxMissing.value());
			cpa.setUseHets(!noHets.value());
			cpa.setMaxDifference(maxDifference.value());
			cpa.setMinUsedClusterSize(minHaplotypeCluster.value());
			cpa.setOverlap(overlap.value());
			fireProgress(10);
			DataSet cpaResult = cpa.performFunction(input);
			
			fireProgress(30);
			
			ViterbiAlgorithmPlugin vap = new ViterbiAlgorithmPlugin(null);
			vap.setFillGapsInAlignment(fillgaps.value());
			vap.setProbHeterozygous(probHeterozygous.value());
			DataSet vapResult = vap.performFunction(cpaResult);
			
			fireProgress(60);
			
			WritePopulationAlignmentPlugin writePap = new WritePopulationAlignmentPlugin(null);
			writePap.setMergeAlignments(mergeAlignments.value());
			writePap.setWriteParentCalls(outParentCalls.value());
			writePap.setWriteNucleotides(outNucleotides.value());
			writePap.setOutputDiploid(!outIUPAC.value());
			DataSet writeResult = writePap.performFunction(vapResult);
			
			List<Datum> allData = new ArrayList<>(writeResult.getDataSet());
			allData.addAll(cpaResult.getDataOfType(TableReport.class));
			
			fireProgress(90); 
//			return writeResult;
			return new DataSet(allData, this);
		} finally {
			fireProgress(100);
		}
	}

	@Override
	public ImageIcon getIcon() {
		return null;
	}

	@Override
	public String getButtonName() {
		return "Impute By FSFHap";
	}

	@Override
	public String getToolTipText() {
		return "Impute variants in full sib families";
	}

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
//     public static void main(String[] args) {
//         GeneratePluginCode.generate(FSFHapImputationPlugin.class);
//     }

     /**
      * Convenience method to run plugin with one return object.
      */
     // TODO: Replace <Type> with specific type.
     public GenotypeTable runPlugin(DataSet input) {
         return (GenotypeTable) performFunction(input).getData(0).getData();
     }

     /**
      * the pedigree file name
      *
      * @return Pedigrees
      */
     public String pedigrees() {
         return pedigreeFilename.value();
     }

     /**
      * Set Pedigrees. the pedigree file name
      *
      * @param value Pedigrees
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin pedigrees(String value) {
         pedigreeFilename = new PluginParameter<>(pedigreeFilename, value);
         return this;
     }

     /**
      * use the cluster algorithm
      *
      * @return Cluster
      */
     public Boolean cluster() {
         return useClusterAlgorithm.value();
     }

     /**
      * Set Cluster. use the cluster algorithm
      *
      * @param value Cluster
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin cluster(Boolean value) {
         useClusterAlgorithm = new PluginParameter<>(useClusterAlgorithm, value);
         return this;
     }

     /**
      * use the windowLD algorithm
      *
      * @return Window L D
      */
     public Boolean windowLD() {
         return useWindowLD.value();
     }

     /**
      * Set Window L D. use the windowLD algorithm
      *
      * @param value Window L D
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin windowLD(Boolean value) {
         useWindowLD = new PluginParameter<>(useWindowLD, value);
         return this;
     }

     /**
      * use the single backcross algorithm
      *
      * @return Bc
      */
     public Boolean bc() {
         return useBCFilter.value();
     }

     /**
      * Set Bc. use the single backcross algorithm
      *
      * @param value Bc
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin bc(Boolean value) {
         useBCFilter = new PluginParameter<>(useBCFilter, value);
         return this;
     }

     /**
      * use the multiple backcross algorithm
      *
      * @return Multbc
      */
     public Boolean multbc() {
         return useMultipleBCFilter.value();
     }

     /**
      * Set Multbc. use the multiple backcross algorithm
      *
      * @param value Multbc
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin multbc(Boolean value) {
         useMultipleBCFilter = new PluginParameter<>(useMultipleBCFilter, value);
         return this;
     }

     /**
      * filter out sites with less than minimumMinorAlleleFrequency
      *
      * @return Min Maf
      */
     public Double minMaf() {
         return minMinorAlleleFreq.value();
     }

     /**
      * Set Min Maf. filter out sites with less than minimumMinorAlleleFrequency
      *
      * @param value Min Maf
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin minMaf(Double value) {
         minMinorAlleleFreq = new PluginParameter<>(minMinorAlleleFreq, value);
         return this;
     }

     /**
      * filter out sites with less than minimumMinorAlleleFrequency
      *
      * @return Window
      */
     public Integer window() {
         return windowSize.value();
     }

     /**
      * Set Window. filter out sites with less than minimumMinorAlleleFrequency
      *
      * @param value Window
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin window(Integer value) {
         windowSize = new PluginParameter<>(windowSize, value);
         return this;
     }

     /**
      * filter out sites not correlated with neighboring sites
      *
      * @return Min R
      */
     public Double minR() {
         return minRforSnps.value();
     }

     /**
      * Set Min R. filter out sites not correlated with neighboring
      * sites
      *
      * @param value Min R
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin minR(Double value) {
         minRforSnps = new PluginParameter<>(minRforSnps, value);
         return this;
     }

     /**
      * filter out sites with proportion missing > maxMissing
      *
      * @return Max Missing
      */
     public Double maxMissing() {
         return maxMissing.value();
     }

     /**
      * Set Max Missing. filter out sites with proportion missing
      * > maxMissing
      *
      * @param value Max Missing
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin maxMissing(Double value) {
         maxMissing = new PluginParameter<>(maxMissing, value);
         return this;
     }

     /**
      * delete heterozygous calls before imputing
      *
      * @return Nohets
      */
     public Boolean nohets() {
         return noHets.value();
     }

     /**
      * Set Nohets. delete heterozygous calls before imputing
      *
      * @param value Nohets
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin nohets(Boolean value) {
         noHets = new PluginParameter<>(noHets, value);
         return this;
     }

     /**
      * use to decide if two haplotypes are equivalent
      *
      * @return Max Diff
      */
     public Integer maxDiff() {
         return maxDifference.value();
     }

     /**
      * Set Max Diff. use to decide if two haplotypes are equivalent
      *
      * @param value Max Diff
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin maxDiff(Integer value) {
         maxDifference = new PluginParameter<>(maxDifference, value);
         return this;
     }

     /**
      * haplotype must be observed at least this often
      *
      * @return Min Hap
      */
     public Integer minHap() {
         return minHaplotypeCluster.value();
     }

     /**
      * Set Min Hap. haplotype must be observed at least this
      * often
      *
      * @param value Min Hap
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin minHap(Integer value) {
         minHaplotypeCluster = new PluginParameter<>(minHaplotypeCluster, value);
         return this;
     }

     /**
      * overlap between adjacent windows
      *
      * @return Overlap
      */
     public Integer overlap() {
         return overlap.value();
     }

     /**
      * Set Overlap. overlap between adjacent windows
      *
      * @param value Overlap
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin overlap(Integer value) {
         overlap = new PluginParameter<>(overlap, value);
         return this;
     }

     /**
      * replace missing values with flanking values if equal
      *
      * @return Fillgaps
      */
     public Boolean fillgaps() {
         return fillgaps.value();
     }

     /**
      * Set Fillgaps. replace missing values with flanking
      * values if equal
      *
      * @param value Fillgaps
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin fillgaps(Boolean value) {
         fillgaps = new PluginParameter<>(fillgaps, value);
         return this;
     }

     /**
      * proportion of sites that are heterozygous
      *
      * @return Phet
      */
     public Double phet() {
         return probHeterozygous.value();
     }

     /**
      * Set Phet. proportion of sites that are heterozygous
      *
      * @param value Phet
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin phet(Double value) {
         probHeterozygous = new PluginParameter<>(probHeterozygous, value);
         return this;
     }

     /**
      * merge families and chromosomes
      *
      * @return Merge
      */
     public Boolean merge() {
         return mergeAlignments.value();
     }

     /**
      * Set Merge. merge families and chromosomes
      *
      * @param value Merge
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin merge(Boolean value) {
         mergeAlignments = new PluginParameter<>(mergeAlignments, value);
         return this;
     }

     /**
      * replace missing values with flanking values if equal
      *
      * @return Out Parents
      */
     public Boolean outParents() {
         return outParentCalls.value();
     }

     /**
      * Set Out Parents. replace missing values with flanking
      * values if equal
      *
      * @param value Out Parents
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin outParents(Boolean value) {
         outParentCalls = new PluginParameter<>(outParentCalls, value);
         return this;
     }

     /**
      * replace missing values with flanking values if equal
      *
      * @return Out Nuc
      */
     public Boolean outNuc() {
         return outNucleotides.value();
     }

     /**
      * Set Out Nuc. replace missing values with flanking values
      * if equal
      *
      * @param value Out Nuc
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin outNuc(Boolean value) {
         outNucleotides = new PluginParameter<>(outNucleotides, value);
         return this;
     }

     /**
      * use IUPAC ambiguity codes for output
      *
      * @return Out I U P A C
      */
     public Boolean outIUPAC() {
         return outIUPAC.value();
     }

     /**
      * Set Out I U P A C. use IUPAC ambiguity codes for output
      *
      * @param value Out I U P A C
      *
      * @return this plugin
      */
     public FSFHapImputationPlugin outIUPAC(Boolean value) {
         outIUPAC = new PluginParameter<>(outIUPAC, value);
         return this;
     }


}
