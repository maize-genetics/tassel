package net.maizegenetics.analysis.imputation;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.*;

public class ImputeProgenyStatesPlugin extends AbstractPlugin {

	private PluginParameter<Boolean> rephase = new PluginParameter.Builder<>("rephase", false, Boolean.class)
			.guiName("Rephase Parents First")
			.description("If true, rephase parents before imputing. If false, input haplotypes from file.")
			.build();
	private PluginParameter<String> input = PluginParameter.getLabelInstance("Input --------------------");
	private PluginParameter<String> parentageFile = new PluginParameter.Builder<>("parentage", null, String.class)
			.guiName("Parentage Input File")
			.inFile()
			.description("The input file containing the parentage which lists the parents of each progeny and whether they were derived by self or outcross.")
			.required(true)
			.build();
	private PluginParameter<String> parentHaplotypeFilename = new PluginParameter.Builder<>("parentHap", null, String.class)
			.guiName("Parent Haplotypes Input File")
			.inFile()
			.description("The input file containing the parent haplotypes expressed as nucleotides. ")
			.build();
	private PluginParameter<String> progenyFile = new PluginParameter.Builder<>("progeny", null, String.class)
			.dependentOnParameter(rephase)
			.guiName("Progeny States Input File")
			.inFile()
			.description("The input file containing the progeny states (parentcalls). Needed for rephasing using haplotype probabilities or for writing breakpoints.")
			.build();
	private PluginParameter<String> output = PluginParameter.getLabelInstance("Output --------------------");
	private PluginParameter<String> imputedFile = new PluginParameter.Builder<>("imputedOut", null, String.class)
			.guiName("Imputed Genotypes Output File")
			.outFile()
			.description("The output file containing the imputed progeny genotypes in hapmap format.")
			.build();
	private PluginParameter<String> statesFile = new PluginParameter.Builder<>("statesOut", null, String.class)
			.guiName("Progeny States Output File")
			.outFile()
			.description("The output file containing the new progeny states (parentcalls) in hapmap format")
			.build();
	private PluginParameter<String> hapProbFile = new PluginParameter.Builder<>("probOut", null, String.class)
			.guiName("Updated Haplotype Probabilities Output File")
			.outFile()
			.description("The output file containing the new parent haplotype probabilities, binary format. A .bin extension will be appended if not present.")
			.build();
	private PluginParameter<Boolean> writeBreakpoints = new PluginParameter.Builder<>("writebp", false, Boolean.class)
			.guiName("Write Breakpoints")
			.description("Write breakpoint file from imputed genotypes and parentage. Progeny states will not be imputed.")
			.build();
	private PluginParameter<String> breakpointOutFile = new PluginParameter.Builder<>("bpOut", null, String.class)
			.guiName("Breakpoint Output File")
			.outFile()
			.description("The name of the output file for writing the breakpoints. If the file exists, it will be overwritten.")
			.build();

	public ImputeProgenyStatesPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
	}
	
	@Override
	public DataSet processData(DataSet input) {

		List<Datum> resultList = new ArrayList<>();
		GenotypeTable inputGenotype = (GenotypeTable) input.getDataOfType(GenotypeTable.class).get(0).getData();
		ImputeCrossProgeny icp = new ImputeCrossProgeny();
		icp.setParentage(parentageFile.value()); //input

		if (writeBreakpoints()) {
			icp.writeBreakpointFile(inputGenotype, breakpointOutFile());
			return null;
		}

		icp.setMyGenotype(inputGenotype); //input
		icp.setParentage(parentageFile.value()); //input
		icp.setHaplotypeMap(parentHaplotypeFilename.value());  //input
		icp.setImputedGenotypeOutFilename(imputedFile.value());  //output
		icp.setParentcallOutFilename(statesFile.value());  //output
		icp.setPhasedParentOutFilename(hapProbFile.value());  //output
		
		if (rephase.value()) {
			
			icp.setParentCallInputFilename(progenyFile.value());
			icp.improveImputedProgenyStates();
		} else  {
			icp.imputeAll();
		}
		
		return new DataSet(resultList, this);
	}

	@Override
	public ImageIcon getIcon() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getButtonName() {
		return "Progeny States";
	}

	@Override
	public String getToolTipText() {
		return "Impute progeny states from parent haplotypes.";
	}

	@Override
	public String pluginDescription() {
		return "Impute progeny states using parental haplotypes either represented as nucleotides or as the probability that "
				+ "a haplotype carries the major allele at a site. The plugin also provides a method for estimating haplotype "
				+ "probabilities from progeny states.";
	}

//	public static void main(String[] args) {
//		GeneratePluginCode.generate(ImputeProgenyStatesPlugin.class);
//	}

	/**
	 * If true, rephase parents before imputing. If false,
	 * input haplotypes from file.
	 *
	 * @return Rephase Parents First
	 */
	public Boolean rephase() {
		return rephase.value();
	}

	/**
	 * Set Rephase Parents First. If true, rephase parents
	 * before imputing. If false, input haplotypes from file.
	 *
	 * @param value Rephase Parents First
	 *
	 * @return this plugin
	 */
	public ImputeProgenyStatesPlugin rephase(Boolean value) {
		rephase = new PluginParameter<>(rephase, value);
		return this;
	}


	/**
	 * The input file containing the parentage which lists
	 * the parents of each progeny and whether they were derived
	 * by self or outcross.
	 *
	 * @return Parentage Input File
	 */
	public String parentageFile() {
		return parentageFile.value();
	}

	/**
	 * Set Parentage Input File. The input file containing
	 * the parentage which lists the parents of each progeny
	 * and whether they were derived by self or outcross.
	 *
	 * @param value Parentage Input File
	 *
	 * @return this plugin
	 */
	public ImputeProgenyStatesPlugin parentageFile(String value) {
		parentageFile = new PluginParameter<>(parentageFile, value);
		return this;
	}

	/**
	 * The input file containing the parent haplotypes expressed
	 * as nucleotides.
	 *
	 * @return Parent Haplotypes Input File
	 */
	public String parentHaplotypeFilename() {
		return parentHaplotypeFilename.value();
	}

	/**
	 * Set Parent Haplotypes Input File. The input file containing
	 * the parent haplotypes expressed as nucleotides.
	 *
	 * @param value Parent Haplotypes Input File
	 *
	 * @return this plugin
	 */
	public ImputeProgenyStatesPlugin parentHaplotypeFilename(String value) {
		parentHaplotypeFilename = new PluginParameter<>(parentHaplotypeFilename, value);
		return this;
	}

	/**
	 * The input file containing the progeny states (parentcalls).
	 * Needed for rephasing using haplotype probabilities
	 * or for writing breakpoints.
	 *
	 * @return Progeny States Input File
	 */
	public String progenyFile() {
		return progenyFile.value();
	}

	/**
	 * Set Progeny States Input File. The input file containing
	 * the progeny states (parentcalls). Needed for rephasing
	 * using haplotype probabilities or for writing breakpoints.
	 *
	 * @param value Progeny States Input File
	 *
	 * @return this plugin
	 */
	public ImputeProgenyStatesPlugin progenyFile(String value) {
		progenyFile = new PluginParameter<>(progenyFile, value);
		return this;
	}

	/**
	 * The output file containing the imputed progeny genotypes
	 * in hapmap format.
	 *
	 * @return Imputed Genotypes Output File
	 */
	public String imputedFile() {
		return imputedFile.value();
	}

	/**
	 * Set Imputed Genotypes Output File. The output file
	 * containing the imputed progeny genotypes in hapmap
	 * format.
	 *
	 * @param value Imputed Genotypes Output File
	 *
	 * @return this plugin
	 */
	public ImputeProgenyStatesPlugin imputedFile(String value) {
		imputedFile = new PluginParameter<>(imputedFile, value);
		return this;
	}

	/**
	 * The output file containing the new progeny states (parentcalls)
	 * in hapmap format
	 *
	 * @return Progeny States Output File
	 */
	public String statesFile() {
		return statesFile.value();
	}

	/**
	 * Set Progeny States Output File. The output file containing
	 * the new progeny states (parentcalls) in hapmap format
	 *
	 * @param value Progeny States Output File
	 *
	 * @return this plugin
	 */
	public ImputeProgenyStatesPlugin statesFile(String value) {
		statesFile = new PluginParameter<>(statesFile, value);
		return this;
	}

	/**
	 * The output file containing the new parent haplotype
	 * probabilities, binary format. A .bin extension will
	 * be appended if not present.
	 *
	 * @return Updated Haplotype Probabilities Output File
	 */
	public String hapProbFile() {
		return hapProbFile.value();
	}

	/**
	 * Set Updated Haplotype Probabilities Output File. The
	 * output file containing the new parent haplotype probabilities,
	 * binary format. A .bin extension will be appended if
	 * not present.
	 *
	 * @param value Updated Haplotype Probabilities Output File
	 *
	 * @return this plugin
	 */
	public ImputeProgenyStatesPlugin hapProbFile(String value) {
		hapProbFile = new PluginParameter<>(hapProbFile, value);
		return this;
	}

	/**
	 * Write breakpoint file from imputed genotypes and parentage.
	 *
	 * @return Write Breakpoints
	 */
	public Boolean writeBreakpoints() {
		return writeBreakpoints.value();
	}

	/**
	 * Set Write Breakpoints. Write breakpoint file from imputed
	 * genotypes and parentage.
	 *
	 * @param value Write Breakpoints
	 *
	 * @return this plugin
	 */
	public ImputeProgenyStatesPlugin writeBreakpoints(Boolean value) {
		writeBreakpoints = new PluginParameter<>(writeBreakpoints, value);
		return this;
	}

	/**
	 * The name of the output file for writing the breakpoints.
	 * If the file exists, it will be overwritten.
	 *
	 * @return Breakpoint Output File
	 */
	public String breakpointOutFile() {
		return breakpointOutFile.value();
	}

	/**
	 * Set Breakpoint Output File. The name of the output
	 * file for writing the breakpoints. If the file exists,
	 * it will be overwritten.
	 *
	 * @param value Breakpoint Output File
	 *
	 * @return this plugin
	 */
	public ImputeProgenyStatesPlugin breakpointOutFile(String value) {
		breakpointOutFile = new PluginParameter<>(breakpointOutFile, value);
		return this;
	}
}
