package net.maizegenetics.analysis.imputation;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.TreeSet;

import javax.swing.ImageIcon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.data.FileLoadPlugin.TasselFileType;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.SimpleTableReport;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportBuilder;

public class ParentPhasingPlugin extends AbstractPlugin {
	private static Logger myLogger = LogManager.getLogger(ParentPhasingPlugin.class);
	private static final byte N = GenotypeTable.UNKNOWN_ALLELE;
	
	private PluginParameter<Boolean> phaseParents = new PluginParameter.Builder<>("phase", false, Boolean.class)
			.guiName("Phase Parents")
			.description("Phase parents")
			.build();
	private PluginParameter<Boolean> rephaseParents = new PluginParameter.Builder<>("rephase", false, Boolean.class)
			.guiName("Re-Phase Parents")
			.description("Rephase parents")
			.build();
	private PluginParameter<String> parentageFile = new PluginParameter.Builder<>("parentage", null, String.class)
			.guiName("Parentage File")
			.inFile()
			.description("The file containing the parentage which lists the parents of each progeny and whether they were derived by self or outcross.")
			.build();

	//dependent on phaseParents
	private PluginParameter<Boolean> selfonly = new PluginParameter.Builder<>("self", false, Boolean.class)
			.dependentOnParameter(phaseParents)
			.guiName("Self Families Only")
			.description("If checked or true, phases using only families created by selfing a single parent. The alternative is to use all families. "
					+ "Compared to using all families, self only generally does a better job of phasing with a lot less missing data.")
			.build();
	private PluginParameter<Integer> windowSize = new PluginParameter.Builder<>("window", 50, Integer.class)
			.dependentOnParameter(selfonly)
			.guiName("Window Size")
			.description("")
			.build();
	private PluginParameter<Double> maxMissing = new PluginParameter.Builder<>("maxMissing", 0.7, Double.class)
			.dependentOnParameter(selfonly)
			.guiName("Max Proportion Missing Data")
			.description("Maximum allowable proportion of missing data for a site.")
			.build();
	
	private PluginParameter<String> outputFile = new PluginParameter.Builder<>("out", null, String.class)
			.dependentOnParameter(phaseParents)
			.guiName("Output File")
			.outFile()
			.description("The file contain the phased parent haplotypes in binary format. A .bin extension will be added to the filename.")
			.build();
	
	
	private PluginParameter<String> separator1 = PluginParameter.getLabelInstance("Files  ------------------------------");
	
	//dependent on rephase
	//need progeny parentcalls, parent haplotypes, and output filename
	private PluginParameter<String> parentCallFilename = new PluginParameter.Builder<>("parentcalls", null, String.class)
			.dependentOnParameter(rephaseParents)
			.guiName("Progeny States (parentcalls)")
			.inFile()
			.description("The genotype file containing the imputed progeny states, probably identified as parentcalls.")
			.build();
	private PluginParameter<String> parentHaplotypeFilename = new PluginParameter.Builder<>("parentHaplotypes", null, String.class)
			.dependentOnParameter(rephaseParents)
			.guiName("Parent Haplotype File")
			.inFile()
			.description("The file containing the haplotypes of the parents. Rephasing parents uses the progeny states to improve these haplotype calls.")
			.build();
	private PluginParameter<String> rephaseOutFile = new PluginParameter.Builder<>("rephaseOut", null, String.class)
			.dependentOnParameter(rephaseParents)
			.guiName("Rephase Output File")
			.outFile()
			.description("The file that will contain the rephased, improved parent haplotype calls.")
			.build();
	
	//dependent on combine
	private PluginParameter<String> separator2 = PluginParameter.getLabelInstance("Combine phased data ------------------------------");
	private PluginParameter<Boolean> combine = new PluginParameter.Builder<>("combine", false, Boolean.class)
			.guiName("Combine Phasing")
			.description("Combines two methods of phasing parents. Any sites that disagree are set to missing. Uses binary files as input. A report comparing the sites in common between the files will be generated.")
			.build();

	private PluginParameter<String> phased1 = new PluginParameter.Builder<>("phased1", null, String.class)
			.dependentOnParameter(combine)
			.guiName("First Phased File")
			.inFile()
			.description("One of the binary phased parent haplotype files.")
			.build();
	private PluginParameter<String> phased2 = new PluginParameter.Builder<>("phased2", null, String.class)
			.dependentOnParameter(combine)
			.guiName("Second Phased File")
			.inFile()
			.description("The other binary phased parent haplotype file.")
			.build();
	private PluginParameter<String> combineOut = new PluginParameter.Builder<>("combineout", null, String.class)
			.dependentOnParameter(combine)
			.guiName("Combine Output Filename")
			.outFile()
			.description("The name of the combined output file. A .bin extension will be added. If this is blank, the files will not be combined and only a report comparing the two files will be generated.")
			.build();


	private PluginParameter<String> separator3 = PluginParameter.getLabelInstance("Convert binary haplotypes to text ------------------------------");
	private PluginParameter<Boolean> convert = new PluginParameter.Builder<>("convert", false, Boolean.class)
			.guiName("Convert Phasing")
			.description("Converts a binary phasing results file to text.")
			.build();
	//dependent on convert
	private PluginParameter<String> phasedIn = new PluginParameter.Builder<>("binaryinput", null, String.class)
			.dependentOnParameter(convert)
			.guiName("File to convert")
			.inFile()
			.description("The name of the binary file to be converted.")
			.build();
	private PluginParameter<String> convertOut = new PluginParameter.Builder<>("convertout", null, String.class)
			.dependentOnParameter(convert)
			.guiName("Combine Output Filename")
			.outFile()
			.description("The name of the converted output file. A .txt extension will be added if not present.")
			.build();
	
	private PluginParameter<String> separator4 = PluginParameter.getLabelInstance("Correct self phasing ------------------------------");
	private PluginParameter<Boolean> correct = new PluginParameter.Builder<>("correct", false, Boolean.class)
			.guiName("Correct Self Phasing")
			.description("Exchanges sections of chromosomes so that self phase matches cross phase.")
			.build();
	private PluginParameter<String> parentname = new PluginParameter.Builder<>("parent", null, String.class)
			.dependentOnParameter(correct)
			.guiName("Parent Name")
			.description("The parent for which phase will be corrected. Required for phase correction.")
			.build();
	private PluginParameter<String> chrname = new PluginParameter.Builder<>("chrom", null, String.class)
			.dependentOnParameter(correct)
			.guiName("Chromosome Name")
			.description("The chromosome for which phase will be corrected. Required for phase correction.")
			.build();
	private PluginParameter<String> selfPhased = new PluginParameter.Builder<>("selfPhased", null, String.class)
			.dependentOnParameter(correct)
			.guiName("Self Phased File")
			.inFile()
			.description("Phased parent haplotype files from self progeny.")
			.build();
	private PluginParameter<String> crossPhased = new PluginParameter.Builder<>("crossPhased", null, String.class)
			.dependentOnParameter(correct)
			.guiName("Cross Phased File")
			.inFile()
			.description("Phased parent haplotype files from self progeny.")
			.build();
	private PluginParameter<String> selfOut = new PluginParameter.Builder<>("selfOut", null, String.class)
			.dependentOnParameter(correct)
			.guiName("Corrected Output Filename")
			.outFile()
			.description("The name of the corrected self phased output file. A .bin extension will be added. If this is blank, the files will not be combined and only a report comparing the two files will be generated.")
			.build();

	
	public ParentPhasingPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
	}
	
	
	@Override
	protected void preProcessParameters(DataSet input) {
		
		if (phaseParents.value()) {
			if (rephaseParents.value()) throw new IllegalArgumentException("Both phase and rephase cannot be chosen. Check at most one of those.");
			if (input.getDataOfType(GenotypeTable.class).size() != 1) throw new IllegalArgumentException("Phasing parents requires exactly one genotype dataset as input.");
			if (outputFile.value() == null || outputFile.value().length() < 1) throw new IllegalArgumentException("No output file name for phaseParents.");
			if (parentageFile.value() == null || parentageFile.value().length() < 1) throw new IllegalArgumentException("No parentage file name for phaseParents.");
		}
		
		if (combine.value()) {
			if (phased1.value() == null || phased1.value().trim().length() < 1) throw new IllegalArgumentException("Combining haplotypes phased1 input filename is missing.");
			if (phased2.value() == null || phased2.value().trim().length() < 1) throw new IllegalArgumentException("Combining haplotypes phased2 input filename is missing.");
			if (combineOut.value() == null || combineOut.value().trim().length() < 1) throw new IllegalArgumentException("Combining haplotypes output filename is missing.");
		}
		
		if (convert.value()) {
			if (phasedIn.value() == null || phasedIn.value().trim().length() < 1) throw new IllegalArgumentException("Converting binary to text requires exactly one input file.");
			if (convertOut.value() == null || convertOut.value().trim().length() < 1) throw new IllegalArgumentException("Converting binary to text requires exactly one output file.");
		}
		
		if (rephaseParents.value()) {
			if (parentageFile.value() == null || parentageFile.value().length() < 1) throw new IllegalArgumentException("No parentage file name for rephaseParents.");
			if (parentCallFilename.value() == null || parentCallFilename.value().length() < 1) throw new IllegalArgumentException("No parent call file name for rephaseParents.");
			if (parentHaplotypeFilename.value() == null || parentHaplotypeFilename.value().length() < 1) throw new IllegalArgumentException("No parent haplotype input file name for rephaseParents.");
			if (rephaseOutFile.value() == null || rephaseOutFile.value().length() < 1) throw new IllegalArgumentException("No parent haplotype output file name for rephaseParents.");
		}
		
		if (correct.value()) {
			if (selfPhased.value() == null || selfPhased.value().trim().length() < 1) throw new IllegalArgumentException("Correcting haplotypes self input filename is missing.");
			if (crossPhased.value() == null || crossPhased.value().trim().length() < 1) throw new IllegalArgumentException("Correcting haplotypes cross input filename is missing.");
			if (selfOut.value() == null || selfOut.value().trim().length() < 1) throw new IllegalArgumentException("Correcting haplotypes output filename is missing.");
			if (parentname.value() == null || parentname.value().trim().length() < 1) throw new IllegalArgumentException("Correcting haplotypes requires a parent name.");
			if (chrname.value() == null || chrname.value().trim().length() < 1) throw new IllegalArgumentException("Correcting haplotypes requires a chromosome.");
		}

	}

	@Override
	public DataSet processData(DataSet input) {
		List<Datum> resultList = new ArrayList<>();
		GenotypeTable myGeno = (GenotypeTable) input.getDataOfType(GenotypeTable.class).get(0).getData();
		
		if (phaseParents.value()) {
			Path parentPath = Paths.get(parentageFile.value());
			Path savepath = Paths.get(appendbin(outputFile.value()));
			
			if (selfonly.value()) {
				int window = windowSize.value();
				double minNotMiss = 1 - maxMissing.value();
				SelfedHaplotypeFinder shf = new SelfedHaplotypeFinder(window, minNotMiss);
				shf.setGenotype(myGeno);
				shf.phaseSelfedParents(parentPath, savepath);
			} else {
				PhaseHighCoverage phc = new PhaseHighCoverage(myGeno);
				phc.setParentage(parentageFile.value());
				phc.phaseParentsUsingAllAvailableProgeny(2.0, savepath);
			}
		}
		
		if(rephaseParents.value()) {
			RephaseParents rp = new RephaseParents(myGeno, parentCallFilename.value(),parentageFile.value(), parentHaplotypeFilename.value());
			  ImputationUtils.serializePhasedHaplotypes(rp.rephaseUsingCrossProgeny(), rephaseOutFile.value());
		}
		
		if (combine.value()) {
			String comment = String.format("Comparison of phasing for:\n%s\n%s.", phased1.value(), phased2.value());
			Datum report = new Datum("Phase Comparison Report", comparePhasing(myGeno), comment);
			resultList.add(report);
			if (combineOut.value() != null && combineOut.value().trim().length() > 0) {
				mergePhasedHaplotypes(myGeno);
			}
		}
		
		if (convert.value()) {
			if (convertOut.value() != null && convertOut.value().trim().length() > 0) {
				formatPhasedDataAsText(myGeno);
			}
		}
		
		if (correct.value()) {
			SelfedHaplotypeFinder shf = new SelfedHaplotypeFinder(50, 1); //parameters are not used for this, so anything works
			shf.setGenotype(myGeno);
			shf.correctSelfPhaseUsingCross(selfPhased.value(), crossPhased.value(), parentname.value(), Integer.parseInt(chrname.value()), selfOut.value());
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
		return "Phase Parents";
	}

	@Override
	public String getToolTipText() {
		return "Phase parents using progeny genotypes";
	}


	@Override
	public String pluginDescription() {
		return "This plugin phases parents of full-sib families using parent and progeny genotypes and a parentage file as input."
				+ "The method can use either selfed families only or all families.";
	}

	private TableReport comparePhasing(GenotypeTable myGenotype) {
		//String head = "parent\tchr\tphased_same\tphased_diff\tmonoPhased\tphasedMono\tsameMono\tdiffMono\n";
		String[] columnNames = new String[]{"parent","chr","phased_same","phased_diff","monoPhased","phasedMono","sameMono","diffMono"};
		TableReportBuilder reportBuilder = TableReportBuilder.getInstance("Phasing Comparison", columnNames);
		
		byte N = GenotypeTable.UNKNOWN_ALLELE;
		Map<String, byte[][]> hapmap1 = loadHaplotypes(phased1.value());
		Map<String, byte[][]> hapmap2 = loadHaplotypes(phased2.value());
		List<String> crossParents = new ArrayList<>(hapmap1.keySet());
		Collections.sort(crossParents);
		int numberOfChromomsomes = myGenotype.numChromosomes();
		
		for (String parent : crossParents) {
			byte[][] firstHaps = hapmap1.get(parent);
			byte[][] secondHaps = hapmap2.get(parent);
			if (secondHaps != null) {
				myLogger.info(String.format("Comparing phasing for %s\n", parent));
				int[] countPhaseSame = new int[numberOfChromomsomes];
				int[] countPhaseDifferent = new int[numberOfChromomsomes];
				int[] countMonoPhased = new int[numberOfChromomsomes]; //self is monomorphic, cross is polymorphic
				int[] countPhasedMono = new int[numberOfChromomsomes]; //selfed is polymorphic, cross is monomorphic
				int[] sameMono = new int[numberOfChromomsomes];
				int[] differentMono = new int[numberOfChromomsomes];
				
				int n = firstHaps[0].length;
				
				for (int s = 0; s < n; s++) {
					if (firstHaps[0][s] == N) continue;
					if (firstHaps[1][s] == N) continue;
					if (secondHaps[0][s] == N) continue;
					if (secondHaps[1][s] == N) continue;

					int chr = myGenotype.chromosome(s).getChromosomeNumber() - 1;
					if (firstHaps[0][s] == firstHaps[1][s]) {
						if (secondHaps[0][s] == secondHaps[1][s]) {
							if (firstHaps[0][s] == secondHaps[0][s]) sameMono[chr]++;
							else differentMono[chr]++;
						} else {
							countMonoPhased[chr]++;
						}
					} else if (secondHaps[0][s] == secondHaps[1][s]) {
						countPhasedMono[chr]++;
					} else {
						if (firstHaps[0][s] == secondHaps[0][s]) countPhaseSame[chr]++;
						else countPhaseDifferent[chr]++;
					}
				}
				
				for (int c = 0; c < numberOfChromomsomes; c++) {
					int same, diff;
					if (countPhaseSame[c] >= countPhaseDifferent[c]) {
						same = countPhaseSame[c];
						diff = countPhaseDifferent[c];
					} else {
						diff = countPhaseSame[c];
						same = countPhaseDifferent[c];
					}

					Object[] row = new Object[8];
					row[0] = parent;
					row[1] = myGenotype.chromosomes()[c].getName();
					row[2] = same;
					row[3] = diff;
					row[4] = countMonoPhased[c];
					row[5] = countPhasedMono[c];
					row[6] = sameMono[c];
					row[7] = differentMono[c];
					reportBuilder.add(row);
				}
			}
		}

		return reportBuilder.build();
	}
	
    private void mergePhasedHaplotypes(GenotypeTable myGeno) {
    	System.out.println("merging self and cross haplotypes");
    	int[] chrstart = myGeno.chromosomesOffsets();
    	int numberOfChromosomes = chrstart.length;
    	int[] chrend = new int[numberOfChromosomes];
    	System.arraycopy(chrstart, 1, chrend, 0, numberOfChromosomes - 1);
    	chrend[numberOfChromosomes - 1] = myGeno.numberOfSites();
    	
		Map<String, byte[][]> hapmap1 = loadHaplotypes(phased1.value());
		Map<String, byte[][]> hapmap2 = loadHaplotypes(phased2.value());
		Map<String, byte[][]> combinedHapmap = new HashMap<>();
		
		//combine the parent lists
		//for each parent if only one type is present use that
		//if both types are present use consensus
		
		TreeSet<String> parentSet = new TreeSet<>();
		parentSet.addAll(hapmap1.keySet());
		parentSet.addAll(hapmap2.keySet());
		
		for (String parent : parentSet) {
			byte[][] hap1 = hapmap1.get(parent);
			byte[][] hap2 = hapmap2.get(parent);
			Optional<byte[][]> combinedhap = smashHap(hap1, hap2, chrstart, chrend, parent);
			if (combinedhap.isPresent()) combinedHapmap.put(parent, combinedhap.get());
		}

		storeHaplotypes(combinedHapmap, appendbin(combineOut.value()));
		System.out.println("Finished merging self and cross haplotypes");
    }

    private Optional<byte[][]> smashHap(byte[][] hap0, byte[][] hap1, int[] chrstart, int[] chrend, String parent) {
    	//hap0 will be the one used if there is any conflict
    	
    	if (hap0 == null && hap1 == null) return Optional.empty();
    	if (hap0 == null) return Optional.of(hap1);
    	if (hap1 == null) return Optional.of(hap0);
    	int nsites = hap0[0].length;
    	byte[][] combined = new byte[2][];
    	combined[0] = Arrays.copyOf(hap0[0], nsites);
    	combined[1] = Arrays.copyOf(hap0[1], nsites);
    	int numberOfChromosomes = chrstart.length;
    	
    	for (int c = 0; c < numberOfChromosomes; c++) {
    		byte[] chap00 = Arrays.copyOfRange(hap0[0], chrstart[c], chrend[c]);
    		byte[] chap01 = Arrays.copyOfRange(hap0[1], chrstart[c], chrend[c]);
    		byte[] chap10 = Arrays.copyOfRange(hap1[0], chrstart[c], chrend[c]);
    		byte[] chap11 = Arrays.copyOfRange(hap1[1], chrstart[c], chrend[c]);
    		int[] comp00 = compareHaplotypes(chap00, chap10);
    		int[] comp01 = compareHaplotypes(chap00, chap11);
    		int[] comp10 = compareHaplotypes(chap01, chap10);
    		int[] comp11 = compareHaplotypes(chap01, chap11);
    		
    		//test whether the haplotype orders match
    		double originalMatch = (comp00[0] + comp11[0]) / (double)(comp00[1] + comp11[1]);
    		double reverseMatch = (comp01[0] + comp10[0]) / (double)(comp01[1] + comp10[1]);
    		System.out.printf("%s, chr %d: original order match = %1.3f, reverse match = %1.3f\n", parent, c + 1, originalMatch, reverseMatch);
    		if (originalMatch > reverseMatch && originalMatch > 0.9) {
    			//they do. If more than 0.9 similar combine them, else return
    			for (int s = chrstart[c]; s < chrend[c]; s++) {
    				if (hap0[0][s] == N && hap1[0][s] != N) {
    					combined[0][s] = hap1[0][s];
    					combined[1][s] = hap1[1][s];
    				} else if (hap0[0][s] != N && hap1[0][s] != N) {
    					if (hap0[0][s] != hap1[0][s] || hap0[1][s] != hap1[1][s]) {
    						combined[0][s] = N;
    						combined[1][s] = N;
    					} 
    				}
    			}
    		} else if (reverseMatch > originalMatch && reverseMatch > 0.9){
    			//no, reverse them. If more than 0.9 similar, combine them
    			for (int s = chrstart[c]; s < chrend[c]; s++) {
    				if (hap0[0][s] == N && hap1[0][s] != N) {
    					combined[0][s] = hap1[1][s];
    					combined[1][s] = hap1[0][s];
    				} else if (hap0[0][s] != N && hap1[0][s] != N) {
    					if (hap0[0][s] != hap1[1][s] || hap0[1][s] != hap1[0][s]) {
    						combined[0][s] = N;
    						combined[1][s] = N;
    					} 
    				}
    			}
    		} 
    	}
    	return Optional.of(combined);
    }

    private int[] compareHaplotypes(byte[] h0, byte[] h1) {
    	int totalCount = 0;
    	int sameCount = 0;
    	for (int i = 0; i < h0.length; i++) {
    		if (h0[i] != N & h1[i] != N) {
    			totalCount++;
    			if (h0[i] == h1[i]) sameCount++;
    		}
    	}
    	return new int[]{sameCount, totalCount};
    }

    private void formatPhasedDataAsText(GenotypeTable myGeno) {
    	Map<String, byte[][]> myHaps = loadHaplotypes(phasedIn.value());
    	List<String> taxa = new ArrayList<>(myHaps.keySet());
    	Collections.sort(taxa);
    	int nsites = myHaps.get(taxa.get(0))[0].length;
    	int ntaxa = taxa.size();
    	
    	try(BufferedWriter bw = Files.newBufferedWriter(Paths.get(appendtxt(convertOut.value())))) {
    		bw.write("snpid\tchr\tpos");
    		for (String name : taxa) bw.write(String.format("\t%s\t%s", name + "_hap1", name + "_hap2"));
    		bw.write("\n");
    		for (int s = 0; s < nsites; s++) {
    			String snpname = myGeno.siteName(s);
    			String chrname = myGeno.chromosomeName(s);
    			int pos = myGeno.chromosomalPosition(s);
    			bw.write(String.format("%s\t%s\t%d", snpname, chrname, pos));
    			for (String taxonName : taxa) {
    				byte[][] haps = myHaps.get(taxonName);
    				bw.write(String.format("\t%s\t%s", NucleotideAlignmentConstants.getHaplotypeNucleotide(haps[0][s]), 
    						NucleotideAlignmentConstants.getHaplotypeNucleotide(haps[1][s])));
    			}
        		bw.write("\n");
    		}
    	} catch(IOException e) {
    		e.printStackTrace();
    	}
    	System.out.println("Finished writing formatted data to file.");
    }

	private GenotypeTable loadGenotype(String filename) {
		FileLoadPlugin flp = new FileLoadPlugin(null, false);
		flp.setTheFileType(TasselFileType.Unknown);
		flp.setOpenFiles(new File[]{new File(filename)});
		return (GenotypeTable) flp.performFunction(null).getData(0).getData();
	}
	
	private Map<String, byte[][]> loadHaplotypes(String filename) {
    	try {
    		FileInputStream fis = new FileInputStream(new File(filename));
            ObjectInputStream ois = new ObjectInputStream(fis);
            Map<String, byte[][]> phasedHaps = (Map<String, byte[][]>) ois.readObject();
            ois.close();
            return phasedHaps;
		} catch (IOException | ClassNotFoundException e) {
			throw new RuntimeException(e);
		}
	}
	
	private void storeHaplotypes(Map<String, byte[][]> hapmap, String name) {
		try {
			FileOutputStream fos = new FileOutputStream(name);
			ObjectOutputStream oos = new ObjectOutputStream(fos);
			oos.writeObject(hapmap);
			oos.close();
		} catch (IOException e) {
			throw new RuntimeException("Unable to save phased haplotypes.", e);
		}

	}
	
	private String appendbin(String original) {
		if (original.endsWith(".bin")) return original;
		else return original + ".bin";
	}

	private String appendtxt(String original) {
		if (original.endsWith(".txt")) return original;
		else return original + ".txt";
	}

}
