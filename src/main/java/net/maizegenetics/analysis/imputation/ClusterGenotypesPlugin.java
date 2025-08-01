package net.maizegenetics.analysis.imputation;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.swing.ImageIcon;

import net.maizegenetics.analysis.clustering.Haplotype;
import net.maizegenetics.analysis.clustering.HaplotypeCluster;
import net.maizegenetics.analysis.clustering.HaplotypeClusterer;
import net.maizegenetics.analysis.data.GetTaxaListPlugin;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.OpenBitSet;

public class ClusterGenotypesPlugin extends AbstractPlugin {
	private GenotypeTable myGenotype;
	private String dataName;
	private final byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
	private static final Pattern comma_space = Pattern.compile("[,\\s]+");

	PluginParameter<Boolean> useSiteList = new PluginParameter.Builder<>("useSiteList", false, Boolean.class)
			.description("If true, use the list of sites. If false, use start site and number of sites. (Default = false)")
			.guiName("Use Site List")
			.build();

	PluginParameter<Integer> startSiteIndex = new PluginParameter.Builder<>("startSite", 0, Integer.class)
			.description("Start the cluster at this site")
			.guiName("Start Site")
			.dependentOnParameter(useSiteList, false)
			.build();

	PluginParameter<Integer> numberOfSites = new PluginParameter.Builder<>("nsites", 30, Integer.class)
			.description("Cluster this many sites. (Default = 30)")
			.guiName("Number of Sites")
			.dependentOnParameter(useSiteList, false)
			.build();

	PluginParameter<String> siteList = new PluginParameter.Builder<>("sites", null, String.class)
			.description("Comma separated list of sites. Example: 1-5, 8, 10-12.")
			.guiName("Site List")
			.dependentOnParameter(useSiteList)
			.build();

	PluginParameter<Integer> maxDiff = new PluginParameter.Builder<>("maxDiff",
			0, Integer.class)
			.description("If the distance between two haplotypes is less than or equal to maxDiff, "
							+ "the haplotypes can be placed in the same cluster. (Default = 0)")
			.guiName("Max Diff")
			.build();

	PluginParameter<Integer> minHap = new PluginParameter.Builder<>("minHap",1, Integer.class)
			.description("Only keep clusters that have minHap members or more. (Default = 1)")
			.guiName("Min Hap")
			.build();

	PluginParameter<Boolean> clusterOnce = new PluginParameter.Builder<>("unique", false, Boolean.class)
			.description("Place each haplotype in only one cluster. (Default = false)")
			.guiName("One cluster per Haplotype")
			.build();
	
	PluginParameter<Boolean> noMissingValues = new PluginParameter.Builder<>("nomiss", false, Boolean.class)
			.description("Use only haplotypes with no missing values. (Default = false)")
			.guiName("No Missing Values")
			.build();

	PluginParameter<Boolean> makePhenotype = new PluginParameter.Builder<>("convert", false, Boolean.class)
			.description("Should the two biggest clusters be converted to a Phenotype? (Default = false)")
			.guiName("Make Phenotype")
			.build();

	public ClusterGenotypesPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
	}

	@Override
	protected void preProcessParameters(DataSet input) {
		List<Datum> genotypes = input.getDataOfType(GenotypeTable.class);
		if (genotypes.size() != 1)
			throw new IllegalArgumentException(
					"Exactly one genotype data set must be selected");
		myGenotype = (GenotypeTable) genotypes.get(0).getData();
		dataName = genotypes.get(0).getName();
	}

	@Override
	public DataSet processData(DataSet input) {
		List<Position> filterPositions;
		ArrayList<Haplotype> hapList = new ArrayList<>();
		int ntaxa = myGenotype.numberOfTaxa();

		if (useSiteList.value()) {
			filterPositions = new ArrayList<>();
			List<Integer> siteNumberList = parseSiteList();
			for (Integer sitenumber : siteNumberList) {
				filterPositions.add(myGenotype.positions().get(sitenumber));
			}
			int nsites = siteNumberList.size();
			for (int t = 0; t < ntaxa; t++) {
				byte[] haplotype = new byte[nsites];
				for (int s = 0; s < nsites; s++) {
					haplotype[s] = myGenotype
							.genotype(t, siteNumberList.get(s));
				}
				if (!noMissingValues.value()) {
					hapList.add(new Haplotype(haplotype));
				} else if (hasNoMissingSites(haplotype)) {
					hapList.add(new Haplotype(haplotype));
				}
			}
		} else {
			int startSite = startSiteIndex.value();
			int endSite = startSite + numberOfSites.value();

			filterPositions = myGenotype.positions()
					.subList(startSite, endSite);
			for (int t = 0; t < ntaxa; t++) {
				byte[] haplotype = myGenotype.genotypeRange(t, startSite,
						endSite);
				if (!noMissingValues.value()) {
					hapList.add(new Haplotype(haplotype));
				} else if (hasNoMissingSites(haplotype)) {
					hapList.add(new Haplotype(haplotype));
				}
			}
		}

		PositionList filterPosList = PositionListBuilder
				.getInstance(filterPositions);

		if (hapList.size() < 2)
			throw new RuntimeException("No haplotypes to cluster");
		HaplotypeClusterer clusterMaker = new HaplotypeClusterer(hapList);
		clusterMaker.makeClusters();
		
		if (clusterOnce.value()) 
		    clusterMaker.moveAllHaplotypesToBiggestCluster(maxDiff.value());
		else 
		    clusterMaker.mergeClusters(maxDiff.value());
		
		clusterMaker.sortClusters();
		ArrayList<HaplotypeCluster> clusterList = clusterMaker.getClusterList();
		GenotypeTableBuilder genoBuilder = GenotypeTableBuilder
				.getTaxaIncremental(filterPosList);
		int count = 0;
		int minSize = minHap.value();
		for (HaplotypeCluster cluster : clusterList) {
			String name = String.format("Cluster%d:%d(%1.1f)", count++,
					cluster.getSize(), cluster.getScore());
			if (cluster.getSize() >= minSize) 
			    genoBuilder.addTaxon(new Taxon(name), cluster.getMajorityHaplotype());
		}

		List<Datum> resultList = new ArrayList<>();
		String comment = String.format(
				"Clusters built from %s\nusing maxDiff = %d and minHap = %d.",
				dataName, maxDiff.value(), minHap.value());
		Datum result = new Datum("Clusters_" + dataName, genoBuilder.build(),
				comment);
		resultList.add(result);

		// TODO provide a Fishers exact test for SNPs
		if (makePhenotype.value()) {
			int nrows = clusterMaker.getClusterList().get(0).getSize();
			nrows += clusterMaker.getClusterList().get(1).getSize();
			List<Taxon> taxa = new ArrayList<>(nrows);
			float[] code = new float[nrows];
			OpenBitSet missing = new OpenBitSet(nrows);
			count = 0;
			for (Haplotype hap : clusterMaker.getClusterList().get(0)
					.getHaplotypeList()) {
				int taxonIndex = hap.taxonIndex;
				taxa.add(myGenotype.taxa().get(taxonIndex));
				code[count++] = 0;
			}
			for (Haplotype hap : clusterMaker.getClusterList().get(1)
					.getHaplotypeList()) {
				int taxonIndex = hap.taxonIndex;
				taxa.add(myGenotype.taxa().get(taxonIndex));
				code[count++] = 1;
			}

			List<PhenotypeAttribute> attrList = new ArrayList<>();
			List<ATTRIBUTE_TYPE> typeList = new ArrayList<>();
			attrList.add(new TaxaAttribute(taxa));
			attrList.add(new NumericAttribute("genotype", code, missing));
			typeList.add(ATTRIBUTE_TYPE.taxa);
			typeList.add(ATTRIBUTE_TYPE.data);
			Phenotype clusterCode = new PhenotypeBuilder()
					.fromAttributeList(attrList, typeList).build().get(0);
			comment = String.format(
					"Phenotype of two largest clusters \nfrom %s", dataName);
			result = new Datum("Coded_Cluster_" + dataName, clusterCode,
					comment);
			resultList.add(result);
		}

		return new DataSet(resultList, this);
	}

	private boolean hasNoMissingSites(byte[] geno) {
		boolean noMissing = true;
		for (byte b : geno) {
			if (b == NN)
				noMissing = false;
		}
		return noMissing;
	}

	private List<Integer> parseSiteList() {

		List<Integer> siteNumbers = new ArrayList<>();
		String siteListString = siteList.value();
		try {
			String[] items = comma_space.split(siteListString);
			for (String str : items) {
				str = str.trim();
				if (str.contains("-")) {
					String[] strints = str.split("-");
					int start = Integer.parseInt(strints[0]);
					int end = Integer.parseInt(strints[1]);
					for (int i = start; i <= end; i++)
						siteNumbers.add(i);
				} else {
					if (str.length() > 0)
						siteNumbers.add(Integer.parseInt(str));
				}
			}
		} catch (Exception e) {
			throw new RuntimeException("Malformed site list: " + siteListString);
		}

		return siteNumbers;
	}

	@Override
	public ImageIcon getIcon() {
        URL imageURL = GetTaxaListPlugin.class.getResource("/net/maizegenetics/analysis/images/pca.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
	}

	@Override
	public String getButtonName() {
		return "Cluster Genotypes";
	}

	@Override
	public String getToolTipText() {
		return "Cluster the genotypes in a dataset";
	}

}
