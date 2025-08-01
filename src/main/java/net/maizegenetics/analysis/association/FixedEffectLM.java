package net.maizegenetics.analysis.association;

import java.util.List;

import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.TableReport;

public interface FixedEffectLM {

	/**
	 * Initializes builders for the site and allele reports. The report names and column headers for a specific analysis are created in this method.
	 */
	public void initializeReportBuilders();

	/**
	 * Solves the linear model
	 */
	public void solve();

	/**
	 * @return		the site report, which contains tests of significance for each site
	 */
	public TableReport siteReport();

	/**
	 * @return		the allele report, which contains estimates and observation numbers for each allele at each site
	 */
	public TableReport alleleReport();

	/**
	 * @return		the site and allele reports as a list of Datum
	 */
	public List<Datum> datumList();
	
	/**
	 * @param permute	if true, conduct an experiment-wise permutation test
	 * @param nperm		the number of permutations to be run
	 */
	public void permutationTest(boolean permute, int nperm);
	
	/**
	 * @param maxP		test results with p > maxP will not be reported.
	 */
	public void maxP(double maxP);
	
	/**
	 * @param savefile	results will be saved to this file instead of memory
	 */
	public void siteReportFilepath(String savefile);
	
	/**
	 * @param savefile	results will be saved to this file instead of memory
	 */
	public void alleleReportFilepath(String savefile);

	/**
	 * @param biallelic		If true, only biallelic sites will be included in the analysis
	 */
	public void biallelicOnly(boolean biallelic);
	
	/**
	 * @param minsize		If a genotype class has fewer than minsize observations, it will not be tested.
	 */
	public void minimumClassSize(int minsize);
	
	/**
	 * @param siteStats		If true, site statistics will be output as a file.
	 */
	public void saveSiteStats(boolean siteStats);
	
	/**
	 * @param filename		The filename to which the site statistics will be written
	 */
	public void siteStatsFile(String filename);
	
	/**
	 * @param append If true, additive and dominance effects will be add to the stats report for bi-allelic loci.
	 */
	public void appendAddDom(boolean append);
}