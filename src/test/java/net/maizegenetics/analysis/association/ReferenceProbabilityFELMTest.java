package net.maizegenetics.analysis.association;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.score.ReferenceProbabilityBuilder;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.GenotypePhenotypeBuilder;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportTestUtils;
import net.maizegenetics.util.TableReportUtils;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class ReferenceProbabilityFELMTest {
	GenotypeTable myGenotype;
//	static {LoggingUtils.setupDebugLogging();}
	
	@Before
	public void setUp() throws Exception {
		String genotypeFilename = TutorialConstants.HAPMAP_CHR_9_10_FILENAME;
		myGenotype = ImportUtils.readFromHapmap(genotypeFilename);
		
		//convert scores to reference probability
		int nsites = 100;
		int ntaxa = myGenotype.numberOfTaxa();
		ArrayList<float[]> probValueList = new ArrayList<float[]>();
		for (int t = 0; t < ntaxa; t++) probValueList.add(new float[nsites]);
			
		ReferenceProbabilityBuilder refBuilder = ReferenceProbabilityBuilder.getInstance(ntaxa, nsites, myGenotype.taxa());
		PositionListBuilder positionBuilder = new PositionListBuilder();
		for (int s = 0; s < 100; s++) {
			positionBuilder.add(myGenotype.positions().get(s));
			byte minor = myGenotype.minorAllele(s);
			byte[] geno = myGenotype.genotypeAllTaxa(s);
			for (int t = 0; t < ntaxa; t++) {
				if (geno[t] == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) probValueList.get(t)[s] = Float.NaN;
				else {
					byte[] alleles = GenotypeTableUtils.getDiploidValues(geno[t]);
					if (alleles[0] == minor) probValueList.get(t)[s] += 0.5;
					if (alleles[1] == minor) probValueList.get(t)[s] += 0.5;
				}
			}
		}
		
		int count = 0;
		for (float[] probs : probValueList) refBuilder.addTaxon(count++, probs);
		
		myGenotype = GenotypeTableBuilder.getInstance(null, positionBuilder.build(), myGenotype.taxa(), null, null, refBuilder.build(), null, null);
	}

	@After
	public void tearDown() throws Exception {
		
	}

	@Test
	public void testWithoutReplication() {
		String filename = TutorialConstants.TRAITS_FILENAME;
		Phenotype myPhenotype = new PhenotypeBuilder().fromFile(filename).keepAttributes(new int[]{0,1}).build().get(0);
		GenotypePhenotype myGenoPheno = new GenotypePhenotypeBuilder().genotype(myGenotype).phenotype(myPhenotype).intersect().build();

		ReferenceProbabilityFELM rplm = new ReferenceProbabilityFELM(new Datum("genopheno",myGenoPheno,"genopheno"), null);
		rplm.solve();
		TableReport siteReport = rplm.siteReport();
		String expectedResults = GeneralConstants.EXPECTED_RESULTS_DIR + "referenceProbFELM.txt";

		//to save the results to use for validation, uncomment the following line
//		TableReportUtils.saveDelimitedTableReport(siteReport, "\t", new File(expectedResults));

		TableReport expectedSiteReport = TableReportUtils.readDelimitedTableReport(expectedResults, "\t");

		File tempFile; 
		try {
			tempFile = File.createTempFile("test", ".txt");
			TableReportUtils.saveDelimitedTableReport(siteReport, "\t", tempFile);
			TableReport observedSiteReport = TableReportUtils.readDelimitedTableReport(tempFile.getPath(), "\t");
			TableReportTestUtils.compareTableReports(expectedSiteReport, observedSiteReport);
		} catch(IOException e) {
			e.printStackTrace();
		}
	}

//	@Test
	public void testWithReplication() {
		String filename = GeneralConstants.DATA_DIR + "CandidateTests/earht_with_rep.txt";
		Phenotype myPhenotype = new PhenotypeBuilder().fromFile(filename).build().get(0);
		GenotypePhenotype myGenoPheno = new GenotypePhenotypeBuilder().genotype(myGenotype).phenotype(myPhenotype).intersect().build();

		ReferenceProbabilityFELM rplm = new ReferenceProbabilityFELM(new Datum("genopheno",myGenoPheno,"genopheno"), null);
		rplm.solve();
		TableReport siteReport = rplm.siteReport();
		
		String expectedResults = GeneralConstants.EXPECTED_RESULTS_DIR + "referenceProbReplicatedTaxa.txt";

		//to save the results to use for validation, uncomment the following line
//		TableReportUtils.saveDelimitedTableReport(siteReport, "\t", new File(expectedResults));
		
		TableReport expectedSiteReport = TableReportUtils.readDelimitedTableReport(expectedResults, "\t");

		File tempFile; 
		try {
			tempFile = File.createTempFile("test", null);
			TableReportUtils.saveDelimitedTableReport(siteReport, "\t", tempFile);
			TableReport observedSiteReport = TableReportUtils.readDelimitedTableReport(tempFile.getPath(), "\t");
			TableReportTestUtils.compareTableReports(expectedSiteReport, observedSiteReport);
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
}
