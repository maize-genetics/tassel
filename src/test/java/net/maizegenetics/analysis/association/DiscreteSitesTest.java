package net.maizegenetics.analysis.association;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import javax.swing.JOptionPane;

import net.maizegenetics.analysis.association.MLMPlugin.CompressionType;
import net.maizegenetics.analysis.distance.KinshipPlugin;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.GenotypePhenotypeBuilder;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportTestUtils;
import net.maizegenetics.util.TableReportUtils;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class DiscreteSitesTest {
	GenotypeTable myGenotype;
//	static {LoggingUtils.setupDebugLogging();}

	@Before
	public void setUp() throws Exception {
		String genotypeFilename = TutorialConstants.HAPMAP_CHR_9_10_FILENAME;
		myGenotype = ImportUtils.readFromHapmap(genotypeFilename);
	}

	@After
	public void tearDown() throws Exception {
	}

//	@Test
	public void testReplicatedTaxa() {
		String filename = GeneralConstants.DATA_DIR + "CandidateTests/earht_with_rep.txt";
		Phenotype myPhenotype = new PhenotypeBuilder().fromFile(filename).build().get(0);
		GenotypePhenotype myGenoPheno = new GenotypePhenotypeBuilder().genotype(myGenotype).phenotype(myPhenotype).intersect().build();

		DiscreteSitesFELM dslm = new DiscreteSitesFELM(new Datum("name", myGenoPheno, "comment"), null);
		
		dslm.solve();
		
		TableReport siteReport = dslm.siteReport();
		String expectedResults = GeneralConstants.EXPECTED_RESULTS_DIR + "discreteSitesReplicatedTaxa.txt";

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

	@Test
	public void testUnreplicatedTaxa() {
		String filename = TutorialConstants.TRAITS_FILENAME;
		Phenotype myPhenotype = new PhenotypeBuilder().fromFile(filename).keepAttributes(new int[]{0,1}).build().get(0);
		GenotypePhenotype myGenoPheno = new GenotypePhenotypeBuilder().genotype(myGenotype).phenotype(myPhenotype).intersect().build();

		DiscreteSitesFELM dslm = new DiscreteSitesFELM(new Datum("name", myGenoPheno, "comment"), null);
		dslm.solve();
		TableReport siteReport = dslm.siteReport();
		String expectedResults = GeneralConstants.EXPECTED_RESULTS_DIR + "discreteSitesFELM.txt";

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
	
	@Test
	public void testHeterozygotes() {
	       // AC and CA should be a single genotype class not two different genotype classes
        // As a result the effect table should contain 3 * nsite rows not 4 * nsite rows
        
        //create a taxa list for twenty taxa
        int ntaxa = 20;
        TaxaList myTaxa = TaxaListBuilder.getInstance(ntaxa);
        
        //create genotypes with 10 SNPs and twenty taxa
        int nsites = 10;
        Chromosome chr = Chromosome.instance(1);
        GenotypeTableBuilder gtBuilder = GenotypeTableBuilder.getSiteIncremental(myTaxa);
        for (int i = 0; i < nsites; i++) gtBuilder.addSite(Position.builder(1, i).build(), MLMTest.randomGenotypes(ntaxa, "A", "C"));
        GenotypeTable gt = gtBuilder.build();
        
        //create phenotype for twenty taxa
        List<PhenotypeAttribute> attrList = new ArrayList<>();
        attrList.add(new TaxaAttribute(myTaxa));
        float[] values = new float[] {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
        attrList.add(new NumericAttribute("trait", values));
        List<Phenotype.ATTRIBUTE_TYPE> myTypes = Arrays.asList(Phenotype.ATTRIBUTE_TYPE.taxa, Phenotype.ATTRIBUTE_TYPE.data);
        Phenotype myPheno = new PhenotypeBuilder().fromAttributeList(attrList, myTypes).build().get(0);
        
        GenotypePhenotype genoPheno = new GenotypePhenotypeBuilder().genotype(gt).phenotype(myPheno).build();
        
        DiscreteSitesFELM dslm = new DiscreteSitesFELM(new Datum("name", genoPheno, "comment"), null);
        dslm.solve();
        TableReport effectTestReport = dslm.alleleReport();
        System.out.println("effectTestReport for DiscreteSitesTest (GLM)");
        String colnames = Arrays.stream(effectTestReport.getTableColumnNames()).map(obj -> obj.toString())
                .collect(Collectors.joining(","));
        System.out.println(colnames);
        for (long row = 0; row < effectTestReport.getRowCount(); row++) {
            String rowString = Arrays.stream(effectTestReport.getRow(row)).map(obj -> obj.toString())
                    .collect(Collectors.joining(","));
            System.out.println(rowString);
        }
        
        assertEquals(nsites * 3L, effectTestReport.getRowCount());

	}
	
}
