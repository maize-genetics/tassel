package net.maizegenetics.analysis.association;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import net.maizegenetics.analysis.association.MLMPlugin.CompressionType;
import net.maizegenetics.analysis.distance.EndelmanDistanceMatrix;
import net.maizegenetics.analysis.distance.KinshipPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
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
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportTestUtils;
import net.maizegenetics.util.TableReportUtils;

import static org.junit.Assert.assertEquals;

import java.util.*;
import java.util.stream.Collectors;

import org.junit.After;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

public class MLMTest {
	private final static String statsNocompNop3d = GeneralConstants.EXPECTED_RESULTS_DIR + "MLM_statistics_nocomp_nop3d.txt";
	private final static String effectsNocompNop3d = GeneralConstants.EXPECTED_RESULTS_DIR + "MLM_effects_nocomp_nop3d.txt";
	private final static String residualsNocompNop3d = GeneralConstants.EXPECTED_RESULTS_DIR + "MLM_residuals_nocomp_nop3d.txt";
	private final static String statsCompp3d = GeneralConstants.EXPECTED_RESULTS_DIR + "MLM_statistics_comp_p3d.txt";
	private final static String effectsCompp3d = GeneralConstants.EXPECTED_RESULTS_DIR + "MLM_effects_comp_p3d.txt";
	private final static String residualsCompp3d = GeneralConstants.EXPECTED_RESULTS_DIR + "MLM_residuals_comp_p3d.txt";
	private final static String compressionCompp3d = GeneralConstants.EXPECTED_RESULTS_DIR + "MLM_compression_comp_p3d.txt";
	private GenotypeTable myGenotype;
	private Phenotype myPhenotype;
	private DistanceMatrix myK;
	private DataSet myDataset;
	
	@Before
	public void setUp() throws Exception {
		GenotypeTable geno = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_CHR_9_10_FILENAME);
		myK = EndelmanDistanceMatrix.getInstance(geno, 6, null);
		geno = FilterGenotypeTable.getInstance(geno, 0, 50);
		myGenotype = new FilterSiteBuilderPlugin().siteMinCount(100).siteMinAlleleFreq(0.1).runPlugin(geno);
		myPhenotype = new PhenotypeBuilder().fromFile(TutorialConstants.TRAITS_FILENAME).keepAttributes(new int[]{0,1}).build().get(0);
		GenotypePhenotype genopheno = new GenotypePhenotypeBuilder().genotype(myGenotype).phenotype(myPhenotype).intersect().build();
		myDataset = new DataSet(new Datum[]{
				new Datum("genotypePhenotype", genopheno, "no comment"),
				new Datum("kinship", myK, "no comment")
				}, null);
		
	}

	@Test
	public void testNoCompressionNoP3D() {
		MLMPlugin myPlugin = new MLMPlugin(null, false);
		myPlugin.setCompressionType(CompressionType.None);
		myPlugin.setUseP3D(false);
		DataSet resultSet = myPlugin.performFunction(myDataset);
		
		TableReport statReferenceReport = TableReportUtils.readDelimitedTableReport(statsNocompNop3d, "\t");
		TableReport effectReferenceReport = TableReportUtils.readDelimitedTableReport(effectsNocompNop3d, "\t");
		TableReport residualReferenceReport = TableReportUtils.readDelimitedTableReport(residualsNocompNop3d, "\t");
		TableReport statTestReport = (TableReport) resultSet.getData(1).getData();
		TableReport effectTestReport  = (TableReport) resultSet.getData(2).getData();
		TableReport residualTestReport  = (TableReport) resultSet.getData(0).getData();
		
		TableReportTestUtils.compareTableReports(statReferenceReport, statTestReport,.01);
		TableReportTestUtils.compareTableReports(effectReferenceReport, effectTestReport,.01);
		TableReportTestUtils.compareTableReports(residualReferenceReport, residualTestReport,.01);
	}

	@Ignore
	@Test
	public void testCompressionP3D() {
		MLMPlugin myPlugin = new MLMPlugin(null, false);
		myPlugin.setCompressionType(CompressionType.Optimum);
		myPlugin.setUseP3D(true);
		DataSet resultSet = myPlugin.performFunction(myDataset);
		
		TableReport statReferenceReport = TableReportUtils.readDelimitedTableReport(statsCompp3d, "\t");
		TableReport effectReferenceReport = TableReportUtils.readDelimitedTableReport(effectsCompp3d, "\t");
		TableReport residualReferenceReport = TableReportUtils.readDelimitedTableReport(residualsCompp3d, "\t");
		TableReport compressionReferenceReport = TableReportUtils.readDelimitedTableReport(compressionCompp3d, "\t"); 
		TableReport statTestReport = (TableReport) resultSet.getData(1).getData();
		TableReport effectTestReport  = (TableReport) resultSet.getData(2).getData();
		TableReport residualTestReport  = (TableReport) resultSet.getData(0).getData();
		TableReport compressionTestReport  = (TableReport) resultSet.getData(3).getData();
		
		TableReportTestUtils.compareTableReports(statReferenceReport, statTestReport,.01);
		TableReportTestUtils.compareTableReports(effectReferenceReport, effectTestReport,.01);
		TableReportTestUtils.compareTableReports(residualReferenceReport, residualTestReport,.01);
		TableReportTestUtils.compareTableReports(compressionReferenceReport, compressionTestReport, .01);
	}

	@Test
	public void testHeterozygotes() {
	    // AC and CA should be a single genotype class not two different genotype classes
	    // As a result the effect table should contain 3 * nsite rows not 4 * nsite rows
	    
	    //create a taxa list for twenty taxa
	    int ntaxa = 20;
	    TaxaList myTaxa = TaxaListBuilder.getInstance(ntaxa);
	    
	    //create genotypes with 10 SNPs and twenty taxa
		// all sites must have 3 genotypes (2 homozygous, 1 heterozygous)
	    int nsites = 10;
	    Chromosome chr = Chromosome.instance(1);
	    GenotypeTableBuilder gtBuilder = GenotypeTableBuilder.getSiteIncremental(myTaxa);
	    for (int i = 0; i < nsites; i++) {
	    	int ngeno = 0;
	    	byte[] geno = new byte[0];
	    	while (ngeno != 4) {
				geno = randomGenotypes(ntaxa, "A", "C");
				Multiset<Integer> genoMultiset = HashMultiset.create();
				for (byte genoValue : geno) genoMultiset.add(Byte.toUnsignedInt(genoValue));
				ngeno = genoMultiset.elementSet().size();
			}

	    	gtBuilder.addSite(Position.builder(1, i).build(), geno);
		}
	    GenotypeTable gt = gtBuilder.build();
	    
	    DataSet kinshipDataset = new KinshipPlugin(null, false).performFunction(DataSet.getDataSet(gt));
	    DistanceMatrix dm = (DistanceMatrix) kinshipDataset.getData(0).getData();
	    
	    //create phenotype for twenty taxa
	    List<PhenotypeAttribute> attrList = new ArrayList<>();
	    attrList.add(new TaxaAttribute(myTaxa));
	    float[] values = new float[] {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
	    attrList.add(new NumericAttribute("trait", values));
	    List<Phenotype.ATTRIBUTE_TYPE> myTypes = Arrays.asList(Phenotype.ATTRIBUTE_TYPE.taxa, Phenotype.ATTRIBUTE_TYPE.data);
	    Phenotype myPheno = new PhenotypeBuilder().fromAttributeList(attrList, myTypes).build().get(0);
	    
	    GenotypePhenotype genoPheno = new GenotypePhenotypeBuilder().genotype(gt).phenotype(myPheno).build();
	    
	    DataSet mlmDataset = new DataSet(new Datum[] {new Datum("genoPheno", genoPheno, "none"), new Datum("kinship", dm, "none")}, null);
	    MLMPlugin myPlugin = new MLMPlugin(null, false);
	    myPlugin.setCompressionType(CompressionType.None);
	    myPlugin.setUseP3D(false);
	    DataSet resultSet = myPlugin.performFunction(mlmDataset);

	    TableReport effectTestReport  = (TableReport) resultSet.getData(2).getData();
	    System.out.println("effectTestReport");
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
	
	public static byte[] randomGenotypes(int n, String nuc1, String nuc2) {
	    Random ran = new Random();
	    byte[] geno = new byte[4];
	    geno[0] = NucleotideAlignmentConstants.getNucleotideDiploidByte(nuc1 + nuc1);
	    geno[1] = NucleotideAlignmentConstants.getNucleotideDiploidByte(nuc1 + nuc2);
	    geno[2] = NucleotideAlignmentConstants.getNucleotideDiploidByte(nuc2 + nuc1);
	    geno[3] = NucleotideAlignmentConstants.getNucleotideDiploidByte(nuc2 + nuc2);
	    byte[] genotype = new byte[n];
	    for (int i = 0; i < n; i++) {
	        genotype[i] = geno[ran.nextInt(4)];
	    }
	    return genotype;
	}
}
