package net.maizegenetics.analysis.data;

import com.google.common.collect.Table;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class GenotypeSummaryPluginTest {
    private GenotypeTable theTBA;

    @Before
    public void setUp() throws Exception {
        String[] seqs = {
                "AACCGGTT",
                "ANCNGNTN",
                "AAACRRTT",
                "CANNRRNN"};
        TaxaListBuilder tLB = new TaxaListBuilder();
        for (int i = 0; i < seqs.length; i++) {
            tLB.add(new Taxon("T" + i));
        }
        PositionListBuilder pLB = new PositionListBuilder();
        for (int i = 0; i < seqs[0].length(); i++) {
            pLB.add(new GeneralPosition.Builder(Chromosome.UNKNOWN, i).build());
        }
        GenotypeCallTable g = GenotypeCallTableBuilder.getInstance(seqs.length, seqs[0].length()).setBases(seqs).build();
        theTBA = GenotypeTableBuilder.getInstance(g, pLB.build(), tLB.build());

    }

    @Test
    public void testTaxaSummaryStats() throws Exception {
        GenotypeSummaryPlugin gsp=new GenotypeSummaryPlugin(null,false);
        DataSet ds=gsp.performFunction(new DataSet(new Datum("TestGenotypeTable",theTBA,""),null));
        TableReport tr=(TableReport)ds.getData(3).getData(); //
        Table<Integer,String,Object> t=TableReportUtils.convertTableReportToGuavaTable(tr);
        Object id=t.get(0, "Number of Sites");
        Assert.assertTrue("Number of Sites wrong", (Integer) t.get(0, "Number of Sites")==8);
        Assert.assertTrue("Number of Sites wrong",(Integer)t.get(1,"Number of Sites")==8);
        Assert.assertTrue("Proportion Heterozygous wrong",(Double)t.get(0,"Proportion Heterozygous")==0.0);
        Assert.assertTrue("Proportion Heterozygous wrong",(Double)t.get(1,"Proportion Heterozygous")==0.0);
        Assert.assertEquals("Proportion Heterozygous wrong",(double)t.get(2,"Proportion Heterozygous"),0.25,0.01);
        Assert.assertEquals("Proportion Heterozygous wrong",(double)t.get(3,"Proportion Heterozygous"),0.50,0.01);
        Assert.assertEquals("Proportion Missing wrong",(double)t.get(1,"Proportion Missing"),0.50,0.01);
        Assert.assertEquals("Proportion Missing wrong",(double)t.get(3,"Proportion Missing"),0.50,0.01);
    }

    @Test
    public void testSiteSummaryStats() throws Exception {
        GenotypeSummaryPlugin gsp=new GenotypeSummaryPlugin(null,false);
        DataSet ds=gsp.performFunction(new DataSet(new Datum("TestGenotypeTable",theTBA,""),null));
        TableReport tr=(TableReport)ds.getData(2).getData(); //
        Table<Integer,String,Object> t=TableReportUtils.convertTableReportToGuavaTable(tr);
//        System.out.println(t.toString());
        Assert.assertEquals("Minor Allele Frequency wrong", 0.25,(Double)t.get(0, "Minor Allele Frequency"),0.01);
        Assert.assertEquals("Minor Allele Frequency wrong", 0.00,(Double)t.get(1, "Minor Allele Frequency"),0.01);
        Assert.assertEquals("Minor Allele Frequency wrong",
                0.333,(Double)t.get(2, "Minor Allele Frequency"),0.01); //test MAF with missing
        Assert.assertEquals("Minor Allele Frequency wrong",
                0.25,(Double)t.get(4, "Minor Allele Frequency"),0.01); //test MAF with hets
        Assert.assertEquals("Proportion Missing wrong",
                0.00,(Double)t.get(0, "Proportion Missing"),0.01); //test missing
        Assert.assertEquals("Proportion Missing wrong",
                0.5,(Double)t.get(7, "Proportion Missing"),0.01); //test missing
    }


}
