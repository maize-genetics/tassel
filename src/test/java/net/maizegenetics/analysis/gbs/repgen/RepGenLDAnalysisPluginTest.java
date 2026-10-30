/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.junit.Test;

import com.google.common.collect.Multimap;

import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;

/**
 * Correlates the per-taxon depths of the {@link RepGenSimData} tags. The two alleles of an amplicon are
 * never seen in the same taxon, so their presence/absence correlation is -1 and r-squared is 1.
 *
 * @author lcj34
 *
 */
public class RepGenLDAnalysisPluginTest {

    private static final double DELTA = 1e-5;

    @Test
    public void testRepGenLDAnalysisPlugin() throws Exception {
        RepGenSimData data = RepGenSimData.createUnder("ldAnalysis");
        String db = RepGenLoadSeqToDBPluginTest.loadDB(data, 5);

        new RepGenLDAnalysisPlugin()
                .inputDB(db)
                .minTaxa(RepGenSimData.TAXA.length)
                .performFunction(null);

        RepGenSQLite repGenData = new RepGenSQLite(db);
        try {
            List<Tag> tagList = new ArrayList<>(repGenData.getTags());
            assertEquals(data.expectedTags.size(), tagList.size());
            Multimap<Tag, TagCorrelationInfo> correlations = repGenData.getCorrelationsForTags(tagList);
            assertEquals(tagList.size(), correlations.keySet().size());
            for (Tag tag : tagList) {
                Collection<TagCorrelationInfo> forTag = correlations.get(tag);
                assertEquals("Each tag is correlated with every other tag", tagList.size() - 1, forTag.size());
                RepGenSimData.ExpectedTag expected1 = data.expectedTags.get(tag.sequence());
                for (TagCorrelationInfo tci : forTag) {
                    RepGenSimData.ExpectedTag expected2 = data.expectedTags.get(tci.tag2().sequence());
                    double[] depths1 = Arrays.stream(expected1.depths).asDoubleStream().toArray();
                    double[] depths2 = Arrays.stream(expected2.depths).asDoubleStream().toArray();
                    assertEquals(new PearsonsCorrelation().correlation(depths1, depths2), tci.t1t2_pearson(), DELTA);
                    assertEquals(new SpearmansCorrelation().correlation(depths1, depths2), tci.t1t2_spearman(), DELTA);
                    assertEquals(new PearsonsCorrelation().correlation(presence(depths1), presence(depths2)), tci.pres_abs_pearson(), DELTA);
                    if (expected1.locus == expected2.locus) {
                        assertEquals("Alleles of one amplicon never share a taxon", -1.0, tci.pres_abs_pearson(), DELTA);
                        assertEquals(1.0, tci.r2(), DELTA);
                    }
                }
            }
        } finally {
            repGenData.close();
        }
    }

    @Test
    public void testRSquaredNeedsMinTaxa() throws Exception {
        RepGenSimData data = RepGenSimData.createUnder("ldAnalysisMinTaxa");
        String db = RepGenLoadSeqToDBPluginTest.loadDB(data, 5);

        new RepGenLDAnalysisPlugin()
                .inputDB(db)
                .minTaxa(RepGenSimData.TAXA.length + 1)
                .performFunction(null);

        RepGenSQLite repGenData = new RepGenSQLite(db);
        try {
            List<Tag> tagList = new ArrayList<>(repGenData.getTags());
            Multimap<Tag, TagCorrelationInfo> correlations = repGenData.getCorrelationsForTags(tagList);
            assertEquals(tagList.size() * (tagList.size() - 1), correlations.size());
            for (TagCorrelationInfo tci : correlations.values()) {
                assertTrue("r-squared is not estimated with fewer than minTaxa taxa", Double.isNaN(tci.r2()) || tci.r2() == 0.0);
            }
        } finally {
            repGenData.close();
        }
    }

    private static double[] presence(double[] depths) {
        return Arrays.stream(depths).map(depth -> depth > 0 ? 1 : 0).toArray();
    }

}
