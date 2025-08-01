/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import java.awt.Frame;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;

import javax.swing.ImageIcon;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;

import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;

/**
 * This class takes a rAmpSeq (formerly RepGen) database and for
 * each tag in the tag table, performs the following tag-tag
 * correlations based on the taxa distribution for each tag
 *   (1) tag-tag Pearson's correlation
 *   (2) tag-tag Spearman's correlation
 *   (3) tag-tag presence/absence correlations
 *   (4) r-squared
 *   
 *   The vectors presented to the analysis methods represent
 *   a list of taxa and the number of times the tag was seen
 *   in that taxa.  The presence/absence vectors have a 1 or 0
 *   as values in each slot.
 * @author lcj34
 *
 */
public class RepGenLDAnalysisPlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(RepGenLDAnalysisPlugin.class);
    
    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with tags and taxa distribution").build();
    private PluginParameter<Integer> minTaxa = new PluginParameter.Builder<>("minTaxa", 20, Integer.class).guiName("Min Taxa for RSquared")
            .description("Minimum number of taxa that must be present for R-squared to be calculated.").build();    
    public RepGenLDAnalysisPlugin() {
        super(null, false);
    }

    public RepGenLDAnalysisPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public RepGenLDAnalysisPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    @Override
    public void postProcessParameters() {

        if (myDBFile.isEmpty() || !Files.exists(Paths.get(inputDB()))) {
            throw new IllegalArgumentException("RepGenLDAnalysisPlugin: postProcessParameters: Input DB not set or found");
        }
    }
    
    @Override
    public DataSet processData(DataSet input) {
        long totalTime = System.nanoTime();
        long time=System.nanoTime();
 
        try {           
            System.out.println("RepGenLDAnalysis:processData begin, get all tags/taxadist from db"); 
            RepGenDataWriter repGenData=new RepGenSQLite(inputDB());

            Map <Tag, TaxaDistribution> tagTaxaMap = repGenData.getAllTagsTaxaMap();
            int tagcount = 0;

            int processedTags = 0;
            System.out.println("TIme to get all tags with taxa from db: " + (System.nanoTime() - totalTime)/1e9 + " seconds.\n");
            time = System.nanoTime();
            Multimap<Tag,TagCorrelationInfo> tagTagCorrelations = Multimaps.synchronizedMultimap(HashMultimap.<Tag,TagCorrelationInfo>create());
            System.out.println("\nStart processing tag correlations.  Number of tags in db: " + tagTaxaMap.keySet().size());
            Set<Tag> tagSet = tagTaxaMap.keySet();
            List<Tag> tagList = new ArrayList<Tag>(tagSet);
            
            for (int tidx = 0; tidx < tagList.size(); tidx++) {
            //for (Tag tag1 : tagTaxaMap.keySet()) {
                Tag tag1 = tagList.get(tidx);
                tagcount++;
                processedTags++;
                // get dist for each taxa
                TaxaDistribution tag1TD = tagTaxaMap.get(tag1);
                if (tag1TD == null) {
                    ((RepGenSQLite)repGenData).close();
                    System.out.println("GetTagTaxaDist: got null tagTD at tagcount " + tagcount);
                    return null; // But this should return an error?
                }
                int[] depths1 = tag1TD.depths(); // gives us the depths for each taxon
                // This needs to be doubles for Pearson
                double[] ddepths1 = new double[depths1.length];
                
                // Apparently no shorter method of casting int to double 
                for (int idx = 0; idx < depths1.length; idx++) {
                    ddepths1[idx] = (double)depths1[idx];
                }

                double[] depthsPrime1 = new double[depths1.length];
                for (int idx = 0; idx < depthsPrime1.length; idx++) {
                    // boolean - presence/absence
                    if (ddepths1[idx] > 0) depthsPrime1[idx] = 1;
                    else depthsPrime1[idx] = 0;
                }
                final int tIdxFinal = tidx; // IntStream forEach must have "final" variable, tidx is not final
                IntStream.range(tidx+1,tagList.size()).parallel().forEach(item -> {
                    calculateCorrelations(tagTagCorrelations, tagTaxaMap, 
                            tagList.get(tIdxFinal), tagList.get(item), ddepths1, depthsPrime1);
                });
                
                if (tagcount > 1000) {
                    // comment out when run for real!
                    System.out.println("FInished processing " + processedTags + " tags, this set took " + (System.nanoTime() - time)/1e9 + " seconds, now load to db ..." );
                    time = System.nanoTime();
                    repGenData.putTagTagCorrelationMatrix(tagTagCorrelations);
                    System.out.println("Loading DB took " + (System.nanoTime() - time)/1e9 + " seconds.\n");
                    tagcount = 0;
                    tagTagCorrelations.clear(); // start fresh with next 1000
                    time = System.nanoTime();
                }
            
            }
            if (tagcount > 0 ) {
                System.out.println("Finished processing last tags, load to DB");
                time = System.nanoTime();
                repGenData.putTagTagCorrelationMatrix(tagTagCorrelations);
                System.out.println("Loading DB took " + (System.nanoTime() - time)/1e9 + " seconds.\n");
                tagcount = 0;
                tagTagCorrelations.clear(); // start fresh with next 1000
            }
            System.out.println("Total number of tags processed: " + processedTags );           
            
            ((RepGenSQLite)repGenData).close();
            
        } catch (Exception exc) {
            System.out.println("RepGenLDAnalysis:process_data:  processing error");
            exc.printStackTrace();
        }
        System.out.println("Process took " + (System.nanoTime() - totalTime)/1e9 + " seconds.\n");
        return null;
    }
    
    public void calculateCorrelations(Multimap<Tag,TagCorrelationInfo> tagTagCorrelations, Map <Tag, TaxaDistribution> tagTaxaMap, 
            Tag tag1, Tag tag2, double[] ddepths1,double[] depthsPrime1) {
        TaxaDistribution tag2TD = tagTaxaMap.get(tag2);
        if (tag2TD == null) {                       
            System.out.println("GetTagTaxaDist: got null tagTD for sequence " + tag2.sequence());
            return ;
        }
        
        // I need doubles below !!
        int[] depths2 = tag2TD.depths(); // gives us the depths for each taxon
        double[] ddepths2 = new double[depths2.length];
        // Apparently no shorter method of casting int to double 
        for (int idx = 0; idx < depths2.length; idx++) {
            ddepths2[idx] = (double)depths2[idx];
        }
        double[] depthsPrime2 = new double[depths2.length];
        for (int idx = 0; idx < depthsPrime1.length; idx++) {
            // boolean - presence/absence
            if (ddepths2[idx] > 0) depthsPrime2[idx] = 1;
            else depthsPrime2[idx] = 0;
        }
        // From analysis.association.GenomicSelectionPlugin.java:
        PearsonsCorrelation Pearsons = new PearsonsCorrelation();
        double p1 = Pearsons.correlation(ddepths1,ddepths2);
        
        SpearmansCorrelation Spearmans = new SpearmansCorrelation();
        double spearval = Spearmans.correlation(ddepths1, ddepths2);
        
        double p2 = Pearsons.correlation(depthsPrime1, depthsPrime2);
        
        // Count number of times both tags appeared in a taxa, number
        // of times neither tag appeared in a taxa, number of times
        // tag1 appeared but not tag2, and number of times tag2 appeared by not tag1
        int t1Nott2 = 0;
        int t2Nott1 = 0;
        int neither = 0;
        int both = 0;
        
        for (int didx = 0; didx < depthsPrime2.length; didx++) {
            if (depthsPrime1[didx] > 0 && depthsPrime2[didx] > 0) both++;
            if (depthsPrime1[didx] == 0 && depthsPrime2[didx] == 0) neither++;
            if (depthsPrime1[didx] > 0 && depthsPrime2[didx] == 0) t1Nott2++;
            if (depthsPrime1[didx] == 0 && depthsPrime2[didx] > 0) t2Nott1++;
        }
        // Calculate r-squared based on presence/absence of tags at each taxa.
        double r2 = LinkageDisequilibrium.calculateRSqr(neither, t1Nott2, t2Nott1, both, minTaxa());
        TagCorrelationInfo tci = new TagCorrelationInfo(tag2,p1,spearval,p2,r2);
        tagTagCorrelations.put(tag1, tci);       
    }
    
    @Override
    public ImageIcon getIcon() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getButtonName() {
        // TODO Auto-generated method stub
        return null;
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
     public static void main(String[] args) {
         GeneratePluginCode.generate(RepGenLDAnalysisPlugin.class);
     }
    @Override
    public String getToolTipText() {
        // TODO Auto-generated method stub
        return null;
    }
    /**
     * Input database file with tags and taxa distribution
     *
     * @return Input DB
     */
    public String inputDB() {
        return myDBFile.value();
    }

    /**
     * Set Input DB. Input database file with tags and taxa
     * distribution
     *
     * @param value Input DB
     *
     * @return this plugin
     */
    public RepGenLDAnalysisPlugin inputDB(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }
    
    /**
     * Minimum number of taxa that must be present for R-squared
     * to be calculated.
     *
     * @return Min Taxa for RSquared
     */
    public Integer minTaxa() {
        return minTaxa.value();
    }

    /**
     * Set Min Taxa for RSquared. Minimum number of taxa that
     * must be present for R-squared to be calculated.
     *
     * @param value Min Taxa for RSquared
     *
     * @return this plugin
     */
    public RepGenLDAnalysisPlugin minTaxa(Integer value) {
        minTaxa = new PluginParameter<>(minTaxa, value);
        return this;
    }
}
