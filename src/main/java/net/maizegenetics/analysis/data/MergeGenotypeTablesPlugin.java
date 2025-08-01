/*
 * MergeGenotypeTablesPlugin
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.BasicGenotypeMergeRule;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;

/**
 * Merge alignments into a single alignment containing all taxa and positions.
 * TODO: Add capacity for depth
 *
 * @author Jason Wallace
 * @author Terry Casstevens
 *
 */
public class MergeGenotypeTablesPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(MergeGenotypeTablesPlugin.class);

    public MergeGenotypeTablesPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        List<Datum> inputs = input.getDataOfType(GenotypeTable.class);

        if ((inputs == null) || (inputs.size() < 2)) {
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), "Must select at least two alignments.");
            } else {
                myLogger.warn("performFunction: Must select at least two alignments.");
            }
            return null;
        }

        try {
            GenotypeTable[] alignments = new GenotypeTable[inputs.size()];
            for (int i = 0; i < inputs.size(); i++) {
                alignments[i] = (GenotypeTable) ((Datum) inputs.get(i)).getData();
            }

            GenotypeTable merged = mergeGenotypeTables(alignments);
            DataSet result = new DataSet(new Datum("Merged Genotype Table", merged, null), this);

            fireDataSetReturned(new PluginEvent(result, MergeGenotypeTablesPlugin.class));

            return result;
        } finally {
            fireProgress(100);
        }

    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = MergeGenotypeTablesPlugin.class.getResource("/net/maizegenetics/analysis/images/Merge.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Merge Genotype Tables";
    }

    @Override
    public String getToolTipText() {
        return "Merge Genotype Tables";
    }

    /**
     * Merge an array of GenotypeTables into a single GenotypeTable
     */
    public static GenotypeTable mergeGenotypeTables(GenotypeTable[] alignments) {
        if ((alignments == null) || (alignments.length == 0)) {
            return null;
        }

        //Create the needed TaxaList and PositionList, including only unique entries
        TaxaList masterTaxa = generateMasterTaxaList(alignments);
        PositionList masterPos = generateMasterPositionList(alignments);

        //Make helper hashmaps of taxa and positions for each input alignment
        //Done to to avoid the O(n) search time to find the index of a Taxon or Position in them
        myLogger.info("Creating helper data structures (to speed merging)");
        HashMap[] taxaHashes = makeTaxaHashes(alignments);
        HashMap[] posHashes = makePositionHashes(alignments);

        //Set up for merging
        double errorRate = 0.01;
        myLogger.info("Merging genotype calls with assumed error rate " + errorRate);
        BasicGenotypeMergeRule mergeRule = new BasicGenotypeMergeRule(errorRate);
        GenotypeTableBuilder genoBuilder = GenotypeTableBuilder.getTaxaIncremental(masterPos, mergeRule);

        //Merge actual calls. Dut to setup of GenotypeBuilder, must get all calls for a taxon and then add as a unit
        for (Taxon t : masterTaxa) {
            //System.out.println("Processing taxon" + t.getName());
            byte[] taxonCalls = new byte[masterPos.size()];
            for (int i_pos = 0; i_pos < masterPos.size(); i_pos++) {
                /*if(i_pos % 5000 == 0){    //For debugging
                 System.out.println("\tProcessed" +  i_pos + "sites");
                 }*/
                Position p = masterPos.get(i_pos);

                //Go through each individual GenotypeTable and assemble all the calls
                ArrayList<Byte> callList = new ArrayList(alignments.length);
                for (int i_align = 0; i_align < alignments.length; i_align++) {
                    GenotypeTable a = alignments[i_align];
                    //Check that this alignment actually has this site-taxon combination
                    if (taxaHashes[i_align].containsKey(t) && posHashes[i_align].containsKey(p)) {
                        int taxonnum = (Integer) taxaHashes[i_align].get(t);
                        int sitenum = (Integer) posHashes[i_align].get(p);
                        byte genotype = a.genotype(taxonnum, sitenum);
                        if (genotype != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) { // Skip missing genotypes
                            callList.add(genotype);
                        }
                    }
                }
                Byte[] calls = callList.toArray(new Byte[0]);
                taxonCalls[i_pos] = mergeCalls(calls); // Where to implement merging rules
            }
            genoBuilder.addTaxon(t, taxonCalls);
        }

        myLogger.info("Finalizing genotype table with " + masterTaxa.size() + " taxa and " + masterPos.size() + "sites");
        GenotypeTable genos = genoBuilder.build();
        return genos;
    }

    /**
     * Generate a list of unique taxa given an array of alignments
     */
    public static TaxaList generateMasterTaxaList(GenotypeTable[] alignments) {
        myLogger.info("Creating unified taxa list");
        //Make a HashSet of taxa, adding each in to get the final union of them all
        HashSet<Taxon> taxaSet = new HashSet();
        for (GenotypeTable a : alignments) {
            taxaSet.addAll(a.taxa());
        }

        //Sort unique taxa, since not done automatically
        Taxon[] taxaArray = taxaSet.toArray(new Taxon[0]);
        Arrays.sort(taxaArray);

        //Make into TaxaList
        TaxaListBuilder taxaBuilder = new TaxaListBuilder();
        taxaBuilder.addAll(taxaArray);
        return taxaBuilder.build();
    }

    /**
     * Generate a list of unique positions given an array of alignments
     */
    public static PositionList generateMasterPositionList(GenotypeTable[] alignments) {
        myLogger.info("Creating unified position list");
        //Make a HashSet of positions, adding each in to get the final union of them all
        HashSet posSet = new HashSet();
        for (GenotypeTable a : alignments) {
            posSet.addAll(a.positions());
        }
        PositionListBuilder posBuilder = new PositionListBuilder();
        posBuilder.addAll(posSet);
        return posBuilder.build();
    }

    /**
     * Takes an array of GenotypeTables and returns a matching array of HashMaps
     * linking each Taxon to its index
     */
    public static HashMap[] makeTaxaHashes(GenotypeTable[] alignments) {
        HashMap[] taxaHashes = new HashMap[alignments.length];
        for (int i = 0; i < alignments.length; i++) {
            taxaHashes[i] = new HashMap<Taxon, Integer>();
            GenotypeTable a = alignments[i];
            Taxon[] taxaArray = a.taxa().toArray(new Taxon[0]);
            for (int t = 0; t < taxaArray.length; t++) {
                taxaHashes[i].put(taxaArray[t], t);
            }
        }
        return taxaHashes;
    }

    /**
     * Takes an array of GenotypeTables and returns a matching array of HashMaps
     * linking each Position to its index
     */
    public static HashMap[] makePositionHashes(GenotypeTable[] alignments) {
        HashMap[] posHashes = new HashMap[alignments.length];
        for (int i = 0; i < alignments.length; i++) {
            posHashes[i] = new HashMap<Position, Integer>();
            GenotypeTable a = alignments[i];
            Position[] posArray = a.positions().toArray(new Position[0]);
            for (int t = 0; t < posArray.length; t++) {
                posHashes[i].put(posArray[t], t);
            }
        }
        return posHashes;
    }

    /**
     * Take an array of Byte genotype calls and merge TODO: Make this smarter;
     * currently just uses the last call (to be same as in TASSEL4)
     */
    public static byte mergeCalls(Byte[] calls) {
        if (calls.length > 0) {
            return (byte) calls[calls.length - 1];
        }
        return GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
    }
}
