/*
 * ConvertAlignmentCoordinatesPlugin
 */
package net.maizegenetics.analysis.data;

import java.awt.*;
import java.io.BufferedReader;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Pattern;
import javax.swing.*;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class ConvertAlignmentCoordinatesPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(ConvertAlignmentCoordinatesPlugin.class);

    private PluginParameter<String> myMapFilename = new PluginParameter.Builder<>("mapFile", null, String.class)
            .description("")
            .inFile()
            .required(true)
            .build();

    private PluginParameter<String> myGenomeVersion = new PluginParameter.Builder<>("genomeVersion", "AGPv3", String.class)
            .description("")
            .build();

    private final HashMap<String, Chromosome> myLociMap = new HashMap<>();
    private final HashMap<String, Chromosome> myAlignmentLociMap = new HashMap<>();

    public ConvertAlignmentCoordinatesPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);

        if (alignInList.size() != 1) {
            throw new IllegalArgumentException("Invalid selection.  Please select one genotype table.");
        }
        Datum datum = alignInList.get(0);
        GenotypeTable alignment = (GenotypeTable) datum.getData();

        Chromosome[] loci = alignment.chromosomes();
        myLociMap.clear();
        for (int i = 0; i < loci.length; i++) {
            myLociMap.put(loci[i].getName(), loci[i]);
            myAlignmentLociMap.put(loci[i].getName(), loci[i]);
        }

        int numSites = alignment.numberOfSites();
        String[] snpIDs = new String[numSites];
        for (int i = 0; i < numSites; i++) {
            snpIDs[i] = alignment.siteName(i);
        }

        int count = 1;
        PositionListBuilder posBuilder = new PositionListBuilder().addAll(alignment.positions()).genomeVersion(genomeVersion());
        try (BufferedReader br = Utils.getBufferedReader(mapFilename())) {

            Pattern sep = Pattern.compile("\\s+");

            // Ignore head, then get first line.
            String inputline = br.readLine();
            inputline = br.readLine();
            int numChanges = 0;
            while (inputline != null) {

                count++;
                inputline = inputline.trim();
                // ID	v1_chr	v1_pos	v2_chr	v2_pos	strand
                String[] parsedline = sep.split(inputline);
                inputline = br.readLine();

                String locus1 = getLocusName(parsedline[1]);

                if (myAlignmentLociMap.get(locus1) != null) {

                    String snpID = parsedline[0];
                    int pos1 = Integer.valueOf(parsedline[2]);
                    String locus2 = getLocusName(parsedline[3]);
                    int pos2 = Integer.valueOf(parsedline[4]);

                    if ((!locus1.equals(locus2)) || (pos1 != pos2)) {

                        int site = getSiteOfSNPID(snpID, snpIDs);
                        if (site < 0) {
                            continue;
                        }

                        if ((pos1 != alignment.chromosomalPosition(site)) || (!locus1.equals(alignment.chromosome(site).getName()))) {
                            myLogger.warn("map file line: " + count + "  SNP ID: " + snpID + "  position: " + pos1 + "  locus: " + locus1 + " position and locus do not match alignment.");
                            myLogger.warn("Alignment SNP ID: " + alignment.siteName(site) + "  position: " + alignment.chromosomalPosition(site) + "  locus: " + alignment.chromosomeName(site));
                            continue;
                        }

                        numChanges++;
                        GeneralPosition.Builder newPos = new GeneralPosition.Builder(alignment.positions().get(site));
                        newPos.chromosome(getLocusObj(locus2)).position(pos2).snpName(snpID);
                        posBuilder.set(site, newPos.build());

                    }

                }

            }

            myLogger.info("Number Changes: " + numChanges);

            Tuple<PositionList, int[]> temp = posBuilder.buildWithSiteRedirect();
            GenotypeTable newGenotypeTable = FilterGenotypeTable.getInstance(alignment, temp.x, temp.y);

            String comment = alignment.positions().genomeVersion() + " converted to " + genomeVersion();
            Datum result = new Datum(datum.getName() + "_NewCoordinates", newGenotypeTable, comment);

            return new DataSet(result, this);

        } catch (Exception e) {
            throw new IllegalStateException("ConvertAlignmentCoordinatesPlugin: processDatum: problem converting alignment: line: " + count + "  message: " + e.getMessage());
        }

    }

    private Chromosome getLocusObj(String locus) {
        Chromosome locusObj = myLociMap.get(locus);
        if (locusObj == null) {
            locusObj = new Chromosome(locus);
            myLociMap.put(locus, locusObj);
        }
        return locusObj;
    }

    private String getLocusName(String input) {
        if (input.startsWith("chr")) {
            return input.substring(3);
        } else {
            return input;
        }
    }

    private int getSiteOfSNPID(String id, String[] snpIDs) {
        for (int i = 0, n = snpIDs.length; i < n; i++) {
            if (id.equals(snpIDs[i])) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Map File
     *
     * @return Map File
     */
    public String mapFilename() {
        return myMapFilename.value();
    }

    /**
     * Set Map File. Map File
     *
     * @param value Map File
     *
     * @return this plugin
     */
    public ConvertAlignmentCoordinatesPlugin mapFilename(String value) {
        myMapFilename = new PluginParameter<>(myMapFilename, value);
        return this;
    }

    /**
     * Genome Version
     *
     * @return Genome Version
     */
    public String genomeVersion() {
        return myGenomeVersion.value();
    }

    /**
     * Set Genome Version. Genome Version
     *
     * @param value Genome Version
     *
     * @return this plugin
     */
    public ConvertAlignmentCoordinatesPlugin genomeVersion(String value) {
        myGenomeVersion = new PluginParameter<>(myGenomeVersion, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Convert Genotype Table Coordinates";
    }

    @Override
    public String getToolTipText() {
        return "Convert Genotype Table Coordinates";
    }
}
