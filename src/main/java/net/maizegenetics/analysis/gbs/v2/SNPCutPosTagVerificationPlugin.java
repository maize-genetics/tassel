/**
 * 
 */
package net.maizegenetics.analysis.gbs.v2;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TagDataWriter;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.google.common.collect.Multimap;

/**
 * This class allows a user to specify a Cut or SNP position for which they
 * would like data printed.  For a Cut Position, the tags associated with that
 * position are printed along with the number of times it appears in each taxa.
 * 
 * For a SNP Position, each allele and the tags associated with that allele are
 * printed along with the number of times the tag appears in each taxa.  The tag
 * is shown both as it is stored in the db, and as a forward strand.  The SNP alignments
 * are based on forward strand.
 * 
 * @author lcj34
 *
 */
public class SNPCutPosTagVerificationPlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(UpdateSNPPositionQualityPlugin.class);

    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with SNP positions stored").build();
    private PluginParameter<String> myChrom = new PluginParameter.Builder<String>("chr", null, String.class).guiName("Chromosome").required(true)
            .description("Chromsome containing the positions").build();
    private PluginParameter<Integer> myPosition = new PluginParameter.Builder<Integer>("pos", null, Integer.class).guiName("Cut or SNP Position").required(true)
            .description("A cut or SNP position number").build();
    private PluginParameter<Byte> myStrand = new PluginParameter.Builder<Byte>("strand", null, Byte.class).guiName("Strand").required(true)
            .description("The strand - 0 for reverse, 1 for forward").build();
    private PluginParameter<String> myPositionType = new PluginParameter.Builder<String>("type", null, String.class).guiName("Type of Position").required(true)
            .description("Type of Position - either snp or cut - for which the TaxaDistribution will be presented").build();
    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<String>("outFile", null, String.class).guiName("Output file").required(true).outFile()
            .description("File name to which tab-delimited output will be written").build();

    private TagDataWriter tdw = null;

    public SNPCutPosTagVerificationPlugin() {
        super(null, false);
    }

    public SNPCutPosTagVerificationPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public SNPCutPosTagVerificationPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        // Need tag depths for each tag that aligns to the specified position
        // Cut position gives us position where tag aligns. 
        // SNP position gives us position within tag where alleles differ
        tdw = new TagDataSQLite(inputDB()); 
        Map<Tag, TaxaDistribution> cutPositionMap = null;
        Multimap<Allele, Map<Tag, TaxaDistribution>> snpPositionMap = null;
        Map<Tag, Position> tagCutPosMap = null;
        try {
            Chromosome myChr = new Chromosome(chrom());
            Position pos = new GeneralPosition.Builder(myChr, cutOrSnpPosition()).strand(strand()).build();      	        	
            TaxaList taxaList = tdw.getTaxaList(); // Used for printing taxon column headers
            if (positionType().equals("cut")) {
                // Get tag/taxon map, print to tab-delimited file
                try {
                    
                } catch (Exception exc) {
                    System.out.println("\nError attempting to grab Cut position data for chrom " + chrom() + ", pos " + cutOrSnpPosition());
                    System.out.println("Please verify, via a database manager tool or other method, that the position requested is stored as a Cut position in the input database.\n");
                    return null;
                }
                cutPositionMap = tdw.getTagsTaxaMap(pos);
                writeCutPositionTagTaxonFile(taxaList, cutPositionMap);
            } else if (positionType().equals("snp")) {
                // create map with alleles, tag and taxon, print to tab-delimited file
                try {
                    snpPositionMap = tdw.getAllelesTagTaxaDistForSNP(pos);
                } catch (Exception exc) {
                    System.out.println("\nError attempting to grab SNP position data for chrom " + chrom() + ", pos " + cutOrSnpPosition());
                    System.out.println("Please verify the position is stored as a SNP position in the input database.");
                    System.out.println("Use a database manager tool to check your db, or run the SNPQualityProfilerPlugin to check for SNP positions.\n");
                    return null;
                }
                
                //            	snpPositionMap.entries().forEach( entry-> {
                //            		System.out.println("LCJ - Allele call is: " + entry.getKey().alleleAsString());
                //            	});
                // get list of tags, send to db to get cut position/strand
                Set<Tag> fullTagList = new HashSet<Tag>();
                snpPositionMap.entries().forEach(entry -> {
                    Set<Tag> tags= entry.getValue().keySet();
                    fullTagList.addAll(tags);
                });
                tagCutPosMap = tdw.getTagCutPosition(fullTagList);
                writeSNPPositionTagTaxonFile(taxaList, snpPositionMap, tagCutPosMap);
            } else {
                myLogger.error("Position type must be specified as either snp or cut\n");
                return null;
            }       	
            ((TagDataSQLite)tdw).close();  
            myLogger.info("SNPCutPosTagVerificationPlugin: Finished writing TaxaDistribution to file for position " + positionType() + ".\n");
        } catch (Exception exc) {
            myLogger.error("SNPCutPosTagVerificationPlugin: caught error " + exc);
            exc.printStackTrace();
        }
        return null;
    }

    private void writeCutPositionTagTaxonFile(TaxaList taxaList, Map<Tag, TaxaDistribution> cutPositionMap) throws Exception{
        BufferedWriter fileWriter = null;
        StringBuilder strB = new StringBuilder();
        if(outputFile()!=null) {
            // taxanumber from TaxaDistribution is in the depths - they are ordered
            // by the taxalist numbers.  Is the TaxaList order alphabetically ???
            // first write the headers, which is a list of the taxa
            strB.append("Chr\tPos\tTag");
            taxaList.stream().forEach(item -> { // column names are the taxon names
                strB.append("\t");
                strB.append(item.getName());
            });
            strB.append("\n");

            cutPositionMap.entrySet().stream().forEach(entry -> {  
                strB.append(chrom());
                strB.append("\t");
                strB.append(cutOrSnpPosition());
                strB.append("\t");
                Tag curTag = entry.getKey();
                strB.append(curTag.sequence()); // add tag sequence in first column

                // This is CUT position - no ALLELEs here
                TaxaDistribution tagTD = entry.getValue();
                int[] depths = tagTD.depths(); // gives us the depths for each taxon
                for (int idx = 0; idx < depths.length; idx++) {
                    strB.append("\t"); 
                    strB.append(depths[idx]);  // add tag depth         		
                }
                strB.append("\n"); // end of line - start next tag           	
            });
            try {  
                fileWriter = new BufferedWriter(new FileWriter(outputFile()));
                fileWriter.write(strB.toString());
            }
            catch(IOException e) {
                myLogger.error("Caught Exception in writeCutPositionTagTaxonFile");
                System.out.println(e);
            }
            fileWriter.close();
        } else {
            myLogger.warn("Outputfile is null - nothing happening here");
        }
    }

    private void writeSNPPositionTagTaxonFile(TaxaList taxaList, 
            Multimap<Allele, Map<Tag, TaxaDistribution>> snpPositionMap, Map<Tag, Position>tagPosMap) throws Exception{
        BufferedWriter fileWriter = null;
        StringBuilder strB = new StringBuilder();
        if(outputFile()!=null) {
            // taxanumber from TaxaDistribution is in the depths - they are ordered
            // by the taxalist numbers.  Is the TaxaList order alphabetically ???
            strB.append("Chr\tSNPPos\tAllele\tTag\tForwardStrand\tTagAsForwardStrand\tCutPos:SNPOffset");
           // strB.append("Chr\tSNPPos\tAllele\tTag\tForwardStrand\tTagAsForwardStrand\tCutPos-SNPOffset"); // first column, ie row header
            taxaList.stream().forEach(item -> { // column names are the taxon names
                strB.append("\t");
                strB.append(item.getName());
            });
            strB.append("\n");

            snpPositionMap.entries().stream().forEach(entry -> {
                Allele curAllele = entry.getKey();
                strB.append(chrom());
                strB.append("\t");
                strB.append(cutOrSnpPosition());
                strB.append("\t");
                strB.append(NucleotideAlignmentConstants.getHaplotypeNucleotide(curAllele.allele()));
                strB.append("\t");
                Map<Tag, TaxaDistribution> curTagTaxa = entry.getValue();
                // Loop through the tag/taxa for this tag. 
                for ( Map.Entry<Tag, TaxaDistribution> tagTaxaMap : curTagTaxa.entrySet()) {
                    Tag curTag = tagTaxaMap.getKey();
                    Position cutPos = tagPosMap.get(curTag);
                    TaxaDistribution tagTD = tagTaxaMap.getValue();
                    strB.append(curTag.sequence()); 
                    strB.append("\t");
                    boolean isForward = cutPos.getAnnotation().getTextAnnotation("forward")[0].equals("true") ? true: false;
                    strB.append(cutPos.getAnnotation().getTextAnnotation("forward")[0]);
                    strB.append("\t");
                    if (isForward) {
                        strB.append(curTag.sequence());
                    } else { // alignments are based on forward strand, create and add for easier SNP verification
                        strB.append(curTag.toReverseComplement());
                    }
                    strB.append("\t");
                    strB.append(cutPos.getPosition());
                    //strB.append("-"); 
                    strB.append(":");
                    //int offSetVal = isForward ? (cutOrSnpPosition() - cutPos.getPosition()) : ()
                    int offSetVal = cutPos.getPosition() - cutOrSnpPosition(); // why is this a minus ??
                    System.out.println("LCJ: writeSNPPosition: cutPos: " + cutPos.getPosition()
                    + ", cutOrSNPPosition: " + cutOrSnpPosition() + " offsetVal: " + offSetVal + " strand:" + isForward);
                    strB.append(offSetVal);

                    int[] depths = tagTD.depths(); // gives us the depths for each taxon
                    for (int idx = 0; idx < depths.length; idx++) {
                        // write the tag depth to each column
                        strB.append("\t");
                        strB.append(depths[idx]);           		
                    }
                    strB.append("\n"); // end of line - start next tag           	
                }
                // new line already added to file
            });
            try {  
                fileWriter = new BufferedWriter(new FileWriter(outputFile()));
                fileWriter.write(strB.toString());
            }
            catch(IOException e) {
                myLogger.error("Caught exception in writeSNPPositionTagTaxonFile");
                System.out.println(e);
            }
            fileWriter.close();
        } else {
            myLogger.warn("Outputfile is null - nothing happening here");
        }
    }


    @Override
    public String getToolTipText() {
        return "Debug tool: Verify which Tags in which taxon map to the specified cut position.  Verify which tags have a SNP at the specified position.";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "SNP/Cut Position Verification";
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public String runPlugin(DataSet input) {
        return (String) performFunction(input).getData(0).getData();
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
    public SNPCutPosTagVerificationPlugin inputDB(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }

    /**
     * Chromosome as String
     *
     * @return Chromosome name
     */
    public String chrom() {
        return myChrom.value();
    }

    /**
     * Set outputDir path.  Directory file path where
     * output files will be writen
     *
     * @param value Directory path
     *
     * @return this plugin
     */
    public SNPCutPosTagVerificationPlugin chrom(String value) {
        myChrom = new PluginParameter<>(myChrom, value);
        return this;
    }

    /**
     * Cut Position
     *
     * @return Cut position
     */
    public Integer cutOrSnpPosition() {
        return myPosition.value();
    }

    /**
     * Set cut position. 
     *
     * @param value cut position
     *
     * @return this plugin
     */
    public SNPCutPosTagVerificationPlugin cutOrSnpPosition(Integer value) {
        myPosition = new PluginParameter<>(myPosition, value);
        return this;
    }

    /**
     * Strand direction
     *
     * @return Strand direction
     */
    public Byte strand() {
        return myStrand.value();
    }

    /**
     * Set strand direction. 
     *
     * @param value strand direction
     *
     * @return this plugin
     */
    public SNPCutPosTagVerificationPlugin strand(Byte value) {
        myStrand = new PluginParameter<>(myStrand, value);
        return this;
    }
    /**
     * SNP Position
     *
     * @return SNP position
     */
    public String positionType() {
        return myPositionType.value();
    }

    /**
     * Set SNP position. 
     *
     * @param value SNP position
     *
     * @return this plugin
     */
    public SNPCutPosTagVerificationPlugin positionType(String value) {
        myPositionType = new PluginParameter<>(myPositionType, value);
        return this;
    }
    /**
     * output directory path
     *
     * @return outPutDir directory string
     */
    public String outputFile() {
        return myOutputFile.value();
    }

    /**
     * Set outputDir path.  Directory file path where
     * output files will be writen
     *
     * @param value Directory path
     *
     * @return this plugin
     */
    public SNPCutPosTagVerificationPlugin outputFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
        return this;
    }
}
