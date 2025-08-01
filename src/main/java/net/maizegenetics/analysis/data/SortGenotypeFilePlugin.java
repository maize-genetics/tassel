package net.maizegenetics.analysis.data;

import java.awt.*;
import java.net.URL;
import javax.swing.*;
import static net.maizegenetics.analysis.data.FileLoadPlugin.FILE_EXT_PLINK_MAP;
import static net.maizegenetics.analysis.data.FileLoadPlugin.FILE_EXT_PLINK_PED;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.io.BuilderFromHapMap;
import net.maizegenetics.dna.snp.io.BuilderFromVCF;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Created by jgw87 on 6/5/14. This plugin takes a Hapmap or VCF genotype file
 * and sorts it according to TASSEL's conventions which rely on the position,
 * locus (chromosome), strand, and SNP name (to facilitate searching).
 */
public class SortGenotypeFilePlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(SortGenotypeFilePlugin.class);

    public static enum SupportedFileTypes {
        Hapmap, VCF, Plink
    }

    private PluginParameter<String> infile
            = new PluginParameter.Builder<>("inputFile", null, String.class)
            .required(true)
            .inFile()
            .guiName("Input file")
            .description("Input file")
            .build();
    private PluginParameter<String> outfile
            = new PluginParameter.Builder<>("outputFile", null, String.class)
            .required(true)
            .outFile()
            .guiName("Output file")
            .description("Output file")
            .build();
    private PluginParameter<SupportedFileTypes> fileType
            = new PluginParameter.Builder<>("fileType", null, SupportedFileTypes.class)
            .required(false)
            .guiName("File type")
            .description("Input/output file type (if not obvious from file name)")
            .build();

    public SortGenotypeFilePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        switch (fileType()) {
            case Hapmap:
                GenotypeTable myHapmap = BuilderFromHapMap.getBuilder(inputFile(), null).sortPositions().build();
                ExportUtils.writeToHapmap(myHapmap, outputFile());
                break;
            case VCF:
                GenotypeTable myVCF = BuilderFromVCF.getBuilder(inputFile()).keepDepth().buildAndSortInMemory();
                ExportUtils.writeToVCF(myVCF, outputFile(), true);
                break;
            case Plink:
                if (inputFile().endsWith(FILE_EXT_PLINK_PED) || inputFile().endsWith(FILE_EXT_PLINK_PED + ".gz")) {
                    String theMapFile = inputFile().replaceFirst(FILE_EXT_PLINK_PED, FILE_EXT_PLINK_MAP);
                    GenotypeTable plink = ImportUtils.readFromPLink(inputFile(), theMapFile, this, true);
                    ExportUtils.writeToPlink(plink, outputFile(), '\t');
                } else if (inputFile().endsWith(FILE_EXT_PLINK_MAP) || inputFile().endsWith(FILE_EXT_PLINK_MAP + ".gz")) {
                    String thePedFile = inputFile().replaceFirst(FILE_EXT_PLINK_MAP, FILE_EXT_PLINK_PED);
                    GenotypeTable plink = ImportUtils.readFromPLink(thePedFile, inputFile(), this, true);
                    ExportUtils.writeToPlink(plink, outputFile(), '\t');
                }
                break;
            default:
                throw new IllegalArgumentException("SortGenotypeFilePlugin: Identified data type does not conform to known types (Hapmap, VCF)");
        }

        return null;
    }

    @Override
    public String pluginDescription() {
        return "This plugin takes a Hapmap, VCF, or Plink genotype file and sorts it according to TASSEL's conventions, "
                + "which rely on the position, locus (chromosome), physical position, and SNP name (to facilitate searching).";
    }

    @Override
    protected void postProcessParameters() {
        // If file type not provided, try to guess from file name
        if (fileType() == null) {
            if (inputFile().toLowerCase().endsWith(".hmp.txt") || inputFile().toLowerCase().endsWith(".hmp.txt.gz")) {
                fileType(SupportedFileTypes.Hapmap);
            } else if (inputFile().toLowerCase().endsWith(".vcf") || inputFile().toLowerCase().endsWith(".vcf.gz")) {
                fileType(SupportedFileTypes.VCF);
            } else if (inputFile().endsWith(FILE_EXT_PLINK_PED) || inputFile().endsWith(FILE_EXT_PLINK_PED + ".gz")
                    || inputFile().endsWith(FILE_EXT_PLINK_MAP) || inputFile().endsWith(FILE_EXT_PLINK_MAP + ".gz")) {
                fileType(SupportedFileTypes.Plink);
            } else {
                throw new UnsupportedOperationException("Unable to guess file type from input file name. Please rename to end in .hmp.txt, .hmp.txt.gz, .vcf, .vcf.gz, FILE_EXT_PLINK_PED, or FILE_EXT_PLINK_MAP");
            }
        }
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = SortGenotypeFilePlugin.class.getResource("/net/maizegenetics/analysis/images/sort.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Sort Genotype File";
    }

    @Override
    public String getToolTipText() {
        return "Sort Genotype File by Positions";
    }

    /**
     * Input file
     *
     * @return Input file
     */
    public String inputFile() {
        return infile.value();
    }

    /**
     * Set Input file. Input file
     *
     * @param value Input file
     * @return this plugin
     */
    public SortGenotypeFilePlugin inputFile(String value) {
        infile = new PluginParameter<>(infile, value);
        return this;
    }

    /**
     * Output file
     *
     * @return Output file
     */
    public String outputFile() {
        return outfile.value();
    }

    /**
     * Set Output file. Output file
     *
     * @param value Output file
     * @return this plugin
     */
    public SortGenotypeFilePlugin outputFile(String value) {
        outfile = new PluginParameter<>(outfile, value);
        return this;
    }

    /**
     * Input/output file type (if not obvious from file name)
     *
     * @return File type
     */
    public SupportedFileTypes fileType() {
        return fileType.value();
    }

    /**
     * Set File type. Input/output file type (if not obvious from file name)
     *
     * @param value File type
     * @return this plugin
     */
    public SortGenotypeFilePlugin fileType(SupportedFileTypes value) {
        fileType = new PluginParameter<>(fileType, value);
        return this;
    }
}
