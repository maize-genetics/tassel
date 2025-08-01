/*
 * AlignmentTableCellRenderer
 */
package net.maizegenetics.gui;

import java.awt.Color;
import java.awt.Component;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.table.DefaultTableCellRenderer;
import net.maizegenetics.dna.map.DonorHaplotypes;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.genotypecall.ProjectionGenotypeCallTable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * @author Terry Casstevens
 */
public class AlignmentTableCellRenderer extends DefaultTableCellRenderer {

    private static final Logger myLogger = LogManager.getLogger(AlignmentTableCellRenderer.class);

    private static int[] NUCLEOTIDE_COLORS = new int[16];

    static {
        NUCLEOTIDE_COLORS[NucleotideAlignmentConstants.A_ALLELE] = 0xFA0000;
        NUCLEOTIDE_COLORS[NucleotideAlignmentConstants.C_ALLELE] = 0x9BCD9B;
        NUCLEOTIDE_COLORS[NucleotideAlignmentConstants.G_ALLELE] = 0x4876FF;
        NUCLEOTIDE_COLORS[NucleotideAlignmentConstants.T_ALLELE] = 0x777D7E;
        NUCLEOTIDE_COLORS[NucleotideAlignmentConstants.GAP_ALLELE] = 0xFFF68F;
        NUCLEOTIDE_COLORS[NucleotideAlignmentConstants.INSERT_ALLELE] = 0xFF8C00;
    }

    private static final Color MAJOR_ALLELE_COLOR = new Color(0xe6cf45);
    private static final Color MINOR_ALLELE_COLOR = new Color(0x45a1e6);
    private static final Color HETEROZYGOUS_COLOR = new Color(0xe64545);
    private static final Color MAJOR_MINOR_ALLELE_COLOR = new Color((MAJOR_ALLELE_COLOR.getRGB() + MINOR_ALLELE_COLOR.getRGB()) % 0xFFFFFF);
    private static final Color[] COLORS_256 = generateColors(256);
    private static final Map<String, Color> COLORS_NUCLEOTIDES = generateNucleotideColors();

    public enum RENDERING_TYPE {

        Nucleotide, NucleotideHeterozygous, MajorAllele, MinorAllele, MajorMinorAllele,
        Heterozygous, ReferenceMasks, GeneticDistanceMasks, Depth, None, TOPM, SNPs,
        NumericGenotype, Projection
    }

    private static final RENDERING_TYPE[] GENOTYPE_RENDERING_TYPES = new RENDERING_TYPE[]{RENDERING_TYPE.MajorMinorAllele,
            RENDERING_TYPE.MajorAllele, RENDERING_TYPE.MinorAllele, RENDERING_TYPE.Heterozygous,
            RENDERING_TYPE.ReferenceMasks, RENDERING_TYPE.GeneticDistanceMasks, RENDERING_TYPE.None};
    protected final AlignmentTableModel myAlignmentTableModel;
    private final GenotypeTable myAlignment;
    private GenotypeTableMask[] myMasks;
    private final RENDERING_TYPE[] mySupportedRenderingTypes;
    private RENDERING_TYPE myRenderingType = RENDERING_TYPE.MajorMinorAllele;
    private final Map<Integer, byte[]> myCachedAlleles = new LinkedHashMap<Integer, byte[]>() {
        protected boolean removeEldestEntry(Map.Entry eldest) {
            return size() > 100;
        }
    };

    public AlignmentTableCellRenderer(AlignmentTableModel model, GenotypeTable alignment, GenotypeTableMask[] masks) {
        myAlignmentTableModel = model;
        myAlignment = alignment;
        myMasks = masks;

        List<RENDERING_TYPE> temp = new ArrayList<>();

        if (myAlignment.hasReferenceProbablity()) {
            temp.add(RENDERING_TYPE.NumericGenotype);
        }

        if (myAlignment.hasGenotype()) {

            try {
                if (myAlignment.genotypeMatrix() instanceof ProjectionGenotypeCallTable) {
                    temp.add(RENDERING_TYPE.Projection);
                }
            } catch (Exception e) {
                // do nothing
            }

            for (int i = 0; i < GENOTYPE_RENDERING_TYPES.length; i++) {
                temp.add(GENOTYPE_RENDERING_TYPES[i]);
            }

            try {
                if (NucleotideAlignmentConstants.isNucleotideEncodings(myAlignment.alleleDefinitions())) {
                    temp.add(RENDERING_TYPE.Nucleotide);
                    temp.add(RENDERING_TYPE.NucleotideHeterozygous);
                }
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
            }
        }

        if (myAlignment.hasDepth()) {
            temp.add(RENDERING_TYPE.Depth);
        }

        mySupportedRenderingTypes = new RENDERING_TYPE[temp.size()];
        for (int i = 0; i < temp.size(); i++) {
            mySupportedRenderingTypes[i] = temp.get(i);
        }

        if (mySupportedRenderingTypes.length != 0) {
            setRenderingType(mySupportedRenderingTypes[0]);
        }

    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        switch (myRenderingType) {
            case Nucleotide:
                return getNucleotideRendering(table, value, isSelected, hasFocus, row, col);
            case NucleotideHeterozygous:
                return getNucleotideHeterozygousRendering(table, value, isSelected, hasFocus, row, col);
            case MajorAllele:
                return getMajorAlleleRendering(table, value, isSelected, hasFocus, row, col);
            case MinorAllele:
                return getMinorAlleleRendering(table, value, isSelected, hasFocus, row, col);
            case Heterozygous:
                return getHeterozygousRendering(table, value, isSelected, hasFocus, row, col);
            case MajorMinorAllele:
                return getMajorMinorAlleleRendering(table, value, isSelected, hasFocus, row, col);
            case ReferenceMasks:
                return getReferenceMasksRendering(table, value, isSelected, hasFocus, row, col);
            case NumericGenotype:
                return getReferenceProbabilityRendering(table, value, isSelected, hasFocus, row, col);
            case Projection:
                return getProjectionRendering(table, value, isSelected, hasFocus, row, col);
            case GeneticDistanceMasks:
                return getGeneticDistanceMasksRendering(table, value, isSelected, hasFocus, row, col);
            case Depth:
                return getDepthMasksRendering(table, value, isSelected, hasFocus, row, col);
            case None:
                return getDefaultRendering(table, value, isSelected, hasFocus, row, col);
            default:
                return getDefaultRendering(table, value, isSelected, hasFocus, row, col);
        }

    }

    protected Component getDefaultRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else {
            comp.setBackground(null);
        }

        return comp;

    }

    protected Component getNucleotideRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        String alleles = myAlignment.genotypeAsString(row, myAlignmentTableModel.getRealColumnIndex(col));

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if (alleles.equals(GenotypeTable.UNKNOWN_ALLELE_STR)) {
            comp.setBackground(null);
        } else {
            comp.setBackground(COLORS_NUCLEOTIDES.get(alleles));
        }

        return comp;

    }

    private Component getNucleotideHeterozygousRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        int site = myAlignmentTableModel.getRealColumnIndex(col);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if (myAlignment.isHeterozygous(row, site)) {
            String alleles = myAlignment.genotypeAsString(row, myAlignmentTableModel.getRealColumnIndex(col));
            comp.setBackground(COLORS_NUCLEOTIDES.get(alleles));
        } else {
            comp.setBackground(null);
        }

        return comp;

    }

    private Component getMajorAlleleRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        int site = myAlignmentTableModel.getRealColumnIndex(col);
        byte[] alleles = myCachedAlleles.get(site);
        if (alleles == null) {
            alleles = myAlignment.alleles(site);
            myCachedAlleles.put(site, alleles);
        }
        byte major = GenotypeTable.UNKNOWN_ALLELE;
        if (alleles.length > 0) {
            major = alleles[0];
        }
        byte[] diploidValues = myAlignment.genotypeArray(row, site);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if ((diploidValues[0] == major) || (diploidValues[1] == major)) {
            comp.setBackground(MAJOR_ALLELE_COLOR);
        } else {
            comp.setBackground(null);
        }

        return comp;

    }

    private Component getHeterozygousRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        int site = myAlignmentTableModel.getRealColumnIndex(col);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if (myAlignment.isHeterozygous(row, site)) {
            comp.setBackground(HETEROZYGOUS_COLOR);
        } else {
            comp.setBackground(null);
        }

        return comp;

    }

    private Component getMinorAlleleRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        int site = myAlignmentTableModel.getRealColumnIndex(col);
        byte[] alleles = myCachedAlleles.get(site);
        if (alleles == null) {
            alleles = myAlignment.alleles(site);
            myCachedAlleles.put(site, alleles);
        }
        byte minor = GenotypeTable.UNKNOWN_ALLELE;
        if (alleles.length > 1) {
            minor = alleles[1];
        }
        byte[] diploidValues = myAlignment.genotypeArray(row, site);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if ((diploidValues[0] == minor) || (diploidValues[1] == minor)) {
            comp.setBackground(MINOR_ALLELE_COLOR);
        } else {
            comp.setBackground(null);
        }

        return comp;

    }

    private Component getMajorMinorAlleleRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        int site = myAlignmentTableModel.getRealColumnIndex(col);
        byte[] alleles = myCachedAlleles.get(site);
        if (alleles == null) {
            alleles = myAlignment.alleles(site);
            myCachedAlleles.put(site, alleles);
        }

        byte[] diploidValues = myAlignment.genotypeArray(row, site);
        if (alleles.length > 1) {
            byte major = alleles[0];
            byte minor = alleles[1];
            if (isSelected) {
                comp.setBackground(Color.DARK_GRAY);
            } else if (((diploidValues[0] == major) && (diploidValues[1] == minor))
                    || ((diploidValues[0] == minor) && (diploidValues[1] == major))) {
                comp.setBackground(MAJOR_MINOR_ALLELE_COLOR);
            } else if ((diploidValues[0] == major) || (diploidValues[1] == major)) {
                comp.setBackground(MAJOR_ALLELE_COLOR);
            } else if ((diploidValues[0] == minor) || (diploidValues[1] == minor)) {
                comp.setBackground(MINOR_ALLELE_COLOR);
            } else {
                comp.setBackground(null);
            }
        } else if (alleles.length == 1) {
            byte major = alleles[0];
            if (isSelected) {
                comp.setBackground(Color.DARK_GRAY);
            } else if ((diploidValues[0] == major) || (diploidValues[1] == major)) {
                comp.setBackground(MAJOR_ALLELE_COLOR);
            } else {
                comp.setBackground(null);
            }
        } else {
            if (isSelected) {
                comp.setBackground(Color.DARK_GRAY);
            } else {
                comp.setBackground(null);
            }
        }

        return comp;

    }

    public Component getReferenceMasksRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        GenotypeTableMask[] masks = getAlignmentMasksOfClass(GenotypeTableMaskReference.class);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if ((masks == null) || (masks.length == 0)) {
            comp.setBackground(null);
        } else if (masks.length == 1) {

            if (masks[0].getMask(row, myAlignmentTableModel.getRealColumnIndex(col)) == 0x1) {
                comp.setBackground(masks[0].getColor());
            } else {
                comp.setBackground(null);
            }

        } else {

            int red = 0;
            int green = 0;
            int blue = 0;

            boolean changed = false;
            for (int i = 0; i < masks.length; i++) {
                if (masks[i].getMask(row, myAlignmentTableModel.getRealColumnIndex(col)) == 0x1) {
                    red = red + masks[i].getColor().getRed();
                    green = green + masks[i].getColor().getGreen();
                    blue = blue + masks[i].getColor().getBlue();
                    changed = true;
                }
            }

            if (changed) {
                red = red % 256;
                green = green % 256;
                blue = blue % 256;
                comp.setBackground(new Color(red, green, blue));
            } else {
                comp.setBackground(null);
            }

        }

        return comp;

    }

    private static final int[] PROJECTION_MASKS = new int[]{0xFF, 0xFF00, 0xFF0000};

    private Component getProjectionRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        int site = myAlignmentTableModel.getRealColumnIndex(col);

        Optional<DonorHaplotypes> optional;
        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if ((optional = donorHaplotype(row, site)).isPresent()) {
            DonorHaplotypes dh = optional.get();
            int color = PROJECTION_MASKS[row % 3];
            color = color | dh.getStartPosition() | dh.getEndPosition();
            comp.setBackground(new Color(color));
        } else {
            comp.setBackground(null);
        }

        return comp;

    }

    private Optional<DonorHaplotypes> donorHaplotype(int taxon, int site) {
        Position position = myAlignment.positions().get(site);
        int physical = position.getPosition();
        return ((ProjectionGenotypeCallTable) myAlignment.genotypeMatrix()).getDonorHaplotypes(taxon).stream()
                .filter(dh -> dh.getChromosome().equals(position.getChromosome()))
                .filter(dh -> (physical >= dh.getStartPosition() && physical <= dh.getEndPosition()))
                .findFirst();
    }

    public Component getReferenceProbabilityRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else {
            int site = myAlignmentTableModel.getRealColumnIndex(col);
            float probability = myAlignment.referenceProbability().value(row, site);
            if (Float.isNaN(probability)) {
                comp.setBackground(Color.GRAY);
            } else {
                comp.setBackground(COLORS_256[255 - Math.round(probability * 255.0f)]);
            }
        }

        return comp;

    }

    public Component getGeneticDistanceMasksRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        GenotypeTableMask[] masks = getAlignmentMasksOfClass(GenotypeTableMaskGeneticDistance.class);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if ((masks == null) || (masks.length == 0)) {
            comp.setBackground(null);
        } else {

            int aveDistance = 0;
            for (int i = 0; i < masks.length; i++) {
                aveDistance += 0xFF & masks[i].getMask(row, myAlignmentTableModel.getRealColumnIndex(col));
            }
            aveDistance /= masks.length;

            comp.setBackground(COLORS_256[aveDistance]);

        }

        return comp;

    }

    public Component getDepthMasksRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if (!myAlignment.hasDepth()) {
            comp.setBackground(null);
        } else {

            int depth = 0;
            int site = myAlignmentTableModel.getRealColumnIndex(col);
            byte[] diploidValues = myAlignment.genotypeArray(row, site);
            int[] depths = myAlignment.depthForAlleles(row, site);
            if (diploidValues[0] < 6) {
                depth = depths[diploidValues[0]];
            }
            if ((diploidValues[1] < 6) && (diploidValues[1] != diploidValues[0])) {
                depth += depths[diploidValues[1]];
            }

            if (depth == 0) {
                comp.setBackground(null);
            } else {
                //depth = (int) Math.rint((double) depth / (double) count * 100.0);
                depth *= 30;
                if (depth > 255) {
                    depth = 255;
                }
                comp.setBackground(COLORS_256[255 - depth]);
            }

        }

        return comp;

    }

    private static Color[] generateColors(int n) {
        Color[] cols = new Color[n];
        for (int i = 0; i < n; i++) {
            cols[i] = Color.getHSBColor(((float) i / (float) n * 0.6f), 0.85f, 0.9f);
        }
        return cols;
    }

    private static Map<String, Color> generateNucleotideColors() {
        Map<String, Color> result = new HashMap<>();
        Iterator itr = NucleotideAlignmentConstants.NUCLEOTIDE_IUPAC_HASH.entrySet().iterator();
        while (itr.hasNext()) {
            Map.Entry current = (Map.Entry) itr.next();
            result.put((String) current.getValue(), Color.BLACK);
        }
        int numCodes = result.size();
        itr = result.entrySet().iterator();
        int count = 0;
        while (itr.hasNext()) {
            Map.Entry current = (Map.Entry) itr.next();
            Color color = Color.getHSBColor((float) count / (float) (numCodes - 1), 0.7f, 0.9f);
            result.put((String) current.getKey(), color);
            count++;
        }
        return result;
    }

    private GenotypeTableMask[] getAlignmentMasksOfClass(Class type) {
        if ((myMasks == null) || (myMasks.length == 0)) {
            return null;
        }
        List<GenotypeTableMask> masks = new ArrayList<>();
        for (int i = 0; i < myMasks.length; i++) {
            if (type.isInstance(myMasks[i])) {
                masks.add(myMasks[i]);
            }
        }
        GenotypeTableMask[] result = new GenotypeTableMask[masks.size()];
        masks.toArray(result);
        return result;
    }

    public RENDERING_TYPE getRenderingType() {
        return myRenderingType;
    }

    final public void setRenderingType(RENDERING_TYPE type) {
        myRenderingType = type;
        myAlignmentTableModel.setRenderingType(type);
    }

    public void setMasks(GenotypeTableMask[] masks) {
        myMasks = masks;
    }

    public RENDERING_TYPE[] getRenderingTypes() {
        return mySupportedRenderingTypes;
    }

}
