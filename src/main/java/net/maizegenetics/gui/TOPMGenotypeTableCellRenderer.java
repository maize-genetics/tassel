/*
 *  TOPMGenotypeTableCellRenderer
 */
package net.maizegenetics.gui;

import java.awt.Color;
import java.awt.Component;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import net.maizegenetics.dna.snp.TOPMGenotypeTable;

/**
 *
 * @author Terry Casstevens
 */
public class TOPMGenotypeTableCellRenderer extends AlignmentTableCellRenderer {

    private static RENDERING_TYPE[] SUPPORTED_RENDERING_TYPES = new RENDERING_TYPE[]{RENDERING_TYPE.Nucleotide,
        RENDERING_TYPE.None, RENDERING_TYPE.TOPM, RENDERING_TYPE.SNPs};

    private final TOPMGenotypeTable myGenotypeTable;

    public TOPMGenotypeTableCellRenderer(AlignmentTableModel model, TOPMGenotypeTable genotypeTable) {
        super(model, genotypeTable, null);
        myGenotypeTable = genotypeTable;
        setRenderingType(RENDERING_TYPE.TOPM);
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        switch (getRenderingType()) {
            case TOPM:
                return getTOPMRendering(table, value, isSelected, hasFocus, row, col);
            case SNPs:
                return getSNPsRendering(table, value, isSelected, hasFocus, row, col);
            default:
                return super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);
        }

    }

    private Component getTOPMRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        setHorizontalAlignment(SwingConstants.CENTER);

        if (isSelected) {
            Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);
            comp.setBackground(Color.DARK_GRAY);
            return comp;
        } else if (myGenotypeTable.isSNP(row, myAlignmentTableModel.getRealColumnIndex(col))) {
            Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);
            comp.setBackground(Color.GREEN);
            return comp;
        } else {
            return super.getNucleotideRendering(table, value, isSelected, hasFocus, row, col);
        }

    }
    
    private Component getSNPsRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        setHorizontalAlignment(SwingConstants.CENTER);

        if (isSelected) {
            Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);
            comp.setBackground(Color.DARK_GRAY);
            return comp;
        } else if (myGenotypeTable.isSNP(row, myAlignmentTableModel.getRealColumnIndex(col))) {
            Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);
            comp.setBackground(Color.GREEN);
            return comp;
        } else {
            return super.getDefaultRendering(table, value, isSelected, hasFocus, row, col);
        }

    }

    public AlignmentTableCellRenderer.RENDERING_TYPE[] getRenderingTypes() {
        return SUPPORTED_RENDERING_TYPES;
    }

}
