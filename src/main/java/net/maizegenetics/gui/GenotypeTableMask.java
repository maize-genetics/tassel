/*
 * GenotypeTableMask
 */
package net.maizegenetics.gui;

import java.awt.Color;

import java.io.Serializable;
import net.maizegenetics.dna.snp.GenotypeTable;

/**
 * @author Terry Casstevens
 */
public interface GenotypeTableMask extends Serializable {

    /**
     * This defines the type of mask.
     *
     * reference - Mask created using taxon as reference.
     * imputed - Mask created to identify imputed values.
     * compared = Mask created to identify differences between two alignments.
     */
    public enum MaskType {

        reference, imputed, compared
    };

    public byte getMask(int taxon, int site);

    public Color getColor();

    public void setColor(Color color);

    public GenotypeTable getAlignment();

    public MaskType getMaskType();

}
