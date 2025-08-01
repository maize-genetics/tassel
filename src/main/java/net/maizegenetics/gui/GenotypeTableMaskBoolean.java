/*
 * GenotypeTableMaskBoolean
 */
package net.maizegenetics.gui;

import net.maizegenetics.taxa.Taxon;

import java.awt.*;

import net.maizegenetics.dna.snp.GenotypeTable;

/**
 *
 * @author Terry Casstevens
 */
public class GenotypeTableMaskBoolean extends AbstractGenotypeTableMask {

    private static final long serialVersionUID = -5197800047652332969L;
    private final byte[][] myMask;

    public GenotypeTableMaskBoolean(GenotypeTable align, byte[][] mask, String name, MaskType type) {
        this(align, mask, name, getNextColor(), type);
    }

    public GenotypeTableMaskBoolean(GenotypeTable align, byte[][] mask, String name, Color color, MaskType type) {

        super(align, name, color, type);

        if (mask.length != align.numberOfTaxa()) {
            throw new IllegalArgumentException("GenotypeTableMaskBoolean: init: number of mask rows should equal number of sequences.");
        }

        int numBytesNeeded = getNumMaskColumns(align.numberOfSites());
        if (numBytesNeeded != mask[0].length) {
            throw new IllegalArgumentException("GenotypeTableMaskBoolean: init: incorrect number of mask columns: " + mask[0].length + "  should be: " + numBytesNeeded);
        }

        myMask = mask;

    }

    public static GenotypeTableMaskBoolean getInstanceCompareReference(GenotypeTable align) {
        String name = align.taxaName(0) + " Reference";
        return getInstanceCompareReference(align, align.genotypeRange(0, 0, align.numberOfSites()), name);
    }

    public static GenotypeTableMaskBoolean getInstanceCompareReference(GenotypeTable align, Taxon id) {
        int index = align.taxa().indexOf(id);
        if (index < 0) {
            throw new IllegalArgumentException("GenotypeTableMaskBoolean: getInstanceCompareReference: unknown id: " + id);
        }
        String name = id.getName() + " Reference";
        return getInstanceCompareReference(align, align.genotypeRange(index, 0, align.numberOfSites()), name);
    }

    public static GenotypeTableMaskBoolean getInstanceCompareReference(GenotypeTable align, int index) {
        if ((index < 0) || (index >= align.numberOfTaxa())) {
            throw new IllegalArgumentException("GenotypeTableMaskBoolean: getInstanceCompareReference: unknown index: " + index);
        }
        String name = align.taxaName(index) + " Reference";
        return getInstanceCompareReference(align, align.genotypeRange(index, 0, align.numberOfSites()), name);
    }

    public static GenotypeTableMaskBoolean getInstanceCompareReference(GenotypeTable align, String id) {
        int index = align.taxa().indexOf(id);
        if (index < 0) {
            throw new IllegalArgumentException("GenotypeTableMaskBoolean: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, align.genotypeRange(index, 0, align.numberOfSites()), id + " Reference");
    }

    public static GenotypeTableMaskBoolean getInstanceCompareReference(GenotypeTable align, byte[] ref, String name) {

        if ((align == null) || (ref == null)) {
            throw new IllegalArgumentException("GenotypeTableMaskBoolean: getInstanceCompareReference: alignment or reference can not be null.");
        }

        if (align.numberOfSites() != ref.length) {
            throw new IllegalArgumentException("GenotypeTableMaskBoolean: getInstanceCompareReference: ref length should equal alignment site count.");
        }

        int numMaskColumns = getNumMaskColumns(ref.length);
        byte[][] mask = new byte[align.numberOfTaxa()][numMaskColumns];

        for (int c = 0, n = align.numberOfSites(); c < n; c++) {

            int currentByteCol = c / 8;
            byte currentColMask = (byte) (0x80 >>> (c % 8));

            for (int r = 0, m = align.numberOfTaxa(); r < m; r++) {

                // TERRY - If allele order switch it should still probably match?
                if (align.genotype(r, c) != ref[c]) {
                    mask[r][currentByteCol] = (byte) (mask[r][currentByteCol] | currentColMask);
                }

            }

        }

        return new GenotypeTableMaskBoolean(align, mask, name, MaskType.reference);

    }

    public static GenotypeTableMaskBoolean getInstanceCompareAlignments(GenotypeTable align1, GenotypeTable align2, String name, MaskType type) {

        if ((align1.numberOfTaxa() != align2.numberOfTaxa())
                || (align1.numberOfSites() != align2.numberOfSites())) {
            throw new IllegalArgumentException("GenotypeTableMaskBoolean: getInstanceCompareAlignments: both alignments should have same number of sequences and sites.");
        }

        int numMaskColumns = getNumMaskColumns(align1.numberOfSites());
        byte[][] mask = new byte[align1.numberOfTaxa()][numMaskColumns];

        for (int c = 0, n = align1.numberOfSites(); c < n; c++) {

            int currentByteCol = c / 8;
            byte currentColMask = (byte) (0x80 >>> (c % 8));

            for (int r = 0, m = align1.numberOfTaxa(); r < m; r++) {

                // TERRY - If allele order switch it should still probably match?
                if (align1.genotype(r, c) != align2.genotype(r, c)) {
                    mask[r][currentByteCol] = (byte) (mask[r][currentByteCol] | currentColMask);
                }

            }

        }

        return new GenotypeTableMaskBoolean(align1, mask, name, type);

    }

    private static int getNumMaskColumns(int numSites) {
        int numMaskColumns = numSites / 8;
        if (numSites % 8 > 0) {
            numMaskColumns++;
        }
        return numMaskColumns;
    }

    @Override
    public byte getMask(int taxon, int site) {
        int maskColumn = site / 8;
        int shift = 7 - (site % 8);
        return (byte) ((myMask[taxon][maskColumn] >>> shift) & 0x1);
    }
}
