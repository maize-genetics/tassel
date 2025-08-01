/*
 * GenotypeTableMaskGeneticDistance
 */
package net.maizegenetics.gui;

import java.awt.*;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.taxa.Taxon;

/**
 *
 * @author Terry Casstevens
 */
public class GenotypeTableMaskGeneticDistance extends AbstractGenotypeTableMask {

    private final int myTaxonReference;
    private final GenotypeTable myAlignment;
    private final double[] myDistances;

    private GenotypeTableMaskGeneticDistance(GenotypeTable align, int taxonReference, double[] distances, String name, Color color) {
        super(align, name, color, GenotypeTableMask.MaskType.reference);
        myTaxonReference = taxonReference;
        myAlignment = align;
        myDistances = distances;
    }

    public static GenotypeTableMaskGeneticDistance getInstanceCompareReference(GenotypeTable align, Taxon id, double[] distances) {
        int index = align.taxa().indexOf(id);
        if (index < 0) {
            throw new IllegalArgumentException("GenotypeTableMaskGeneticDistance: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index, distances);
    }

    public static GenotypeTableMaskGeneticDistance getInstanceCompareReference(GenotypeTable align, String id, double[] distances) {
        int index = align.taxa().indexOf(id);
        if (index < 0) {
            throw new IllegalArgumentException("GenotypeTableMaskGeneticDistance: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index, distances);
    }

    public static GenotypeTableMaskGeneticDistance getInstanceCompareReference(GenotypeTable align, int index, double[] distances) {
        if ((index < 0) || (index >= align.numberOfTaxa())) {
            throw new IllegalArgumentException("GenotypeTableMaskGeneticDistance: getInstanceCompareReference: unknown index: " + index);
        }
        String name = align.taxaName(index) + " Genetic Distance";
        return new GenotypeTableMaskGeneticDistance(align, index, distances, name, null);
    }

    @Override
    public byte getMask(int taxon, int site) {
        return (byte) (myDistances[taxon] * 255.0);
    }
}
