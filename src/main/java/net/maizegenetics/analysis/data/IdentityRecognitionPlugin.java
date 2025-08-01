/*
 *  IdentityRecognitionPlugin
 * 
 *  Created on Dec 20, 2016
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.analysis.distance.IBSDistanceMatrix3Alleles;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.gui.GenotypeTableMaskGeneticDistance;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.DistanceMatrixUtils;
import net.maizegenetics.util.Tuple;

import javax.swing.*;
import java.awt.*;
import java.util.List;

/**
 * @author Terry Casstevens
 */
public class IdentityRecognitionPlugin extends AbstractPlugin {

    public IdentityRecognitionPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        List<Datum> genotypes = input.getDataOfType(GenotypeTable.class);
        if (genotypes.size() > 1) {
            throw new IllegalArgumentException("IdentityRecognitionPlugin: performFunction: must specify only one genotype table.");
        }

        if (genotypes.size() == 0) {
            List<Datum> distanceMatrix = input.getDataOfType(DistanceMatrix.class);
            if (distanceMatrix.size() != 1) {
                throw new IllegalArgumentException("IdentityRecognitionPlugin: performFunction: must specify only one genotype table or distance matrix.");
            }
            return orderByMostRelated((DistanceMatrix) distanceMatrix.get(0).getData(), distanceMatrix.get(0).getName());
        } else {
            return orderByMostRelated((GenotypeTable) genotypes.get(0).getData(), genotypes.get(0).getName());
        }

    }

    private DataSet orderByMostRelated(DistanceMatrix matrix, String name) {
        DistanceMatrix result = DistanceMatrixUtils.clusterBySmallestDistance(matrix);
        Datum datum = new Datum(name + " Sorted", result, null);
        return new DataSet(datum, this);
    }

    private DataSet orderByMostRelated(GenotypeTable genotype, String name) {

        DistanceMatrix distances = IBSDistanceMatrix3Alleles.getInstance(genotype);

        int numTaxa = genotype.numberOfTaxa();

        double lowestDistance = Double.POSITIVE_INFINITY;
        int taxon1 = -1;
        int taxon2 = -1;
        for (int t1 = 0; t1 < numTaxa; t1++) {
            for (int t2 = t1 + 1; t2 < numTaxa; t2++) {
                if (distances.getDistance(t1, t2) < lowestDistance) {
                    lowestDistance = distances.getDistance(t1, t2);
                    taxon1 = t1;
                    taxon2 = t2;
                }
            }
        }

        Tuple<GenotypeTable, double[]> sortedByGeneticDistance = FilterGenotypeTable.getInstanceTaxaOrderedByGeneticDistance(genotype, taxon1);
        Datum genotypeDatum = new Datum("Genetic Distance with " + genotype.taxaName(taxon1), sortedByGeneticDistance.x, null);

        GenotypeTableMaskGeneticDistance mask = GenotypeTableMaskGeneticDistance.getInstanceCompareReference(sortedByGeneticDistance.x, 0, sortedByGeneticDistance.y);
        Datum maskDatum = new Datum(mask.toString(), mask, null);

        Datum distanceDatum = new Datum("Genetic Distance of " + name, distances, null);

        return new DataSet(new Datum[]{genotypeDatum, maskDatum, distanceDatum}, this);

    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Taxa Identity Recognition";
    }

    @Override
    public String getToolTipText() {
        return "Taxa Identity Recognition";
    }

}
