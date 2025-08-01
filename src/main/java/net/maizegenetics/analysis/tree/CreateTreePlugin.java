/*
 * CreateTreePlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.analysis.tree;

import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.tree.NeighborJoiningTree;
import net.maizegenetics.taxa.tree.Tree;
import net.maizegenetics.taxa.tree.UPGMATree;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public class CreateTreePlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(CreateTreePlugin.class);

    public static enum CLUSTERING_METHOD {

        Neighbor_Joining, UPGMA
    }

    private PluginParameter<CLUSTERING_METHOD> myClusteringMethod = new PluginParameter.Builder<>("clusteringMethod", CLUSTERING_METHOD.Neighbor_Joining, CLUSTERING_METHOD.class)
            .description("")
            .range(CLUSTERING_METHOD.values())
            .build();
    private PluginParameter<Boolean> mySaveDistanceMatrix = new PluginParameter.Builder<>("saveDistanceMatrix", true, Boolean.class)
            .description("")
            .build();

    private Datum myGenotypeTable = null;
    private Datum myDistanceMatrix = null;

    /**
     * Creates a new instance of CreateTreePlugin
     */
    public CreateTreePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {

        myGenotypeTable = null;

        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
        if (alignInList.size() >= 1) {
            myGenotypeTable = alignInList.get(0);
            return;
        }

        List<Datum> matrixList = input.getDataOfType(DistanceMatrix.class);
        if (matrixList.size() >= 1) {
            myDistanceMatrix = matrixList.get(0);
            return;
        }

        throw new IllegalArgumentException("CreateTreePlugin: Invalid selection.  Please select a Genotype Table.");

    }

    @Override
    public DataSet processData(DataSet input) {

        Datum datum = null;
        DistanceMatrix distanceMatrix = null;

        if (myGenotypeTable != null) {
            datum = myGenotypeTable;
            distanceMatrix = IBSDistanceMatrix.getInstance((GenotypeTable) myGenotypeTable.getData(), this);
        } else if (myDistanceMatrix != null) {
            datum = myDistanceMatrix;
            distanceMatrix = (DistanceMatrix) myDistanceMatrix.getData();
        } else {
            throw new IllegalArgumentException("CreateTreePlugin: Invalid selection.  Please select a Genotype Table or Distance Matrix.");
        }

        List<Datum> results = new ArrayList<>();

        if (clusteringMethod() == CLUSTERING_METHOD.Neighbor_Joining) {
            Tree theTree = new NeighborJoiningTree(distanceMatrix);
            results.add(new Datum("Tree:" + datum.getName(), theTree, "NJ Tree"));
        } else if (clusteringMethod() == CLUSTERING_METHOD.UPGMA) {
            Tree theTree = new UPGMATree(distanceMatrix);
            results.add(new Datum("Tree:" + datum.getName(), theTree, "UPGMA Tree"));
        } else {
            throw new IllegalArgumentException("CreateTreePlugin: processData: Unknown clustering method: " + clusteringMethod());
        }

        if (saveDistanceMatrix()) {
            results.add(new Datum("Matrix:" + datum.getName(), distanceMatrix, "Distance Matrix"));
        }

        return new DataSet(results, this);

    }

    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Cladogram/Cladogram";
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    @Override
    public ImageIcon getIcon() {
        URL imageURL = CreateTreePlugin.class.getResource("/net/maizegenetics/analysis/images/Tree.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    @Override
    public String getButtonName() {
        return "Create Tree";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Create a tree";
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    public Tree runPlugin(DataSet input) {
        return (Tree) performFunction(input).getData(0).getData();
    }

    /**
     * Clustering Method
     *
     * @return Clustering Method
     */
    public CLUSTERING_METHOD clusteringMethod() {
        return myClusteringMethod.value();
    }

    /**
     * Set Clustering Method. Clustering Method
     *
     * @param value Clustering Method
     *
     * @return this plugin
     */
    public CreateTreePlugin clusteringMethod(CLUSTERING_METHOD value) {
        myClusteringMethod = new PluginParameter<>(myClusteringMethod, value);
        return this;
    }

    /**
     * Save Distance Matrix
     *
     * @return Save Distance Matrix
     */
    public Boolean saveDistanceMatrix() {
        return mySaveDistanceMatrix.value();
    }

    /**
     * Set Save Distance Matrix. Save Distance Matrix
     *
     * @param value Save Distance Matrix
     *
     * @return this plugin
     */
    public CreateTreePlugin saveDistanceMatrix(Boolean value) {
        mySaveDistanceMatrix = new PluginParameter<>(mySaveDistanceMatrix, value);
        return this;
    }

}
