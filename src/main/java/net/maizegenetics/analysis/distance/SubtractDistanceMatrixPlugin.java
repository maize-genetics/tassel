/*
 *  SubtractDistanceMatrixPlugin
 * 
 *  Created on Nov 23, 2015
 */
package net.maizegenetics.analysis.distance;

import java.awt.Frame;
import java.util.List;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.DistanceMatrixBuilder;
import net.maizegenetics.taxa.distance.DistanceMatrixWithCounts;
import net.maizegenetics.taxa.distance.ReadDistanceMatrix;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class SubtractDistanceMatrixPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(SubtractDistanceMatrixPlugin.class);

    private PluginParameter<String> myWholeMatrix = new PluginParameter.Builder<>("wholeMatrix", null, String.class)
            .description("The filename of the whole matrix which will be used to subtract input sub-matrices.")
            .inFile()
            .required(true)
            .build();

    private DistanceMatrix myCurrentlyLoadedMatrix = null;

    public SubtractDistanceMatrixPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> matrices = input.getDataOfType(DistanceMatrix.class);
        if (matrices.isEmpty()) {
            throw new IllegalArgumentException("SubtractDistanceMatrixPlugin: preProcessParameters: must input at least one Distance Matrix to subtract.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        if (myCurrentlyLoadedMatrix == null) {
            myCurrentlyLoadedMatrix = ReadDistanceMatrix.readDistanceMatrix(wholeMatrix());
        }

        GeneralAnnotation annotations = myCurrentlyLoadedMatrix.annotations();

        String matrixType = null;
        try {
            matrixType = annotations.getTextAnnotation(DistanceMatrixBuilder.MATRIX_TYPE)[0];
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("SubtractDistanceMatrixPlugin: processData: the whole matrix: "
                    + wholeMatrix() + " doesn't have annotation: " + DistanceMatrixBuilder.MATRIX_TYPE
                    + ". The matrix must be exported with a more recent build of Tassel.");
        }
        if (matrixType.equals(KinshipPlugin.KINSHIP_METHOD.Centered_IBS.toString())) {
            List<Datum> matricesList = input.getDataOfType(DistanceMatrix.class);
            DistanceMatrix[] matrices = new DistanceMatrix[matricesList.size()];
            for (int i = 0; i < matricesList.size(); i++) {
                matrices[i] = (DistanceMatrix) matricesList.get(i).getData();
            }
            DistanceMatrix result = EndelmanDistanceMatrix.subtractEndelmanDistance(matrices, myCurrentlyLoadedMatrix, this);
            return new DataSet(new Datum(Utils.getFilename(wholeMatrix()) + "_Rest", result, null), this);
        } else if (matrixType.equals(KinshipPlugin.KINSHIP_METHOD.Normalized_IBS.toString())) {
            List<Datum> matricesList = input.getDataOfType(DistanceMatrixWithCounts.class);
            if (matricesList.isEmpty()) {
                throw new IllegalArgumentException("SubtractDistanceMatrixPlugin: processData: must input at least one Distance Matrix with counts to subtract.");
            }
            DistanceMatrixWithCounts[] matrices = new DistanceMatrixWithCounts[matricesList.size()];
            for (int i = 0; i < matricesList.size(); i++) {
                matrices[i] = (DistanceMatrixWithCounts) matricesList.get(i).getData();
            }
            DistanceMatrixWithCounts tempMatrix = null;
            try {
                tempMatrix = (DistanceMatrixWithCounts) myCurrentlyLoadedMatrix;
            } catch (ClassCastException ex) {
                myLogger.debug(ex.getMessage(), ex);
                throw new IllegalArgumentException("SubtractDistanceMatrixPlugin: processData: whole matrix must have counts.");
            }
            DistanceMatrix result = GCTADistanceMatrix.subtractGCTADistance(matrices, tempMatrix, this);
            return new DataSet(new Datum(Utils.getFilename(wholeMatrix()) + "_Rest", result, null), this);
        } else {
            throw new UnsupportedOperationException("SubstractDistanceMatrixPlugin: processData: unsupported matrix type: " + matrixType);
        }

    }

    /**
     * The filename of the whole matrix which will be used to subtract input
     * sub-matrices.
     *
     * @return Whole Matrix
     */
    public String wholeMatrix() {
        return myWholeMatrix.value();
    }

    /**
     * Set Whole Matrix. The filename of the whole matrix which will be used to
     * subtract input sub-matrices.
     *
     * @param value Whole Matrix
     *
     * @return this plugin
     */
    public SubtractDistanceMatrixPlugin wholeMatrix(String value) {
        myWholeMatrix = new PluginParameter<>(myWholeMatrix, value);
        myCurrentlyLoadedMatrix = null;
        return this;
    }

    public SubtractDistanceMatrixPlugin wholeMatrix(DistanceMatrix matrix) {
        myWholeMatrix = new PluginParameter<>(myWholeMatrix, matrix.getTableTitle());
        myCurrentlyLoadedMatrix = matrix;
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Subtract Distance Matrix";
    }

    @Override
    public String getToolTipText() {
        return "Subtract Distance Matrix";
    }

}
