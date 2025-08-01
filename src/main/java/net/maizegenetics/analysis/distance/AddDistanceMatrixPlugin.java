/*
 *  AddDistanceMatrixPlugin
 * 
 *  Created on Feb 22, 2016
 */
package net.maizegenetics.analysis.distance;

import java.awt.Frame;
import java.util.List;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.DistanceMatrixBuilder;
import net.maizegenetics.taxa.distance.DistanceMatrixWithCounts;
import net.maizegenetics.util.GeneralAnnotation;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class AddDistanceMatrixPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(AddDistanceMatrixPlugin.class);

    private KinshipPlugin.KINSHIP_METHOD myMethod = null;

    public AddDistanceMatrixPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {

        List<Datum> matrices = input.getDataOfType(DistanceMatrix.class);
        if (matrices.size() < 2) {
            throw new IllegalArgumentException("AddDistanceMatrixPlugin: preProcessParameters: must input at least two Distance Matrices to add.");
        }

        String matrixType = null;
        GeneralAnnotation annotations = ((DistanceMatrix) matrices.get(0).getData()).annotations();
        try {
            matrixType = annotations.getTextAnnotation(DistanceMatrixBuilder.MATRIX_TYPE)[0];
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("AddDistanceMatrixPlugin: preProcessParameters: all matrices "
                    + " must have annotation: " + DistanceMatrixBuilder.MATRIX_TYPE
                    + ". Matrices must be exported with a more recent build of Tassel.");
        }

        for (int i = 1; i < matrices.size(); i++) {
            String currentType = null;
            GeneralAnnotation currentAnno = ((DistanceMatrix) matrices.get(i).getData()).annotations();
            try {
                currentType = currentAnno.getTextAnnotation(DistanceMatrixBuilder.MATRIX_TYPE)[0];
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("AddDistanceMatrixPlugin: preProcessParameters: all matrices "
                        + " must have annotation: " + DistanceMatrixBuilder.MATRIX_TYPE
                        + ". Matrices must be exported with a more recent build of Tassel.");
            }
            if (!matrixType.equals(currentType)) {
                throw new IllegalStateException("AddDistanceMatrixPlugin: preProcessParameters: all matrices must have same matrix type.");
            }
        }

        myMethod = KinshipPlugin.KINSHIP_METHOD.valueOf(matrixType);

    }

    @Override
    public DataSet processData(DataSet input) {

        if (myMethod == KinshipPlugin.KINSHIP_METHOD.Centered_IBS) {
            List<Datum> matricesList = input.getDataOfType(DistanceMatrix.class);
            DistanceMatrix[] matrices = new DistanceMatrix[matricesList.size()];
            for (int i = 0; i < matricesList.size(); i++) {
                matrices[i] = (DistanceMatrix) matricesList.get(i).getData();
            }
            DistanceMatrix result = EndelmanDistanceMatrix.addEndelmanDistance(matrices, this);
            return new DataSet(new Datum("SumDistanceMatrix", result, null), this);
        } else if (myMethod == KinshipPlugin.KINSHIP_METHOD.Normalized_IBS) {
            List<Datum> matricesList = input.getDataOfType(DistanceMatrixWithCounts.class);
            DistanceMatrixWithCounts[] matrices = new DistanceMatrixWithCounts[matricesList.size()];
            for (int i = 0; i < matricesList.size(); i++) {
                matrices[i] = (DistanceMatrixWithCounts) matricesList.get(i).getData();
            }
            DistanceMatrix result = GCTADistanceMatrix.addGCTADistance(matrices, this);
            return new DataSet(new Datum("SumDistanceMatrix", result, null), this);
        } else {
            throw new UnsupportedOperationException("AddDistanceMatrixPlugin: processData: unsupported matrix type: " + myMethod);
        }

    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Add Distance Matrix";
    }

    @Override
    public String getToolTipText() {
        return "Add Distance Matrix";
    }

}
