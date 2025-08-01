package net.maizegenetics.analysis.distance;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.analysis.association.AssociationUtils;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.stats.PCA.ClassicMds;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.TableReportBuilder;

public class MultiDimensionalScalingPlugin extends AbstractPlugin {

    private PluginParameter<Integer> numberOfAxes = new PluginParameter.Builder<>("axes", 5, Integer.class)
            .description("The number of axes or dimensions and associated eigenvalues to be returned by the analysis.")
            .guiName("Number of Axes")
            .build();

    private PluginParameter<Boolean> removeNaN = new PluginParameter.Builder<>("removeNaN", true, Boolean.class)
            .description("Remove NaNs from matrix before performing MDS")
            .guiName("Remove NaNs")
            .build();

    public MultiDimensionalScalingPlugin(Frame parent, boolean interactive) {
        super(parent, interactive);
    }

    @Override
    public DataSet processData(DataSet input) {
        List<Datum> resultList = new ArrayList<>();
        List<Datum> myDatumList = input.getDataOfType(DistanceMatrix.class);
        if (myDatumList.size() < 1) {
            throw new RuntimeException("MDS requires a Distance Matrix as input.");
        }
        int ndata = myDatumList.size();

        fireProgress(10);
        int counter = 0;
        for (Datum myDatum : myDatumList) {
            DistanceMatrix myDistanceMatrix = (DistanceMatrix) myDatum.getData();
            if (removeNaN()) {
                myDistanceMatrix = RemoveNaNFromDistanceMatrixPlugin.runPlugin(myDistanceMatrix);
            }
            ClassicMds myMDS = new ClassicMds(myDistanceMatrix);
            int numberOfAxesToReport = numberOfAxes.value();

            //get requested number of axes and package as a Phenotype (covariates)
            List<PhenotypeAttribute> attrList = new ArrayList<>();
            List<ATTRIBUTE_TYPE> typeList = new ArrayList<>();
            attrList.add(new TaxaAttribute(myDistanceMatrix.getTaxaList()));
            typeList.add(ATTRIBUTE_TYPE.taxa);
            int ntaxa = myDistanceMatrix.getSize();
            for (int i = 0; i < numberOfAxesToReport; i++) {
                typeList.add(ATTRIBUTE_TYPE.covariate);
                float[] floatValues = AssociationUtils.convertDoubleArrayToFloat(myMDS.getPrincipalCoordinate(i));
                NumericAttribute pcAttr = new NumericAttribute("PC" + (i + 1), floatValues, new OpenBitSet(ntaxa));
                attrList.add(pcAttr);
            }
            Phenotype myPhenotype = new PhenotypeBuilder().fromAttributeList(attrList, typeList).build().get(0);
            String name = "MDS_PCs_" + myDatum.getName();
            String comment = "Principal Coordinates from MDS analysis of " + myDatum.getName();
            resultList.add(new Datum(name, myPhenotype, comment));

            //get requested number of eigenvalues and package as a TableReport
            String[] colnames = new String[]{"PC", "eigenvalue"};
            TableReportBuilder reportBuilder = TableReportBuilder.getInstance("", colnames);
            for (int i = 0; i < numberOfAxesToReport; i++) {
                reportBuilder.add(new Object[]{new Integer(i + 1), new Double(myMDS.getEigenvalue(i))});
            }
            name = "MDS_Eigenvalues_" + myDatum.getName();
            comment = "Eigenvalues from MDS analysis of " + myDatum.getName();
            resultList.add(new Datum(name, reportBuilder.build(), comment));

            counter++;
            fireProgress(Math.min(99, counter * 100 / ndata));
        }

        return new DataSet(resultList, this);
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = MultiDimensionalScalingPlugin.class.getResource("/net/maizegenetics/analysis/images/pca.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "MDS";
    }

    @Override
    public String getToolTipText() {
        return "Perform classic multidimensional scaling (principal coordinate analysis)";
    }

    @Override
    public String pluginDescription() {
        return "MultiDimensionalScalingPlugin takes a DistanceMatrix as input and performs classic multidimensional scaling, which is also know as principal coordinate analysis (PCO).";
    }
    
    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/MDS/MDS";
    }

    /**
     * The number of axes or dimensions and associated eigenvalues to be
     * returned by the analysis.
     *
     * @return Number of Axes
     */
    public Integer numberOfAxes() {
        return numberOfAxes.value();
    }

    /**
     * Set Number of Axes. The number of axes or dimensions and associated
     * eigenvalues to be returned by the analysis.
     *
     * @param value Number of Axes
     *
     * @return this plugin
     */
    public MultiDimensionalScalingPlugin numberOfAxes(Integer value) {
        numberOfAxes = new PluginParameter<>(numberOfAxes, value);
        return this;
    }

    /**
     * Remove NaNs from matrix before performing MDS
     *
     * @return Remove NaNs
     */
    public Boolean removeNaN() {
        return removeNaN.value();
    }

    /**
     * Set Remove NaNs. Remove NaNs from matrix before performing MDS
     *
     * @param value Remove NaNs
     *
     * @return this plugin
     */
    public MultiDimensionalScalingPlugin removeNaN(Boolean value) {
        removeNaN = new PluginParameter<>(removeNaN, value);
        return this;
    }

}
