/*
 * UnionAlignmentPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.CombineGenotypeTable;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.GenotypePhenotypeBuilder;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeBuilder;

import java.awt.Frame;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import javax.swing.*;
import java.io.StringWriter;
import java.net.URL;
import java.util.List;

/**
 *
 * @author Ed Buckler
 */
public class UnionAlignmentPlugin extends AbstractPlugin {

    /**
     * Creates a new instance of UnionAlignmentPlugin
     */
    public UnionAlignmentPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {
        Datum joinedDatum = processData(input, true);
        if (joinedDatum == null) {
            return null;
        }
        DataSet output = new DataSet(joinedDatum, this);
        //I am setting the firing class as the metadata - so that the control panel know where the event is coming from
        fireDataSetReturned(new PluginEvent(output, UnionAlignmentPlugin.class));
        return output;
    }

    protected Datum processData(DataSet input, boolean isUnion) {

        try {
            Datum outDatum = null;
            String userMessage = "This action requires multiple items be simultaneously selected from the data tree "
                    + "(Ctrl + mouse click).  Please select genotype, trait, and population structure data to isUnion "
                    + "from the data tree.";

            GenotypeTable aa = null;
            Phenotype ca = null;
            List<Datum> aaVector = input.getDataOfType(GenotypeTable.class);
            List<Datum> caVector = input.getDataOfType(Phenotype.class);
            if ((aaVector.size() + caVector.size()) < 2) {
                JOptionPane.showMessageDialog(getParentFrame(), userMessage);
                return null;
            }

            StringWriter sw = new StringWriter();
            Object result = null;
            if (aaVector.size() == 1) {
                aa = (GenotypeTable) aaVector.get(0).getData();
            } else if (aaVector.size() > 1) {
                GenotypeTable[] temp = new GenotypeTable[aaVector.size()];
                for (int i = 0; i < aaVector.size(); i++) {
                    temp[i] = (GenotypeTable) aaVector.get(i).getData();
                }
                aa = CombineGenotypeTable.getInstance(temp, isUnion);
                result = aa;
            }
            if (caVector.size() > 0) {
                ca = (Phenotype) caVector.get(0).getData();
            }
            for (int i = 1; i < caVector.size(); i++) {
                Phenotype ta = (Phenotype) caVector.get(i).getData();
                if (isUnion) {
                    ca = new PhenotypeBuilder().fromPhenotype(ca).fromPhenotype(ta).unionJoin().build().get(0);
                } else {
                    ca = new PhenotypeBuilder().fromPhenotype(ca).fromPhenotype(ta).intersectJoin().build().get(0);
                }
                result = ca;
            }
            if ((ca != null) && (aa != null)) {
                //then make a concatenated alignment
                GenotypePhenotype aac = null;
                if (isUnion) {
                    aac = new GenotypePhenotypeBuilder().genotype(aa).phenotype(ca).union().build();
                } else {
                    aac = new GenotypePhenotypeBuilder().genotype(aa).phenotype(ca).intersect().build();
                }
                result = aac;
            }
            String theName = this.getConcatenatedName(input);
            if (isUnion) {
                sw.append("Union Join\n");
            } else {
                sw.append("Intersect Join\n");
            }
            String theComment = sw.toString();
            outDatum = new Datum(theName, result, theComment);
            return outDatum;
        } finally {
            fireProgress(100);
        }
    }

    protected String getConcatenatedName(DataSet theTDS) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < theTDS.getSize(); i++) {
            sb.append(theTDS.getData(i).getName());
            if (i + 1 != theTDS.getSize()) {
                sb.append(" + ");
            }
        }
        return sb.toString();
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    @Override
    public ImageIcon getIcon() {
        URL imageURL = UnionAlignmentPlugin.class.getResource("/net/maizegenetics/analysis/images/UnionJoin.gif");
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
        return "Union Join";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Join Datasets by Union of Taxa";
    }
}
