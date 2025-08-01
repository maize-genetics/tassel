/*
 * NumericalGenotypePlugin.java
 *
 * Created on May 8, 2008
 *
 */
package net.maizegenetics.analysis.numericaltransform;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.score.ReferenceProbabilityBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.List;

/**
 * Numerical Genotype Plugin.
 *
 * @author Terry Casstevens
 * @author Peter Bradbury
 */
public class NumericalGenotypePlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(NumericalGenotypePlugin.class);

    public static enum TRANSFORM_TYPE {
        as_minor, as_major, as_missing, use_all_alleles, collapse, separated
    };
    
    private TRANSFORM_TYPE myTransformType = TRANSFORM_TYPE.as_minor;

    /** Creates a new instance of the NumericalGenotypePlugin */
    public NumericalGenotypePlugin(Frame parentFrame, boolean isInteractive){
        super(parentFrame, isInteractive);
    }

    public NumericalGenotypePlugin(){
        super(null, false);
    }
    
    @Override
    public DataSet processData(DataSet input) {
 
        List<Datum> datumList = input.getDataOfType(GenotypeTable.class);

        //check size of datumList, throw error if not equal to one
        if (datumList.size() != 1){
            throw new IllegalArgumentException("NumericalGenotypePlugin: select exactly one genotype dataset to transform.");
        }

    	//only the as_minor method is implemented so just display an information Dialog
        String msg = "Genotypes will be converted to the probability that an allele selected at random at a site is "
                + "\nthe major allele. In other words, homozygous major is 1.0, homozygous minor is 0.0, "
                + "\nand heterozygous is 0.5. The results are equivalent to the TASSEL version 4 \"collapse\" method.";
        
        if (isInteractive()) {
            int convertGenotype = JOptionPane.showConfirmDialog(getParentFrame(), msg, "ReferenceProbability from Genotype", JOptionPane.OK_CANCEL_OPTION, JOptionPane.INFORMATION_MESSAGE);
            if (convertGenotype == JOptionPane.CANCEL_OPTION || convertGenotype == JOptionPane.CLOSED_OPTION) return null;
        } else {
        	myLogger.info(msg);
        }
        
        
        //load the GenotypeTable.
        GenotypeTable myGenotype = (GenotypeTable) datumList.get(0).getData();
        GenotypeTable myNewGenotype;
        StringBuilder nameBuilder = new StringBuilder(datumList.get(0).getName());
        nameBuilder.append("_with_Probability");
        StringBuilder commentBuilder = new StringBuilder();
        if (datumList.get(0).getComment() != null) commentBuilder.append(datumList.get(0).getComment());
        
        switch(myTransformType) {
        case as_minor:
        	myNewGenotype = setAlternateMinorAllelesToMinor(myGenotype);
        	commentBuilder.append("\nReference Probability computed by setting alternate minor alleles");
        	commentBuilder.append("\nto the most common minor allele.");
        	break;
        case as_major:
        	myNewGenotype = setAlternateMinorAllelesToMajor(myGenotype);
        	commentBuilder.append("\nReference Probability computed by setting alternate minor alleles");
        	commentBuilder.append("\nto the major allele.");
        	break;
        case as_missing:
        	myNewGenotype = setAlternateMinorAllelesToMissing(myGenotype);
        	commentBuilder.append("\nReference Probability computed by setting alternate minor alleles");
        	commentBuilder.append("\nto missing.");
        	break;
        case use_all_alleles:
        	myNewGenotype = transformToAlleleProbability(myGenotype);
        	commentBuilder.append("\nWith allele probabilities computed.");
        	break;
        default:
        	throw new IllegalStateException("NumericalGenotypePlugin: processData: unknown transform type: " + myTransformType.toString());
        }
        
        Datum newDatum = new Datum(nameBuilder.toString(), myNewGenotype, commentBuilder.toString());
        return new DataSet(newDatum, this);

    }

    public static GenotypeTable setAlternateMinorAllelesToMinor(GenotypeTable myGenotype) {
        int nsites = myGenotype.numberOfSites();
        int ntaxa = myGenotype.numberOfTaxa();

        float[][] data = GenotypeTableUtils.convertGenotypeToFloatProbability(myGenotype, false);

        //build new ReferenceProbability
        ReferenceProbabilityBuilder refBuilder = ReferenceProbabilityBuilder.getInstance(ntaxa, nsites, myGenotype.taxa());
        for (int t = 0; t < ntaxa; t++) {
            refBuilder.addTaxon(t, data[t]);
        }

        //build new GenotypeTable
        return GenotypeTableBuilder.getInstance(myGenotype.genotypeMatrix(), myGenotype.positions(),
                myGenotype.taxa(), myGenotype.depth(), myGenotype.alleleProbability(), refBuilder.build(), myGenotype.dosage(),
                myGenotype.annotations());
    }
    
    public static GenotypeTable setAlternateMinorAllelesToMajor(GenotypeTable myGenotype) {
    	//TODO implement
    	throw new UnsupportedOperationException();
    }

    public static GenotypeTable setAlternateMinorAllelesToMissing(GenotypeTable myGenotype) {
    	//TODO implement
        throw new UnsupportedOperationException();
    }

    public static GenotypeTable transformToAlleleProbability(GenotypeTable myGenotype) {
    	//TODO implement
        throw new UnsupportedOperationException();
    }

    @Override
    public String pluginDescription(){
        return "This plugin creates a numerical genotype table for input genotype data";
    }

    @Override
    public  ImageIcon getIcon(){
        URL imageURL = NumericalGenotypePlugin.class.getResource("/net/maizegenetics/analysis/images/NumericalGenotype.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName(){
        return "Numerical Genotype";
    }


    @Override
    public String getToolTipText(){
        return "Numerical Genotype";
    }

    // Only one type currently supported
    private void setTransformType(TRANSFORM_TYPE type) {
    	myTransformType = type;
    }
}
