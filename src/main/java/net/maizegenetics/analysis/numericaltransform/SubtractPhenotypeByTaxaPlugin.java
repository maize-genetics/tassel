package net.maizegenetics.analysis.numericaltransform;

import java.awt.Frame;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;

public class SubtractPhenotypeByTaxaPlugin extends AbstractPlugin{
    private PluginParameter<Boolean> useAbsoluteDifference =
            new PluginParameter.Builder<>("useAbsDiff", false, Boolean.class)
            .description("Boolean to allow the user to use Absolute Difference when doing the subtraction.")
            .build();
    
    public SubtractPhenotypeByTaxaPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    public DataSet processData(DataSet input) {
        //do whatever your plugin does
        List<Datum> datumList = input.getDataOfType(Phenotype.class);

        //check size of datumList, throw error if not equal to one
        if (datumList.size() != 2){
            throw new IllegalArgumentException("SubtractPhenotypeByTaxaPlugin: select exactly two phenotype dataset to subtract.");
        }
        
        SubtractPhenotype subPheno = new SubtractPhenotype();
        Phenotype myPhenotype1 = (Phenotype) datumList.get(0).getData();
        Phenotype myPhenotype2 = (Phenotype) datumList.get(1).getData();
        
        try {
        Phenotype subtractedPhenotype = subPheno.subtractPhenotype(myPhenotype1, myPhenotype2,useAbsoluteDifference.value());
        if (myPhenotype1 != null) {
            String name = subtractedPhenotype.name();

            Datum td = new Datum(name, subtractedPhenotype, null);
            
            //todo need to add logic of directories.
            DataSet tds = new DataSet(td, this);
            return tds;
        }
        }
        catch(Exception e) {
//            throw new IllegalArgumentException("Error building Phenotype Object.  "
//                    + "Check to make sure the dimensions of your data match the number of taxa and the number of attributes.");
            throw new IllegalArgumentException(e.getMessage());
        }
        return null;  // Note: this can return null
    }
    
    @Override
    public ImageIcon getIcon() {
        return null;
    }
    @Override
    public String getButtonName() {
        return "Subtract Phenotypes By Taxa";
    }
    @Override
    public String getToolTipText() {
        return "Subtract Phenotypes for Matching Taxa.";
    }
    
    @Override
    public String getCitation() {
        return "Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES. "
                + "(2007) TASSEL: Software for association mapping of complex traits in diverse "
                + "samples. Bioinformatics 23:2633­2635.";
    }
    
    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel­5­source/wiki/UserManual/Kinship/Kinship"; 
    }   
    
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(SubtractPhenotypeByTaxaPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public Phenotype runPlugin(DataSet input) {
        return (Phenotype) performFunction(input).getData(0).getData();
    }

    /**
     * Use Abs Diff
     *
     * @return Use Abs Diff
     */
    public Boolean useAbsoluteDifference() {
        return useAbsoluteDifference.value();
    }

    /**
     * Set Use Abs Diff. Use Abs Diff
     *
     * @param value Use Abs Diff
     *
     * @return this plugin
     */
    public SubtractPhenotypeByTaxaPlugin useAbsoluteDifference(Boolean value) {
        useAbsoluteDifference = new PluginParameter<>(useAbsoluteDifference, value);
        return this;
    }


}
