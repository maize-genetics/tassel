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

public class BoxCoxTransformationPlugin extends AbstractPlugin{
    
    private PluginParameter<Double> startLambda = 
            new PluginParameter.Builder<>("startLambda",-5.0,Double.class)
            .description("Parameter to set the starting point for the Lambda search in Box Cox.  "
                    + "The algorithm will start at this value and iterate by Lambda Step until it reaches End Lambda.")
            .build();
    private PluginParameter<Double> endLambda = 
            new PluginParameter.Builder<>("endLambda",5.0,Double.class)
            .description("Parameter to set the ending point for the Lambda search in Box Cox. "
                    + "The algorithm will start at Start Lambda and iterate by Lambda Step until it reaches this value.")
            .build();
    private PluginParameter<Double> stepLambda = 
            new PluginParameter.Builder<>("LambdaStep",0.2,Double.class)
            .description("Parameter to set the iteration step for the Lambda search in Box Cox. "
                    + "The algorithm will start at Start Lambda and iterate by this value until it reaches the End Lambda.")
            .build();
    private PluginParameter<Boolean> addSmallValue =
            new PluginParameter.Builder<>("addSmallValue", true, Boolean.class)
            .description("Boolean to allow for a small random value to be added to each observed value. The value is calculated by the following:\n"
                    + "rand(0:1) * .5 * minObservedValue")
            .build();
    private PluginParameter<Long> randomSeed = 
            new PluginParameter.Builder<>("randomSeed",12345l,Long.class)
            .description("Random Seed to be used in the Random Small Number Generation")
            .dependentOnParameter(addSmallValue)
            .build();
    
    
    public BoxCoxTransformationPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    public DataSet processData(DataSet input) {
        //do whatever your plugin does
        List<Datum> datumList = input.getDataOfType(Phenotype.class);

        //check size of datumList, throw error if not equal to one
        if (datumList.size() != 1){
            throw new IllegalArgumentException("BoxCoxTransformationPluging: select exactly one phenotype dataset to average.");
        }
        if(startLambda.value()>endLambda.value()) {
            throw new IllegalArgumentException("Start Lambda must be smaller than End Lambda");
        }
        
        BoxCoxTransformation boxCoxPheno = new BoxCoxTransformation();
        try {
        Phenotype myPhenotype = boxCoxPheno.applyBoxCox((Phenotype) datumList.get(0).getData(),addSmallValue.value(),randomSeed.value(),
                startLambda.value(),endLambda.value(),stepLambda.value());
        
        if (myPhenotype != null) {
            String name = myPhenotype.name();

            Datum td = new Datum(name, myPhenotype, null);
            
            //todo need to add logic of directories.
            DataSet tds = new DataSet(td, this);
            return tds;
        }
        }
        catch(Exception e) {
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
        return "Box Cox Transformation";
    }
    @Override
    public String getToolTipText() {
        return "Box Cox Transformation.";
    }
    
    @Override
    public String getCitation() {
        return "Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES. "
                + "(2007) TASSEL: Software for association mapping of complex traits in diverse "
                + "samples. Bioinformatics 23:2633­2635.";
    }
    
    @Override
    public String pluginUserManualURL() {
        //TODO add in the documentation link once it is written up
        return "https://bitbucket.org/tasseladmin/tassel­5­source/wiki/UserManual/Kinship/Kinship"; 
    } 
    
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(BoxCoxTransformationPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public Phenotype runPlugin(DataSet input) {
        return (Phenotype) performFunction(input).getData(0).getData();
    }

    /**
     * Add Small Value
     *
     * @return Add Small Value
     */
    public Boolean addSmallValue() {
        return addSmallValue.value();
    }

    /**
     * Set Add Small Value. Add Small Value
     *
     * @param value Add Small Value
     *
     * @return this plugin
     */
    public BoxCoxTransformationPlugin addSmallValue(Boolean value) {
        addSmallValue = new PluginParameter<>(addSmallValue, value);
        return this;
    }

    /**
     * Random Seed
     *
     * @return Random Seed
     */
    public Long randomSeed() {
        return randomSeed.value();
    }

    /**
     * Set Random Seed. Random Seed
     *
     * @param value Random Seed
     *
     * @return this plugin
     */
    public BoxCoxTransformationPlugin randomSeed(Long value) {
        randomSeed = new PluginParameter<>(randomSeed, value);
        return this;
    }

    /**
     * Start Lambda
     *
     * @return Start Lambda
     */
    public Double startLambda() {
        return startLambda.value();
    }

    /**
     * Set Start Lambda. Start Lambda
     *
     * @param value Start Lambda
     *
     * @return this plugin
     */
    public BoxCoxTransformationPlugin startLambda(Double value) {
        startLambda = new PluginParameter<>(startLambda, value);
        return this;
    }

    /**
     * End Lambda
     *
     * @return End Lambda
     */
    public Double endLambda() {
        return endLambda.value();
    }

    /**
     * Set End Lambda. End Lambda
     *
     * @param value End Lambda
     *
     * @return this plugin
     */
    public BoxCoxTransformationPlugin endLambda(Double value) {
        endLambda = new PluginParameter<>(endLambda, value);
        return this;
    }

    /**
     * Step Lambda
     *
     * @return Step Lambda
     */
    public Double stepLambda() {
        return stepLambda.value();
    }

    /**
     * Set Step Lambda. Step Lambda
     *
     * @param value Step Lambda
     *
     * @return this plugin
     */
    public BoxCoxTransformationPlugin stepLambda(Double value) {
        stepLambda = new PluginParameter<>(stepLambda, value);
        return this;
    }

}
