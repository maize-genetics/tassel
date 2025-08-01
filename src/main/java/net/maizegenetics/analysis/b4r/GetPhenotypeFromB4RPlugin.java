package net.maizegenetics.analysis.b4r;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.stream.Collectors;

import javax.swing.ImageIcon;

import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

public class GetPhenotypeFromB4RPlugin extends net.maizegenetics.plugindef.AbstractPlugin {
    private PluginParameter<String> myB4RDBName = new PluginParameter.Builder<>("b4rdb", null, String.class)
            .required(true)
            .guiName("B4R Database Name")
            .description("B4R Database Name")
            .build();

    private PluginParameter<String> myB4RHost = new PluginParameter.Builder<>("b4rhost", "localhost", String.class)
            .guiName("B4R Host Name")
            .description("B4R Host Name")
            .build();

    private PluginParameter<String> myB4RUser = new PluginParameter.Builder<>("b4ruser", null, String.class)
            .required(true)
            .guiName("B4R User")
            .description("B4R User Name")
            .build();

    private PluginParameter<String> myB4RPassword = new PluginParameter.Builder<>("b4rpassword", "", String.class)
            .guiName("B4R Password")
            .description("Password")
            .password()
            .build();
    
    
    private PluginParameter<String> myStudyName = new PluginParameter.Builder<>("studyName","",
            String.class)
                        .description("Name of the study")
                        .guiName("Study Name")
                        .build();
    private PluginParameter<String> myTaxaName = new PluginParameter.Builder<>("taxaName","",
            String.class)
                        .description("Name of the taxa")
                        .guiName("Taxa Name")
                        .build();
    private PluginParameter<String> myVariableName = new PluginParameter.Builder<>("varName","",
            String.class)
                        .description("Name of the variable")
                        .guiName("Variable Name")
                        .build();
    
    public GetPhenotypeFromB4RPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    public DataSet processData(DataSet input) {
        //do whatever your plugin does
        Object result = null;
        try {
        result = B4RPhenotypeUtils.getPhenotypeFromB4RWithDuplicates(b4RHost(), b4RDBName(), b4RUser(), b4RPassword(), studyName(), taxaName(), variableName());
        
        }
        catch(Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Could not create Phenotype Object");
        }
        if (result != null) {
            
            ArrayList<String> outputNameList = new ArrayList<String>();
            outputNameList.add("B4R_Phenotype");
            if(!myStudyName.value().equals("")) {
                outputNameList.add("Study_"+myStudyName.value());
            }
            if(!myTaxaName.value().equals("")) {
                outputNameList.add("Taxa_"+myTaxaName.value());
            }
            if(!myVariableName.value().equals("")) {
                outputNameList.add("Var_"+myVariableName.value());
            }
            
            String name = outputNameList.stream().collect(Collectors.joining("_"));

            Datum td = new Datum(name, result, null);
            
            //todo need to add logic of directories.
            DataSet tds = new DataSet(td, this);
            return tds;
        }
        return null;  // Note: this can return null
    }
    
    @Override
    public ImageIcon getIcon() {
        return null;
    }
    @Override
    public String getButtonName() {
        return "Phenotype From B4R";
    }
    @Override
    public String getToolTipText() {
        return "Pull a Study from B4R";
    }
    
    @Override
    public String getCitation() {
        return "Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES. "
                + "(2007) TASSEL: Software for association mapping of complex traits in diverse "
                + "samples. Bioinformatics 23:2633­2635.";
    }
    
    
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(GetPhenotypeFromB4RPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public Phenotype runPlugin(DataSet input) {
        return (Phenotype) performFunction(input).getData(0).getData();
    }

    /**
     * B4R Database Name
     *
     * @return B4R Database Name
     */
    public String b4RDBName() {
        return myB4RDBName.value();
    }

    /**
     * Set B4R Database Name. B4R Database Name
     *
     * @param value B4R Database Name
     *
     * @return this plugin
     */
    public GetPhenotypeFromB4RPlugin b4RDBName(String value) {
        myB4RDBName = new PluginParameter<>(myB4RDBName, value);
        return this;
    }

    /**
     * B4R Host Name
     *
     * @return B4R Host Name
     */
    public String b4RHost() {
        return myB4RHost.value();
    }

    /**
     * Set B4R Host Name. B4R Host Name
     *
     * @param value B4R Host Name
     *
     * @return this plugin
     */
    public GetPhenotypeFromB4RPlugin b4RHost(String value) {
        myB4RHost = new PluginParameter<>(myB4RHost, value);
        return this;
    }

    /**
     * B4R User Name
     *
     * @return B4R User
     */
    public String b4RUser() {
        return myB4RUser.value();
    }

    /**
     * Set B4R User. B4R User Name
     *
     * @param value B4R User
     *
     * @return this plugin
     */
    public GetPhenotypeFromB4RPlugin b4RUser(String value) {
        myB4RUser = new PluginParameter<>(myB4RUser, value);
        return this;
    }

    /**
     * Password
     *
     * @return B4R Password
     */
    public String b4RPassword() {
        return myB4RPassword.value();
    }

    /**
     * Set B4R Password. Password
     *
     * @param value B4R Password
     *
     * @return this plugin
     */
    public GetPhenotypeFromB4RPlugin b4RPassword(String value) {
        myB4RPassword = new PluginParameter<>(myB4RPassword, value);
        return this;
    }

    /**
     * Name of the study
     *
     * @return Study Name
     */
    public String studyName() {
        return myStudyName.value();
    }

    /**
     * Set Study Name. Name of the study
     *
     * @param value Study Name
     *
     * @return this plugin
     */
    public GetPhenotypeFromB4RPlugin studyName(String value) {
        myStudyName = new PluginParameter<>(myStudyName, value);
        return this;
    }

    /**
     * Name of the taxa
     *
     * @return Taxa Name
     */
    public String taxaName() {
        return myTaxaName.value();
    }

    /**
     * Set Taxa Name. Name of the taxa
     *
     * @param value Taxa Name
     *
     * @return this plugin
     */
    public GetPhenotypeFromB4RPlugin taxaName(String value) {
        myTaxaName = new PluginParameter<>(myTaxaName, value);
        return this;
    }

    /**
     * Name of the variable
     *
     * @return Variable Name
     */
    public String variableName() {
        return myVariableName.value();
    }

    /**
     * Set Variable Name. Name of the variable
     *
     * @param value Variable Name
     *
     * @return this plugin
     */
    public GetPhenotypeFromB4RPlugin variableName(String value) {
        myVariableName = new PluginParameter<>(myVariableName, value);
        return this;
    }

//
//    
//    public static void main(String[] args) {
//        GeneratePluginCode.generate(GetPhenotypeFromB4RPlugin.class);
//    }
}
