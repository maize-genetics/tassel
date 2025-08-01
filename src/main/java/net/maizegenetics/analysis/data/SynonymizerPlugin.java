/*
 * SynonymizerPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.taxa.IdentifierSynonymizer;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.gui.TableReportNoPagingTableModel;
import net.maizegenetics.tassel.TASSELMainFrame;
import net.maizegenetics.util.GeneralAnnotation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import javax.swing.table.TableModel;
import javax.swing.tree.DefaultMutableTreeNode;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URL;
import java.text.Collator;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author Ed Buckler, Zack Miller
 */
public class SynonymizerPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(SynonymizerPlugin.class);
    private String delimiter;
    /**
     * Creates a new instance of SynonymizerPlugin
     */
    public SynonymizerPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {
        try {
            List<Datum> data = new ArrayList<Datum>();
            for (int i = 0, n = input.getSize(); i < n; i++) {
                Datum current = input.getData(i);
                Object currentData = current.getData();
                if (currentData instanceof GenotypeTable) {
                    TaxaList idGroup = ((GenotypeTable) currentData).taxa();
                    Datum idGroupDatum = new Datum(current.getName(), idGroup, current.getComment());
                    data.add(idGroupDatum);
                } else if (currentData instanceof Phenotype) {
                    TaxaList idGroup = ((Phenotype) currentData).taxa();
                    Datum idGroupDatum = new Datum(current.getName(), idGroup, current.getComment());
                    data.add(idGroupDatum);
                } else {
                    data.add(current);
                }
            }
            DataSet newInput = new DataSet(data, this);
            if(isInteractive()) {
                
                TASSELMainFrame frame = (TASSELMainFrame)getParentFrame();
                Map map = frame.getDataTreePanel().getDataList();
                Object[] datumArray = map.keySet().toArray(); 
                
                //check to see if the user has loaded at least 2 files
                if(datumArray.length<2) {
                    String msg = "Error:  Make sure at least 2 files are loaded into TASSEL before attempting to Synonymize.";
                    JOptionPane.showMessageDialog(getParentFrame(), msg);
                }
                else {
                    int response = JOptionPane.showConfirmDialog(frame,
                                    "Would you like to run the full Synonymizer Pipeline?",
                                    "Full Synonymize",
                                    JOptionPane.YES_NO_OPTION);
                   
                    if(response==0){
                        int[] fileOptions = getFileChoice(data,datumArray);
                        DataSet step1Result = runStep1(data,datumArray,fileOptions);
                        Datum step1DatumResult = (Datum)((ArrayList)step1Result.getDataSet()).get(0);
                        runStep2(step1DatumResult,true);
                        runStep3(datumArray,step1DatumResult,fileOptions);
                    }
                    else {
                        boolean[] menuArray = {false,false,false};
                        SynMenuDialog menuDiag = new SynMenuDialog(menuArray,getParentFrame());
                        menuDiag.setLocationRelativeTo(getParentFrame());
                        menuDiag.setVisible(true);
                        
                        if(menuArray[0] == true) {
                            int[] fileOptions = getFileChoice(data,datumArray);
                            return runStep1(data,datumArray,fileOptions);
                        }
                        else if(menuArray[1] == true) {
                            //Do a quick check to make sure there is at least one Synonym loaded
                            boolean hasSynonymFile = false;
                            for(int i = 0;i<datumArray.length;i++) {
                                Datum currentDatum = (Datum)datumArray[i];
                                if(currentDatum.getDataType().equals(IdentifierSynonymizer.class)) {
                                    hasSynonymFile = true;
                                }
                            }
                            if(!hasSynonymFile) {
                                String msg = "Error:  No Synonymize files have been found.  Make sure at least one Synonymize File and at least one Standard File are loaded."
                                        + "\nPlease run the first step.";
                                JOptionPane.showMessageDialog(getParentFrame(), msg);
                            }
                            else {
                              //Check to see if user has selected files like before
                                String[] initialSelections = new String[1];
                                boolean errorOut = false;
                                boolean validSelection = false;
                                //If there is a selection
                                if(data.size()>0) {
                                    //Set comboboxes to match
                                    if(data.get(0).getDataType().equals(IdentifierSynonymizer.class)) {
                                        initialSelections[0] = data.get(0).getName();
                                        validSelection = true;
                                    }
                                }
                                
                                if(!validSelection) {
                                  //First file which is an IdentifierSynonymizer
                                    initialSelections[0] = "";
                                    for(int i = 0; i<datumArray.length;i++) {
                                        Datum currentDatum = (Datum)datumArray[i];
                                        if(currentDatum.getDataType().equals(IdentifierSynonymizer.class)) {
                                            initialSelections[0] = currentDatum.getName();
                                            break;
                                        }
                                    }
                                    if(initialSelections[0].equals("")) {
                                        String msg = "Error:  No Synonymize files have been found.  Please run the first step.";
                                        JOptionPane.showMessageDialog(getParentFrame(), msg);
                                        errorOut = true;
                                    }
                                
                                }
                               
                                if(!errorOut) {
                                    int[] fileOptions = new int[1];
                                    fileOptions[0] = -1;
                                    SynonymizerFileChooser fileChooseDiag= new SynonymizerFileChooser(getParentFrame(),datumArray,fileOptions,initialSelections,"Step2");
                                    fileChooseDiag.setLocationRelativeTo(getParentFrame());
                                    fileChooseDiag.setVisible(true);
                                    if(fileOptions[0]!=-1) {
                                        Datum current = (Datum)datumArray[fileOptions[0]];
                                        IdentifierSynonymizer is = (IdentifierSynonymizer) current.getData();
                                        SynonymizerDialog theSD = new SynonymizerDialog(is, getParentFrame());
                                        theSD.setLocationRelativeTo(getParentFrame());
                                        theSD.setVisible(true);
                                    }
                                }
                            }                    
                        }
                        else if(menuArray[2]==true) {
                            //Do a quick check to make sure there is at least one Synonym loaded
                            boolean hasSynonymFile = false;
                            for(int i = 0;i<datumArray.length;i++) {
                                Datum currentDatum = (Datum)datumArray[i];
                                if(currentDatum.getDataType().equals(IdentifierSynonymizer.class)) {
                                    hasSynonymFile = true;
                                }
                            }
                            if(!hasSynonymFile) {
                                String msg = "Error:  No Synonymize files have been found.  Make sure at least one Synonymize File and at least one Standard File are loaded."
                                        + "\nPlease run the first step.";
                                JOptionPane.showMessageDialog(getParentFrame(), msg);
                            }
                            else {
                              //Check to see if user has selected files like before
                                //initSelections[0] is the SynonymizerFile
                                //initSelections[1] is the File to be Synonymized
                                String[] initialSelections = new String[2];
                                boolean errorOut = false;
                                boolean validSelection = false;
                                //If there is a selection
                                if(data.size()>1) {
                                    //Set comboboxes to match
                                    if(data.get(0).getDataType().equals(IdentifierSynonymizer.class)) {
                                        initialSelections[0] = data.get(0).getName();
                                        initialSelections[1] = data.get(1).getName();
                                        validSelection = true;
                                    }
                                }
                                if(!validSelection) {
                                  //First file which is an IdentifierSynonymizer
                                    initialSelections[0] = "";
                                    initialSelections[1] = "";
                                    for(int i = 0; i<datumArray.length;i++) {
                                        Datum currentDatum = (Datum)datumArray[i];
                                        if(initialSelections[1].equals("") && !currentDatum.getDataType().equals(IdentifierSynonymizer.class)) {
                                            initialSelections[1] = currentDatum.getName();
                                        }
                                        if(currentDatum.getDataType().equals(IdentifierSynonymizer.class)) {
                                            initialSelections[0] = currentDatum.getName();
                                            break;
                                        }
                                    }
                                    if(initialSelections[0].equals("")) {
                                        String msg = "Error:  No Synonymize files have been found.  Please run the first step.";
                                        JOptionPane.showMessageDialog(getParentFrame(), msg);
                                        errorOut = true;
                                    }
                                }
                               
                                if(!errorOut) {
                                    int[] fileOptions = new int[2];
                                    SynonymizerFileChooser fileChooseDiag= new SynonymizerFileChooser(getParentFrame(),datumArray,fileOptions,initialSelections,"Step3");
                                    fileChooseDiag.setLocationRelativeTo(getParentFrame());
                                    fileChooseDiag.setVisible(true);
                                    
                                    if(fileOptions[0]!=-1) {
                                        Datum firstDatum = (Datum)datumArray[fileOptions[0]];
                                        IdentifierSynonymizer first = (IdentifierSynonymizer)firstDatum.getData();
                                        
                                        //Figure out what original Data Structure is
                                        Datum current = (Datum)datumArray[fileOptions[1]];
                                        
                                        Object currentData = current.getData();
                                        
                                        if (currentData instanceof GenotypeTable) {
                                            //Build a new Genotype table
                                            //TaxaList idGroup = ((GenotypeTable) currentData).taxa();
                                            //Datum idGroupDatum = new Datum(current.getName(), idGroup, current.getComment());
                                            if(first.checkSynForDups()) {
                                                String msg = "Error:  Duplicate Taxa to Be Synonymized.  Please Manually Update the Synonym List.";
                                                JOptionPane.showMessageDialog(getParentFrame(), msg);
                                                errorOut = true; 
                                            }
                                            else {
                                                //First tell Synonymizer to swap synonyms with real names
                                                ArrayList<TaxaList> swapped = first.swapSynonyms();
                                                
                                                GenotypeTable currentGenotypeTable = (GenotypeTable)currentData;
                                                GenotypeTable a=GenotypeTableBuilder.getInstance(currentGenotypeTable.genotypeMatrix(), 
                                                                                                currentGenotypeTable.positions(), swapped.get(0)); 
                                                Datum idGroupDatum = new Datum(current.getName()+"_Synonymized",a,current.getComment());
                                                ArrayList<Datum> listOfDatum = new ArrayList<>();
                                                listOfDatum.add(idGroupDatum);
                                                //Fire off new genotype table
                                                DataSet output = new DataSet(listOfDatum, this);
                                                fireDataSetReturned(new PluginEvent(output, SynonymizerPlugin.class));
                                            }
                                        } else if (currentData instanceof Phenotype) {
                                            Phenotype currentPhenotype = (Phenotype)currentData;
                                            //Get all information in Phenotype
                                            List<PhenotypeAttribute> attributeList = currentPhenotype.attributeListCopy();
                                            ArrayList<PhenotypeAttribute> newAttributeList = new ArrayList<>();
                                            ArrayList<Taxon> listOfTaxon = new ArrayList<Taxon>();
                                            //Figure out the Types of each attribute
                                            ArrayList<Phenotype.ATTRIBUTE_TYPE> types = new ArrayList<>();
                                            for(int i = 0; i<attributeList.size();i++) {
                                                types.add(currentPhenotype.attributeType(i));
                                                if(currentPhenotype.attributeType(i).equals(Phenotype.ATTRIBUTE_TYPE.taxa)) {
                                                    //For the Taxa attribute go through and update the names
                                                    List<Taxon> oldTaxa = currentPhenotype.taxaAttribute().allTaxaAsList();
                                                    TaxaList synTl = first.getTaxaList();
                                                    HashMap<String,Taxon> taxonMap = new HashMap<String,Taxon>();
                                                    for(Taxon taxon:synTl) {
                                                        taxonMap.put(taxon.getName(),taxon);
                                                    }
                                                    for(Taxon taxon:oldTaxa) {
                                                        //if taxonMap has the taxonName as a key
                                                        if(taxonMap.containsKey(taxon.getName())){
                                                            Taxon mappedTaxon = taxonMap.get(taxon.getName());
                                                            //Check the SynonymList
                                                            GeneralAnnotation ga = mappedTaxon.getAnnotation();
                                                            String[] values = ga.getTextAnnotation(Taxon.SynonymKey);
                                                            
                                                            if(values.length==0 || values==null) {
                                                                //Add taxon to listOfTaxon
                                                                listOfTaxon.add(taxon);
                                                            }
                                                            else {
                                                                //else make a new Taxon
                                                                Taxon.Builder tb = new Taxon.Builder(ga.getTextAnnotation(Taxon.SynonymKey)[0]);
                                                                String[] keys = ga.getTextAnnotation(Taxon.SynonymKey);
                                                                //Copy Keys
                                                                for(String key:keys) {
                                                                    if(!key.equals(Taxon.SynonymKey)) {
                                                                        String[] values2 = ga.getTextAnnotation(key);
                                                                        for(String value:values2) {
                                                                            tb.addAnno(key, value);
                                                                        }
                                                                    }
                                                                }
                                                                String[] synVals = ga.getTextAnnotation(Taxon.SynonymKey);
                                                                tb.addAnno(Taxon.SynonymKey, taxon.getName());
                                                                for(int j = 1; j<synVals.length;j++) {
                                                                    tb.addAnno(Taxon.SynonymKey, synVals[j]);
                                                                }
                                                                listOfTaxon.add(tb.build());
                                                            }
                                                        } 
                                                        else {
                                                            //Add taxon to listOfTaxon
                                                            listOfTaxon.add(taxon);
                                                        }
                                                    }
                                                    //Add the new Taxon to the the type structure
                                                    newAttributeList.add(new TaxaAttribute(listOfTaxon));
                                                }
                                                else {
                                                    newAttributeList.add(attributeList.get(i));
                                                }
                                            }
                                            //Use public PhenotypeBuilder fromAttributeList(List<PhenotypeAttribute> attributes, List<ATTRIBUTE_TYPE> types) 
                                            PhenotypeBuilder phb = new PhenotypeBuilder();
                                            phb.fromAttributeList(newAttributeList,types);
                                            //Build Phenotype
                                            Phenotype buildPhenotype = phb.build().get(0);
                                            
                                            Datum idGroupDatum = new Datum(current.getName()+"_Synonymized",buildPhenotype,current.getComment());
                                            ArrayList<Datum> listOfDatum = new ArrayList<>();
                                            listOfDatum.add(idGroupDatum);
                                            //Fire off new genotype table
                                            DataSet output = new DataSet(listOfDatum, this);
                                            fireDataSetReturned(new PluginEvent(output, SynonymizerPlugin.class));
                                           
                                        } else {
                                            //Dont do anything
                                            String msg = "Error:  Unsupported File Type. Please use a GenotypeTable or Phenotype.";
                                            JOptionPane.showMessageDialog(getParentFrame(), msg);
                                            errorOut = true; 
                                        }
                                    }
                                }                          
                            }
                        }
                    }
                }
                return null;
            }
            else {
                int alignCnt = newInput.getDataOfType(TaxaList.class).size();
                int synCnt = newInput.getDataOfType(IdentifierSynonymizer.class).size();
                if ((synCnt == 0) && (alignCnt > 1)) {  //create a new synonymizer Step 1
                    Datum td = createSynonymizer(newInput);
                    DataSet output = new DataSet(td, this);
                    fireDataSetReturned(new PluginEvent(output, SynonymizerPlugin.class));
                    return output;
                } else if ((synCnt == 1) && (alignCnt > 0)) {   //apply synonymizer to alignments Step 3
                    //applySynonymsToIdGroups(newInput);
                } else if ((synCnt == 1) && (alignCnt == 0)) {  //Step 2
                    if (isInteractive()) {
                        Datum inputDatum = newInput.getDataOfType(IdentifierSynonymizer.class).get(0);
                        IdentifierSynonymizer is = (IdentifierSynonymizer) inputDatum.getData();
                        SynonymizerDialog theSD = new SynonymizerDialog(is, getParentFrame());
                        theSD.setLocationRelativeTo(getParentFrame());
                        theSD.setVisible(true);
                    }
                } else {
                    String msg = "To create a synonym list:\n Please first select the reference taxa names and then the synonym taxa names (use Ctrl key)\n"
                            + "To apply a synonym list to a dataset:\n Select a synonym list and then the taxa names to be changed (use Ctrl key)";
                    if (isInteractive()) {
                        JOptionPane.showMessageDialog(getParentFrame(), msg);
                    } else {
                        myLogger.error(msg);
                    }
                }
    
                return null;
            }
        } finally {
            fireProgress(100);
        }
    }

    private Datum createSynonymizer(DataSet input) {
        Datum td = null;
        StringBuilder synonymSets = new StringBuilder();
        for (int i = 1; i < input.getSize(); i++) {
            synonymSets.append(input.getData(i).getName());
            synonymSets.append("\n");
        }
        boolean performFunction = true;
        //String msg = "You have selected to apply synonym list " + input.getData(0).getName() + " to the following dataset:\n"
        //        + synonymSets.toString();
        String msg = "You have selected to generate a synonym list from " + input.getData(0).getName() + " to be applied to the following dataset:\n"
                + synonymSets.toString();
        
        if (isInteractive()) {
            int response = JOptionPane.showOptionDialog(getParentFrame(), msg, "Verify Selection",
                    JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
            if (response == JOptionPane.CANCEL_OPTION) {
                performFunction = false;
            }
        } else {
            myLogger.info(msg);
        }
        if (performFunction) {
            List<Datum> idList = input.getDataOfType(TaxaList.class);
            TaxaList[] aa = new TaxaList[idList.size() - 1];
            for (int i = 1; i < idList.size(); i++) {
                aa[i - 1] = (TaxaList) idList.get(i).getData();
            }
            IdentifierSynonymizer ts = new IdentifierSynonymizer((TaxaList) idList.get(0).getData(), aa);
            StringWriter sw = new StringWriter();
            ts.report(new PrintWriter(sw));
            td = new Datum(input.getData(0).getName() + " Synonyms", ts, "Taxa synonyms\n" + sw.toString());
        }
        return td;
    }
    
    private Datum createSynonymizer(DataSet input,int technique) {
        Datum td = null;
        StringBuilder synonymSets = new StringBuilder();
        for (int i = 1; i < input.getSize(); i++) {
            synonymSets.append(input.getData(i).getName());
            synonymSets.append("\n");
        }
        boolean performFunction = true;
        //String msg = "You have selected to apply synonym list " + input.getData(0).getName() + " to the following dataset:\n"
        //        + synonymSets.toString();
        String msg = "You have selected to generate a synonym list from " + input.getData(0).getName() + " to be applied to the following dataset:\n"
                + synonymSets.toString();
        
        if (isInteractive()) {
            int response = JOptionPane.showOptionDialog(getParentFrame(), msg, "Verify Selection",
                    JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
            if (response == JOptionPane.CANCEL_OPTION) {
                performFunction = false;
            }
        } else {
            myLogger.info(msg);
        }
        if (performFunction) {
            List<Datum> idList = input.getDataOfType(TaxaList.class);
            TaxaList[] aa = new TaxaList[idList.size() - 1];
            for (int i = 1; i < idList.size(); i++) {
                aa[i - 1] = (TaxaList) idList.get(i).getData();
            }
            IdentifierSynonymizer ts = new IdentifierSynonymizer((TaxaList) idList.get(0).getData(), aa,technique);
            StringWriter sw = new StringWriter();
            ts.report(new PrintWriter(sw));
            //td = new Datum(input.getData(0).getName() + " Synonyms", ts, "Taxa synonyms\n" + sw.toString());
            td = new Datum(input.getData(1).getName() + " Synonymizer Object", ts, "Taxa synonyms\n"+sw.toString());
        }
        return td;
    }
    private Datum createSynonymizer(DataSet input,int technique,String delimiter) {
        Datum td = null;
        StringBuilder synonymSets = new StringBuilder();
        for (int i = 1; i < input.getSize(); i++) {
            synonymSets.append(input.getData(i).getName());
            synonymSets.append("\n");
        }
        boolean performFunction = true;
        //String msg = "You have selected to apply synonym list " + input.getData(0).getName() + " to the following dataset:\n"
        //        + synonymSets.toString();
        String msg = "You have selected to generate a synonym list from " + input.getData(0).getName() + " to be applied to the following dataset:\n"
                + synonymSets.toString();
        
        if (isInteractive()) {
            int response = JOptionPane.showOptionDialog(getParentFrame(), msg, "Verify Selection",
                    JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
            if (response == JOptionPane.CANCEL_OPTION) {
                performFunction = false;
            }
        } else {
            myLogger.info(msg);
        }
        if (performFunction) {
            List<Datum> idList = input.getDataOfType(TaxaList.class);
            TaxaList[] aa = new TaxaList[idList.size() - 1];
            for (int i = 1; i < idList.size(); i++) {
                aa[i - 1] = (TaxaList) idList.get(i).getData();
            }
            IdentifierSynonymizer ts = new IdentifierSynonymizer((TaxaList) idList.get(0).getData(), aa,technique,delimiter);
            StringWriter sw = new StringWriter();
            ts.report(new PrintWriter(sw));
            //td = new Datum(input.getData(0).getName() + " Synonyms", ts, "Taxa synonyms\n" + sw.toString());
            td = new Datum(input.getData(1).getName() + " Synonymizer Object", ts, "Taxa synonyms\n"+sw.toString());
        }
        return td;
    }

    /*
    private void applySynonymsToIdGroups(DataSet input) {
        StringBuilder synonymSets = new StringBuilder();
        for (int i = 1; i < input.getSize(); i++) {
            synonymSets.append(input.getData(i).getName());
            synonymSets.append("\n");
        }
        boolean performFunction = true;
        //String msg = "You have selected " + input.getData(0).getName() + " as the reference name dataset.\n"
        //        + "The synonyms will be extracted from the following: \n" + synonymSets.toString();
        String msg = "You have selected " + input.getData(0).getName() + " as the list of Synonyms.\n"
                + "The synonyms will be written to the following: \n" + synonymSets.toString();
        if (isInteractive()) {
            int response = JOptionPane.showOptionDialog(getParentFrame(), msg, "Verify Selection",
                    JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
            if (response == JOptionPane.CANCEL_OPTION) {
                performFunction = false;
            }
        } else {
            myLogger.info(msg);
        }
        if (performFunction) {
            IdentifierSynonymizer is = (IdentifierSynonymizer) input.getDataOfType(IdentifierSynonymizer.class).get(0).getData();
            List<Datum> idList = input.getDataOfType(TaxaList.class);
            TaxaList[] aa = new TaxaList[idList.size()];
            for (int i = 0; i < idList.size(); i++) {
                aa[i] = (TaxaList) idList.get(i).getData();
            }
            is.changeAlignmentIdentifiers(aa);
        }
    }
    */

    private void applySynonymsToOldFile(Datum current, IdentifierSynonymizer first){
        Object currentData = current.getData();
        if (currentData instanceof GenotypeTable) {
            //Build a new Genotype table
            //TaxaList idGroup = ((GenotypeTable) currentData).taxa();
            //Datum idGroupDatum = new Datum(current.getName(), idGroup, current.getComment());
            if(first.checkSynForDups()) {
                String msg = "Error:  Duplicate Taxa to Be Synonymized.  Please Manually Update the Synonym List.";
                JOptionPane.showMessageDialog(getParentFrame(), msg);
            }
            else {
                //First tell Synonymizer to swap synonyms with real names
                ArrayList<TaxaList> swapped = first.swapSynonyms();
                
                GenotypeTable currentGenotypeTable = (GenotypeTable)currentData;
                GenotypeTable a=GenotypeTableBuilder.getInstance(currentGenotypeTable.genotypeMatrix(), 
                                                                currentGenotypeTable.positions(), swapped.get(0)); 
                Datum idGroupDatum = new Datum(current.getName()+"_Synonymized",a,current.getComment());
                ArrayList<Datum> listOfDatum = new ArrayList<>();
                listOfDatum.add(idGroupDatum);
                //Fire off new genotype table
                DataSet output = new DataSet(listOfDatum, this);
                fireDataSetReturned(new PluginEvent(output, SynonymizerPlugin.class));
            }
        } else if (currentData instanceof Phenotype) {
            Phenotype currentPhenotype = (Phenotype)currentData;
            //Get all information in Phenotype
            List<PhenotypeAttribute> attributeList = currentPhenotype.attributeListCopy();
            ArrayList<PhenotypeAttribute> newAttributeList = new ArrayList<>();
            ArrayList<Taxon> listOfTaxon = new ArrayList<Taxon>();
            //Figure out the Types of each attribute
            ArrayList<Phenotype.ATTRIBUTE_TYPE> types = new ArrayList<>();
            for(int i = 0; i<attributeList.size();i++) {
                types.add(currentPhenotype.attributeType(i));
                if(currentPhenotype.attributeType(i).equals(Phenotype.ATTRIBUTE_TYPE.taxa)) {
                    //For the Taxa attribute go through and update the names
                    List<Taxon> oldTaxa = currentPhenotype.taxaAttribute().allTaxaAsList();
                    TaxaList synTl = first.getTaxaList();
                    HashMap<String,Taxon> taxonMap = new HashMap<String,Taxon>();
                    for(Taxon taxon:synTl) {
                        taxonMap.put(taxon.getName(),taxon);
                    }
                    for(Taxon taxon:oldTaxa) {
                        //if taxonMap has the taxonName as a key
                        if(taxonMap.containsKey(taxon.getName())){
                            Taxon mappedTaxon = taxonMap.get(taxon.getName());
                            //Check the SynonymList
                            GeneralAnnotation ga = mappedTaxon.getAnnotation();
                            String[] values = ga.getTextAnnotation(Taxon.SynonymKey);
                            
                            if(values.length==0 || values==null) {
                                //Add taxon to listOfTaxon
                                listOfTaxon.add(taxon);
                            }
                            else {
                                //else make a new Taxon
                                Taxon.Builder tb = new Taxon.Builder(ga.getTextAnnotation(Taxon.SynonymKey)[0]);
                                String[] keys = ga.getTextAnnotation(Taxon.SynonymKey);
                                //Copy Keys
                                for(String key:keys) {
                                    if(!key.equals(Taxon.SynonymKey)) {
                                        String[] values2 = ga.getTextAnnotation(key);
                                        for(String value:values2) {
                                            tb.addAnno(key, value);
                                        }
                                    }
                                }
                                String[] synVals = ga.getTextAnnotation(Taxon.SynonymKey);
                                tb.addAnno(Taxon.SynonymKey, taxon.getName());
                                for(int j = 1; j<synVals.length;j++) {
                                    tb.addAnno(Taxon.SynonymKey, synVals[j]);
                                }
                                listOfTaxon.add(tb.build());
                            }
                        } 
                        else {
                            //Add taxon to listOfTaxon
                            listOfTaxon.add(taxon);
                        }
                    }
                    //Add the new Taxon to the the type structure
                    newAttributeList.add(new TaxaAttribute(listOfTaxon));
                }
                else {
                    newAttributeList.add(attributeList.get(i));
                }
            }
            //Use public PhenotypeBuilder fromAttributeList(List<PhenotypeAttribute> attributes, List<ATTRIBUTE_TYPE> types) 
            PhenotypeBuilder phb = new PhenotypeBuilder();
            phb.fromAttributeList(newAttributeList,types);
            //Build Phenotype
            Phenotype buildPhenotype = phb.build().get(0);
            //Build a new Phenotype object
            
            /*
            TaxaList idGroup = ((Phenotype) currentData).taxa();
            Datum idGroupDatum = new Datum(current.getName(), idGroup, current.getComment());
            */
            
            Datum idGroupDatum = new Datum(current.getName()+"_Synonymized",buildPhenotype,current.getComment());
            ArrayList<Datum> listOfDatum = new ArrayList<>();
            listOfDatum.add(idGroupDatum);
            //Fire off new genotype table
            DataSet output = new DataSet(listOfDatum, this);
            fireDataSetReturned(new PluginEvent(output, SynonymizerPlugin.class));
           
        } else {
            //Dont do anything
            String msg = "Error:  Unsupported File Type. Please use a GenotypeTable or Phenotype.";
            JOptionPane.showMessageDialog(getParentFrame(), msg);
        }     
    }
    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = SynonymizerPlugin.class.getResource("/net/maizegenetics/analysis/images/Synonymizer.gif");
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
    public String getButtonName() {
        return "Synonymizer";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Unify Taxa Names";
    }

    public int[] getFileChoice(List<Datum> data, Object[] datumArray) {
      //Check to see if user has selected files like before
        String[] initialSelections = new String[2];
        if(data.size()>1) {
            //Set comboboxes to match
            initialSelections[0] = data.get(0).getName();
            initialSelections[1] = data.get(1).getName();
        }
        else {
            //First and Second file
            Datum firstDatum = (Datum)datumArray[0];
            initialSelections[0] = firstDatum.getName();
            Datum secondDatum = (Datum)datumArray[1];
            initialSelections[1] = secondDatum.getName();
        }

        int[] fileOptions = new int[3];
        SynonymizerFileChooser fileChooseDiag= new SynonymizerFileChooser(getParentFrame(),datumArray,fileOptions,initialSelections,"Step1");
        fileChooseDiag.setLocationRelativeTo(getParentFrame());
        fileChooseDiag.setVisible(true);
        if(fileOptions[2]==7) {
            delimiter = fileChooseDiag.getDelim();
        }
        System.out.println("Delim: "+fileChooseDiag.getDelim());
        
        return fileOptions;
    }
    public DataSet runStep1(List<Datum> data, Object[] datumArray, int[] fileOptions ) {
        
        if(fileOptions[0]!=-1) {
            ArrayList<Datum> datumList = new ArrayList<Datum>();
            for(int i = 0; i<fileOptions.length-1;i++) {
                Datum current = (Datum)datumArray[fileOptions[i]];
                Object currentData = current.getData();
                if (currentData instanceof GenotypeTable) {
                    TaxaList idGroup = ((GenotypeTable) currentData).taxa();
                    Datum idGroupDatum = new Datum(current.getName(), idGroup, current.getComment());
                    datumList.add(idGroupDatum);
                } else if (currentData instanceof Phenotype) {
                    TaxaList idGroup = ((Phenotype) currentData).taxa();
                    Datum idGroupDatum = new Datum(current.getName(), idGroup, current.getComment());
                    datumList.add(idGroupDatum);
                } else {
                    datumList.add(current);
                }
            }
            DataSet newInputDataSet = new DataSet(datumList, this);
            Datum td = null;
            if(fileOptions[2] == 7) {
                td = createSynonymizer(newInputDataSet,fileOptions[2],delimiter);
            }
            else {
                td = createSynonymizer(newInputDataSet,fileOptions[2]);
            }
            //Datum td = createSynonymizer(newInputDataSet,fileOptions[2]);
            DataSet output = new DataSet(td, this);
            fireDataSetReturned(new PluginEvent(output, SynonymizerPlugin.class));
            return output;
        
        }
        else {
           return null;
        }
        
    }

    public void runStep2(Datum currentData, boolean end2end) {
        IdentifierSynonymizer is = (IdentifierSynonymizer) currentData.getData();
        SynonymizerDialog theSD = new SynonymizerDialog(is, getParentFrame());
        theSD.setLocationRelativeTo(getParentFrame());
        theSD.setVisible(true);
    }

    public void runStep3(Object[] data, Datum step1Result, int[] fileOptions) {
        IdentifierSynonymizer first = (IdentifierSynonymizer)step1Result.getData();
        Datum current = (Datum)data[fileOptions[1]];
        Object currentData = current.getData();
        applySynonymsToOldFile(current,first);
    }
    
}

/**
 */
class SynonymizerDialog extends JDialog {

    private JPanel jPanel1 = new JPanel();
    private JTextField ThresholdTextField = new JTextField();
    private JButton setThresholdButton = new JButton();
    private JButton CancelButton = new JButton();
    private Frame theFrame;
    private double threshold = 1.0;
    boolean isCanceled;
    JButton okButton = new JButton();
    JList matchList = new JList();
    JScrollPane newRealNameScrollPane1 = new JScrollPane();
    JTable synTable = new JTable();
    JButton selectSynButton = new JButton();
    JButton setNoSynButton = new JButton();
    JLabel jLabel1 = new JLabel();
    JPanel jPanel2 = new JPanel();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JScrollPane theATP;
    JTable theNameTable;
    IdentifierSynonymizer theTS;
    JCheckBox cbxSortAlphabetically = new JCheckBox();
    
    JLabel addTaxa = new JLabel("Add a new Synonym");
    JTextField addTaxaField = new JTextField();
    JButton addTaxaButton = new JButton();

    public SynonymizerDialog(IdentifierSynonymizer ts, Frame theFrame) {
        super((Frame) theFrame, true);
        this.theFrame = theFrame;
        this.theTS = ts;
        try {
            theNameTable = new JTable(new TableReportNoPagingTableModel(this.theTS));
            theNameTable.setAutoCreateRowSorter(true);
            theNameTable.setCellEditor(null);
            theNameTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
            
            matchList.setAutoscrolls(true);
            theATP = new JScrollPane(theNameTable);
            jbInit();
            pack();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public SynonymizerDialog() {
        try {
            jbInit();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void jbInit() throws Exception {
        jPanel1.setLayout(gridBagLayout1);
        setThresholdButton.setText("Apply threshold");
        setThresholdButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                setThresholdButton_actionPerformed(e);
            }
        });
        CancelButton.setText("Cancel");
        CancelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                CancelButton_actionPerformed(e);
            }
        });
        okButton.setText("OK");
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                okButton_actionPerformed(e);
            }
        });
        ThresholdTextField.setPreferredSize(new Dimension(30, 30));
        selectSynButton.setFont(new java.awt.Font("Dialog", Font.BOLD, 14));
        selectSynButton.setToolTipText("Set synonym to selected taxon");
        selectSynButton.setText("<");
        selectSynButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                selectSynButton_actionPerformed(e);
            }
        });
        // setNoSynButton.setFont(new java.awt.Font("Dialog", Font.BOLD, 14));
        setNoSynButton.setToolTipText("Set selected taxon to no synonym");
        setNoSynButton.setText("No Synonym");
        setNoSynButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                setNoSynButton_actionPerformed(e);
            }
        });
        theNameTable.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                synTable_mouseClicked(e);
            }
        });
        jLabel1.setFont(new java.awt.Font("Dialog", Font.BOLD, 14));
        jLabel1.setText("Synonymizer");
        newRealNameScrollPane1.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        cbxSortAlphabetically.setActionCommand("jCheckBox1");
        cbxSortAlphabetically.setHorizontalAlignment(SwingConstants.CENTER);
        cbxSortAlphabetically.setHorizontalTextPosition(SwingConstants.TRAILING);
        cbxSortAlphabetically.setSelectedIcon(null);
        cbxSortAlphabetically.setText("Sort Alphabetically");
        cbxSortAlphabetically.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cbxSortAlphabetically_actionPerformed(e);
            }
        });

        addTaxaButton.setText("Add Synonym");
        addTaxaField.setPreferredSize(new Dimension(60, 30));
        
        addTaxaButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                addTaxaButton_actionPerformed(e);
            }
        });
        
        this.getContentPane().add(jPanel1, BorderLayout.CENTER);

        jPanel2.add(theATP);
        newRealNameScrollPane1.getViewport().add(matchList);
        jPanel1.add(ThresholdTextField, new GridBagConstraints(1, 2, 1, 2, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 50, 0));
        jPanel1.add(setThresholdButton, new GridBagConstraints(2, 3, 3, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 10, 7));
        jPanel1.add(newRealNameScrollPane1, new GridBagConstraints(2, 1, 2, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 100, 300));
        jPanel1.add(selectSynButton, new GridBagConstraints(1, 1, 1, 1, 0.5, 0.5, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 10, 4));
        jPanel1.add(setNoSynButton, new GridBagConstraints(0, 2, 1, 2, 0.5, 0.5, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 10, 4));
        jPanel1.add(CancelButton, new GridBagConstraints(2, 5, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(19, 0, 10, 25), 15, 7));
        jPanel1.add(jLabel1, new GridBagConstraints(0, 0, 3, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(8, 20, 0, 0), 184, 5));
        jPanel1.add(theATP, new GridBagConstraints(0, 1, 1, 1, 2.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 200, 300));
        jPanel1.add(okButton, new GridBagConstraints(0, 5, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(19, 69, 10, 49), 33, 10));
        jPanel1.add(addTaxa, new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0, GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0,0,0,0), 33, 10));
        jPanel1.add(addTaxaField, new GridBagConstraints(1, 4, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 70, 7));
        jPanel1.add(addTaxaButton, new GridBagConstraints(2, 4, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0,0,0,0), 37, 7));
        
        jPanel1.add(cbxSortAlphabetically, new GridBagConstraints(3, 2, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTHWEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        this.setSize(800, 600);
        this.setTitle(" Threshold for synonymizer ");
    }

    public IdentifierSynonymizer getIdentifierSynonymizer() {
        return theTS;
    }

    void deleteByThreshold(double threshold) {
        theTS.deleteByThreshold(threshold);
        TableReportNoPagingTableModel model = (TableReportNoPagingTableModel)theNameTable.getModel();
        model.fireTableChanged();        
    }

    void setThresholdButton_actionPerformed(ActionEvent e) {
        threshold = getMatchThreshold();
        deleteByThreshold(threshold);
    }

    void CancelButton_actionPerformed(ActionEvent e) {
        theTS.resetTempTaxaList();
        isCanceled = true;
        this.setVisible(false);
    }

    double getMatchThreshold() {
        double th = threshold;
        try {
            th = Double.parseDouble(ThresholdTextField.getText().trim());
            if ((th < 0) || (th > 1)) {
                throw new NumberFormatException();
            }
        } catch (NumberFormatException nfe) {
            JOptionPane.showMessageDialog(theFrame, "Please enter an double between 0 and 1.");
        }
        return th;
    }

    public boolean isCanceled() {
        return isCanceled;
    }

    void selectSynButton_actionPerformed(ActionEvent e) {
        try {
            String newRealName = (String) matchList.getSelectedValue();
            int theRow = theNameTable.getSelectedRow();
            theTS.updateSynonym(theRow,newRealName);
            TableReportNoPagingTableModel model = (TableReportNoPagingTableModel)theNameTable.getModel();
            model.fireTableChanged();
        } catch (Exception ex) {
            System.out.println("Make sure both a row and a new name are selected");
        }
    }
    void addTaxaButton_actionPerformed(ActionEvent e) {
        try {
            String newTaxa = addTaxaField.getText();
            int theRow = theNameTable.getSelectedRow();
            theTS.updateSynonym(theRow,newTaxa);
            TableReportNoPagingTableModel model = (TableReportNoPagingTableModel)theNameTable.getModel();
            model.fireTableChanged();
        } catch (Exception ex) {
            System.out.println("Make sure both a row and a new name are selected");
        }
    }

    void setNoSynButton_actionPerformed(ActionEvent e) {
        try {
            int theRow = theNameTable.getSelectedRow();
            theTS.removeSynonyms(theRow);
            TableReportNoPagingTableModel model = (TableReportNoPagingTableModel)theNameTable.getModel();
            model.fireTableChanged();
        } catch(Exception ex) {
            System.out.println("Make Sure a row is selected");
        }
    }

    void synTable_mouseClicked(MouseEvent e) {
        if (cbxSortAlphabetically.isSelected()) {
            sortListAlphabetically();
        } else {
            sortListByMatchScore();
        }
    }

    void okButton_actionPerformed(ActionEvent e) {
        theTS.saveTempTaxaList();
        isCanceled = false;
        this.setVisible(false);
    }

    void sortListByMatchScore() {
        Object theSynonym = theNameTable.getModel().getValueAt(theNameTable.getSelectedRow(), 0);
        ArrayList findOrderedMatches = theTS.findOrderedMatches((String) theSynonym, 4);
        DefaultListModel dlm = new DefaultListModel();
        Object[] a = findOrderedMatches.toArray();
        for (int i = 0; i < a.length; i++) {
            dlm.insertElementAt(a[i], i);
        }
        matchList.setModel(dlm);
    }

    private void cbxSortAlphabetically_actionPerformed(ActionEvent e) {

        if (cbxSortAlphabetically.isSelected()) {
            sortListAlphabetically();
        } else {
            sortListByMatchScore();
        }
    }

    private void sortListAlphabetically() {
        DefaultListModel listModel = (DefaultListModel) matchList.getModel();

        int itemCount = listModel.getSize();
        String[] a = new String[itemCount];

        listModel.copyInto(a);

        sortArray(Collator.getInstance(), a);

        for (int i = 0; i < itemCount; i++) {
            listModel.setElementAt(a[i], i);
        }
    }

    private void sortArray(Collator collator, String[] strArray) {
        String tmp;
        if (strArray.length == 1) {
            return;
        }
        for (int i = 0; i < strArray.length; i++) {
            for (int j = i + 1; j < strArray.length; j++) {
                if (collator.compare(strArray[i], strArray[j]) > 0) {
                    tmp = strArray[i];
                    strArray[i] = strArray[j];
                    strArray[j] = tmp;
                }
            }
        }
    }
}
class SynMenuDialog extends JDialog {
    private JFrame frmSynonymizerOperationMode;
    private boolean[] option;

    /**
     * Create the application.
     */
    public SynMenuDialog(boolean[] option,Frame frame) {
        super((Frame) frame, true);
        this.option = option;
        this.frmSynonymizerOperationMode = (JFrame)frame;
        initialize();
        this.pack();
    }

    /**
     * Initialize the contents of the frame.
     */
    private void initialize() {
        BorderLayout borderLayout = (BorderLayout) frmSynonymizerOperationMode.getContentPane().getLayout();
        borderLayout.setVgap(5);        
        
        JPanel basePanel = new JPanel();
        basePanel.setPreferredSize(new Dimension(500,300));
        this.getContentPane().add(basePanel, BorderLayout.CENTER);
        basePanel.setLayout(new GridLayout(3, 1, 0, 0));
        
        JButton btnChooseOpt1 = new JButton("Choose");
        btnChooseOpt1.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    menuButtonPressed(0);    
                    setVisible(false);//TODO Change to dispose with the rest of them
                }
        });
        basePanel.add(btnChooseOpt1);
        
        JLabel lblGenerateSynonymList = new JLabel("<html>1. Generate Synonym List From 2 Files Containing Taxa.</html>");
        basePanel.add(lblGenerateSynonymList);
        
        JButton btnChooseOpt2 = new JButton("Choose");
        btnChooseOpt2.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    menuButtonPressed(1);
                    setVisible(false);
                }
        });
        basePanel.add(btnChooseOpt2);
        
        JLabel lblManuallyEditSynonym = new JLabel("<html>2. Manually Edit Synonym List</html>");
        basePanel.add(lblManuallyEditSynonym);
        
        JButton btnChooseOpt3 = new JButton("Choose");
        btnChooseOpt3.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    menuButtonPressed(2);
                    setVisible(false);
                }
        });
        basePanel.add(btnChooseOpt3);
        
        JLabel lblApplySynonymList = new JLabel("<html>3. Apply a Synonym List to a Target File<html>");
        basePanel.add(lblApplySynonymList);
        
        JPanel cancelButtonPanel = new JPanel();
        this.getContentPane().add(cancelButtonPanel, BorderLayout.SOUTH);
        
        JButton btnCancel = new JButton("Cancel");
        btnCancel.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                        cancelButtonPressed();
                }
        });
        cancelButtonPanel.add(btnCancel);
        
        JLabel lblNewLabel = new JLabel("Choose Synonymizer Mode");
        lblNewLabel.setVerticalAlignment(SwingConstants.CENTER);
        lblNewLabel.setHorizontalAlignment(SwingConstants.CENTER);
        lblNewLabel.setEnabled(true);
        this.getContentPane().add(lblNewLabel, BorderLayout.NORTH);
        
        this.setTitle("Synonymizer Operation Mode");
        this.getContentPane().setSize(500, 300);
    }
    
    void cancelButtonPressed() {
            this.setVisible(false);
    }
    void menuButtonPressed(int index) {
        if(index == option.length) {
            for(int i = 0; i < option.length; i++) {
                option[i] = true;
            }
        }
        else {
            option[index] = true;
        }
    }

}

class SynonymizerFileChooser extends JDialog {
    
    JFrame theFrame;
    Object[] datum;
    int[] fileOptions;
    String[] initialSelections;
    String step;
    String delimiter;
    /**
     * Create the application.
     */
    public SynonymizerFileChooser(Frame theFrame,Object[] datum,int[] fileOptions,String[] initialSelections,String step) {
        super((Frame) theFrame, true);
        this.theFrame = (JFrame)theFrame;
        this.datum = datum;
        this.fileOptions = fileOptions;
        this.initialSelections = initialSelections;
        this.step = step;
        initialize();
        pack();
    }

    /**
     * Initialize the contents of the frame.
     */
    private void initialize() {
        if(step.equals("Step1")) {
            init_Step1();
        }
        else if(step.equals("Step2")) {
            init_Step2();
        }
        else if(step.equals("Step3")) {
            init_Step3();
        } 
    }

    private void init_Step1() {
        ArrayList<JComboBox> comboBoxes = new ArrayList<JComboBox>();
        
        JPanel panel = new JPanel();
        panel.setPreferredSize(new Dimension(500,300));
        this.getContentPane().add(panel, BorderLayout.CENTER);
        GridBagLayout gbl_panel = new GridBagLayout();
        gbl_panel.columnWidths = new int[]{0, 0, 0, 0};
        gbl_panel.rowHeights = new int[]{49, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        gbl_panel.columnWeights = new double[]{0.0, 1.0, 0.0, Double.MIN_VALUE};
        gbl_panel.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Double.MIN_VALUE};
        panel.setLayout(gbl_panel);
                
        JLabel lblSelectAReference = new JLabel("Select a File to Generate Synonyms");
   
        GridBagConstraints gbc_lblSelectAReference = getConstraints(1,1,new Insets(0, 10, 5, 10));
        gbc_lblSelectAReference.anchor = GridBagConstraints.WEST;
        panel.add(lblSelectAReference, gbc_lblSelectAReference);
        
        String[] listOfFileNames = getFileNamesFromDatum(datum);
        comboBoxes.add(getComboBoxWithSelection(listOfFileNames,initialSelections[0]));
        
        GridBagConstraints gbc_comboBoxReference = getConstraints(1,2,new Insets(0, 10, 5, 10));
        gbc_comboBoxReference.fill = GridBagConstraints.HORIZONTAL;
        panel.add(comboBoxes.get(0), gbc_comboBoxReference);
       
        JLabel lblSelectASynonymize = new JLabel("Select A File to Be Synonymized");
        GridBagConstraints gbc_lblSelectASynonymize = getConstraints(1,4,new Insets(0, 10, 5, 10));
        gbc_lblSelectASynonymize.anchor = GridBagConstraints.WEST;
        panel.add(lblSelectASynonymize, gbc_lblSelectASynonymize);
        
        comboBoxes.add(getComboBoxWithSelection(listOfFileNames,initialSelections[1]));
          
        GridBagConstraints gbc_comboBoxSynonym = getConstraints(1,5,new Insets(0, 10, 5, 10));
        gbc_comboBoxSynonym.fill = GridBagConstraints.HORIZONTAL;       
        panel.add(comboBoxes.get(1), gbc_comboBoxSynonym);
        
        JLabel lblSynonymizeTechnique = new JLabel("Similarity Metric");
        GridBagConstraints gbc_lblSynonymizeTechnique = getConstraints(1,7, new Insets(0, 10, 5, 10));
        gbc_lblSynonymizeTechnique.anchor = GridBagConstraints.WEST;
        panel.add(lblSynonymizeTechnique, gbc_lblSynonymizeTechnique);
        
        comboBoxes.add(getComboBoxWithSelection(new String[] {"Dice's Coefficient(Default Technique)",
                                                                "Edit Distance",
                                                                "Dynamic Time Warping using Hamming Distance",
                                                                "Dynamic Time Warping using Keyboard Distance",
                                                                "Hamming Distance using Soundex Encoding",
                                                                "Dice's Coefficient using Metaphone Encoding",
                                                                "Edit Distance using Metaphone Encoding",
                                                                "Delimiter Based Coefficient"
                                                                },"Dice's Coefficient(Default Technique)"));

        GridBagConstraints gbc_comboBoxTechniques = getConstraints(1,8,new Insets(0, 10, 0, 10));
        gbc_comboBoxTechniques.fill = GridBagConstraints.HORIZONTAL;
       
        
        
        panel.add(comboBoxes.get(2), gbc_comboBoxTechniques);
        
        
        JLabel lblSynonymizeDelim = new JLabel("Delimiter(Not Required)");
        lblSynonymizeDelim.setToolTipText("For use with concatenated Taxa names");
        GridBagConstraints gbc_lblSynonymizeDelim = getConstraints(1,10, new Insets(0, 10, 0, 10));
        gbc_lblSynonymizeDelim.anchor = GridBagConstraints.WEST;
        panel.add(lblSynonymizeDelim, gbc_lblSynonymizeDelim);
        
        JTextField fldSynonymizeDelim = new JTextField();
        fldSynonymizeDelim.setPreferredSize(new Dimension(60,30));
        GridBagConstraints gbc_fldSynonymizeDelim = getConstraints(1,11,new Insets(0,10,20,10));
        gbc_fldSynonymizeDelim.anchor = GridBagConstraints.WEST;
        gbc_fldSynonymizeDelim.fill = GridBagConstraints.HORIZONTAL;
        panel.add(fldSynonymizeDelim, gbc_fldSynonymizeDelim);
        
        lblSynonymizeDelim.setVisible(false);
        fldSynonymizeDelim.setVisible(false);
        
        comboBoxes.get(2).addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if(comboBoxes.get(2).getSelectedIndex()==7) {
                    lblSynonymizeDelim.setVisible(true);
                    fldSynonymizeDelim.setVisible(true);
                }
                else {
                    lblSynonymizeDelim.setVisible(false);
                    fldSynonymizeDelim.setVisible(false);
                }
            }
        });
        JLabel lblSelectFilesTo = new JLabel("Select Files to Generate Synonym List");
        lblSelectFilesTo.setHorizontalAlignment(SwingConstants.CENTER);
        this.getContentPane().add(lblSelectFilesTo, BorderLayout.NORTH);
        
        JPanel bottomPanel = new JPanel();
        this.getContentPane().add(bottomPanel, BorderLayout.SOUTH);
        GridBagLayout gbl_bottomPanel = new GridBagLayout();
        gbl_bottomPanel.columnWidths = new int[]{0, 0, 0, 43, 19, 0, 0};
        gbl_bottomPanel.rowHeights = new int[] {15, 0, 15, 0};
        gbl_bottomPanel.columnWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Double.MIN_VALUE};
        gbl_bottomPanel.rowWeights = new double[]{0.0, 0.0, 0.0, Double.MIN_VALUE};
        bottomPanel.setLayout(gbl_bottomPanel);
        
        JButton btnCancel = new JButton("Cancel");
        btnCancel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButtonPressed(comboBoxes);
            }
        });
   
        GridBagConstraints gbc_btnCancel = getConstraints(4,1,new Insets(0,0,5,5));
        bottomPanel.add(btnCancel, gbc_btnCancel);
        
        JButton btnSubmit = new JButton("Submit");
        btnSubmit.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                submitButtonPressed(comboBoxes,fldSynonymizeDelim.getText());
            }
        });
      
        GridBagConstraints gbc_btnSubmit = getConstraints(5,1,new Insets(0,0,5,0));        
        bottomPanel.add(btnSubmit, gbc_btnSubmit);
        
        this.getContentPane().setSize(500,300);
        this.setTitle("Choose Files and the Synonymize Technique");
    }
    
    private void init_Step2() {
        ArrayList<JComboBox> comboBoxes = new ArrayList<JComboBox>();
        ArrayList<Integer> comboBoxLabel = new ArrayList<Integer>();
        
        JPanel panel = new JPanel();
        panel.setPreferredSize(new Dimension(500,300));
        this.getContentPane().add(panel, BorderLayout.CENTER);
        GridBagLayout gbl_panel = new GridBagLayout();
        gbl_panel.columnWidths = new int[]{0, 0, 0, 0};
        gbl_panel.rowHeights = new int[]{49, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        gbl_panel.columnWeights = new double[]{0.0, 1.0, 0.0, Double.MIN_VALUE};
        gbl_panel.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Double.MIN_VALUE};
        panel.setLayout(gbl_panel);
        
        JLabel lblSelectAReference = new JLabel("Select a File to Manually Edited:");

        GridBagConstraints gbc_lblSelectAReference = getConstraints(1,1,new Insets(0, 10, 5, 10));
        gbc_lblSelectAReference.anchor = GridBagConstraints.WEST;
        panel.add(lblSelectAReference, gbc_lblSelectAReference);
        
        String[] listOfFileNames = getFileNamesFromDatumSyn(datum);
        
        comboBoxes.add(getComboBoxWithSelection(listOfFileNames,initialSelections[0]));

        GridBagConstraints gbc_comboBoxReference = getConstraints(1, 2, new Insets(0, 10, 5, 10));
        gbc_comboBoxReference.fill = GridBagConstraints.HORIZONTAL;
        panel.add(comboBoxes.get(0), gbc_comboBoxReference);
        comboBoxLabel.add(1);
        
        JPanel bottomPanel = new JPanel();
        this.getContentPane().add(bottomPanel, BorderLayout.SOUTH);
        GridBagLayout gbl_bottomPanel = new GridBagLayout();
        gbl_bottomPanel.columnWidths = new int[]{0, 0, 0, 43, 19, 0, 0};
        gbl_bottomPanel.rowHeights = new int[] {15, 0, 15, 0};
        gbl_bottomPanel.columnWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Double.MIN_VALUE};
        gbl_bottomPanel.rowWeights = new double[]{0.0, 0.0, 0.0, Double.MIN_VALUE};
        bottomPanel.setLayout(gbl_bottomPanel);
        
        JButton btnCancel = new JButton("Cancel");
        btnCancel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButtonPressed(comboBoxes);
            }
        });

        GridBagConstraints gbc_btnCancel = getConstraints(4,1, new Insets(0, 0, 5, 5));
        bottomPanel.add(btnCancel, gbc_btnCancel);
        
        JButton btnSubmit = new JButton("Submit");
        btnSubmit.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                submitButtonPressed(comboBoxes,getSyn_NonSynMapping(datum),comboBoxLabel);
            }
        });

        GridBagConstraints gbc_btnSubmit = getConstraints(5, 1, new Insets(0, 0, 5, 5));
        bottomPanel.add(btnSubmit, gbc_btnSubmit);
        
        this.getContentPane().setSize(500,300);
        this.setTitle("Choose File for Manual Update");
    }
    
    private void init_Step3() {
        ArrayList<JComboBox> comboBoxes = new ArrayList<JComboBox>();
        ArrayList<Integer> comboBoxLabel = new ArrayList<Integer>();
        
        JPanel panel = new JPanel();
        panel.setPreferredSize(new Dimension(500,300));
        this.getContentPane().add(panel, BorderLayout.CENTER);
        GridBagLayout gbl_panel = new GridBagLayout();
        gbl_panel.columnWidths = new int[]{0, 0, 0, 0};
        gbl_panel.rowHeights = new int[]{49, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        gbl_panel.columnWeights = new double[]{0.0, 1.0, 0.0, Double.MIN_VALUE};
        gbl_panel.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Double.MIN_VALUE};
        panel.setLayout(gbl_panel);
        
        JLabel lblSelectAReference = new JLabel("Select a Synonym File");

        GridBagConstraints gbc_lblSelectAReference = getConstraints(1, 1, new Insets(0, 10, 5, 10));
        gbc_lblSelectAReference.anchor = GridBagConstraints.WEST;
        
        panel.add(lblSelectAReference, gbc_lblSelectAReference);
        
        String[] listOfSynFileNames = getFileNamesFromDatumSyn(datum);
        
        comboBoxes.add(getComboBoxWithSelection(listOfSynFileNames,initialSelections[0]));
             
        GridBagConstraints gbc_comboBoxReference = getConstraints(1, 2, new Insets(0, 10, 5, 10));
        gbc_comboBoxReference.fill = GridBagConstraints.HORIZONTAL;
        
        panel.add(comboBoxes.get(0),gbc_comboBoxReference);
        
        comboBoxLabel.add(1);
        
        JLabel lblSelectASynonymize = new JLabel("Select A File to Be Updated with Synonyms");
 
        GridBagConstraints gbc_lblSelectASynonymize = getConstraints(1, 4, new Insets(0, 10, 5, 10));
        gbc_lblSelectASynonymize.anchor = GridBagConstraints.WEST;
        panel.add(lblSelectASynonymize, gbc_lblSelectASynonymize);
        
        String[] listOfFileNames = getFileNamesFromDatumNonSyn(datum);

        comboBoxes.add(getComboBoxWithSelection(listOfFileNames,initialSelections[1]));
        
        GridBagConstraints gbc_comboBoxSynonym = getConstraints(1, 5, new Insets(0, 10, 5, 10));
        gbc_comboBoxSynonym.fill = GridBagConstraints.HORIZONTAL;
        panel.add(comboBoxes.get(1),gbc_comboBoxSynonym);
        comboBoxLabel.add(0);
        
        JLabel lblSelectFilesTo = new JLabel("Select Files to Generate Synonym List");
        lblSelectFilesTo.setHorizontalAlignment(SwingConstants.CENTER);
        this.getContentPane().add(lblSelectFilesTo, BorderLayout.NORTH);
        
        JPanel bottomPanel = new JPanel();
        this.getContentPane().add(bottomPanel, BorderLayout.SOUTH);
        GridBagLayout gbl_bottomPanel = new GridBagLayout();
        gbl_bottomPanel.columnWidths = new int[]{0, 0, 0, 43, 19, 0, 0};
        gbl_bottomPanel.rowHeights = new int[] {15, 0, 15, 0};
        gbl_bottomPanel.columnWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Double.MIN_VALUE};
        gbl_bottomPanel.rowWeights = new double[]{0.0, 0.0, 0.0, Double.MIN_VALUE};
        bottomPanel.setLayout(gbl_bottomPanel);
        
        JButton btnCancel = new JButton("Cancel");
        btnCancel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButtonPressed(comboBoxes);
            }
        });

        GridBagConstraints gbc_btnCancel = getConstraints(4, 1, new Insets(0, 0, 5, 5));
        bottomPanel.add(btnCancel, gbc_btnCancel);
        
        JButton btnSubmit = new JButton("Submit");
        btnSubmit.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                submitButtonPressed(comboBoxes,getSyn_NonSynMapping(datum),comboBoxLabel);
            }
        });

        GridBagConstraints gbc_btnSubmit = getConstraints(5, 1, new Insets(0, 0, 5, 5));
        bottomPanel.add(btnSubmit, gbc_btnSubmit);
        
        this.getContentPane().setSize(500,300);
        this.setTitle("Choose Files");
    }
    
    private JComboBox getComboBoxWithSelection(String[] model, String selection) {
        JComboBox comboBox = new JComboBox();
        comboBox.setModel(new DefaultComboBoxModel(model));
        for(int i = 0; i<model.length; i++) {
            if(model[i].equals(selection)) {
                comboBox.setSelectedIndex(i);
            }
        }
        return comboBox;
    }
    
    private GridBagConstraints getConstraints(int gridx, int gridy, Insets inset) {
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = inset;
        gbc.gridx = gridx;
        gbc.gridy = gridy;
        return gbc;
   
    }
    private String[] getFileNamesFromDatum(Object[] datum) {
        String[] names = new String[datum.length];
        for(int i = 0; i<datum.length;i++) {
            Datum singleRecord = (Datum)datum[i];
            names[i] = singleRecord.getName();
        }
        return names;
    }
    private String[] getFileNamesFromDatumSyn(Object[] datum) {
        ArrayList<String> names = new ArrayList<String>();
        
        for(Object record:datum) {
            Datum singleRecord = (Datum)record;
            if(singleRecord.getDataType().equals(IdentifierSynonymizer.class)){
                names.add(singleRecord.getName());
            }
        }
        
        String[] namesArray = new String[names.size()];
        for(int i = 0; i<namesArray.length;i++) {
            namesArray[i] = names.get(i);
        }
        return namesArray;
    }
    private String[] getFileNamesFromDatumSyn(Object[] datum, ArrayList<Integer> indexArray) {
        ArrayList<String> names = new ArrayList<String>();
        
        for(int i = 0; i<datum.length;i++) {
            Datum singleRecord = (Datum)datum[i];
            if(singleRecord.getDataType().equals(IdentifierSynonymizer.class)){
                indexArray.add(i);
                names.add(singleRecord.getName());
            }
        }
        
        String[] namesArray = new String[names.size()];
        for(int i = 0; i<namesArray.length;i++) {
            namesArray[i] = names.get(i);
        }
        return namesArray;
    }
    private ArrayList<ArrayList<Integer>> getSyn_NonSynMapping(Object[] datum) {
        ArrayList<Integer> synList = new ArrayList<Integer>();
        ArrayList<Integer> nonSynList = new ArrayList<Integer>();
        
        for(int i = 0; i<datum.length; i++) {
            Datum singleRecord = (Datum)datum[i];
            if(singleRecord.getDataType().equals(IdentifierSynonymizer.class)){
                synList.add(i);
            }
            else {
                nonSynList.add(i);
            }
        }
        ArrayList<ArrayList<Integer>> combinedList = new ArrayList<ArrayList<Integer>>();
        combinedList.add(nonSynList);
        combinedList.add(synList);
        
        return combinedList;
    }
    private String[] getFileNamesFromDatumNonSyn(Object[] datum) {
        ArrayList<String> names = new ArrayList<String>();
        
        for(Object record:datum) {
            Datum singleRecord = (Datum)record;
            if(!singleRecord.getDataType().equals(IdentifierSynonymizer.class)){
                names.add(singleRecord.getName());
            }
        }
        
        String[] namesArray = new String[names.size()];
        for(int i = 0; i<namesArray.length;i++) {
            namesArray[i] = names.get(i);
        }
        return namesArray;
    }
    void cancelButtonPressed(ArrayList<JComboBox> comboBoxes) {
        for(int i = 0; i<fileOptions.length;i++) {
            fileOptions[i] = -1;
        }
        setVisible(false);
    }
    void submitButtonPressed(ArrayList<JComboBox> comboBoxes) {
        for(int i = 0; i<fileOptions.length;i++) {
            fileOptions[i] = comboBoxes.get(i).getSelectedIndex();
        }
        setVisible(false);
    }
    void submitButtonPressed(ArrayList<JComboBox> comboBoxes,String delim) {
        for(int i = 0; i<fileOptions.length;i++) {
            fileOptions[i] = comboBoxes.get(i).getSelectedIndex();
        }
        this.delimiter = delim;
        setVisible(false);
    }
    void submitButtonPressed(ArrayList<JComboBox> comboBoxes, ArrayList<ArrayList<Integer>> synListMap, ArrayList<Integer> comboBoxLabel) {
        for(int i = 0; i<fileOptions.length; i++) {
            fileOptions[i] = synListMap.get(comboBoxLabel.get(i))
                                       .get(comboBoxes.get(i).getSelectedIndex());
        }
        setVisible(false);
    }
    public String getDelim() {
        return delimiter;
    }
}
