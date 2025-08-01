package net.maizegenetics.analysis.numericaltransform;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

import org.apache.commons.math3.stat.descriptive.moment.Mean;

import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeUtils;
import net.maizegenetics.taxa.Taxon;

public class AvgPhenotype {

    public Phenotype averagePheno(Phenotype pheno, boolean addSmallValue, long randomSeed) throws Exception {
      
        
        //Get the variable names and put it into a list
        ArrayList<PhenotypeAttribute> attributes = (ArrayList<PhenotypeAttribute>)pheno.attributeListOfType(ATTRIBUTE_TYPE.data);

        //get the covariates
        attributes.addAll((ArrayList<PhenotypeAttribute>)pheno.attributeListOfType(ATTRIBUTE_TYPE.covariate));
        
        ArrayList<String> phenoNames = (ArrayList<String>)attributes.stream().map((phenoAttr) -> phenoAttr.name()).collect(Collectors.toList());
        
        
        //Loop through the List<Taxa> get the unique names
        ArrayList<Taxon> taxonList = (ArrayList<Taxon>)pheno.taxaAttribute().allTaxaAsList();
        ArrayList<String> uniqueTaxaNames =(ArrayList<String>) pheno.taxaAttribute().allTaxaAsList().stream()
                                                                .map((taxon) -> taxon.getName()).distinct().collect(Collectors.toList());
        
        HashMap<String, ArrayList<ArrayList<Double>>> rowMapByTaxa= new HashMap<String,ArrayList<ArrayList<Double>>>();
        
        //Create an ArrayList object for each unique taxaName
        for(String taxaName : uniqueTaxaNames) {
            ArrayList<ArrayList<Double>> list = new ArrayList<ArrayList<Double>>();
            rowMapByTaxa.put(taxaName, list);
        }
        double minValue = Double.MAX_VALUE;
        
        //Go through all taxon and use the name to associate duplicates together
        for(int i = 0; i < taxonList.size(); i++) {
            //Make a Container to hold the values
            ArrayList<Double> currentRow = new ArrayList<Double>();
            for(String varName : phenoNames) {
                //Add this to the duplicate list
                currentRow.add(((Float)pheno.value(i,pheno.attributeIndexForName(varName))).doubleValue());
                
                if(((Float)pheno.value(i,pheno.attributeIndexForName(varName))).doubleValue() < minValue && 
                        ((Float)pheno.value(i,pheno.attributeIndexForName(varName))).doubleValue() > 0) {
                    minValue = ((Float)pheno.value(i,pheno.attributeIndexForName(varName))).doubleValue();
                }
            }
            //Add the duplicate row(by taxa name) to the mapping
            rowMapByTaxa.get(taxonList.get(i)).add(currentRow);
        }
        
        //Select a random small value where 0<randVal<.5*minValue
        //Doing this by finding a random double between 0 and 1 then multiplying against .5*minValue
        Random rand = new Random(randomSeed);
        
        //Create 2-D ArrayList to store all of the averaged values
        ArrayList<ArrayList<Double>> avgValues = new ArrayList<ArrayList<Double>>();
        //for each taxa sum up all values and divide
        //Use Commons Math in case we want to do a geometric mean
        for(String taxaName : uniqueTaxaNames) {
            //Get the list of all the duplicate values for the taxaName
            ArrayList<ArrayList<Double>> currentValues = rowMapByTaxa.get(taxaName);
            
            //Make the container to hold the current average values
            ArrayList<Double> currentAverages = new ArrayList<Double>();
            //Loop through this 2-D array to and apply the average
            for(int col = 0; col < currentValues.get(0).size(); col++) {
                double[] values = new double[currentValues.size()];
                Mean meanObj = new Mean();
                for(int row = 0; row < currentValues.size(); row++) {
                   values[row] = currentValues.get(row).get(col);
                   //if we need to add in small values do it here.  will shift the average by the small val amount
                   if(addSmallValue) {
                       values[row] = values[row] + rand.nextDouble() * 0.5 * minValue; 
                   }
                }
                //Add the averaged value to the finished current row
                currentAverages.add(meanObj.evaluate(values));
            }
            //Add the averaged current row to the avgValues 2-D array
            avgValues.add(currentAverages);
        }
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        for(int i = 0; i < pheno.attributeListCopy().size(); i++) {
            originalTypes.add(pheno.attributeType(i));
        }
        try {
            //Take the taxa List, the variable name list and the matrix and send to PhenoUtils package to get a phenotype back
            return PhenotypeUtils.createPhenotypeFromTransform(uniqueTaxaNames, phenoNames, avgValues,pheno.name()+"_AveragedByTaxa",originalTypes);
        }
        catch(Exception e) {
            throw e;
        }
    }
}
