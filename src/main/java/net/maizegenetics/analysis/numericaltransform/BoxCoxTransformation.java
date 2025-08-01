package net.maizegenetics.analysis.numericaltransform;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.stream.Collectors;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeUtils;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.taxa.Taxon;

public class BoxCoxTransformation {

    public Phenotype applyBoxCox(Phenotype pheno, boolean addSmallVal, long randomSeed,double startLambda, double endLambda, double stepLambda) throws Exception{
      //Get the variable names and put it into a list
        ArrayList<PhenotypeAttribute> attributes = (ArrayList<PhenotypeAttribute>)pheno.attributeListOfType(ATTRIBUTE_TYPE.data);
        
        //get the covariates
        attributes.addAll((ArrayList<PhenotypeAttribute>)pheno.attributeListOfType(ATTRIBUTE_TYPE.covariate));
        
        ArrayList<String> phenoNames = (ArrayList<String>)attributes.stream().map((phenoAttr) -> phenoAttr.name()).collect(Collectors.toList());
        
        //Loop through the List<Taxa> get the unique names
        ArrayList<Taxon> taxonList = (ArrayList<Taxon>)pheno.taxaAttribute().allTaxaAsList();
        ArrayList<String> uniqueTaxaNames =(ArrayList<String>) pheno.taxaAttribute().allTaxaAsList().stream()
                                                                .map((taxon) -> taxon.getName()).collect(Collectors.toList());
        
        //loop through the all the values and find the minimum value
        double minValue = Double.MAX_VALUE;
        
        HashMap<String,ArrayList<Double>> phenoToDataList = new HashMap<String,ArrayList<Double>>();
        
        for(String varName : phenoNames) {
            ArrayList<Double> currentPhenoValues = new ArrayList<Double>();
            for(int i = 0; i < taxonList.size(); i++) {
                double currentValue = ((Float)pheno.value(i,pheno.attributeIndexForName(varName))).doubleValue();
                currentPhenoValues.add(currentValue);
                if( currentValue < minValue && currentValue > 0) {
                    minValue = currentValue;
                }
            }
            phenoToDataList.put(varName,currentPhenoValues);
        }
        
        if(!addSmallVal) {
            minValue = 0.0;
        }
        else {
            System.out.println(minValue);
        }
        
       //Loop through the Phenotype Attributes and Box Cox them
        double[][] fullTransformedVals = new double[uniqueTaxaNames.size()][phenoNames.size()];

        for(int phenoNameIndex = 0; phenoNameIndex < phenoNames.size(); phenoNameIndex++){
            double[] boxCoxValues = computeBoxCox(phenoToDataList.get(phenoNames.get(phenoNameIndex)),minValue,randomSeed,startLambda, endLambda, stepLambda);
            for(int taxaIndex = 0; taxaIndex < uniqueTaxaNames.size(); taxaIndex++) {
                fullTransformedVals[taxaIndex][phenoNameIndex] = boxCoxValues[taxaIndex];
            }
        }
        
        //Create 2-D ArrayList to store all of the transformed values
        ArrayList<ArrayList<Double>> avgValues = new ArrayList<ArrayList<Double>>();
        for(int taxaCounter = 0; taxaCounter < fullTransformedVals.length; taxaCounter++) {
            ArrayList<Double> currentTaxaVals = new ArrayList<Double>();
            for(int phenoCounter = 0; phenoCounter < fullTransformedVals[taxaCounter].length; phenoCounter++) {
                currentTaxaVals.add(fullTransformedVals[taxaCounter][phenoCounter]);
            }
            avgValues.add(currentTaxaVals);
        }
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        for(int i = 0; i < pheno.attributeListCopy().size(); i++) {
            originalTypes.add(pheno.attributeType(i));
        }
        
        try {
            //Take the taxa List, the variable name list and the matrix and send to PhenoUtils package to get a phenotype back
            return PhenotypeUtils.createPhenotypeFromTransform(uniqueTaxaNames, phenoNames, avgValues,pheno.name()+"_BoxCoxTransformed",originalTypes);
        }
        catch(Exception e) {
            throw e;
        }
        
    }
    
    public static double[] computeBoxCox(ArrayList<Double> phenotypeValues,double minimumValue,long randomSeed, double startLambda, double endLambda, double stepLambda) {
        //add in the small value if needed
        //Select a random small value where 0<randVal<.5*minValue
        //Doing this by finding a random double between 0 and 1 then multiplying against .5*minValue
        Random rand = new Random(randomSeed);
        ArrayList<Double> addedSmallPhenoValues = (ArrayList<Double>) phenotypeValues.stream().map((currentVal) -> {
            return rand.nextDouble() * 0.5 * minimumValue + currentVal;
        }).collect(Collectors.toList());
       
        //Check attribute to make sure we have variation
        if(!attributeHasVariation(addedSmallPhenoValues)) {
            //if we dont we need to return the original value as box cox will not do anything
            double[] valuesToReturn = new double[addedSmallPhenoValues.size()];
            for(int i = 0; i < addedSmallPhenoValues.size(); i++) {
                valuesToReturn[i] = addedSmallPhenoValues.get(i);
            }
            return valuesToReturn;
        }
        else {
            //Loop through lambdas and compute the KST
            double currentLambda = startLambda;
            double[] currentValues = new double[addedSmallPhenoValues.size()];
            double bestTestStat = Double.MAX_VALUE;
            
            double bestMean = 0.0;
            double bestStDev = 0.0;
            
            //Loop through lambdas -5 to 5 increasing by .2 each time
            for(double lambda = startLambda; lambda < endLambda; lambda += stepLambda) {
                //create a new array to hold the values
                double[] transformedValues = new double[addedSmallPhenoValues.size()];
                
                //get the transformed values
                for(int i = 0; i < addedSmallPhenoValues.size(); i++) {
                    transformedValues[i] = boxCoxTransform(addedSmallPhenoValues.get(i), lambda);
                }
                
                //Check to see how well it fits a normal distribution
                //Use the following method
                //
                //Get Mean and standard deviation of the transformed data
                //Create a Normal Distribution using these input params
                //Use KolmogorovSmirnovTest to see how well the data fits the normal distribution
                //Pick the one with the lowest values
                StandardDeviation sdev = new StandardDeviation();
                double sampleStandardDev = sdev.evaluate(transformedValues);

                //calculate the mean
                Mean meanObj = new Mean();
                double meanVal = meanObj.evaluate(transformedValues);
                
                //Skip this lambda if we have a negative sample standard deviation
                if(sampleStandardDev<=0.0 || Double.isNaN(sampleStandardDev) || Double.isNaN(meanVal)) {
                    continue;
                }
                //Compute kst
                NormalDistribution normDist = new NormalDistribution(meanVal, sampleStandardDev);
                KolmogorovSmirnovTest kst = new KolmogorovSmirnovTest();
                
                //Get the statistic
                double testStat = kst.kolmogorovSmirnovStatistic(normDist, transformedValues);
                
                //If it is less than our current best, we update the statistics
                if(testStat < bestTestStat) {
                    currentLambda = lambda;
                    currentValues = transformedValues;
                    bestTestStat = testStat;
                    bestMean = meanVal;
                    bestStDev = sampleStandardDev;
                }
            }
            return currentValues;
        }
    }
    
    
    private static double boxCoxTransform(double value, double lambda) {
        if(lambda == 0) {
            //do the natural log
            return Math.log(value);
        }
        else {
            //Compute (value^lamdba - 1)/lambda
            return (Math.pow(value, lambda) - 1.0)/lambda;
        }
    }
    
    private static boolean attributeHasVariation(ArrayList<Double> phenoAttribute) {
        double prevVal = phenoAttribute.get(0);
        
        for(int i = 1; i < phenoAttribute.size(); i++) { 
            if(phenoAttribute.get(i) != prevVal) {
                return true;
            }
        }
        return false;
    }
}
