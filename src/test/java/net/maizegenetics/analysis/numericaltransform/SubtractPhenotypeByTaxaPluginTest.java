package net.maizegenetics.analysis.numericaltransform;

import java.util.ArrayList;
import java.util.Arrays;

import junit.framework.Assert;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeUtils;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;

import org.junit.Test;

public class SubtractPhenotypeByTaxaPluginTest {

    //With Same Taxa and Variables
    @Test
    public void testSubtractPhenoTest() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNames2 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNames2 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> data2 = new ArrayList<ArrayList<Double>>();
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,100.0,100.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,300.0,300.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> dataKnownSub = new ArrayList<ArrayList<Double>>();
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{200.0,100.0,0.0})));
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{-200.0,-100.0,0.0})));
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{0.0,0.0,0.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"DataSet1",originalTypes);
            Phenotype pheno2 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames2, varNames2, data2,"DataSet2",originalTypes);
            Phenotype phenoKnownSub = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownSub, varNamesKnownSub, dataKnownSub,"KnownSubtraction",originalTypes);
            
            SubtractPhenotype subPhenotypeObject = new SubtractPhenotype();
            Phenotype phenoSub = subPhenotypeObject.subtractPhenotype(pheno1, pheno2, false);
            
            boolean isSame = true;
            for(int row = 0; row < phenoSub.getRowCount(); row++){
                for(int col = 0; col < phenoSub.getColumnCount(); col++) {
                    if(!phenoKnownSub.getValueAt(row, col).equals(phenoSub.getValueAt(row, col))) {
                        System.out.println(phenoKnownSub.getValueAt(row, col) + " "+phenoSub.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Subtract, Same Taxa, Same Variables Test Fails",isSame);
        }
        catch(Exception e) {
            Assert.assertTrue(false);
            e.printStackTrace();
        }
    }
    
  //With Different Taxa and same Variables
    @Test
    public void testSubtractPhenoTestDifferentTaxa() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNames2 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML248","B37"}));
        ArrayList<String> taxaNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"B73","B37"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNames2 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> data2 = new ArrayList<ArrayList<Double>>();
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,100.0,100.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,300.0,300.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> dataKnownSub = new ArrayList<ArrayList<Double>>();
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{200.0,100.0,0.0})));
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{0.0,0.0,0.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"Dataset1",originalTypes);
            Phenotype pheno2 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames2, varNames2, data2,"Dataset2",originalTypes);
            Phenotype phenoKnownSub = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownSub, varNamesKnownSub, dataKnownSub,"DatasetSubtracted",originalTypes);
            
            SubtractPhenotype subPhenotypeObject = new SubtractPhenotype();
            Phenotype phenoSub = subPhenotypeObject.subtractPhenotype(pheno1, pheno2, false);
            
            boolean isSame = true;
            for(int row = 0; row < phenoSub.getRowCount(); row++){
                for(int col = 0; col < phenoSub.getColumnCount(); col++) {
                    if(!phenoKnownSub.getValueAt(row, col).equals(phenoSub.getValueAt(row, col))) {
                        System.out.println(phenoKnownSub.getValueAt(row, col) + " "+phenoSub.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Subtract, Different Taxa, Same Variables Test Fails",isSame);
        }
        catch(Exception e) {
            Assert.assertTrue(false);
            e.printStackTrace();
        }
    }
    
  //With same Taxa and different Variables
    @Test
    public void testSubtractPhenoTestDifferentVariable() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNames2 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNames2 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var4", "var3"}));
        ArrayList<String> varNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> data2 = new ArrayList<ArrayList<Double>>();
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,100.0,100.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,300.0,300.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> dataKnownSub = new ArrayList<ArrayList<Double>>();
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{200.0,0.0})));
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{-200.0,0.0})));
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{0.0,0.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        
        ArrayList<ATTRIBUTE_TYPE> newTypes = new ArrayList<>();
        newTypes.add(ATTRIBUTE_TYPE.taxa);
        newTypes.add(ATTRIBUTE_TYPE.data);
        newTypes.add(ATTRIBUTE_TYPE.data);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"Dataset1",originalTypes);
            Phenotype pheno2 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames2, varNames2, data2,"Dataset2",originalTypes);
            Phenotype phenoKnownSub = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownSub, varNamesKnownSub, dataKnownSub,"Dataset3",newTypes);
            
            SubtractPhenotype subPhenotypeObject = new SubtractPhenotype();
            Phenotype phenoSub = subPhenotypeObject.subtractPhenotype(pheno1, pheno2, false);
            
            boolean isSame = true;
            for(int row = 0; row < phenoSub.getRowCount(); row++){
                for(int col = 0; col < phenoSub.getColumnCount(); col++) {
                    if(!phenoKnownSub.getValueAt(row, col).equals(phenoSub.getValueAt(row, col))) {
                        System.out.println(phenoKnownSub.getValueAt(row, col) + " "+phenoSub.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Subtract, Same Taxa, Different Variables Test Fails",isSame);
        }
        catch(Exception e) {
            e.printStackTrace();
            Assert.assertTrue(false);
            
        }
    }
    //With different Taxa and different Variables
    @Test
    public void testSubtractPhenoTestDifferentTaxaAndVariable() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNames2 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML248","B37"}));
        ArrayList<String> taxaNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"B73","B37"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNames2 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var4", "var3"}));
        ArrayList<String> varNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> data2 = new ArrayList<ArrayList<Double>>();
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,100.0,100.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,300.0,300.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> dataKnownSub = new ArrayList<ArrayList<Double>>();
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{200.0,0.0})));
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{0.0,0.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        
        ArrayList<ATTRIBUTE_TYPE> newTypes = new ArrayList<>();
        newTypes.add(ATTRIBUTE_TYPE.taxa);
        newTypes.add(ATTRIBUTE_TYPE.data);
        newTypes.add(ATTRIBUTE_TYPE.data);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"Dataset1",originalTypes);
            Phenotype pheno2 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames2, varNames2, data2,"Dataset2",originalTypes);
            Phenotype phenoKnownSub = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownSub, varNamesKnownSub, dataKnownSub,"DatasetSubtracted",newTypes);
            
            SubtractPhenotype subPhenotypeObject = new SubtractPhenotype();
            Phenotype phenoSub = subPhenotypeObject.subtractPhenotype(pheno1, pheno2, false);
            
            boolean isSame = true;
            for(int row = 0; row < phenoSub.getRowCount(); row++){
                for(int col = 0; col < phenoSub.getColumnCount(); col++) {
                    if(!phenoKnownSub.getValueAt(row, col).equals(phenoSub.getValueAt(row, col))) {
                        System.out.println(phenoKnownSub.getValueAt(row, col) + " "+phenoSub.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Subtract, Different Taxa, Different Variables Test Fails",isSame);
        }
        catch(Exception e) {
            Assert.assertTrue(false);
            e.printStackTrace();
        }
    }
    
  //With Same Taxa and Variables
    @Test
    public void testSubtractPhenoTestWithCov() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNames2 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNames2 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> data2 = new ArrayList<ArrayList<Double>>();
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,100.0,100.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,300.0,300.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> dataKnownSub = new ArrayList<ArrayList<Double>>();
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{200.0,100.0,0.0})));
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{-200.0,-100.0,0.0})));
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{0.0,0.0,0.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.covariate);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"DataSet1",originalTypes);
            Phenotype pheno2 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames2, varNames2, data2,"DataSet2",originalTypes);
            Phenotype phenoKnownSub = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownSub, varNamesKnownSub, dataKnownSub,"KnownSubtraction",originalTypes);
            
            SubtractPhenotype subPhenotypeObject = new SubtractPhenotype();
            Phenotype phenoSub = subPhenotypeObject.subtractPhenotype(pheno1, pheno2, false);
            
            boolean isSame = true;
            for(int row = 0; row < phenoSub.getRowCount(); row++){
                for(int col = 0; col < phenoSub.getColumnCount(); col++) {
                    if(!phenoKnownSub.getValueAt(row, col).equals(phenoSub.getValueAt(row, col))) {
                        System.out.println(phenoKnownSub.getValueAt(row, col) + " "+phenoSub.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Subtract, Same Taxa, Same Variables With Covariates Test Fails",isSame);
            
            boolean attributesSame = true;
            for(int attrCounter = 0; attrCounter < originalTypes.size(); attrCounter++) {
                if(!originalTypes.get(attrCounter).equals(phenoSub.attributeType(attrCounter))) {
                    attributesSame = false;
                }
            }
            Assert.assertTrue("Phenotype Subtract, Same Taxa, Same Variables With Covariates Test Fails. Attributes are not matching",attributesSame);
            
        }
        catch(Exception e) {
            Assert.assertTrue(false);
            e.printStackTrace();
        }
    }
    
  //With Different Taxa and same Variables
    @Test
    public void testSubtractPhenoTestDifferentTaxaWithCov() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNames2 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML248","B37"}));
        ArrayList<String> taxaNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"B73","B37"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNames2 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> data2 = new ArrayList<ArrayList<Double>>();
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,100.0,100.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,300.0,300.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> dataKnownSub = new ArrayList<ArrayList<Double>>();
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{200.0,100.0,0.0})));
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{0.0,0.0,0.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.covariate);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"Dataset1",originalTypes);
            Phenotype pheno2 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames2, varNames2, data2,"Dataset2",originalTypes);
            Phenotype phenoKnownSub = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownSub, varNamesKnownSub, dataKnownSub,"DatasetSubtracted",originalTypes);
            
            SubtractPhenotype subPhenotypeObject = new SubtractPhenotype();
            Phenotype phenoSub = subPhenotypeObject.subtractPhenotype(pheno1, pheno2, false);
            
            boolean isSame = true;
            for(int row = 0; row < phenoSub.getRowCount(); row++){
                for(int col = 0; col < phenoSub.getColumnCount(); col++) {
                    if(!phenoKnownSub.getValueAt(row, col).equals(phenoSub.getValueAt(row, col))) {
                        System.out.println(phenoKnownSub.getValueAt(row, col) + " "+phenoSub.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Subtract, Different Taxa, Same Variables With Covariates Test Fails",isSame);
            
            boolean attributesSame = true;
            for(int attrCounter = 0; attrCounter < originalTypes.size(); attrCounter++) {
                if(!originalTypes.get(attrCounter).equals(phenoSub.attributeType(attrCounter))) {
                    attributesSame = false;
                }
            }
            Assert.assertTrue("Phenotype Subtract, Different Taxa, Same Variables With Covariates Test Fails. Attributes are not matching",attributesSame);
            
        }
        catch(Exception e) {
            Assert.assertTrue(false);
            e.printStackTrace();
        }
    }
    
  //With same Taxa and different Variables
    @Test
    public void testSubtractPhenoTestDifferentVariableWithCovariate() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNames2 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNames2 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var4", "var3"}));
        ArrayList<String> varNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> data2 = new ArrayList<ArrayList<Double>>();
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,100.0,100.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,300.0,300.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> dataKnownSub = new ArrayList<ArrayList<Double>>();
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{200.0,0.0})));
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{-200.0,0.0})));
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{0.0,0.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.covariate);
        
        ArrayList<ATTRIBUTE_TYPE> newTypes = new ArrayList<>();
        newTypes.add(ATTRIBUTE_TYPE.taxa);
        newTypes.add(ATTRIBUTE_TYPE.data);
        newTypes.add(ATTRIBUTE_TYPE.covariate);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"Dataset1",originalTypes);
            Phenotype pheno2 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames2, varNames2, data2,"Dataset2",originalTypes);
            Phenotype phenoKnownSub = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownSub, varNamesKnownSub, dataKnownSub,"Dataset3",newTypes);
            
            SubtractPhenotype subPhenotypeObject = new SubtractPhenotype();
            Phenotype phenoSub = subPhenotypeObject.subtractPhenotype(pheno1, pheno2, false);
            
            boolean isSame = true;
            for(int row = 0; row < phenoSub.getRowCount(); row++){
                for(int col = 0; col < phenoSub.getColumnCount(); col++) {
                    if(!phenoKnownSub.getValueAt(row, col).equals(phenoSub.getValueAt(row, col))) {
                        System.out.println(phenoKnownSub.getValueAt(row, col) + " "+phenoSub.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Subtract, Same Taxa, Different Variables With Covariates Test Fails",isSame);
            
            boolean attributesSame = true;
            for(int attrCounter = 0; attrCounter < newTypes.size(); attrCounter++) {
                if(!newTypes.get(attrCounter).equals(phenoSub.attributeType(attrCounter))) {
                    attributesSame = false;
                }
            }
            Assert.assertTrue("Phenotype Subtract, Same Taxa, Different Variables With Covariates Test Fails. Attributes are not matching",attributesSame);
            
        }
        catch(Exception e) {
            e.printStackTrace();
            Assert.assertTrue(false);
            
        }
    }
    
    //With different Taxa and different Variables
    @Test
    public void testSubtractPhenoTestDifferentTaxaAndVariableWithCov() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNames2 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML248","B37"}));
        ArrayList<String> taxaNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"B73","B37"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNames2 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var4", "var3"}));
        ArrayList<String> varNamesKnownSub = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> data2 = new ArrayList<ArrayList<Double>>();
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,100.0,100.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,300.0,300.0})));
        data2.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ArrayList<Double>> dataKnownSub = new ArrayList<ArrayList<Double>>();
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{200.0,0.0})));
        dataKnownSub.add(new ArrayList<Double>(Arrays.asList(new Double[]{0.0,0.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.covariate);
        
        ArrayList<ATTRIBUTE_TYPE> newTypes = new ArrayList<>();
        newTypes.add(ATTRIBUTE_TYPE.taxa);
        newTypes.add(ATTRIBUTE_TYPE.data);
        newTypes.add(ATTRIBUTE_TYPE.covariate);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"Dataset1",originalTypes);
            Phenotype pheno2 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames2, varNames2, data2,"Dataset2",originalTypes);
            Phenotype phenoKnownSub = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownSub, varNamesKnownSub, dataKnownSub,"DatasetSubtracted",newTypes);
            
            SubtractPhenotype subPhenotypeObject = new SubtractPhenotype();
            Phenotype phenoSub = subPhenotypeObject.subtractPhenotype(pheno1, pheno2, false);
            
            boolean isSame = true;
            for(int row = 0; row < phenoSub.getRowCount(); row++){
                for(int col = 0; col < phenoSub.getColumnCount(); col++) {
                    if(!phenoKnownSub.getValueAt(row, col).equals(phenoSub.getValueAt(row, col))) {
                        System.out.println(phenoKnownSub.getValueAt(row, col) + " "+phenoSub.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Subtract, Different Taxa, Different Variables With Covariates Test Fails",isSame);
            
            boolean attributesSame = true;
            for(int attrCounter = 0; attrCounter < newTypes.size(); attrCounter++) {
                if(!newTypes.get(attrCounter).equals(phenoSub.attributeType(attrCounter))) {
                    attributesSame = false;
                }
            }
            Assert.assertTrue("Phenotype Subtract, Different Taxa, Different Variables With Covariates Test Fails. Attributes are not matching",attributesSame);
            
        }
        catch(Exception e) {
            Assert.assertTrue(false);
            e.printStackTrace();
        }
    }
    
}
