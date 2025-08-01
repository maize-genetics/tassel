package net.maizegenetics.analysis.numericaltransform;

import java.util.ArrayList;
import java.util.Arrays;

import org.junit.Test;

import junit.framework.Assert;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeUtils;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;

public class AvgPhenotypeByTaxaPluginTest {

    //Unit Test where all taxa are different
    @Test
    public void testAvgPhenoByTaxa() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNamesKnownAvg = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNamesKnownAvg = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        
        ArrayList<ArrayList<Double>> dataKnownAvg = new ArrayList<ArrayList<Double>>();
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"Dataset1",originalTypes);
            Phenotype phenoKnownAvg = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownAvg, varNamesKnownAvg, dataKnownAvg,"AverageDataset",originalTypes);
            
            
            AvgPhenotype avgPhenotypeObject = new AvgPhenotype();
            Phenotype phenoAvg = avgPhenotypeObject.averagePheno(pheno1, false,12345l);
            
           
            boolean isSame = true;
            for(int row = 0; row < phenoAvg.getRowCount(); row++){
                for(int col = 0; col < phenoAvg.getColumnCount(); col++) {
                    if(!phenoKnownAvg.getValueAt(row, col).equals(phenoAvg.getValueAt(row, col))) {
                        System.out.println(phenoKnownAvg.getValueAt(row, col) + " "+phenoAvg.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Average, All Unique Taxa, Test Fails",isSame);
        }
        catch(Exception e) {
            Assert.assertTrue(false);
            e.printStackTrace();
        }
    }
    
    //Unit Test where one taxa is a duplicate
    @Test
    public void testAvgPhenoByTaxaOneDuplicateTaxa() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37","B73"}));
        ArrayList<String> taxaNamesKnownAvg = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNamesKnownAvg = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,100.0,300.0})));
        
        ArrayList<ArrayList<Double>> dataKnownAvg = new ArrayList<ArrayList<Double>>();
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,150.0,200.0})));
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"Dataset1",originalTypes);
            Phenotype phenoKnownAvg = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownAvg, varNamesKnownAvg, dataKnownAvg,"AverageDataset",originalTypes);
            
            
            AvgPhenotype avgPhenotypeObject = new AvgPhenotype();
            Phenotype phenoAvg = avgPhenotypeObject.averagePheno(pheno1, false,12345l);
            
           
            boolean isSame = true;
            for(int row = 0; row < phenoAvg.getRowCount(); row++){
                for(int col = 0; col < phenoAvg.getColumnCount(); col++) {
                    if(!phenoKnownAvg.getValueAt(row, col).equals(phenoAvg.getValueAt(row, col))) {
                        System.out.println(phenoKnownAvg.getValueAt(row, col) + " "+phenoAvg.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Average, One Duplicate Taxa, Test Fails",isSame);
        }
        catch(Exception e) {
            e.printStackTrace();
            Assert.assertTrue(false);
            
        }
    }
    
  //Unit Test where all taxa are duplicate
    @Test
    public void testAvgPhenoByTaxaAllDuplicateTaxa() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","B73","B73","B73"}));
        ArrayList<String> taxaNamesKnownAvg = new ArrayList<String>(Arrays.asList(new String[]{"B73"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNamesKnownAvg = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,400.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,400.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,400.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,100.0,400.0})));
        
        ArrayList<ArrayList<Double>> dataKnownAvg = new ArrayList<ArrayList<Double>>();
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,250.0,400.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"Dataset1",originalTypes);
            Phenotype phenoKnownAvg = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownAvg, varNamesKnownAvg, dataKnownAvg,"AverageDataset",originalTypes);
            
            
            AvgPhenotype avgPhenotypeObject = new AvgPhenotype();
            Phenotype phenoAvg = avgPhenotypeObject.averagePheno(pheno1, false,12345l);
            
           
            boolean isSame = true;
            for(int row = 0; row < phenoAvg.getRowCount(); row++){
                for(int col = 0; col < phenoAvg.getColumnCount(); col++) {
                    if(!phenoKnownAvg.getValueAt(row, col).equals(phenoAvg.getValueAt(row, col))) {
                        System.out.println(phenoKnownAvg.getValueAt(row, col) + " "+phenoAvg.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Average, All Duplicate Taxa",isSame);
        }
        catch(Exception e) {
            Assert.assertTrue(false);
            e.printStackTrace();
        }
    }
    
    //Unit Test for attributes other than data
    @Test
    public void testAvgPhenoByTaxaWithCov() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        ArrayList<String> taxaNamesKnownAvg = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNamesKnownAvg = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        
        ArrayList<ArrayList<Double>> dataKnownAvg = new ArrayList<ArrayList<Double>>();
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.covariate);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"Dataset1",originalTypes);
            Phenotype phenoKnownAvg = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownAvg, varNamesKnownAvg, dataKnownAvg,"AverageDataset",originalTypes);
            
            
            AvgPhenotype avgPhenotypeObject = new AvgPhenotype();
            Phenotype phenoAvg = avgPhenotypeObject.averagePheno(pheno1, false,12345l);
            
           
            boolean isSame = true;
            for(int row = 0; row < phenoAvg.getRowCount(); row++){
                for(int col = 0; col < phenoAvg.getColumnCount(); col++) {
                    if(!phenoKnownAvg.getValueAt(row, col).equals(phenoAvg.getValueAt(row, col))) {
                        System.out.println(phenoKnownAvg.getValueAt(row, col) + " "+phenoAvg.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Average, All Unique Taxa, Test Fails",isSame);
            
            boolean attributesSame = true;
            for(int attrCounter = 0; attrCounter < originalTypes.size(); attrCounter++) {
                if(!originalTypes.get(attrCounter).equals(phenoAvg.attributeType(attrCounter))) {
                    attributesSame = false;
                }
            }
            Assert.assertTrue("Phenotype Average, All Unique Taxa, Test Fails. Attributes are not matching",attributesSame);
            
            
        }
        catch(Exception e) {
            Assert.assertTrue(false);
            e.printStackTrace();
        }
    }
  //Unit Test where one taxa is a duplicate
    @Test
    public void testAvgPhenoByTaxaOneDuplicateTaxaWithCov() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37","B73"}));
        ArrayList<String> taxaNamesKnownAvg = new ArrayList<String>(Arrays.asList(new String[]{"B73","CML246","B37"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNamesKnownAvg = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,100.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,100.0,300.0})));
        
        ArrayList<ArrayList<Double>> dataKnownAvg = new ArrayList<ArrayList<Double>>();
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,150.0,200.0})));
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,300.0})));
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,500.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.covariate);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"Dataset1",originalTypes);
            Phenotype phenoKnownAvg = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownAvg, varNamesKnownAvg, dataKnownAvg,"AverageDataset",originalTypes);
            
            
            AvgPhenotype avgPhenotypeObject = new AvgPhenotype();
            Phenotype phenoAvg = avgPhenotypeObject.averagePheno(pheno1, false,12345l);
            
           
            boolean isSame = true;
            for(int row = 0; row < phenoAvg.getRowCount(); row++){
                for(int col = 0; col < phenoAvg.getColumnCount(); col++) {
                    if(!phenoKnownAvg.getValueAt(row, col).equals(phenoAvg.getValueAt(row, col))) {
                        System.out.println(phenoKnownAvg.getValueAt(row, col) + " "+phenoAvg.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Average, One Duplicate Taxa, Test Fails",isSame);
            
            boolean attributesSame = true;
            for(int attrCounter = 0; attrCounter < originalTypes.size(); attrCounter++) {
                if(!originalTypes.get(attrCounter).equals(phenoAvg.attributeType(attrCounter))) {
                    attributesSame = false;
                }
            }
            Assert.assertTrue("Phenotype Average, One Duplicate Taxa, Test Fails. Attributes are not matching",attributesSame);
            
        }
        catch(Exception e) {
            e.printStackTrace();
            Assert.assertTrue(false);
            
        }
    }
    
  //Unit Test where all taxa are duplicate
    @Test
    public void testAvgPhenoByTaxaAllDuplicateTaxaWithCov() {
        ArrayList<String> taxaNames1 = new ArrayList<String>(Arrays.asList(new String[]{"B73","B73","B73","B73"}));
        ArrayList<String> taxaNamesKnownAvg = new ArrayList<String>(Arrays.asList(new String[]{"B73"}));
        
        ArrayList<String> varNames1 = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        ArrayList<String> varNamesKnownAvg = new ArrayList<String>(Arrays.asList(new String[]{"var1", "var2", "var3"}));
        
        ArrayList<ArrayList<Double>> data1 = new ArrayList<ArrayList<Double>>();
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,200.0,400.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{100.0,200.0,400.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{500.0,500.0,400.0})));
        data1.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,100.0,400.0})));
        
        ArrayList<ArrayList<Double>> dataKnownAvg = new ArrayList<ArrayList<Double>>();
        dataKnownAvg.add(new ArrayList<Double>(Arrays.asList(new Double[]{300.0,250.0,400.0})));
        
        ArrayList<ATTRIBUTE_TYPE> originalTypes = new ArrayList<>();
        originalTypes.add(ATTRIBUTE_TYPE.taxa);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.data);
        originalTypes.add(ATTRIBUTE_TYPE.covariate);
        
        try {
            Phenotype pheno1 = PhenotypeUtils.createPhenotypeFromTransform(taxaNames1, varNames1, data1,"Dataset1",originalTypes);
            Phenotype phenoKnownAvg = PhenotypeUtils.createPhenotypeFromTransform(taxaNamesKnownAvg, varNamesKnownAvg, dataKnownAvg,"AverageDataset",originalTypes);
            
            
            AvgPhenotype avgPhenotypeObject = new AvgPhenotype();
            Phenotype phenoAvg = avgPhenotypeObject.averagePheno(pheno1, false,12345l);
            
           
            boolean isSame = true;
            for(int row = 0; row < phenoAvg.getRowCount(); row++){
                for(int col = 0; col < phenoAvg.getColumnCount(); col++) {
                    if(!phenoKnownAvg.getValueAt(row, col).equals(phenoAvg.getValueAt(row, col))) {
                        System.out.println(phenoKnownAvg.getValueAt(row, col) + " "+phenoAvg.getValueAt(row, col));
                        isSame = false;
                    }
                }
            }
            
            Assert.assertTrue("Phenotype Average, All Duplicate Taxa",isSame);
            

            boolean attributesSame = true;
            for(int attrCounter = 0; attrCounter < originalTypes.size(); attrCounter++) {
                if(!originalTypes.get(attrCounter).equals(phenoAvg.attributeType(attrCounter))) {
                    attributesSame = false;
                }
            }
            Assert.assertTrue("Phenotype Average, All Duplicate Taxa. Attributes are not matching",attributesSame);
        }
        catch(Exception e) {
            Assert.assertTrue(false);
            e.printStackTrace();
        }
    }
}
