/**
 * 
 */
package net.maizegenetics.analysis.data;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import org.junit.Test;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;

/**
 * @author lcj34
 * @author Josh Guy
 *
 */
public class GenosToABHPluginTest{
    
    /**
     * Test not specifying output format
     */
    @Test
    public void GenosToABHPluginCallSmallGenosOutNULL() {
 
        String actualResultsFile = GeneralConstants.TEMP_DIR + "ABHTestCaseJunitOut.csv";
        String expectedResultsFile = GeneralConstants.DATA_DIR + "GenotypeTableTests/genosToABH_expectedOutC.csv";
        String parentA = GeneralConstants.DATA_DIR + "GenotypeTableTests/parentAFile.txt";
        String parentB = GeneralConstants.DATA_DIR + "GenotypeTableTests/parentBFile.txt";
        String inputFile = GeneralConstants.DATA_DIR + "GenotypeTableTests/genosToABH.hmp.txt";
        
        System.out.println("Starting test GenosToABHPluginCallSmallGenosOutNULL");
        
        DataSet result = ImportUtils.readDataSet(inputFile);
        new GenosToABHPlugin()
        .outfile(actualResultsFile)
        .parentA(parentA)
        .parentB(parentB)
        .performFunction(result);
        
        
        Multimap<String, Tuple<String,String>> actualResults = parseResultsFile(actualResultsFile);
        Multimap<String, Tuple<String,String>> expectedResults = parseResultsFile( expectedResultsFile);
        if (actualResults == null || expectedResults == null) {
            System.out.println("GenosToABHPluginCallSmallGenosOutNULL FAILED - bad files");
            return;
        }       
        assertEquals(actualResults.size(), expectedResults.size());
        
        // Get the key sets for the actual and expected results                 
        Collection aTaxonColl = (Collection) actualResults.keySet();
        Collection eTaxonColl = (Collection) expectedResults.keySet();
        assertEquals(aTaxonColl.size(),eTaxonColl.size());
        
        System.out.println("aTaxonColl size: "+ aTaxonColl.size() + " eTaxonColl size: " + eTaxonColl.size());
        // For each taxon in the actualResults file, verify it appears in the expectedResults file

        actualResults.keySet().stream().forEach(aTaxon -> {
            assertTrue(eTaxonColl.contains(aTaxon));
            // Get the values for this taxon as stored in the actual and expected results files
            // Verify each taxon contains the same sites/allele calls
            Collection<Tuple<String, String>> aValues = (Collection) actualResults.get(aTaxon);
            Collection<Tuple<String, String>> eValues = (Collection) expectedResults.get(aTaxon);
            for (Tuple<String, String> myTuple : aValues) {
                boolean containsTuple = eValues.contains(myTuple);
                System.out.println("Contained in Expected List: " + containsTuple + " ActualTuple X: " + myTuple.x + " Y: " + myTuple.y);                
                assertTrue(eValues.contains(myTuple));
            }           
        });
        
        System.out.println("\n");
        // Verify all taxon in expectedResults occurs in actualResults
        expectedResults.keySet().stream().forEach(eTaxon -> {
            assertTrue(eTaxonColl.contains(eTaxon));
            // Get the values for this taxon as stored in the actual and expected results files
            // Verify each taxon contains the same sites/allele calls
            Collection<Tuple<String, String>> aValues = (Collection) actualResults.get(eTaxon);
            Collection<Tuple<String, String>> eValues = (Collection) expectedResults.get(eTaxon);
            for (Tuple<String, String> myTuple : eValues) {
                boolean containsTuple = aValues.contains(myTuple);
                System.out.println("Contained in Actual List: " + containsTuple + " ExpectedTuple X: " + myTuple.x + " Y: " + myTuple.y);                
                assertTrue(aValues.contains(myTuple));
            }           
        });
        
        System.out.println("Test Complete !!");
    }
    
    /**
     * Test given ABH output format
     */
    @Test
    public void GenosToABHPluginCallSmallGenosOutC() {
 
        String actualResultsFile = GeneralConstants.TEMP_DIR + "ABHTestCaseJunitOut.csv";
        String expectedResultsFile = GeneralConstants.DATA_DIR + "GenotypeTableTests/genosToABH_expectedOutC.csv";
        String parentA = GeneralConstants.DATA_DIR + "GenotypeTableTests/parentAFile.txt";
        String parentB = GeneralConstants.DATA_DIR + "GenotypeTableTests/parentBFile.txt";
        String inputFile = GeneralConstants.DATA_DIR + "GenotypeTableTests/genosToABH.hmp.txt";
        
        System.out.println("Starting test GenosToABHPluginCallSmallGenosOutC");
        
        DataSet result = ImportUtils.readDataSet(inputFile);
        new GenosToABHPlugin()
        .outfile(actualResultsFile)
        .parentA(parentA)
        .parentB(parentB)
        .outputFormat(GenosToABHPlugin.OUTPUT_CHECK.c)
        .performFunction(result);
        
        
        Multimap<String, Tuple<String,String>> actualResults = parseResultsFile(actualResultsFile);
        Multimap<String, Tuple<String,String>> expectedResults = parseResultsFile( expectedResultsFile);
        if (actualResults == null || expectedResults == null) {
            System.out.println("GenosToABHPluginCallSmallGenosOutC FAILED - bad files");
            return;
        }       
        assertEquals(actualResults.size(), expectedResults.size());
        
        // Get the key sets for the actual and expected results                 
        Collection aTaxonColl = (Collection) actualResults.keySet();
        Collection eTaxonColl = (Collection) expectedResults.keySet();
        assertEquals(aTaxonColl.size(),eTaxonColl.size());
        
        System.out.println("aTaxonColl size: "+ aTaxonColl.size() + " eTaxonColl size: " + eTaxonColl.size());
        // For each taxon in the actualResults file, verify it appears in the expectedResults file

        actualResults.keySet().stream().forEach(aTaxon -> {
            assertTrue(eTaxonColl.contains(aTaxon));
            // Get the values for this taxon as stored in the actual and expected results files
            // Verify each taxon contains the same sites/allele calls
            Collection<Tuple<String, String>> aValues = (Collection) actualResults.get(aTaxon);
            Collection<Tuple<String, String>> eValues = (Collection) expectedResults.get(aTaxon);
            for (Tuple<String, String> myTuple : aValues) {
                boolean containsTuple = eValues.contains(myTuple);
                System.out.println("Contained in Expected List: " + containsTuple + " ActualTuple X: " + myTuple.x + " Y: " + myTuple.y);                
                assertTrue(eValues.contains(myTuple));
            }           
        });
        
        System.out.println("\n");
        // Verify all taxon in expectedResults occurs in actualResults
        expectedResults.keySet().stream().forEach(eTaxon -> {
            assertTrue(eTaxonColl.contains(eTaxon));
            // Get the values for this taxon as stored in the actual and expected results files
            // Verify each taxon contains the same sites/allele calls
            Collection<Tuple<String, String>> aValues = (Collection) actualResults.get(eTaxon);
            Collection<Tuple<String, String>> eValues = (Collection) expectedResults.get(eTaxon);
            for (Tuple<String, String> myTuple : eValues) {
                boolean containsTuple = aValues.contains(myTuple);
                System.out.println("Contained in Actual List: " + containsTuple + " ExpectedTuple X: " + myTuple.x + " Y: " + myTuple.y);                
                assertTrue(aValues.contains(myTuple));
            }           
        });
        
        System.out.println("Test Complete !!");
    }

    /**
     * Test given integer output format
     */
    @Test
    public void GenosToABHPluginCallSmallGenosOutI() {
        
        String actualResultsFile = GeneralConstants.TEMP_DIR + "ABHTestCaseJunitOut.csv";
        String expectedResultsFile = GeneralConstants.DATA_DIR + "GenotypeTableTests/genosToABH_expectedOutI.csv";
        String parentA = GeneralConstants.DATA_DIR + "GenotypeTableTests/parentAFile.txt";
        String parentB = GeneralConstants.DATA_DIR + "GenotypeTableTests/parentBFile.txt";
        String inputFile = GeneralConstants.DATA_DIR + "GenotypeTableTests/genosToABH.hmp.txt";
        
        System.out.println("Starting test GenosToABHPluginCallSmallGenosOutI");
        
        DataSet result = ImportUtils.readDataSet(inputFile);
        new GenosToABHPlugin()
        .outfile(actualResultsFile)
        .parentA(parentA)
        .parentB(parentB)
        .outputFormat(GenosToABHPlugin.OUTPUT_CHECK.i)
        .performFunction(result);
        
        
        Multimap<String, Tuple<String,String>> actualResults = parseResultsFile(actualResultsFile);
        Multimap<String, Tuple<String,String>> expectedResults = parseResultsFile( expectedResultsFile);
        if (actualResults == null || expectedResults == null) {
            System.out.println("GenosToABHPluginCallSmallGenosOutI FAILED - bad files");
            return;
        }       
        assertEquals(actualResults.size(), expectedResults.size());
        
        // Get the key sets for the actual and expected results                 
        Collection aTaxonColl = (Collection) actualResults.keySet();
        Collection eTaxonColl = (Collection) expectedResults.keySet();
        assertEquals(aTaxonColl.size(),eTaxonColl.size());
        
        System.out.println("aTaxonColl size: "+ aTaxonColl.size() + " eTaxonColl size: " + eTaxonColl.size());
        // For each taxon in the actualResults file, verify it appears in the expectedResults file

        actualResults.keySet().stream().forEach(aTaxon -> {
            assertTrue(eTaxonColl.contains(aTaxon));
            // Get the values for this taxon as stored in the actual and expected results files
            // Verify each taxon contains the same sites/allele calls
            Collection<Tuple<String, String>> aValues = (Collection) actualResults.get(aTaxon);
            Collection<Tuple<String, String>> eValues = (Collection) expectedResults.get(aTaxon);
            for (Tuple<String, String> myTuple : aValues) {
                boolean containsTuple = eValues.contains(myTuple);
                System.out.println("Contained in Expected List: " + containsTuple + " ActualTuple X: " + myTuple.x + " Y: " + myTuple.y);                
                assertTrue(eValues.contains(myTuple));
            }           
        });
        
        System.out.println("\n");
        // Verify all taxon in expectedResults occurs in actualResults
        expectedResults.keySet().stream().forEach(eTaxon -> {
            assertTrue(eTaxonColl.contains(eTaxon));
            // Get the values for this taxon as stored in the actual and expected results files
            // Verify each taxon contains the same sites/allele calls
            Collection<Tuple<String, String>> aValues = (Collection) actualResults.get(eTaxon);
            Collection<Tuple<String, String>> eValues = (Collection) expectedResults.get(eTaxon);
            for (Tuple<String, String> myTuple : eValues) {
                boolean containsTuple = aValues.contains(myTuple);
                System.out.println("Contained in Actual List: " + containsTuple + " ExpectedTuple X: " + myTuple.x + " Y: " + myTuple.y);                
                assertTrue(aValues.contains(myTuple));
            }           
        });
        
        System.out.println("Test Complete !!");
    }
   
    /**
     * Test given real output format
     */
    @Test
    public void GenosToABHPluginCallSmallGenosOutR() {
        
        String actualResultsFile = GeneralConstants.TEMP_DIR + "ABHTestCaseJunitOut.csv";
        String expectedResultsFile = GeneralConstants.DATA_DIR + "GenotypeTableTests/genosToABH_expectedOutR.csv";
        String parentA = GeneralConstants.DATA_DIR + "GenotypeTableTests/parentAFile.txt";
        String parentB = GeneralConstants.DATA_DIR + "GenotypeTableTests/parentBFile.txt";
        String inputFile = GeneralConstants.DATA_DIR + "GenotypeTableTests/genosToABH.hmp.txt";
        
        System.out.println("Starting test GenosToABHPluginCallSmallGenosOutR");
        
        DataSet result = ImportUtils.readDataSet(inputFile);
        new GenosToABHPlugin()
        .outfile(actualResultsFile)
        .parentA(parentA)
        .parentB(parentB)
        .outputFormat(GenosToABHPlugin.OUTPUT_CHECK.r)
        .performFunction(result);
        
        
        Multimap<String, Tuple<String,String>> actualResults = parseResultsFile(actualResultsFile);
        Multimap<String, Tuple<String,String>> expectedResults = parseResultsFile( expectedResultsFile);
        if (actualResults == null || expectedResults == null) {
            System.out.println("GenosToABHPluginCallSmallGenosOutR FAILED - bad files");
            return;
        }       
        assertEquals(actualResults.size(), expectedResults.size());
        
        // Get the key sets for the actual and expected results                 
        Collection aTaxonColl = (Collection) actualResults.keySet();
        Collection eTaxonColl = (Collection) expectedResults.keySet();
        assertEquals(aTaxonColl.size(),eTaxonColl.size());
        
        System.out.println("aTaxonColl size: "+ aTaxonColl.size() + " eTaxonColl size: " + eTaxonColl.size());
        // For each taxon in the actualResults file, verify it appears in the expectedResults file

        actualResults.keySet().stream().forEach(aTaxon -> {
            assertTrue(eTaxonColl.contains(aTaxon));
            // Get the values for this taxon as stored in the actual and expected results files
            // Verify each taxon contains the same sites/allele calls
            Collection<Tuple<String, String>> aValues = (Collection) actualResults.get(aTaxon);
            Collection<Tuple<String, String>> eValues = (Collection) expectedResults.get(aTaxon);
            for (Tuple<String, String> myTuple : aValues) {
                boolean containsTuple = eValues.contains(myTuple);
                System.out.println("Contained in Expected List: " + containsTuple + " ActualTuple X: " + myTuple.x + " Y: " + myTuple.y);                
                assertTrue(eValues.contains(myTuple));
            }           
        });
        
        System.out.println("\n");
        // Verify all taxon in expectedResults occurs in actualResults
        expectedResults.keySet().stream().forEach(eTaxon -> {
            assertTrue(eTaxonColl.contains(eTaxon));
            // Get the values for this taxon as stored in the actual and expected results files
            // Verify each taxon contains the same sites/allele calls
            Collection<Tuple<String, String>> aValues = (Collection) actualResults.get(eTaxon);
            Collection<Tuple<String, String>> eValues = (Collection) expectedResults.get(eTaxon);
            for (Tuple<String, String> myTuple : eValues) {
                boolean containsTuple = aValues.contains(myTuple);
                System.out.println("Contained in Actual List: " + containsTuple + " ExpectedTuple X: " + myTuple.x + " Y: " + myTuple.y);                
                assertTrue(aValues.contains(myTuple));
            }           
        });
        
        System.out.println("Test Complete !!");
    }
       
    private Multimap<String, Tuple<String,String>> parseResultsFile(String resultsFile) {
        // put the data into a hashmap with taxon as the key, then a tuple of site/allele 
        Multimap<String, Tuple<String,String>> resultsMap = HashMultimap.create();
        
        ArrayList<String> aSiteName = new ArrayList<String>();
        
        try (BufferedReader reader = Utils.getBufferedReader(resultsFile)) {
            String line;
            // Read each line of the parentFile, add the parent taxon
            // to the ArrayList of parent names
            String line1 = reader.readLine(); // first line should be "id" with site names
            if (line1 == null) {
                System.out.println("No lines in output file !!");
                return null;                
            }
            String line2 = reader.readLine(); // second line should be chromosomes
            if (line2 == null) {
                System.out.println("No second line in output file !!");
                return null;                
            }
            String[] siteTokens = line1.split(",");          
            for (int idx = 1; idx < siteTokens.length; idx++) {
                aSiteName.add(siteTokens[idx]); // skip first element, which is header "id"
            }
 
            // Remaining lines are the taxon names with allele calls at sites
            while ((line = reader.readLine()) != null) {              
               // Add the taxon with site/allele value to the map
                String[] tokens = line.split(",");
                String taxaName = tokens[0];
                for (int idx = 1; idx < tokens.length; idx++) {
                    String siteName = aSiteName.get(idx-1); // -1 because of the header skip above
                    Tuple <String,String> siteValue = new Tuple<String, String>(siteName, tokens[idx]);
                    resultsMap.put(taxaName, siteValue);
                }
                System.out.println("parseResultsFile: Number of sites actualResults: " + aSiteName.size() + 
                        " size of ActualResults: " + resultsMap.size());
            }
            
        } catch (IOException exc) {
            exc.printStackTrace();
        }   
        return resultsMap;
    }
    
//  @Test
//  public void GenosToABHPluginCallLargeData() {
//      new GenosToABHPlugin()
//      .infile("/Users/lcj34/development/hackathon_2015/filt_informsites.hmp.txt")
//      .outfile("/Users/lcj34/development/hackathon_2015/genosToABH_noParents.csv")
//      .parentA("/Users/lcj34/development/hackathon_2015/parentAFile.txt")
//      .parentB("/Users/lcj34/development/hackathon_2015/parentBFile.txt")
//      .performFunction(null);
//  }
}
