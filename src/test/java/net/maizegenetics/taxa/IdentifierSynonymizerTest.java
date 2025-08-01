package net.maizegenetics.taxa;

import static org.junit.Assert.*;

import org.junit.Ignore;
import org.junit.Test;

import net.maizegenetics.taxa.IdentifierSynonymizer;

public class IdentifierSynonymizerTest {

    @Test
    public void diceTest() {
        int technique = 0;//For dice's coefficient
        String methodName = "Dice";
        //Initialize String vars;
        String testString1 = "";
        String testString2 = "";
            
        //NOTE Higher numbers = lower similarity.  We are dealing with distances here.  Expected values calculated by hand.
        //Check to see if both strings are empty
        assertEquals(methodName+" Test Fail for null string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
        //Check to see if one string is empty
        testString1 = "Hello World";
        testString2 = "";
        assertEquals(methodName+" Test Fail for single null string comparison: Expected: 1.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),1.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
        
        //Check to make sure if both strings are the same we get 1.0
        testString1 = "Hello World";
        testString2 = "Hello World";
        assertEquals(methodName+" Test Fail for same string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
    
        //Check to see if the expected value is returned from different strings
        testString1 = "Hello World";
        testString2 = "Hello Tassel";
        assertEquals(methodName+" Test Fail for different string comparison: Expected: 0.5789, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.5789,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
    }
    
    @Test
    public void editDistTest() {
        int technique = 1;//For edit Distance
        String methodName = "EditDistance";
        //Initialize String vars;
        String testString1 = "";
        String testString2 = "";
            
        //NOTE Higher numbers = lower similarity.  We are dealing with distances here. Expected Values are calculated by hand
        //Check to see if both strings are empty
        assertEquals(methodName+" Test Fail for null string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
     
        //Check to see if one string is empty
        testString1 = "Hello World";
        testString2 = "";
        assertEquals(methodName+" Test Fail for single null string comparison: Expected: 10.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),10.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
       
        
        //Check to make sure if both strings are the same we get 1.0
        testString1 = "Hello World";
        testString2 = "Hello World";
        assertEquals(methodName+" Test Fail for same string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
    
        //Check to see if the expected value is returned from different strings
        testString1 = "Hello World";
        testString2 = "Hello Tassel";
        assertEquals(methodName+" Test Fail for different string comparison: Expected: 6.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),6.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
    }
    
    @Test
    public void dtwHammingDistTest() {
        int technique = 2;//For DTW with hamming
        String methodName = "Dynamic Time Warping With Hamming Distance";
        //Initialize String vars;
        String testString1 = "";
        String testString2 = "";
            
        //NOTE Higher numbers = lower similarity.  We are dealing with distances here. Expected Values are calculated by hand
        //Check to see if both strings are empty
        assertEquals(methodName+" Test Fail for null string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
        //Check to see if one string is empty
        testString1 = "Hello World";
        testString2 = "";
        
        try {
            IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique);
            fail(methodName+" Test Fail for single null string comparison: Expected: ArrayIndexOutOfBoundsException, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique));
        }
        catch(IndexOutOfBoundsException e) {
            assertTrue(true);
        }
        
        
        //Check to make sure if both strings are the same we get 1.0
        testString1 = "Hello World";
        testString2 = "Hello World";
        assertEquals(methodName+" Test Fail for same string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
    
        //Check to see if the expected value is returned from different strings
        testString1 = "Hello World";
        testString2 = "Hello Tassel";
        assertEquals(methodName+" Test Fail for different string comparison: Expected: 6.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),6.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
    }
    
    @Test
    public void dtwKeyboardDistTest() {
        int technique = 3;//For DTW with Keyboard
        String methodName = "Dynamic Time Warping With Keyboard Distance";
        //Initialize String vars;
        String testString1 = "";
        String testString2 = "";
            
        //NOTE Higher numbers = lower similarity.  We are dealing with distances here. Expected Values are calculated by hand
        //Check to see if both strings are empty
        assertEquals(methodName+" Test Fail for null string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
        //Check to see if one string is empty
        testString1 = "Hello World";
        testString2 = "";
        
        try {
            IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique);
            fail(methodName+" Test Fail for single null string comparison: Expected: ArrayIndexOutOfBoundsException, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique));
        }
        catch(IndexOutOfBoundsException e) {
            assertTrue(true);
        }
        
        
        //Check to make sure if both strings are the same we get 1.0
        testString1 = "Hello World";
        testString2 = "Hello World";
        assertEquals(methodName+" Test Fail for same string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
    
        //Check to see if the expected value is returned from different strings
        testString1 = "Hello World";
        testString2 = "Hello Tassel";
        assertEquals(methodName+" Test Fail for different string comparison: Expected: 20.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),20.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
    }
    
    
    @Test
    public void hammingSoundexTest() {
        int technique = 4;//For DTW with Keyboard
        String methodName = "Hamming Distance Using Soundex Encoding";
        //Initialize String vars;
        String testString1 = "";
        String testString2 = "";
            
        //NOTE Higher numbers = lower similarity.  We are dealing with distances here. Expected Values are calculated by hand
        //Check to see if both strings are empty
        assertEquals(methodName+" Test Fail for null string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
        //Check to see if one string is empty
        testString1 = "Hello World";
        testString2 = "";
        assertEquals(methodName+" Test Fail for null string comparison: Expected: 4.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),4.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
        //Check to make sure if both strings are the same we get 1.0
        testString1 = "Hello World";
        testString2 = "Hello World";
        assertEquals(methodName+" Test Fail for same string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
    
        //Check to see if the expected value is returned from different strings
        testString1 = "Hello World";
        testString2 = "Hello Tassel";
        //System.out.println(IdentifierSynonymizer.soundex2(testString1,true,true,true));
        //System.out.println(IdentifierSynonymizer.soundex2(testString2, true, true, true));
        assertEquals(methodName+" Test Fail for different string comparison: Expected: 2, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),2.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
    }
    
    @Test
    public void diceWithMetaphoneTest() {
        int technique = 5;//For dice's coefficient with Metaphone
        String methodName = "Dice With Metaphone";
        //Initialize String vars;
        String testString1 = "";
        String testString2 = "";
            
        //NOTE Higher numbers = lower similarity.  We are dealing with distances here.  Expected values calculated by hand.
        //Check to see if both strings are empty
        assertEquals(methodName+" Test Fail for null string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
        //Check to see if one string is empty
        testString1 = "Hello World";
        testString2 = "";
        assertEquals(methodName+" Test Fail for single null string comparison: Expected: 1.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),1.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
        
        //Check to make sure if both strings are the same we get 0.0
        testString1 = "Hello World";
        testString2 = "Hello World";
        assertEquals(methodName+" Test Fail for same string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
    
        //Check to see if the expected value is returned from different strings
        testString1 = "Hello World";
        testString2 = "Hello Tassel";
        assertEquals(methodName+" Test Fail for different string comparison: Expected: 0.6666, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.6666,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
    }
    
    @Test
    public void editWithMetaphoneTest() {
        int technique = 6;//For dice's coefficient with Metaphone
        String methodName = "Edit Distance With Metaphone";
        //Initialize String vars;
        String testString1 = "";
        String testString2 = "";
            
        //NOTE Higher numbers = lower similarity.  We are dealing with distances here.  Expected values calculated by hand.
        //Check to see if both strings are empty
        assertEquals(methodName+" Test Fail for null string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
        //Check to see if one string is empty
        testString1 = "Hello World";
        testString2 = "";
        //assertEquals(methodName+" Test Fail for single null string comparison: Expected: 1.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),1.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        assertEquals(methodName+" Test Fail for single null string comparison: Expected: 4.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),4.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
        
        
        //Check to make sure if both strings are the same we get 0.0
        testString1 = "Hello World";
        testString2 = "Hello World";
        assertEquals(methodName+" Test Fail for same string comparison: Expected: 0.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),0.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001); 
    
        //Check to see if the expected value is returned from different strings
        testString1 = "Hello World";
        testString2 = "Hello Tassel";
        assertEquals(methodName+" Test Fail for different string comparison: Expected: 2.0, got:"+IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),2.0,IdentifierSynonymizer.getScore(testString1, testString2, true, true, true, technique),.0001);   
    }
    

}
