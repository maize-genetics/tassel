package net.maizegenetics.util;

import java.util.ArrayList;

public class BuilderFromVCFUtil {

    //This method is designed to import multiple indels on a single VCF line
    public static char[][] createSmallCallTable(ArrayList<String> variants) {
        
        //Calculate the longest variant sequence.  This is how big we will need each char array to be
        int maxVarLength = 0;
        for(String variant : variants) {
            maxVarLength = (variant.length() > maxVarLength)?variant.length() : maxVarLength;
        }
        //Create a call table
        char[][] alleleCallTable = new char[variants.size()][maxVarLength];
       
        for(int varCounter = 0; varCounter < alleleCallTable.length; varCounter++ ) {
            for(int posCounter = 0; posCounter < alleleCallTable[varCounter].length; posCounter++) {
                if(variants.get(varCounter).length()<=posCounter) {
                    alleleCallTable[varCounter][posCounter] = '-';
                }
                else {
                    alleleCallTable[varCounter][posCounter] = variants.get(varCounter).charAt(posCounter);
                }
            }
        }
        return alleleCallTable;
    }
}
