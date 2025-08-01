package net.maizegenetics.analysis.b4r;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;

import net.maizegenetics.analysis.gobii.GOBIIPostgresConnection;
import net.maizegenetics.phenotype.CorePhenotype;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeUtils;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;

public class B4RPhenotypeUtils {
    public static Phenotype getPhenotypeFromB4RMergeDuplicates(String dbURL, String db, String userName, String password, 
                                            String studyName, String taxaName, String variableName) throws Exception {
        //Set up the connection
        Connection conn = setupB4RConnection(dbURL, db, userName, password);
        //Build the query
        String query = buildSQLQuery(studyName, taxaName,variableName);
        System.out.println(query);
        //Scan through the result set and build datasets
        ArrayList<ArrayList<String>> phenotypesArrayList = pullPhenotypes(conn, query, studyName, taxaName, variableName);
        System.out.println("Create the Phenotype Object");
        //printPhenos(phenotypesArrayList);
        //From datasets, build CorePhenotypeObject and return it
        Phenotype pheno = PhenotypeUtils.createPhenotypeObjectFromDB(phenotypesArrayList);
        return pheno;
    }
    public static Phenotype getPhenotypeFromB4RWithDuplicates(String dbURL, String db, String userName, String password, 
            String studyName, String taxaName, String variableName) throws Exception {
        //Set up the connection
        Connection conn = setupB4RConnection(dbURL, db, userName, password);
        //Build the query
        String query = buildSQLQuery(studyName, taxaName,variableName);
        System.out.println(query);
        //Scan through the result set and build datasets
        ArrayList<ArrayList<String>> phenotypesArrayList = pullPhenotypes(conn, query, studyName, taxaName, variableName);
        System.out.println("Create the Phenotype Object");
        //printPhenos(phenotypesArrayList);
        //From datasets, build CorePhenotypeObject and return it
        Phenotype pheno = PhenotypeUtils.createPhenotypeObjectFromDB2(phenotypesArrayList);
        return pheno;
    }
    
    public static Phenotype getPhenotypeFromB4RWithRepNumber(String dbURL, String db, String userName, String password, 
            String studyName, String taxaName, String variableName, String repNumber) throws Exception {
        //Set up the connection
        Connection conn = setupB4RConnection(dbURL, db, userName, password);
        //Build the query
        String query = buildSQLQueryWithRep(studyName, taxaName,variableName,repNumber);
        //Scan through the result set and build datasets
        ArrayList<ArrayList<String>> phenotypesArrayList = pullPhenotypes(conn, query, studyName, taxaName, variableName, repNumber);
        System.out.println("Create the Phenotype Object");
        //printPhenos(phenotypesArrayList);
        //From datasets, build CorePhenotypeObject and return it
        Phenotype pheno = PhenotypeUtils.createPhenotypeObjectFromDB2(phenotypesArrayList);
        return pheno;
    }
    
    //Get the B4R connection
    private static Connection setupB4RConnection(String dbURL, String db, String userName, String password) {
        return GOBIIPostgresConnection.connection(dbURL, userName, password, db); 
    }
    
    private static String buildSQLQuery(String studyName, String taxaName, String variableName) {
        String query = "SELECT S.title, E.product_name, P.plotno, Var.name, PD.value\n" +
                        "FROM operational.study AS S\n"+
                        "INNER JOIN operational.entry AS E ON S.id = E.study_id\n" +
                        "INNER JOIN operational.plot AS P ON S.id = P.study_id AND E.id = P.entry_id\n" +
                        "INNER JOIN operational.plot_data AS PD ON S.id = PD.study_id AND E.id = PD.entry_id AND P.id = PD.plot_id\n"+
                        "INNER JOIN master.variable AS Var ON PD.variable_id = Var.id\n";
        
        
        ArrayList<String> whereTerms = new ArrayList<String>();
        
        if(studyName.length() > 0) {
            whereTerms.add("S.title = ? ");
        }
        if(taxaName.length() > 0) {
            whereTerms.add("E.product_name = ? ");
        }
        if(variableName.length() > 0) {
            whereTerms.add("var.name = ? ");
        }
        
        if(whereTerms.size() == 0) {
             return query;
        }
        else {
            String whereComponent = whereTerms.stream().collect(Collectors.joining(" AND "));
            return query + " WHERE "+whereComponent;
        } 
    }
    
    private static String buildSQLQueryWithRep(String studyName, String taxaName, String variableName, String repNum) {
        String query = "SELECT S.title, E.product_name, P.plotno, Var.name, PD.value\n" +
                        "FROM operational.study AS S\n"+
                        "INNER JOIN operational.entry AS E ON S.id = E.study_id\n" +
                        "INNER JOIN operational.plot AS P ON S.id = P.study_id AND E.id = P.entry_id\n" +
                        "INNER JOIN operational.plot_data AS PD ON S.id = PD.study_id AND E.id = PD.entry_id AND P.id = PD.plot_id\n"+
                        "INNER JOIN master.variable AS Var ON PD.variable_id = Var.id\n";
        
        
        ArrayList<String> whereTerms = new ArrayList<String>();
        
        if(studyName.length() > 0) {
            whereTerms.add("S.title = ? ");
        }
        if(taxaName.length() > 0) {
            whereTerms.add("E.product_name = ? ");
        }
        if(variableName.length() > 0) {
            whereTerms.add("var.name = ? ");
        }
        if(repNum.length() > 0) {
            whereTerms.add("P.rep = ? ");
        }
        
        if(whereTerms.size() == 0) {
             return query;
        }
        else {
            String whereComponent = whereTerms.stream().collect(Collectors.joining(" AND "));
            return query + " WHERE "+whereComponent;
        } 
    }
    private static ArrayList<ArrayList<String>> pullPhenotypes(Connection connection, String query, 
                                                            String studyName, String taxaName, String variableName) throws Exception{
        ArrayList<ArrayList<String>> phenotypeObject = new ArrayList<ArrayList<String>>();
        //new list for title
        phenotypeObject.add(new ArrayList<String>());
        //new list for product_name
        phenotypeObject.add(new ArrayList<String>());
        //new list for plotno
        phenotypeObject.add(new ArrayList<String>());
        //new list for name
        phenotypeObject.add(new ArrayList<String>());
        //new list for value
        phenotypeObject.add(new ArrayList<String>());
        
        try {
            PreparedStatement phenoPS = connection.prepareStatement(query);
            //Prepare the statement
            int psCounter = 1;
            if(studyName.length() > 0) {
                phenoPS.setString(psCounter, studyName);
                psCounter++;
            }
            if(taxaName.length() > 0) {
                phenoPS.setString(psCounter, taxaName);
                psCounter++;
            }
            if(variableName.length() > 0) {
                phenoPS.setString(psCounter, variableName);
                psCounter++;
            }
            
            ResultSet phenoResults = phenoPS.executeQuery();
            System.out.println("Retrieved Result Set from DB");
            while(phenoResults.next()) {
                for(int i = 0; i < phenotypeObject.size(); i++) {
                    phenotypeObject.get(i).add(phenoResults.getString(i+1));
                }
            }
            
        }
        catch(SQLException e) {
            throw e;
        }
        
        return phenotypeObject;
    }
    
    private static ArrayList<ArrayList<String>> pullPhenotypes(Connection connection, String query, 
            String studyName, String taxaName, String variableName,String repNo) throws Exception{
        ArrayList<ArrayList<String>> phenotypeObject = new ArrayList<ArrayList<String>>();
        //new list for title
        phenotypeObject.add(new ArrayList<String>());
        //new list for product_name
        phenotypeObject.add(new ArrayList<String>());
        //new list for plotno
        phenotypeObject.add(new ArrayList<String>());
        //new list for name
        phenotypeObject.add(new ArrayList<String>());
        //new list for value
        phenotypeObject.add(new ArrayList<String>());
        
        try {
            PreparedStatement phenoPS = connection.prepareStatement(query);
            //Prepare the statement
            int psCounter = 1;
            if(studyName.length() > 0) {
                phenoPS.setString(psCounter, studyName);
                psCounter++;
            }
            if(taxaName.length() > 0) {
                phenoPS.setString(psCounter, taxaName);
                psCounter++;
            }
            if(variableName.length() > 0) {
                phenoPS.setString(psCounter, variableName);
                psCounter++;
            }
            if(repNo.length() > 0) {
                phenoPS.setInt(psCounter, Integer.parseInt(repNo));
                psCounter++;
            }
            
            ResultSet phenoResults = phenoPS.executeQuery();
            System.out.println("Retrieved Result Set from DB");
            while(phenoResults.next()) {
                for(int i = 0; i < phenotypeObject.size(); i++) {
                    phenotypeObject.get(i).add(phenoResults.getString(i+1));
                }
            }
            
        }
        catch(SQLException e) {
            throw e;
        }
        
            return phenotypeObject;
    }
    
    private static void printPhenos(ArrayList<ArrayList<String>> phenos) {
        System.out.println("tissueName\ttaxaName\tplotno\tvarName\tvalue");
        for(int i = 0; i < phenos.get(0).size(); i++) {
            for(int j = 0; j < phenos.size(); j++) {
                System.out.print(phenos.get(j).get(i) + "\t");
            }
            System.out.println();
        }
    }
    
    
}
