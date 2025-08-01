/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Properties;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.postgresql.copy.CopyManager;
import org.postgresql.core.BaseConnection;

import net.maizegenetics.analysis.data.GenomeAnnosDBQueryToPositionListPlugin;
import net.maizegenetics.util.Utils;

/**
 * These methods were largely borrowed from Jeff's privatemaizegenetics GenomeAnnosDB.java, 
 * then modified as needed for GOBII postgres access.  
 * 
 * Currently, this file is ONLY connecting to postgres, not monetd. if plugins are
 * created that need to connect to monetdb, a new connectToDB() method for monetdb
 * must be created, or the existing one will need to be changed to indicate the
 * database type for connection.
 * 
 * For connecting to a monetdb instance, see example in 
 * privatemaizegenetics.lynn.MonetdbFileProcessing.MonetDBQtoPosList.connectToDBOrDie
 * 
 * @author lcj34
 *
 */
public class GOBIIDbUtils {
    private static final Logger myLogger = LogManager.getLogger(GenomeAnnosDBQueryToPositionListPlugin.class);
    private static String errorMessage;
    public static Connection connectToDB(String configFile) {
        Properties props = new Properties();
        try {
            BufferedInputStream inputStream = new BufferedInputStream(new FileInputStream(configFile));
            props.load(inputStream);
        } catch (IOException e) {
            errorMessage = "Problem reading DB connection config file (" + configFile + "):\n\t" + e;
            myLogger.error(errorMessage);
            return null;
        }
        
        String user = props.getProperty("user");
        String password = props.getProperty("password");
        String dbName = props.getProperty("DB");
       // System.out.println("\nLCJ - connectToDB, props values are: " + props.toString()); // NOTE - this prints out the password
//        String port = props.getProperty("port");
        if (user == null) {
            errorMessage = "ERROR: Please provide a line with the user name (user=<userName>) in the DB connection config file (" + configFile + ")";
            myLogger.error(errorMessage);
            return null;
        }
        if (password == null) {
            errorMessage = "ERROR: Please provide a line with the password (password=<yourPassword>) in the DB connection config file (" + configFile + ")";
            myLogger.error(errorMessage);
            return null;
        }
        if (dbName == null) {
            errorMessage = "ERROR: Please provide a line with the DB name (DB=<dbName>) in the DB connection config file (" + configFile + ")";
            myLogger.error(errorMessage);
            return null;
        }
        return connectToDatabaseOrDie(dbName, props);
    }
    public static Connection connectToDatabaseOrDie(String dbName, Properties props) {
        Connection conn = null;
        String url = "not connected yet";
        String host = props.getProperty("host");
        //String host = null;
        try {
            Class.forName("org.postgresql.Driver");
            if (host == null) {
                url = "jdbc:postgresql://localhost:5432/" + dbName;
            } else {
                url = "jdbc:postgresql://"+host+":5432/" + dbName;
            }
            String user = props.getProperty("user");
            System.out.println("Attempting connection with user " + user + " and url " + url);
            conn = DriverManager.getConnection(url, props);
        } catch (ClassNotFoundException e) {
            errorMessage = e.getMessage();
            myLogger.error(errorMessage);
            return null;
        } catch (SQLException e) {
            errorMessage = e.getMessage();
            myLogger.error(errorMessage);
            return null;
        }
        myLogger.info("\nUsing DB:  " + url + "\n");
        return conn;
    }
    public static String[] readNextLine(BufferedReader reader, String inFile) {
        String line;
        try {
            line = reader.readLine();
        } catch (IOException e) {
            System.err.println("\n\nProblem reading data file (" + inFile + "):\n\t" + e);
            System.exit(1);
            return null;
        }
        if (line == null) {
            return null;
        }
        return line.split("\t", -1);  // negative limit (-1) means no limit on the number of fields and trailing blank fields are not discarded
    }
    
    // The CopyManager is the API for PostgreSQL COPY bulk data transfer.
    // using "," as the delimiter (csv files)
    public static void postgreSQLCopyFromReader(Connection conn, String table, BufferedReader reader) {
        try {
            CopyManager copyManager = new CopyManager((BaseConnection) conn);            
            copyManager.copyIn("COPY " + table + " FROM STDIN", reader);
         // Specify the comma as delimiter.  Default is tab-delimited
            //copyManager.copyIn("COPY " + table + " FROM STDIN WITH DELIMITER ','", reader);
        } catch (SQLException sqle) {
            System.err.println("\n\nProblem populating table from file (COPY " + table + " FROM STDIN):\n\t" + sqle);
            System.exit(1);
        } catch (FileNotFoundException fnfe) {
            System.err.println("\n\n" + fnfe);
            System.exit(1);
        } catch (IOException ioe) {
            System.err.println("\n\n" + ioe);
            System.exit(1);
        }
    }
    public static ResultSet executePostgreSQLQuery(Statement st, String query, boolean echoQuery) {
        if (echoQuery) {
            System.out.println("\n\n" + query);
        }
        try {
            return st.executeQuery(query);
        } catch (SQLException e) {
            System.err.println("\n\nProblem executing query (" + query + "):\n\t" + e);
            System.exit(1);
            return null;
        }
    }
    public static void populateTableFromFile(Connection conn, String table, String sourceFile, boolean header) {
        System.out.println("\n\nPopulating the "+table+" table from the tab-delimited text file:\n   " + sourceFile);
        BufferedReader reader = Utils.getBufferedReader(sourceFile, 524288);
        if (header) {
            readNextLine(reader, sourceFile);
        }
        postgreSQLCopyFromReader(conn, table, reader);
        printPostgreSQLResultSet(executePostgreSQLQuery(createPostgreSQLStatement(conn), "SELECT count(*) from "+table, true));
    }
    
    public static Statement createPostgreSQLStatement(Connection conn) {
        try {
            return conn.createStatement();
        } catch (SQLException e) {
            System.err.println("\n\nProblem creating statement:\n\t" + e);
            System.exit(1);
            return null;
        }
    }

    public static void printPostgreSQLResultSet(ResultSet rs) {
        printPostgreSQLResultSet(rs, 1);
    }
    public static void printPostgreSQLResultSet(ResultSet rs, int stride) {
        try {
            ResultSetMetaData rsmd = rs.getMetaData();
            int[] type = new int[rsmd.getColumnCount()+1];
            for (int col = 1; col <= rsmd.getColumnCount(); col++) {
                type[col] = rsmd.getColumnType(col);
                System.out.print(rsmd.getColumnLabel(col));
                if (col < rsmd.getColumnCount()) {
                    System.out.print("\t");
                }
            }
            System.out.print("\n");
            int nResults = 0;
            while (rs.next()) {
                if (nResults % stride == 0) {
                    for (int col = 1; col <= rsmd.getColumnCount(); col++) {
                        if (type[col] == java.sql.Types.REAL || type[col] == java.sql.Types.FLOAT) {
                            System.out.print(rs.getFloat(col));
                        } else if (type[col] == java.sql.Types.DOUBLE) {
                            System.out.print(rs.getDouble(col));
                        } else {
                            System.out.print(rs.getString(col));
                        }
                        if (col < rsmd.getColumnCount()) {
                            System.out.print("\t");
                        }
                    }
                    System.out.print("\n");
                }
                nResults++;
            }
            System.out.print("\n");
            rs.close();
        } catch (SQLException se) {
            System.err.println(se.getMessage());
        }
    }
    public static void printPostgreSQLResultSetToFile(ResultSet rs, String outFile) {
        BufferedWriter writer = Utils.getBufferedWriter(outFile);
        StringBuilder sb = new StringBuilder();
        int nRows = 0;
        try {
            ResultSetMetaData rsmd = rs.getMetaData();
            int[] type = new int[rsmd.getColumnCount()+1];
            for (int col = 1; col <= rsmd.getColumnCount(); col++) {
                type[col] = rsmd.getColumnType(col);
                sb.append(rsmd.getColumnLabel(col));
                if (col < rsmd.getColumnCount()) {
                    sb.append("\t");
                }
            }
            sb.append("\n");
            while (rs.next()) {
                for (int col = 1; col <= rsmd.getColumnCount(); col++) {
                    if (type[col] == java.sql.Types.REAL || type[col] == java.sql.Types.FLOAT) {
                        sb.append(rs.getFloat(col));
                    } else if (type[col] == java.sql.Types.DOUBLE) {
                        sb.append(rs.getDouble(col));
                    } else {
                        sb.append(rs.getString(col));
                    }
                    if (col < rsmd.getColumnCount()) {
                        sb.append("\t");
                    }
                }
                sb.append("\n");
                nRows++;
                if (nRows % 1000 == 0) {
                    writer.append(sb.toString());
                    sb = new StringBuilder();
                }
            }
            writer.append(sb.toString());
            writer.close();
            rs.close();
        } catch (SQLException se) {
            System.err.println(se.getMessage());
        } catch (IOException ioe) {
            System.err.println(ioe.getMessage());
        }
    }

}
