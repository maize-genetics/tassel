/*
 *  BMSConnection
 * 
 *  Created on Aug 5, 2016
 */
package net.maizegenetics.analysis.gobii;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Properties;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class BMSConnection {

    private static final Logger myLogger = LogManager.getLogger(BMSConnection.class);

    private BMSConnection() {
        //utility
    }

    public static Connection connection(String propertiesFile) {

        Properties properties = new Properties();
        try {
            properties.load(Utils.getBufferedReader(propertiesFile));
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalArgumentException("BMSConnection: connection: problem reading properties file: " + propertiesFile);
        }

        String host = properties.getProperty("host");
        String user = properties.getProperty("user");
        String password = properties.getProperty("password");
        String dbName = properties.getProperty("DB");

        return connection(host, user, password, dbName);

    }

    public static Connection connection(String host, String user, String password, String dbName) {

        Connection connection = null;
        String url = "jdbc:mysql://" + host + "/" + dbName;
        try {
            Class.forName("com.mysql.jdbc.Driver");
            connection = DriverManager.getConnection(url, user, password);
        } catch (ClassNotFoundException e) {
            myLogger.error(e.getMessage(), e);
            throw new IllegalStateException("BMSConnection: connection: org.mysql.Driver can't be found");
        } catch (SQLException e) {
            myLogger.error(e.getMessage(), e);
            throw new IllegalStateException("BMSConnection: connection: problem connecting to database: " + e.getMessage());
        }
        myLogger.info("Connected to database:  " + url + "\n");
        return connection;

    }

    public static void taxaList(Connection connection, Map<String, Taxon.Builder> gids) {

        if (connection == null) {
            throw new IllegalArgumentException("BMSConnection: taxaList: Must specify database connection.");
        }

        //
        // germplasm.external_code (GOBII Postgres) and names.gid (BMS) is GID
        //
        // select gid, nval from names
        // where gid in (1, 2, 3);
        //
        StringBuilder builder = new StringBuilder();
        builder.append("select gid, nval from names ");
        builder.append("where gid in (");
        boolean first = true;
        for (Map.Entry<String, Taxon.Builder> current : gids.entrySet()) {
            if (first) {
                first = false;
            } else {
                builder.append(", ");
            }
            builder.append(current.getKey());
        }
        builder.append(");");

        String query = builder.toString();
        myLogger.info("taxaList: query statement: " + query);

        try (ResultSet rs = connection.createStatement().executeQuery(query)) {
            while (rs.next()) {
                String gid = rs.getString("gid");
                if (gids.get(gid) != null) {
                    throw new IllegalStateException("BMSConnection: gid is duplicated: " + gid);
                }
                gids.put(gid, new Taxon.Builder(rs.getString("nval")));
            }

        } catch (Exception se) {
            myLogger.debug(se.getMessage(), se);
            throw new IllegalStateException("BMSConnection: taxaList: Problem querying the database: " + se.getMessage());
        }

    }

    public static void main(String[] args) {
        
        LoggingUtils.setupDebugLogging();

        String[] gidsStrings = new String[]{"42000018",
            "42000250",
            "42000118",
            "42000139",
            "42000179",
            "42000038",
            "42000215",
            "42000227",
            "42000273",
            "42000186",
            "42000057"};

        Connection connection = connection("/home/tmc46/bms_panzea.txt");
        Map<String, Taxon.Builder> gids = new LinkedHashMap<>();
        for (String current : gidsStrings) {
            gids.put(current, null);
        }
        taxaList(connection, gids);
        
        TaxaListBuilder builder = new TaxaListBuilder();
        for (String key : gids.keySet()) {
            builder.add(gids.get(key).build());
        }
        
        TaxaList taxa = builder.build();
        
        System.out.println("number of taxa: " + taxa.numberOfTaxa());
        
        for (Taxon current : taxa) {
            System.out.println(current.getName());
        }

    }

}
