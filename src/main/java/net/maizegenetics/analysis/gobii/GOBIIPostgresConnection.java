/*
 *  GOBIIPostgresConnection
 * 
 *  Created on Jun 15, 2016
 */
package net.maizegenetics.analysis.gobii;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Properties;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class GOBIIPostgresConnection {

    private static final Logger myLogger = LogManager.getLogger(GOBIIPostgresConnection.class);

    private GOBIIPostgresConnection() {
        //utility
    }

    /**
     * Creates a Postgres database connection given a properties file
     *
     * @param propertiesFile properties file
     *
     * @return Postgres database connection
     */
    public static Connection connection(String propertiesFile) {

        Properties properties = new Properties();
        try {
            properties.load(Utils.getBufferedReader(propertiesFile));
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalArgumentException("GOBIIPostgresConnection: connection: problem reading properties file: " + propertiesFile);
        }

        String host = properties.getProperty("host");
        String user = properties.getProperty("user");
        String password = properties.getProperty("password");
        String dbName = properties.getProperty("DB");

        return connection(host, user, password, dbName);

    }

    /**
     * Creates a Postgres database connection.
     *
     * @param host hostname
     * @param user user id
     * @param password password
     * @param dbName database name
     *
     * @return Postgres database connection
     */
    public static Connection connection(String host, String user, String password, String dbName) {

        Connection connection = null;
        String url = "jdbc:postgresql://" + host + "/" + dbName;
        try {
            Class.forName("org.postgresql.Driver");
            connection = DriverManager.getConnection(url, user, password);
        } catch (ClassNotFoundException e) {
            myLogger.error(e.getMessage(), e);
            throw new IllegalStateException("GOBIIPostgresConnection: connection: org.postgresql.Driver can't be found");
        } catch (SQLException e) {
            myLogger.error(e.getMessage(), e);
            throw new IllegalStateException("GOBIIPostgresConnection: connection: problem connecting to database: " + e.getMessage());
        }
        myLogger.info("Connected to database:  " + url + "\n");
        return connection;

    }

    public static TaxaList taxaList(Connection postgres, String datasetName, Connection bms) {

        if (postgres == null) {
            throw new IllegalArgumentException("GOBIIPostgresConnection: taxaList: Must specify database connection.");
        }

        //
        // germplasm.external_code is GID that maps to BMS
        //
        // select distinct(germplasm.external_code) from dataset, dataset_dnarun, dnarun,
        // dnasample, germplasm
        // where dataset.name='maize282_raw_AGPv2'
        // and dataset_dnarun.dataset_id = dataset.dataset_id
        // and dnarun.dnarun_id = dataset_dnarun.dnarun_id
        // and dnarun.dnasample_id = dnasample.dnasample_id
        // and dnasample.germplasm_id = germplasm.germplasm_id
        // order by dataset_dnarun.dnarun_idx;
        //
        StringBuilder builder = new StringBuilder();
        builder.append("select germplasm.external_code from dataset, dataset_dnarun, dnarun, dnasample, germplasm ");
        builder.append("where dataset.name='");
        builder.append(datasetName);
        builder.append("'");
        builder.append(" and dataset_dnarun.dataset_id = dataset.dataset_id");
        builder.append(" and dnarun.dnarun_id = dataset_dnarun.dnarun_id");
        builder.append(" and dnarun.dnasample_id = dnasample.dnasample_id");
        builder.append(" and dnasample.germplasm_id = germplasm.germplasm_id");
        builder.append(" order by dataset_dnarun.dnarun_idx;");

        String query = builder.toString();
        myLogger.info("taxaList: query statement: " + query);

        Map<String, Taxon.Builder> gids = new LinkedHashMap<>();
        try (ResultSet rs = postgres.createStatement().executeQuery(query)) {

            while (rs.next()) {
                String externalCode = rs.getString("external_code");
                gids.put(externalCode, null);
            }

        } catch (Exception se) {
            myLogger.debug(se.getMessage(), se);
            throw new IllegalStateException("GOBIIPostgresConnection: taxaList: Problem querying the database: " + se.getMessage());
        }

        BMSConnection.taxaList(bms, gids);
        TaxaListBuilder taxa = new TaxaListBuilder();
        for (Map.Entry<String, Taxon.Builder> current : gids.entrySet()) {
            taxa.add(current.getValue().build());
        }
        return taxa.build();

    }

    public static PositionList positionList(Connection connection, String datasetName) {

        if (connection == null) {
            throw new IllegalArgumentException("GOBIIPostgresConnection: positionList: Must specify database connection.");
        }

        //
        // select marker.name, marker_linkage_group.start, marker_linkage_group.stop, linkage_group.name
        // from dataset, dataset_marker, marker, marker_linkage_group, linkage_group
        // where dataset.name='maize282_raw_AGPv2'
        // and dataset.dataset_id=dataset_marker.dataset_id
        // and dataset_marker.marker_id=marker.marker_id
        // and marker.marker_id=marker_linkage_group.marker_id
        // and marker_linkage_group.linkage_group_id=linkage_group.linkage_group_id
        //  order by dataset_marker.marker_idx;
        //
        StringBuilder builder = new StringBuilder();
        builder.append("select marker.name, marker_linkage_group.start, marker_linkage_group.stop, linkage_group.name ");
        builder.append("from dataset, dataset_marker, marker, marker_linkage_group, linkage_group ");
        builder.append("where dataset.name='");
        builder.append(datasetName);
        builder.append("'");
        builder.append(" and dataset.dataset_id=dataset_marker.dataset_id");
        builder.append(" and dataset_marker.marker_id=marker.marker_id");
        builder.append(" and marker.marker_id=marker_linkage_group.marker_id");
        builder.append(" and marker_linkage_group.linkage_group_id=linkage_group.linkage_group_id");
        builder.append(" order by dataset_marker.marker_idx;");

        String query = builder.toString();
        myLogger.info("positionList: query statement: " + query);

        Map<String, Chromosome> chromosomes = new HashMap<>();
        try (ResultSet rs = connection.createStatement().executeQuery(query)) {
            PositionListBuilder positions = new PositionListBuilder();
            while (rs.next()) {

                // marker_linkage_group.start
                int start = rs.getInt("start");
                // marker_linkage_group.stop
                int end = rs.getInt("stop");
                if (start != end) {
                    throw new IllegalArgumentException("GOBIIPostgresConnection: positionList: start position: " + start + " and end position: " + end + " should be the same.");
                }

                // marker.name
                String snpName = rs.getString(1);

                // linkage_group.name
                String chrStr = rs.getString(4);
                Chromosome chr = chromosomes.get(chrStr);
                if (chr == null) {
                    chr = new Chromosome(chrStr);
                    chromosomes.put(chrStr, chr);
                }

                GeneralPosition.Builder current = new GeneralPosition.Builder(chr, start);
                current.snpName(snpName);
                positions.add(current.build());

            }
            return positions.build();
        } catch (Exception se) {
            myLogger.debug(se.getMessage(), se);
            throw new IllegalStateException("GOBIIPostgresConnection: taxaList: Problem querying the database: " + se.getMessage());
        }

    }

    public static String hdf5Filename(Connection connection, String datasetName) {

        if (connection == null) {
            throw new IllegalArgumentException("GOBIIPostgresConnection: hdf5Filename: Must specify database connection.");
        }

        //
        // select data_file from dataset where dataset.name='name';
        //
        StringBuilder builder = new StringBuilder();
        builder.append("select data_file ");
        builder.append("from dataset ");
        builder.append("where dataset.name='");
        builder.append(datasetName);
        builder.append("';");

        String query = builder.toString();
        myLogger.info("hdf5Filename: query statement: " + query);

        try (ResultSet rs = connection.createStatement().executeQuery(query)) {
            if (!rs.next()) {
                throw new IllegalStateException("GOBIIPostgresConnection: hdf5Filename: dataset name: " + datasetName + " doesn't map to any hdf5 filename.");
            }
            String result = rs.getString("data_file");
            if (rs.next()) {
                throw new IllegalStateException("GOBIIPostgresConnection: hdf5Filename: dataset name: " + datasetName + " maps to more than one hdf5 filename.");
            }
            if ((result == null) || (result.length() == 0)) {
                throw new IllegalStateException("GOBIIPostgresConnection: hdf5Filename: dataset name: " + datasetName + " doesn't map to any hdf5 filename.");
            }
            myLogger.info("Dataset Name: " + datasetName + " HDF5 Filename: " + result);
            return result;
        } catch (Exception se) {
            myLogger.debug(se.getMessage(), se);
            throw new IllegalStateException("GOBIIPostgresConnection: hdf5Filename: Problem querying the database: " + se.getMessage());
        }

    }

    public static void printAvailableDatasets(Connection connection) {

        if (connection == null) {
            throw new IllegalArgumentException("GOBIIPostgresConnection: printAvailableDatasets: Must specify database connection.");
        }

        //
        // select name from dataset;
        //
        StringBuilder builder = new StringBuilder();
        builder.append("select name from dataset;");

        String query = builder.toString();
        myLogger.info("printAvailableDatasets: query statement: " + query);

        myLogger.info("Avaliable Datasets...");
        try (ResultSet rs = connection.createStatement().executeQuery(query)) {
            while (rs.next()) {
                System.out.println(rs.getString("name"));
            }
        } catch (Exception se) {
            myLogger.debug(se.getMessage(), se);
            throw new IllegalStateException("GOBIIPostgresConnection: printAvailableDatasets: Problem querying the database: " + se.getMessage());
        }

    }

}
