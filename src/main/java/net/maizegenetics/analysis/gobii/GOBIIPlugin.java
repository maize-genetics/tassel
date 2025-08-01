/*
 *  GOBIIPlugin
 * 
 *  Created on July 25, 2016
 */
package net.maizegenetics.analysis.gobii;

import java.awt.Frame;
import java.sql.Connection;
import javax.swing.ImageIcon;
import net.maizegenetics.analysis.data.GenotypeSummaryPlugin;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.GOBIIGenotypeCallTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.taxa.TaxaList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class GOBIIPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(GOBIIPlugin.class);

    private PluginParameter<String> myDataset = new PluginParameter.Builder<>("dataset", null, String.class)
            .required(true)
            .description("Dataset Name")
            .guiName("Dataset Name")
            .build();

    private final PluginParameter<String> myDatabaseLabel = PluginParameter.getLabelInstance("GOBII Postgres Database Properties");

    private PluginParameter<String> myDBName = new PluginParameter.Builder<>("db", null, String.class)
            .required(true)
            .guiName("Database Name")
            .description("Database Name")
            .build();

    private PluginParameter<String> myUser = new PluginParameter.Builder<>("user", null, String.class)
            .required(true)
            .description("User Name")
            .build();

    private PluginParameter<String> myPassword = new PluginParameter.Builder<>("password", "", String.class)
            .description("Password")
            .password()
            .build();

    private final PluginParameter<String> myBMSLabel = PluginParameter.getLabelInstance("BMS Database Properties");

    private PluginParameter<String> myBMSDBName = new PluginParameter.Builder<>("bmsdb", null, String.class)
            .required(true)
            .guiName("BMS Database Name")
            .description("BMS Database Name")
            .build();

    private PluginParameter<String> myBMSHost = new PluginParameter.Builder<>("bmshost", "localhost", String.class)
            .guiName("BMS Host Name")
            .description("BMS Host Name")
            .build();

    private PluginParameter<String> myBMSUser = new PluginParameter.Builder<>("bmsuser", null, String.class)
            .required(true)
            .guiName("BMS User")
            .description("BMS User Name")
            .build();

    private PluginParameter<String> myBMSPassword = new PluginParameter.Builder<>("bmspassword", "", String.class)
            .guiName("BMS Password")
            .description("Password")
            .password()
            .build();

    public GOBIIPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        if ((dBName() == null) || (dBName().length() == 0)) {
            dBName(TasselPrefs.getGOBIIDB());
        }
        if ((user() == null) || (user().length() == 0)) {
            user(TasselPrefs.getGOBIIUser());
        }
    }

    @Override
    protected void postProcessParameters() {
        TasselPrefs.putGOBIIDB(dBName());
        TasselPrefs.putGOBIIUser(user());
        TasselPrefs.putBMSHost(bmsHost());
        TasselPrefs.putBMSDB(bmsDBName());
        TasselPrefs.putBMSUser(bmsUser());
    }

    @Override
    public DataSet processData(DataSet input) {

        Connection postgres = GOBIIPostgresConnection.connection("localhost", user(), myPassword.value(), dBName());

        Connection bms = BMSConnection.connection(bmsHost(), bmsUser(), myBMSPassword.value(), bmsDBName());
        String hdf5Filename = null;
        try {
            hdf5Filename = GOBIIPostgresConnection.hdf5Filename(postgres, dataset());
        } catch (Exception e) {
            GOBIIPostgresConnection.printAvailableDatasets(postgres);
            throw e;
        }

        TaxaList taxa = GOBIIPostgresConnection.taxaList(postgres, dataset(), bms);
        myLogger.info("Number of Taxa: " + taxa.numberOfTaxa());

        PositionList positions = GOBIIPostgresConnection.positionList(postgres, dataset());
        myLogger.info("Number of Positions: " + positions.numberOfSites());

        GOBIIGenotypeCallTable genotypes = GOBIIGenotypeCallTable.getInstance(taxa.numberOfTaxa(), positions.numberOfSites(), false, hdf5Filename);

        GenotypeTable genotypeTable = GenotypeTableBuilder.getInstance(genotypes, positions, taxa);
        DataSet result = new DataSet(new Datum("gobii:" + dataset(), genotypeTable, null), this);
        GenotypeSummaryPlugin.printSimpleSummary(result);
        return result;

    }

    /**
     * Dataset Name
     *
     * @return Dataset Name
     */
    public String dataset() {
        return myDataset.value();
    }

    /**
     * Set Dataset Name. Dataset Name
     *
     * @param value Dataset Name
     *
     * @return this plugin
     */
    public GOBIIPlugin dataset(String value) {
        myDataset = new PluginParameter<>(myDataset, value);
        return this;
    }

    /**
     * Database Name
     *
     * @return Database Name
     */
    public String dBName() {
        return myDBName.value();
    }

    /**
     * Set Database Name. Database Name
     *
     * @param value Database Name
     *
     * @return this plugin
     */
    public GOBIIPlugin dBName(String value) {
        myDBName = new PluginParameter<>(myDBName, value);
        return this;
    }

    /**
     * User Name
     *
     * @return User
     */
    public String user() {
        return myUser.value();
    }

    /**
     * Set User. User Name
     *
     * @param value User
     *
     * @return this plugin
     */
    public GOBIIPlugin user(String value) {
        myUser = new PluginParameter<>(myUser, value);
        return this;
    }

    /**
     * Set Password. Password
     *
     * @param value Password
     *
     * @return this plugin
     */
    public GOBIIPlugin password(String value) {
        myPassword = new PluginParameter<>(myPassword, value);
        return this;
    }

    /**
     * Database Name
     *
     * @return Database Name
     */
    public String bmsDBName() {
        return myBMSDBName.value();
    }

    /**
     * Set Database Name. Database Name
     *
     * @param value Database Name
     *
     * @return this plugin
     */
    public GOBIIPlugin bmsDBName(String value) {
        myBMSDBName = new PluginParameter<>(myBMSDBName, value);
        return this;
    }

    /**
     * BMS Host Name
     *
     * @return BMS Host Name
     */
    public String bmsHost() {
        return myBMSHost.value();
    }

    /**
     * Set BMS Host Name. BMS Host Name
     *
     * @param value BMS Host Name
     *
     * @return this plugin
     */
    public GOBIIPlugin bmsHost(String value) {
        myBMSHost = new PluginParameter<>(myBMSHost, value);
        return this;
    }

    /**
     * User Name
     *
     * @return BMS User
     */
    public String bmsUser() {
        return myBMSUser.value();
    }

    /**
     * Set BMS User. User Name
     *
     * @param value BMS user
     *
     * @return this plugin
     */
    public GOBIIPlugin bmsUser(String value) {
        myBMSUser = new PluginParameter<>(myBMSUser, value);
        return this;
    }

    /**
     * Set BMS Password. Password
     *
     * @param value BMS Password
     *
     * @return this plugin
     */
    public GOBIIPlugin bmsPassword(String value) {
        myBMSPassword = new PluginParameter<>(myBMSPassword, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "GOBII Plugin";
    }

    @Override
    public String getToolTipText() {
        return "GOBII Plugin";
    }

}
