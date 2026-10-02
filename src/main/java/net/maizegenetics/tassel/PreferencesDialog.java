package net.maizegenetics.tassel;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.ParameterCache;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.prefs.TasselPrefs;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import net.maizegenetics.plugindef.GeneratedGuiBoilerplate;

/**
 * @author Terry Casstevens
 */
public class PreferencesDialog extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(PreferencesDialog.class);

    private PluginParameter<Boolean> myRetainRareAlleles = new PluginParameter.Builder<>("retainRareAlleles", TasselPrefs.ALIGNMENT_RETAIN_RARE_ALLELES_DEFAULT, Boolean.class)
            .description("True if rare alleles should be retained.  This has no effect on Nucleotide Data as all alleles will be retained regardless.")
            .build();

    private PluginParameter<Boolean> mySendLogToConsole = new PluginParameter.Builder<>("sendLogToConsole", TasselPrefs.TASSEL_LOG_SEND_TO_CONSOLE_DEFAULT, Boolean.class)
            .description("Flag whether to send logging to the console.")
            .build();

    private PluginParameter<Boolean> myDarkTheme = new PluginParameter.Builder<>("darkTheme", TasselPrefs.TASSEL_DARK_THEME_DEFAULT, Boolean.class)
            .description("Use the dark color theme (unchecked = light theme).")
            .build();

    private PluginParameter<String> myConfigFile = new PluginParameter.Builder<>("configFile", TasselPrefs.TASSEL_CONFIG_FILE_DEFAULT, String.class)
            .description("Global configuration file")
            .required(false)
            .inFile()
            .build();

    private PluginParameter<LocaleWrapper> myLocale = new PluginParameter.Builder<>("locale", null, LocaleWrapper.class)
            .description("Default locale for TASSEL")
            .required(false)
            .objectListSingleSelect()
            .build();

    public PreferencesDialog(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
        // Don't Load global parameters here
        // It will cause GUI to fail startup if error occurs
        // ParameterCache.load(TasselPrefs.getConfigFile());
        Locale locale = TasselPrefs.getLocale();
        Locale.setDefault(locale);
        myLogger.info("TasselPrefs: default locale set: " + locale.getDisplayName());
    }

    @Override
    protected void preProcessParameters(DataSet input) {

        setParameter(myRetainRareAlleles, TasselPrefs.getAlignmentRetainRareAlleles());
        setParameter(mySendLogToConsole, TasselPrefs.getLogSendToConsole());
        setParameter(myDarkTheme, TasselPrefs.getDarkTheme());
        setParameter(myConfigFile, TasselPrefs.getConfigFile());

        List<LocaleWrapper> temp = new ArrayList<>();
        temp.add(new LocaleWrapper(TasselPrefs.getLocale()));
        temp.add(new LocaleWrapper(Locale.US));
        for (Locale current : Locale.getAvailableLocales()) {
            if (current.toString() != null && !current.toString().isEmpty())
                temp.add(new LocaleWrapper(current));
        }
        myLocale = new PluginParameter<>(myLocale, temp);

    }

    @Override
    public DataSet processData(DataSet input) {

        TasselPrefs.putAlignmentRetainRareAlleles(retainRareAlleles());

        TasselPrefs.putLogSendToConsole(sendLogToConsole());
        TasselLogging.updateLoggingLocation();

        if (TasselPrefs.getDarkTheme() != darkTheme()) {
            TasselPrefs.putDarkTheme(darkTheme());
            TASSELMainApp.setupLookAndFeel();
            com.formdev.flatlaf.FlatLaf.updateUI();
        }

        TasselPrefs.putConfigFile(configFile());
        ParameterCache.load(TasselPrefs.getConfigFile());

        TasselPrefs.putLocale(locale().myLocale);
        Locale.setDefault(locale().myLocale);

        ((TASSELMainFrame) getParentFrame()).updatePluginsWithGlobalConfigParameters();

        return null;

    }

    /**
     * Retain Rare Alleles
     *
     * @return Retain Rare Alleles
     */
    @GeneratedGuiBoilerplate
    public Boolean retainRareAlleles() {
        return myRetainRareAlleles.value();
    }

    /**
     * Set Retain Rare Alleles. Retain Rare Alleles
     *
     * @param value Retain Rare Alleles
     *
     * @return this plugin
     */
    @GeneratedGuiBoilerplate
    public PreferencesDialog retainRareAlleles(Boolean value) {
        myRetainRareAlleles = new PluginParameter<>(myRetainRareAlleles, value);
        return this;
    }

    /**
     * Flag whether to send logging to the console.
     *
     * @return Send Log To Console
     */
    @GeneratedGuiBoilerplate
    public Boolean sendLogToConsole() {
        return mySendLogToConsole.value();
    }

    /**
     * Set Send Log To Console. Flag whether to send logging to the console.
     *
     * @param value Send Log To Console
     *
     * @return this plugin
     */
    @GeneratedGuiBoilerplate
    public PreferencesDialog sendLogToConsole(Boolean value) {
        mySendLogToConsole = new PluginParameter<>(mySendLogToConsole, value);
        return this;
    }

    /**
     * Use the dark color theme.
     *
     * @return Dark Theme
     */
    public Boolean darkTheme() {
        return myDarkTheme.value();
    }

    /**
     * Set Dark Theme. Use the dark color theme (unchecked = light theme).
     *
     * @param value Dark Theme
     *
     * @return this plugin
     */
    public PreferencesDialog darkTheme(Boolean value) {
        myDarkTheme = new PluginParameter<>(myDarkTheme, value);
        return this;
    }

    /**
     * Global configuration file
     *
     * @return Config File
     */
    @GeneratedGuiBoilerplate
    public String configFile() {
        return myConfigFile.value();
    }

    /**
     * Set Config File. Global configuration file
     *
     * @param value Config File
     *
     * @return this plugin
     */
    @GeneratedGuiBoilerplate
    public PreferencesDialog configFile(String value) {
        myConfigFile = new PluginParameter<>(myConfigFile, value);
        return this;
    }

    /**
     * Default locale for TASSEL
     *
     * @return Locale
     */
    @GeneratedGuiBoilerplate
    public LocaleWrapper locale() {
        return myLocale.value();
    }

    /**
     * Set Locale. Default locale for TASSEL
     *
     * @param value Locale
     *
     * @return this plugin
     */
    @GeneratedGuiBoilerplate
    public PreferencesDialog locale(LocaleWrapper value) {
        myLocale = new PluginParameter<>(myLocale, value);
        return this;
    }

    @Override
    @GeneratedGuiBoilerplate
    public ImageIcon getIcon() {
        URL imageURL = TasselLogging.class.getResource("/net/maizegenetics/analysis/images/preferences.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    @GeneratedGuiBoilerplate
    public String getButtonName() {
        return "Preferences";
    }

    @Override
    @GeneratedGuiBoilerplate
    public String getToolTipText() {
        return "Preferences";
    }

    private class LocaleWrapper {
        public final Locale myLocale;

        private LocaleWrapper(Locale locale) {
            myLocale = locale;
        }

        @Override
        public String toString() {
            return myLocale.getDisplayName();
        }
    }

}
