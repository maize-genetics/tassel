package net.maizegenetics.plugindef;

import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.util.Enumeration;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Properties;

/**
 * This class is for storing a cache of parameters retrieved from a {@link java.util.Properties} file that is available
 * to all plugins on the command line. This is used when run_pipeline.pl -configParameters config.txt... is used.
 *
 * config.txt example entries as follow... First five are using by GetDBConnectionPlugin -config config.txt but that
 * should be refactored to individual PluginParameters. Then entry example would be GetDBConnectionPlugin.DBType=sqlite
 * KinshipPlugin.method=TYPE is the standard.  PLUGIN_NAME.COMMAND_LINE_NAME=VALUE
 *
 * host=localHost user=sqlite password=sqlite DB=/tempFileDir/outputDir/phgTestDB.db DBtype=sqlite
 * KinshipPlugin.method=Normalized_IBS
 *
 * @author Terry Casstevens Created January 05, 2018
 */
public class ParameterCache {

    private static final Logger myLogger = LogManager.getLogger(ParameterCache.class);

    private static Properties CACHE = null;

    private ParameterCache() {
    }

    /**
     * Loads the parameter cache with the values in the given {@link java.util.Properties} file.
     *
     * @param filename file name
     */
    public static void load(String filename) {

        if (filename == null || filename.trim().isEmpty()) {
            CACHE = null;
            return;
        }

        myLogger.info("load: loading parameter cache with: " + filename);

        CACHE = new Properties();
        try (BufferedReader reader = Utils.getBufferedReader(filename)) {
            CACHE.load(reader);
            CACHE.put("configFile", filename);
        } catch (Exception e) {
            CACHE = null;
            myLogger.debug(e.getMessage(), e);
            throw new IllegalArgumentException("ParameterCache: load: problem reading properties file: " + filename + "\n" + e.getMessage());
        }

        for (String key : CACHE.stringPropertyNames()) {
            if (key.toLowerCase().contains("password") || key.toLowerCase().contains("passwd")) {
                myLogger.info("ParameterCache: key: " + key + " value: ********");
            } else {
                myLogger.info("ParameterCache: key: " + key + " value: " + CACHE.getProperty(key));
            }
        }

        verify(filename);

    }

    /**
     * Verifies specified parameters exist on given plugin
     */
    private static void verify(String filename) {

        for (Map.Entry<Object, Object> current : CACHE.entrySet()) {
            int dot = ((String) current.getKey()).indexOf('.');
            if (dot != -1) {
                String possibleClassName = ((String) current.getKey()).substring(0, dot);
                Plugin plugin = null;
                List<String> matches = Utils.getFullyQualifiedClassNames(possibleClassName);
                for (String match : matches) {
                    plugin = Plugin.getPluginInstance(match, null);
                    if (plugin != null) {
                        break;
                    }
                }
                if (plugin == null) {
                    throw new IllegalArgumentException("ParameterCache: verify: unknown plugin: " + possibleClassName + " in file: " + filename);
                } else {
                    String parameterName = ((String) current.getKey()).substring(dot + 1);
                    if (!((AbstractPlugin) plugin).isPluginParameter(parameterName)) {
                        throw new IllegalArgumentException("ParameterCache: verify: unknown parameter: " + parameterName + " for plugin: " + possibleClassName + " in file: " + filename);
                    }
                }
            }

        }

    }

    public static Optional<String> value(Plugin plugin, String parameter) {
        return value(plugin.getClass().getName(), parameter);
    }

    /**
     * Returns the value if any for the given plugin and parameter.  Value for PluginClassName.parameter will be
     * returned if it exists.  If not, value for parameter will be returned.  Otherwise, an empty optional.
     *
     * @param plugin plugin
     * @param parameter parameter
     *
     * @return value if exists
     */
    public static Optional<String> value(String plugin, String parameter) {

        if (CACHE == null) {
            return Optional.empty();
        }

        String value = CACHE.getProperty(Utils.getBasename(plugin) + "." + parameter);
        if (value != null) {
            return Optional.of(value);
        }

        value = CACHE.getProperty(parameter);
        return Optional.ofNullable(value);

    }

    /**
     * Returns the value if any for the given key.
     *
     * @param key key
     *
     * @return value if exists
     */
    public static Optional<String> value(String key) {
        if (CACHE == null) {
            return Optional.empty();
        }
        return Optional.ofNullable(CACHE.getProperty(key));
    }

    /**
     * Returns true if this cache has been loaded with values.
     *
     * @return true if this cache has been loaded with values.
     */
    public static boolean hasValues() {
        return CACHE != null;
    }

    public static Enumeration<?> keys() {
        if (CACHE == null) {
            return null;
        }
        return CACHE.propertyNames();
    }

}
