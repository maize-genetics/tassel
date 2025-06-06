/*
 * AbstractPlugin.java
 *
 * Created on December 22, 2006, 5:03 PM
 *
 */
package net.maizegenetics.plugindef;

import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.io.JSONUtils;
import net.maizegenetics.dna.snp.io.ReadBedfile;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.gui.SelectFromAvailableDialog;
import net.maizegenetics.gui.SiteNamesAvailableListModel;
import net.maizegenetics.gui.TaxaAvailableListModel;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.ReadDistanceMatrix;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.File;
import java.lang.reflect.Field;
import java.math.BigDecimal;
import java.net.URI;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.List;
import java.util.regex.Pattern;

/**
 * @author Terry Casstevens
 */
abstract public class AbstractPlugin implements Plugin {

    private static final Logger myLogger = LogManager.getLogger(AbstractPlugin.class);

    public static final String DEFAULT_CITATION = "Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES. (2007) TASSEL: Software for association mapping of complex traits in diverse samples. Bioinformatics 23:2633-2635.";

    public static final String POSITION_LIST_NONE = "None";
    public static final String TAXA_LIST_NONE = "None";

    private final List<PluginListener> myListeners = new ArrayList<>();
    private final List<Plugin> myInputs = new ArrayList<>();
    private DataSet myCurrentInputData = null;
    private final Frame myParentFrame;
    private final boolean myIsInteractive;
    private boolean myTrace = false;
    private boolean myThreaded = false;
    protected boolean myWasCancelled = false;

    /**
     * Creates a new instance of AbstractPlugin
     */
    public AbstractPlugin() {
        this(null, true);
    }

    /**
     * Creates a new instance of AbstractPlugin
     */
    public AbstractPlugin(Frame parentFrame, boolean isInteractive) {
        myParentFrame = parentFrame;
        myIsInteractive = isInteractive;
    }

    @Override
    public DataSet performFunction(DataSet input) {

        LocalDateTime time = LocalDateTime.now();
        String timeStr = time.format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
        myLogger.info("Starting " + getClass().getName() + ": time: " + timeStr);

        myCurrentInputData = input;

        try {

            preProcessParameters(input);

            if (isInteractive()) {
                if (!setParametersViaGUI()) {
                    return null;
                }
            }

            checkRequiredParameters();
            postProcessParameters();
            printParameterValues();
            checkParameters();

            DataSet output = processData(input);
            time = LocalDateTime.now();
            timeStr = time.format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
            myLogger.info("Finished " + getClass().getName() + ": time: " + timeStr);
            fireDataSetReturned(new PluginEvent(output, getClass()));
            return output;

        } catch (Exception e) {

            if (isInteractive()) {
                myLogger.debug(e.getMessage(), e);
                DialogUtils.showError(e.getMessage() + "\n", getParentFrame());
            } else {
                myLogger.debug(e.getMessage(), e);
                printUsage();
                myLogger.error(e.getMessage());
                if (!getListeners().isEmpty()) {
                    System.exit(1);
                }
            }
            return null;

        } finally {
            fireProgress(100);
        }

    }

    protected void preProcessParameters(DataSet input) {
        // do nothing
    }

    protected void postProcessParameters() {
        // do nothing
    }

    @Override
    public DataSet processData(DataSet input) {
        throw new UnsupportedOperationException();
    }

    protected List<Field> getParameterFields() {

        List<Field> result = new ArrayList<>();
        Field[] fields = getClass().getDeclaredFields();
        for (Field current : fields) {
            if (current.getType().isAssignableFrom(PluginParameter.class)) {
                current.setAccessible(true);
                result.add(current);
            }
        }

        return result;
    }

    private List<PluginParameter<?>> getParameterInstances() {

        List<PluginParameter<?>> result = new ArrayList<>();
        Field[] fields = getClass().getDeclaredFields();
        for (Field current : fields) {
            if (current.getType().isAssignableFrom(PluginParameter.class)) {
                current.setAccessible(true);
                try {
                    PluginParameter<?> parameter = (PluginParameter) current.get(this);
                    result.add(parameter);
                } catch (Exception e) {
                    e.printStackTrace();
                    throw new IllegalArgumentException("AbstractPlugin: getParameterInstances: problem getting parameter instances");
                }

            }
        }

        return result;
    }

    private Field getParameterField(String key) {

        Field[] fields = getClass().getDeclaredFields();
        for (Field current : fields) {
            if (current.getType().isAssignableFrom(PluginParameter.class)) {
                try {
                    current.setAccessible(true);
                    PluginParameter<?> parameter = (PluginParameter) current.get(this);
                    if (parameter.cmdLineName().equals(key)) {
                        return current;
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                    throw new IllegalArgumentException("AbstractPlugin: getParameterField: problem with key: " + key);
                }

            }
        }

        throw new IllegalArgumentException("AbstractPlugin: getParameterField: unknown key: " + key);
    }

    private PluginParameter<?> getParameterInstance(String key) {
        try {
            Field field = getParameterField(key);
            return (PluginParameter) field.get(this);
        } catch (Exception e) {
            return null;
        }
    }

    public boolean isPluginParameter(String key) {
        if (getParameterInstance(key) == null) {
            return false;
        } else {
            return true;
        }
    }

    public static <T> T convert(String input, Class<T> outputClass) {
        try {
            if ((input == null) || (input.length() == 0)) {
                return null;
            } else if (outputClass.isEnum()) {
                return (T) Enum.valueOf((Class<Enum>) outputClass, input);
            } else if (outputClass.isAssignableFrom(String.class)) {
                return (T) input;
            } else if (outputClass.isAssignableFrom(Integer.class)) {
                char groupingSeparator = DecimalFormatSymbols.getInstance(Locale.getDefault()).getGroupingSeparator();
                input = input.replace(String.valueOf(groupingSeparator), "");
                char decimalSeparator = DecimalFormatSymbols.getInstance(Locale.getDefault()).getDecimalSeparator();
                input = input.replace(String.valueOf(decimalSeparator), ".");
                return (T) new Integer(new BigDecimal(input).intValueExact());
            } else if (outputClass.isAssignableFrom(Double.class)) {
                char groupingSeparator = DecimalFormatSymbols.getInstance(Locale.getDefault()).getGroupingSeparator();
                input = input.replace(String.valueOf(groupingSeparator), "");
                char decimalSeparator = DecimalFormatSymbols.getInstance(Locale.getDefault()).getDecimalSeparator();
                input = input.replace(String.valueOf(decimalSeparator), ".");
                return (T) new Double(new BigDecimal(input).doubleValue());
            } else if (outputClass.isAssignableFrom(Float.class)) {
                char groupingSeparator = DecimalFormatSymbols.getInstance(Locale.getDefault()).getGroupingSeparator();
                input = input.replace(String.valueOf(groupingSeparator), "");
                char decimalSeparator = DecimalFormatSymbols.getInstance(Locale.getDefault()).getDecimalSeparator();
                input = input.replace(String.valueOf(decimalSeparator), ".");
                return (T) new Double(new BigDecimal(input).floatValue());
            } else if (outputClass.isAssignableFrom(List.class)) {
                return (T) getListFromString(input);
            } else if (outputClass.isAssignableFrom(SortedSet.class)) {
                return (T) getSortedSet(input);
            } else if (outputClass.isAssignableFrom(PositionList.class)) {
                String test = input.trim().substring(Math.max(0, input.length() - 8)).toLowerCase();
                if ((test.endsWith(".bed")) || (test.endsWith(".bed.gz"))) {
                    return (T) ReadBedfile.getPositionList(input);
                } else if (test.endsWith(".json") || test.endsWith(".json.gz")) {
                    return (T) JSONUtils.importPositionListFromJSON(input);
                } else {
                    try {
                        GenotypeTable temp = ImportUtils.read(input);
                        return (T) temp.positions();
                    } catch (Exception e) {
                        myLogger.debug(e.getMessage(), e);
                        throw new IllegalArgumentException("AbstractPlugin: convert: don't know who to covert: " + input + " to postion list");
                    }
                }
            } else if (outputClass.isAssignableFrom(TaxaList.class)) {
                String test = input.trim().substring(Math.max(0, input.length() - 8)).toLowerCase();
                if (test.endsWith(".json") || test.endsWith(".json.gz")) {
                    return (T) JSONUtils.importTaxaListFromJSON(input);
                } else if (test.endsWith(".txt")) {
                    TaxaListBuilder builder = new TaxaListBuilder();
                    try (BufferedReader br = Utils.getBufferedReader(input)) {
                        String line = br.readLine();
                        Pattern sep = Pattern.compile("[\\s,]+");

                        while (line != null) {
                            line = line.trim();
                            String[] parsedline = sep.split(line);
                            for (int i = 0; i < parsedline.length; i++) {
                                if ((parsedline[i] != null) || (parsedline[i].length() != 0)) {
                                    builder.add(parsedline[i]);
                                }
                            }
                            line = br.readLine();
                        }
                    }

                    return (T) builder.build();
                } else {
                    String[] taxa = input.trim().split(",");
                    return (T) new TaxaListBuilder().addAll(taxa).build();
                }
            } else if (outputClass.isAssignableFrom(DistanceMatrix.class)) {
                return (T) ReadDistanceMatrix.readDistanceMatrix(input);
            } else if (outputClass.isAssignableFrom(Character.class)) {
                if (input.length() != 1) {
                    throw new IllegalArgumentException("Should be one character");
                }
                return (T) new Character(input.charAt(0));
            } else {
                return outputClass.getConstructor(String.class).newInstance(input);
            }
        } catch (Exception nfe) {
            myLogger.debug(nfe.getMessage(), nfe);
            String message = nfe.getMessage();
            if (message == null) {
                throw new IllegalArgumentException("Problem converting: " + input + " to " + Utils.getBasename(outputClass.getName()));
            } else {
                throw new IllegalArgumentException(message + " Problem converting: " + input + " to " + Utils.getBasename(outputClass.getName()));
            }
        }
    }

    private static List<String> getListFromString(String str) {

        if ((str == null) || (str.length() == 0) || str.equalsIgnoreCase("null")) {
            return null;
        }
        String[] tokens = str.split(",");
        List<String> result = new ArrayList<>();
        for (String current : tokens) {
            current = current.trim();
            if (current.length() != 0) {
                result.add(current);
            }
        }
        return result;

    }

    private static SortedSet<Integer> getSortedSet(String str) {

        if ((str == null) || (str.length() == 0) || str.equalsIgnoreCase("null")) {
            return null;
        }
        String[] tokens = str.split(",");
        SortedSet<Integer> result = new TreeSet<>();
        for (String current : tokens) {
            current = current.trim();
            if (current.length() != 0) {
                result.add(Integer.valueOf(current));
            }
        }
        return result;

    }

    @Override
    public void setParameters(String[] args) {

        if (args != null) {

            for (int i = 0; i < args.length; i++) {
                String arg = args[i];
                if (args[i].startsWith("-")) {
                    arg = arg.substring(1);
                    PluginParameter<?> parameter = getParameterInstance(arg);
                    if (parameter == null) {
                        myLogger.error("Unrecognized argument: " + args[i]);
                        printUsage();
                        System.exit(1);
                    }
                    if ((i == args.length - 1) || (args[i + 1]).startsWith("-")) {
                        if (parameter.valueType().isAssignableFrom(Boolean.class)) {
                            setParameter(arg, Boolean.TRUE);
                        } else if (Number.class.isAssignableFrom(parameter.valueType())) {
                            setParameter(arg, args[i + 1]);
                            i++;
                        } else {
                            myLogger.error("Parameter requires a value: " + args[i]);
                            printUsage();
                            System.exit(1);
                        }
                    } else {
                        setParameter(arg, args[i + 1]);
                        i++;
                    }
                } else {
                    myLogger.error("Argument expected to start with dash(-): " + args[i]);
                    printUsage();
                    System.exit(1);
                }
            }

        }

    }

    private void setFieldsToConfigParameters(Map<String, JComponent> parameterFields) {

        final List<PluginParameter<?>> parameterInstances = getParameterInstances();
        if (parameterInstances.isEmpty()) {
            return;
        }

        for (final PluginParameter<?> current : parameterInstances) {
            JComponent component = parameterFields.get(current.cmdLineName());
            setFieldToConfigParameters(component, current);
        }

    }

    private void setFieldToConfigParameters(JComponent component, PluginParameter<?> parameter) {

        Optional<String> configValue = ParameterCache.value(this, parameter.cmdLineName());
        if (!configValue.isPresent()) {
            return;
        }
        try {
            if (component instanceof JTextField) {
                ((JTextField) component).setText(configValue.get());
            } else if (component instanceof JCheckBox) {
                Boolean value = convert(configValue.get(), Boolean.class);
                ((JCheckBox) component).setSelected(value);
            } else if (component instanceof JComboBox) {
                Object value = convert(configValue.get(), parameter.valueType());
                ((JComboBox) component).setSelectedItem(value);
            }
        } catch (Exception e) {
            myLogger.warn("setFieldToConfigParameters: problem with configuration key: " + this.getClass().getName() + "." + parameter.cmdLineName() + "  value: " + configValue.get() + "\n" + e.getMessage());
        }

    }

    public void setConfigParameters() {

        if (ParameterCache.hasValues()) {

            for (PluginParameter<?> parameter : getParameterInstances()) {
                Optional<String> value = ParameterCache.value(this, parameter.cmdLineName());
                if (value.isPresent()) {
                    try {
                        setParameter(parameter, value.get());
                    } catch (Exception e) {
                        myLogger.warn("setConfigParameters: problem with configuration key: " + this.getClass().getName() + "." + parameter.cmdLineName() + "  value: " + value.get() + "\n" + e.getMessage());
                    }
                }
            }

        }

    }

    private void checkRequiredParameters() {

        List<String> cmdLineNames = new ArrayList<>();
        for (PluginParameter<?> current : getParameterInstances()) {
            if (cmdLineNames.contains(current.cmdLineName())) {
                if (isInteractive()) {
                    throw new IllegalStateException(current.cmdLineName() + " exist multiple times for this plugin.");
                } else {
                    myLogger.error("-" + current.cmdLineName() + " exist multiple times for this plugin.\n");
                    printUsage();
                    System.exit(1);
                }
            } else {
                cmdLineNames.add(current.cmdLineName());
            }

            if (current.required()) {
                if (current.isEmpty()) {
                    if (isInteractive()) {
                        throw new IllegalStateException(current.guiName() + " must be defined.");
                    } else {
                        myLogger.error("-" + current.cmdLineName() + " is required.\n");
                        printUsage();
                        System.exit(1);
                    }
                }
            }
        }

    }

    /**
     * Verification checks of parameters.
     */
    private void checkParameters() {

        for (PluginParameter<?> current : getParameterInstances()) {

            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.IN_FILE) {
                if (!current.isEmpty()) {
                    String filename = current.value().toString();
                    File theFile = new File(filename);
                    if (!theFile.exists()) {
                        if (isInteractive()) {
                            throw new IllegalStateException(current.guiName() + ": " + filename + " doesn't exist.");
                        } else {
                            myLogger.error("-" + current.cmdLineName() + ": " + filename + " doesn't exist\n");
                            printUsage();
                            System.exit(1);
                        }
                    }
                    if (!theFile.isFile()) {
                        if (isInteractive()) {
                            throw new IllegalStateException(current.guiName() + ": " + filename + " isn't a file.");
                        } else {
                            myLogger.error("-" + current.cmdLineName() + ": " + filename + " isn't a file\n");
                            printUsage();
                            System.exit(1);
                        }
                    }
                }
            }

            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.OUT_FILE) {
                if (!current.isEmpty()) {
                    String filename = current.value().toString();
                    String outFolder = Utils.getDirectory(filename);
                    File outDir = new File(outFolder);
                    if (!outDir.isDirectory()) {
                        if (isInteractive()) {
                            throw new IllegalStateException(current.guiName() + ": Output Directory: " + outFolder + " doesn't exist.");
                        } else {
                            myLogger.error("-" + current.cmdLineName() + ": Output Directory: " + outFolder + " doesn't exist\n");
                            printUsage();
                            System.exit(1);
                        }
                    }
                }
            }

            if ((current.parameterType() == PluginParameter.PARAMETER_TYPE.IN_DIR)
                    || (current.parameterType() == PluginParameter.PARAMETER_TYPE.OUT_DIR)) {
                if (!current.isEmpty()) {
                    String dirname = current.value().toString();
                    File directory = new File(dirname);
                    if (!directory.isDirectory()) {
                        if (isInteractive()) {
                            throw new IllegalStateException(current.guiName() + ": Directory: " + dirname + " doesn't exist.");
                        } else {
                            myLogger.error("-" + current.cmdLineName() + ": Directory: " + dirname + " doesn't exist\n");
                            printUsage();
                            System.exit(1);
                        }
                    }
                }
            }

        }

    }

    protected void printParameterValues() {
        List<PluginParameter<?>> parameters = getParameterInstances();
        if ((parameters == null) || (parameters.isEmpty())) {
            return;
        }
        StringBuilder builder = new StringBuilder();
        builder.append("\n");
        builder.append(Utils.getBasename(getClass().getName()));
        builder.append(" Parameters\n");
        for (PluginParameter<?> current : parameters) {
            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                continue;
            }
            builder.append(current.cmdLineName());
            builder.append(": ");
            Object value = current.value();
            if (value instanceof PositionList) {
                builder.append(((PositionList) value).numberOfSites());
                builder.append(" positions");
            } else if (value instanceof List) {
                builder.append((Arrays.toString(((List) value).toArray())));
            } else if (value instanceof DistanceMatrix) {
                DistanceMatrix temp = (DistanceMatrix) value;
                builder.append(temp.getColumnCount());
                builder.append(" columns x ");
                builder.append(temp.getRowCount());
                builder.append(" rows");
            } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.PASSWORD) {
                builder.append("?????");
            } else {
                builder.append(value);
            }
            builder.append("\n");
        }
        myLogger.info(builder.toString());
    }

    private void printUsage() {

        StringBuilder builder = new StringBuilder();
        String description = pluginDescription();
        if (description != null) {
            builder.append("\n");
            builder.append(Utils.getBasename(getClass().getName())).append(" Description...\n");
            builder.append(description);
            builder.append("\n");
        }
        builder.append("\nUsage:\n");
        builder.append(Utils.getBasename(getClass().getName())).append(" <options>\n");
        for (PluginParameter<?> current : getParameterInstances()) {
            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                continue;
            }
            builder.append("-");
            builder.append(current.cmdLineName());
            builder.append(" ");
            if (current.valueType().isAssignableFrom(Boolean.class)) {
                builder.append("<true | false>");
            } else {
                builder.append("<");
                builder.append(current.guiName());
                builder.append(">");
            }
            builder.append(" : ");
            builder.append(current.description());
            if (current.hasRange()) {
                builder.append(" ");
                builder.append(current.rangeToString());
            }
            if (current.defaultValue() != null) {
                builder.append(" (Default: ");
                builder.append(current.defaultValue());
                builder.append(")");
            }
            if (current.required()) {
                builder.append(" (required)");
            }
            builder.append("\n");
        }

        myLogger.info(builder.toString());
    }

    @Override
    public String getUsage() {

        StringBuilder builder = new StringBuilder();
        builder.append(Utils.getBasename(getClass().getName()));
        builder.append("\n");
        String description = pluginDescription();
        if (description != null) {
            builder.append("\nDescription: ");
            builder.append(description);
            builder.append("\n\n");
        }
        for (PluginParameter<?> current : getParameterInstances()) {
            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                continue;
            }
            builder.append("\n");
            builder.append(current.guiName());
            builder.append(" : ");
            builder.append(current.description());
            if (current.hasRange()) {
                builder.append(" ");
                builder.append(current.rangeToString());
            }
            if (current.defaultValue() != null) {
                builder.append(" (Default: ");
                builder.append(current.defaultValue());
                builder.append(")");
            }
            if (current.required()) {
                builder.append(" (required)");
            }
            builder.append("\n");
        }

        return builder.toString();
    }

    public Map<String, Map<PARAMETER_PROPERTIES, String>> usageParameters() {

        Map<String, Map<PARAMETER_PROPERTIES, String>> result = new HashMap<>();
        for (PluginParameter<?> current : getParameterInstances()) {
            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                continue;
            }
            Map<PARAMETER_PROPERTIES, String> temp = new HashMap<>();
            temp.put(PARAMETER_PROPERTIES.Required, Boolean.toString(current.required()));
            String defaultValue = (current.defaultValue() == null) ? null : current.defaultValue().toString();
            temp.put(PARAMETER_PROPERTIES.Default, defaultValue);
            temp.put(PARAMETER_PROPERTIES.Description, current.description());
            result.put(current.cmdLineName(), temp);
        }

        return result;

    }

    @Override
    public Map<String, String> pluginParameters() {

        Map<String, String> result = new HashMap<>();

        for (PluginParameter<?> current : getParameterInstances()) {
            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                continue;
            }
            result.put(current.cmdLineName(), (current.value() == null)? null : current.value().toString());
        }

        return result;

    }

    @Override
    public Object getParameter(Enum key) {
        return getParameterInstance(key.toString()).value();
    }

    @Override
    public Object getParameter(String key) {
        return getParameterInstance(key).value();
    }

    @Override
    public Plugin setParameter(PluginParameter<?> param, Object value) {
        if (value == null) {
            setParameter(param.cmdLineName(), value);
        } else if (value instanceof String) {
            setParameter(param.cmdLineName(), (String) value);
        } else {
            setParameter(param.cmdLineName(), value);
        }
        return this;
    }

    @Override
    public Plugin setParameter(String key, Object value) {

        try {

            Field field = getParameterField(key);
            PluginParameter parameter = (PluginParameter) field.get(this);
            if (parameter == null) {
                throw new IllegalArgumentException("setParameter: Unknown Parameter: " + key);
            }
            if ((parameter.hasRange()) && (!parameter.acceptsValue(value))) {
                throw new IllegalArgumentException("setParameter: " + parameter.cmdLineName() + " value: " + value.toString() + " outside range: " + parameter.rangeToString());
            }
            PluginParameter newParameter = new PluginParameter<>(parameter, value);
            field.set(this, newParameter);

        } catch (Exception e) {
            if (isInteractive()) {
                try {
                    throw e;
                } catch (IllegalAccessException ex) {
                    myLogger.error(ex.getMessage(), ex);
                }
            } else {
                myLogger.error(key + ": " + e.getMessage());
                printUsage();
                myLogger.debug(e.getMessage(), e);
                System.exit(1);
            }
        }

        return this;
    }

    @Override
    public Plugin setParameter(String key, String value) {

        try {
            PluginParameter parameter = getParameterInstance(key);
            return setParameter(key, convert(value, parameter.valueType()));
        } catch (Exception e) {
            if (isInteractive()) {
                throw new IllegalArgumentException(getParameterInstance(key).guiName() + ": " + e.getMessage());
            } else {
                myLogger.error(key + ": " + e.getMessage());
                printUsage();
                System.exit(1);
            }
        }
        return this;
    }

    private static final int TEXT_FIELD_WIDTH = 25;

    private boolean parametersAreSet = true;

    /**
     * Generates dialog based on this plugins define parameters.
     *
     * @return true if OK clicked, false if canceled
     */
    private boolean setParametersViaGUI() {

        final List<PluginParameter<?>> parameterInstances = getParameterInstances();
        if (parameterInstances.isEmpty()) {
            return true;
        }

        final JDialog dialog = new JDialog(getParentFrame(), null, true);
        dialog.addWindowListener(new java.awt.event.WindowAdapter() {
            @Override
            public void windowClosing(WindowEvent e) {
                parametersAreSet = false;
                dialog.setVisible(false);
            }
        });

        final Map<String, JComponent> parameterFields = new HashMap<>();

        parametersAreSet = true;

        JButton okButton = new JButton();
        okButton.setActionCommand("Ok");
        okButton.setText("Ok");
        okButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    for (final PluginParameter<?> current : parameterInstances) {
                        JComponent component = parameterFields.get(current.cmdLineName());
                        if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                            // do nothing
                        } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.GENOTYPE_TABLE) {
                            GenotypeWrapper input = (GenotypeWrapper) ((JComboBox) component).getSelectedItem();
                            if (input != null) {
                                setParameter(current.cmdLineName(), input.myObj);
                            }
                        } else if (PositionList.class.isAssignableFrom(current.valueType())) {
                            if (component instanceof JComboBox) {
                                Object temp = ((JComboBox) component).getSelectedItem();
                                if (temp == POSITION_LIST_NONE) {
                                    setParameter(current.cmdLineName(), null);
                                } else {
                                    setParameter(current.cmdLineName(), ((Datum) temp).getData());
                                }
                            } else {
                                String input = ((JTextField) component).getText().trim();
                                setParameter(current.cmdLineName(), input);
                            }
                        } else if (TaxaList.class.isAssignableFrom(current.valueType())) {
                            if (component instanceof JComboBox) {
                                Object temp = ((JComboBox) component).getSelectedItem();
                                if (temp == TAXA_LIST_NONE) {
                                    setParameter(current.cmdLineName(), null);
                                } else {
                                    setParameter(current.cmdLineName(), ((Datum) temp).getData());
                                }
                            } else {
                                String input = ((JTextField) component).getText().trim();
                                setParameter(current.cmdLineName(), input);
                            }
                        } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.DISTANCE_MATRIX) {
                            if (component instanceof JComboBox) {
                                Object temp = ((JComboBox) component).getSelectedItem();
                                if (temp == null) {
                                    throw new IllegalArgumentException("setParametersViaGUI: must specify a distance matrix.");
                                } else {
                                    setParameter(current.cmdLineName(), ((Datum) temp).getData());
                                }
                            } else {
                                String input = ((JTextField) component).getText().trim();
                                setParameter(current.cmdLineName(), input);
                            }
                        } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.OBJECT_LIST_SINGLE_SELECT) {
                            Object selectedObjects = ((JComboBox<?>) component).getSelectedItem();
                            setParameter(current.cmdLineName(), selectedObjects);
                        } else if (component instanceof JTextField) {
                            String input = ((JTextField) component).getText().trim();
                            setParameter(current.cmdLineName(), input);
                        } else if (component instanceof JCheckBox) {
                            if (((JCheckBox) component).isSelected()) {
                                setParameter(current.cmdLineName(), Boolean.TRUE);
                            } else {
                                setParameter(current.cmdLineName(), Boolean.FALSE);
                            }
                        } else if (component instanceof JComboBox) {
                            Enum temp = (Enum) ((JComboBox) component).getSelectedItem();
                            setParameter(current.cmdLineName(), temp);
                        }
                    }
                } catch (Exception ex) {
                    myLogger.debug(ex.getMessage(), ex);
                    StringBuilder builder = new StringBuilder();
                    builder.append("Problem Setting Parameters: ");
                    builder.append("\n");
                    builder.append(Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(ex), 50));
                    String str = builder.toString();
                    DialogUtils.showError(str, getParentFrame());
                    return;
                }
                dialog.setVisible(false);
            }
        });

        JButton cancelButton = new JButton();
        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                parametersAreSet = false;
                dialog.setVisible(false);
            }
        });

        JButton defaultsButton = new JButton();
        defaultsButton.setText("Defaults");
        defaultsButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                setFieldsToDefault(parameterFields);
                setFieldsToConfigParameters(parameterFields);
            }
        });

        JButton userManualButton = new JButton();
        userManualButton.setText("User Manual");
        userManualButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    Desktop desktop = Desktop.getDesktop();
                    URI uri = new URI(pluginUserManualURL());
                    desktop.browse(uri);
                } catch (Exception ex) {
                    myLogger.debug(ex.getMessage(), ex);
                }
            }
        });

        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        boolean show_citation = !DEFAULT_CITATION.equals(getCitation());
        JTextPane citationText = null;
        if (show_citation) {
            citationText = new JTextPane();
            citationText.setContentType("text/html");
            citationText.setMargin(new Insets(5, 5, 5, 5));
            citationText.setEditable(false);
            JScrollPane scroll = new JScrollPane(citationText);
            scroll.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
            scroll.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);
            scroll.setPreferredSize(new Dimension(scroll.getWidth(), 45));
            panel.add(scroll);
        }

        for (final PluginParameter<?> current : getParameterInstances()) {
            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.GENOTYPE_TABLE) {
                Datum datum = getGenotypeTable();
                JComboBox menu = new JComboBox();
                if (datum != null) {
                    String name = datum.getName();
                    GenotypeTable table;
                    if (datum.getData() instanceof GenotypeTable) {
                        table = (GenotypeTable) datum.getData();
                    } else if (datum.getData() instanceof GenotypePhenotype) {
                        table = ((GenotypePhenotype) datum.getData()).genotypeTable();
                    } else {
                        throw new IllegalStateException("AbstractPlugin: setParametersViaGUI: unknown GenotypeTable type: " + datum.getData().getClass().getName());
                    }

                    if (current.acceptsValue(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype) && table.hasGenotype()) {
                        menu.addItem(new GenotypeWrapper(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, "Genotype (" + name + ")"));
                    }
                    if (current.acceptsValue(GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability) && table.hasReferenceProbablity()) {
                        menu.addItem(new GenotypeWrapper(GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability, "Reference Probability (" + name + ")"));
                    }
                    if (current.acceptsValue(GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability) && table.hasAlleleProbabilities()) {
                        menu.addItem(new GenotypeWrapper(GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability, "Allele Probability (" + name + ")"));
                    }
                    if (current.acceptsValue(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Depth) && table.hasDepth()) {
                        menu.addItem(new GenotypeWrapper(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Depth, "Depth (" + name + ")"));
                    }
                    if (current.acceptsValue(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Dosage) && table.hasDosage()) {
                        menu.addItem(new GenotypeWrapper(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Dosage, "Dosage (" + name + ")"));
                    }
                    menu.setSelectedIndex(0);
                }
                createEnableDisableAction(current, parameterFields, menu);
                JPanel temp = new JPanel(new FlowLayout(FlowLayout.CENTER));
                temp.add(new JLabel(current.guiName()));
                temp.add(menu);
                temp.setToolTipText(getToolTip(current));
                panel.add(temp);
                parameterFields.put(current.cmdLineName(), menu);
            } else if (PositionList.class.isAssignableFrom(current.valueType())) {
                Datum datum = getPositionList();
                if (datum != null) {
                    JComboBox menu = new JComboBox();
                    menu.addItem(POSITION_LIST_NONE);
                    menu.addItem(datum);
                    menu.setSelectedIndex(0);
                    createEnableDisableAction(current, parameterFields, menu);
                    JPanel temp = new JPanel(new FlowLayout(FlowLayout.CENTER));
                    temp.add(new JLabel(current.guiName()));
                    temp.add(menu);
                    temp.setToolTipText(getToolTip(current));
                    panel.add(temp);
                    parameterFields.put(current.cmdLineName(), menu);
                } else {
                    JTextField field = new JTextField(TEXT_FIELD_WIDTH - 8);
                    JButton browse = getOpenFile(dialog, field);
                    JPanel line = getLine(current.guiName(), field, browse, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, new JComponent[]{field, browse}, field);
                    panel.add(line);
                    parameterFields.put(current.cmdLineName(), field);
                }
            } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.DISTANCE_MATRIX) {
                List<Datum> matrices = getDistanceMatrices();
                if (!matrices.isEmpty()) {
                    JComboBox menu = new JComboBox();
                    for (Datum matrix : matrices) {
                        menu.addItem(matrix);
                    }
                    menu.setSelectedIndex(0);
                    createEnableDisableAction(current, parameterFields, menu);
                    JPanel temp = new JPanel(new FlowLayout(FlowLayout.CENTER));
                    temp.add(new JLabel(current.guiName()));
                    temp.add(menu);
                    temp.setToolTipText(getToolTip(current));
                    panel.add(temp);
                    parameterFields.put(current.cmdLineName(), menu);
                } else {
                    JTextField field = new JTextField(TEXT_FIELD_WIDTH - 8);
                    JButton browse = getOpenFile(dialog, field);
                    JPanel line = getLine(current.guiName(), field, browse, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, new JComponent[]{field, browse}, field);
                    panel.add(line);
                    parameterFields.put(current.cmdLineName(), field);
                }
            } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.OBJECT_LIST_SINGLE_SELECT) {
                JPanel listPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
                listPanel.add(new JLabel(current.guiName()));
                JComboBox<?> list = new JComboBox<>(current.possibleValues().toArray());
                list.setSelectedIndex(0);
                createEnableDisableAction(current, parameterFields, list);
                listPanel.add(list);
                listPanel.setToolTipText(getToolTip(current));
                panel.add(listPanel);
                parameterFields.put(current.cmdLineName(), list);
            } else if (current.valueType().isEnum()) {
                JComboBox menu = new JComboBox();
                Object[] values = current.valueType().getEnumConstants();
                for (Object item : values) {
                    menu.addItem(item);
                }
                menu.setSelectedItem(current.value());
                createEnableDisableAction(current, parameterFields, menu);
                JPanel temp = new JPanel(new FlowLayout(FlowLayout.CENTER));
                temp.add(new JLabel(current.guiName()));
                temp.add(menu);
                temp.setToolTipText(getToolTip(current));
                panel.add(temp);
                parameterFields.put(current.cmdLineName(), menu);
            } else if (Boolean.class.isAssignableFrom(current.valueType())) {
                JCheckBox check = new JCheckBox(current.guiName());
                check.setToolTipText(getToolTip(current));
                if (current.value() == Boolean.TRUE) {
                    check.setSelected(true);
                } else {
                    check.setSelected(false);
                }
                createEnableDisableAction(current, parameterFields, check);
                JPanel temp = new JPanel(new FlowLayout(FlowLayout.CENTER));
                temp.add(check);
                panel.add(temp);
                parameterFields.put(current.cmdLineName(), check);
            } else if (TaxaList.class.isAssignableFrom(current.valueType())) {

                List<Datum> datums = getTaxaListDatum();
                if (datums != null) {
                    JComboBox menu = new JComboBox();
                    menu.addItem(TAXA_LIST_NONE);
                    for (Datum datum : datums) {
                        menu.addItem(datum);
                    }
                    menu.setSelectedIndex(0);
                    createEnableDisableAction(current, parameterFields, menu);
                    JPanel temp = new JPanel(new FlowLayout(FlowLayout.CENTER));
                    temp.add(new JLabel(current.guiName()));
                    temp.add(menu);
                    temp.setToolTipText(getToolTip(current));
                    panel.add(temp);
                    parameterFields.put(current.cmdLineName(), menu);
                } else {
                    TaxaList taxa = getTaxaList();
                    JTextField field;
                    if (taxa == null) {
                        field = new JTextField(TEXT_FIELD_WIDTH);
                    } else {
                        field = new JTextField(TEXT_FIELD_WIDTH - 7);
                    }
                    if (current.value() != null) {
                        field.setText(current.value().toString());
                    }
                    JPanel taxaPanel = getTaxaListPanel(current.guiName(), field, current.description(), dialog, taxa);
                    panel.add(taxaPanel);
                    parameterFields.put(current.cmdLineName(), field);
                }

            } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.SITE_NAME_LIST) {
                PositionList positions = getSiteNameList();
                JTextField field;
                if (positions == null) {
                    field = new JTextField(TEXT_FIELD_WIDTH);
                } else {
                    field = new JTextField(TEXT_FIELD_WIDTH - 7);
                }
                if (current.value() != null) {
                    field.setText(current.value().toString());
                }
                JPanel positionsPanel = getPositionListPanel(current.guiName(), field, current.description(), dialog, positions);
                List<JComponent> componentList = new ArrayList<>();
                for (Component component : positionsPanel.getComponents()) {
                    if (component instanceof JComponent) {
                        componentList.add((JComponent) component);
                    }
                }
                createEnableDisableAction(current, parameterFields, componentList.toArray(new JComponent[0]), field);
                panel.add(positionsPanel);
                parameterFields.put(current.cmdLineName(), field);
            } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                JPanel temp = new JPanel(new FlowLayout(FlowLayout.CENTER));
                JLabel label = new JLabel(current.guiName());
                label.setFont(new Font("Dialog", Font.BOLD, 14));
                temp.add(label);
                panel.add(temp);
            } else {
                final JTextField field;
                if (current.parameterType() == PluginParameter.PARAMETER_TYPE.PASSWORD) {
                    field = new JPasswordField(TEXT_FIELD_WIDTH);
                } else if (current.parameterType() != PluginParameter.PARAMETER_TYPE.NA) {
                    field = new JTextField(TEXT_FIELD_WIDTH - 8);
                } else {
                    field = new JTextField(TEXT_FIELD_WIDTH);
                }

                if (current.value() != null) {
                    setFieldToValue(field, current, current.value());
                }

                field.addFocusListener(new FocusAdapter() {
                    @Override
                    public void focusLost(FocusEvent e) {
                        String input = field.getText().trim();
                        try {
                            PluginParameter parameter = getParameterInstance(current.cmdLineName());
                            if (!current.acceptsValue(input)) {
                                JOptionPane.showMessageDialog(dialog, current.guiName() + " range: " + current.rangeToString());
                                setFieldToValue(field, parameter, parameter.value());
                            } else {
                                setFieldToValue(field, parameter, convert(input, parameter.valueType()));
                            }
                        } catch (Exception ex) {
                            myLogger.debug(ex.getMessage(), ex);
                            JOptionPane.showMessageDialog(dialog, current.guiName() + ": " + ex.getMessage());
                            field.setText(getParameterInstance(current.cmdLineName()).value().toString());
                        }
                    }
                });

                String label = null;
                if (current.required()) {
                    label = current.guiName() + "*";
                } else {
                    label = current.guiName();
                }

                JPanel line = null;
                if (current.parameterType() == PluginParameter.PARAMETER_TYPE.IN_FILE) {
                    JButton browse = getOpenFile(dialog, field);
                    line = getLine(label, field, browse, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, new JComponent[]{field, browse}, field);
                } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.OUT_FILE) {
                    JButton browse = getSaveFile(dialog, field);
                    line = getLine(label, field, browse, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, new JComponent[]{field, browse}, field);
                } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.IN_DIR) {
                    JButton browse = getOpenDir(dialog, field);
                    line = getLine(label, field, browse, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, new JComponent[]{field, browse}, field);
                } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.OUT_DIR) {
                    JButton browse = getSaveDir(dialog, field);
                    line = getLine(label, field, browse, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, new JComponent[]{field, browse}, field);
                } else {
                    line = getLine(label, field, null, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, field);
                }
                panel.add(line);

                parameterFields.put(current.cmdLineName(), field);
            }
        }

        JTabbedPane tabbedPane = new JTabbedPane();
        tabbedPane.add(new JScrollPane(panel, ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED, ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER), getButtonName());

        JPanel pnlButtons = new JPanel();
        pnlButtons.setLayout(new FlowLayout());
        pnlButtons.add(okButton);
        pnlButtons.add(cancelButton);
        pnlButtons.add(defaultsButton);
        pnlButtons.add(userManualButton);
        dialog.getContentPane().add(tabbedPane, BorderLayout.CENTER);
        dialog.getContentPane().add(pnlButtons, BorderLayout.SOUTH);

        JTextPane helpText = new JTextPane();
        helpText.setMargin(new Insets(10, 10, 10, 10));
        helpText.setEditable(false);
        helpText.setContentType("text/html");
        tabbedPane.add(new JScrollPane(helpText), "Help");
        dialog.pack();
        if (show_citation) {
            citationText.setText(getCitationHTML(dialog.getWidth() / 9));
            dialog.setMinimumSize(null);
            dialog.pack();
        }
        helpText.setText(getUsageHTML());
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        if (screenSize.getHeight() - 125 < dialog.getHeight()) {
            dialog.setSize(Math.max(dialog.getWidth(), 550), (int) screenSize.getHeight() - 125);
        } else {
            dialog.setSize(Math.max(dialog.getWidth(), 550), Math.max(dialog.getHeight(), 250));
        }
        dialog.setResizable(false);
        dialog.setLocationRelativeTo(getParentFrame());
        dialog.setVisible(true);
        return parametersAreSet;

    }

    /**
     * This sets the GUI field text to appropriate formatting for integers, doubles, floats, strings, etc. for the
     * user's default locale.
     *
     * @param field text field
     * @param parameter associated plugin parameter
     * @param value value to format
     */
    private void setFieldToValue(JTextField field, PluginParameter<?> parameter, Object value) {
        if (value == null) {
            field.setText(null);
        } else if (Integer.class.isAssignableFrom(parameter.valueType())) {
            field.setText(NumberFormat.getInstance().format(value));
        } else if ((Double.class.isAssignableFrom(parameter.valueType())) ||
                (Float.class.isAssignableFrom(parameter.valueType()))) {
            DecimalFormat temp = new DecimalFormat();
            temp.setMaximumFractionDigits(5);
            temp.setMinimumFractionDigits(1);
            field.setText(temp.format(value));
        } else {
            field.setText(value.toString());
        }
    }

    private void setFieldsToDefault(Map<String, JComponent> parameterFields) {

        final List<PluginParameter<?>> parameterInstances = getParameterInstances();
        if (parameterInstances.isEmpty()) {
            return;
        }

        for (final PluginParameter<?> current : parameterInstances) {
            JComponent component = parameterFields.get(current.cmdLineName());
            setFieldToDefault(component, current);
        }

    }

    private void setFieldToDefault(JComponent component, PluginParameter<?> current) {
        if (component instanceof JTextField) {
            Object defaultValue = current.defaultValue();
            setFieldToValue((JTextField) component, current, defaultValue);
            setParameter(current.cmdLineName(), defaultValue);
        } else if (component instanceof JCheckBox) {
            Boolean value = (Boolean) current.defaultValue();
            ((JCheckBox) component).setSelected(value);
            setParameter(current.cmdLineName(), value);
        } else if (component instanceof JComboBox) {
            ((JComboBox) component).setSelectedItem(current.defaultValue());
            setParameter(current.cmdLineName(), current.defaultValue());
        }
    }

    @Override
    public void setParametersToDefault() {

        final List<PluginParameter<?>> parameterInstances = getParameterInstances();
        if (parameterInstances.isEmpty()) {
            return;
        }

        for (final PluginParameter<?> current : parameterInstances) {
            setParameter(current.cmdLineName(), current.defaultValue());
        }

    }

    private String getCitationHTML(int lineWidth) {
        String citation = getCitation();
        int count = 10;
        StringBuilder builder = new StringBuilder();
        builder.append("<html><center>Citation: ");
        for (int i = 0, n = citation.length(); i < n; i++) {
            count++;
            if (citation.charAt(i) == '\n') {
                builder.append("<br>");
                count = 0;
            } else if ((count > lineWidth) && (citation.charAt(i) == ' ')) {
                builder.append("<br>");
                count = 0;
            } else {
                builder.append(citation.charAt(i));
            }
        }
        builder.append("</center></html>");
        return builder.toString();
    }

    public String getUsageHTML() {

        StringBuilder builder = new StringBuilder();
        builder.append("<html><center><strong>");
        builder.append(Utils.getBasename(getClass().getName()));
        builder.append("</strong>");
        String description = pluginDescription();
        if (description != null) {
            builder.append("<br><br><strong>Description:</strong> ");
            builder.append(description);
        }
        builder.append("<br><br>");
        builder.append("<table border='1'>");
        builder.append("<tr><th>Parameter</th><th>Description</th><th>Values</th><th>Default</th></tr>");
        for (PluginParameter<?> current : getParameterInstances()) {
            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                continue;
            }
            builder.append("<tr>");

            builder.append("<th>");
            if (current.required()) {
                builder.append("<font color='red'>");
            }
            builder.append(current.guiName());
            if (current.required()) {
                builder.append("</font>");
            }
            builder.append("</th>");

            builder.append("<td>");
            builder.append(current.description());
            builder.append("</td>");

            builder.append("<td>");
            if (current.hasPossibleValues()) {
                String range = current.possibleValuesString(true);
                if ((range.charAt(0) == '[') && (range.charAt(range.length() - 1) == ']')) {
                    range = range.substring(1, range.length() - 1);
                }
                StringBuilder buildRange = new StringBuilder();
                for (char rangeChr : range.toCharArray()) {
                    if (rangeChr == '_') {
                        buildRange.append("\n");
                    }
                    buildRange.append(rangeChr);
                }
                builder.append(buildRange.toString());
            } else if (current.hasRange()) {
                String range = current.rangeToString(true);
                if ((range.charAt(0) == '[') && (range.charAt(range.length() - 1) == ']')) {
                    range = range.substring(1, range.length() - 1);
                }
                StringBuilder buildRange = new StringBuilder();
                for (char rangeChr : range.toCharArray()) {
                    if (rangeChr == '_') {
                        buildRange.append("\n");
                    }
                    buildRange.append(rangeChr);
                }
                builder.append(buildRange.toString());
            } else if (current.valueType().isAssignableFrom(Boolean.class)) {
                builder.append("true, false");
            }
            builder.append("</td>");

            builder.append("<td>");
            if (current.defaultValue() != null) {
                String defaultValue = current.defaultValue().toString();
                StringBuilder buildDefault = new StringBuilder();
                for (char defaultChr : defaultValue.toCharArray()) {
                    if (defaultChr == '_') {
                        buildDefault.append("\n");
                    }
                    buildDefault.append(defaultChr);
                }
                builder.append(buildDefault.toString());
            }
            builder.append("</td>");

            builder.append("</tr>");
        }
        builder.append("</table>");

        builder.append("</center>");

        builder.append("<br><font color='red'>* parameters in red are required</font>");

        builder.append("</html>");

        return builder.toString();
    }

    private void createEnableDisableAction(PluginParameter<?> current, Map<String, JComponent> parameterFields, final JComponent component) {
        createEnableDisableAction(current, parameterFields, new JComponent[]{component}, component);
    }

    private void createEnableDisableAction(final PluginParameter<?> current, Map<String, JComponent> parameterFields, final JComponent[] components, final JComponent input) {

        if (current.dependentOnParameter() != null) {
            JComponent depends = parameterFields.get(current.dependentOnParameter().cmdLineName());
            if (depends instanceof JCheckBox) {
                final JCheckBox checkBox = (JCheckBox) depends;

                for (JComponent component : components) {
                    if (checkBox.isSelected() == (Boolean) current.dependentOnParameterValue()[0]) {
                        component.setEnabled(true);
                    } else {
                        component.setEnabled(false);
                    }
                }

                checkBox.addItemListener(new ItemListener() {

                    @Override
                    public void itemStateChanged(ItemEvent e) {
                        for (JComponent component : components) {
                            if (checkBox.isSelected() == (Boolean) current.dependentOnParameterValue()[0]) {
                                component.setEnabled(true);
                            } else {
                                component.setEnabled(false);
                            }
                        }
                    }
                });

            } else if (depends instanceof JComboBox) {
                final JComboBox comboBox = (JComboBox) depends;

                for (JComponent component : components) {
                    Object[] values = current.dependentOnParameterValue();
                    component.setEnabled(false);
                    for (Object value : values) {
                        if (comboBox.getSelectedItem() == value) {
                            component.setEnabled(true);
                            break;
                        }
                    }
                }

                comboBox.addItemListener((ItemEvent e) -> {
                    for (JComponent component : components) {
                        Object[] values = current.dependentOnParameterValue();
                        component.setEnabled(false);
                        for (Object value : values) {
                            if (comboBox.getSelectedItem() == value) {
                                component.setEnabled(true);
                                break;
                            }
                        }
                    }
                });
            }
        }

    }

    private static final int DEFAULT_TOOL_TIP_LINE_LENGTH = 50;

    private String getToolTip(PluginParameter<?> parameter) {
        String description = parameter.description();
        int count = 0;
        StringBuilder builder = new StringBuilder();
        builder.append("<html>");
        for (int i = 0, n = description.length(); i < n; i++) {
            count++;
            if (description.charAt(i) == '\n') {
                builder.append("<br>");
                count = 0;
            } else if ((count > DEFAULT_TOOL_TIP_LINE_LENGTH) && (description.charAt(i) == ' ')) {
                builder.append("<br>");
                count = 0;
            } else {
                builder.append(description.charAt(i));
            }
        }
        builder.append("</html>");
        return builder.toString();
    }

    private JPanel getLine(String label, JTextField ref, JButton button, String description) {

        JPanel result = new JPanel(new FlowLayout(FlowLayout.RIGHT));
        result.setToolTipText(description);

        result.add(new JLabel(label));
        ref.setEditable(true);
        ref.setHorizontalAlignment(JTextField.LEFT);
        ref.setAlignmentX(JTextField.CENTER_ALIGNMENT);
        ref.setAlignmentY(JTextField.CENTER_ALIGNMENT);
        ref.setMaximumSize(ref.getPreferredSize());
        result.add(ref);
        if (button != null) {
            result.add(button);
        }

        return result;

    }

    private JPanel getTaxaListPanel(String label, final JTextField ref, String description, final JDialog parent, final TaxaList taxa) {

        JPanel result = new JPanel(new FlowLayout(FlowLayout.RIGHT));
        result.setToolTipText(description);

        result.add(new JLabel(label));
        ref.setEditable(true);
        ref.setHorizontalAlignment(JTextField.LEFT);
        ref.setAlignmentX(JTextField.CENTER_ALIGNMENT);
        ref.setAlignmentY(JTextField.CENTER_ALIGNMENT);
        ref.setMaximumSize(ref.getPreferredSize());
        result.add(ref);

        if (taxa != null) {
            final SelectFromAvailableDialog dialog = new SelectFromAvailableDialog(getParentFrame(), "Taxa Filter", new TaxaAvailableListModel(taxa));
            JButton taxaButton = new JButton(new AbstractAction() {

                @Override
                public void actionPerformed(ActionEvent e) {
                    dialog.setLocationRelativeTo(parent);
                    dialog.setVisible(true);
                    if (!dialog.isCanceled()) {
                        int[] indicesToKeep = dialog.getDesiredIndices();
                        StringBuilder builder = new StringBuilder();
                        for (int i = 0; i < indicesToKeep.length; i++) {
                            if (i != 0) {
                                builder.append(",");
                            }
                            builder.append(taxa.taxaName(indicesToKeep[i]));
                        }
                        ref.setText(builder.toString());
                    }
                    dialog.setVisible(false);
                }
            });
            taxaButton.setText("Select");
            result.add(taxaButton);
        }

        JButton browse = getOpenFile(parent, ref);
        result.add(browse);

        return result;

    }

    private JPanel getPositionListPanel(String label, final JTextField ref, String description, final JDialog parent, final PositionList positions) {

        JPanel result = new JPanel(new FlowLayout(FlowLayout.RIGHT));
        result.setToolTipText(description);

        result.add(new JLabel(label));
        ref.setEditable(true);
        ref.setHorizontalAlignment(JTextField.LEFT);
        ref.setAlignmentX(JTextField.CENTER_ALIGNMENT);
        ref.setAlignmentY(JTextField.CENTER_ALIGNMENT);
        ref.setMaximumSize(ref.getPreferredSize());
        result.add(ref);

        if (positions != null) {
            final SelectFromAvailableDialog dialog = new SelectFromAvailableDialog(getParentFrame(), "Site Name Filter", new SiteNamesAvailableListModel(positions));
            JButton siteNamesButton = new JButton(new AbstractAction() {

                @Override
                public void actionPerformed(ActionEvent e) {
                    dialog.setLocationRelativeTo(parent);
                    dialog.setVisible(true);
                    if (!dialog.isCanceled()) {
                        int[] indicesToKeep = dialog.getDesiredIndices();
                        StringBuilder builder = new StringBuilder();
                        for (int i = 0; i < indicesToKeep.length; i++) {
                            if (i != 0) {
                                builder.append(",");
                            }
                            builder.append(positions.siteName(indicesToKeep[i]));
                        }
                        ref.setText(builder.toString());
                    }
                    dialog.setVisible(false);
                }
            });
            siteNamesButton.setText("Select");
            result.add(siteNamesButton);
        }

        return result;

    }

    private JButton getOpenFile(final JDialog parent, final JTextField textField) {

        final JFileChooser fileChooser = new JFileChooser(TasselPrefs.getOpenDir());

        JButton result = new JButton("Browse");

        result.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                if (fileChooser.showOpenDialog(parent) == JFileChooser.APPROVE_OPTION) {
                    File file = fileChooser.getSelectedFile();
                    textField.setText(file.getPath());
                    TasselPrefs.putOpenDir(fileChooser.getCurrentDirectory().getPath());
                }
            }

        });

        return result;
    }

    private JButton getSaveFile(final JDialog parent, final JTextField textField) {

        final JFileChooser fileChooser = new JFileChooser(TasselPrefs.getSaveDir());

        JButton result = new JButton("Browse");

        result.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                if (fileChooser.showSaveDialog(parent) == JFileChooser.APPROVE_OPTION) {
                    File file = fileChooser.getSelectedFile();
                    textField.setText(file.getPath());
                    TasselPrefs.putSaveDir(fileChooser.getCurrentDirectory().getPath());
                }
            }

        });

        return result;
    }

    private JButton getOpenDir(final JDialog parent, final JTextField textField) {

        final JFileChooser fileChooser = new JFileChooser(Utils.getDirectory(TasselPrefs.getOpenDir()));
        fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

        JButton result = new JButton("Browse");

        result.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                if (fileChooser.showOpenDialog(parent) == JFileChooser.APPROVE_OPTION) {
                    File file = fileChooser.getSelectedFile();
                    textField.setText(file.getPath());
                    TasselPrefs.putOpenDir(file.getPath());
                }
            }

        });

        return result;
    }

    private JButton getSaveDir(final JDialog parent, final JTextField textField) {

        final JFileChooser fileChooser = new JFileChooser(Utils.getDirectory(TasselPrefs.getSaveDir()));
        fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

        JButton result = new JButton("Browse");

        result.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                if (fileChooser.showOpenDialog(parent) == JFileChooser.APPROVE_OPTION) {
                    File file = fileChooser.getSelectedFile();
                    textField.setText(file.getPath());
                    TasselPrefs.putSaveDir(file.getPath());
                }
            }

        });

        return result;
    }

    private class GenotypeWrapper {

        private final Object myObj;
        private final String myName;

        public GenotypeWrapper(Object obj, String name) {
            myObj = obj;
            myName = name;
        }

        @Override
        public String toString() {
            return myName;
        }

        public Object getObject() {
            return myObj;
        }

    }

    private Datum getGenotypeTable() {

        if (myCurrentInputData == null) {
            return null;
        }

        List<Datum> genotypeTables = myCurrentInputData.getDataOfType(GenotypeTable.class);
        if (!genotypeTables.isEmpty()) {
            return genotypeTables.get(0);
        }

        genotypeTables = myCurrentInputData.getDataOfType(GenotypePhenotype.class);
        if (!genotypeTables.isEmpty()) {
            return genotypeTables.get(0);
        }

        return null;
    }

    private TaxaList getTaxaList() {

        if (myCurrentInputData == null) {
            return null;
        }

        List<Datum> taxaList = myCurrentInputData.getDataOfType(GenotypeTable.class);
        if (!taxaList.isEmpty()) {
            return ((GenotypeTable) taxaList.get(0).getData()).taxa();
        }

        taxaList = myCurrentInputData.getDataOfType(TaxaList.class);
        if (!taxaList.isEmpty()) {
            return (TaxaList) taxaList.get(0).getData();
        }

        return null;
    }

    private List<Datum> getTaxaListDatum() {

        if (myCurrentInputData == null) {
            return null;
        }

        List<Datum> taxaList = myCurrentInputData.getDataOfType(TaxaList.class);
        if (!taxaList.isEmpty()) {
            return taxaList;
        }

        return null;

    }

    private PositionList getSiteNameList() {

        if (myCurrentInputData == null) {
            return null;
        }

        List<Datum> positionList = myCurrentInputData.getDataOfType(GenotypeTable.class);
        if (!positionList.isEmpty()) {
            return ((GenotypeTable) positionList.get(0).getData()).positions();
        }

        positionList = myCurrentInputData.getDataOfType(PositionList.class);
        if (!positionList.isEmpty()) {
            return (PositionList) positionList.get(0).getData();
        }

        return null;
    }

    private Datum getPositionList() {

        if (myCurrentInputData == null) {
            return null;
        }

        List<Datum> positionList = myCurrentInputData.getDataOfType(PositionList.class);
        if (!positionList.isEmpty()) {
            return positionList.get(0);
        }

        return null;
    }

    private List<Datum> getDistanceMatrices() {
        if (myCurrentInputData != null) {
            return myCurrentInputData.getDataOfType(DistanceMatrix.class);
        } else {
            return null;
        }
    }

    /**
     * Returns menu that can be added to main menu bar.
     *
     * @return menu
     */
    @Override
    public JMenu getMenu() {
        return null;
    }

    /**
     * Sets up this plugin to receive input from another plugin.
     *
     * @param input input
     */
    @Override
    public void receiveInput(Plugin input) {

        if (input == null) {
            throw new IllegalArgumentException("AbstractPlugin: receiveInput: input can not be null.");
        }

        if (!myInputs.contains(input)) {
            myInputs.add(input);
        }

        input.addListener(this);

    }

    /**
     * GUI Panel for this plugin.
     *
     * @return panel
     */
    @Override
    public JPanel getPanel() {
        return null;
    }

    /**
     * If interactive = true, the plugin will create dialogs and panels to interacts with the user
     *
     * @return boolean
     */
    @Override
    public boolean isInteractive() {
        return myIsInteractive;
    }

    /**
     * Parent Frame for this plugin. Can be null.
     *
     * @return frame
     */
    @Override
    public Frame getParentFrame() {
        return myParentFrame;
    }

    /**
     * Adds listener to this plugin.
     *
     * @param listener listener to add
     */
    @Override
    public void addListener(PluginListener listener) {

        synchronized (myListeners) {
            if ((listener != null) && (!myListeners.contains(listener))) {
                myListeners.add(listener);
            }
        }

    }

    public List<PluginListener> getListeners() {
        return Collections.unmodifiableList(myListeners);
    }

    public boolean hasListeners() {
        return !myListeners.isEmpty();
    }

    public List<Plugin> getInputs() {
        return myInputs;
    }

    /**
     * Returns data set after complete.
     *
     * @param event event
     */
    protected void fireDataSetReturned(PluginEvent event) {

        synchronized (myListeners) {
            Iterator<PluginListener> itr = myListeners.iterator();
            while (itr.hasNext()) {
                try {
                    if (myThreaded) {
                        PluginListener current = itr.next();
                        ThreadedPluginListener thread = new ThreadedPluginListener(current, event);
                        thread.start();
                    } else {
                        PluginListener current = itr.next();
                        current.dataSetReturned(event);
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

    }

    /**
     * Returns data set after complete.
     *
     * @param data data set
     */
    protected void fireDataSetReturned(DataSet data) {
        fireDataSetReturned(new PluginEvent(data));
    }

    private static final List<String> myPrintedCitations = new ArrayList<>();

    /**
     * Returns progress of execution.
     *
     * @param event event
     */
    protected void fireProgress(PluginEvent event) {

        synchronized (myListeners) {
            Iterator<PluginListener> itr = myListeners.iterator();
            while (itr.hasNext()) {
                PluginListener current = itr.next();
                current.progress(event);
            }
        }

        DataSet ds = (DataSet) event.getSource();
        if (ds != null) {
            List<Datum> percentage = ds.getDataOfType(Integer.class);

            if (percentage.size() > 0) {
                Datum datum = percentage.get(0);
                Integer percent = (Integer) datum.getData();
                if (percent == 100) {
                    if (!myPrintedCitations.contains(getCitation())) {
                        myLogger.info(getClass().getName() + "  Citation: " + getCitation());
                        myPrintedCitations.add(getCitation());
                    }
                }
            }
        }

    }

    /**
     * Returns progress of execution.
     *
     * @param percent percentage between 0 and 100 inclusive.
     */
    protected void fireProgress(Integer percent) {

        if ((percent < 0) || (percent > 100)) {
            throw new IllegalArgumentException("AbstractPlugin: fireProgress: percent must be between 0 and 100 inclusive.  arg: " + percent);
        }

        Datum percentage = new Datum("Percent", percent, null);
        fireProgress(new PluginEvent(new DataSet(percentage, this)));

    }

    @Override
    public String getCitation() {
        return DEFAULT_CITATION;
    }

    @Override
    public String pluginDescription() {
        return null;
    }

    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual";
    }

    //
    // Methods for PluginListener.
    //

    /**
     * Returns data set after complete.
     *
     * @param event event
     */
    @Override
    public void dataSetReturned(PluginEvent event) {

        DataSet input = (DataSet) event.getSource();

        performFunction(input);

    }

    /**
     * No operation for this abstract class.
     */
    @Override
    public void progress(PluginEvent event) {
        // The default action of a plugin is to do
        // nothing when another plugin reports its
        // progress.   This is intended to be implemented
        // by GUI applications to show the user the
        // progress of an interactive action.
    }

    public void reverseTrace(int indent) {

        if (myTrace) {
            return;
        }

        indent(indent);
        System.out.println(getClass().getName());

        Iterator<Plugin> itr = myInputs.iterator();
        while (itr.hasNext()) {
            try {
                AbstractPlugin current = (AbstractPlugin) itr.next();
                current.reverseTrace(indent + 3);
            } catch (Exception e) {
                // do nothing
            }
        }

        myTrace = true;

    }

    public void trace(int indent) {

        if (myTrace) {
            return;
        }

        indent(indent);
        System.out.println(getClass().getName());

        Iterator<PluginListener> itr = myListeners.iterator();
        while (itr.hasNext()) {
            try {
                AbstractPlugin current = (AbstractPlugin) itr.next();
                current.trace(indent + 3);
            } catch (Exception e) {
                // do nothing
            }
        }

        myTrace = true;

    }

    private void indent(int indent) {

        for (int i = 0; i < indent; i++) {
            System.out.print(" ");
        }

    }

    @Override
    public void setThreaded(boolean threaded) {
        myThreaded = threaded;
    }

    @Override
    public boolean cancel() {
        return false;
    }

    @Override
    public void run() {
        performFunction(null);
    }

    @Override
    public void progress(int percent, Object meta) {
        fireProgress(percent);
    }

    @Override
    public boolean wasCancelled() {
        return myWasCancelled;
    }

}
