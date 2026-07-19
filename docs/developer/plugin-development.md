# Developing Plugins

Almost everything a user can do in TASSEL — loading data, filtering, running an
analysis, exporting results — is implemented as a **plugin**. Plugins are
*self-describing*: a single set of parameter declarations drives both the
graphical dialog and the command-line flags, so you write your logic once and get
both front-ends for free.

This page walks through writing a plugin from scratch. For a real, complete
example, read
[`KinshipPlugin`](https://github.com/maize-genetics/tassel/blob/main/src/main/java/net/maizegenetics/analysis/distance/KinshipPlugin.java).

## Anatomy of a plugin

Every plugin should:

1. **Extend `AbstractPlugin`** (`net.maizegenetics.plugindef.AbstractPlugin`).
2. **Provide the standard constructor:**

    ```java
    public MyPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    ```

3. **Declare its inputs as `PluginParameter` fields** (see below).
4. **Implement `processData(DataSet input)`** to do the work:

    ```java
    @Override
    public DataSet processData(DataSet input) {
        // ...do whatever your plugin does...
        return result; // may be null
    }
    ```

5. **Provide the GUI hooks** so the plugin appears correctly in the interface.
   Give `getButtonName()` and `getToolTipText()` meaningful values:

    ```java
    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Example Function";
    }

    @Override
    public String getToolTipText() {
        return "Performs Function by Algorithm";
    }
    ```

A plugin **may** also override:

- `preProcessParameters(DataSet input)` — validation/setup that must run *before*
  the user is prompted (for example, checking that a genotype table is selected).
- `postProcessParameters()` — logic that runs *after* parameters are set.
- `pluginDescription()` — a human-readable description.
- `getCitation()` — the citation users should reference.
- `pluginUserManualURL()` — a link to the relevant user-manual page.

A plugin **should not**:

- Call `System.exit()`.
- Implement `performFunction()` (that belongs to `AbstractPlugin`).
- Handle errors with dialogs or manual logging. Instead, **throw an exception**
  and let `AbstractPlugin` present it appropriately for GUI or CLI:

    ```java
    // Don't do this:
    if (alignInList.size() != 1) {
        String msg = "Invalid selection. Please select one genotype alignment.";
        if (isInteractive()) {
            JOptionPane.showMessageDialog(getParentFrame(), msg);
        } else {
            myLogger.error(msg);
        }
        return null;
    }

    // Do this instead:
    if (alignInList.size() != 1) {
        throw new IllegalArgumentException(
            "Invalid selection. Please select one genotype alignment.");
    }
    ```

- Keep a `main()` method in committed code. A temporary `main()` is used to
  auto-generate getters/setters (see [below](#generating-getters-and-setters)),
  but comment it out or remove it before committing.

## Declaring parameters

Declare each parameter as a private `PluginParameter` field using
`PluginParameter.Builder`. The declaration order determines the order fields
appear in the GUI dialog, so group related parameters together.

The builder constructor takes three arguments:

- **command-line name** — no spaces; use `camelCase` for multi-word names.
- **default value** — or `null` for no default.
- **class type** — e.g. `String.class`, `Double.class`, an enum class.

```java
private PluginParameter<String> inputFile =
    new PluginParameter.Builder<>("inputFile", null, String.class)
        .inFile()
        .required(true)
        .description("The genotype file to read.")
        .build();
```

Always finish the chain with `.build()`.

### Common builder methods

| Method | Effect |
| ------ | ------ |
| `.description("…")` | Short description of the parameter. |
| `.guiName("…")` | Override the GUI label (defaults to a title-cased version of the CLI name, e.g. `inputFile` → "Input File"). |
| `.inFile()` / `.outFile()` | Parameter is an input/output **file** path. |
| `.inDir()` / `.outDir()` | Parameter is an input/output **directory** path. |
| `.required(true)` | Mark as required. A required parameter cannot also have a default value. |
| `.range(Range.closed(0.0, 1.0))` | Restrict to an inclusive numeric range (Guava `Range`). |
| `.range(MyEnum.values())` | Restrict to a set of enum values. |
| `.units("centimorgans")` | Document the parameter's units. |
| `.dependentOnParameter(other)` | Enable this parameter only when a prior boolean parameter is `true`. |
| `.dependentOnParameter(other, value)` | Enable only when `other` equals `value`. |
| `.dependentOnParameter(other, new Object[]{a, b})` | Enable when `other` is one of the listed values. |
| `.genotypeTable()` | Let the user pick which component of a genotype table to use. |
| `.distanceMatrix()` | Let the user pick among selected distance matrices. |
| `.taxaNameList()` / `.siteNameList()` | Searchable taxa/site name selection. |
| `.positionList()` | Accept a `.json.gz` file, imported as a `PositionList`. |

### Getter/setter methods

For each parameter, provide a paired accessor: a *getter* that returns the value
and a *setter* that takes a value and returns the plugin (for fluent chaining).
By convention these methods are named after the parameter (no literal `get`/`set`
prefix):

```java
// getter
public String inputFile() {
    return inputFile.value();
}

// setter
public MyPlugin inputFile(String value) {
    inputFile = new PluginParameter<>(inputFile, value);
    return this;
}
```

These wrap the inherited `getParameter(...)` / `setParameter(...)`, but providing
them makes your plugin much easier to call from other code.

#### Generating getters and setters

The `GeneratePluginCode` class writes all the accessors (with Javadoc) for you.
Temporarily add a `main` method, run it, copy the console output into your class,
then remove the `main`:

```java
public static void main(String[] args) {
    GeneratePluginCode.generate(MyPlugin.class);
}
```

## Data flow: `DataSet` and `Datum`

Plugins exchange data as `DataSet` objects. A `DataSet` is a collection of
`Datum` items, each wrapping a typed payload (such as a `GenotypeTable`,
`Phenotype`, or `DistanceMatrix`) along with a name and comment.

Inside `processData`, pull out the inputs you need by type:

```java
List<Datum> genotypes = input.getDataOfType(GenotypeTable.class);
if (genotypes == null || genotypes.isEmpty()) {
    throw new IllegalArgumentException(
        "MyPlugin: Nothing selected. Please select a genotype.");
}
```

Build results as new `Datum` objects and return them in a `DataSet`:

```java
Datum result = new Datum(resultName, resultObject, comment);
return new DataSet(result, this);
```

## Logging

Use Log4j 2 for logging:

```java
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

private static final Logger myLogger = LogManager.getLogger(MyPlugin.class);
```

## A convenience `runPlugin` method

It is common to add a typed convenience method so callers can run the plugin and
get a single result object back directly:

```java
public DistanceMatrix runPlugin(DataSet input) {
    return (DistanceMatrix) performFunction(input).getData(0).getData();
}
```

## Wiring a plugin into the pipeline and GUI

Because the plugin declares its own parameters, the command-line
[pipeline](../pipelines/tassel5-pipeline-cli.md) and the GUI can both drive it
without any per-plugin parsing code. Once your plugin is on the classpath, its
parameters are exposed as CLI flags and as dialog fields automatically.
