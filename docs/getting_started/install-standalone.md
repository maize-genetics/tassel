---
hide:
  - navigation
---

# Install the Standalone Distribution

The **standalone** distribution is a self-contained bundle of the TASSEL JAR, its
libraries, and helper launch scripts. It is ideal for running TASSEL from the
command line, on a server, or inside an analysis pipeline.

Unlike the native GUI installers, the standalone distribution does **not** bundle
a Java runtime, so you need Java installed yourself.

## Prerequisites

- **Java 21** (JDK or JRE) on your `PATH`. Verify with:

    ```bash
    java -version
    ```

- **Perl** to use the `.pl` launch scripts (pre-installed on macOS and most
  Linux distributions). On Windows, use the provided `.bat` scripts instead.
- *(Optional, recommended)* **OpenBLAS** for faster native matrix operations.

## Download

Grab the latest standalone archive from the [Download page](../download/index.md)
(select **Standalone**), or directly:

- [`tassel-5-standalone-v5.2.97.tar.gz`](https://github.com/maize-genetics/tassel/releases/download/v5.2.97/tassel-5-standalone-v5.2.97.tar.gz)

A `.zip` archive is also published on the
[releases page](https://github.com/maize-genetics/tassel/releases) for each
version.

## Extract

=== "macOS / Linux"

    ```bash
    tar -xzf tassel-5-standalone-v5.2.97.tar.gz
    cd tassel-5-standalone
    ```

=== "Windows"

    Extract the `.zip` with File Explorer (right-click &rarr; **Extract All**),
    then open a Command Prompt in the extracted `tassel-5-standalone` folder.

## Run

The distribution ships two entry points:

- **`start_tassel`** &mdash; launches the graphical interface.
- **`run_pipeline`** &mdash; runs the command-line pipeline.

=== "macOS / Linux"

    ```bash
    # Launch the GUI
    ./start_tassel.pl

    # Run a pipeline (example: print the help/usage)
    ./run_pipeline.pl
    ```

=== "Windows"

    ```bat
    :: Launch the GUI
    start_tassel.bat

    :: Run a pipeline
    run_pipeline.bat
    ```

### Allocating memory

The launch scripts accept a maximum heap size. For large datasets, increase it,
for example to 20&nbsp;GB:

```bash
./run_pipeline.pl -Xmx20g ...
```

## Next steps

- Learn the command-line workflow in the [Pipelines](../pipelines/index.md)
  section and the
  [Pipeline Tutorial](../pipelines/pipeline-tutorial.md).
- Prefer a point-and-click app? See the [GUI installation guide](install-gui.md).
