# Executing TASSEL

This page describes how to launch TASSEL, adjust the Java heap size for the GUI
and command line, and troubleshoot common startup errors. For the full
command-line flag reference, see the
[Pipeline Reference](../pipelines/tassel5-pipeline-cli.md) and the
[GBSv2 Pipeline](../gbsv2_pipeline/index.md).

## Change Heap Size (GUI)

TASSEL runs on the Java Virtual Machine, which reserves a fixed range of memory
called the heap. Increase the maximum heap (`-Xmx`) if you run out of memory on
large datasets.

### macOS

1. In the Applications folder, double-click `TASSEL 5`.
2. Right-click the `TASSEL 5` application and choose **Show Package Contents**.
3. Open the `Contents` folder.
4. Edit the file `vmoptions.txt`.
5. Change `-Xms` (minimum heap, e.g. `-Xms512m`) and `-Xmx` (maximum heap, e.g. `-Xmx5g`) as appropriate.

### Windows

1. Go to the TASSEL installation directory (default: `C:\Program Files\TASSEL5`).
2. Edit the file `Tassel 5.vmoptions`.
3. Change `-Xms` (minimum heap, e.g. `-Xms512m`) and `-Xmx` (maximum heap, e.g. `-Xmx5g`) as appropriate.

!!! note

    If User Account Control (UAC) is enabled, start your text editor with
    **Run as administrator** before editing `Tassel 5.vmoptions`.

## Command Line (macOS / Linux / Windows Bash)

Add the TASSEL directory to your `PATH` so you can run the scripts from any
location, for example in `.bash_profile`:

```bash
PATH=$PATH:/analysis_tools/tassel5.0_standalone
```

Then use these commands:

```bash
./start_tassel.pl -Xmx4g   # Start the GUI with a 4 GB max heap size
./run_pipeline.pl -Xmx5g   # Run the command line with a 5 GB max heap size
```

## Command Line (Windows)

```bat
start_tassel.bat   :: Start the GUI
run_pipeline.bat   :: Run the command line
```

!!! note

    To increase the heap size, edit the `.bat` file you are using and change the
    `-Xmx` parameter.

## Troubleshooting

!!! warning "Could not create the Java Virtual Machine"

    If you get an error message containing `Error: Could not create the Java
    Virtual Machine`, `The specified size exceeds the maximum representable
    size`, or `Could not reserve enough space for 2097152KB object heap`, you
    most likely have a mismatch between your Java installation and operating
    system: one is 32-bit and the other is 64-bit.

    Run `java -version` to check whether your Java matches your operating system.
    Note that there may be multiple Java installations on your machine, and your
    web browser may use a different one than your command line. Run
    `java -d64 -version` to see if your installation supports 64-bit.

    See the
    [Oracle HotSpot heap FAQ](https://www.oracle.com/technetwork/java/hotspotfaq-138619.html#gc_heap_32bit)
    for details.

!!! warning "Invalid maximum heap size"

    ```text
    Invalid maximum heap size: -Xmx4g
    The specified size exceeds the maximum representable size.
    Error: Could not create the Java Virtual Machine.
    ```

    This may be due to a limit on the maximum heap, which is system dependent.
    Try reducing the `-Xmx` setting.

!!! warning "UnsupportedClassVersionError"

    ```text
    Exception in thread "main" java.lang.UnsupportedClassVersionError:
    net/maizegenetics/pipeline/TasselPipeline : Unsupported major.minor version 51.0
    ```

    You need a higher version of Java. TASSEL 5 requires Java 1.8.

!!! warning "Could not find or load main class"

    If you get `Error: Could not find or load main class
    net.maizegenetics.tassel.TASSELMainApp` while running `start_tassel.pl` or
    `run_pipeline.pl` in the Bash shell on Windows, change the following line to
    use a `;` instead of a `:`.

    ```perl
    my $CP = join(":", @fl);
    ```

!!! note "macOS Gatekeeper"

    TASSEL does not have an Apple signature, so you may need to adjust Gatekeeper
    to install the program. See
    [Apple's instructions](https://support.apple.com/en-us/HT202491).
