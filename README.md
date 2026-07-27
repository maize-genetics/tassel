### TASSEL

TASSEL (**T**rait **A**nalysis by a**SS**ociation, **E**volution and **L**inkage) is a
software package for evaluating trait associations, evolutionary patterns, and
linkage disequilibrium in genetic data. It is designed to handle the diversity of
data types and sizes common in modern genomics, including a variety of genotype and
phenotype file formats.

* [Latest Builds](http://www.maizegenetics.net/tassel/)

## Overview

TASSEL implements a broad range of statistical genetics and population genomics
methods, including:

* **Association analysis** - general linear model (GLM), mixed linear model (MLM),
  fast/multithreaded association, stepwise model fitting, and genomic selection.
* **Kinship and distance** - centered/normalized IBS, additive (A) and hybrid (H)
  relationship matrices, and multidimensional scaling.
* **Population genetics** - linkage disequilibrium, sequence diversity, and
  principal component analysis.
* **Data handling** - import, filter, merge, transform, and export genotype and
  phenotype data across many formats (HapMap, VCF, PLINK, HDF5, and more).

TASSEL can be used through its Swing-based graphical user interface or through a
command-line pipeline for scripted and batch workflows. Functionality is organized
into composable *plugins* that can be chained together into analysis pipelines.

## Requirements

* **Java 21** (the build targets and compiles against Java 21 bytecode).
* **Gradle** - a Gradle wrapper (`./gradlew`) is included, so a local Gradle
  installation is not required.
* **OpenBLAS** (optional, recommended) - a native BLAS library speeds up matrix
  operations and is required by some statistical tests. The build auto-detects
  common install locations, or you can set the `BLAS_LIB_PATH` environment variable
  to point at the directory containing the native library.
  * **macOS** - `brew install openblas`
  * **Linux (Debian/Ubuntu)** - `apt-get install libopenblas-dev`

## How to Run

### Build

Build the project and assemble the runnable JAR (`build/libs/sTASSEL.jar`) along
with its runtime dependencies (`build/libs/lib/`):

```bash
./gradlew build
```

### Launch the graphical interface

Run the application through Gradle:

```bash
./gradlew run
```

Or launch the built JAR directly (the manifest sets the main class and classpath):

```bash
java -jar build/libs/sTASSEL.jar
```

### Run the command-line pipeline

TASSEL analyses can be scripted through the pipeline entry point
(`net.maizegenetics.pipeline.TasselPipeline`), for example:

```bash
java -classpath 'build/libs/sTASSEL.jar:build/libs/lib/*' \
  net.maizegenetics.pipeline.TasselPipeline -h
```

Increase the heap for large datasets by adding a `-Xmx` flag, e.g. `-Xmx10g`.

## Standalone Version

Prebuilt standalone distributions of TASSEL are available from the
[Latest Builds](http://www.maizegenetics.net/tassel/) page. A standalone
distribution bundles `dist/sTASSEL.jar` together with all runtime dependencies in
`lib/`, and includes convenience launcher scripts in the same directory:

* `start_tassel.pl` - launches the graphical interface.
* `run_pipeline.pl` / `run_pipeline.bat` - runs the command-line pipeline
  (`net.maizegenetics.pipeline.TasselPipeline`) for scripting and batch analyses.

These scripts build the Java classpath automatically from `lib/` and
`dist/sTASSEL.jar`, so you do not need to construct it by hand.

### Using `run_pipeline` in scripts

`run_pipeline.pl` (macOS/Linux) and `run_pipeline.bat` (Windows) accept TASSEL
pipeline arguments directly. Any `-Xms`/`-Xmx` arguments are pulled out and passed
to the JVM, while all other arguments are forwarded to the pipeline. If no memory
flags are given, the scripts default to `-Xms512m -Xmx1536m`.

Show pipeline help:

```bash
./run_pipeline.pl
```

A typical analysis - load a HapMap genotype file and a phenotype file, run an MLM
association, and write the results - chains plugins together with `-fork`,
`-input`, and `-combine` directives:

```bash
./run_pipeline.pl -Xmx8g \
  -fork1 -h mydata.hmp.txt \
  -fork2 -r mytraits.txt \
  -combine3 -input1 -input2 -intersect \
  -mlm -export mlm_results
```

On Windows, invoke the batch launcher instead:

```bat
run_pipeline.bat -Xmx8g -fork1 -h mydata.hmp.txt -fork2 -r mytraits.txt -combine3 -input1 -input2 -intersect -mlm -export mlm_results
```

Because `run_pipeline` is a plain command-line invocation, it can be embedded in
shell scripts, cron jobs, or workflow managers to batch-process many datasets. For
example, to run the same pipeline over several genotype files:

```bash
#!/bin/bash
for geno in genotypes/*.hmp.txt; do
  name=$(basename "$geno" .hmp.txt)
  ./run_pipeline.pl -Xmx8g \
    -fork1 -h "$geno" \
    -fork2 -r mytraits.txt \
    -combine3 -input1 -input2 -intersect \
    -mlm -export "results/${name}_mlm"
done
```

## Development

### Project layout

* `src/main/java/net/maizegenetics/` - main source (Java and Kotlin).
* `src/test/java/net/maizegenetics/` - unit and integration tests.
* `build.gradle.kts` / `settings.gradle.kts` - Gradle build configuration.

### Test data

Many tests depend on a shared test-data archive that is downloaded into the
git-ignored `dataFiles/` directory. Fetch it once after a clean checkout:

```bash
./gradlew fetchTestData
```

### Running tests

The statistical-correctness suite is the enforced CI gate and is wired into
`check`:

```bash
./gradlew statisticsTest
```

Run the full (non-blocking) test suite:

```bash
./gradlew test
```

### Coverage

Coverage is measured through Kover (using the JaCoCo engine), focusing on branch
coverage of the analysis and pipeline logic while excluding GUI code:

```bash
./gradlew koverHtmlReport
```

### Documentation

The user and developer documentation lives in `docs/` and is published with
MkDocs. Preview it locally with:

```bash
pip install mkdocs-material
mkdocs serve
```

API documentation is generated with Dokka:

```bash
./gradlew dokkaHtml
```

### Contributing

Start from the branch that matches your change: `develop` for features and
non-urgent fixes, `main` for hotfixes and documentation-only changes.
Documentation changes take a shortcut — they skip the test suite and publish
without cutting a release. See
[Contributing](docs/CONTRIBUTING.md#the-git-workflow) for the full workflow.
