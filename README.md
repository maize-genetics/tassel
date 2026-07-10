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

API documentation is generated with Dokka:

```bash
./gradlew dokkaHtml
```
