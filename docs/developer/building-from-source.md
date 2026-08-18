# Building from Source

TASSEL is built with [Gradle](https://gradle.org/). A Gradle **wrapper**
(`./gradlew`) is committed to the repository, so you do **not** need a local
Gradle installation — the wrapper downloads the correct version automatically.

## Prerequisites

| Requirement | Details |
| ----------- | ------- |
| **Java 21 (JDK)** | The build uses a Gradle toolchain targeting Java 21 and compiles to Java 21 bytecode. Verify with `java -version`. |
| **Git** | Needed to clone the repository. Verify with `git --version`. |
| **OpenBLAS** (optional, recommended) | A native BLAS library that accelerates matrix operations and is required by some statistical tests. |

### Installing OpenBLAS

Some analyses (kinship, MLM, PCA, and the statistical test suite) call into a
native BLAS library through JNI. Installing OpenBLAS is optional for a plain
build, but recommended for full functionality and required to run the
statistics test gate.

=== "macOS"

    ```bash
    brew install openblas
    ```

=== "Linux (Debian/Ubuntu)"

    ```bash
    sudo apt-get update
    sudo apt-get install -y libopenblas-dev
    ```

The build and test tasks auto-detect common install locations
(`/opt/homebrew/opt/openblas/lib` and `/usr/local/opt/openblas/lib` on macOS,
`/usr/lib/x86_64-linux-gnu` on Linux). If your library lives elsewhere, point at
it with the `BLAS_LIB_PATH` environment variable:

```bash
export BLAS_LIB_PATH=/path/to/dir/containing/libopenblas
```

## Getting the source

```bash
git clone https://github.com/maize-genetics/tassel.git
cd tassel
```

To update an existing checkout to the latest code:

```bash
git pull
```

## Building

Build the project and assemble the runnable JAR:

```bash
./gradlew build
```

This produces:

- `build/libs/sTASSEL.jar` — the TASSEL application JAR. Its manifest declares
  the main class and a `Class-Path` that points at `lib/`.
- `build/libs/lib/` — all runtime dependencies, copied next to the JAR so the
  manifest classpath resolves.

A self-contained JAR that needs no `lib/` directory is available too, but it is
built only on demand because it is over 70 MB:

```bash
./gradlew shadowJar   # build/libs/tassel-<version>-jar-with-dependencies.jar
```

To do a clean rebuild:

```bash
./gradlew clean build
```

!!! tip "Skipping tests during a build"
    The full `build` runs the enforced `statisticsTest` gate (wired into
    `check`). To build the JAR without running tests — for example, when you
    just need an artifact quickly — skip the test tasks:

    ```bash
    ./gradlew clean build -x test -x statisticsTest
    ```

## Running

### Graphical interface

Run through Gradle:

```bash
./gradlew run
```

Or launch the built JAR directly (the manifest supplies the main class and
classpath):

```bash
java -jar build/libs/sTASSEL.jar
```

### Command-line pipeline

TASSEL analyses can be scripted through the pipeline entry point,
`net.maizegenetics.pipeline.TasselPipeline`:

```bash
java -classpath 'build/libs/sTASSEL.jar:build/libs/lib/*' \
  net.maizegenetics.pipeline.TasselPipeline -h
```

Increase the heap for large datasets with a `-Xmx` flag, e.g. `-Xmx10g`. For
the full pipeline command language, see the
[TASSEL 5 Pipeline (CLI)](../pipelines/tassel5-pipeline-cli.md) documentation.

## Useful Gradle tasks

| Task | Purpose |
| ---- | ------- |
| `./gradlew build` | Compile, test, and assemble `sTASSEL.jar` + `lib/`. |
| `./gradlew run` | Launch the Swing GUI. |
| `./gradlew test` | Run the full (non-blocking) test suite. |
| `./gradlew statisticsTest` | Run the enforced statistical-correctness gate. |
| `./gradlew koverHtmlReport` | Generate an HTML code-coverage report. |
| `./gradlew dokkaGenerate` | Generate API documentation with Dokka. |
| `./gradlew shadowJar` | Assemble the self-contained `jar-with-dependencies` JAR. |
| `./gradlew fetchTestData` | Download the shared test-data archive into `dataFiles/`. |
| `./gradlew printVersion` | Print the current project version. |

Run `./gradlew tasks` to see the full list.

## Common issues

- **`UnsatisfiedLinkError` / BLAS failures** — OpenBLAS is not installed or not
  found. Install it (see above) or set `BLAS_LIB_PATH`.
- **Wrong Java version** — the toolchain targets Java 21. Make sure a JDK 21 is
  installed and discoverable; Gradle can auto-provision a toolchain, but a local
  JDK 21 avoids surprises.
- **Out-of-memory during large analyses** — increase the JVM heap with `-Xmx`.
