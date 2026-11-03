# Building from Source

TASSEL is built with [Gradle](https://gradle.org/). A Gradle **wrapper**
(`./gradlew`) is committed to the repository, so you do **not** need a local
Gradle installation — the wrapper downloads the correct version automatically.

## Prerequisites

| Requirement | Details |
| ----------- | ------- |
| **Java 21 (JDK)** | The build uses a Gradle toolchain targeting Java 21 and compiles to Java 21 bytecode. Verify with `java -version`. |
| **Git** | Needed to clone the repository. Verify with `git --version`. |
| **C compiler and BLAS** (optional) | Only needed to build the native BLAS backend described below. |

### Native BLAS backend

TASSEL's matrix code runs on [EJML](https://ejml.org/), a pure-Java library,
unless the optional **TasselBlas** JNI library is on `java.library.path`.
TasselBlas calls the platform's BLAS and LAPACK instead. You don't need it to
build, run, or test TASSEL; it's needed only for the native BLAS variant of the
statistics gate.

The `buildTasselBlas` task compiles it into `build/native/tasselblas/` with the
system C compiler (`cc`):

=== "macOS"

    Links Apple's Accelerate framework, which ships with macOS. If `cc` is
    missing, install the Xcode command-line tools:

    ```bash
    xcode-select --install
    ```

=== "Linux (Debian/Ubuntu)"

    Links OpenBLAS:

    ```bash
    sudo apt-get update
    sudo apt-get install -y build-essential libopenblas-dev
    ```

Windows is not supported.

```bash
./gradlew buildTasselBlas       # compile the library only
./gradlew statisticsTestBlas    # statistics gate on BLAS (builds the library first)
./gradlew test -PtasselBlas     # full suite, with DoubleMatrixTest's BLAS case required
```

When TasselBlas loads, `DoubleMatrixFactory` uses BLAS by default. To pick the
backend explicitly, set the `tassel.matrix.factory` system property to `ejml`,
`colt`, or `blas`. To use the library outside Gradle, pass its directory to the
JVM:

```bash
java -Djava.library.path=build/native/tasselblas -jar build/libs/sTASSEL.jar
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
| `./gradlew statisticsTestBlas` | Run the same gate on the native BLAS backend. |
| `./gradlew buildTasselBlas` | Compile the optional TasselBlas JNI library. |
| `./gradlew koverHtmlReport` | Generate an HTML code-coverage report. |
| `./gradlew dokkaGenerate` | Generate API documentation with Dokka. |
| `./gradlew shadowJar` | Assemble the self-contained `jar-with-dependencies` JAR. |
| `./gradlew fetchTestData` | Download the shared test-data archive into `dataFiles/`. |
| `./gradlew printVersion` | Print the current project version. |

Run `./gradlew tasks` to see the full list.

## Common issues

- **`buildTasselBlas` can't find `cblas.h` or `-lopenblas`** — on Linux,
  install `libopenblas-dev` (see [Native BLAS backend](#native-blas-backend)).
- **`TasselBlas failed to load` in `DoubleMatrixTest`** — the library was built
  but the JVM couldn't load it. The assertion message carries the underlying
  `UnsatisfiedLinkError`, which names the missing file or symbol.
- **Wrong Java version** — the toolchain targets Java 21. Make sure a JDK 21 is
  installed and discoverable; Gradle can auto-provision a toolchain, but a local
  JDK 21 avoids surprises.
- **Out-of-memory during large analyses** — increase the JVM heap with `-Xmx`.
