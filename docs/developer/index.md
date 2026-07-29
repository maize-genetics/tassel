# Developer Guide

This guide is for developers who want to build TASSEL from source, understand
its internal architecture, add new analysis functionality, or contribute code
back to the project.

TASSEL is written primarily in **Java** (with some **Kotlin**) and built with
**Gradle**. The source lives on GitHub at
[maize-genetics/tassel](https://github.com/maize-genetics/tassel).

## Where to start

| If you want to…                                            | Read                                            |
| ---------------------------------------------------------- | ----------------------------------------------- |
| Compile the code and produce a runnable `sTASSEL.jar`      | [Building from Source](building-from-source.md) |
| Understand how the code is organized                       | [Project Architecture](architecture.md)         |
| Add a new analysis, filter, or data operation             | [Developing Plugins](plugin-development.md)      |
| Run the test suite and check statistical correctness       | [Testing](testing.md)                           |
| Submit a change through a pull request                     | [Contributing](../CONTRIBUTING.md)              |
| Cut a release and publish artifacts                        | [Releasing](releasing.md)                       |

!!! note "Branching model"
    TASSEL uses a `develop` integration branch: normal work flows
    `feature/* -> develop -> main`, and `main` holds released code. Branch off
    `develop` for your change (see [Contributing](../CONTRIBUTING.md)); releases
    are promoted from `develop` to `main` (see [Releasing](releasing.md)).

## Toolchain at a glance

| Tool          | Version / Notes                                                        |
| ------------- | --------------------------------------------------------------------- |
| Java (JDK)    | **21** — the build compiles to Java 21 bytecode via a Gradle toolchain |
| Kotlin        | 2.1.x (configured by the Gradle Kotlin plugin)                        |
| Build system  | **Gradle** — use the bundled wrapper (`./gradlew`); no local install  |
| Native BLAS   | **OpenBLAS** (optional, recommended) for fast matrix operations       |
| VCS           | **Git** / GitHub                                                       |
| Recommended IDE | IntelliJ IDEA (imports the Gradle project directly)                 |

## Quick start

```bash
# Clone the repository
git clone https://github.com/maize-genetics/tassel.git
cd tassel

# Build the runnable JAR (output in build/libs/)
./gradlew build

# Launch the graphical interface
./gradlew run
```

See [Building from Source](building-from-source.md) for prerequisites (including
OpenBLAS setup) and more detail.
