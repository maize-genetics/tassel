---
hide:
  - navigation
render_macros: true
---

# Use TASSEL as a Library (Maven / Gradle)

TASSEL is published to **Maven Central**, so you can add it to any JVM project
(Java or Kotlin) and call its APIs directly instead of shelling out to the
command line.

## Coordinates

| Coordinate | Value |
| --- | --- |
| **Group** | `net.maizegenetics` |
| **Artifact** | `tassel` |
| **Latest version** | `{{ version }}` |

Each release publishes four artifacts:

| Artifact | Classifier | Use |
| --- | --- | --- |
| `tassel-<version>.jar` | *(none)* | The library. Dependencies resolve through the POM. |
| `tassel-<version>-jar-with-dependencies.jar` | `jar-with-dependencies` | Self-contained "fat" JAR with every dependency bundled, for environments without dependency resolution. |
| `tassel-<version>-sources.jar` | `sources` | Sources, for IDE navigation and debugging. |
| `tassel-<version>-javadoc.jar` | `javadoc` | API documentation. |

!!! warning "Skip 5.2.97"
    The 5.2.97 artifacts on Maven Central are defective: the fat JAR is missing, and
    the sources and javadoc JARs are copies of the main JAR. Maven Central is
    immutable, so 5.2.97 cannot be repaired — use 5.2.98 or later.

## Add the dependency

=== "Gradle (Kotlin DSL)"

    ```kotlin
    // build.gradle.kts
    repositories {
        mavenCentral()
    }

    dependencies {
        implementation("net.maizegenetics:tassel:{{ version }}")
    }
    ```

=== "Gradle (Groovy DSL)"

    ```groovy
    // build.gradle
    repositories {
        mavenCentral()
    }

    dependencies {
        implementation 'net.maizegenetics:tassel:{{ version }}'
    }
    ```

=== "Maven"

    ```xml
    <!-- pom.xml -->
    <dependency>
      <groupId>net.maizegenetics</groupId>
      <artifactId>tassel</artifactId>
      <version>{{ version }}</version>
    </dependency>
    ```

!!! tip "Always up to date"
    Check the latest published version on
    [Maven Central](https://central.sonatype.com/artifact/net.maizegenetics/tassel)
    or the [Download page](../download/index.md).

## Requirements

- **Java 21** or newer. TASSEL is compiled to Java 21 bytecode.
- *(Optional, recommended)* **OpenBLAS** installed on the host for fast native
  matrix operations used by some analyses.

## Verify your setup

After syncing dependencies, confirm TASSEL classes resolve. For example, print
the version from the `net.maizegenetics` package on your classpath, or run a
minimal program that imports a TASSEL class such as
`net.maizegenetics.tassel.TASSELMainFrame`.

## Next steps

- Building TASSEL yourself instead of consuming the release? See
  [Building from Source](../developer/building-from-source.md).
- Extending TASSEL with new analyses? See
  [Developing Plugins](../developer/plugin-development.md).
