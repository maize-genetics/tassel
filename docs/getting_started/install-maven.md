---
hide:
  - navigation
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
| **Latest version** | `5.2.97` |

## Add the dependency

=== "Gradle (Kotlin DSL)"

    ```kotlin
    // build.gradle.kts
    repositories {
        mavenCentral()
    }

    dependencies {
        implementation("net.maizegenetics:tassel:5.2.97")
    }
    ```

=== "Gradle (Groovy DSL)"

    ```groovy
    // build.gradle
    repositories {
        mavenCentral()
    }

    dependencies {
        implementation 'net.maizegenetics:tassel:5.2.97'
    }
    ```

=== "Maven"

    ```xml
    <!-- pom.xml -->
    <dependency>
      <groupId>net.maizegenetics</groupId>
      <artifactId>tassel</artifactId>
      <version>5.2.97</version>
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
