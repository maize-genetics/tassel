# Testing

TASSEL uses [JUnit](https://junit.org/) tests run through Gradle. Tests live
under `src/test/java/net/maizegenetics/`, mirroring the main source layout.

There are two test entry points with very different roles:

| Task | Blocking? | Purpose |
| ---- | --------- | ------- |
| `./gradlew statisticsTest` | **Yes** (CI gate) | Verifies TASSEL's numeric results against R-validated expected values. |
| `./gradlew test` | No | Runs the full suite, including pipeline/IO tests that are still being stabilized. |

## Fetching test data

Many tests read from a shared test-data archive that is downloaded into the
**git-ignored `dataFiles/`** directory. Fetch it once after a clean checkout:

```bash
./gradlew fetchTestData
```

This downloads and extracts the `tassel_test_data` release archive from
[maize-genetics/tassel_test_data](https://github.com/maize-genetics/tassel_test_data).
The task is a no-op if `dataFiles/` already exists and is non-empty.

!!! note "OpenBLAS required for statistics tests"
    The statistical tests exercise native BLAS routines. Install OpenBLAS (see
    [Building from Source](building-from-source.md#installing-openblas)) or set
    `BLAS_LIB_PATH` before running them.

## The statistics gate (required)

`statisticsTest` is the **enforced CI gate**. It runs only the classes that
verify TASSEL's statistical correctness — kinship, MLM, GLM, PCA, linkage
disequilibrium, distance matrices, linear models, and related numeric transforms
— with `ignoreFailures = false`, so any failure is visible and blocks a merge.

```bash
./gradlew statisticsTest
```

It is wired into `check`, so `./gradlew check` also enforces it. Reports are
written to:

- `build/reports/tests/statisticsTest/` (HTML)
- `build/test-results/statisticsTest/` (JUnit XML)

If you change any analysis or statistics code, run this gate locally before
opening a pull request.

## The full suite (non-blocking)

```bash
./gradlew test
```

The broad `test` task runs everything but is currently **non-blocking**
(`ignoreFailures = true`) while some pipeline and I/O tests are being fixed, and
it excludes a set of environment-sensitive tests (certain GBS, HDF5, and
hard-coded-path tests). Failures here are informative but do not block CI.

## Coverage

Coverage is measured with [Kover](https://github.com/Kotlin/kotlinx-kover) using
the JaCoCo engine, focused on **branch** coverage of the analysis and pipeline
logic. Pure GUI/Swing code is excluded from the report so coverage reflects
exercised analytical logic.

```bash
# Human-readable HTML report
./gradlew koverHtmlReport

# XML report (used by CI / Codecov)
./gradlew koverXmlReport
```

The HTML report is written under `build/reports/kover/`.

!!! tip "Excluding generated GUI boilerplate"
    Auto-generated plugin accessors and GUI hook methods are annotated with
    `@GeneratedGuiBoilerplate`. Because its name contains "Generated", JaCoCo
    automatically drops those methods from coverage.

## What CI runs

The GitHub Actions workflow (`.github/workflows/coverage.yml`) runs on pull
requests that touch `src/**`, on Ubuntu with JDK 21 and OpenBLAS installed. It
has two jobs:

1. **Statistics gate (required)** — installs OpenBLAS, downloads the test data,
   and runs `./gradlew statisticsTest`. Must be green to merge.
2. **Full suite & coverage (non-blocking)** — runs
   `./gradlew test koverXmlReport` and uploads coverage to Codecov.

## Writing tests

- Place tests in `src/test/java/net/maizegenetics/`, mirroring the package of the
  code under test.
- For statistical code, prefer asserting against known-good expected values
  (e.g. R-validated results stored in the test-data archive) rather than
  re-deriving them in the test.
- If a new statistical test should be enforced, add its fully-qualified class
  name to the `statisticsClasses` list in `build.gradle.kts` so it becomes part
  of the required gate.
