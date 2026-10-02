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

## The GBSv2 suite

The GBSv2 tests self-generate deterministic inputs through `GBSSimData`, so they
run as part of the broad `test` task without any extra setup. Two additional
tasks exercise the same classes against the real Chr9 datasets in `dataFiles/`:

```bash
./gradlew gbsTestSmall   # Chr9_10-200000, ~200 KB — fast, use while iterating
./gradlew gbsTestLarge   # Chr9_10-20000000, ~20 MB — slow, full validation
```

`gbsTest` is a back-compat alias for `gbsTestLarge`. The dataset reaches the
test JVM as `-Dgbs.test.dataset` and is read by
`GBSConstants.RAW_SEQ_CURRENT_TEST`, which must be set before that class loads.
A handful of tests need inputs that only ship with the 20 MB dataset and
self-skip on the small one through JUnit `Assume` guards.

Both tasks run with `ignoreFailures = true`, so check the reports under
`build/reports/tests/` rather than relying on the build's exit status.

Background on how these tests were rehabilitated is in
`docs/llm_notes/gbs-tests/` (internal notes, excluded from the docs site).

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

The GitHub Actions workflow (`.github/workflows/coverage.yml`) runs on Ubuntu
with JDK 21 and OpenBLAS installed. It has three jobs:

1. **Detect source changes** — classifies the pull request with
   `dorny/paths-filter`. The two heavy jobs below `needs:` this one and skip
   when nothing under `src/**` changed. They are skipped via `if:` rather than a
   workflow-level `paths:` filter, because a job skipped through `if:` reports
   success to branch protection while a workflow skipped by `paths:` stays
   pending forever and would block docs-only pull requests.
2. **Statistics gate (required)** — installs OpenBLAS, downloads the test data,
   and runs `./gradlew statisticsTest`. Must be green to merge.
3. **TASSEL 5 CI** — runs
   `./gradlew clean test koverXmlReport koverVerify --continue` and uploads
   coverage to Codecov. The suite itself is non-blocking
   (`ignoreFailures = true`), so this job gates on `koverVerify`.

!!! note "Release builds skip the gate"
    `statisticsTest` is wired into `check`, and `build` depends on `check`, so
    the packaging steps in `jdeploy.yml` and `nightly.yml` pass
    `-x statisticsTest` alongside the `-x test` and `-x koverVerify` they
    already used. Those jobs only package; the gate is enforced on the pull
    request and on pushes to `develop`, before anything is promoted to `main`.

## Writing tests

- Place tests in `src/test/java/net/maizegenetics/`, mirroring the package of the
  code under test.
- For statistical code, prefer asserting against known-good expected values
  (e.g. R-validated results stored in the test-data archive) rather than
  re-deriving them in the test.
- If a new statistical test should be enforced, add its fully-qualified class
  name to the `statisticsClasses` list in `build.gradle.kts` so it becomes part
  of the required gate.
