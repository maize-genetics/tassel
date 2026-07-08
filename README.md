### TASSEL

* [Latest Builds](http://www.maizegenetics.net/tassel/)

---

## Development

### Running the test suite locally

**1. Fetch test data (one-time, after a clean checkout)**

Tests require genotype and expected-result files that are distributed separately
(they are too large to commit and are git-ignored via `dataFiles/`).
Download and unpack them with the bundled Gradle task:

```bash
./gradlew fetchTestData
```

This downloads the `tassel_test_data_v1` release archive from GitHub, extracts
it into `dataFiles/`, and deletes the archive.  The task is a no-op if
`dataFiles/` already contains files.

> **GBS tests self-generate their data.** The GBSv2 tests (e.g.
> `GBSSeqToTagDBPluginTest`, `DiscoverySNPCallerPluginV2Test`,
> `SNPQualityProfilerPluginTest`, `SNPCutPosTagVerificationPluginTest`,
> `GetTagSequenceFromDBPluginTest`, `GetTagTaxaDistFromDBTest`,
> `ParseBarcodeReadTest`) build their own deterministic reference FASTA, key
> file, FASTQ, and (aligner-free) SAM at runtime via
> `net.maizegenetics.analysis.gbs.v2.GBSSimData`. They need no `fetchTestData`
> download and no external aligner (BWA/bowtie/PEAR). A small set of legacy GBS
> tests remain excluded in `build.gradle.kts` (each with a one-line rationale):
> they require the jhdf5 native library, byte-exact golden fixtures, or
> hardcoded dev-machine paths.

**2. Run the statistical-correctness gate**

The statistics gate exercises kinship, MLM, GLM, PCA, LD, distance matrices,
and linear models against R-validated expected results.  It runs with
`ignoreFailures = false`, so failures are visible rather than silently swallowed:

```bash
./gradlew statisticsTest
# HTML report: build/reports/tests/statisticsTest/index.html
```

**3. Run the full test suite**

The broad suite includes all non-excluded tests with `ignoreFailures = true`
(failures are reported but do not fail the build while IO/pipeline tests are
being fixed):

```bash
./gradlew test
# HTML report: build/reports/tests/test/index.html
```

**4. OpenBLAS (needed for some numeric tests)**

Set `BLAS_LIB_PATH` if OpenBLAS is in a non-standard location:

```bash
BLAS_LIB_PATH=/opt/homebrew/opt/openblas/lib ./gradlew statisticsTest
```

On macOS with Homebrew, `brew install openblas` places the library at
`/opt/homebrew/opt/openblas/lib` (Apple Silicon) or
`/usr/local/opt/openblas/lib` (Intel), which the build script detects
automatically.