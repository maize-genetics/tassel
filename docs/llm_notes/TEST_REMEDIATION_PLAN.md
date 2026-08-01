# TASSEL Test Suite — Issues & Remediation Plan

**Audience:** developers joining to fix the test suite.
**Status:** the build passes green, but that green is not trustworthy. This document
explains why, inventories the problems, and lays out a phased plan to restore
enforceable test coverage — with special attention to the **statistical-correctness
tests**, which are the highest priority to protect (TASSEL is the analysis core for
`rTASSEL` and is used to teach GWAS to ~25,000 students/year).

---

## TL;DR

- **148** JUnit test classes exist under `src/test/java`.
- **47** are hard-excluded from execution in `build.gradle.kts`.
- **89** of the 148 depend on external data that is **not in the repository**
  (`dataFiles/` is git-ignored) or on **hardcoded developer-machine paths**
  (`/Users/pbradbury/...`, `/Volumes/Macintosh HD 2/...` — 32 files).
- **`ignoreFailures = true`** in the test task means **any failure is swallowed**
  and the build still reports `BUILD SUCCESSFUL`.
- Net effect: the statistical checks (kinship, MLM/GLM, PCA, LD, distance) are
  present and mathematically meaningful, but on a clean checkout they **cannot find
  their data, fail, and are silently ignored**. They are effectively advisory today.

**Proof (reproducible now):**
```
./gradlew test --tests "net.maizegenetics.analysis.distance.KinshipTest"
# KinshipTest > testKinship FAILED  (dataFiles/Tutorial/mdp_genotype.hmp.txt doesn't exist)
# KinshipTest > testEndelman FAILED
# ...
# BUILD SUCCESSFUL   <-- failures ignored
```

---

## Root causes

### 1. `ignoreFailures = true`
[`build.gradle.kts`](build.gradle.kts) test block:
```kotlin
ignoreFailures = true // currently setting this to 'true' until we figure out failing tests
```
This is the single most damaging setting: it converts every red test into a
non-event. Nothing in CI or local builds can regress and be noticed.

### 2. Test data lives outside the repo
Test inputs and R-validated expected results are referenced through constants in
`src/test/java/net/maizegenetics/constants/`:
```java
// GeneralConstants.java
public static final String DATA_DIR = "dataFiles/";
public static final String EXPECTED_RESULTS_DIR = DATA_DIR + "ExpectedResults/";
// TutorialConstants.java
public static final String TUTORIAL_DIR = DATA_DIR + "Tutorial/";
public static final String HAPMAP_FILENAME = TUTORIAL_DIR + "mdp_genotype.hmp.txt";
```
…but `dataFiles/` (and `tempDir/`) are git-ignored ([`.gitignore`](.gitignore) lines 46–47).
A fresh clone therefore has **no** test data, so all 89 data-dependent tests fail
at the "input doesn't exist" check.

### 3. Hardcoded absolute paths (32 files)
Many tests point at a specific developer's disk, e.g.:
```
/Users/pbradbury/Documents/projects/tassel/validation/kinshipTestData.hmp.txt
/Volumes/Macintosh HD 2/temp/fast_assoc_test.txt
```
These can never pass on CI or another developer's machine and must be routed
through the `dataFiles/`-relative constants instead.

### 4. Native / external dependencies not provisioned
- **OpenBLAS** is required by some numeric tests (`build.gradle.kts` wires
  `-Djava.library.path` to `/opt/homebrew/opt/openblas/lib` etc.). Not installed →
  those tests fail.
- **HDF5 (jhdf5 native lib)**, **PostgreSQL/GOBII**, **MonetDB**, and GBS aligners
  are required by whole categories of excluded tests.

### 5. Tests are skipped entirely in the distribution build
`./gradlew installDist -x test` (and the packaging path) skip tests, so the app can
ship without any test ever executing.

---

## Inventory of excluded tests (47 classes)

Excluded via `exclude(...)` in the `test { }` block of `build.gradle.kts`.

| Category | Count | Reason | Statistical correctness? |
|---|---:|---|---|
| GBS pipeline — `analysis/gbs`, `gbs/v2`, `gbs/repgen` | 26 | Need FASTQ/BAM, SQLite/HDF5 DBs, external aligners, large fixtures | No — sequencing/bioinformatics infra |
| External DBs — `analysis/gobii`, `analysis/monetdb` | 8 | Need live Postgres / MonetDB | No |
| HDF5-native — `LowLevelCopyOfHDF5Test`, `TagsOnPhysMapHDF5Test`, `SplitHDF5ByChromosomePluginTest`, `BuildUnfinishedHDF5GenotypesPluginTest`, `DistanceMatrixHDF5Test` | 5 | Need jhdf5 native lib + big files | 1 is distance **I/O** (math covered elsewhere) |
| RNA — `analysis/rna` | 2 | Sequencing DB fixtures | No |
| Misc plugin/data — `GenomeFeatureBuilderTest`, `BasicGenotypeMergeRuleTest`, `ThinSitesByPositionPluginTest`, `LDKNNiImputationPluginTest`, `GenomeAnnosDBQueryToPositionListPluginTest` | 5 | Hardcoded paths / heavy fixtures | `LDKNNiImputation` = imputation |
| `FastMultithreadedAssociationPluginTest` | 1 | Performance variant, hardcoded `/Volumes/...` path | **Yes (association)** — but math overlaps `MLMTest`/`EqtlAssociationPluginTest` |

**Key point:** only ~2 of the 47 excluded classes are statistical-correctness tests,
and both are I/O- or performance-variants whose underlying math is exercised by tests
that are *not* excluded. **The exclusion list is not the main threat to statistical
coverage — the silent-failure problem (root causes 1–3) is.**

---

## Statistical-correctness tests that SHOULD be the protected core

These are present and **not** on the exclude list. They are the ones to get running
and enforced first. All currently fail on a clean checkout for lack of `dataFiles/`.

| Area | Test classes (package `net.maizegenetics.*`) |
|---|---|
| Association / GWAS | `analysis.association.MLMTest`, `ReferenceProbabilityFELMTest`, `PhenotypeLMTest`, `GenomicSelectionPluginTest`, `EqtlAssociationPluginTest`, `DiscreteSitesTest` |
| Kinship / distance | `analysis.distance.KinshipTest`, `CenteredIBSTest`, `NormalizedIBSTest`, `DominanceCenteredIBSTest`, `DominanceNormalizedIBSTest`, `IBSDistanceMatrixTest`, `AMatrixPluginTest` |
| Linear models | `stats.linearmodels.ModelEffectTest`, `SolveByOrtholgonalizingTest` |
| PCA | `stats.PCA.PrinCompTest` |
| Statistics | `stats.statistics.FisherExactTest` *(passes today — no external data)* |
| Linkage disequilibrium | `popgen.LinkageDisequilibriumTest` |
| Model fitting | `analysis.modelfitter.StepwiseAdditiveModelFitterTest`, `AdditiveSiteTest` |
| Matrix algebra | `matrixalgebra.Matrix.DoubleMatrixTest` *(passes today)* |
| Numeric transforms | `analysis.numericaltransform.*` (5) |

---

## Remediation plan (phased)

### Phase 0 — Make the current state visible (½ day)
- Run the full suite once **with `ignoreFailures = true` left as-is** and publish the
  HTML/XML report (`build/reports/tests/test/index.html`) so everyone sees the true
  pass/fail/error counts. This is the baseline.
- Record which failures are "missing data" vs "assertion/logic" vs "missing native lib".

### Phase 1 — Restore test data & remove hardcoded paths (owner: 1 dev, ~1 week)
1. Assemble the `dataFiles/` tree the tests expect (`Tutorial/`, `ExpectedResults/`,
   `CandidateTests/`, `GenomeFeatures/`, …). Source from the existing developer copies.
2. Decide storage: commit a **small** fixture subset into `src/test/resources`
   (preferred for the statistical core — a few hundred KB of genotype + expected
   matrices), and/or fetch larger data in CI via Git LFS or a download step.
3. Replace the **32** hardcoded `/Users/...` and `/Volumes/...` paths with the
   `GeneralConstants` / `TutorialConstants` relative paths.
4. Verify each statistical-core test finds its input and expected-results file.

### Phase 2 — Split a "must-pass" statistics group and stop ignoring failures (owner: 1 dev, ~3 days)
1. Tag the statistical-correctness classes (JUnit 4 `@Category`, or migrate those to
   JUnit 5 `@Tag("statistics")`).
2. Add a Gradle task, e.g. `statisticsTest`, that runs **only** that group with
   **`ignoreFailures = false`** — this becomes a required CI gate.
3. Keep the broad `test` task non-blocking for now (pipeline/IO tests stay
   `ignoreFailures = true`) so the team can burn those down incrementally without
   blocking releases.
4. Wire `statisticsTest` into `.github/workflows` as a required check.

### Phase 3 — Provision native/external deps in CI (owner: devops, ~2–3 days)
- Install **OpenBLAS** on CI runners (already parameterized via `BLAS_LIB_PATH` /
  the OS detection in `build.gradle.kts`).
- Provision **HDF5** so the HDF5 storage tests can be un-excluded where they cover
  real behavior (esp. `DistanceMatrixHDF5Test`).
- Stand up disposable **Postgres**/**MonetDB** (containers) for the GOBII/MonetDB
  tests, or formally mark those as integration-only and document that decision.

### Phase 4 — Burn down the broad suite & raise coverage (ongoing, multiple devs)
- Work category-by-category (GBS v2 first — largest block) to fix or explicitly
  quarantine each excluded test with a documented reason.
- Once the broad `test` task is green, flip its `ignoreFailures` to `false` too.
- Raise the Kover line-coverage minimum (currently 15% in `build.gradle.kts`) in
  steps as coverage improves.

---

## Definition of done
- [ ] `dataFiles/` fixtures available to every clone/CI run; **zero** hardcoded
      absolute paths remain in `src/test`.
- [ ] A `statisticsTest` CI gate runs the statistical-correctness classes with
      `ignoreFailures = false` and is **green because they pass**, not because
      failures are ignored.
- [ ] Every remaining exclusion in `build.gradle.kts` has a one-line comment stating
      why (native dep / external service / integration-only).
- [ ] Distribution/release builds run at least the `statisticsTest` gate (not
      `-x test` blanket skip).
- [ ] Kover minimum coverage raised above the current 15% floor.

---

## Reference — key files
- Test task config & exclusions: [`build.gradle.kts`](build.gradle.kts) (`test { }` block, `ignoreFailures`, `exclude(...)`)
- Path constants: `src/test/java/net/maizegenetics/constants/GeneralConstants.java`, `TutorialConstants.java`, `GenomeFeatureConstants.java`
- Git-ignored data dirs: [`.gitignore`](.gitignore) lines 46–47 (`dataFiles/`, `tempDir/`)
- Representative statistical test: `src/test/java/net/maizegenetics/analysis/distance/KinshipTest.java` (R-validated kinship, fails on missing data)

## Reference — useful commands
```bash
# Full suite (report shows true state even though build stays green)
./gradlew test ; open build/reports/tests/test/index.html

# A single class or package
./gradlew test --tests "net.maizegenetics.analysis.distance.*"

# Point BLAS explicitly if OpenBLAS is in a non-standard location
BLAS_LIB_PATH=/opt/homebrew/opt/openblas/lib ./gradlew test
```
