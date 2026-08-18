# GBS test run #2 (complete 20 MB dataset) — results, root causes, timing

Branch: `gbs-new-test`. Task: `./gradlew gbsTest --rerun-tasks`.
Dataset: 20 MB (`Chr9_10-20000000`) with **both** raw fastqs now valid
(`C05` 166 MB, `C08` 156 MB). Raw log tail: `notes/gbs-tests/logs/gbsTest-run2-complete-dataset.tail.txt`.

## Result: 26 tests, 8 failed, 1 skipped — run time 1h 21m

Adding the valid `C08` fixed only **one** test vs run #1 (9 → 8 failures). The GBSv2
NPEs did **not** clear, which disproves the "empty DB from missing fastq" theory — the
v2 failures are a real code/test issue, not a data gap.

### Passing (17)
Legacy: `DiscoverySNPCallerPluginTest`, `FastqToTagCountPluginTest`,
`MergeMultipleTagCountPluginTest`, `ModifyTBTHDF5PluginTest`, `ParseBarcodeReadTest`,
`SAMConverterPluginTest`, `SeqToTBTHDF5PluginTest`, `TagCountToFastqPluginTest`.
GBSv2: `BarcodeTrieTest`, `DiscoverySNPCallerPluginV2Test` (`testFullSNPCaller`,
`testFilterAlignedTags`, `testCreateReferenceTag`, `testChromPosition`,
`testSinglePositionFiltering`, `testAlignTags`), `GBSSeqToTagDBPluginTest.testTagExportPlugin`,
`EvaluateSNPCallQualityOfPipelineTest.testBiologyOfDiscoveryTBT`.

### Skipped (1)
`GBSSeqToTagDBPluginTest.testSAMImportPlugin` (`@Ignore`).

### Failed (8)
| test | error | line |
|---|---|---|
| `ProductionPipelineMainTest.testProductionPipelineMain` | AssertionError (checksum) | 76 |
| `ProductionSNPCallerPluginTest.testPerformFunction` | AssertionError (genotypes) | 55 |
| `DiscoverySNPCallerPluginV2Test.testNumTagsWithReference` | NullPointerException | 377 |
| `EvaluateSNPCallQualityOfPipelineTest.pipeLineWithVariableSitesTest` | NullPointerException | 340 |
| `EvaluateSNPCallQualityOfPipelineTest.agpv2WithVariableSitesTest` | NullPointerException | 340 |
| `EvaluateSNPCallQualityOfPipelineTest.agpv2IncludingInvariantSites` | NullPointerException | 340 |
| `EvaluateSNPCallQualityOfPipelineTest.pipelineIncludingInvariantSites` | NullPointerException | 340 |
| `GBSSeqToTagDBPluginTest.testKeepOldData` | NullPointerException | 114 |

## Root causes

### 1. `System.exit(1)` in legacy plugins crashes the test executor
`FastqToTagCountPlugin`, `ProductionSNPCallerPlugin` (×5), `ProductionPipelineMain`, and
`SimpleGenotypeSBit` call `System.exit(1)` on error paths. When one fires inside a test,
it kills the worker JVM → `Process 'Gradle Test Executor 4' finished with non-zero exit
value 1`. This:
- fails the whole `gbsTest` task even though `ignoreFailures = true`, and
- **prevents JUnit XML/HTML report generation** (only `build/test-results/gbsTest/binary/`
  is written) — so we have no per-test timings and must scrape the console log.

Fix: replace `System.exit` in these code paths with thrown exceptions so failures stay
inside the test JVM.

### 2. GBSv2 NPEs at `TagDataSQLite` reads — not data-related
All six v2 failures NPE at `new TagDataSQLite(GBS_GBS2DB_FILE)` / `getTags()` /
`getSNPPositions()` (lines 114, 340, 377). Valid `C08` did not change this. Prime
suspects, in order:
- The suite shares **one** DB path (`GBS_GBS2DB_FILE`) across all v2 tests with no
  isolation; a crashed/half-written DB from an earlier test (see #1) leaves later reads
  returning null. `testKeepOldData` deletes then rebuilds the shared DB.
- `EvaluateSNPCallQualityOfPipelineTest` failures are all at the same line (340) reached
  via helper `evaluateConservedSites()` — likely the SNP-position query returns null when
  the DB has no SNPs written.
Needs a focused single-test debug run (isolate one v2 test, inspect the DB it builds).

### 3. Legacy assertion failures (2)
`ProductionPipelineMain` (checksum of `Pipeline_Testing.hmp.txt.gz`) and
`ProductionSNPCaller` (genotype comparison) mismatch the 20 MB expected results. Could be
genuine output drift or fallout from #1 (partial outputs). Re-check in isolation.

## Why it is so slow, and how to fix

### Where the time goes
- **Data volume:** 20 MB dataset = ~306 MB of fastq (158M + 148M). The 200 KB dataset is
  ~308 KB — about **1000× smaller**. Run #1 (only `C05`, 158 MB) took 19 min; run #2 (both,
  306 MB) took 1h 21m. Non-linear jump also reflects machine contention during the run.
- **Redundant rebuilds:** `GBSSeqToTagDBPlugin`/`FastqToTagCountPlugin` is invoked **~40+
  times** across the suite (14× `GBSSeqToTagDBPluginTest`, 10× `EvaluateSNPCallQuality`,
  6× `DiscoverySNPCallerPluginV2`, …). There is **no shared setup** — every test re-parses
  the full fastqs and rebuilds the tag DB from scratch, then writes/compares large
  HDF5/DB files.
- **No parallelism:** all tests write the same `GBS_GBS2DB_FILE` and `tempDir` paths, so
  they cannot run in parallel without colliding.

### Recommendations (highest impact first)
1. **Run the suite on the 200 KB dataset.** Point `GBSConstants.RAW_SEQ_CURRENT_TEST` at
   `Chr9_10-200000/` (fastqs already staged, `ExpectedResults/GBS/Chr9_10-200000/` present).
   Regenerate the expected values that are dataset-specific (legacy MD5s; the few hardcoded
   v2 tag counts/sequences). Expected effect: minutes → **seconds**. Keep the 20 MB dataset
   as an opt-in "slow"/nightly profile.
2. **Isolate per-test state.** Give each test its own temp DB (e.g. a `@Rule TemporaryFolder`
   or unique db filename) instead of the shared `GBS_GBS2DB_FILE`. Removes the cross-test
   corruption cascade (#2) *and* unlocks parallelism.
3. **Enable parallel forks** once state is isolated: `maxParallelForks = N` on the test
   task. With `-Xmx10g` per fork, size N to available RAM.
4. **Share expensive fixtures** where a test only *reads* the DB: build it once in
   `@BeforeClass` rather than per method.
5. **Replace `System.exit` with exceptions** (root cause #1) so a single failure doesn't
   abort the run and reports are always produced.
6. **Don't use `--rerun-tasks` for iterative work** — rely on Gradle up-to-date checks;
   only the test logic needs to re-run, not a full recompile.

## Next actions
- Prototype the 200 KB switch on a couple of representative tests (one legacy MD5, one v2)
  to confirm the speedup and quantify expected-value regeneration effort.
- Debug the shared-DB NPE cascade by running one `EvaluateSNPCallQuality` test in isolation.
- File the `System.exit` removal as a small, self-contained fix.
