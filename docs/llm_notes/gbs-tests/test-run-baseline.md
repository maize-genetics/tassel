# GBS test run — baseline (2026-07-10)

Branch: `gbs-new-test`. Task: `./gradlew gbsTest` (new task added to `build.gradle.kts`).
Dataset: 20 MB (`Chr9_10-20000000`), as pinned by `GBSConstants.RAW_SEQ_CURRENT_TEST`.
Input staged: valid `C05F2ACXX_5_fastq.gz` only — the 20 MB `C08L7ACXX_6_fastq.gz` from the
recovered tarball is corrupt (truncated), so it was left out.

## Build fixes needed to run at all
- `openchart:openchart:1.4.2` (transitive of biojava `forester`) is unresolvable in every
  repo → added a global `configurations.all { exclude(... openchart ...) }`.
- GBS tests are excluded from the main `test` task → added a dedicated `gbsTest` task
  including `**/analysis/gbs/*Test.class` and `**/analysis/gbs/v2/*Test.class`.
- Symlinked `dataFiles/` (bitbucket `tassel-5-test` + staged fastqs) into the project root,
  since tests use relative `dataFiles/` / `tempDir/` paths.

## Result: 26 tests, 9 failed, 1 skipped (run time ~19 min)

### Legacy GBS — passing where data is intact
Verified from the run log (exact MD5 matches against `ExpectedResults/`):
- `FastqToTagCountPlugin` → `C05F2ACXX_5.cnt` MD5 `e31d4c16…` matches. 188,975 tags.
- `MergeMultipleTagCountPlugin` → `Merged_Tag_Counts.cnt` MD5 `3245c24d…` matches.
- `SAMConverterPlugin`, `DiscoverySNPCallerPlugin`, `ModifyTBTHDF5Plugin`, `SeqToTBTHDF5Plugin`
  all executed against 20 MB expected results.
- One legacy path hit a null `InputStream` / HDF5 read error on
  `tempDir/GBS/Pipeline_Testing.hmp.txt.gz` (ProductionSNPCaller-style flow needing a file
  not produced) — needs follow-up.

### GBSv2 — 9 failures, all the same shape
All failures are `NullPointerException` at `new TagDataSQLite(GBS_GBS2DB_FILE)` /
`tagData.getTags()` / `getSNPPositions()`:
- `DiscoverySNPCallerPluginV2Test.testNumTagsWithReference` (line 377)
- `EvaluateSNPCallQualityOfPipelineTest` ×4 (`pipelineIncludingInvariantSites`,
  `pipeLineWithVariableSitesTest`, `agpv2WithVariableSitesTest`, `agpv2IncludingInvariantSites`) (line 340)
- `GBSSeqToTagDBPluginTest.testKeepOldData` (line 114)
- (plus a few more not captured in the truncated console; total 9)

Passing v2 tests: `testChromPosition`, `testFullSNPCaller`, `testAlignTags`,
`testSinglePositionFiltering`, `GBSSeqToTagDBPluginTest.testTagExportPlugin`,
`EvaluateSNPCallQualityOfPipelineTest.testBiologyOfDiscoveryTBT`.
Skipped: `GBSSeqToTagDBPluginTest.testSAMImportPlugin` (`@Ignore`).

### Root cause of the v2 NPEs
Run log shows the v2 `GBSSeqToTagDBPlugin` produced a near-empty DB: `t=0 tagsRead=2 outCnt=0`.
An empty tag DB → `getTags()` returns empty/null → NPE in the assertions that read it.
The same `C05` fastq yields 188,975 tags via the legacy plugin, so the fastq is fine — the
empty v2 DB is consistent with the 20 MB dataset missing its valid second fastq (`C08`).

### Red herring
`EnzymeList - ERROR! Cannot load Enzyme List -- …/ini4j/…/enzymes.ini` is only a warning:
`EnzymeList` (v2) looks for `enzymes.ini` next to the ini4j jar, fails, and falls back to
`loadDefaults()`. Enzymes still load. Not a cause of the failures (but the misleading path
lookup is worth cleaning up).

## Next steps
1. Drop a **valid 20 MB `C08L7ACXX_6_fastq.gz`** into `dataFiles/GBS/Chr9_10-20000000/`
   (from a better tarball / lab backup) and re-run — expected to clear most v2 NPEs.
2. If a valid 20 MB `C08` can't be found, re-point `RAW_SEQ_CURRENT_TEST` to the fully-intact
   200 KB dataset and regenerate v2 expected tag-sequence values.
3. Investigate the legacy `Pipeline_Testing.hmp.txt.gz` null-stream path separately.
4. Optional cleanup: make `EnzymeList` load `enzymes.ini` from classpath resources.
