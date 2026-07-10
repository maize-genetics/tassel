# GBS tests — handoff

Branch: `gbs-new-test` (off `main` @ `a5ea4cad`)
Goal: get the GBS tests into a sane, runnable state.

## Where we are

- Cloned `maize-genetics/tassel`, created branch `gbs-new-test`.
- Did a user-reachability call analysis of the GBS code — see
  [gbs-accessibility-analysis.md](gbs-accessibility-analysis.md).
- No code changes yet. This is planning/notes only.

## Key facts (see analysis doc for detail)

- Legacy GBS is deprecated in the GUI (struck-through menu) and mostly CLI-invisible.
  GBSv2 is the active pipeline.
- All GBS tests are blocked on missing `dataFiles/GBS/...` inputs — nothing under
  `dataFiles/` is in the repo, and there's no LFS/download mechanism.
- Only 3 tests are `@Ignore`d (2 `SAMToGBSdbPluginTest`, 1 `RepGenLoadSeqToDBPluginTest`).

## Test data source (resolved)

Most data is in `https://bitbucket.org/tasseladmin/tassel-5-test` (branch `master`,
`dataFiles/` ≈ 484 MB): reference genome, key files, and **all expected results**.
The **one gap** is the raw input fastqs (`C05F2ACXX_5_fastq.gz`, `C08L7ACXX_6_fastq.gz`)
— never in git; the README's download URL (`maizegenetics.net/tassel/GBSTestData.tar`)
is **dead** (now returns a Wix HTML page). See the analysis doc for full detail.

## Decisions needed before coding

1. **Raw fastqs** — locate `GBSTestData.tar` (lab backup?), OR commit the tiny
   200 KB-dataset fastqs and re-point `RAW_SEQ_CURRENT_TEST`, OR build synthetic
   fixtures + regenerate expected results.
2. **Build-time data wiring** — vendor `dataFiles/` vs env-var/CWD vs Gradle fetch.
3. **Legacy tier disposition** — fix / quarantine / delete tests targeting
   deprecated or non-reachable classes.

## Suggested plan (once decisions are made)

1. Establish how test data is provided (fetch script or bundled fixtures).
2. Prioritize GBSv2 tests (active pipeline): make them green with available data.
3. Legacy tests: quarantine with an explicit `@Ignore("deprecated GBS pipeline — see notes/gbs-tests")`
   or delete, tier-by-tier per the analysis doc.
4. Verify `./gradlew test` selection runs the intended GBS suite; capture what still
   needs external data.

## Useful pointers

- Test data constants: `src/test/java/net/maizegenetics/constants/GBSConstants.java`
- CLI plugin resolution: `src/main/java/net/maizegenetics/pipeline/TasselPipeline.java` (~1594)
- GUI menus: `src/main/java/net/maizegenetics/tassel/TASSELMainFrame.java`
  (`getGBSv2Menu()`, `getGBSMenu()`)
- GBS tests: `src/test/java/net/maizegenetics/analysis/gbs/{,v2,repgen}/`

## Non-git data bundle (for machine handoff)

All large non-git files live in `/Users/esb33/Developer/gbs-test-data/` (≈ 1.0 GB):
`dataFiles/` (complete corpus + all 4 raw fastqs, gzip-verified) and
`GBSTestData_full.tar` (archival). Copy that whole folder to the laptop and symlink its
`dataFiles/` into the repo root — see `gbs-test-data/README.md`. The repo's
`tassel/dataFiles` currently symlinks there.

## Log

- 2026-07-10 — branch created, accessibility analysis + notes written.
- 2026-07-10 — recovered raw fastqs; ran `gbsTest` baseline (26 tests, 9 v2 NPEs from an
  empty tag DB caused by the then-missing 20 MB `C08`). See `test-run-baseline.md`.
- 2026-07-10 — obtained complete `C08` (+ full tar) from lab file server; assembled the
  `gbs-test-data/` handoff bundle with all 4 fastqs valid.
- 2026-07-10 — ran `gbsTest` on the complete 20 MB dataset (1h 21m): 26 tests, 8 failed,
  1 skipped. Valid `C08` fixed only 1 test — the v2 NPEs are a real code/test issue, not a
  data gap. Full results + root causes + timing/perf recommendations in
  `test-run-2-results-and-timing.md`. Headline: switch the suite to the 200 KB dataset
  (minutes→seconds), isolate per-test DBs, and remove `System.exit` from legacy plugins.
