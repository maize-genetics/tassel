# GBS tests — handoff

Branch: `gbs-new-test` (off `main` @ `a5ea4cad`)
Goal: get the GBS tests into a sane, runnable state.

## Where we are

- Cloned `maize-genetics/tassel`, created branch `gbs-new-test`.
- Did a user-reachability call analysis of the GBS code — see
  [gbs-accessibility-analysis.md](gbs-accessibility-analysis.md).
- First code changes landed 2026-07-10 (see Log): `System.exit` removal + small/large
  dataset selection mechanism. Data-dependent validation still pending the data download.

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
- 2026-07-10 — **code change: removed `System.exit` from the test-hit legacy plugins**
  (root cause #1). Replaced `System.exit(1)` error paths with thrown `IllegalStateException`
  (carrying the original message/cause) in `FastqToTagCountPlugin` (2),
  `ProductionSNPCallerPlugin` (5), `ProductionPipelineMain` (2), `SimpleGenotypeSBit` (2).
  Several were the `try { throw ISE } catch { System.exit }` anti-pattern, simplified to a
  direct throw. Failures now stay inside the test JVM (no more "Gradle Test Executor
  finished with non-zero exit value 1"), so JUnit reports generate. `main()` methods and the
  untested legacy/`pana` classes were intentionally left alone. `compileJava` clean.
- 2026-07-10 — **code change: small/large dataset selection.** `GBSConstants.RAW_SEQ_CURRENT_TEST`
  now reads `-Dgbs.test.dataset` (default = 20 MB `Chr9_10-20000000/`, so existing behavior
  is unchanged). Added Gradle tasks `gbsTestSmall` (`Chr9_10-200000/`), `gbsTestLarge`
  (`Chr9_10-20000000/`), and kept `gbsTest` as an alias for the large set. `compileTestJava`
  clean; all three tasks register.
- 2026-07-10 — **data bundle arrived**; wired `dataFiles/` symlink → `gbs-test-data/dataFiles`.
  Ran the profiles and drove the real root causes to ground:
  - **TRUE root cause of the v2 NPE cascade = a build regression, not shared-DB corruption.**
    The modernized Gradle build only copies `src/main/resources`, but ~93 runtime resources
    (notably `dna/tag/TagSchema.sql`, plus GUI icons/xml/html) live co-located under
    `src/main/java` and were **absent from the classpath at runtime — in tests *and*
    production**. So `new TagDataSQLite()` built a *schemaless* DB → "no such table" →
    `myTaxaList` NPE. Fixed in `build.gradle.kts` by adding `src/main/java` as a resource root
    (excluding `*.java/*.kt/*.c/*.h`). This is the fix the run-#2 notes were missing — their
    "shared DB with no isolation" theory was wrong.
  - **System.exit, the deeper source:** the killer isn't only the GBS plugins — `AbstractPlugin`
    (base class of every plugin) does `printUsage();System.exit(1)` on any non-interactive
    parameter/`processData` error. Per Ed's decision we did **not** change `AbstractPlugin`
    (broad blast radius; `TasselPipeline` already exits non-zero via its own top-level catch).
    Instead we removed the stray `System.exit`s from the GBS plugins actually hit
    (added `ProductionSNPCallerPluginV2` ×1 and `DiscoverySNPCallerPlugin` ×6 to the earlier
    four) and made the tests never feed a plugin missing inputs (below).
  - **Redundant DB rebuilds (Ed's efficiency point):** ~10 raw-fastq→tag-DB builds across the
    suite, zero `@BeforeClass` sharing. With the resource bug fixed, sharing is now *safe*
    (the thing that looked like cross-test corruption was the resource bug). Added a per-class
    `@BeforeClass` shared build to `GBSSeqToTagDBPluginTest` (6 builds → 1); read-only tests
    reuse it, the mutating append/keepOld tests still rebuild (that's what they verify).
  - **Small vs large split, expressed in the tests:** the two Gradle profiles now differ *only*
    by `-Dgbs.test.dataset`. 20 MB-only tests self-skip on the small dataset via
    `Assume.assumeTrue(RAW_SEQ_CURRENT_TEST == Chr9_10-20000000)` —
    `EvaluateSNPCallQualityOfPipelineTest` (needs the bowtie SAM), `DiscoverySNPCallerPluginV2Test
    .testNumTagsWithReference`, `ProductionSNPCallerPluginTest`, `ProductionPipelineMainTest`.
    Added dir-creation `@Before`s so plugins never see a missing `-db` parent dir.
  - **`@Ignore`d two genuinely unrunnable tests:** `GBSSeqToTagDBPluginTest.testSAMImportPlugin`
    (needs an external bowtie run) and all of `RNADeMultiPlexSeqToDBPluginTest` (hardcoded
    `/Users/lcj34/...` paths that exist on no machine).
  - **Result — `./gradlew gbsTestSmall` is GREEN: 38 passed, 13 skipped, 0 failed, ~1m37s,
    no JVM kill.** `gbsTestLarge` (full 20 MB) running to confirm; expected to clear the run-#2
    v2 NPEs thanks to the resource fix.

- 2026-07-10 (night) — **`gbsTestLarge` (20 MB) finished: 54m17s, 44 pass / 2 fail / 5 skip.**
  Both failures are legacy GBSv1 Production assertion tests (not chased per Ed's GBSv2 focus);
  **all six run-#2 v2 NPEs are gone** and 22 JUnit XMLs generated (both confirm the resource +
  System.exit fixes on the large dataset). Then produced the overnight roadmap + a JFR profile —
  see **[roadmap-and-profiling.md](roadmap-and-profiling.md)**. Headlines:
  - **Timing:** legacy GBSv1 = 88% of the run (SeqToTBTHDF5 28min + ModifyTBTHDF5 14min +
    FastqToTagCount 4min); GBSv2 = 12%. **Deleting GBSv1 alone takes the suite 54min → ~6.5min.**
  - **Profile (DB-build):** #1 hot spot (~40%) is `BarcodeTrie.TrieNode.containsKey`, which
    allocates a fresh `ArrayList` + linear-scans per read char — replaceable with
    `children[c-'A'] != null`. `seqDifferences` is NOT hot here (it's an alignment-phase method).
  - Roadmap also contains the full **legacy-deletion plan** (only `Barcode` must be preserved)
    and the **GBSv2 test-speedup plan**.

  **Still open (for review):** commit the test-remediation changes; execute the GBSv1 deletion as
  its own PR; try the BarcodeTrie fix (4a) + re-profile; `enzymes.ini` warning still cosmetic.
