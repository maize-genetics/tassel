# GBS user-accessibility & test analysis

Branch: `gbs-new-test`
Date: 2026-07-10
Scope: understand whether the legacy GBS pipeline is reachable by users before
deciding how much to invest in fixing the GBS tests.

## TL;DR

- There are exactly **two user entry points**: the Swing GUI (`TASSELMainFrame`)
  and the CLI (`TasselPipeline`).
- The **legacy GBS pipeline is deliberately being retired**: its entire GUI menu
  is rendered struck-through (deprecated), and most legacy classes are not in any
  menu at all.
- **GBSv2 is the active, shipped pipeline** and is where test effort is worthwhile.
- The real blocker on the tests is **missing large data files**, not `@Ignore`:
  every GBS test (legacy *and* v2) resolves inputs through `GBSConstants` to
  `dataFiles/GBS/...`, and **none of `dataFiles/` is in the repo**.

## Entry points

| Entry point | File | How a class becomes reachable |
|---|---|---|
| GUI | `src/main/java/net/maizegenetics/tassel/TASSELMainFrame.java` | Must be wired into a `JMenu`. Two GBS menus: `getGBSv2Menu()` (active) and `getGBSMenu()` (whole menu struck through via `addMenuItemDeprecated`). |
| CLI | `src/main/java/net/maizegenetics/pipeline/TasselPipeline.java` (~line 1594) | Resolves `-AnyPluginClassName` by reflection (`Utils.getFullyQualifiedClassNames` + `Plugin.getPluginInstance`). Any `Plugin` subclass on the classpath is callable by name, menu or not. |

## Reachability tiers (legacy `analysis/gbs/`)

1. **GBSv2 — active** (GUI menu + CLI), 9 plugins:
   `GBSSeqToTagDBPlugin`, `DiscoverySNPCallerPluginV2`, `ProductionSNPCallerPluginV2`,
   `SAMToGBSdbPlugin`, `TagExportToFastqPlugin`, `SNPQualityProfilerPlugin`,
   `UpdateSNPPositionQualityPlugin`, `GetTagSequenceFromDBPlugin`,
   `SNPCutPosTagVerificationPlugin`.

2. **Legacy GBS — in GUI but deprecated** (strikethrough menu + CLI), 11 plugins:
   `FastqToTagCountPlugin`, `MergeMultipleTagCountPlugin`, `TagCountToFastqPlugin`,
   `SAMConverterPlugin`, `SeqToTBTHDF5Plugin`, `ModifyTBTHDF5Plugin`,
   `DiscoverySNPCallerPlugin`, `ProductionSNPCallerPlugin`, `BinaryToTextPlugin`,
   `UTagCountToTagPairPlugin`, `UTagPairToTOPMPlugin`.

3. **Legacy GBS — CLI-only** (no menu; reachable only if the user knows the exact
   class name), ~22 plugins:
   `FastqToTBTPlugin`, `QseqToTagCountPlugin`, `QseqToTBTPlugin`, `KmerToTagCountPlugin`,
   `KmerToTBTPlugin`, `MergeDuplicateSNPsPlugin`, `MergeMultipleTOPMPlugin`,
   `MergeTagsByTaxaFilesPlugin`, `MergeTagsByTaxaFilesByRowPlugin`, `TOPMSummaryPlugin`,
   `TagAgainstAnchorPlugin`, `TextToBinaryPlugin`, `AnnotateTOPMwSAMPlugin`,
   `SAMWGMapConverterPlugin`, `KeepSpecifiedSitesInTOPMPlugin`,
   `KeepSpecifiedReadsinFastqPlugin`, `CompareGenosBetweenHapMapFilesPlugin`,
   `ContigPETagCountPlugin`, `FastqToPETagCountPlugin`, `MergePETagCountPlugin`,
   `QseqToPETagCountPlugin`, `ProductionPipeline`.

4. **Not user-reachable** — non-plugins (internal helpers, `main()`-only scripts,
   or dead code), ~23 classes: `Barcode *`, `ParseBarcodeRead *`, `ReadBarcodeResult *`,
   `ProductionPipelineMain (main)`, `TerryPipelines (main)`, `TagMatchFinder (main)`,
   `TagAgainstAnchor`, `TagAgainstAnchorHypothesis`, `PolymorphismFinder`,
   `SmithWaterman`, `Clusters`, `UNetworkFilter`, `UTagPairFinder`, `SimpleGenotypeSBit`,
   `AnnotateTOPM`, `SNPLogging`, `TagBlockPosition`, `PEParseBarcodeRead`,
   `PEReadBarcodeResult`, `ShortReadBarcodeResult`.

   `*` = **live dependency of GBSv2/RNA.** `Barcode` and the barcode-parsing helpers
   are imported by `RNADeMultiPlexSeqToDBPlugin` and
   `ConvertOldFastqToModernFormatPlugin`, so the legacy package cannot be deleted
   wholesale even though its own plugins are retired.

## Test situation

### `@Ignore` is small
Only 3 ignored tests total across all GBS tests:
- `analysis/gbs/v2/SAMToGBSdbPluginTest.java` — 2
- `analysis/gbs/repgen/RepGenLoadSeqToDBPluginTest.java` — 1

### Missing large files is the real blocker
- Inputs resolve through `src/test/java/net/maizegenetics/constants/GBSConstants.java`:
  - `GBS_DATA_DIR = "dataFiles/GBS/"`
  - `GBS_INPUT_DIR = dataFiles/GBS/Chr9_10-20000000/`
  - fastq input `C05F2ACXX_5_fastq.gz`, reference `ZmB73_RefGen_v2_chr9_10_1st20MB.fasta`,
    expected TOPM/HDF5 results under `dataFiles/ExpectedResults/GBS/...`.
- `GeneralConstants.DATA_DIR = "dataFiles/"` — a **relative** path resolved against the
  test working directory.
- **Nothing under `dataFiles/` is committed.** No `.gitattributes`, no git-LFS, no
  download step in `build.gradle.kts`.
- Both legacy and GBSv2 tests use the same `GBSConstants.GBS_INPUT_DIR`, so the whole
  GBS test corpus is blocked on the same missing inputs.

### Legacy test shape
The legacy GBS tests (e.g. `FastqToTagCountPluginTest`) are not `@Ignore`d but hard-fail
on missing inputs — they `mkdirs()` an output dir then read `GBS_INPUT_DIR`, and compare
MD5s against `dataFiles/ExpectedResults/...` that don't exist.

## Implication for effort

- Fixing **legacy** GBS tests means sourcing large old B73 GBS data for a pipeline the
  UI already marks deprecated and that is mostly CLI-invisible — low value.
  Defensible options: quarantine (`@Ignore("<reason>")`) or delete the tests whose
  target class is in the CLI-only / not-reachable tiers.
- Fixing **GBSv2** tests is worthwhile (active pipeline), but step one is locating the
  `dataFiles/GBS/` corpus — decide between *restore the data* vs *shrink tests to small
  bundled fixtures*.

## Test-data source (resolved 2026-07-10)

Data repo: `https://bitbucket.org/tasseladmin/tassel-5-test` (branch `master`).
It is laid out as a working directory — `dataFiles/`, `tempDir/`, plus an old
snapshot of `src/` and `run_gbs_suite.pl`. Tests were historically run with this
repo as the CWD, which is why `GBSConstants` paths are relative to `dataFiles/`.

What it provides (`dataFiles/` ≈ 484 MB total):
- `dataFiles/GBS/` (41 MB): `ZmB73_RefGen_v2_chr9_10_1st20MB.fasta` (42 MB, real file,
  no LFS), `Pipeline_Testing_key.txt`, `Pipeline_Testing_Mirror_key.txt`, `README.txt`.
- `dataFiles/ExpectedResults/GBS/` (99 MB): **all expected outputs** for both the
  200 KB and 20 MB datasets — `.cnt`, `.topm`, `.h5`, `.sam`, hmp, donor haplos, etc.

What it does **not** provide — the raw input fastqs:
- `C05F2ACXX_5_fastq.gz` / `C08L7ACXX_6_fastq.gz`, 20 MB dataset (166 MB + 156 MB) and
  200 KB dataset (180 KB + 135 KB). Never committed to git; `README.txt` points to the
  now-dead `maizegenetics.net/tassel/GBSTestData.tar`.

## Raw fastqs recovered from the Wayback Machine (2026-07-10)

Retrieved `GBSTestData.tar` (209 MB) — saved at `/Users/esb33/Developer/GBSTestData.tar`.
The archive is **truncated** (Wayback capture cut off inside the last file). Integrity
after extraction (`gzip -t`):

| Dataset | C05F2ACXX_5 | C08L7ACXX_6 |
|---|---|---|
| `Chr9_10-200000` (200 KB) | OK (180 KB) | OK (135 KB) |
| `Chr9_10-20000000` (20 MB) | OK (166 MB) | **CORRUPT** (only ~43 of 156 MB) |

Net: the **200 KB dataset is fully usable**; the 20 MB dataset is missing a valid
second fastq. Valid fastqs have been staged into the local
`tassel-5-test/dataFiles/GBS/<dataset>/` layout the tests expect, giving a complete
runnable `dataFiles/` locally (the corrupt 20 MB `C08` was intentionally not copied).

### 200 KB vs 20 MB feasibility

`ExpectedResults/GBS/Chr9_10-200000/` covers 8 plugins: `FastqToTagCountPlugin`,
`MergeMultipleTagCountPlugin`, `TagCountToFastqPlugin`, `SAMConverterPlugin`,
`SeqToTBTHDF5Plugin`, `ModifyTBTHDF5Plugin`, `DiscoverySNPCallerPlugin`, `bowtie2.1`.
It is **missing** `ProductionSNPCallerPlugin` and `maizeTestDataDonorHaplos` (imputation)
— those exist only under the 20 MB tree, whose second fastq is corrupt.

Legacy tests compare output-file MD5s against `ExpectedResults/`, so they can be
re-pointed to the 200 KB dataset (except the Production/imputation ones). GBSv2 tests
(`GBSSeqToTagDBPluginTest`) assert relational invariants plus a few **dataset-specific
expected tag sequences** — switching them to 200 KB would require regenerating those
expected values. `GBSConstants.RAW_SEQ_CURRENT_TEST` currently = `Chr9_10-20000000/`.

### Still to source
- A valid 20 MB `C08L7ACXX_6_fastq.gz` (a Buckler-lab backup, or an earlier/complete
  Wayback capture), if the 20 MB pipeline is to be kept.
- Fix the dead download URL in `dataFiles/GBS/README.txt` regardless.

Consequence: `GBSConstants.RAW_SEQ_CURRENT_TEST = Chr9_10-20000000/`, so the tests
currently expect the 20 MB raw fastqs. Every test that starts from raw fastq — legacy
(`FastqToTagCountPlugin`, `SeqToTBTHDF5Plugin`, …) **and** GBSv2
(`GBSSeqToTagDBPluginTest` → `GBS_INPUT_DIR`) — is blocked on those missing fastqs.
Tests that start from an intermediate artifact present in `ExpectedResults/` could run
without the raw fastqs.

## Open questions

1. **Raw fastqs** — where can we get `C05F2ACXX_5_fastq.gz` / `C08L7ACXX_6_fastq.gz`?
   Options: (a) a Buckler-lab server/backup copy of `GBSTestData.tar`; (b) commit the
   tiny 200 KB-dataset fastqs and re-point `RAW_SEQ_CURRENT_TEST` to `Chr9_10-200000/`;
   (c) build small synthetic fixtures and regenerate expected results. The dead download
   URL should also be fixed in `README.txt` regardless.
2. **Where the data lives at build time** — copy `dataFiles/` into the tassel repo,
   reference it via env var/CWD, or add a Gradle fetch step. Decide before wiring CI.
3. **Legacy tier** — fix, quarantine, or delete the deprecated-pipeline tests?
