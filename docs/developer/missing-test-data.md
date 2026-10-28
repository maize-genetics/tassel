# Tests waiting on missing test data

Some tests stay in `src/test` but are marked `@Ignore` because they need files
that aren't in the [`tassel_test_data`](https://github.com/maize-genetics/tassel_test_data)
archive. They used to read those files from developers' own machines (see
[#24](https://github.com/maize-genetics/tassel/issues/24)). Their paths now point
to where each file would go under `dataFiles/`, so a test can be switched back on
once its files are in the archive.

Each `@Ignore` names what the test is waiting on and links
[#36](https://github.com/maize-genetics/tassel/issues/36).

## Which tests need which files

All of these classes are in `net.maizegenetics.dna.snp.io`.

| Test | Missing files |
| ---- | ------------- |
| `BuilderFromGVCFTest.pullRealV4Sequence` | W22 GVCF, full AGPv4 genome |
| `BuilderFromGVCFTest.pullW22FASTASequence` | W22 GVCF, full AGPv4 genome, gene-anchor CSV |
| `BuilderFromGVCFTest.pullB73FASTASequence` | Full AGPv4 genome, gene-anchor CSV |
| `BuilderFromGVCFTest.pullB73FASTASequenceFullGVCF` | Full AGPv4 genome, merged-anchor CSV |
| `BuilderFromGVCFTest.getChr8FromB73V4` | Full AGPv4 genome |
| `BuilderFromGVCFTest.getChr8FromB73V3` | Full AGPv3 genome |
| `BuilderFromGVCFTest.pullAnchor2kFromW22` | W22 GVCF, full AGPv4 genome |
| `BuilderFromGenotypeHDF5Test.testDeflationBySite` | HapMap v2 chromosome 8 HDF5 |
| `BuilderFromGenotypeHDF5Test.testDeflation` | HapMap v2 chromosome 8 HDF5 |

## The missing files

Paths are relative to `dataFiles/`.

### W22 GVCF

- **Path:** `CandidateTests/VCFFiles/W22_haplotype_caller_output.g.vcf.gz`
- **What it is:** a genomic VCF for the W22 inbred, aligned to B73 AGPv4. The
  original folder name (`W22SentieonTest`) suggests it was called with Sentieon's
  haplotype caller.
- **Where to get it:** unknown. It's the W22 counterpart of
  `B73_haplotype_caller_output.g.vcf.gz`, which is already in the archive (105 MB
  compressed). Both came from the pan-genome work of the tests' author (`zrm22`).
  If no copy survives, it has to be regenerated from W22 reads.

### Full B73 AGPv4 genome

- **Path:** `CandidateTests/Zea_mays.AGPv4.dna.toplevel.fa.gz`
- **What it is:** the whole B73 RefGen_v4 reference, all chromosomes.
- **Where to get it:** most likely Ensembl Plants, whose naming scheme this file
  name follows.
- The archive already has chromosome 8 on its own, as
  `CandidateTests/CHR8_Zea_mays.AGPv4.dna.toplevel.fa.gz`. `getChr8FromB73V4` is
  the code that made it.

### Full B73 AGPv3 genome

- **Path:** `CandidateTests/Zea_mays.AGPv3.29.dna.genome.fa.gz`
- **What it is:** the whole B73 RefGen_v3 reference. The `.29.` in the name is
  most likely the Ensembl Plants release number.
- Only `getChr8FromB73V3` needs it. That method made the chromosome 8 extract
  `CHR8_Zea_mays.AGPv3.dna.toplevel.fa.gz` that's already in the archive.

### Chromosome 8 anchor CSVs

- **Paths:** `CandidateTests/anchorsFile_chr8_genesExactCoordinates.csv` and
  `CandidateTests/anchorsFile_chr8_MergedPlus1000orGapDiffcoordinates.csv`
- **What they are:** lists of chromosome 8 regions. The tests skip the first line
  as a header, then read `chromosome,start,end` from the first three columns of
  each row. Going by the names, the first holds exact gene coordinates and the
  second holds merged gene regions padded by 1,000 bp.
- **Where to get them:** unknown; same origin as the W22 GVCF.

### Maize HapMap v2, chromosome 8, as HDF5

- **Path:** `CandidateTests/maizeHapMapV2_B73RefGenV2_201203028_chr8.hmp.h5`
- **What it is:** chromosome 8 of the public maize HapMap v2 release, converted
  from HapMap text format to TASSEL 5 HDF5.
- **Where to get it:** download `maizeHapMapV2_B73RefGenV2_201203028_chr8.hmp.txt.gz`
  from the Panzea HapMap v2 release and convert it with TASSEL. HDF5 needs unique
  taxon names, and the method that originally did this conversion
  (`testHMPtoHDF5`, since removed) renamed duplicate taxa before writing. Its code
  is still in the history:

```bash
git show 198c06f1:src/test/java/net/maizegenetics/dna/snp/io/BuilderFromGenotypeHDF5Test.java
```

## Hurdles before re-enabling

### Archive size and CI time

Adding a file means publishing a new `tassel_test_data` release and updating its
URL in three places: the `fetchTestData` task in `build.gradle.kts`, both test
jobs in `.github/workflows/coverage.yml`, and `.github/workflows/nightly.yml`.

Every CI run downloads and unpacks the whole archive in each of those jobs. The
archive is about 484 MB unpacked today. A full maize genome is about 2.1 billion
bases, and going by the 48 MB chromosome 8 extract, each compressed genome is
likely 500 MB or more on its own. Options include a separate, optional archive
for large files that only manual runs fetch, or the smaller substitutes
described below.

### Memory and run time

`GenomeSequenceBuilder` holds every chromosome of the reference in memory, and
several of these methods load the reference twice (once through the GVCF
builder, once directly). The test JVM is capped at 10 GB (`-Xmx10g` in
`build.gradle.kts`). `BuilderFromGVCFTest` used to be one of the slowest classes
in the suite ([#37](https://github.com/maize-genetics/tassel/issues/37)),
because each method reloaded its GVCF and reference. The active methods that use
the B73 GVCF now share one loaded copy; re-enabled methods should share their
inputs the same way.

### The GVCF tests don't check anything yet

The `BuilderFromGVCFTest` methods above write FASTA or CSV files to `tempDir/`
or print to the console. None of them asserts anything, and most of them catch
every exception and only print it. With the data in place they would pass even
if the code under test were broken. Before removing an `@Ignore`:

- Let exceptions fail the test, either by removing the `try`/`catch` or by
  calling `Assert.fail` in the `catch` as the active tests in the class do.
- Compare the output with an expected file or a known sequence.

### Region statistics are switched off

`pullB73FASTASequenceFullGVCF` writes per-region statistics, but the counters
behind them are commented out in `GVCFGenomeSequenceBuilder`:
`getPreviousRegionStats()` always returns an empty map, and `resetCounters()`
does nothing. Even with the data in place, every statistic in the output would
be `null`. Either restore the counters or drop the statistics from the test.

### The chromosome 8 extractors aren't tests

`getChr8FromB73V4` and `getChr8FromB73V3` produced the `CHR8_*.fa.gz` files that
are already in the archive. They're worth keeping as a record of where those
files came from, but they belong in a script rather than the test suite.

### The deflation methods are benchmarks

`testDeflationBySite` and `testDeflation` time Java's built-in `Deflater` on
genotype data. They don't test TASSEL code, so their `@Ignore` is labelled
`benchmark:` and they shouldn't run in CI even once the file exists.
`testDeflationBySite` also reads 10 blocks of 65,536 sites (655,360 sites in
total). Chromosome 8 of HapMap v2 is large enough, but a smaller substitute
would fail with an index error unless the loop is limited to the file's site
count.

## Smaller substitutes

- **Chromosome 8 reference instead of the full AGPv4 genome.**
  `pullW22FASTASequence`, `pullB73FASTASequence` and
  `pullB73FASTASequenceFullGVCF` only pull chromosome 8 regions, as long as the
  anchor CSVs list only chromosome 8 as their names suggest. They could use the
  `CHR8_Zea_mays.AGPv4.dna.toplevel.fa.gz` extract already in the archive, the
  same way the active GVCF tests pair it with the B73 GVCF. That would leave
  only `pullRealV4Sequence`, `pullAnchor2kFromW22` (which reads chromosome 1)
  and `getChr8FromB73V4` needing the full genome.
- **Generated regions instead of the original CSVs.** Any CSV with a header row
  and `chromosome,start,end` rows lets the anchor tests run. The output can't
  then be compared with what the original author saw, but these tests have no
  expected output yet anyway. A test could also write its own small region list
  to `tempDir/` and need no new data at all.

## Other ignored tests from the same cleanup

- **`BuilderFromVCFTest.convertToHDF5KeepDepthTest`** is ignored because of a
  bug, not missing data. Converting a VCF straight to HDF5 with `keepDepth()`
  writes no depth: `GenotypeTableBuilder.addSiteBlock` receives the depth block
  but only writes genotype calls. `convertToHDF5KeepDepthGenotypesTest` covers
  the genotype half and runs normally. Remove the `@Ignore` once `addSiteBlock`
  writes depth.
- **`ImportUtilsTest.speedFileTests`** is a timing benchmark and is labelled
  `benchmark:`.
