# GBS test remediation — roadmap & profiling

Date: 2026-07-10 (night). Branch: `gbs-new-test`.
Context: the `System.exit`/data/build-regression work is done and `gbsTestSmall` is green
(38 pass / 13 skip / 0 fail). This doc plans the next four items Ed asked for.

Focus decision (Ed): **prioritize GBSv2; GBSv1 (legacy) may be deleted outright.**

---

## 1. Timing report — GBS legacy vs GBSv2

`gbsTestLarge` (20 MB) completed: **54m17s, 51 tests, 44 passed, 2 failed, 5 skipped.** The 2
failures are both **legacy GBSv1** Production assertion tests (`ProductionPipelineMainTest`,
`ProductionSNPCallerPluginTest`) — we agreed not to chase these. **All six of run-#2's GBSv2 NPEs
are gone** (the `TagSchema.sql` resource fix), and 22 JUnit XMLs generated (was zero before — the
`System.exit` fix). No JVM kill.

### Per-class wall time (from `build/test-results/gbsTestLarge/*.xml`)

| sec | tier | test class |
|----:|:----:|------------|
| **1693.0** | v1 | `SeqToTBTHDF5PluginTest` |
| **815.6** | v1 | `ModifyTBTHDF5PluginTest` |
| **252.2** | v1 | `FastqToTagCountPluginTest` |
| 179.2 | v2 | `EvaluateSNPCallQualityOfPipelineTest` |
| 70.0 | v2 | `GBSSeqToTagDBPluginTest` |
| 64.1 | v2 | `SNPQualityProfilerPluginTest` |
| 43.4 | v1 | `DiscoverySNPCallerPluginTest` |
| 41.2 | v2 | `GBSv2BiologyCompareTest` |
| 33.0 | v2 | `DiscoverySNPCallerPluginV2Test` |
| 28.4 | v1 | `ProductionPipelineMainTest` |
| 16.3 | v1 | `ProductionSNPCallerPluginTest` |
| 3.0 | v1 | `MergeMultipleTagCountPluginTest` |
| <1 | — | the remaining 10 classes |

### Headline
- **Legacy GBSv1 = 2852.6 s (~47.5 min, 88% of the run). GBSv2 = 389.6 s (~6.5 min, 12%).**
- The **top 3 sinks are all legacy HDF5 tests** — `SeqToTBTHDF5` (28 min!) + `ModifyTBTHDF5`
  (13.6 min) + `FastqToTagCount` (4 min) = **~46 of the 54 minutes.**
- **⇒ Deleting GBSv1 (§2), which is planned anyway, cuts the suite 54 min → ~6.5 min by itself.**
  It is by far the biggest test-speed lever — bigger than any fixture/parallelism change.
- Within v2, the long poles are `Evaluate` (179 s) and `GBSSeqToTagDB` (70 s) — the DB-build
  tests, which §3.2 (shared fixture) and §4 (algorithm speedups) target.

---

## 2. Plan — delete GBS legacy (GBSv1) code + tests

### Key dependency finding (makes this tractable)
The **only** legacy `analysis/gbs` class the active v2/RNA code depends on is **`Barcode`**
(imported by `BarcodeTrie`, `GBSUtils`, `GBSSeqToTagDBPlugin`, `ProductionSNPCallerPluginV2`).
`Barcode` itself imports nothing from the legacy package and only references `BaseEncoder`
(which lives in `net.maizegenetics.dna`, shared infra — stays). Legacy `SmithWaterman` is
**not** used by v2.

So the deletion reduces to: **preserve `Barcode`, delete the rest of `analysis/gbs/*` (legacy),
and delete their tests.**

### Proposed order (each step compiles + `gbsTestSmall` green before the next)
1. **Move `Barcode` out of the legacy package** into `analysis/gbs/v2` (or a neutral
   `analysis/gbs/common`), update the 4 v2 imports. Keeps the one live dependency and lets us
   delete the legacy package wholesale.
   - Check first: does `Barcode` reference other same-package legacy types
     (`ParseBarcodeRead`, `ReadBarcodeResult`)? Grep showed only `BaseEncoder`. Confirm before moving.
2. **Delete legacy test classes** (`src/test/.../analysis/gbs/*Test.java`, non-v2): `FastqToTagCountPluginTest`,
   `MergeMultipleTagCountPluginTest`, `TagCountToFastqPluginTest`, `SAMConverterPluginTest`,
   `SeqToTBTHDF5PluginTest`, `ModifyTBTHDF5PluginTest`, `DiscoverySNPCallerPluginTest`,
   `ProductionSNPCallerPluginTest`, `ProductionPipelineMainTest`, `ParseBarcodeReadTest`.
   - These are the tests whose `System.exit` I already patched; that patch effort is sunk but
     harmless — the classes go away with the tests.
3. **Delete legacy main classes** `analysis/gbs/*.java` (the 53 non-`Barcode` files), plus the
   `pana/` subpackage if it is also legacy-only (verify no v2 use — earlier analysis put pana in
   the CLI-only/unreachable tier).
   - Watch for shared utility types that v2 or other packages import (grep each class name across
     `src/main` before deleting). Candidates that might be referenced elsewhere: `FilePacking`,
     `TagsByTaxa*`, `BaseEncoder` (already in `dna/`, safe). Delete in dependency order, compiling
     between batches.
4. **Remove the legacy include from the gbs test tasks**: change the `gbsTest*` `include(...)` from
   `**/analysis/gbs/*Test.class` + `v2` to **v2 (+ repgen) only**, since legacy tests are gone.
5. **Update the GUI**: `TASSELMainFrame.getGBSMenu()` (the struck-through legacy menu) and its
   `addMenuItemDeprecated` wiring can be removed once the legacy plugins are gone.

### Risks / gates
- The `pana/` subpackage and CLI-only legacy plugins (`TagAgainstAnchor`, `AnnotateTOPM`, …) may
  have their own web of intra-package deps; delete as a connected component, not file-by-file.
- Anything still referenced by `TasselPipeline` reflection-by-name would only fail at runtime, not
  compile — grep `TasselPipeline` and the GUI for the deprecated class names before removing.
- Keep this a **separate PR/commit series** from the test-remediation work.

---

## 3. Plan — speed up the remaining GBSv2 tests

Highest-impact first (the first two are the big ones):

1. **Make the 200 KB dataset the default `gbsTest` target.** Already wired: `gbsTestSmall`
   (~seconds/minutes) vs `gbsTestLarge` (~1h). Point CI / the default at small; keep large as a
   nightly/opt-in. Biggest practical win for day-to-day.
2. **Extend the `@BeforeClass` shared-DB fixture** (done for `GBSSeqToTagDBPluginTest`, 6→1) to
   **`EvaluateSNPCallQualityOfPipelineTest`**: build DB + export + SAM-import once in
   `@BeforeClass`; the two pipeline methods differ only by MAF in the Discovery step (2 full
   rebuilds → 1). This is the biggest single v2 build cost on the large profile.
3. **Isolate per-test DB paths** (`@Rule TemporaryFolder` or a unique db filename per test) so
   tests never collide on `GBS_GBS2DB_FILE`. Removes the last ordering coupling and **unlocks
   parallelism**.
4. **Enable parallel forks** once (3) lands: `maxParallelForks` on the test task, sized to
   RAM (each fork is `-Xmx10g`). With per-test DBs there are no shared-file collisions.
5. **Trim redundant rebuilds in the mutating tests** where the assertion doesn't actually need a
   fresh parse (audit `testKeepOldData` / `AppendTest` — they legitimately rebuild, so likely leave).
6. **Drop kover instrumentation for routine local runs** (the `-javaagent` kover agent adds
   overhead); keep it only for the coverage task.

Expected: large profile from ~1h → a few minutes once small is the default and the DB is shared +
parallel; keep large as a full-fidelity nightly.

---

## 4. Read-mapping algorithm — optimization candidates (PROFILED)

**JFR profile taken** of `GBSSeqToTagDBPluginTest` on the 20 MB fastqs (`settings=profile`,
`/tmp/gbs-21473.jfr`). The profile **overturned the static-analysis guess** — the win is the
barcode trie, and `seqDifferences` is not even hot in this phase.

### Hot leaf methods (execution samples, DB-build phase)
| samples | method | note |
|---:|---|---|
| **2429** | `BarcodeTrie.longestPrefix` → `HashMap`/`ArrayList` | **~40% of all samples — #1** |
| 775 | `ConcurrentHashMap$Traverser.advance` | tag-map iteration (stats/purge) |
| 739 | `BufferedReader.implReadLine` | fastq I/O |
| 734 | `BaseEncoder.getLongFromSeq` | via `AbstractTag.getLongArrayFromSeq` (tag build) |
| 476 | `GBSSeqToTagDBPlugin.removeSecondCutSiteIndexOf` | per-read `substring` |
| 424 | `AbstractTag.getLongArrayFromSeq` | 2-bit encode |
| 274 | `BaseEncoder.getFirstLowQualityPos` | quality scan |
Top allocated types: `byte[]`, **`Object[]` (the per-char ArrayList, see 4a)**, `String`
(substrings). Real GC pressure (16.6k allocation samples).

### ✅ 4a DONE (2026-07-11) — measured result
Replaced `crawl.containsKey(ch)` with a null-check on the already-computed `child` in
`longestPrefix`. Re-profiled the same 20 MB `GBSSeqToTagDBPluginTest` workload:
- **Total execution samples 6749 → 2777 (~59% less CPU).**
- **`BarcodeTrie.longestPrefix` 2673 → 389 samples (−85%);** `HashMap.getNode` 2429 → 32.
- All 5 tests pass on 20 MB incl. tag-count assertions (functional equivalence).
- **Residual (389):** now the `result += ch` string concat + `barcodeInformation.get(result)`
  (`String.hashCode` newly visible). **Follow-up:** stop building `result`; store the `Barcode`
  on the terminal `TrieNode` and return it directly — removes the per-read String + final map get.
- Also moved (Ed's request) `SmithWaterman` → `analysis/gbs/v2` (was legacy `gbs`, self-contained,
  only user was legacy `TagMatchFinder` which got an import); preserves it through the GBSv1 delete.

### 4a (original note). `BarcodeTrie.TrieNode.containsKey(char)` → array null-check  ★★ THE win (profiled #1)
`v2/BarcodeTrie.java:211`. Called for every char of every read (millions of times), it currently:
```java
List followers = new ArrayList();                 // fresh allocation EVERY call
for (TrieNode x : children) if (x!=null) followers.add(x.character); // scan 26 slots, autobox
return followers.contains(c);                      // linear scan, autobox c
```
It is semantically identical to `children[c-'A'] != null` — an O(1) array check with **zero
allocation**. Also `longestPrefix` calls both `getNode(ch)` and `containsKey(ch)` (double work);
collapse to one null-checked `getNode`. Expected: gut most of the ~40% #1 hot spot **and** a large
chunk of the allocation/GC pressure. Trivial, safe, unit-testable (same barcode assignments).

### 4b. Per-read `String.substring` allocation (`removeSecondCutSiteIndexOf`, encode) — profiled #4/#6
`GBSSeqToTagDBPlugin.java:336/371` and `AbstractTag.getLongArrayFromSeq`. `String` is the #3
allocated type. Operate on the read `String` + start offset with `indexOf(str,fromIndex)` /
`regionMatches`, materializing only the final tag substring; encode chunks with index math instead
of `substring`. Do after 4a and re-profile.

### 4c. Tag-map churn — `ConcurrentHashMap$Traverser.advance` (profiled #2)
775 samples iterating the tag map (`calcTagMapStats`, `removeTagsWithoutReplication`,
`removeTagByCount`). Audit how often the full map is walked; the get-then-put at
`GBSSeqToTagDBPlugin.java:341-347` is also racy across the file-parallel streams — switch to
`merge`/`compute`. Secondary to 4a.

### 4d. `BaseEncoder.seqDifferences` — bit loop → POPCNT (NOT hot here; alignment phase only)
`dna/BaseEncoder.java:463/488` still worth doing (32-iteration loop →
`Long.bitCount((diff|(diff>>>1)) & 0x5555…)`), but it did **not** appear in the DB-build profile —
it's used in tag **alignment** (`DiscoverySNPCallerPluginV2`). Profile `Evaluate`/Discovery
separately before investing; likely relevant to the 179 s `Evaluate` class.

### Profiling method (reproducible)
```
JAVA_TOOL_OPTIONS="-XX:StartFlightRecording=filename=/tmp/gbs-%p.jfr,settings=profile,dumponexit=true" \
  ./gradlew gbsTestLarge --tests "*GBSSeqToTagDBPluginTest"
# largest /tmp/gbs-<pid>.jfr is the test worker; analyze with `jfr print --events jdk.ExecutionSample`
```

### Recommended order to try tomorrow
1. **4a** (BarcodeTrie) — biggest, trivial, safe. Add a unit test asserting identical barcode
   assignment on a fixture, then re-profile to measure the drop.
2. **4b** (substring) if the re-profile still shows `String`/GC high.
3. Profile the **alignment** phase (`Evaluate`, 179 s) → decide on **4d** (POPCNT) there.

---

## 5. Alignment-phase profile + `seqDifferences` (2026-07-11)

**Done:** `BaseEncoder.seqDifferences` family rewritten to the correct 2-bit-parallel popcount
(`(diff | diff>>>1) & 0x5555…`, then `Long.bitCount`) with randomized equivalence tests — commit
`ef86af8f`. Note it is **legacy-only today** (callers `TagMatchFinder`, `ParseBarcodeRead`); kept in
`dna/BaseEncoder` as a correct/fast primitive for the aligner idea below.

### Alignment-phase profile (Evaluate `pipelineIncludingInvariantSites`, 20 MB, 2553 samples)
Per-sample attribution:
| share | area |
|---:|---|
| **59.2%** | **BioJava ClustalW MSA** (`Alignments.getMultipleSequenceAlignment` in `DiscoverySNPCallerPluginV2.alignTags`) |
| 20.3% | GBS DB build (fastq parse / BarcodeTrie) |
| 8.9% | Discovery (non-biojava) |
| 7.4% | TagDataSQLite I/O |

**Key finding inside the 59%:** the single hottest leaf (~28% of the *whole* pipeline) is
`SimpleSubstitutionMatrix.getValue → getIndexOfCompound → ArrayList.indexOf` — BioJava does a
**linear list scan to map each nucleotide compound to its matrix row/col index, on every DP cell**.
So roughly half the alignment cost is this lookup, not the DP recursion (`setScorePoint`,
`getSubstitutionScoreVector` ≈ the other ~31%).

### Two levers (recommendation)
- **(L1) Low-risk, same output — fix the substitution-matrix lookup.** Give the aligner a
  substitution matrix / compound-set with O(1) compound→index (array/`EnumMap`-style) instead of
  `List.indexOf`, or a DNA-specialized scorer. Produces **identical scores → identical alignments →
  identical SNP calls**, so it's verifiable against `ExpectedResults` with low risk, and could
  reclaim a large chunk of the ~28%. Feasibility depends on what BioJava's
  `Alignments.getMultipleSequenceAlignment(lst)` lets us inject (it's currently called with all
  defaults in `alignTags`). **Recommended first step.**
- **(L2) High-value, gated — 2-bit fast path replacing ClustalW for the common case.** Tags at a cut
  site that are equal-length and indel-free need only columnar/Hamming comparison (the new
  `seqDifferences` is the primitive); fall back to ClustalW only when lengths differ / indels are
  needed. Requires: (a) instrument `alignTags` to measure how often a cut site's tags are uniform
  length (fast-path hit rate), and (b) a SNP-call equivalence harness vs `ExpectedResults` before any
  swap. Bigger win, higher risk — do after L1 and only with sign-off.

### Next action
Prototype **L1** (custom O(1) substitution scorer for `alignTags`) and re-profile + diff SNP calls;
report before touching L2.

### ✅ L1 DONE (2026-07-11) — commit `513ef9fb`
`IndexedNucleotideSubstitutionMatrix` wraps nuc-4.4 with a precomputed dense score table +
`IdentityHashMap` compound index (IdentityHashMap avoids `NucleotideCompound.hashCode()`/`equals()`,
which are `toString()`-based — a `HashMap` first attempt just moved the cost into `String.hashCode`).
Passed to `Alignments.getMultipleSequenceAlignment(lst, FAST_NUC_4_4)`.
- **Correctness:** scores byte-identical → SNP calls unchanged; `DiscoverySNPCallerPluginV2Test`
  (incl. `testAlignTags`, `testFullSNPCaller`) and `gbsTestSmall` green.
- **Measured (Evaluate, 20 MB, JFR):** linear compound lookup 709→0; alignment 59%→49% of pipeline;
  total samples 2553→2046 (~20% overall, ~33% off the alignment phase).
- **Honest caveat:** end-to-end test **wall time flat** (38s→37s) — alignment CPU is not this
  test's wall-clock bottleneck. Remaining alignment cost is BioJava's O(n²) profile-profile DP
  (`AlignerHelper.setScorePoint`, the per-cell compound loop) — only removable by **L2** (replace the
  aligner), which is high-risk (SNP-call equivalence) for an uncertain wall-time payoff on this test.
- **Recommendation:** keep L1 (correct, low-risk CPU win). **Hold L2** unless production profiling on
  a many-cut-site dataset shows alignment dominates wall time.

---

## Suggested tomorrow order
1. Fill in §1 timing + §4 JFR profile (data-gathering, machine now free).
2. Try **4a** (`seqDifferences` POPCNT) — smallest, safest, add an equivalence unit test.
3. Try **4b** (substring allocation) if the profile shows GC/`substring` dominates.
4. Start the **legacy-deletion** series (§2) as its own PR, `Barcode` move first.
5. Make **200 KB the default** test target (§3.1) + extend the fixture to Evaluate (§3.2).
