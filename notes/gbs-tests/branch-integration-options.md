# Integrating the two GBS test-remediation branches

Two overlapping efforts must combine for the next version:

- **Ed** — `gbs-new-test` (test-infra fixes + GBSv2 read-mapping optimizations) → `remove-gbs-legacy`
  (deletes the legacy GBSv1 pipeline, PanA, and the dead `RNADeMultiPlexSeqToDBPlugin`; moves
  `Barcode` to `gbs/v2`). Verified green: `gbsTestSmall` 1m42s, `gbsTestLarge` 6m30s (was 54 min).
- **Colleague** — `test-remediation` (12 commits): HDF5/jhdf5 native support, an independent
  ~220-line `build.gradle.kts` overhaul, coverage excludes, v2 sim data, GUI codegen annotations,
  "purge hard-coded paths", "enforce statistical correctness".

They conflict on ~13 files (build.gradle.kts, several GBS plugins + tests) — including a
**delete-vs-modify** on `RNADeMultiPlexSeqToDBPlugin` and on legacy classes/tests (e.g.
`ProductionSNPCallerPlugin`, `DiscoverySNPCallerPlugin`, `FastqToTagCountPlugin`,
`ProductionPipelineMainTest`) that Ed deletes but the colleague edited.

## Guiding fact
Ed's removals are the correct baseline: the legacy GBSv1 pipeline, PanA, and the never-deployed RNA
plugin are dead/deprecated. **You should not remediate, HDF5-fix, or annotate code that is being
deleted.** A large share of the overlap is exactly that — so the *order* of integration decides how
much reconciliation work is wasted.

## Options

### A. Ed's branches first, then colleague rebases onto the result  ★ recommended
Merge `remove-gbs-legacy` (which carries `gbs-new-test`) into the release/integration branch first;
the colleague then rebases/re-applies `test-remediation` on top.
- **Strengths:** deletions land first, so the codebase the colleague reconciles against is much
  smaller; every colleague edit to a now-deleted legacy/RNA file simply drops out (correctly —
  that work was remediating code slated for removal). Foundational change (what code exists) settles
  before incremental change (how it's tested/built). Clear, correct baseline.
- **Weaknesses:** the colleague must redo/rebase and will lose edits to deleted files (intended);
  requires sequencing/cooperation. Any of the colleague's *non-legacy* work (jhdf5, sim data,
  non-GBS build/coverage) still needs a genuine merge.

### B. Colleague's `test-remediation` first, then Ed rebases on top
- **Strengths:** the colleague's HDF5/build work lands first with no rebase for them.
- **Weaknesses:** backwards — Ed then re-deletes files the colleague just modified (every
  delete-vs-modify resolved manually), and the colleague's effort spent remediating legacy tests is
  wasted anyway once Ed deletes them. Maximizes conflict and rework. Not recommended.

### C. Reconcile both into one combined branch in a single pass
- **Strengths:** one clean integrated result; conflicts resolved once.
- **Weaknesses:** highest-effort manual merge across ~13 files + two divergent build.gradle
  overhauls; needs deep knowledge of *both* efforts at once; easy to silently drop intent from
  either side. Overkill given (A) makes most conflicts vanish.

### D. Cherry-pick only the non-overlapping pieces from one branch onto the other
- **Strengths:** avoids wholesale conflict; keep only what's wanted.
- **Weaknesses:** tedious and error-prone; hard to be sure nothing is missed, especially in the two
  build.gradle overhauls. Fine as a supplement to (A), not a primary strategy.

## Recommendation
**Option A — Ed's deletions first, colleague rebases on top.** It puts the correct, foundational
change (removing dead code) first, which dissolves most of the overlap instead of forcing it to be
resolved by hand, and prevents wasted effort remediating code that is being deleted. Concretely:
1. Land `gbs-new-test`, then `remove-gbs-legacy`, on the integration branch.
2. Colleague rebases `test-remediation`; edits to deleted legacy/RNA files drop out.
3. Reconcile only what remains: the two `build.gradle.kts` overhauls, jhdf5/HDF5 support, sim data,
   and coverage excludes — for GBSv2/kept code only.
