# TASSEL Web UI Architecture — One UI, Two Compute Tiers

**Audience:** developers planning the post-Swing front-end.
**Purpose:** capture the target architecture — a single web UI that runs small jobs
**client-side in WebAssembly** and delegates big jobs to the **same JVM analysis core
server-side** — plus the migration path to get there without disrupting `rTASSEL`.

---

## Why this shape

Three constraints drive the design:

1. **Teaching is the dominant use case.** ~25,000 students/year use the TASSEL GUI to
   learn GWAS, on mixed OSes and browsers, often on lab machines. The biggest UX win
   is **zero-install, runs in a browser** — ideally with **no per-student server cost**.
2. **`rTASSEL` depends on the JVM core.** The analysis engine must stay a JVM library,
   callable in-process from R (via rJava). That surface must not change.
3. **Swing has a low ceiling.** Modernizing "each main graphic component" in Swing is
   diminishing returns; any real refresh is a UI rewrite regardless of target toolkit.

The insight that reconciles these: **teaching-scale GWAS fits inside the browser**, while
research-scale does not. So build **one web UI** and route each job to the cheapest place
it can correctly run — client Wasm for small/tutorial data, server JVM for large data —
over the **same core algorithms**.

---

## Guiding principle: separate the headless analysis core first

This is the prerequisite for everything below, and TASSEL is already ~90% there
(`TasselPipeline` runs the full engine headless as a CLI).

- Carve the analysis engine into a **UI-free JVM library module** (no `javax.swing`
  dependency anywhere in it).
- Give it **two consumers**:
  - **in-process** — `rTASSEL` (unchanged) and any future JVM desktop shell;
  - **service API** — a thin server that exposes the core to the web UI.

Done well, this de-risks the entire UI effort and is valuable on its own even if the
front-end decision changes.

---

## The architecture

```
                         ┌───────────────────────────────────────┐
                         │              Web UI (one codebase)     │
                         │   screens, wizards, tables, plots      │
                         │   Compose Multiplatform/Wasm  OR  web  │
                         └───────────────┬───────────────────────┘
                                         │  submit job(spec, data)
                                         ▼
                         ┌───────────────────────────────────────┐
                         │            Compute Router              │
                         │  decides client vs server per job,     │
                         │  based on: algorithm, data size,       │
                         │  memory estimate, server availability  │
                         └───────┬───────────────────────┬───────┘
                 small/teaching  │                       │  large/research
                                 ▼                       ▼
        ┌────────────────────────────────┐   ┌────────────────────────────────┐
        │   Client Compute Engine (Wasm) │   │        Server API (thin)        │
        │  ported teaching algorithms:   │   │   REST/gRPC over the JVM core   │
        │  GLM, MLM(Q+K), kinship, PCA,  │   │   job queue, result storage     │
        │  LD, basic filters             │   └───────────────┬────────────────┘
        │  (SIMD + Wasm threads)         │                   ▼
        └────────────────────────────────┘   ┌────────────────────────────────┐
                                              │     JVM Analysis Core (shared)  │
                                              │  full TASSEL engine, native     │
                                              │  BLAS/HDF5/SQLite, big-data      │
                                              │  ← also used in-process by      │
                                              │    rTASSEL and desktop shell    │
                                              └────────────────────────────────┘
```

**Key property:** the client Wasm engine and the server JVM core implement the **same
algorithms** and must produce **numerically identical results** (within tolerance) — so a
student's browser run and a researcher's server run agree. This is enforced by shared
correctness tests (see `TEST_REMEDIATION_PLAN.md` — the statistical must-pass group
becomes the cross-tier oracle).

---

## Compute routing policy (starting heuristic)

| Job | Client (Wasm) | Server (JVM core) |
|---|---|---|
| GLM association | Any teaching size; up to ~mid-size | Large marker sets, memory-bound |
| MLM (Q+K): kinship + per-marker | Small n (hundreds of taxa) — O(n³) eig is cheap | Thousands of taxa (O(n³) + RAM) |
| Kinship / distance matrices | Small–medium n | Large n, HDF5-backed |
| PCA | Small–medium | Large |
| LD, basic filters, format conversion | Yes | If data is huge |
| GBS pipeline, HDF5/DB-backed, imputation | No (native deps) | **Always server** |

Route by a cheap up-front estimate (taxa × sites × algorithm cost + memory estimate vs
the ~4 GB Wasm ceiling), with a manual override and graceful fallback to server if a
client job exceeds budget.

---

## File & data handling (see browser capability notes)

- **Load (all browsers):** `<input type=file>` / drag-and-drop, with `File.stream()` for
  large genotype files. Covers the tutorial "open genotype + phenotype" flow everywhere.
- **Working data (all browsers):** OPFS (sandboxed, fast, worker-sync) for staging
  converted matrices / a wasm-SQLite scratch DB / intermediate results.
- **Project folder read/write (Chromium only):** File System Access API for a
  desktop-like "open a project directory" experience — treat as progressive enhancement,
  not a baseline requirement.
- **Results (all browsers):** Blob download; `showSaveFilePicker` save-in-place on Chromium.
- **Server jobs:** stream/upload the input (or reference server-side storage); return a
  job id, poll/stream progress, download results.

---

## Deployment models (all from the same UI codebase)

1. **Static teaching site — client-only, no backend.** Pure Wasm compute on tutorial
   data. Zero install for students *and* zero hosting/compute cost — a static host (CDN)
   serves it. This is the headline win for the 25k-student use case.
2. **Hybrid — client + backend.** Adds the server JVM core for research-scale jobs; the
   router sends big jobs there. For labs/institutions running real datasets.
3. **Desktop shell — in-process core.** If a native desktop app is still wanted (power
   users, offline), a Compose Multiplatform/JVM shell can embed the core directly (no
   IPC), reusing the same UI screens.

---

## What gets built vs reused

**Reused (unchanged):**
- The JVM analysis core and its algorithms.
- `rTASSEL`'s integration surface.
- The statistical-correctness test suite (becomes the cross-tier numerical oracle).

**Built:**
- The headless-core module extraction.
- The web UI (screens, wizards, data tables, plots — plots replace the current
  Swing/JFreeChart charts).
- The **client Wasm compute engine** — a port of the *teaching subset* of algorithms
  (GLM, MLM, kinship, PCA, LD, basic filters). Candidate toolchain: **Kotlin/Wasm**
  (aligns with the ongoing Java→Kotlin migration) or a JS+Wasm numeric stack.
- The thin **server API** over the JVM core (job submit/status/result).
- The **compute router**.

**The hard part — native dependencies.** The full JVM core uses native **BLAS
(OpenBLAS), jhdf5, native SQLite, snappy**; these do **not** cross to Wasm. So the client
engine cannot be "the JVM core compiled to Wasm" — it is a *reimplementation of the
teaching algorithms* on a Wasm-friendly linear-algebra stack. The full core stays on the
server for anything native-backed.

---

## Phased roadmap

| Phase | Goal | Notes |
|---|---|---|
| **0 (done)** | FlatLaf stopgap | Current Swing GUI looks modern for this year's cohort |
| **1** | Extract the headless JVM core module | Protects `rTASSEL`; unblocks all UI options |
| **2** | Server API + web UI shell; **all compute server-side first** | Fastest path to a modern web UI; proves the screens/plots without the Wasm engine |
| **3** | Client Wasm engine for the teaching algorithm set + the router | Port GLM/MLM/kinship/PCA/LD; validate numerically against the server via shared tests |
| **4** | Static teaching deployment (client-only) | Zero-install, zero-cost student build |
| **(opt)** | Compose MP desktop shell (in-process core) | Only if a native desktop app is still wanted |

Sequencing rationale: getting a **modern web UI with server compute (Phase 2)** delivers
value early and de-risks the UI itself; the **client Wasm engine (Phase 3)** is the
optimization that unlocks the zero-cost teaching deployment, and can be validated against
the already-running server implementation.

---

## Risks & open questions

- **Compose Web/Wasm maturity** for heavy, interactive data grids and genomics plots
  (Manhattan/QQ) — de-risk with a Phase-2 spike on the genotype-table + Manhattan-plot
  screen before committing.
- **Which algorithms make the "teaching set"** ported to Wasm — needs a short scoping
  pass against the actual tutorials/courses.
- **Client engine toolchain** — Kotlin/Wasm vs a JS+Wasm numeric stack; prototype both on
  MLM before deciding.
- **Cross-tier numerical agreement** — needs the statistical test suite fixed first
  (`TEST_REMEDIATION_PLAN.md`) so it can serve as the shared oracle.
- **Cross-origin isolation** (COOP/COEP headers) required for Wasm threads/SharedArrayBuffer
  — a hosting/config constraint for the static site.
- **Charting** — replacing JFreeChart with a web/Compose charting approach.

---

## One-line summary

Separate the JVM core once; put one modern web UI on top; **run small/teaching GWAS in the
browser (Wasm) and send big/research GWAS to the same core on a server** — with a static,
zero-cost browser build for the 25,000 students and the full engine intact for `rTASSEL`.
