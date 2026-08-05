# TASSEL Modernization Plan

## Context

TASSEL (v5.2.96) is a ~245K-line Swing desktop application for genetic trait
association, evolution, and linkage analysis. The codebase is 963 Java files and
20 Kotlin files. It builds with Gradle (Kotlin DSL) on a JDK 21 toolchain and is
already Kotlin-enabled (Kotlin 2.1 plugin, mixed Java/Kotlin in a single
`src/main/java` source tree).

This plan covers three things:

1. **Build/run fix** (already applied) so the project compiles and launches.
2. **Appearance refresh** using a modern Swing Look-and-Feel (FlatLaf) — low risk,
   high visual impact, stays 100% Swing.
3. **Incremental Java → Kotlin migration** strategy — opportunistic, file-by-file,
   leaving the performance-sensitive numeric core in Java for now.

All work happens on the `modernization` branch; `main` is left untouched.

---

## 0. Build/run fix — DONE

**Problem:** `org.biojava.thirdparty:forester:1.039` pulls in the transitive
dependency `openchart:openchart:1.4.2`, which no longer exists in Maven Central,
scijava, jitpack, or any other public repo — so dependency resolution (and thus
the whole build) fails.

**Fix applied** in [`build.gradle.kts`](build.gradle.kts): a global exclusion,
since `openchart` is not referenced anywhere in TASSEL source (`forester` itself
is used only by `ArchaeopteryxPlugin.java`):

```kotlin
configurations.all {
    exclude(group = "openchart", module = "openchart")
}
```

**Verified:** `./gradlew installDist -x test` → BUILD SUCCESSFUL, and
`build/install/tassel/bin/tassel` launches the GUI cleanly.

> Follow-up worth doing later: confirm the Archaeopteryx phylogenetic-tree viewer
> plugin still works at runtime without `openchart`, or drop the plugin if it's
> dead. `openchart` is a charting sub-dependency; the tree viewer may not need it.

---

## 1. Appearance refresh — FlatLaf drop-in

### Why FlatLaf
The app currently sets **no Look-and-Feel**, so it renders in Java's dated
cross-platform "Metal" theme. [FlatLaf](https://www.formdev.com/flatlaf/) is the
de-facto modern flat L&F for Swing: actively maintained, light/dark themes,
native macOS window integration, and proper HiDPI/Retina scaling — all without
touching the 243 Swing files individually.

### Changes

**a) Add the dependency** in `build.gradle.kts`:
```kotlin
implementation("com.formdev:flatlaf:3.5.4")          // check for latest at build time
implementation("com.formdev:flatlaf-extras:3.5.4")   // optional: SVG icons, UI inspector
```

**b) Install the L&F before any Swing component is created.** The single
injection point is `TASSELMainApp.main()`
([TASSELMainApp.java:37](src/main/java/net/maizegenetics/tassel/TASSELMainApp.java))
— set it as the very first thing in `main`, before `new TASSELMainFrame()`:
```java
// macOS: use the system menu bar and proper app name
System.setProperty("apple.laf.useScreenMenuBar", "true");
System.setProperty("apple.awt.application.name", "TASSEL");
try {
    UIManager.setLookAndFeel(new com.formdev.flatlaf.FlatLightLaf());
} catch (Exception e) {
    myLogger.warn("Could not set FlatLaf look and feel", e);
}
```
Persist the user's light/dark choice via the existing `TasselPrefs` mechanism so
it survives restarts (add a `Preferences` menu toggle in `PreferencesDialog`).

**c) Make the handful of hardcoded colors/fonts theme-aware.** Only ~6 files
hardcode `new Color(...)` / `new Font(...)`. Replace literals like
`Color.lightGray` and `Color.white` with `UIManager.getColor("...")` keys so they
follow the active theme. Representative spots:
- [TASSELMainFrame.java:177](src/main/java/net/maizegenetics/tassel/TASSELMainFrame.java) — `mainPanelTextArea` forced to `Monospaced 12`
- [TASSELMainFrame.java:180](src/main/java/net/maizegenetics/tassel/TASSELMainFrame.java) — status field `Color.lightGray`
- [TASSELMainFrame.java:435,458,471](src/main/java/net/maizegenetics/tassel/TASSELMainFrame.java) — menu items forced `Color.white`

**d) HiDPI / macOS integration** comes largely for free with FlatLaf. Optionally
enable unified window decorations and a light/dark accent color.

### Explicitly out of scope for this phase
SVG icon migration (208 files use bitmap `ImageIcon`), toolbar/menu/table
restyling, and any full JavaFX/Compose rewrite. These can be follow-ups; the
drop-in above delivers the biggest visual win at minimal risk.

### Effort
~1 day. Isolated, reversible, no functional-logic changes.

---

## 2. Java → Kotlin migration — incremental / opportunistic

### Principles
- The build already compiles mixed Java/Kotlin from the same source tree, and
  Kotlin↔Java interop is seamless, so **no build changes are needed** — just add
  `.kt` files and delete the corresponding `.java`.
- **Convert file-by-file**, using IntelliJ IDEA's *Convert Java File to Kotlin*
  (or Android Studio) followed by manual cleanup of the auto-generated code
  (nullability, `!!`, `companion object`, idiomatic collections).
- **Never big-bang.** Test coverage is low (Kover minimum is 15%, and many test
  classes are excluded in `build.gradle.kts`), so wide automated conversion is
  too risky to verify.

### Priority order (safest → riskiest)
1. **New code** — write all new plugins/utilities in Kotlin (already the team's
   habit: see the 20 existing `.kt` files, e.g. the `analysis/tree/*` plugins,
   `util/ChangeLog.kt`, `matrixalgebra/*`).
2. **Leaf utilities** with no/simple dependencies — `util/`, small helpers,
   self-contained data holders. Low blast radius, easy to eyeball.
3. **GUI / Swing classes** — convert cleanly and benefit from Kotlin's
   conciseness (listeners, builders). Verify visually by launching the app.
4. **Plugin classes** (`plugindef`, `analysis/*Plugin`) — well-bounded units with
   a clear `processData` contract.
5. **Data model** (`taxa`, `phenotype`, parts of `dna`) — convert carefully;
   watch for primitive-array performance and public-API compatibility.

### Leave in Java for now
The **performance-critical numeric hot paths** — `dna/snp` genotype iteration,
`matrixalgebra`, and inner loops in `analysis` (association, distance, LD).
Kotlin handles these fine, but primitive-array micro-optimizations and
allocation behavior need per-file profiling; defer until there's test coverage
to catch regressions.

### Guardrails per conversion
- Keep the public API (class/method signatures, package) identical so Java
  callers are unaffected; use `@JvmStatic`/`@JvmField`/`@JvmOverloads` where Java
  code depends on static access or overloads.
- One class per PR/commit where practical; run `./gradlew build` + launch the app
  after each batch.
- Add or preserve tests for anything in the numeric core before touching it.

### Effort
Ongoing / continuous — not a single deliverable. Expect a long tail; the aim is
steady drift toward Kotlin, not a deadline.

---

## Verification

For the build fix (already validated) and any future refresh/migration work:

```bash
# Compile + assemble a runnable distribution (skip the heavy, partially-excluded tests)
./gradlew installDist -x test

# Launch the GUI
./build/install/tassel/bin/tassel

# Full build incl. tests (tests currently run with ignoreFailures = true)
./gradlew build
```

- **UI refresh:** launch the app and visually confirm the FlatLaf theme (and the
  light/dark toggle) render correctly; check macOS menu-bar integration and
  Retina sharpness.
- **Kotlin migration:** after each converted batch, `./gradlew build` must pass
  and the app must launch and exercise the affected screen/plugin. Prefer
  converting classes that already have tests, and add tests before converting
  numeric code.

## Suggested sequencing
1. Land the build fix (done).
2. Ship the FlatLaf drop-in (Section 1) as its own PR — small, self-contained,
   instantly visible.
3. Begin opportunistic Kotlin conversion (Section 2) as a rolling effort,
   starting with new code and leaf utilities.
