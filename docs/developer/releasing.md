# Releasing

TASSEL releases are largely **automated through GitHub Actions**. Merging code to
`main` builds the application, publishes a GitHub release with installers and a
standalone distribution, and redeploys the documentation site.

## Branching model

TASSEL uses a `develop` integration branch and reserves `main` for released
code:

```
feature/*  --PR-->  develop  --promotion PR-->  main  --> release
hotfix/*   ------------------------PR------------------->  main  --> release (then synced to develop)
docs/*     ------------------------PR------------------->  main  --> site deploy only
```

* `develop` — the default branch and integration target for all normal work
  (features, enhancements, non-urgent bug fixes). See
  [Contributing](../CONTRIBUTING.md).
* `main` — released code only. A push of code here triggers the release
  automation described below.

A release is cut by opening a **promotion PR** from `develop` into `main`.
Because `main` is where the automation fires, everything on this page keys off
that promotion merge.

Documentation is the one thing that reaches `main` without producing a release —
see [Documentation-only changes](#documentation-only-changes).

## Versioning

The project version is declared in `build.gradle.kts`:

```kotlin
version = "5.2.98"
```

The same value is duplicated in `TASSELMainFrame.version`, which is what the GUI and
pipeline banner report. Bump both together.

Bump this value on `develop` before opening the promotion PR. You can confirm
the current version with:

```bash
./gradlew printVersion
```

The build workflow reads this version to name the GitHub release
(`v<version>`), the standalone archives, and the installer bundles.

## What happens on merge to `main`

Two workflows chain together automatically.

### 1. Build & publish release (`jdeploy.yml`)

Triggered on push to `main`, this "jDeploy CI with Gradle" workflow:

1. Sets up JDK 21 and reads the version via `./gradlew printVersion`.
2. Builds the app: `./gradlew clean build -x test -x koverVerify`.
3. Assembles the **standalone distribution** in `dist/tassel-5-standalone/` —
   `sTASSEL.jar`, the `lib/` dependencies, and the launcher scripts
   (`start_tassel.pl`, `run_pipeline.pl`, `run_pipeline.bat`) with their paths
   rewritten for the flattened layout.
4. Archives the standalone as both `.zip` and `.tar.gz`.
5. Extracts changelog content from the latest merged PR description and renders
   release notes from `.github/release_template.md`.
6. Creates a **GitHub Release** tagged `v<version>` and uploads the standalone
   archives.
7. Runs [jDeploy](https://www.jdeploy.com/) to build native **installer bundles**
   (Linux x64, macOS arm64/x64, Windows x64) and attaches them to the release.
8. Updates the download links in the release notes and in
   `docs/overrides/partials/featured_downloads.html`, then commits the updated
   `docs/changelog.md` and download links back to the repository.

### 2. Deploy the documentation site (`deploy_project_site.yml`)

After the build workflow completes on `main`, the "Site Builder and Deployment"
workflow builds this MkDocs site with `mkdocs gh-deploy --force`, publishing it to
GitHub Pages. It also runs directly on any push to `main` that touches `docs/**`
or `mkdocs.yml`, which is what publishes documentation-only merges.

## Documentation-only changes

A push to `main` where **every** changed file is documentation (`docs/**`, any
`*.md`, or `mkdocs.yml`) skips the release automation entirely: both
`jdeploy.yml` and `run_publish_maven.yml` carry a `paths-ignore` filter for those
paths. No application build, no GitHub release, no Maven Central publish — only
the site deploy above.

Two consequences worth knowing:

* The filter skips a run only when *all* changed files are documentation. A PR
  mixing docs with a workflow or source edit still triggers the full release
  path, which is why the documentation track refuses mixed changes.
* Because these commits land on `main` first, they need to come back down to
  `develop`. That is what
  [`sync_main_to_develop.yml`](../CONTRIBUTING.md#keeping-develop-in-sync-with-main)
  is for.

See [Documentation track](../CONTRIBUTING.md#documentation-track) for how to open
one of these PRs.

## Publishing to Maven Central (`run_publish_maven.yml`)

Publishing the library to Maven Central runs on pushes of code to `main` and on
version tags, and can also be triggered manually. Documentation-only pushes are
excluded, and nightly `dev-*` tags never publish. It relies on the
`maven-publish`, `signing`, `shadow`, and `jreleaser` configuration in
`build.gradle.kts`:

- The publication artifact id is `tassel` under group `net.maizegenetics`.
- Artifacts are GPG-signed. The signing key and passphrase are supplied via the
  `JRELEASER_GPG_SECRET_KEY` and `JRELEASER_GPG_PASSPHRASE` environment
  variables/secrets.
- API docs are generated with Dokka (`dokkaJar`) and attached to the publication
  alongside the sources JAR.
- JReleaser deploys the staged repository (`build/staging-deploy`) to Sonatype /
  Maven Central.

### Published artifacts

Every release must publish four JARs:

| Artifact | Produced by |
| -------- | ----------- |
| `tassel-<version>.jar` | `jar` (also written locally as `sTASSEL.jar`) |
| `tassel-<version>-jar-with-dependencies.jar` | `shadowJar` |
| `tassel-<version>-sources.jar` | `sourcesJar` |
| `tassel-<version>-javadoc.jar` | `dokkaJar` |

The fat JAR carries the `jar-with-dependencies` classifier that
`maven-assembly-plugin` used through 5.2.96, so consumers of that coordinate are
unaffected by the move to Gradle. It is deliberately excluded from `assemble`, so
`./gradlew build` does not pay for it; only `publish` builds it.

### Artifact verification

Maven Central is immutable — a bad upload can only be corrected by cutting a new
version. Release 5.2.97 shipped without a fat JAR, and its sources and javadoc JARs
were copies of the main JAR, because every `Jar` task was writing to the same
`sTASSEL.jar` file.

To prevent a repeat, `.github/scripts/verify-staged-artifacts.sh` runs between
staging and deployment and fails the workflow unless all four JARs exist, are
signed, have distinct checksums, and contain what their classifier claims. Run it
yourself against a local staging directory:

```bash
./gradlew clean publish
.github/scripts/verify-staged-artifacts.sh "$(./gradlew printVersion -q | tail -n1)"
```

## Nightly dev builds

Every night, `.github/workflows/nightly.yml` builds `develop` and publishes an
unstable standalone **prerelease** tagged `dev-YYYYMMDD` — but only when
`develop` has new commits since the last nightly. These builds are for testing
only and never publish to Maven Central. The nightly also runs the full test
suite, so a red nightly means `develop` is broken.

## Cutting a release: checklist

Releases are promoted out of `develop`:

1. Ensure `develop` is green (the statistics gate passes).
2. On `develop`, update `version` in `build.gradle.kts`.
3. Open a **promotion PR** from `develop` into `main`. Make sure the PR
   description contains the changelog block (between the
   `<!-- BEGIN CHANGELOG -->` / `<!-- END CHANGELOG -->` markers) so release
   notes and `docs/changelog.md` are generated correctly — the automation resolves
   the PR from the commit being built.
4. Merge the promotion PR. Pushing to `main` runs the build/release and
   site-deploy workflows automatically.
5. Verify the new GitHub Release, its attached standalone archives and
   installers, and the updated [Version History](../changelog.md).
6. If publishing to Maven Central, run the `run_publish_maven.yml` workflow and
   confirm the artifacts appear on Central.

## Hotfixes

For a critical bug in an already-released version, do **not** wait for the
normal `develop` cycle:

1. Branch from `main`: `git switch main && git pull && git switch -c hotfix/short-description`.
2. Bump the **patch** portion of `version` in `build.gradle.kts` and open a PR
   targeting `main` using the hotfix PR template.
3. After the PR merges and the release publishes, the fix has to reach `develop`
   too, or the next promotion will reintroduce the bug. You do not need to
   cherry-pick it: an automated `main` → `develop`
   [back-merge PR](../CONTRIBUTING.md#keeping-develop-in-sync-with-main) is
   opened for you. Confirming that PR merges is part of the hotfix checklist,
   and when resolving its conflicts you keep **`develop`'s** `version`.

## Local dry runs

You can reproduce the key build steps locally to sanity-check a release before
merging:

```bash
# Build exactly as the release workflow does
./gradlew clean build -x test -x koverVerify

# Confirm the version string
./gradlew printVersion

# Build the docs site locally (requires mkdocs-material)
mkdocs serve
```
