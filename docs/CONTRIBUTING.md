# Contributing to TASSEL
Thank you for your interest in contributing to TASSEL! We welcome contributions from anyone, and are grateful for even the smallest of fixes!

TASSEL is developed on GitHub at [maize-genetics/tassel](https://github.com/maize-genetics/tassel). This page covers everything you need to contribute: setting up your environment, reporting issues, and the mechanics of the Git workflow, pull requests, testing, and code review.

## Table of Contents
* [Code of Conduct](#code-of-conduct)
* [Getting Started](#getting-started)
* [How to Contribute](#how-to-contribute)
     * [Reporting Bugs](#reporting-bugs)
     * [Suggesting Enhancements and New Features](#suggesting-enhancements-and-new-features)
     * [Submitting Code Changes](#submitting-code-changes)
* [The Git Workflow](#the-git-workflow)
     * [Which track am I on?](#which-track-am-i-on)
     * [Feature track (the default)](#feature-track-the-default)
     * [Hotfix track](#hotfix-track)
     * [Documentation track](#documentation-track)
     * [Keeping `develop` in sync with `main`](#keeping-develop-in-sync-with-main)
* [Opening a Pull Request](#opening-a-pull-request)
* [Testing and Continuous Integration](#testing-and-continuous-integration)
* [Code Review](#code-review)
* [Coding Tips](#coding-tips)

## Code of Conduct
Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## Getting Started
The TASSEL project is written in Java (with some Kotlin) and uses the Gradle build system. To get started, you will need to install the following:

* [Java 21](https://www.oracle.com/java/technologies/downloads/#java17)
* [Git](https://git-scm.com/downloads)

It is recommended to use an IDE to make any code changes. Our group prefers using [IntelliJ IDEA](https://www.jetbrains.com/idea/download/).

Before writing code, install the toolchain and confirm you can build the project — see [Building from Source](developer/building-from-source.md). For anything beyond a trivial fix, open (or find) a GitHub [issue](https://github.com/maize-genetics/tassel/issues) describing the bug or enhancement first, so the work can be discussed and coordinated.

Contributing **documentation only**? You do not need Java, Gradle, or the test data at all — just Python and MkDocs to preview your changes:

```bash
pip install mkdocs-material
mkdocs serve
```

## How to Contribute
For any code changes, you will need to fork the TASSEL repository and create a pull request. For more information on how to do this, please see [this guide](https://docs.github.com/en/get-started/quickstart/fork-a-repo).

### Reporting Bugs
If you find a bug, please first check the [TASSEL issue tracker](https://github.com/maize-genetics/tassel/issues) to see if the bug has already been reported. If it has not, please create a new issue with the bug report template.

### Suggesting Enhancements and New Features
If you have an idea for an enhancement, please create a new issue with the enhancement template in the [TASSEL issue tracker](https://github.com/maize-genetics/tassel/issues).

### Submitting Code Changes
To submit a code change, you first will need to fork the TASSEL repository, make your changes on a branch of your fork, and then submit a Pull Request to the TASSEL repository. The sections below describe the Git workflow, pull request expectations, and the checks your change must pass.

## The Git Workflow
TASSEL uses a branch-and-pull-request model with two long-lived branches:

* **`main`** is the released branch. Every merge of code to `main` produces a new build and release, so `main` should always reflect what users have.
* **`develop`** is the integration branch where normal work accumulates. A [nightly build](developer/releasing.md) publishes an unstable prerelease from it, and it is promoted to `main` when a release is cut.

Which branch you start from depends on the kind of change you are making. There are three tracks, each with its own pull request template.

### Which track am I on?

| Your change | Branch from | PR into | Template | Runs the test suite? |
| --- | --- | --- | --- | --- |
| Feature, enhancement, refactor, non-urgent bug fix, CI work | `develop` | `develop` | feature (default) | Yes |
| Critical fix for a bug already affecting a released version | `main` | `main` | hotfix | Yes |
| Documentation only (`docs/**`, any `*.md`, `mkdocs.yml`) | `main` | `main` | docs | No — site build only |

When in doubt, use the feature track. It is never *wrong*, only slower.

Common to all three tracks:

1. **Fork** the repository (external contributors) or create a branch (team members).
2. **Branch** off the base branch for your track. Branches are cheap — use one per logical piece of work.
3. **Commit** focused, well-described changes.

    ```bash
    git add path/to/changed/files
    git commit -m "Short, imperative summary of the change"
    ```

4. **Push** your branch and open a pull request against the base branch for your track.

Pull the base branch into your branch periodically to reduce merge conflicts:

```bash
git checkout my-branch
git fetch origin
git merge origin/develop   # or origin/main, depending on your track
```

### Feature track (the default)
Nearly all work belongs here, including changes that mix documentation with code.

```bash
git checkout develop
git pull
git checkout -b my-feature
```

Open the PR against `develop`. The default PR template applies. Your change ships to users at the next promotion of `develop` to `main`.

### Hotfix track
Only for an urgent fix to a bug that is already affecting a released version. Branch from `main`, keep the change as small as possible, bump the patch version, and open the PR against `main` with the hotfix template:

```bash
git checkout main
git pull
git checkout -b hotfix/short-description
```

Because this bypasses `develop`, the fix has to reach `develop` too, or the next promotion will reintroduce the bug. You do not need to cherry-pick it: merging to `main` opens an automated `main` → `develop` [back-merge PR](#keeping-develop-in-sync-with-main). Confirming that PR merges is part of the hotfix checklist.

When resolving conflicts in that PR, keep **`develop`'s** `version` in `build.gradle.kts` — a hotfix bumps the patch version of the released line and must not overwrite an in-progress version on `develop`.

### Documentation track
Use this when the change touches **nothing but** documentation: files under `docs/`, any `*.md` file (including `README.md`), and `mkdocs.yml`. Documentation is neither a feature that should wait for the next release nor an emergency, so it gets a shorter path: straight to `main`, with no test suite and no release.

Branch from `main` and prefix the branch name with `docs/`:

```bash
git checkout main
git pull
git checkout -b docs/fix-mlm-example
```

Open the PR against `main` using the docs template (append `?template=docs.md` to the compare URL). The `docs/` branch prefix — or a `documentation` label on the PR — is what marks the PR as being on this track.

What is different about this track:

* **No test suite.** The Java build and tests are skipped; the only check is a fast MkDocs site build that catches a broken `nav` entry or an unbuildable page.
* **No version bump.** Do not change `version` in `build.gradle.kts`.
* **No changelog block.** Documentation merges never reach the release-notes automation, so the `CHANGELOG` markers are not needed.
* **No release.** Merging documentation to `main` does not build the application, cut a GitHub release, or publish to Maven Central. It only redeploys the documentation site, so your change is live within a few minutes.
* **`develop` is synced for you.** Because these commits land on `main` first, an automated [back-merge PR](#keeping-develop-in-sync-with-main) brings them back down.

If a "documentation" change turns out to also need a code edit, move it to the feature track — the **Docs track guard** check fails any `docs/`-prefixed PR that touches files outside the documentation paths.

### Keeping `develop` in sync with `main`
Three things land on `main` without going through `develop`: documentation merges, hotfixes, and the release automation's own commits (it writes `docs/changelog.md` and the download links after each release). Left alone, each one becomes a conflict to untangle at the next promotion.

So a workflow opens a single `main` → `develop` pull request whenever `main` gains content that `develop` does not have, and keeps that one PR up to date as more commits land. A person still merges it — a back-merge can be conflict-free and still be semantically wrong, which is exactly what review is for.

Two things to know when you merge one:

* **Prefer `main` when resolving conflicts**, since everything in the PR is already released — except `build.gradle.kts`, where you keep `develop`'s `version`.
* **No checks run on the sync PR itself**, because it is opened by the CI token. If it carries code, the test suite runs against `develop` once you merge, and the [nightly build](developer/releasing.md#nightly-development-builds-nightlyyml) is the backstop.

## Opening a Pull Request
When you open a PR:

- Confirm the base branch matches your track (`develop` for features, `main` for hotfixes and documentation).
- Fill out the PR template with a clear description of *what* changed and *why*.
- Reference any related issue (e.g. "Closes #123").
- Keep PRs focused. Smaller, single-purpose PRs are reviewed faster.
- Add reviewers from the TASSEL team. If you are unsure who should review, add **@zrm22** and additional reviewers will be assigned.

After you have submitted your Pull Request, verify that all of the automated checks have passed. If any of the checks have failed, review the error message and make any necessary changes. If you are unsure how to fix the error, reach out to the TASSEL team for assistance.

### Changelog notes
The release automation extracts changelog content from the merged PR's description (the text between the `<!-- BEGIN CHANGELOG -->` and `<!-- END CHANGELOG -->` markers in the template). Fill this in so your change is reflected in the published [Version History](changelog.md). Documentation-only PRs are the exception — they do not trigger a release, so they have no changelog block.

## Testing and Continuous Integration
TASSEL uses [JUnit](https://junit.org/) tests run through Gradle. Please add or update tests for your change and make sure the required checks pass before opening a Pull Request.

Which checks run depends on what you changed. CI inspects the changed paths rather than the branch, so a PR that touches no compiled sources skips the Java jobs, and a PR that touches no documentation skips the site build:

| Changed paths | Checks that run |
| --- | --- |
| `src/**` | Full Java CI (build, tests, coverage) |
| `docs/**`, `*.md`, `mkdocs.yml` only | MkDocs site build (about a minute) |
| Both | Both |

A skipped job reports success, so it never blocks a merge.

The full Java CI also runs on pushes to `develop` that touch `src/**` or the build files. That is what verifies a hotfix once it has been back-merged, since no checks run on the sync PR itself.

If you are on the documentation track, you can preview exactly what CI builds with:

```bash
mkdocs serve
```

Fetch the shared test-data archive once after a clean checkout (it is downloaded into the git-ignored `dataFiles/` directory):

```bash
./gradlew fetchTestData
```

Opening or updating a PR that touches `src/**` triggers the CI workflow (JDK 21 with OpenBLAS installed). It has two jobs, matching the two local test entry points:

* **Statistics gate (required):** verifies TASSEL's numeric results against R-validated expected values. This **must pass** for a PR to be mergeable. Run it locally before pushing, especially if you touched analysis or statistics code:

    ```bash
    ./gradlew statisticsTest
    ```

* **Full suite & coverage (non-blocking):** runs everything, including pipeline and I/O tests that are still being stabilized, plus coverage reporting. Failures here are informative but do not block merging:

    ```bash
    ./gradlew test
    ```

The statistical tests exercise native BLAS routines, so you will need OpenBLAS installed (see [Building from Source](developer/building-from-source.md#installing-openblas)).

For full details on the test layout, coverage reports, what CI runs, and how to add an enforced test, see the [Testing guide](developer/testing.md).

## Code Review
A member of the TASSEL team will review your Pull Request and may request changes. Push follow-up commits to the same branch to update the PR. Once approved, the change is merged into its base branch.

Merges of code to `main` trigger the build-and-release automation, which produces a new build and standalone distribution — so a merged hotfix, or a promotion of `develop`, generally results in a new released build. Documentation-only merges are excluded from that automation and instead redeploy the documentation site. See [Releasing](developer/releasing.md).

## Coding Tips
- Match the style and structure of the surrounding code.
- Implement new user-facing functionality as a plugin — see [Developing Plugins](developer/plugin-development.md).
- Add or update tests for your change; for statistical code, wire enforced tests into the statistics gate (see [Testing](developer/testing.md)).
- Prefer throwing exceptions over handling errors inline in plugins; let `AbstractPlugin` present them.
