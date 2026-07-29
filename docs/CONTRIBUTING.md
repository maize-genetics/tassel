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

## How to Contribute
For any code changes, you will need to fork the TASSEL repository and create a pull request. For more information on how to do this, please see [this guide](https://docs.github.com/en/get-started/quickstart/fork-a-repo).

### Reporting Bugs
If you find a bug, please first check the [TASSEL issue tracker](https://github.com/maize-genetics/tassel/issues) to see if the bug has already been reported. If it has not, please create a new issue with the bug report template.

### Suggesting Enhancements and New Features
If you have an idea for an enhancement, please create a new issue with the enhancement template in the [TASSEL issue tracker](https://github.com/maize-genetics/tassel/issues).

### Submitting Code Changes
To submit a code change, you first will need to fork the TASSEL repository, make your changes on a branch of your fork, and then submit a Pull Request to the TASSEL repository. The sections below describe the Git workflow, pull request expectations, and the checks your change must pass.

## The Git Workflow
TASSEL uses a `develop` integration branch: everyday work is merged into
`develop` first, and `main` always reflects the latest released version. Normal
contributions branch off `develop` and are merged back into `develop`. In short:

1. **Fork** the repository (external contributors) or create a branch (team members).
2. **Branch** off `develop` for your change, using a `feature/*` name. Branches
   are cheap — use one per logical piece of work.

    ```bash
    git switch develop
    git pull
    git switch -c feature/my-change
    ```

3. **Commit** focused, well-described changes.

    ```bash
    git add path/to/changed/files
    git commit -m "Short, imperative summary of the change"
    ```

4. **Push** your branch.

    ```bash
    git push -u origin feature/my-change
    ```

5. **Open a pull request** against `develop`.

### Keeping your branch current
Pull the latest `develop` into your branch periodically to reduce merge conflicts:

```bash
git switch develop
git pull
git switch feature/my-change
git merge develop
```

## Opening a Pull Request
When you open a PR:

- Fill out the PR template with a clear description of *what* changed and *why*.
  Normal work uses the default (feature) template with `develop` as the base
  branch. Critical fixes to an already-released version instead use the **hotfix**
  template with `main` as the base branch (see [Releasing](developer/releasing.md#hotfixes)).
- Reference any related issue (e.g. "Closes #123").
- Keep PRs focused. Smaller, single-purpose PRs are reviewed faster.
- Add reviewers from the TASSEL team. If you are unsure who should review, add **@zrm22** and additional reviewers will be assigned.

After you have submitted your Pull Request, verify that all of the automated checks have passed. If any of the checks have failed, review the error message and make any necessary changes. If you are unsure how to fix the error, reach out to the TASSEL team for assistance.

### Changelog notes
The release automation extracts changelog content from the merged PR's description (the text between the `<!-- BEGIN CHANGELOG -->` and `<!-- END CHANGELOG -->` markers in the template). Fill this in so your change is reflected in the published [Version History](changelog.md).

## Testing and Continuous Integration
TASSEL uses [JUnit](https://junit.org/) tests run through Gradle. Please add or update tests for your change and make sure the required checks pass before opening a Pull Request.

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
A member of the TASSEL team will review your Pull Request and may request changes. Push follow-up commits to the same branch to update the PR. Once approved, the change is merged into `develop`.

Merging to `develop` does **not** publish a release. Releases happen when `develop` is promoted to `main` through a separate promotion PR, and merges to `main` trigger the build-and-release automation. See [Releasing](developer/releasing.md).

## Coding Tips
- Match the style and structure of the surrounding code.
- Implement new user-facing functionality as a plugin — see [Developing Plugins](developer/plugin-development.md).
- Add or update tests for your change; for statistical code, wire enforced tests into the statistics gate (see [Testing](developer/testing.md)).
- Prefer throwing exceptions over handling errors inline in plugins; let `AbstractPlugin` present them.
