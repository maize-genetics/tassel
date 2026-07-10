# Contributing Code

This page covers the mechanics of contributing code to TASSEL: the Git workflow,
pull requests, and code review. For the higher-level contribution policy, see the
project's [Contributing guide](../CONTRIBUTING.md) and
[Code of Conduct](../CODE_OF_CONDUCT.md).

TASSEL is developed on GitHub at
[maize-genetics/tassel](https://github.com/maize-genetics/tassel).

## Before you start

- Install the toolchain and confirm you can build the project — see
  [Building from Source](building-from-source.md).
- For anything beyond a trivial fix, open (or find) a GitHub
  [issue](https://github.com/maize-genetics/tassel/issues) describing the bug or
  enhancement first, so the work can be discussed and coordinated.

## The Git workflow

TASSEL uses a branch-and-pull-request model. In short:

1. **Fork** the repository (external contributors) or create a branch (team
   members).
2. **Branch** off `main` for your change. Branches are cheap — use one per
   logical piece of work and keep `main` clean.

    ```bash
    git checkout main
    git pull
    git checkout -b my-feature
    ```

3. **Commit** focused, well-described changes.

    ```bash
    git add path/to/changed/files
    git commit -m "Short, imperative summary of the change"
    ```

4. **Push** your branch.

    ```bash
    git push -u origin my-feature
    ```

5. **Open a pull request** against `main`.

### Keeping your branch current

Pull the latest `main` into your branch periodically to reduce merge conflicts:

```bash
git checkout main
git pull
git checkout my-feature
git merge main
```

## Opening a pull request

When you open a PR:

- Fill out the PR template with a clear description of *what* changed and *why*.
- Reference any related issue (e.g. "Closes #123").
- Keep PRs focused. Smaller, single-purpose PRs are reviewed faster.
- Add reviewers from the TASSEL team. If you are unsure who should review, add
  **@zrm22** and additional reviewers will be assigned.

### Changelog notes

The release automation extracts changelog content from the merged PR's
description (the text between the `<!-- BEGIN CHANGELOG -->` and
`<!-- END CHANGELOG -->` markers in the template). Fill this in so your change is
reflected in the published [Version History](../changelog.md).

## Continuous integration

Opening or updating a PR that touches `src/**` triggers the CI workflow. It runs
on JDK 21 with OpenBLAS installed and has two jobs:

- **Statistics gate (required)** — `./gradlew statisticsTest`. This **must pass**
  for the PR to be mergeable.
- **Full suite & coverage (non-blocking)** — the broader `test` task plus
  coverage reporting. Failures are visible but do not block merging while those
  tests are being stabilized.

Run the statistics gate locally before pushing to catch problems early:

```bash
./gradlew statisticsTest
```

See [Testing](testing.md) for details.

## Code review

A member of the TASSEL team will review your PR and may request changes. Push
follow-up commits to the same branch to update the PR. Once approved, the change
is merged into `main`.

Merges to `main` trigger the build-and-release automation, which produces a new
build and standalone distribution — so a merged PR generally results in a new
released build. See [Releasing](releasing.md).

## Coding tips

- Match the style and structure of the surrounding code.
- Implement new user-facing functionality as a plugin — see
  [Developing Plugins](plugin-development.md).
- Add or update tests for your change; for statistical code, wire enforced tests
  into the statistics gate (see [Testing](testing.md)).
- Prefer throwing exceptions over handling errors inline in plugins; let
  `AbstractPlugin` present them.
