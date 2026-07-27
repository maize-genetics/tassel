<!--
  DOCUMENTATION ONLY
  Base branch: main  <-- set the PR base to `main`.

  Use this template when the change touches nothing but documentation:
  `docs/**`, any `*.md` file (including README.md), and `mkdocs.yml`.

  Branch from `main` and name the branch `docs/short-description`. The `docs/`
  prefix is what tells CI to treat this as a documentation change: the Java test
  suite is skipped, only the MkDocs site build runs, and merging publishes the
  site without cutting a software release.

  If this PR also changes source, the build, or workflows, CI will fail the
  "Docs track guard" check. Use the feature template against `develop` instead:
  https://github.com/maize-genetics/tassel/compare/develop...HEAD?template=feature.md
-->

## Description

_What documentation changed, and why? Link any related issue._

## Pages affected

_List the pages or sections touched, so a reviewer knows where to look._

-

## Checklist:

- [ ] The base branch of this PR is `main`, and I branched from `main`
- [ ] My branch is named `docs/...` (or the PR carries the `documentation` label)
- [ ] **Only** documentation files changed — nothing under `src/**`, no `build.gradle.kts`, no workflow changes
- [ ] I did **not** bump `version` in `build.gradle.kts` (documentation merges do not cut a release)
- [ ] I previewed the rendered result locally with `mkdocs serve`
- [ ] Any new page is wired into the `nav` in `mkdocs.yml`
- [ ] Links and images resolve, and code samples are accurate for the current release

<!--
  After merge: the site redeploys automatically, and a `main` -> `develop` sync
  PR is opened (or updated) so `develop` picks these changes up. No changelog
  block is needed here — documentation merges never reach the release-notes
  automation.
-->
