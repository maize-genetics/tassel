<!--
  HOTFIX - CRITICAL BUG IN A RELEASED VERSION
  Base branch: main  <-- set the PR base to `main`.

  Use this template ONLY for urgent fixes to a bug that is already affecting
  a released version. Branch from `main` (e.g. `hotfix/short-description`).

  For features, enhancements, or non-urgent bug fixes, use the feature
  template instead (base branch `develop`):
  https://github.com/maize-genetics/tassel/compare/develop...HEAD?template=feature.md

  For a documentation-only change, use the docs template. It also targets
  `main`, but skips the test suite and needs no version bump:
  https://github.com/maize-genetics/tassel/compare/main...HEAD?template=docs.md

  IMPORTANT: this fix must reach `develop` as well, or the next promotion will
  reintroduce the bug. An automated `main` -> `develop` sync PR is opened for
  you after merge — you own confirming that it merges (see the checklist).
-->

## Summary of the critical issue

_What is broken in the released version, and what is the user-facing impact?
Link the issue if one exists._

- Affected release version:
- Impact / severity:

## The fix

_Describe the minimal change that resolves the issue. Keep hotfixes as small
and targeted as possible._

<!-- BEGIN CHANGELOG -->

<!-- END CHANGELOG -->


## Type of change

- [x] `BUGFIX` (non-breaking change which fixes an issue)
- [ ] `CHANGE` (fix that would cause existing functionality to not work as expected)

## Hotfix checklist:

- [ ] The base branch of this PR is `main`
- [ ] I branched from `main` (not `develop`)
- [ ] I bumped the **patch** portion of `version` in `build.gradle.kts` (e.g. `5.2.97` -> `5.2.98`)
- [ ] The `CHANGELOG` tags above contain the user-facing change note
- [ ] The change is minimal and scoped to the critical fix
- [ ] **Backport confirmed:** I will check that the automated `main` → `develop` sync PR merges after this lands, and resolve any conflicts (keeping `develop`'s `version` in `build.gradle.kts`)
