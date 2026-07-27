<!--
  HOTFIX - CRITICAL BUG IN A RELEASED VERSION
  Base branch: main  <-- set the PR base to `main`.

  Use this template ONLY for urgent fixes to a bug that is already affecting
  a released version. Branch from `main` (e.g. `hotfix/short-description`).

  For features, enhancements, or non-urgent bug fixes, use the feature
  template instead (base branch `develop`):
  https://github.com/maize-genetics/tassel/compare/develop...HEAD?template=feature.md

  IMPORTANT: after this hotfix merges to `main`, it MUST be backported to
  `develop` (see the checklist) so the next promotion does not reintroduce
  the bug.
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
- [ ] **Backport planned:** I will open a follow-up PR (or cherry-pick) to bring this fix into `develop` immediately after merge
