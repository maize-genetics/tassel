<!--
  FEATURE / ENHANCEMENT / NON-CRITICAL BUGFIX
  Base branch: develop  <-- set the PR base to `develop`, NOT `main`.

  Use this template for normal work: new features, enhancements, refactors,
  CI, and bug fixes that are NOT urgent production issues. It also covers
  changes that mix documentation with code.

  For a critical bug already affecting a released version, use the hotfix
  template instead:
  https://github.com/maize-genetics/tassel/compare/main...HEAD?template=hotfix.md

  For a documentation-only change, use the docs template (base branch `main`),
  which skips the test suite and goes live without a release:
  https://github.com/maize-genetics/tassel/compare/main...HEAD?template=docs.md
-->

## Description

_Provide a summary of your changes including motivation, context,
and a bullet list of **concise**, **human-readable changes**. The
bullet list should be in the following `CHANGELOG` tags. If these
changes fix a bug or resolves a feature request, be sure to link to
that issue._

<!-- BEGIN CHANGELOG -->

<!-- END CHANGELOG -->


## Type of change

_What type of changes does your code introduce? Put an `x` in boxes that apply._

- [ ] `CHANGE` (fix or feature that would cause existing functionality to not work as expected)
- [ ] `FEATURE` (non-breaking change which adds functionality)
- [ ] `BUGFIX` (non-breaking change which fixes an issue)
- [ ] `ENHANCEMENT` (non-breaking change which improves existing functionality)
- [ ] `NONE` (if none of the other choices apply. Example, tooling, build system, CI, docs, etc.)

## Checklist:

- [ ] The base branch of this PR is `develop`
- [ ] I have updated the `version` variable in `build.gradle.kts`
- [ ] I have performed a self-review of my code
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] I have updated relevant documentation
