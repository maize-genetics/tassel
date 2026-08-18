<!--
  This template was modified from the PhenoApps/fieldbook pull_request_template.md template.

  This is the DEFAULT template, used for features, enhancements, and
  non-urgent bug fixes. The base branch for this work should be `develop`.

  Working on a CRITICAL fix for an already-released version? Use the hotfix
  template (base branch `main`) by opening the PR with `?template=hotfix.md`
  appended to the compare URL, e.g.:
  https://github.com/maize-genetics/tassel/compare/main...HEAD?template=hotfix.md

  Changing ONLY documentation (docs/**, *.md, mkdocs.yml)? Use the docs
  template (base branch `main`) — it skips the test suite and publishes the
  site without cutting a release:
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

- [ ] The base branch of this PR is `develop` (use the hotfix template if this targets `main`)
- [ ] I have updated the `version` variable in `build.gradle.kts`
- [ ] I have performed a self-review of my code
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] I have updated relevant documentation
