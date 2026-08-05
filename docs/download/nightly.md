---
title: Nightly builds
render_macros: true
---

# Nightly builds

Nightly builds are automated snapshots of the
[`develop` branch](https://github.com/maize-genetics/tassel/tree/develop),
published so that fixes and new features can be tested before they are released.

!!! warning "These are not releases"
    A nightly build is unstable, has no version number of its own, and has not
    been through the release process. It is never published to Maven Central.
    For analyses you intend to publish, use the
    [latest stable release](index.md) instead.

## Latest nightly

The `dev-latest` tag always points at the most recent nightly, so these links
never go stale:

<div class="feature-grid">
  <div class="feature-cards">
    <a href="https://github.com/maize-genetics/tassel/releases/download/dev-latest/tassel-5-standalone-nightly.tar.gz" class="feature-card">
      <i class="fa fa-terminal fa-2x"></i>
      <strong>Standalone (.tar.gz)</strong>
    </a>
    <a href="https://github.com/maize-genetics/tassel/releases/download/dev-latest/tassel-5-standalone-nightly.zip" class="feature-card">
      <i class="fa fa-terminal fa-2x"></i>
      <strong>Standalone (.zip)</strong>
    </a>
  </div>
</div>

The version and commit each build came from are recorded in the
[`dev-latest` release notes](https://github.com/maize-genetics/tassel/releases/tag/dev-latest).
Because the filenames are constant, the nightly can also be fetched from a
script:

```bash
curl -LO https://github.com/maize-genetics/tassel/releases/download/dev-latest/tassel-5-standalone-nightly.tar.gz
tar -xzf tassel-5-standalone-nightly.tar.gz
```

The same release also carries a dated copy of each archive
(`tassel-5-standalone-v<version>-dev.<date>.tar.gz`), so you can tell which
build you downloaded. The constant-named archive is a copy of that dated one, so
it unpacks into a dated folder such as `tassel-5-standalone-v{{ version }}-dev.20260801/`
— which is the other way to identify a build fetched from the stable URL above.
Running it is otherwise identical to the released standalone distribution — see
the [standalone guide](../getting_started/install-standalone.md).

## All nightly builds { #all }

Every nightly is also archived under its own dated tag, `dev-YYYYMMDD`. The
newest 14 are kept; older ones are deleted automatically.

<div id="nightly-list" data-state="loading">
  <p id="nightly-status">Loading nightly builds&hellip;</p>
</div>

<div id="nightly-fallback">
  <p>
    If the list above is unavailable, browse the
    <a href="https://github.com/maize-genetics/tassel/releases">GitHub releases
    page</a> and look for prereleases tagged <code>dev-</code>.
  </p>
</div>

## How they are built

The [nightly workflow](https://github.com/maize-genetics/tassel/blob/main/.github/workflows/nightly.yml)
runs at 05:00 UTC and builds only when `develop` has new commits since the last
nightly, so there is no build on quiet days. It runs the full test suite first:
a nightly that fails to appear means `develop` is broken, not that nothing
changed. For the full picture, see
[Releasing](../developer/releasing.md#nightly-dev-builds).
