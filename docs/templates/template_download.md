---
title: Download TASSEL 5
---

# Download TASSEL 5

Get the latest release of TASSEL 5. If you want to pick another platform,
distribution, or past version, [see the "Pick a platform &amp; version" section](#choose).

<div class="dl-heroes">
  <div class="dl-hero dl-hero--desktop" id="dl-recommended" data-state="loading">
    <div class="dl-hero__body">
      <p class="dl-hero__eyebrow">Desktop application</p>
      <h2 class="dl-hero__title" id="dl-rec-title">Detecting your platform&hellip;</h2>
      <p class="dl-hero__meta" id="dl-rec-meta"></p>
    </div>
    <a class="md-button md-button--primary dl-hero__btn" id="dl-rec-btn" href="#choose">Choose below</a>
  </div>

  <div class="dl-hero dl-hero--cli" id="dl-standalone" data-state="loading">
    <div class="dl-hero__body">
      <p class="dl-hero__eyebrow">Standalone CLI</p>
      <h2 class="dl-hero__title" id="dl-cli-title">Loading latest version&hellip;</h2>
      <p class="dl-hero__meta" id="dl-cli-meta"></p>
    </div>
    <a class="md-button md-button--primary dl-hero__btn" id="dl-cli-btn" href="#choose">Choose below</a>
  </div>
</div>

## Install as a library (Maven / Gradle)

Add the TASSEL library to your JVM project from Maven Central:

=== "Gradle (Kotlin DSL)"

    ```kotlin
    // build.gradle.kts
    repositories {
        mavenCentral()
    }

    dependencies {
        implementation("net.maizegenetics:tassel:${VERSION}")
    }
    ```

=== "Gradle (Groovy DSL)"

    ```groovy
    // build.gradle
    repositories {
        mavenCentral()
    }

    dependencies {
        implementation 'net.maizegenetics:tassel:${VERSION}'
    }
    ```

=== "Maven"

    ```xml
    <!-- pom.xml -->
    <dependency>
      <groupId>net.maizegenetics</groupId>
      <artifactId>tassel</artifactId>
      <version>${VERSION}</version>
    </dependency>
    ```

For a full walkthrough, see the
[Maven / JVM Library guide](../getting_started/install-maven.md).

## Pick a platform &amp; version { #choose }

<div class="dl-picker">
  <label class="dl-field">
    <span>Distribution</span>
    <select id="dl-dist">
      <option value="installer">GUI Installer (native app)</option>
      <option value="standalone">Standalone (command line)</option>
    </select>
  </label>
  <label class="dl-field" id="dl-os-field">
    <span>Operating system</span>
    <select id="dl-os"></select>
  </label>
  <label class="dl-field">
    <span>Version</span>
    <select id="dl-version"></select>
  </label>
</div>

<div class="dl-result" id="dl-result" aria-live="polite"></div>

<div class="dl-fallback" id="dl-fallback">
  <h3>Direct downloads (latest release)</h3>
  <p>If the live selector above is unavailable, use these direct links:</p>
  <div class="feature-grid">
    <div class="feature-cards">
      <a href="https://github.com/maize-genetics/tassel/releases/download/v${VERSION}/${REL_MAC_SILICON}" class="feature-card">
        <i class="fa fa-apple fa-2x"></i>
        <strong>macOS (Apple Silicon)</strong>
      </a>
      <a href="https://github.com/maize-genetics/tassel/releases/download/v${VERSION}/${REL_MAC_INTEL}" class="feature-card">
        <i class="fa fa-apple fa-2x"></i>
        <strong>macOS (Intel)</strong>
      </a>
      <a href="https://github.com/maize-genetics/tassel/releases/download/v${VERSION}/${REL_LINUX}" class="feature-card">
        <i class="fa fa-linux fa-2x"></i>
        <strong>Linux/Unix</strong>
      </a>
      <a href="https://github.com/maize-genetics/tassel/releases/download/v${VERSION}/${REL_WINDOWS}" class="feature-card">
        <i class="fa fa-windows fa-2x"></i>
        <strong>Windows (64-bit)</strong>
      </a>
      <a href="https://github.com/maize-genetics/tassel/releases/download/v${VERSION}/tassel-5-standalone-v${VERSION}.tar.gz" class="feature-card">
        <i class="fa fa-terminal fa-2x"></i>
        <strong>Standalone</strong>
      </a>
    </div>
  </div>
</div>

!!! tip "Which distribution should I choose?"
    - **GUI Installer (native app):** the easiest way to run the TASSEL desktop
      application. Native installers are built with
      [jDeploy](https://www.jdeploy.com/) and bundle everything you need. See the
      [GUI installation guide](../getting_started/install-gui.md).
    - **Standalone (command line):** a self-contained distribution for running
      TASSEL from the terminal or in pipelines. See the
      [standalone guide](../getting_started/install-standalone.md).
    - **Library (Maven / Gradle):** add TASSEL to your own JVM project.

!!! note
    Native GUI installers, the standalone distribution, and the Maven library are
    all published per version, so older versions remain available in the selector
    and on the
    [GitHub releases page](https://github.com/maize-genetics/tassel/releases).

!!! info "Want to test an unreleased fix?"
    Nightly builds of the `develop` branch are published as unstable
    prereleases. See [Nightly builds](nightly.md) for the latest one and the
    archive of past builds.

!!! warning "Upgrading from an install named TASSEL 5 main?"
    Installers from 5.2.97 and earlier registered the desktop application as
    **TASSEL 5 main**; it is now plain **TASSEL 5**. Uninstall the old copy once
    before installing — see
    [Updating](../getting_started/install-gui.md#updating).
