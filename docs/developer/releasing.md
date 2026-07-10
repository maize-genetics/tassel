# Releasing

TASSEL releases are largely **automated through GitHub Actions**. Merging to
`main` builds the application, publishes a GitHub release with installers and a
standalone distribution, and redeploys the documentation site. Publishing to
Maven Central is a separate, manually-triggered workflow.

## Versioning

The project version is declared in `build.gradle.kts`:

```kotlin
version = "5.2.96"
```

Bump this value when cutting a new release. You can confirm the current version
with:

```bash
./gradlew printVersion
```

The build workflow reads this version to name the GitHub release
(`v<version>`), the standalone archives, and the installer bundles.

## What happens on merge to `main`

Two workflows chain together automatically.

### 1. Build & publish release (`jdeploy.yml`)

Triggered on push to `main`, this "jDeploy CI with Gradle" workflow:

1. Sets up JDK 21 and reads the version via `./gradlew printVersion`.
2. Builds the app: `./gradlew clean build -x test -x koverVerify`.
3. Assembles the **standalone distribution** in `dist/tassel-5-standalone/` —
   `sTASSEL.jar`, the `lib/` dependencies, and the launcher scripts
   (`start_tassel.pl`, `run_pipeline.pl`, `run_pipeline.bat`) with their paths
   rewritten for the flattened layout.
4. Archives the standalone as both `.zip` and `.tar.gz`.
5. Extracts changelog content from the latest merged PR description and renders
   release notes from `.github/release_template.md`.
6. Creates a **GitHub Release** tagged `v<version>` and uploads the standalone
   archives.
7. Runs [jDeploy](https://www.jdeploy.com/) to build native **installer bundles**
   (Linux x64, macOS arm64/x64, Windows x64) and attaches them to the release.
8. Updates the download links in the release notes and in
   `docs/overrides/partials/featured_downloads.html`, then commits the updated
   `docs/changelog.md` and download links back to the repository.

### 2. Deploy the documentation site (`deploy_project_site.yml`)

After the build workflow completes on `main`, the "Site Builder and Deployment"
workflow builds this MkDocs site with `mkdocs gh-deploy --force`, publishing it to
GitHub Pages.

## Publishing to Maven Central (`run_publish_maven.yml`)

Publishing the library to Maven Central is handled by a separate workflow and is
**not** triggered automatically by every merge. It relies on the `maven-publish`,
`signing`, and `jreleaser` configuration in `build.gradle.kts`:

- The publication artifact id is `tassel` under group `net.maizegenetics`.
- Artifacts are GPG-signed. The signing key and passphrase are supplied via the
  `JRELEASER_GPG_SECRET_KEY` and `JRELEASER_GPG_PASSPHRASE` environment
  variables/secrets.
- API docs are generated with Dokka (`dokkaJar`) and attached to the publication
  alongside the sources JAR.
- JReleaser deploys the staged repository (`build/staging-deploy`) to Sonatype /
  Maven Central.

## Cutting a release: checklist

1. Ensure `main` is green (the statistics gate passes).
2. Update `version` in `build.gradle.kts`.
3. Make sure the merging PR's description contains the changelog block (between
   the `<!-- BEGIN CHANGELOG -->` / `<!-- END CHANGELOG -->` markers) so release
   notes and `docs/changelog.md` are generated correctly.
4. Merge to `main`. The build/release and site-deploy workflows run
   automatically.
5. Verify the new GitHub Release, its attached standalone archives and
   installers, and the updated [Version History](../changelog.md).
6. If publishing to Maven Central, run the `run_publish_maven.yml` workflow and
   confirm the artifacts appear on Central.

## Local dry runs

You can reproduce the key build steps locally to sanity-check a release before
merging:

```bash
# Build exactly as the release workflow does
./gradlew clean build -x test -x koverVerify

# Confirm the version string
./gradlew printVersion

# Build the docs site locally (requires mkdocs-material)
mkdocs serve
```
