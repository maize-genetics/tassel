---
hide:
  - navigation
---

# Install the Desktop Application

The easiest way to run the TASSEL desktop application is with a **native
installer**. These installers are produced by
[jDeploy](https://www.jdeploy.com/) and bundle a compatible Java runtime, so you
do **not** need to install Java separately.

Installers are attached to each versioned release, alongside the standalone
archives, so the installer you download records the exact TASSEL version it
installs.

!!! tip "Just want the download?"
    Visit the [Download page](../download/index.md). It auto-detects your
    operating system and offers the correct installer, plus every other build.

## Choose your platform

Every installer filename includes the TASSEL version, so pick the asset for your
platform from the
[latest release](https://github.com/maize-genetics/tassel/releases/latest):

| Platform | File type |
| -------- | --------- |
| macOS (Apple Silicon) | `TASSEL.5.Installer-mac-arm64-...tgz` |
| macOS (Intel) | `TASSEL.5.Installer-mac-x64-...tgz` |
| Windows (64-bit) | `TASSEL.5.Installer-win-x64-...exe` |
| Linux / Unix | `TASSEL.5.Installer-linux-x64-...tar.gz` |

## macOS

1. Download the `TASSEL.5.Installer-mac-arm64-...tgz` asset if you have Apple
   Silicon (M1/M2/M3 and newer), or `TASSEL.5.Installer-mac-x64-...tgz` for
   Intel. Not sure which you have? Click the Apple menu &rarr; **About This Mac**
   and check the **Chip** / **Processor** line.
2. Double-click the downloaded `.tgz` to expand it, then open the resulting app.
3. If macOS Gatekeeper blocks the app because it is from an unidentified
   developer, right-click (or Control-click) the app and choose **Open**, then
   confirm. You only need to do this the first time.

## Windows

1. Download the `TASSEL.5.Installer-win-x64-...exe` asset.
2. Run the `.exe`. If Windows SmartScreen appears, choose **More info** &rarr;
   **Run anyway**.
3. Follow the prompts. TASSEL is added to the Start menu when the installer
   finishes.

## Linux

1. Download the `TASSEL.5.Installer-linux-x64-...tar.gz` asset.
2. Extract the archive:

    ```bash
    tar -xzf TASSEL.5.Installer-linux-x64-*.tar.gz
    ```

3. Run the extracted launcher to start TASSEL.

## Updating

The installer asks how you want updates handled, defaulting to tracking stable
releases. With that default, TASSEL updates itself when a newer release is
published, so you do not need to download an installer again. You can always
install a specific build from the
[GitHub releases page](https://github.com/maize-genetics/tassel/releases).

!!! warning "Already have an app named TASSEL 5 main?"
    Installers from 5.2.97 and earlier registered the application as
    **TASSEL 5 main**, because they were built from the `main` branch rather than
    a version. From 5.2.98 onward it is registered as plain **TASSEL 5**. Your
    operating system treats these as two different applications, and the old one
    will not receive further updates, so uninstall **TASSEL 5 main** once and
    then run the current installer.

## Next steps

- New to the interface? Start with the [User Manual](../user_manual/index.md).
- Prefer the command line? See the [Standalone guide](install-standalone.md).
