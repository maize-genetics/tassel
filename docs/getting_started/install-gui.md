---
hide:
  - navigation
---

# Install the Desktop Application

The easiest way to run the TASSEL desktop application is with a **native
installer**. These installers are produced by
[jDeploy](https://www.jdeploy.com/) and bundle a compatible Java runtime, so you
do **not** need to install Java separately.

Installers are published to the rolling `main` release and always reflect the
latest build.

!!! tip "Just want the download?"
    Visit the [Download page](../download/index.md). It auto-detects your
    operating system and offers the correct installer, plus every other build.

## Choose your platform

| Platform | File type |
| -------- | --------- |
| macOS (Apple Silicon) | `TASSEL.5.Installer-mac-arm64-...tgz` |
| macOS (Intel) | `TASSEL.5.Installer-mac-x64-...tgz` |
| Windows (64-bit) | `TASSEL.5.Installer-win-x64-...exe` |
| Linux / Unix | `TASSEL.5.Installer-linux-x64-...tar.gz` |

## macOS

1. Download the installer that matches your chip:
   [Apple Silicon](https://github.com/maize-genetics/tassel/releases/download/main/TASSEL.5.Installer-mac-arm64-@main_26DT.tgz)
   (M1/M2/M3 and newer) or
   [Intel](https://github.com/maize-genetics/tassel/releases/download/main/TASSEL.5.Installer-mac-x64-@main_26DT.tgz).
   Not sure which you have? Click the Apple menu &rarr; **About This Mac** and
   check the **Chip** / **Processor** line.
2. Double-click the downloaded `.tgz` to expand it, then open the resulting app.
3. If macOS Gatekeeper blocks the app because it is from an unidentified
   developer, right-click (or Control-click) the app and choose **Open**, then
   confirm. You only need to do this the first time.

## Windows

1. Download the [Windows installer](https://github.com/maize-genetics/tassel/releases/download/main/TASSEL.5.Installer-win-x64-@main_26DT.exe).
2. Run the `.exe`. If Windows SmartScreen appears, choose **More info** &rarr;
   **Run anyway**.
3. Follow the prompts. TASSEL is added to the Start menu when the installer
   finishes.

## Linux

1. Download the [Linux installer](https://github.com/maize-genetics/tassel/releases/download/main/TASSEL.5.Installer-linux-x64-@main_26DT.tar.gz).
2. Extract the archive:

    ```bash
    tar -xzf TASSEL.5.Installer-linux-x64-*.tar.gz
    ```

3. Run the extracted launcher to start TASSEL.

## Updating

Because native installers track the rolling `main` release, downloading and
reinstalling the latest installer updates you to the newest build. You can always
find previous builds on the
[GitHub releases page](https://github.com/maize-genetics/tassel/releases).

## Next steps

- New to the interface? Start with the [User Manual](../user_manual/index.md).
- Prefer the command line? See the [Standalone guide](install-standalone.md).
