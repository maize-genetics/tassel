<img src="docs/images/tassel_logo.svg" align="right" width="100" alt="TASSEL logo">

## TASSEL

[![TASSEL 5 CI](https://github.com/maize-genetics/tassel/actions/workflows/coverage.yml/badge.svg)](https://github.com/maize-genetics/tassel/actions/workflows/coverage.yml)
[![Maven Central](https://img.shields.io/maven-central/v/net.maizegenetics/tassel?label=Maven%20Central)](https://central.sonatype.com/artifact/net.maizegenetics/tassel)

TASSEL (**T**rait **A**nalysis by a**SS**ociation, **E**volution and **L**inkage) is a
software package for evaluating trait associations, evolutionary patterns, and
linkage disequilibrium in genetic data. It is designed to handle the diversity of
data types and sizes common in modern genomics, including a variety of genotype and
phenotype file formats.

* [**Website**](https://maize-genetics.github.io/tassel/)
* [**Latest Builds**](https://github.com/maize-genetics/tassel/releases/latest)

### Quick start

If you want to build this from the source, TASSEL needs **Java 21**. 
A Gradle wrapper (`./gradlew`) is included, so a local Gradle installation is 
not required. A native BLAS library (**OpenBLAS**) is optional but 
recommended for fast matrix operations.

Build the runnable JAR and launch the graphical interface:

```bash
./gradlew build   # assemble build/libs/sTASSEL.jar
./gradlew run     # launch the GUI
```

For prerequisites (including OpenBLAS setup), command-line pipelines, and the
standalone distribution, see the documentation links below.

### Documentation

Full documentation lives on the [TASSEL docs site](https://maize-genetics.github.io/tassel/):

| Resource                                                                                   | Description                                                                  |
| ------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------- |
| [Getting Started & Installation](https://maize-genetics.github.io/tassel/getting_started/) | Install the desktop app, standalone distribution, JVM library, or rTASSEL.   |
| [User Manual](https://maize-genetics.github.io/tassel/user_manual/)                        | GUI reference for every analysis, filter, and data operation.                |
| [Pipelines](https://maize-genetics.github.io/tassel/pipelines/)                            | Command-line pipeline documentation for scripted and batch workflows.        |
| [Developer Guide](https://maize-genetics.github.io/tassel/developer/)                      | Build from source, architecture, plugin development, testing, and releasing. |
| [Contributing](https://maize-genetics.github.io/tassel/CONTRIBUTING/)                      | Git workflow, pull requests, and code review.                                |
| [Version History](https://maize-genetics.github.io/tassel/changelog/)                      | Release notes and changelog.                                                 |

To preview the site locally, install `mkdocs-material` and run `mkdocs serve`.

### Contributing

Start from the branch that matches your change: `develop` for features and
non-urgent fixes, `main` for hotfixes and documentation-only changes.
Documentation changes take a shortcut — they skip the test suite and publish
without cutting a release. See
[Contributing](https://maize-genetics.github.io/tassel/CONTRIBUTING/#the-git-workflow)
for the full workflow.

### Citation

If you use TASSEL in your research, please cite the overall package:

> Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES. (2007)
> [TASSEL: Software for association mapping of complex traits in diverse samples.](https://academic.oup.com/bioinformatics/article/23/19/2633/185151)
> Bioinformatics 23:2633-2635.

BibTeX entry:

```bibtex
@article{Bradbury2007,
  author  = {Bradbury, Peter J. and Zhang, Zhiwu and Kroon, Dallas E. and Casstevens, Terry M. and Ramdoss, Yogesh and Buckler, Edward S.},
  title   = {{TASSEL}: Software for association mapping of complex traits in diverse samples},
  journal = {Bioinformatics},
  volume  = {23},
  number  = {19},
  pages   = {2633--2635},
  year    = {2007},
  doi     = {10.1093/bioinformatics/btm311},
  url     = {https://academic.oup.com/bioinformatics/article/23/19/2633/185151}
}
```

