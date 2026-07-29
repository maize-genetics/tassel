---
hide:
  - navigation
---

# ![rTASSEL](../images/rtassel_logo.svg){ style="height:1.4em; vertical-align:middle" } Use TASSEL from R (rTASSEL)

[**rTASSEL**](https://rtassel.maizegenetics.net/) is an R-based front-end for
accessing key TASSEL 5 methods and tools. It lets you run powerful TASSEL
analyses within a unified R workflow, without switching between environments.

## Install

rTASSEL is installed from GitHub with [`pak`](https://pak.r-lib.org/).

=== "Stable release"

    ```r
    # install.packages("pak")
    pak::pak("maize-genetics/rTASSEL@v0.12.0")
    ```

=== "Developmental version"

    ```r
    # install.packages("pak")
    pak::pak("maize-genetics/rTASSEL")
    ```

!!! note "Requirements"
    rTASSEL uses [`rJava`](https://www.rforge.net/rJava/), so you need a working
    version of **Java (>= 8)**. Mac/Linux users may need to run
    `R CMD javareconf` if they run into issues installing `rJava` via `pak`.

## Verify your setup

For an overview of the available functions once installed, run:

```r
help(package = "rTASSEL")
```

!!! quote "Citation"
    Monier et al., (2022). rTASSEL: An R interface to TASSEL for analyzing
    genomic diversity. *Journal of Open Source Software*, 7(76), 4530,
    <https://doi.org/10.21105/joss.04530>

    For BibTeX users:

    ```bibtex
    @article{Monier2022,
      author  = {Monier, Brandon and Casstevens, Terry M. and Bradbury, Peter J. and Buckler, Edward S.},
      title   = {{rTASSEL}: An {R} interface to {TASSEL} for analyzing genomic diversity},
      journal = {Journal of Open Source Software},
      volume  = {7},
      number  = {76},
      pages   = {4530},
      year    = {2022},
      doi     = {10.21105/joss.04530},
      url     = {https://doi.org/10.21105/joss.04530}
    }
    ```

## Learn more

- [rTASSEL website](https://rtassel.maizegenetics.net/): documentation and articles.
- [Getting started article](https://rtassel.maizegenetics.net/articles/rTASSEL.html): a guided walkthrough of common pipelines.
- [Interactive Binder demo](https://mybinder.org/v2/gh/maize-genetics/rTASSEL_sandbox/main?labpath=getting_started.ipynb): try rTASSEL without installing it locally.
- [rTASSEL on GitHub](https://github.com/maize-genetics/rTASSEL): source code and issue tracker.
