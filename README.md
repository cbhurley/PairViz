
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PairViz

<!-- badges: start -->
<!-- badges: end -->

The goal of PairViz is to improving graphics by ameliorating order
effects, using Eulerian tours and Hamiltonian decompositions of graphs.

## Installation

You can install the released version of PairViz from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("PairViz")
```

You will also need to install the graph package from Bioconductor.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("graph")
```

To get started, check out the PairVizIntroduction vignette.

## References

C.B. Hurley and R.W. Oldford, Pairwise display of high dimensional
information via Eulerian tours and Hamiltonian decompositions. Journal
of Computational and Graphical Statistics. 19(10), pp. 861–886, 2010.

C.B. Hurley and R.W. Oldford, Eulerian tour algorithms for data
visualization and the PairViz package. Computational Statistics, 26(4),
pp 613–633, 2011.
