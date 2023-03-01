concordance
================

[![License: GPL
v3](https://img.shields.io/badge/License-BSD_3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![](https://img.shields.io/badge/devel%20version-1.0.0-purple.svg)](https://github.com/knrumsey/concordance)

<!-- README.md is generated from README.Rmd. Please edit that file -->

`concordance` is an R package for performing concordance analyses and
for the discovery of active subspaces (Constantine (2015)) in
high-dimensions. The “workhorse” of the package is the `C_bass()`
function, which estimates *Constantine’s* \(C\) matrix for a given
computer model and behaves similarly to the `activegp::C_gp()` function
(described in Wycoff, Binois, and Wild (2021)). The `C_bass()` function,
which relies on a Bayesian MARS emulator (as described in Francom and
Sansó (2020) and implemented in the [BASS
package](https://CRAN.R-project.org/package=BASS)) is likely to be more
efficient and accurate when the dimension of the input space is large
and admits a large class of *measures* for the inputs.

To install this package, use

``` r
# install.packages("devtools")
devtools::install_github("knrumsey/concordance")
```

A paper describing this method will be added once published. A vignette
for this package is also available.

## References

<div id="refs" class="references">

<div id="ref-constantine2015active">

Constantine, Paul G. 2015. *Active Subspaces: Emerging Ideas for
Dimension Reduction in Parameter Studies*. SIAM.

</div>

<div id="ref-francom2020bass">

Francom, Devin, and Bruno Sansó. 2020. “Bass: An R Package for Fitting
and Performing Sensitivity Analysis of Bayesian Adaptive Spline
Surfaces.” *Journal of Statistical Software* 94 (LA-UR-20-23587).

</div>

<div id="ref-wycoff2021sequential">

Wycoff, Nathan, Mickael Binois, and Stefan M Wild. 2021. “Sequential
Learning of Active Subspaces.” *Journal of Computational and Graphical
Statistics* 30 (4): 1224–37.

</div>

</div>
