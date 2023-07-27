concordance
================

[![License: GPL
v3](https://img.shields.io/badge/License-BSD_3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![](https://img.shields.io/badge/devel%20version-1.0.0-purple.svg)](https://github.com/knrumsey/concordance)

<!-- README.md is generated from README.Rmd. Please edit that file -->

<div class="figure">

<img src="inst/logos/CONCORDANCE.png" alt="This logo was designed by Imagine AI Art Studio" width="50%" />
<p class="caption">
This logo was designed by Imagine AI Art Studio
</p>

</div>

### Description

`concordance` is an R package for performing concordance analyses and
for the [discovery of active subspaces in
high-dimensions](https://arxiv.org/pdf/2307.11241.pdf) (described in
Rumsey, Francom, and Wiel (2023)). The “workhorse” of the package is the
`C_bass()` function, which estimates *Constantine’s* $C$ matrix (see
(**constantine2015?**))for a given computer model and behaves similarly
to the `activegp::C_gp()` function (described in Wycoff, Binois, and
Wild (2021)). The `C_bass()` function, which relies on a Bayesian MARS
emulator (as described in Francom and Sansó (2020) and implemented in
the [BASS package](https://CRAN.R-project.org/package=BASS)) is likely
to be more efficient and accurate when the dimension of the input space
is large and admits a large class of *measures* for the inputs.

To install this package, use

``` r
# install.packages("remotes")
remotes::install_github("knrumsey/concordance@main")
```

This R package is based on the method described
[here](https://arxiv.org/pdf/2307.11241.pdf). A vignette for this
package is also available and can be viewed by typing
`browseVignettes("concordance")` in R (after installation).

## Copyright Notice

© *2023. Triad National Security, LLC. All rights reserved.*

*This program was produced under U.S. Government contract
89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
operated by Triad National Security, LLC for the U.S. Department of
Energy/National Nuclear Security Administration. All rights in the
program are reserved by Triad National Security, LLC, and the U.S.
Department of Energy/National Nuclear Security Administration. The
Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to
reproduce, prepare derivative works, distribute copies to the public,
perform publicly and display publicly, and to permit others to do so.*

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-francom2020bass" class="csl-entry">

Francom, Devin, and Bruno Sansó. 2020. “Bass: An r Package for Fitting
and Performing Sensitivity Analysis of Bayesian Adaptive Spline
Surfaces.” *Journal of Statistical Software* 94 (LA-UR-20-23587).

</div>

<div id="ref-rumsey2023discovering" class="csl-entry">

Rumsey, Kellin N, Devin Francom, and Scott Vander Wiel. 2023.
“Discovering Active Subspaces for High-Dimensional Computer Models.”
*arXiv Preprint arXiv:2307.11241*.

</div>

<div id="ref-wycoff2021sequential" class="csl-entry">

Wycoff, Nathan, Mickael Binois, and Stefan M Wild. 2021. “Sequential
Learning of Active Subspaces.” *Journal of Computational and Graphical
Statistics* 30 (4): 1224–37.

</div>

</div>
