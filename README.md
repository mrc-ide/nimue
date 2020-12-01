
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nimue

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R build
status](https://github.com/mrc-ide/nimue/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/nimue/actions)
[![codecov](https://codecov.io/gh/mrc-ide/nimue/branch/master/graph/badge.svg)](https://codecov.io/gh/mrc-ide/nimue)
<!-- badges: end -->

## IMPORTANT NOTES

:warning: This code is released with no support. Please submit any
questions or bugs as [issues](https://github.com/mrc-ide/nimue/issues)
and we will try to address them as quickly as possible.

:warning: This model is in active development and so parameter name and
behaviours, and output file formats will change without notice.

:warning: The model is stochastic. Multiple runs with different seeds
should be undertaken to see average behaviour.

:warning: As with any mathematical model, it is easy to misconfigure
inputs and therefore get meaningless outputs. Please contact the authors
if you intend to publish results using `nimue`.

Nimue is built on the shoulders of
[squire](https://mrc-ide.github.io/squire/) and
[sircovid](https://mrc-ide.github.io/sircovid/). The [squire
website](https://mrc-ide.github.io/squire/) is the best initial
reference to refer to for base model structure and parameterisation
details.

Nimue is named after the [Lady of the
Lake](https://en.wikipedia.org/wiki/Lady_of_the_Lake)

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mrc-ide/nimue")
```
