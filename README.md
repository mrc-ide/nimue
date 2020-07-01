
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

The goal of nimue is to …

Nimue is named after the [Lady of the
Lake](https://en.wikipedia.org/wiki/Lady_of_the_Lake)

## Installation

You can install the released version of nimue from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("nimue")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mrc-ide/nimue")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(nimue)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
