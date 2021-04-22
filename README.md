
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SPR: Tools to Learn Statistical Programming in R

<!-- badges: start -->
<!-- badges: end -->

SPR, an R package to simulate clinical data as part of training in R
statistical programming.

## Description

This R package consists of functions that simulate clinical data. These
functions are intended to be used along with Vignettes to illustrate how
to fit statistical models to clinical data in R. Simulation studies can
be conducted using these functions that allow statistical programmers to
evaluate the performance of models in clinical applications.

## Overview of Vignettes

Vignettes are located under `Articles` at the top of this page.

The Vignettes are the focus of the Training. The R package is intended
to be used to simulate data to illustrate the statistical programming in
R with clinical data.

### Missing Data

-   Vignette 1. Fit a MMRM in R
-   Vignette 2. Simulation study to evaluate MMRM performance, includes
    FWE evaluation, FDR correction
-   Vignette 3. How do I handle missing data?
-   Vignette 4. Continuation of missing data handling: What is the
    relationship between MCAR/MAR/MNAR?
-   Vignette 5. Continuation of Missing Data Handling: What is the
    relationship between MCAR/MAR/Conditional MCAR?

### Different Data Distributions: Cross Sectional Data

-   Vignette 6. Continuous data (Normal error distribution)
-   Vignette 7. Binomial Data
-   Vignette 8. Ordinal Data (probit and logit)
-   Vignette 9. Multinomial Data
-   Vignette 10. Poisson and Negative Binomial models

### Under Construction: Missing Data with Binomial/Ordinal Data

-   Vignette 11. Longitudinal Binomial Data with MCAR and Conditional
    MCAR Drop out
-   Vignette 12. Longitudinal Ordinal Data

## Installation

You can install the released version of SPR from github:

``` r
install.packages('devtools')
library(devtools)
devtools::install_github("CJangelo/SPR")
```

Note that this R package is not intended for CRAN. It is for training
purposes.

## License

This R package is free and open source software (License: GPL (&gt;=
3)).
