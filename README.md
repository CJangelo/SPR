
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SPR: Tools to Learn Statistical Programming in R

<!-- badges: start -->
<!-- badges: end -->

SPR, an R package to simulate clinical data as part of training in R
statistical programming.

## Installation

You can install the released version of SPR from github:

``` r
install.packages('devtools')
library(devtools)
devtools::install_github("CJangelo/SPR")
```

Note that this R package is not intended for CRAN. It is for training
purposes.

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

### Fundamentals

This covers the basics of doing analysis as a member of our team.

### Vignettes: Different Data Distributions, Cross Sectional Data

-   Continuous data (Normal error distribution)
-   Binomial data
-   Ordinal data (probit)
-   Ordinal data (logit)
-   Multinomial data
-   Poisson and Negative Binomial data
-   Beta Distribution
-   Zero-inflated Poisson and Negative Binomial data

### Longitudinal Data

Longitudinal data leads to two complications:

-   Marginal versus Conditional model specification
-   Missing Data

#### Marginal versus Conditional Models

As a resource, please review the terrific explanations provided by
Professor of Biostatistics Dimitris Rizopoulos:
<http://www.drizopoulos.com/courses/EMC/CE08.pdf>

One important note for the Vignettes:

-   Continuous data: marginal and conditional model estimates **are**
    equivalent
-   Binomial/Ordinal data: marginal and conditional model estimates
    **are NOT** equivalent

#### Missing Data

Tons of resources on this - can???t go wrong with the NRC panel???s report:
<https://www.ncbi.nlm.nih.gov/books/NBK209904/>

### Vignettes: Longitudinal Data and Missing Data Handling

#### Continuous Data

-   Fit a MMRM in R
-   Simulation study to evaluate MMRM performance, includes FWE
    evaluation, FDR correction
-   How to Handle Missing Data?
-   Continuation of missing data handling: What is the relationship
    between MCAR/MAR/MNAR?
-   Continuation of Missing Data Handling: What is the relationship
    between MCAR/MAR/Conditional MCAR?
-   MMRM and Mixed Effects Model, Simulation Study - what is the
    relationship between marginal and conditional models with continuous
    data?

#### Categorical Data (Binomial, Ordinal)

-   Longitudinal Binomial Data with the Generalized Linear Mixed Model:
    Simulation Study
-   Longitudinal Ordinal Data, Marginal Specification
-   Longitudinal Ordinal Data, Conditional Specification

#### All Data Types

-   Longitudinal Data, All Data Types, Conditional Specification

## Statistics Resources & References

-   We all need reminders for ourselves, and references to point clients
    to
-   Common statistical tests are linear models (or: how to do basic
    stats)
    -   <https://lindeloev.github.io/tests-as-linear/>
-   Missing Data in clinical trials (and analysis in RCT more broadly)
    -   <https://www.ncbi.nlm.nih.gov/books/NBK209904/>
-   Code Guidelines
    -   <https://style.tidyverse.org/>
-   Jargon
    -   NRC stat report had a section on jargon:
        <https://www.ncbi.nlm.nih.gov/books/NBK209903/>
    -   Searchable resource from FDA Working Group:
        <https://www.ncbi.nlm.nih.gov/books/NBK338448/#IX-B>
-   Bug report
    -   How do we report issues/problems to developers/stackoverflow
        community/colleagues?
        -   <https://github.com/rstudio/rstudio/wiki/Writing-Good-Bug-Reports>
        -   Summary: send fully reproducible code and output from
            ???sessionInfo()???
-   Statistical Misconceptions
    -   <https://discourse.datamethods.org/t/reference-collection-to-push-back-against-common-statistical-myths/1787>

## License

This R package is free and open source software (License: GPL (&gt;=
3)).
