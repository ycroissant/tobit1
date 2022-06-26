## Overview

`tobit1` is an `R` package which fits and tests models for which the
response is truncated. 


## Estimators

A large amount of estimating methods are provided:

- maximum likelihood estimator,
- two-steps estimator,
- symetrically trimmed estimator,
- non-linear least square estimator.

The estimation can be performed either on the censored or the
truncated sample.

Full fonctionalities are provided for the most popular case where the
response is zero-truncated response. Some estimators are also provided
for responses where the left truncation point is not 0, where the
response is right truncated and also for two-limits response.

## Specification tests

Specification tests are of particular interest for tobit models as the
consistency of most estimators relies on the hypothesis of
homoscedasticity and normality. `tobit1` provide an `haustest` for the
Hausman (1978) test and, for conditional moment tests, the `cmtest`
can be used as it provides specific methods for `tobit1` objects.

## Endogeneity

The endogeneity test of Smith and Blundell (1986) is provided as the
`smbltest` function. The estimation of simultaneous-equation Tobit
models can be performed using either a two-steps estimator or maximum
likelihood using the `ivtobit` function.

## Installation

`tobit1` is not currently available on `CRAN`. Use:

`remotes::install_github("ycroissant/tobit1", build_vignettes = TRUE)`
