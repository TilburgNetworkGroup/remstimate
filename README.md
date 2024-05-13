<img align="right" width="185" src='man/figures/remstimate-logo.svg'>

# remstimate

<!-- badges: start -->
[![github-repo-status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-package-version](https://img.shields.io/github/r-package/v/TilburgNetworkGroup/remstimate)](https://www.github.com/TilburgNetworkGroup/remstimate)
[![R-CMD-check](https://github.com/TilburgNetworkGroup/remstimate/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/TilburgNetworkGroup/remstimate/actions/workflows/check-standard.yaml)
[![codecov](https://codecov.io/gh/TilburgNetworkGroup/remstimate/graph/badge.svg?token=8NZ4T6E4N9)](https://codecov.io/gh/TilburgNetworkGroup/remstimate)
[![grand-total-downloads](http://cranlogs.r-pkg.org/badges/grand-total/remstimate)](https://cran.r-project.org/package=remstimate)
<!-- badges: end -->

### Optimization Tools for Relational Event History data
The `remstimate` package provides a set of functions that perform necessary calculations when modeling a Relational Event History. It can perform tie-oriented as well as actor-oriented modeling. The main function is `remstimate::remstimate()` which provides four different estimation methods: 

- `"MLE"`, maximizing the model likelihood
- `"GDADAMAX"`, optimization based on the gradient 
- `"BSIR"`, Bayesian Sampling Importance Resampling
- `"HMC"`, Hamiltonian Monte Carlo


### Installation

Install the package in R from CRAN:

```
install.packages("remstimate")
```


### Author and maintainer
Giuseppe Arena, Tilburg University (Tilburg, The Netherlands). (g.arena@tilburguniversity.edu)


### Funding
The funder of this work is the ERC and the ERC project number is 758791.

