
# remstimate 2.3.9

* minor fix on titles of Schoenfeld's residuals plots;
* correction on learning rate parameter of "GDADAMAX" method.

# remstimate 2.3.8

* DESCRIPTION file, removed quotes from acronyms;
* documentation (return values added where missing);
* tests (removed test on 'ncores' argument);
* vignette (switched to 'html_document', edited seed value and parameters for `"HMC"` method).
.
# remstimate 2.3.7 (initial CRAN release)

* trace plot for posterior draws (only via `"HMC"` method) with highest posterior density intervals;
* histograms of posterior draws (for both `"BSIR"` and `"HMC"` method) with highest posterior density intervals;
* the content of `summary.remstimate()` can be saved by assigning it to a variable;
* Watanabe-Akaike's Information Criterion (WAIC) calculation available. New argument, logical value, `WAIC`, with default `FALSE`. If `WAIC=TRUE`, then an optional argument `nsimWAIC` (number of draws to estimate the WAIC) is set by default to `500` (the user can supply a different value);
* updated parallelization of routines as to the computation of likelihood, gradient and hessian in the actor-oriented modeling framework;
* updated parallelization of routines for computing model residuals and WAIC;
* added new tests, increased coverage.

# remstimate 2.3.6 - 2.3.1

* several minor updates before initial CRAN release `remstimate 2.3.7`.

# remstimate 2.3.0

* Simultaneous events are supported (likelihood functions adapted).

# remstimate 2.2.0:

* major updates of Rcpp functions;
* now using `tinytest` for testing (all tests are converted from `testthat`);
* new examples to all functions in the package;
* remove Rcpp functions (`emp_dist_longest_batch()`, for tuning L parameter in `"HMC"` method);
* vignette created for actor-oriented modeling framework. 

# remstimate 2.1.0:

* new remstimate logo;
* update of functions' description;
* major update of `remstimate()` function:
    - removed input argument `model`;
    - initial controls of input argument;
    - output structure;
    - `"HMC"` for actor-oriented framework;
    - added initial values of `"HMC"` method (experimental);
* major updates of Rcpp functions:
    - new experimental Rcpp function `set_seed`;
    - input argument `edgelist` changed with `dyad`;
    - new input arguments `actor1` and `actor2` (for actor-oriented routines);
    - removed old Rcpp functions: `getDyadIndex()`, `getDyadComposition()`;
    - added dependency of exported C++ routines from `remify`, using now `remify::getDyadIndex()` and `remify::getDyadComposition()`;
* added new tests using `testthat`. Tests are divided in multiple .R scripts depending on what is tested;
* update summary and print of actor-oriented framework in `summary.remstimate`;
* `predict.remstimate` and `plot.remstimate` function are not yet available;
* added dependencies: `remify` and `remstats`;
* added data for tie-oriented and actor-oriented examples;
* new README.md with badges.

# remstimate 2.0.0:

* This version is adapted to the latest changes coming from `remify 2.0.0` and it can estimate a Tie-Oriented model as well as an Actor-Oriented model. Models can be estimated by means of different methods: `"MLE"`, `"GDADAMAX"` (replacing the former `"GD"` and `"GDADAM"`), `"BSIR"` and `"HMC"`. Methods like `"BSIR"` and `"HMC"` are ready-to-use but still under a continuous development in order to improve the user-experience;
* removed experimental "fast method" to compute the likelihood;
* added dependencies `trust` and `parallel`;
* experimental Rcpp functions `posteriorRank()` and `remDerivativesStandard_lambdas()` (added but finally removed in the version 2.0.0);
* added new Rcpp functions: `getDyadIndex()`, `getDyadComposition()`.

# remstimate 1.0.0:

* _11/12/2020_:
    - Methods working with the function `remstimate()` are: `"MLE"`, `"GD"`, `"GDADAM"`, `"BSIR"`, `"HMC"`. However the output lacks of a structure attributes and methods;
* _04/09/2020_ :
    - _messages_ becomes again an Rcpp file. This exstension appears to suit better the intent of the content/aim of error and warning messages;
    - _remstimate.R_ contains the main function `remstimate()` which is aimed to run either a Frequentist or a Bayesian approach by using different optimization/methods. It also includes a switch to the "fast method" to compute the likelihood. The "fast method" is run if the actual improvement (percentage of improvement) is higher than a threshold set by the user (default is `0.5`);
    - _reh.cpp_ : contains `getRisksetMatrix()`, `getRisksetCube()`, `convertInputREH()`, `getBinaryREH()` and `reh()`. This last function is the one that preprocesses the input given by the user which consists in: edgelist, riskset and covariates. intereventTime variable and covariates input still need to be preprocessed via specific utility functions;
    - _remstimate.cpp_ : contains `remDerivatives()` (which returns the value of loglikelihood, gradient, hessian at a specific parameter value), `lpd()` (log-pointwise density), utility functions for the fast method (`cube2matrix()`, `getUniqueVectors()`, `computeTimes()`, `computeOccurrencies()`) which is run with the function `remDerivativesFast()`;
    - Since `compute_stats()` is not an exported function in `remstats`, _getStats.R_ / _compute_stats.cpp_ / _compute_stats.h_ are temporary files so as to calculate statistics, run the estimation and compare estimates with `relevent::rem()`.  `getStats()` is the alias of `remstats::remstats()` with some modifications at the stage of the preprocessing of the network;
* _07/07/2020_ :
    - `reh()` will be the only preprocessing function and it is coded in Rcpp (see _reh.cpp_ file) whereas the R function and the _reh.R_ file are removed;
    - utility functions called inside `reh()` are added inside the _reh.cpp_ file, before the `reh()` function itself;
    - _messages.cpp_ becomes a header file _messages.h_ and the aim/content remains the same;
* _03/07/2020_ :
    - _reh.h_ changed to _reh.cpp_ and it contains utility functions used in _reh.R_ within the R function `reh(...)` ;
    - _messages.cpp_ will contain functions `errorMessage(cond)` and `warningMessage(cond)` that will return appropriate error/warning messages according to the _cond_ argument;
* _11/06/2020_ :
    - Created _reh.h_ were the utility functions to preprocess data will be developed;
    - _remstimateBoost.h_ will contain the routines that speed up the computation of the loglikelihood and its first and second derivatives;
* _20/04/2020_ :
    - created repository with first commit;
    - package only contains three functions: `remCpp(...)` (uses optim to find the maximum likelihood estimates of REM),
    `nllik(...)` (returns the negative log-likelihood value for an observed event sequence, by specifying a vector of parameters and statistics),
    `lpd(...)` (calculates the same as `nllik(...)` but only for a specific time point and without taking the negative of the value).
