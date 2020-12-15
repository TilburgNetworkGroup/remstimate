# NEWS (last update on December 11, 2020)

_20/04/2020_ :
* created repository with first commit;
* package only contains three functions: `remCpp(...)` (uses optim to find the maximum likelihood estimates of REM),
`nllik(...)` (returns the negative log-likelihood value for an observed event sequence, by specifying a vector of parameters and statistics),
`lpd(...)` (calculates the same as `nllik(...)` but only for a specific time point and without taking the negative of the value);

_11/06/2020_ :
* Created _reh.h_ were the utility functions to preprocess data will be developed;
* _remstimateBoost.h_ will contain the routines that speed up the computation of the loglikelihood and its first and second derivatives.

_03/07/2020_ :
*  _reh.h_ changed to _reh.cpp_ and it contains utility functions used in _reh.R_ within the R function `reh(...)` ;
* _messages.cpp_ will contain functions `errorMessage(cond)` and `warningMessage(cond)` that will return appropriate error/warning messages according to the _cond_ argument.

_07/07/2020_ :
*  `reh()` will be the only preprocessing function and it is coded in Rcpp (see _reh.cpp_ file) whereas the R function and the _reh.R_ file are removed;
* utility functions called inside `reh()` are added inside the _reh.cpp_ file, before the `reh()` function itself;
* _messages.cpp_ becomes a header file _messages.h_ and the aim/content remains the same.

_04/09/2020_ :
*  _messages_ becomes again an Rcpp file. This exstension appears to suit better the intent of the content/aim of error and warning messages;
* _remstimate.R_ contains the main function `remstimate()` which is aimed to run either a MLE or a Bayesian approach by using different optimization/methods. It also includes a switch to the fast method to compute the likelihood. The fast method is run if the actual improvement (percentage of improvement) is higher than a threshold set by the user (default 0.5);
* _reh.cpp_ : contains `getRisksetMatrix()`, `getRisksetCube()`, `convertInputREH()`, `getBinaryREH()` and `reh()`. This last function is the one that preprocesses the input given by the user which consists in: edgelist, riskset and covariates. intereventTime variable and covariates input still need to be preprocessed via specific utility functions;
* _remstimate.cpp_ : contains `remDerivatives()` (which returns the value of loglikelihood, gradient, hessian at a specific parameter value), `lpd()` (log-pointwise density), utility functions for the fast method (`cube2matrix()`, `getUniqueVectors()`, `computeTimes()`, `computeOccurrencies()`) which is run with the function `remDerivativesFast()`;
* Since `compute_stats()` is not an exported function in `remstats`, _getStats.R_ / _compute_stats.cpp_ / _compute_stats.h_ are temporary files so as to calculate statistics, run the estimation and compare estimates with `relevent::rem()`.  `getStats()` is the alias of `remstats::remstats()` with some modifications at the stage of the preprocessing of the network;

_11/12/2020_ :
* Methods working on `remstimate()` are: _MLE_, _GD_, _GDADAM_, _BSIR_, _HMC_. However the output lacks of a structure attributes and methods.
