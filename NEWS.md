# NEWS (last update on July 7, 2020)

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
*  `reh()` will be the only preprocessing function and it is coded in Rcpp (see reh.cpp file) whereas the R function and the _reh.R_ file are removed;
* utility functions called inside `reh()` are added inside the _reh.cpp_ file, before the `reh()` function itself;
* _messages.cpp_ becomes an header file _messages.h_ and the aim/content remains the same.
