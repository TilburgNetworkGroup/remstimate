# NEWS (last update on June 11, 2020)

_20/04/2020_ :
* created repository with first commit;
* package only contains three functions: _remCpp_ (uses optim to find the maximum likelihood estimates of REM),
_nllik_ (returns the negative log-likelihood value for an observed event sequence, by specifying a vector of parameters and statistics),
_lpd_ (calculates the same as _nllik_ but only for a specific time point and without taking the negative of the value);

_11/06/2020_ :
* Created _reh.h_ were the utility functions to preprocess data will be developed;
* _remstimateBoost.cpp_ will contain the routines that speed up the computation of the loglikelihood and its first and second derivatives.
