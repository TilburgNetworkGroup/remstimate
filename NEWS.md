# NEWS (last update on April 20, 2020)

_04/20/2020_ :
* created repository with first commit;
* package only contains three functions: _remCpp_ (uses optim to find the maximum likelihood estimates of REM),
_nllik_ (returns the negative log-likelihood value for an observed event sequence, by specifying a vector of parameters and statistics),
_lpd_ (calculates the same as _nllik_ but only for a specific time point and without taking the negative of the value);
