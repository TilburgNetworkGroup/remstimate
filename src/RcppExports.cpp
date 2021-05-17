// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// warningMessage
std::string warningMessage(int cond);
RcppExport SEXP _remstimate_warningMessage(SEXP condSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type cond(condSEXP);
    rcpp_result_gen = Rcpp::wrap(warningMessage(cond));
    return rcpp_result_gen;
END_RCPP
}
// errorMessage
std::string errorMessage(int cond);
RcppExport SEXP _remstimate_errorMessage(SEXP condSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type cond(condSEXP);
    rcpp_result_gen = Rcpp::wrap(errorMessage(cond));
    return rcpp_result_gen;
END_RCPP
}
// getUniqueVectors
arma::mat getUniqueVectors(arma::cube stats);
RcppExport SEXP _remstimate_getUniqueVectors(SEXP statsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type stats(statsSEXP);
    rcpp_result_gen = Rcpp::wrap(getUniqueVectors(stats));
    return rcpp_result_gen;
END_RCPP
}
// computeTimes
arma::vec computeTimes(const arma::mat& unique_vectors_stats, const arma::uword& M, const arma::cube& stats, const arma::vec& intereventTime);
RcppExport SEXP _remstimate_computeTimes(SEXP unique_vectors_statsSEXP, SEXP MSEXP, SEXP statsSEXP, SEXP intereventTimeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type unique_vectors_stats(unique_vectors_statsSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type stats(statsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type intereventTime(intereventTimeSEXP);
    rcpp_result_gen = Rcpp::wrap(computeTimes(unique_vectors_stats, M, stats, intereventTime));
    return rcpp_result_gen;
END_RCPP
}
// computeOccurrencies
arma::vec computeOccurrencies(const Rcpp::DataFrame& edgelist, const arma::ucube& risksetCube, const arma::uword& M, const arma::mat& unique_vectors_stats, const arma::cube& stats);
RcppExport SEXP _remstimate_computeOccurrencies(SEXP edgelistSEXP, SEXP risksetCubeSEXP, SEXP MSEXP, SEXP unique_vectors_statsSEXP, SEXP statsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type edgelist(edgelistSEXP);
    Rcpp::traits::input_parameter< const arma::ucube& >::type risksetCube(risksetCubeSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type unique_vectors_stats(unique_vectors_statsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type stats(statsSEXP);
    rcpp_result_gen = Rcpp::wrap(computeOccurrencies(edgelist, risksetCube, M, unique_vectors_stats, stats));
    return rcpp_result_gen;
END_RCPP
}
// remDerivativesStandard
Rcpp::List remDerivativesStandard(const arma::vec& pars, const arma::cube& stats, const arma::mat& event_binary, const arma::vec& interevent_time, bool ordinal, int ncores, bool gradient, bool hessian);
RcppExport SEXP _remstimate_remDerivativesStandard(SEXP parsSEXP, SEXP statsSEXP, SEXP event_binarySEXP, SEXP interevent_timeSEXP, SEXP ordinalSEXP, SEXP ncoresSEXP, SEXP gradientSEXP, SEXP hessianSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type stats(statsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type event_binary(event_binarySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type interevent_time(interevent_timeSEXP);
    Rcpp::traits::input_parameter< bool >::type ordinal(ordinalSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< bool >::type gradient(gradientSEXP);
    Rcpp::traits::input_parameter< bool >::type hessian(hessianSEXP);
    rcpp_result_gen = Rcpp::wrap(remDerivativesStandard(pars, stats, event_binary, interevent_time, ordinal, ncores, gradient, hessian));
    return rcpp_result_gen;
END_RCPP
}
// remDerivativesFast
Rcpp::List remDerivativesFast(const arma::vec& pars, const arma::vec& times_r, const arma::vec& occurrencies_r, const arma::mat& unique_vectors_stats, bool gradient, bool hessian);
RcppExport SEXP _remstimate_remDerivativesFast(SEXP parsSEXP, SEXP times_rSEXP, SEXP occurrencies_rSEXP, SEXP unique_vectors_statsSEXP, SEXP gradientSEXP, SEXP hessianSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type times_r(times_rSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type occurrencies_r(occurrencies_rSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type unique_vectors_stats(unique_vectors_statsSEXP);
    Rcpp::traits::input_parameter< bool >::type gradient(gradientSEXP);
    Rcpp::traits::input_parameter< bool >::type hessian(hessianSEXP);
    rcpp_result_gen = Rcpp::wrap(remDerivativesFast(pars, times_r, occurrencies_r, unique_vectors_stats, gradient, hessian));
    return rcpp_result_gen;
END_RCPP
}
// remDerivativesSenderRates
Rcpp::List remDerivativesSenderRates(const arma::vec& pars, const arma::cube& stats, const arma::cube& risksetCube, const arma::mat& event_binary, const arma::vec& interevent_time, const arma::mat& edgelist);
RcppExport SEXP _remstimate_remDerivativesSenderRates(SEXP parsSEXP, SEXP statsSEXP, SEXP risksetCubeSEXP, SEXP event_binarySEXP, SEXP interevent_timeSEXP, SEXP edgelistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type stats(statsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type risksetCube(risksetCubeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type event_binary(event_binarySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type interevent_time(interevent_timeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type edgelist(edgelistSEXP);
    rcpp_result_gen = Rcpp::wrap(remDerivativesSenderRates(pars, stats, risksetCube, event_binary, interevent_time, edgelist));
    return rcpp_result_gen;
END_RCPP
}
// remDerivativesReceiverChoice
Rcpp::List remDerivativesReceiverChoice(const arma::vec& pars, const arma::cube& stats, const arma::cube& risksetCube, const arma::mat& event_binary, const arma::vec& interevent_time, const arma::mat& edgelist, int N);
RcppExport SEXP _remstimate_remDerivativesReceiverChoice(SEXP parsSEXP, SEXP statsSEXP, SEXP risksetCubeSEXP, SEXP event_binarySEXP, SEXP interevent_timeSEXP, SEXP edgelistSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type stats(statsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type risksetCube(risksetCubeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type event_binary(event_binarySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type interevent_time(interevent_timeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type edgelist(edgelistSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(remDerivativesReceiverChoice(pars, stats, risksetCube, event_binary, interevent_time, edgelist, N));
    return rcpp_result_gen;
END_RCPP
}
// remDerivatives
Rcpp::List remDerivatives(const arma::vec& pars, const arma::cube& stats, const arma::mat& event_binary, const arma::vec& interevent_time, std::string model, bool ordinal, int ncores, bool fast, bool gradient, bool hessian);
RcppExport SEXP _remstimate_remDerivatives(SEXP parsSEXP, SEXP statsSEXP, SEXP event_binarySEXP, SEXP interevent_timeSEXP, SEXP modelSEXP, SEXP ordinalSEXP, SEXP ncoresSEXP, SEXP fastSEXP, SEXP gradientSEXP, SEXP hessianSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type stats(statsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type event_binary(event_binarySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type interevent_time(interevent_timeSEXP);
    Rcpp::traits::input_parameter< std::string >::type model(modelSEXP);
    Rcpp::traits::input_parameter< bool >::type ordinal(ordinalSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< bool >::type fast(fastSEXP);
    Rcpp::traits::input_parameter< bool >::type gradient(gradientSEXP);
    Rcpp::traits::input_parameter< bool >::type hessian(hessianSEXP);
    rcpp_result_gen = Rcpp::wrap(remDerivatives(pars, stats, event_binary, interevent_time, model, ordinal, ncores, fast, gradient, hessian));
    return rcpp_result_gen;
END_RCPP
}
// GD
Rcpp::List GD(const arma::vec& pars, const arma::cube& stats, const arma::mat& event_binary, const arma::vec& interevent_time, std::string model, bool ordinal, int ncores, bool fast, int epochs, double learning_rate);
RcppExport SEXP _remstimate_GD(SEXP parsSEXP, SEXP statsSEXP, SEXP event_binarySEXP, SEXP interevent_timeSEXP, SEXP modelSEXP, SEXP ordinalSEXP, SEXP ncoresSEXP, SEXP fastSEXP, SEXP epochsSEXP, SEXP learning_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type stats(statsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type event_binary(event_binarySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type interevent_time(interevent_timeSEXP);
    Rcpp::traits::input_parameter< std::string >::type model(modelSEXP);
    Rcpp::traits::input_parameter< bool >::type ordinal(ordinalSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< bool >::type fast(fastSEXP);
    Rcpp::traits::input_parameter< int >::type epochs(epochsSEXP);
    Rcpp::traits::input_parameter< double >::type learning_rate(learning_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(GD(pars, stats, event_binary, interevent_time, model, ordinal, ncores, fast, epochs, learning_rate));
    return rcpp_result_gen;
END_RCPP
}
// GDADAM
Rcpp::List GDADAM(const arma::vec& pars, const arma::cube& stats, const arma::mat& event_binary, const arma::vec& interevent_time, std::string model, bool ordinal, int ncores, bool fast, int epochs, double learning_rate, double beta1, double beta2, double eta);
RcppExport SEXP _remstimate_GDADAM(SEXP parsSEXP, SEXP statsSEXP, SEXP event_binarySEXP, SEXP interevent_timeSEXP, SEXP modelSEXP, SEXP ordinalSEXP, SEXP ncoresSEXP, SEXP fastSEXP, SEXP epochsSEXP, SEXP learning_rateSEXP, SEXP beta1SEXP, SEXP beta2SEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type stats(statsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type event_binary(event_binarySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type interevent_time(interevent_timeSEXP);
    Rcpp::traits::input_parameter< std::string >::type model(modelSEXP);
    Rcpp::traits::input_parameter< bool >::type ordinal(ordinalSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< bool >::type fast(fastSEXP);
    Rcpp::traits::input_parameter< int >::type epochs(epochsSEXP);
    Rcpp::traits::input_parameter< double >::type learning_rate(learning_rateSEXP);
    Rcpp::traits::input_parameter< double >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< double >::type beta2(beta2SEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(GDADAM(pars, stats, event_binary, interevent_time, model, ordinal, ncores, fast, epochs, learning_rate, beta1, beta2, eta));
    return rcpp_result_gen;
END_RCPP
}
// logPostHMC
double logPostHMC(const arma::vec& meanPrior, const arma::mat& sigmaPrior, const arma::vec& pars, const arma::cube& stats, const arma::mat& event_binary, const arma::vec& interevent_time, std::string model, bool ordinal, int ncores, bool fast);
RcppExport SEXP _remstimate_logPostHMC(SEXP meanPriorSEXP, SEXP sigmaPriorSEXP, SEXP parsSEXP, SEXP statsSEXP, SEXP event_binarySEXP, SEXP interevent_timeSEXP, SEXP modelSEXP, SEXP ordinalSEXP, SEXP ncoresSEXP, SEXP fastSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type meanPrior(meanPriorSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigmaPrior(sigmaPriorSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type stats(statsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type event_binary(event_binarySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type interevent_time(interevent_timeSEXP);
    Rcpp::traits::input_parameter< std::string >::type model(modelSEXP);
    Rcpp::traits::input_parameter< bool >::type ordinal(ordinalSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< bool >::type fast(fastSEXP);
    rcpp_result_gen = Rcpp::wrap(logPostHMC(meanPrior, sigmaPrior, pars, stats, event_binary, interevent_time, model, ordinal, ncores, fast));
    return rcpp_result_gen;
END_RCPP
}
// logPostGradientHMC
arma::vec logPostGradientHMC(const arma::vec& meanPrior, const arma::mat& sigmaPrior, const arma::vec& pars, const arma::cube& stats, const arma::mat& event_binary, const arma::vec& interevent_time, std::string model, bool ordinal, int ncores, bool fast);
RcppExport SEXP _remstimate_logPostGradientHMC(SEXP meanPriorSEXP, SEXP sigmaPriorSEXP, SEXP parsSEXP, SEXP statsSEXP, SEXP event_binarySEXP, SEXP interevent_timeSEXP, SEXP modelSEXP, SEXP ordinalSEXP, SEXP ncoresSEXP, SEXP fastSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type meanPrior(meanPriorSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigmaPrior(sigmaPriorSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type stats(statsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type event_binary(event_binarySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type interevent_time(interevent_timeSEXP);
    Rcpp::traits::input_parameter< std::string >::type model(modelSEXP);
    Rcpp::traits::input_parameter< bool >::type ordinal(ordinalSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< bool >::type fast(fastSEXP);
    rcpp_result_gen = Rcpp::wrap(logPostGradientHMC(meanPrior, sigmaPrior, pars, stats, event_binary, interevent_time, model, ordinal, ncores, fast));
    return rcpp_result_gen;
END_RCPP
}
// iterHMC
arma::vec iterHMC(arma::uword L, double epsilon, const arma::vec& meanPrior, const arma::mat& sigmaPrior, const arma::vec& pars, const arma::cube& stats, const arma::mat& event_binary, const arma::vec& interevent_time, std::string model, bool ordinal, int ncores, bool fast);
RcppExport SEXP _remstimate_iterHMC(SEXP LSEXP, SEXP epsilonSEXP, SEXP meanPriorSEXP, SEXP sigmaPriorSEXP, SEXP parsSEXP, SEXP statsSEXP, SEXP event_binarySEXP, SEXP interevent_timeSEXP, SEXP modelSEXP, SEXP ordinalSEXP, SEXP ncoresSEXP, SEXP fastSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type meanPrior(meanPriorSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigmaPrior(sigmaPriorSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type stats(statsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type event_binary(event_binarySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type interevent_time(interevent_timeSEXP);
    Rcpp::traits::input_parameter< std::string >::type model(modelSEXP);
    Rcpp::traits::input_parameter< bool >::type ordinal(ordinalSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< bool >::type fast(fastSEXP);
    rcpp_result_gen = Rcpp::wrap(iterHMC(L, epsilon, meanPrior, sigmaPrior, pars, stats, event_binary, interevent_time, model, ordinal, ncores, fast));
    return rcpp_result_gen;
END_RCPP
}
// burninHMC
arma::cube burninHMC(const arma::cube& samples, arma::uword burnin, arma::uword thin);
RcppExport SEXP _remstimate_burninHMC(SEXP samplesSEXP, SEXP burninSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(burninHMC(samples, burnin, thin));
    return rcpp_result_gen;
END_RCPP
}
// HMC
arma::cube HMC(arma::mat pars_init, arma::uword nsim, arma::uword nchains, arma::uword burnin, const arma::vec& meanPrior, const arma::mat& sigmaPrior, const arma::cube& stats, const arma::mat& event_binary, const arma::vec& interevent_time, std::string model, bool ordinal, int ncores, bool fast, arma::uword thin, arma::uword L, double epsilon);
RcppExport SEXP _remstimate_HMC(SEXP pars_initSEXP, SEXP nsimSEXP, SEXP nchainsSEXP, SEXP burninSEXP, SEXP meanPriorSEXP, SEXP sigmaPriorSEXP, SEXP statsSEXP, SEXP event_binarySEXP, SEXP interevent_timeSEXP, SEXP modelSEXP, SEXP ordinalSEXP, SEXP ncoresSEXP, SEXP fastSEXP, SEXP thinSEXP, SEXP LSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type pars_init(pars_initSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type nchains(nchainsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type meanPrior(meanPriorSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigmaPrior(sigmaPriorSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type stats(statsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type event_binary(event_binarySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type interevent_time(interevent_timeSEXP);
    Rcpp::traits::input_parameter< std::string >::type model(modelSEXP);
    Rcpp::traits::input_parameter< bool >::type ordinal(ordinalSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< bool >::type fast(fastSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(HMC(pars_init, nsim, nchains, burnin, meanPrior, sigmaPrior, stats, event_binary, interevent_time, model, ordinal, ncores, fast, thin, L, epsilon));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_remstimate_warningMessage", (DL_FUNC) &_remstimate_warningMessage, 1},
    {"_remstimate_errorMessage", (DL_FUNC) &_remstimate_errorMessage, 1},
    {"_remstimate_getUniqueVectors", (DL_FUNC) &_remstimate_getUniqueVectors, 1},
    {"_remstimate_computeTimes", (DL_FUNC) &_remstimate_computeTimes, 4},
    {"_remstimate_computeOccurrencies", (DL_FUNC) &_remstimate_computeOccurrencies, 5},
    {"_remstimate_remDerivativesStandard", (DL_FUNC) &_remstimate_remDerivativesStandard, 8},
    {"_remstimate_remDerivativesFast", (DL_FUNC) &_remstimate_remDerivativesFast, 6},
    {"_remstimate_remDerivativesSenderRates", (DL_FUNC) &_remstimate_remDerivativesSenderRates, 6},
    {"_remstimate_remDerivativesReceiverChoice", (DL_FUNC) &_remstimate_remDerivativesReceiverChoice, 7},
    {"_remstimate_remDerivatives", (DL_FUNC) &_remstimate_remDerivatives, 10},
    {"_remstimate_GD", (DL_FUNC) &_remstimate_GD, 10},
    {"_remstimate_GDADAM", (DL_FUNC) &_remstimate_GDADAM, 13},
    {"_remstimate_logPostHMC", (DL_FUNC) &_remstimate_logPostHMC, 10},
    {"_remstimate_logPostGradientHMC", (DL_FUNC) &_remstimate_logPostGradientHMC, 10},
    {"_remstimate_iterHMC", (DL_FUNC) &_remstimate_iterHMC, 12},
    {"_remstimate_burninHMC", (DL_FUNC) &_remstimate_burninHMC, 3},
    {"_remstimate_HMC", (DL_FUNC) &_remstimate_HMC, 16},
    {NULL, NULL, 0}
};

RcppExport void R_init_remstimate(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
