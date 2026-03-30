#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

// remDerivativesSampled
//
// Importance-weighted sampled likelihood for the tie-oriented REM.
//
// INTERVAL timing log-likelihood at event m:
//   log L_m = x_{case} β  -  Δt_m * Σ_{s} exp(x_s β) / π_s
//
// ORDINAL timing log-likelihood at event m:
//   log L_m = x_{case} β  -  log( Σ_{s} exp(x_s β) / π_s )
//
// Both sums are importance-weighted: cases have π_s = 1,
// controls have π_s = samp_prob[m,s].
// When samp_num = D and all π_s = 1, both recover the full likelihood.
//
// [[Rcpp::export]]
Rcpp::List remDerivativesSampled(const arma::vec &pars,
                                  const arma::cube &stats,
                                  const Rcpp::List &case_pos,
                                  const arma::mat &samp_prob,
                                  const arma::vec &interevent_time,
                                  bool ordinal  = false,
                                  bool gradient = true,
                                  bool hessian  = true,
                                  int  ncores   = 1) {

  arma::uword M = stats.n_slices;
  arma::uword S = stats.n_rows;
  arma::uword P = stats.n_cols;
  arma::uword m, s, k, l;

  arma::vec loglik(M, arma::fill::zeros);
  arma::mat grad(P, M, arma::fill::zeros);
  arma::cube hess(P, P, M, arma::fill::zeros);

  #ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(ncores);
  #pragma omp parallel for if(ncores>1) private(m,s,k,l) \
    shared(M,S,P,loglik,grad,hess,pars,stats,case_pos,samp_prob,interevent_time,ordinal,gradient,hessian)
  #endif
  for (m = 0; m < M; m++) {
    arma::mat stats_m = stats.slice(m).t(); // [P x S]
    arma::vec log_lambda = stats_m.t() * pars; // [S]

    arma::uvec cases = Rcpp::as<arma::uvec>(case_pos[m]);
    arma::uword n_cases = cases.n_elem;

    double loglik_m = 0.0;
    arma::vec grad_m(P, arma::fill::zeros);
    arma::mat hess_m(P, P, arma::fill::zeros);

    if (n_cases == 0) { loglik(m) = 0.0; continue; }

    // (1) Numerator: sum log_lambda for cases
    loglik_m = arma::accu(log_lambda(cases));
    if (gradient) grad_m += arma::sum(stats_m.cols(cases), 1);

    // (2) Build importance weights: 1 for cases, 1/pi_s for controls
    arma::vec weights(S, arma::fill::ones);
    for (s = 0; s < S; s++) {
      bool is_case = false;
      for (arma::uword c = 0; c < n_cases; c++) {
        if (s == cases(c)) { is_case = true; break; }
      }
      if (!is_case) {
        double pi_s = samp_prob(m, s);
        weights(s) = (pi_s > 0.0) ? 1.0 / pi_s : 0.0;
      }
    }

    // lambda_s = exp(x_s β) / pi_s
    arma::vec lambda_s = arma::exp(log_lambda) % weights;
    double wtd_sum = arma::accu(lambda_s);
    if (wtd_sum <= 0.0) wtd_sum = 1e-300;
    arma::vec wtd_stats = stats_m * lambda_s; // [P]

    if (!ordinal) {
      // INTERVAL: subtract Δt_m * Σ lambda_s
      double dt_m = interevent_time(m);
      loglik_m -= dt_m * wtd_sum;
      if (gradient) grad_m -= dt_m * wtd_stats;
      if (hessian) {
        for (k = 0; k < P; k++) {
          for (l = k; l < P; l++) {
            double h_kl = 0.0;
            for (s = 0; s < S; s++)
              h_kl -= dt_m * stats_m(k,s) * stats_m(l,s) * lambda_s(s);
            hess_m(k,l) = h_kl;
            hess_m(l,k) = h_kl;
          }
        }
      }
    } else {
      // ORDINAL: subtract log(Σ lambda_s)
      loglik_m -= std::log(wtd_sum);
      if (gradient) grad_m -= wtd_stats / wtd_sum;
      if (hessian) {
        for (k = 0; k < P; k++) {
          for (l = k; l < P; l++) {
            double h_kl = 0.0;
            for (s = 0; s < S; s++)
              h_kl -= stats_m(k,s) * stats_m(l,s) * lambda_s(s);
            h_kl /= wtd_sum;
            h_kl += (wtd_stats(k)/wtd_sum) * (wtd_stats(l)/wtd_sum);
            hess_m(k,l) = h_kl;
            hess_m(l,k) = h_kl;
          }
        }
      }
    }

    loglik(m) = loglik_m;
    grad.col(m) = grad_m;
    hess.slice(m) = hess_m;
  }

  if (gradient && !hessian) {
    return Rcpp::List::create(
      Rcpp::Named("value")    = -arma::accu(loglik),
      Rcpp::Named("gradient") = -arma::sum(grad, 1));
  } else if (!gradient && !hessian) {
    return Rcpp::List::create(Rcpp::Named("value") = -arma::accu(loglik));
  } else {
    arma::cube H = -arma::sum(hess, 2);
    return Rcpp::List::create(
      Rcpp::Named("value")    = -arma::accu(loglik),
      Rcpp::Named("gradient") = -arma::sum(grad, 1),
      Rcpp::Named("hessian")  = H.slice(0));
  }
}
