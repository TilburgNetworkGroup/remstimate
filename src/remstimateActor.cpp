#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;

//' remDerivativesActor
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values
//'
//' @param pars is a vector of parameters The first P_d elements are for the dyadic effect parameters and the next P_s elements are for the sender effect parameters (note: the order must be aligned with the column order in 'stats_d' and 'stats_s' respectively)
//' @param stats_d is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param stats_s is cube of M slices. Each slice is a matrix of dimensions N*U with statistics of interest by column and senders by row.
//' @param risksetCube
//' @param event_binary is a matrix [M*D] of 1/0/-1 : 1 indicating the observed dyad and 0 (-1) the non observed dyads that could have (have not) occurred.
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param edgelist, output from remify, (note: indices of the actors must start from 0)
//'
//'
//' @return list of values: loglik, grad, fisher
//'
// [[Rcpp::export]]
Rcpp::List remDerivativesActor(
    const arma::vec &pars,
    const arma::cube &stats_d,
    const arma::cube &stats_s,
    const arma::cube &risksetCube,
    const arma::mat &event_binary,
    const arma::vec &interevent_time,
	const arma::mat &edgelist
  ){

    arma::uword P_d = stats_d.n_cols; // number of parameters for dyad stats
    arma::uword P_s = stats_s.n_cols; // number of parameters for send stats
    arma::uword D = event_binary.n_cols; // number of dyads
    arma::uword M = edgelist.n_rows; // number of events
    arma::uword N = stats_s.n_rows;//Number of actors
	arma::vec pars_d = pars(arma::span(0,(P_d-1)));
	arma::vec pars_s = pars(arma::span(P_d,(P_d+P_s-1)));

	// variables for internal computation
    arma::vec lambda_d(D) ;
    arma::vec lambda_s(N) ;
	arma::vec weighted_sum_current_event(P_s);
	arma::vec expected_stat_current_event(P_d);
	arma::mat fisher_current_event_s(P_s,P_s);
	arma::mat fisher_current_event_d(P_d,P_d);
	arma::uword dyad =0;
	double denom = 0;

	//output
    double loglik = 0.0;
    arma::mat fisher_d(P_d,P_d,arma::fill::zeros);
	arma::mat fisher_s(P_s,P_s,arma::fill::zeros);
    arma::vec grad_d(P_d,arma::fill::zeros);
    arma::vec grad_s(P_s,arma::fill::zeros);


	//loop through all events
    for(arma::uword m = 0; m < M; m++){
        arma::mat stats_d_m = stats_d.slice(m).t(); // dimensions : [P_d*D] we want to access dyads by column
        arma::mat stats_s_m = stats_s.slice(m).t(); // dimensions: [P_s * N]

		//this is exp(beta^T X) dot product
        lambda_d = arma::exp(stats_d_m.t() * pars_d);
        lambda_s = arma::exp(stats_s_m.t() * pars_s);

		//actors in remify edgelist are 0 indexed hence no -1
		int sender = edgelist(m,1);
        int receiver = edgelist(m,2);

		//reset internal variables
		fisher_current_event_s.zeros();
		fisher_current_event_d.zeros();
		expected_stat_current_event.zeros();
		denom = 0;

		//loop thought all actors
		for(int i = 0;i<N; i++){
		    dyad = risksetCube(sender,i,0);
			if(i!=sender && event_binary(m,dyad)!= -1){
				dyad = risksetCube(sender,i,0);
				//loglik
				denom += lambda_d(dyad); // exp(param_d * X_sender_i) (4)
				fisher_current_event_d += lambda_d(dyad)*(stats_d_m.col(dyad) * stats_d_m.col(dyad).t());

				expected_stat_current_event += lambda_d(dyad) * stats_d_m.col(dyad);
			}
		    fisher_current_event_s += lambda_s(i)*(stats_s_m.col(i) * stats_s_m.col(i).t() );
		}

		dyad = risksetCube(sender,receiver,0);

		loglik += std::log(lambda_s(sender)) + std::log(lambda_d(dyad)); //(1) params_s^T * X_sender + (3) params_d^T * X_sender_receiver

		loglik -=  sum(lambda_s) * interevent_time(m); // (2);
		loglik -= std::log(denom);

		grad_s += stats_s_m(sender); //(5)
		weighted_sum_current_event = stats_s_m * lambda_s; //exp(params_s * X_i)*X_i
		grad_s -= interevent_time(m) * weighted_sum_current_event; //(6)

		grad_d += stats_d_m.col(dyad); //(7)
		expected_stat_current_event /= denom; //(8)
		grad_d -= expected_stat_current_event;

		fisher_s += interevent_time(m) * fisher_current_event_s;

		fisher_current_event_d /= denom;
		fisher_current_event_d -= (expected_stat_current_event * expected_stat_current_event.t());
		fisher_d += fisher_current_event_d;
    }

	arma::vec grad(P_d+P_s,arma::fill::zeros);
	grad(arma::span(0,(P_d-1))) = grad_d;
	grad(arma::span(P_d,(P_d+P_s-1) )) = grad_s;

	arma::mat fisher(P_s+P_d,P_s+P_d,arma::fill::zeros);
	fisher(arma::span(0,P_d-1),arma::span(0,P_d-1)) = fisher_d;
	fisher(arma::span(P_d,(P_d+P_s-1)),arma::span(P_d,(P_d+P_s-1))) = fisher_s;

    return Rcpp::List::create(Rcpp::Named("value") = loglik, Rcpp::Named("gradient") = grad, Rcpp::Named("fisher") = fisher);
}
