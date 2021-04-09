#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;

//' remDerivativesActor
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values
//'
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param event_binary is a matrix [M*D] of 1/0/-1 : 1 indicating the observed dyad and 0 (-1) the non observed dyads that could have (have not) occurred.
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param gradient boolean true/false whether to return gradient value
//' @param hessian boolean true/false whether to return hessian value
//'
//' @return list of values: loglik, gradient, hessian
//'
// [[Rcpp::export]]
Rcpp::List remDerivativesActor(
    const arma::vec &pars_d,
    const arma::vec &pars_s,
    const arma::cube &stats_d,
    const arma::cube &stats_s,
    const arma::cube &risksetCube,
    const arma::mat &event_binary,
    const arma::vec &interevent_time,
	const arma::mat &edgelist
  ){
    arma::uword P_d = pars_d.n_elem; // number of parameters
    arma::uword P_s = pars_s.n_elem; // number of parameters for send stats
    arma::uword D = event_binary.n_cols; // number of dyads
    arma::uword M = edgelist.n_rows; // number of events
    arma::uword N = stats_s.n_rows;//Number of actors

    cout <<"M"<<M<<"D:"<<D<<"N:"<<N<<"P_d:"<<P_d<<"P_s:"<<P_s<<endl;


	// variables for internal computation
    arma::vec lambda_d(D) ;
    arma::vec lambda_s(N) ;
	arma::vec weighted_sum_current_event(P_s);
	arma::vec expected_stat_current_event(P_d);
	arma::mat fisher_current_event_s(P_s,P_s);
	arma::mat fisher_current_event_d(P_d,P_d);
	arma::uword dyad =0; //dyad
	double denom = 0;
	//output
    double loglik = 0.0;
    arma::mat fisher_d(P_d,P_d,arma::fill::zeros);
	arma::mat fisher_s(P_s,P_s,arma::fill::zeros);
    arma::vec grad_d(P_d,arma::fill::zeros);
    arma::vec grad_s(P_s,arma::fill::zeros);

	//return Rcpp::List::create(Rcpp::Named("value") = loglik, Rcpp::Named("gradient_d") = grad_d, Rcpp::Named("fisher_d") = fisher_d,Rcpp::Named("gradient_s") = grad_s, Rcpp::Named("fisher_s") = fisher_s);

	//loop through all events
    //for(arma::uword m = 0; m < M; m++){
        int m = 5;
        cout<<"event:"<<m<<endl;
        arma::mat stats_d_m = stats_d.slice(m).t(); // dimensions : [P_d*D] we want to access dyads by column
        arma::mat stats_s_m = stats_s.slice(m).t(); // dimensions: [P_s * N]

		//std::cout <<"1"<<std::endl;

		//this is exp(beta^T X) dot product
        lambda_d = arma::exp(stats_d_m.t() * pars_d);
        lambda_s = arma::exp(stats_s_m.t() * pars_s);

		std::cout <<"2"<<std::endl;
        cout << size(edgelist)<<endl;

		//actors in remify edgelist are 0 indexed hence no -1
		int sender = edgelist(m,1);
        int receiver = edgelist(m,2);

		//reset internal variables
		fisher_current_event_s.zeros();
		fisher_current_event_d.zeros();
		expected_stat_current_event.zeros();
		denom = 0;

		//loop thought all actors
		cout << "sender: "<<sender <<" receiver: "<<receiver<<endl;

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
		//expected_stat_current_event = stats_d_m * lambda_d;//exp(param_d*X_sender_i)*X_sender_i
		expected_stat_current_event /= denom; //(8)
		grad_d -= expected_stat_current_event;

		//cout << "(7):"<<stats_d_m(0,dyad)<< "-"<<stats_d_m(1,dyad)<<endl;
		//cout <<"expected stat (8)"<<-expected_stat_current_event(0)<<" -  "<<-expected_stat_current_event(1)<<endl;

		fisher_s += interevent_time(m) * fisher_current_event_s;

		//arma::mat temp = expected_stat_current_event * expected_stat_current_event.t();

		fisher_current_event_d /= denom;

		//cout <<"fisher d"<<fisher_current_event_d(0,0)<<" - "<<fisher_current_event_d(0,1)<<" - "<<fisher_current_event_d(1,0)<<" - "<<fisher_current_event_d(1,1)<<endl;

		fisher_current_event_d -= (expected_stat_current_event * expected_stat_current_event.t());
		//cout <<"temp d"<<temp(0,0)<<" - "<<temp(0,1)<<" - "<<temp(1,0)<<" - "<<temp(1,1)<<endl;
		fisher_d += fisher_current_event_d;
    //}

    return Rcpp::List::create(Rcpp::Named("value") = loglik, Rcpp::Named("gradient_d") = grad_d, Rcpp::Named("fisher_d") = fisher_d,Rcpp::Named("gradient_s") = grad_s, Rcpp::Named("fisher_s") = fisher_s);
}

/***R

out <- remDerivativesActor(
    pars_d,
    pars_s,
    stats_d,
    stats_s,
    rsCube,
    reh.out$rehBinary,
    reh.out$intereventTime,
    as.matrix(reh.out$edgelist)
)

 */
