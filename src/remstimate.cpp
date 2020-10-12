#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <omp.h>
#include <RcppArmadilloExtensions/sample.h>
#include <typeinfo>
#include <map>
#include <iterator>
#include <string>
#include <remstats.h>


// compute_statsCpp (description is the same as in remstats)
//
// Calls functions in compute_effects.h to compute statistics, scales 
// statistics, and combines all computed statistics in an array. 
// 
//  effects: vector of length p with integers referring to effect types.
//  edgelist: matrix [time, sender/actor1, receiver/actor2, (type), 
// riskset position]
//  riskset: matrix [sender/actor1, receiver/actor2, (type)]
//  start: integer value referring to the first row + 1in the edgelist 
// for which effects have to be computed
//  stop: integer value referring to the last row + 1 in the edgelist 
// for which effects have to be computed
//  values: list of length p with matrices with exogenous information 
// for the exogenous effects and NULL for all other effects
//  scaling: vector of length p with integer values referring to the 
// type of scaling method that has to be applied to the statistic.
//  memory_value: vector of length p with numeric values referring to 
// the length of the window in which past events are considered for endogenous 
// effects
//  with_type: vector of length p with logical values indicating whether 
// effects have to be computed considering event types a dependent variable
//  event_weights: matrix with p columns where each column refers to the 
// weights of the events for a specific effect
//  equal_val: vector of length p 
//
//  statistics: array [timepoint x riskset position x statistic]
//
//[[Rcpp::export]]
arma::cube compute_statsCpp(const arma::vec& effects, const arma::mat& edgelist, 
    const arma::mat& riskset, int start, int stop, 
    const Rcpp::List& values, const arma::vec& scaling, 
    const arma::vec& memory_value, const arma::vec& with_type, 
    const arma::mat& event_weights, const arma::vec& equal_val) {

arma::cube out = remstats::compute_stats(effects, edgelist, riskset, start, stop, 
        values, scaling, memory_value, with_type, event_weights, equal_val);
        return out;

}

//' nLoglik
//'
//' function that returns the negative of the loglikelihood value at specific parameters' values
//' 
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param event_binary is a matrix [M*D] of 1/0/-1 : 1 indicating the observed dyad and 0 (-1) the non observed dyads that could have (have not) occurred.
//' @param interevent_time the time difference between the current time point and the previous event time.
//'
//' @return list of values: loglik, gradient, hessian
//'
//' @export
// [[Rcpp::export]]
double nLoglik(const arma::vec &pars,const arma::cube &stats, const arma::mat &event_binary, const arma::vec &interevent_time){
    arma::uword D = event_binary.n_cols;
    arma::uword M = event_binary.n_rows; // number of events

    arma::uword d,m; 
    arma::vec log_lambda(D,arma::fill::zeros) ;

    double loglik = 0.0;
    
    for(m = 0; m < M; m++){
        arma::mat stats_m = stats.slice(m); // dimensions : [D*U]
        log_lambda = stats_m.t() * pars;
        
        for(d = 0; d < D; d++){
            if(event_binary(m,d)!=-1){ // ignoring impossible events that are not in risk set                
                if(event_binary(m,d) == 1){ // if event occured
                    loglik += log_lambda.at(d);
                }               
                double dtelp = exp(log_lambda.at(d))*interevent_time.at(m);
                loglik -= dtelp;                
            }            
        }        
    }
    return -loglik;
}

//' remDerivatives
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values
//' 
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param event_binary is a matrix [M*D] of 1/0/-1 : 1 indicating the observed dyad and 0 (-1) the non observed dyads that could have (have not) occurred.
//' @param interevent_time the time difference between the current time point and the previous event time.
//'
//' @return list of values: loglik, gradient, hessian
//'
//' @export
// [[Rcpp::export]]
Rcpp::List remDerivatives(const arma::vec &pars,const arma::cube &stats, const arma::mat &event_binary, const arma::vec &interevent_time){
    arma::uword U = pars.n_elem; // number of parameters
    arma::uword D = event_binary.n_cols;
    arma::uword M = event_binary.n_rows; // number of events

    arma::uword d,m,l,k;
    arma::vec log_lambda(D,arma::fill::zeros) ;

    double loglik = 0.0;
    arma::mat hess(U,U,arma::fill::zeros);
    arma::vec grad(U,arma::fill::zeros);
    
    for(m = 0; m < M; m++){
        arma::mat stats_m = stats.slice(m); // dimensions : [D*U]
        log_lambda = stats_m.t() * pars;
        
        for(d = 0; d < D; d++){
            if(event_binary(m,d)!=-1){ // ignoring impossible events that are not in risk set                
                if(event_binary(m,d) == 1){ // if event occured
                    loglik += log_lambda.at(d);
                    grad += stats_m.col(d);
                }               
                double dtelp = exp(log_lambda.at(d))*interevent_time.at(m);
                loglik -= dtelp;                
                grad -= stats_m.col(d)*dtelp;               
                for(k = 0; k < U; k++){
                    for (l = k; l < U; l++){
                        hess(k,l) -= stats_m.at(l,d)*stats_m.at(k,d)*dtelp;
                        hess(l,k) = hess(k,l);
                    }
                }    
            }            
        }        
    }
    return Rcpp::List::create(Rcpp::Named("value") = -loglik, Rcpp::Named("gradient") = -grad, Rcpp::Named("hessian") = -hess);
}

//' lpd (Log-Pointwise Density of REM - to rewrite according to the 0/1/-1 event vector)
//'
//' @param pars is a vector of parameters (note: the order must be aligned witht the column order in 'stats')
//' @param stats is a matrix of dimensions n_dyads*variables with statistics of interest by column and dyads by row.
//' @param event is a vector of 1/0 : 1 indicating the observed dyad and 0 the non observed dyads.
//' @param interevent_time the time difference between the current time point and the previous event time.
//'
//' @return log-pointwise density value of a specific time point
//'
//' @export
// [[Rcpp::export]]
double lpd(arma::vec pars, arma::mat stats, arma::uvec event, double interevent_time){
        arma::uword n_dyads = event.n_elem;
        arma::uword i;
        arma::vec log_lambda = stats * pars;
        double lpd = 0.0;
        for(i = 0; i < n_dyads; i++){
            if(event(i) == 0){
                lpd -= exp(log_lambda(i))*interevent_time;
            }
            else{
                lpd += log_lambda(i)-exp(log_lambda(i))*interevent_time;
            }
        }
        return lpd;
    }

// /////////////////////////////////////////////////////////////////////////////////
// ///////////(BEGIN)     remstimateFAST routine functions      (BEGIN)///////////// 
// /////////////////////////////////////////////////////////////////////////////////

//' cube2matrix
//'
//' A function to rearrange the cube of statistics into a matrix.
//'
//' @param stats cube structure of dimensions [M*D*U] filled with statistics values. 
//'
//' @return matrix of dimensions [(M*D)*U]
arma::mat cube2matrix(arma::cube stats){
  arma::uword u,m;
  arma::uword new_row = 0;
  arma::mat out(stats.n_rows*stats.n_slices, stats.n_cols, arma::fill::zeros); 
  for(u = 0; u < stats.n_slices; u++){
    for(m = 0; m < stats.n_rows; m++){
      out.row(new_row) = stats.slice(u).row(m);
      new_row += 1;
    }
  }
  return out;
}

//' getUniqueVectors
//'
//' A function to retrieve only the unique vectors of statistics observed throught times points and dyads. This function is based on the result shown by the Appendix C in the paper 'Hierarchical models for relational event sequences', DuBois et al. 2013 (pp. 308-309).
//'
//' @param stats cube of statistics with dimensions [M*D*U]
//'
//' @return matrix with only unique vectors of statistics with dimensions [R*U]
//'
//' @export
// [[Rcpp::export]]
arma::mat getUniqueVectors(arma::cube stats){

  // transform the cube of statistics to matrix first
  arma::mat A = cube2matrix(stats);
  arma::uword a,a2compare; 
  arma::uvec indices(A.n_rows, arma::fill::zeros);
  
  for(a = 0; a < (A.n_rows-1); a++){
    for(a2compare = (a+1); a2compare < A.n_rows; a2compare++){
      if(arma::approx_equal(A.row(a2compare), A.row(a), "absdiff", 0.000001)){
         indices(a2compare) = 1;
         break;
      }
    }
  }

  return A.rows(arma::find(indices == 0));
}

//' computeTimes 
//'
//' A function to compute the sum of interevent times for those vector of statistics that occurre more than once (output of getUniqueVectors()). This function is based on the result shown by the Appendix C in the paper 'Hierarchical models for relational event sequences', DuBois et al. 2013 (pp. 308-309).
//'
//' @param unique_vectors_stats matrix of unique vectors of statistics (output of getUniqueVectors()).
//' @param M number of observed relational events.
//' @param stats array of statistics with dimensons [D*U*M].
//' @param intereventTime vector of time differences between two subsequent time points (i.d., waiting time between t[m] and t[m-1]).
//'
//' @return vector of sum of interevent times per each unique_vector_stats element
//'
//' @export
// [[Rcpp::export]]
arma::vec computeTimes(const arma::mat& unique_vectors_stats, const arma::uword& M, const arma::cube& stats, const arma::vec& intereventTime){
 
  arma::uword R = unique_vectors_stats.n_rows; // number of unique vector of statistics
  arma::uword D = stats.n_rows; // number of dyads

  arma::uword m,r,d; 
  arma::mat mat_counts(R, M, arma::fill::zeros);
  arma::vec out(R, arma::fill::zeros);
  
  for(m = 0; m < M; m++){
    for(r = 0; r < R; r++){
      for(d = 0; d < D; d++){
        if(arma::approx_equal(unique_vectors_stats.row(r), stats.slice(m).row(d), "absdiff", 0.000001)){
          mat_counts(r,m) += 1;
        }
      }
    }
  }

  out = mat_counts * intereventTime; // we could include this step inside the loop above 
  
  return out;
}

//' computeOccurrencies
//'
//' A function to compute how many times each of the unique vector of statistics returned by getUniqueVectors() occurred in the network (as in contributing to the hazard in the likelihood). This function is based on the result shown by the Appendix C in the paper 'Hierarchical models for relational event sequences', DuBois et al. 2013 (pp. 308-309).
//'
//' @param edgelist is the preprocessed edgelist dataframe with information about [time,sender,receiver,type,weight] by row.
//' @param risksetMatrix matrix object inside the output list of the preprocessed relational event history.
//' @param M number of observed relational events.
//' @param unique_vectors_stats matrix of unique vectors of statistics (output of getUniqueVectors()).
//' @param stats array of statistics with dimensons [D*U*M]
//'
//' @return vector of q's
//'
//' @export
// [[Rcpp::export]]
arma::vec computeOccurrencies(const Rcpp::DataFrame& edgelist, const arma::umat& risksetMatrix, const arma::uword& M, const arma::mat& unique_vectors_stats, const arma::cube& stats){

  arma::uword r,m,d;
  arma::uword R = unique_vectors_stats.n_rows;
  arma::uword U = unique_vectors_stats.n_cols;
  Rcpp::IntegerVector sender = edgelist["sender"];
  Rcpp::IntegerVector receiver = edgelist["receiver"];
  arma::rowvec stats_event_m(U,arma::fill::zeros);
  arma::vec out(R);

  for(m = 0; m < M; m++){
    d = risksetMatrix(sender(m),receiver(m));
    stats_event_m = stats.slice(m).row(d);
    for(r = 0; r < R; r++){
      if(arma::approx_equal(stats_event_m, unique_vectors_stats.row(r), "absdiff", 0.000001)){
        out(r) += 1;
      }
    }
  }

  return out;
}


//' remDerivativesFast (a function that returns a list of 0th/1st/2nd order derivatives of loglikelihood evaluated in pars)
//'
//' description of the function here
//'
//' @param pars vector of parameters 
//' @param times_r  former m
//' @param occurrencies_r former q
//' @param unique_vectors_stats former U
//'
//' @return list of value/gradient/hessian in pars
//'
//' @export
// [[Rcpp::export]]
Rcpp::List remDerivativesFast(const arma::vec& pars, const arma::vec& times_r, const arma::vec& occurrencies_r, const arma::mat& unique_vectors_stats){
  
  arma::uword U = pars.n_elem;
  arma::uword u,k,l;

  double loglik = sum(occurrencies_r.t()*U*pars - times_r.t() * exp(unique_vectors_stats*pars));
  
  arma::vec grad(U);
  
  for(u = 0; u < U; u++){
    grad.row(u) =  sum(occurrencies_r % unique_vectors_stats.col(u) - times_r % exp(unique_vectors_stats * pars) % unique_vectors_stats.col(u));
  }
  
  arma::mat hess(U, U);
    
  for(k = 0; k < U; k++){
    for(l = 0; l < U; l++){
      hess(k,l) = sum(- times_r % exp(unique_vectors_stats * pars) % unique_vectors_stats.col(l) % unique_vectors_stats.col(k));
    }
  }
  
  return (Rcpp::List::create(Rcpp::Named("value") = -loglik, 
                             Rcpp::Named("gradient") = -grad,
                             Rcpp::Named("hessian") = -hess));
}

// /////////////////////////////////////////////////////////////////////////////////
// /////////////(END)     remstimateFAST routine functions      (END)/////////////// 
// /////////////////////////////////////////////////////////////////////////////////