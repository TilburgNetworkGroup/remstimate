#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <RcppArmadilloExtensions/sample.h>
#include <typeinfo>
#include <map>
#include <iterator>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef REMSTIMATEBOOST_H
#define REMSTIMATEBOOST_H


// things to do :
// (1) change naming convention of functions and arguments and make it consistent with the whole package (as well as the other packages in remverse);
// (2) some functions will already be defined in other files in the pckg;
// (3) some arguments can be changed with others that would require less memory to process.
    
// further developement to discuss on 12/06 :
// (1) what about the calculation of the percentage of potential improvement?
// (2) gradient and hessian according to the optimization procedure.

//' cube_to_matrix (a function to transform a cube to a matrix)
//'
//' description of the function here
//'
//' @param S cube structure to rearrange in a matrix of dimension [N(N-1)M*statistics]
//'
//' @return matrix
//'
//' @export
// [[Rcpp::export]]  
arma::mat cube_to_matrix(arma::cube S){
  arma::mat aux(S.n_rows*S.n_slices, S.n_cols, arma::fill::zeros); 
  arma::uword num = 0;
  
  for(arma::uword i = 0; i < S.n_slices; ++i){
    for(arma::uword j = 0; j < S.n_rows; ++j){
      aux.row(num) = S.slice(i).row(j);
      num += 1;
    }
  }
  
  return aux;
}


//' get_unique_vectors (a function to )
//'
//' description of the function here
//'
//' @param A  matrix of statistics per dyad and per time point
//'
//' @return matrix with only unique vectors
//'
//' @export
// [[Rcpp::export]]
arma::mat get_unique_vectors(const arma::mat& A){
  arma::uvec indices(A.n_rows, arma::fill::zeros);
  
  for(arma::uword i = 0; i < A.n_rows; ++i){
    for(arma::uword j = i+1; j < A.n_rows; ++j){
      if(arma::approx_equal(A.row(j), A.row(i), "absdiff", 0.000001)){
         indices(j) = 1; 
         break;
      }
    } 
  }  
  return A.rows(arma::find(indices == 0));
}

//' get_events_index (to remove, since a similar function will be written in reh.h)
//'
//' description of the function here
//'
//' @param edgelist matrix of time sender and receiver
//' @param riskset matrix of possible dyads
//'
//' @return vector of indices (indicating the dyad occurred at each time point)
//'
//' @export
// [[Rcpp::export]]
arma::uvec get_events_index(const arma::mat& edgelist, const arma::mat& riskset){
  arma::uvec index(edgelist.n_rows);
  
  for(arma::uword i = 0; i < edgelist.n_rows; ++i){
    
    arma::uword s0 = edgelist(i,1);
    arma::uword r0 = edgelist(i,2);
    
    for(arma::uword j = 0; j < riskset.n_rows; ++j){
      if((riskset(j,0) == s0) & (riskset(j,1) == r0)){
        index(i) = j;
      }
    }
  }
  return index;
}

//' compute_q (one of the most important functions to compute quantity q, see paper)
//'
//' description of the function here
//'
//' @param index get_events_index output
//' @param edgelist (this argument is only used to define a dimension, which is the number of events, then it can be omitted)
//' @param U matrix of unique vectors
//' @param S array of statistics with dimensons [dyads*statistics*time]
//'
//' @return vector of q's
//'
//' @export
// [[Rcpp::export]]
arma::vec compute_q(const arma::vec& index, const arma::mat& edgelist, const arma::mat& U, const arma::cube& S){
  arma::mat sReal(edgelist.n_rows, U.n_cols);
  arma::mat Uindex(U.n_rows, sReal.n_rows, arma::fill::zeros);
  
  for(arma::uword i = 0; i < index.n_elem; ++i){
    sReal.row(i) = S.slice(i).row(index(i));
  }
  
  for(arma::uword i = 0; i < sReal.n_rows; ++i){
    for(arma::uword j = 0; j < U.n_rows; ++j){
      if(arma::approx_equal(sReal.row(i), U.row(j), "absdiff", 0.000001)){
        Uindex(j,i) = 1;
      }
    }
  }
  
  arma::vec q(U.n_rows);
  
  for(arma::uword i = 0; i < q.n_elem; ++i){
    q(i) = arma::sum(Uindex.row(i));
  }
  
  return q;
}

//' compute_m (one of the most important functions to compute quantity m, see paper)
//'
//' description of the function here
//'
//' @param index get_events_index output
//' @param edgelist (this argument is only used to define a dimension, which is the number of events, then it can be omitted)
//' @param U matrix of unique vectors
//' @param S array of statistics with dimensons [dyads*statistics*time]
//'
//' @return vector of m's
//'
//' @export
// [[Rcpp::export]]
arma::vec compute_m(const arma::vec& index, const arma::mat& edgelist, const arma::mat& U, const arma::cube& S){
  arma::mat MM(U.n_rows, edgelist.n_rows, arma::fill::zeros);
  
  for(arma::uword i = 0; i < edgelist.n_rows; ++i){
    for(arma::uword j = 0; j < U.n_rows; ++j){
      for(arma::uword l = 0; l < S.n_rows; ++l){
        
        if(l == index(i)){continue;}
        
        if(arma::approx_equal(U.row(j), S.slice(i).row(l), "absdiff", 0.000001)){
          MM(j,i) += 1;
        }
      }
    }
  }
  
  arma::vec time(edgelist.n_rows+1, arma::fill::zeros);
  
  for(arma::uword i = 0; i < edgelist.n_rows; ++i) time(i+1) = edgelist(i,0);
  
  time = arma::diff(time);
  
  arma::vec m(U.n_rows);
  
  m = MM * time;
  
  return m;
}

//' logLike (computing loglikelihood value)
//'
//' description of the function here
//'
//' @param beta vector of parameters' value where to compute the value of the loglikelihood
//' @param U matrix of unique vectors
//' @param q vector computed with compute_q
//' @param m vector computed with compute_m
//'
//' @return vector of q's
//'
//' @export
// [[Rcpp::export]]
double logLike(arma::vec beta, const arma::mat& U, const arma::vec q, const arma::vec m){
  return(arma::sum(q.t() * (U * beta) - m.t() * exp(U * beta)));
}

#endif