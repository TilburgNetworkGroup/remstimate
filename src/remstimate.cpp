#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <omp.h>
#include <RcppArmadilloExtensions/sample.h>
#include <typeinfo>
#include <map>
#include <iterator>
#include <string>
#include "reh.h"
#include "remstimateBoost.h"

// loglikelihood function

// gradient function

// hessian function

// function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values

//////////////////////////////////////////////////////////////////////////////////
////////////             old functions are BELOW             /////////////////////
//////////////////////////////////////////////////////////////////////////////////

//' lpd (Log-Pointwise Density of REM)
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

//' nllik (Negative Log-Likelihood of REM)
//'
//' @param pars is a vector of parameters (note: the order must be the same as the column order in 'stats') at which to calculate the likelihood value
//' @param stats is a cube of dimensions n_dyads*variables*M with statistics of interest by column and dyads by row.
//' @param event_binary is a matrix of ones and zeros of dimensions M*n_dyads : 1 indicating the observed dyad and 0 the non-observed dyads.
//' @param interevent_time the vector of time differences between the current time point and the previous event time (note: interevent time at t_1 is t_1-0=t_1, since t_0 = 0).
//' @param threads
//'
//' @return negative log likelihood value
//'
//' @export
// [[Rcpp::export]]  
double nllik(arma::vec pars, 
            arma::cube stats, 
            arma::umat event_binary, 
            arma::vec interevent_time,
            int threads){

        arma::uword n_dyads = event_binary.n_cols;
        arma::uword i,m;
        arma::uword M = event_binary.n_rows;
        arma::vec log_lambda(n_dyads,arma::fill::zeros) ;
        arma::vec llik(M,arma::fill::zeros);

        omp_set_dynamic(0);           // disabling dynamic teams 
        omp_set_num_threads(threads); // number of threads for all consecutive parallel regions
        #pragma omp parallel for private(m,i,log_lambda) shared(n_dyads,M,stats,event_binary,interevent_time,llik)
        for(m = 0; m < M; m++)
        {
            log_lambda = stats.slice(m) * pars;
            for(i = 0; i < n_dyads; i++){
                if(event_binary(m,i) == 0){
                    llik(m) -= exp(log_lambda(i))*interevent_time(m);
                }
                else{
                    llik(m) += log_lambda(i)-exp(log_lambda(i))*interevent_time(m);
                }
            }
        }
        return -sum(llik);
    }
