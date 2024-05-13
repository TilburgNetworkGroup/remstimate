#include <RcppArmadillo.h>
#include <iostream>
#include <RcppArmadilloExtensions/sample.h> // used for the sample function inside RcppArmadillo
#include <iterator>
#include <numeric>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <remify/remify.h> 
#ifdef _OPENMP
#include <omp.h>
#endif

// /////////////////////////////////////////////////////////////////////////////////
// //////////(BEGIN)              remDerivatives                  (BEGIN)///////////
// /////////////////////////////////////////////////////////////////////////////////  

// remDerivativesStandard (tie-oriented modeling)
//
// function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values (for tie-oriented modeling).
//
// @param pars is a vector of parameters (note: the order must be aligned with the column order of the array 'stats').
// @param stats is a cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
// @param dyad is a list of M vectors, each one containing the ID of observed dyads from the 'remify' object.
// @param omit_dyad is a list of two objects: a vector "time" and a matrix "riskset". Two objects for handling changing risksets. The object is NULL if no change of the risk set structure is defined.
// @param interevent_time the time difference between the current time point and the previous event time.
// @param ordinal boolean that indicate whether to use the ordinal or interval timing likelihood (default is false).
// @param ncores integer referring to the number of threads for the parallelization (default is 1).
// @param gradient boolean true/false whether to return gradient value (default is true).
// @param hessian boolean true/false whether to return hessian value (default is true).
//
// @return list of values: loglik, gradient and hessian.
Rcpp::List remDerivativesStandard(const arma::vec &pars,
                                          const arma::cube &stats,
                                          const arma::field<arma::uvec> &dyad,
                                          const Rcpp::List &omit_dyad,
                                          const arma::vec &interevent_time,
                                          bool ordinal = false,
                                          int ncores = 1,
                                          bool gradient = true,
                                          bool hessian = true
                                          ){
    arma::uword U = stats.n_cols; // number of parameters 
    arma::uword D = stats.n_rows;
    arma::uword M = stats.n_slices; // number of events

    arma::uword d,m,l,k;
    //arma::vec log_lambda(D,arma::fill::zeros);

    arma::vec loglik(M,arma::fill::zeros);
    arma::cube hess(U,U,M,arma::fill::zeros);
    arma::mat grad(U,M,arma::fill::zeros);

    // initializing two empty objects for the two omit_dyad list objects
    // they will be filled only if omit_dyad has length >0 (i.e. it is not an empy list)
    arma::vec riskset_time_vec(M); 
    arma::mat riskset_mat;
    if(omit_dyad.size()>0){
      riskset_time_vec = Rcpp::as<arma::vec>(omit_dyad["time"]);
      riskset_mat = Rcpp::as<arma::mat>(omit_dyad["riskset"]);
    }
    else{
      riskset_time_vec.fill(-1); // to simplify the ifelse in the loop below
    }
    if(ordinal & hessian) gradient = true; // because in the ordinal likelihood we want to avoid to double compute certain quantities

    if(!ordinal){ // interval likelihood
      #ifdef _OPENMP
      omp_set_dynamic(0);         
      omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
      #pragma omp parallel for if(ncores>1) private(m,d,l,k) shared(M,D,U,loglik,hess,grad,pars,stats,riskset_time_vec,riskset_mat,dyad,gradient,hessian,interevent_time)
      #endif
      for(m = 0; m < M; m++){
          arma::mat stats_m = stats.slice(m).t(); // dimensions : [U*D] we want to access dyads by column
          arma::vec log_lambda = stats_m.t() * pars;
          int riskset_time_m = riskset_time_vec(m);
          double loglik_m = 0.0;
          arma::vec grad_m(U,arma::fill::zeros);
          arma::mat hess_m(U,U,arma::fill::zeros);
          arma::uvec events_occurred = dyad(m)-1; //  -1 because dyads' IDs must range between 0 and D-1
          // for the occured events we comput loglik and gradient contributes
          loglik_m += arma::accu(log_lambda(events_occurred));
          grad_m += arma::sum(stats_m.cols(events_occurred),1);
          // dealing with the risk set
          if(riskset_time_m!=(-1)){ // if the 'riskset_time_m' is different from (-1), then a dynamic riskset is observed
            for(d = 0; d < D; d++){
              if(riskset_mat(riskset_time_m,d)){ // ignoring impossible events that are not in risk set 
                  double ratewt = exp(log_lambda(d))*interevent_time(m);  // `ratewt` means (rate * waiting time)
                  loglik_m -= ratewt;
                  if(gradient){
                      grad_m -= stats_m.col(d)*ratewt;
                  }
                  if(hessian){
                    for(k = 0; k < U; k++){
                      for (l = k; l < U; l++){
                        hess_m(k,l) -= stats_m(l,d)*stats_m(k,d)*ratewt;
                        hess_m(l,k) = hess_m(k,l);
                      }
                    }
                  }

              }
            }
          }
          else{ // loop over all dyad (because for all the time points the riskset is fixed)
            for(d = 0; d < D; d++){
              double ratewt = exp(log_lambda(d))*interevent_time(m);  // `ratewt` means (rate * waiting time)
              loglik_m -= ratewt;
              if(gradient){
                grad_m -= stats_m.col(d)*ratewt;
              }
              if(hessian){
                for(k = 0; k < U; k++){
                  for (l = k; l < U; l++){
                      hess_m(k,l) -= stats_m(l,d)*stats_m(k,d)*ratewt;
                      hess_m(l,k) = hess_m(k,l);
                  }
                }
              }
            }
          }
          loglik(m) = loglik_m;
          grad.col(m) = grad_m;
          hess.slice(m) = hess_m;
      }
    }
    else{ // ordinal likelihood
      #ifdef _OPENMP
      omp_set_dynamic(0);         
      omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
      #pragma omp parallel for if(ncores>1) private(m,d,l,k) shared(M,D,U,loglik,hess,pars,grad,stats,riskset_time_vec,riskset_mat,dyad,gradient,hessian)
      #endif
      for(m = 0; m < M; m++){
        arma::mat stats_m = stats.slice(m).t();
        arma::vec log_lambda = stats_m.t() * pars;
        int riskset_time_m = riskset_time_vec(m);
        double surv = 0.0;
        double loglik_m = 0.0;
        arma::vec grad_m(U,arma::fill::zeros);
        arma::mat hess_m(U,U,arma::fill::zeros);
        arma::uvec events_occurred = dyad(m)-1;  // -1 because dyads' IDs must range between 0 and D-1
        // for the occured events we comput loglik and gradient contributes
        loglik_m += arma::accu(log_lambda(events_occurred));
        if(gradient){
          grad_m += arma::sum(stats_m.cols(events_occurred),1); // first component of the gradient
        }
        if(riskset_time_m!=(-1)){ // if the riskset_time_m is different from (-1), then a dynamic riskset is observed
          arma::vec lambda_only_dyads_at_risk = riskset_mat.row(riskset_time_m).t() % exp(log_lambda); // lambda is zero for events not at risk
          surv = arma::accu(lambda_only_dyads_at_risk); 
          if(gradient){
            grad_m *= surv; // gradient multiplied by surv here
            arma::vec lambda_d_stats_m = stats_m * lambda_only_dyads_at_risk;
            grad_m -= lambda_d_stats_m;
            grad_m /= surv; // gradient divided by surv here, so (first_component*surv-second_component)/surv
          }
         for(d = 0; d < D; d++){
            if(hessian){
              for(k = 0; k < U; k++){
                for (l = k; l < U; l++){
                  hess_m(k,l) -= stats_m.at(l,d)*stats_m.at(k,d)*lambda_only_dyads_at_risk.at(d);
                  hess_m(l,k) = hess_m(k,l);
                }
              }
            }   
          }
        }
        else{ // loop over all the dyads
          arma::vec lambda_only_dyads_at_risk = exp(log_lambda);
          surv = arma::accu(lambda_only_dyads_at_risk); 
          if(gradient){
            grad_m *= surv; // gradient multiplied by surv here
            arma::vec lambda_d_stats_m = stats_m * lambda_only_dyads_at_risk;
            grad_m -= lambda_d_stats_m;
            grad_m /= surv; // gradient divided by surv here, so (first_component*surv-second_component)/surv
          }
         for(d = 0; d < D; d++){
            if(hessian){
              for(k = 0; k < U; k++){
                for (l = k; l < U; l++){
                  hess_m(k,l) -= stats_m.at(l,d)*stats_m.at(k,d)*lambda_only_dyads_at_risk.at(d);
                  hess_m(l,k) = hess_m(k,l);
                }
              }
            }   
          }
        }
        loglik_m -= log(surv);
        if(hessian){
          for(k = 0; k < U; k++){
            for (l = k; l < U; l++){
              hess_m(k,l) /= surv;
              hess_m(k,l) += (grad_m(k) * grad_m(l))/(std::pow(surv,2));
              hess_m(l,k) = hess_m(k,l);
            }
          }
        }
        loglik(m) = loglik_m;
        grad.col(m) = grad_m;
        hess.slice(m) = hess_m;
      }
    }

    if(gradient && !hessian){
        return Rcpp::List::create(Rcpp::Named("value") = -sum(loglik), Rcpp::Named("gradient") = -sum(grad,1));
      }else if(!gradient && !hessian){
        return Rcpp::List::create(Rcpp::Named("value") = -sum(loglik));
      }else{
        arma::cube H = -sum(hess,2);
        return Rcpp::List::create(Rcpp::Named("value") = -sum(loglik), Rcpp::Named("gradient") = -sum(grad,1), Rcpp::Named("hessian") = H.slice(0));
      }
}



// remDerivativesSenderRates (actor-oriented modeling - sender model)
//
// function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values for estimating the sender rate parameters for the actor oriented model.
//
// @param pars is a vector of parameters (note: the order must be aligned with the column order of the array 'stats$sender_stats').
// @param stats is a cube of M slices. Each slice is a matrix of dimensions N*U with statistics of interest by column and actors by row.
// @param actor1 is a list actor1's observed at each time point (attr(reh,"actor1")-1).
// @param omit_dyad is a list of two objects: a vector "time" and a matrix "risksetSender". Two objects for handling changing risksets. The object is NULL if no change of the risk set structure is defined.
// @param interevent_time the time difference between the current time point and the previous event time.
// @param ordinal boolean that indicate whether to use the ordinal or interval timing likelihood (default is false).
// @param ncores number of cores for parallelization
// @param gradient boolean true/false whether to return gradient value (default is true).
// @param hessian boolean true/false whether to return hessian value (default is true).
//
// @return list of values: loglik, gradient and hessian.
Rcpp::List remDerivativesSenderRates(
        const arma::vec &pars,
        const arma::cube &stats,
        const arma::field<arma::uvec> &actor1,     
        const Rcpp::List &omit_dyad,
        const arma::vec &interevent_time,
        bool ordinal = false,
        int ncores = 1,
        bool gradient = true,
        bool hessian  = true){
  arma::uword U = stats.n_cols; // number of parameters 
  int M = stats.n_slices; // number of events
  int N = stats.n_rows; //Number of actors
  int n, m;

  //output
  arma::vec loglik(M,arma::fill::zeros);
  arma::cube fisher(U,U,M,arma::fill::zeros);
  arma::mat grad(U,M,arma::fill::zeros);

  // omit_dyad 
  arma::vec riskset_time_vec(M); 
  arma::mat riskset_mat;

  if(omit_dyad.size()>0){
    riskset_time_vec = Rcpp::as<arma::vec>(omit_dyad["time"]);
    riskset_mat = Rcpp::as<arma::mat>(omit_dyad["risksetSender"]);
  }
  else{
    riskset_time_vec.fill(-1); // to simplify the ifelse in the loop below
  }

  #ifdef _OPENMP
  omp_set_dynamic(0);         // disabling dynamic teams
  omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
  #pragma omp parallel for if(ncores>1) private(m,n) shared(N,M,stats,pars,actor1,interevent_time,riskset_time_vec,riskset_mat,loglik,grad,fisher,ordinal,hessian,gradient)
  #endif
  for(m = 0; m < M; m++){
    arma::mat stats_m = stats.slice(m).t(); 
    //this is exp(beta^T X) dot product
    arma::vec lambda_s = arma::exp(stats_m.t() * pars);
    int riskset_time_m = riskset_time_vec(m); // risk set at time m
    double loglik_m = 0.0;
    arma::vec grad_m(U,arma::fill::zeros);
    arma::vec expected_stat_m(U,arma::fill::zeros);
    arma::mat fisher_m(U,U,arma::fill::zeros);

    // sender(s) interacting at time t[m]
    arma::uvec sender = actor1(m)-1; // -1 because actors' IDs must range between 0 and N-1

    // loglik and gradient (first addend)
    loglik_m += arma::accu(arma::log(lambda_s(sender)));
    if(gradient){
      grad_m += arma::sum(stats_m.cols(sender),1); //grad
    }

    // changes in the riskset at m-th event
    if(riskset_time_m!=(-1)){
      double sum_lambda = sum(riskset_mat.row(riskset_time_m).t() % lambda_s); 
      if(ordinal){
        loglik_m -= std::log(sum_lambda); // loglik second addend
        if(gradient | hessian){
          for(n = 0; n < N; n++){         //loop throughout all actors
            if(riskset_mat(riskset_time_m,n) == 1){
              expected_stat_m += lambda_s(n) * (stats_m.col(n));
              expected_stat_m /= sum_lambda; //exp(params_s * X_i)*X_i / sum_h (exp(params * X_h))
              if(gradient){
                grad_m -= expected_stat_m; //grad
              }
              if(hessian){
                fisher_m += (expected_stat_m * expected_stat_m.t());
                fisher_m -= (lambda_s(n) / sum_lambda) *( stats_m.col(n) * (stats_m.col(n).t()));
                //fisher += fisher_m;
              }   
            }
          }
        }
      } 
      else{
        loglik_m -=  sum_lambda * interevent_time(m);
        if(gradient){
          grad_m -= interevent_time(m) * stats_m * (riskset_mat.row(riskset_time_m).t() % lambda_s); 
        }
        if(hessian){
          for(n = 0; n < N; n++){         //loop throughout all actors
            fisher_m += (riskset_mat(riskset_time_m,n) * lambda_s(n)) * (stats_m.col(n) * (stats_m.col(n).t()) );
          }
          //fisher -= interevent_time(m) * fisher_m;
        }
      }
    }
    else {  //no dynamic riskset
      double sum_lambda = sum(lambda_s);
      if(ordinal){
        loglik_m -= std::log(sum_lambda);
        if(gradient | hessian){
          for(n = 0; n < N; n++){         //loop throughout all actors
            expected_stat_m += lambda_s(n) * (stats_m.col(n));
            expected_stat_m /= sum_lambda; //exp(params_s * X_i)*X_i / sum_h (exp(params * X_h))
            if(gradient){
              grad_m -= expected_stat_m;
            }

            if(hessian){
              fisher_m += (expected_stat_m * expected_stat_m.t());
              fisher_m -= (lambda_s(n) / sum_lambda) * (stats_m.col(n) * (stats_m.col(n).t()));
              //fisher += fisher_m;
            }
          }
        }
      } 
      else{
        loglik_m -=  sum_lambda * interevent_time(m);
        if(gradient){
          grad_m -= interevent_time(m) * stats_m * lambda_s; //(6)
        }
        if(hessian){
          for(n = 0; n < N; n++){         //loop throughout all actors
            fisher_m += lambda_s(n)*(stats_m.col(n) * (stats_m.col(n).t()));
          }
          //fisher -= interevent_time(m) * fisher_m;
        }
      }
    }

    loglik(m) = loglik_m;
    if(gradient){ // assigning gradient value at time m
      grad.col(m) = grad_m;
    }
    if(hessian){ // assigning hessian value at time m (different value between ordinal and interval likelihood)
      if(ordinal){
        fisher.slice(m) = fisher_m;
      }
      else{
        fisher.slice(m) = (-1.0) * interevent_time(m) * fisher_m;
      }
    }
  }

  // processing output
  if(gradient && !hessian){
    return Rcpp::List::create(Rcpp::Named("value") = -sum(loglik), Rcpp::Named("gradient") = -sum(grad,1));
  }else if(!gradient && !hessian){
    return Rcpp::List::create(Rcpp::Named("value") = -sum(loglik));
  }else{
    arma::cube H = -sum(fisher,2);
    return Rcpp::List::create(Rcpp::Named("value") = -sum(loglik), Rcpp::Named("gradient") = -sum(grad,1), Rcpp::Named("hessian") = H.slice(0));
  }
}


// remDerivativesReceiverChoice (actor-oriented modeling - receiver model)
//
// function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values for estimating the receiver choice parameters for the actor oriented model.
//
// @param pars is a vector of parameters (note: the order must be aligned with the column order of the array 'stats$receiver_stats').
// @param stats is cube of M slices. Each slice is a matrix of dimensions N*U with statistics of interest by column and actors by row.
// @param actor1 list of actor1's observed at each time point (attr(reh,"actor1")-1).
// @param actor2 list of actor2's observed at each time point (attr(reh,"actor2")-1).
// @param omit_dyad is a list of two objects: a vector "time" and a matrix "riskset". Two objects for handling changing risksets. The object is NULL if no change of the risk set structure is defined.
// @param interevent_time the time difference between the current time point and the previous event time.
// @param ncores number of cores for parallelization
// @param N the number of actors in the network (reh$N).
// @param gradient boolean true/false whether to return gradient value (default is true).
// @param hessian boolean true/false whether to return hessian value (default is true).
//
// @return list of values: loglik, gradient and hessian.
Rcpp::List remDerivativesReceiverChoice(
        const arma::vec &pars,
        const arma::cube &stats,
        const arma::field<arma::uvec> &actor1,
        const arma::field<arma::uvec> &actor2,
        const Rcpp::List &omit_dyad,
        const arma::vec &interevent_time,
        int N,
        int ncores = 1,
        bool gradient = true,
        bool hessian  = true){

    int U = stats.n_cols; // number of parameters for dyad stats
    int M = stats.n_slices; // number of events
    int m;
    int n;

    //output
    arma::vec loglik(M,arma::fill::zeros);
    arma::cube fisher(U,U,M,arma::fill::zeros);
    arma::mat grad(U,M,arma::fill::zeros);

    // omit_dyad 
    arma::vec riskset_time_vec(M); 
    arma::mat riskset_mat;
    if(omit_dyad.size()>0){
      riskset_time_vec = Rcpp::as<arma::vec>(omit_dyad["time"]);
      riskset_mat = Rcpp::as<arma::mat>(omit_dyad["riskset"]);
    }
    else{
      riskset_time_vec.fill(-1); // to simplify the ifelse in the loop below
    }

    #ifdef _OPENMP
    omp_set_dynamic(0);         // disabling dynamic teams
    omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
    #pragma omp parallel for if(ncores>1) private(m,n) shared(N,M,stats,pars,actor1,actor2,riskset_time_vec,riskset_mat,loglik,grad,fisher,hessian,gradient)
    #endif
    for(m = 0; m < M; m++){
      arma::mat stats_m = stats.slice(m).t(); // dimensions : [P*N] we want to access dyads by column

      //this is exp(beta^T X) dot product
      arma::vec lambda_d = arma::exp(stats_m.t() * pars);

      //actors
      arma::uvec sender = actor1(m)-1; // -1 because actors' IDs must range between 0 and N-1, we will call sender[0] because actor1 is an unlisted list of vectors (so at each m there is only one sender)
      arma::uvec receiver = actor2(m)-1; // -1 because actors' IDs must range between 0 and N-1, we will call receiver[0] because actor2 is an unlisted list of vectors (so at each m there is only one receiver)
      int sender_0_int = static_cast<int>(sender(0));
      // riskset_m
      int riskset_time_m = riskset_time_vec(m);

      double loglik_m = 0.0;
      arma::vec expected_stat_m(U,arma::fill::zeros);
      arma::mat fisher_m(U,U,arma::fill::zeros);

      double denom = 0.0;
      
      // changes in the riskset at m-th event
      if(riskset_time_m!=(-1)){ 
          for(n = 0; n<N; n++){
              arma::uword dyad = remify::getDyadIndex(sender(0),n,0,N,true);
              if((n!=sender_0_int) && (riskset_mat(riskset_time_m,dyad) == 1)){ // dynamic riskset
                  //loglik
                  denom += lambda_d(n); // exp(param_d * X_sender_i) 
                  if(hessian){
                  fisher_m -= lambda_d(n)*(stats_m.col(n) * stats_m.col(n).t());
                  }
                  if(gradient | hessian){
                  expected_stat_m += lambda_d(n) * stats_m.col(n);
                  }
              }
          }
      }
      else{ // no changes in the riskset at m-th event
          for(n = 0;n<N; n++){ //loop throught all actors
          if(n!=sender_0_int){ // 
              // loglik
              denom += lambda_d(n); // exp(param_d * X_sender_i) )
              if(hessian){
                  fisher_m -= lambda_d(n)*(stats_m.col(n) * stats_m.col(n).t());
              }
              if(gradient | hessian){
              expected_stat_m += lambda_d(n) * stats_m.col(n);
              }
          }
          }
      }

      loglik_m +=  std::log(lambda_d(receiver(0))); // params_d^T * X_sender_receiver
      loglik_m -= std::log(denom);
      // assigning likelihood value
      loglik(m) = loglik_m;
      if(gradient | hessian){
          expected_stat_m /= denom; 
      }
      // assigning gradient
      if(gradient){
          grad.col(m) = stats_m.col(receiver(0)) - expected_stat_m; 
      }
      if(hessian){
          fisher_m /= denom;
          fisher_m += (expected_stat_m * expected_stat_m.t());
          fisher.slice(m) = fisher_m;
      } 
    }

  // processing output
  if(gradient && !hessian){
    return Rcpp::List::create(Rcpp::Named("value") = -sum(loglik), Rcpp::Named("gradient") = -sum(grad,1));
  }else if(!gradient && !hessian){
    return Rcpp::List::create(Rcpp::Named("value") = -sum(loglik));
  }else{
    arma::cube H = -sum(fisher,2);
    return Rcpp::List::create(Rcpp::Named("value") = -sum(loglik), Rcpp::Named("gradient") = -sum(grad,1), Rcpp::Named("hessian") = H.slice(0));
  }
}



// remDerivatives (for actor-oriented and tie-oriented modeling)
//
// function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values.
//
// @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats').
// @param stats is cube of M slices. Each slice is a matrix of dimensions D*U (or N*U for the actor-oriented model) with statistics of interest by column and dyads (actors, for the actor-oriented model) by row.
// @param actor1 list of actor1's observed per each time point (attr(reh,"actor1")-1).
// @param actor2 list of actor2's observed per each time point (attr(reh,"actor2")-1).
// @param dyad list of dyads observed per each time point (from the attribute attr(remify,"dyad")-1).
// @param omit_dyad is a list of two objects: a vector "time" and a matrix "riskset" (or "risksetSender" for the sender model). Two objects for handling changing risksets. The object is NULL if no change of the risk set structure is defined.
// @param interevent_time the time difference between the current time point and the previous event time.
// @param model either "actor" or "tie" oriented model.
// @param ordinal boolean that indicates whether to use the ordinal or interval timing likelihood (default is false).
// @param ncores number of threads to use for the parallelization (default is 1).
// @param gradient boolean true/false whether to return gradient value (default is true).
// @param hessian boolean true/false whether to return hessian value (default is true).
// @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false). Default is true.
// @param N number of actors in the network. This argument is used only in the ReceiverChoice likelihood (model = "actor"). Default is NULL.
//
// @return list of values: loglik, gradient, hessian.
//
// [[Rcpp::export]]
Rcpp::List remDerivatives(const arma::vec &pars,
                                  const arma::cube &stats,
                                  const arma::field<arma::uvec> &actor1,
                                  const arma::field<arma::uvec> &actor2,
                                  const arma::field<arma::uvec> &dyad,
                                  const Rcpp::List &omit_dyad,
                                  const arma::vec &interevent_time,
                                  std::string model,
                                  bool ordinal = false,
                                  int ncores = 1,
                                  bool gradient = true,
                                  bool hessian = true,
                                  bool senderRate = true,
                                  Rcpp::Nullable<int> N = R_NilValue){
  Rcpp::List out;
  std::vector<std::string> models = {"tie","actor"};
  std::vector<std::string>::iterator itr = std::find(models.begin(), models.end(), model);
  auto which_model = std::distance(models.begin(), itr);

  switch (which_model)
  {
    case 0: { out = remDerivativesStandard(pars,stats,dyad,omit_dyad,interevent_time,ordinal,ncores,gradient,hessian);
      break;}

    case 1: { 
        switch (senderRate){ // both likelihood miss paralellization
          case 0 : {
                    out = remDerivativesReceiverChoice(pars,stats,actor1,actor2,omit_dyad,interevent_time,Rcpp::as<int>(N),ncores,gradient,hessian);
                    break;
                    }
          case 1 : {
                    out = remDerivativesSenderRates(pars,stats,actor1,omit_dyad,interevent_time,ordinal,ncores,gradient,hessian);
                    break;
                    }
        }
      }
  }
  return out;
}

// /////////////////////////////////////////////////////////////////////////////////
// ///////////(END)                remDerivatives                 (END)/////////////
// /////////////////////////////////////////////////////////////////////////////////



// /////////////////////////////////////////////////////////////////////////////////
// ///////////(BEGIN)             Gradient Descent              (BEGIN)/////////////
// /////////////////////////////////////////////////////////////////////////////////

// GDADAMAX (for actor-oriented and tie-oriented modeling)
//
// function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values.
//
// @param pars vector of model parameters.
// @param stats array of statistics.
// @param actor1 list of actor1's observed per each time point (attr(reh,"actor1")-1).
// @param actor2 list of actor2's observed per each time point (attr(reh,"actor2")-1).
// @param dyad list of dyads observed per each time point (from the attribute attr(remify,"dyad")-1).
// @param omit_dyad is a list of two objects: vector "time" and matrix "riskset" (or "risksetSender" for the sender model). Two object for handling changing risksets. NULL if no change is defined.
// @param interevent_time vector of interevent times (inside the reh object).
// @param model either "actor" or "tie" oriented model.
// @param ordinal boolean that indicate whether to use the ordinal or interval timing likelihood (default is false)
// @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false).
// @param gradient boolean true/false whether to return gradient value.
// @param hessian boolean true/false whether to return hessian value.
// @param N number of actors in the network. This argument is used only in the ReceiverChoice likelihood (model = "actor"). Default is NULL.
// @param ncores number of threads to use for the parallelization (default is 1).
// @param epochs number of epochs (default is 1e03).
// @param learning_rate learning rate (default is 0.002).
// @param beta1 hyperparameter beta1 (default is 0.9).
// @param beta2 hyperparameter beta2 (default is 0.999).
// @param epsilon hyperparameter eta (default is 0.01).
//
// @return optimization with GDADAM
//
// [[Rcpp::export]]
Rcpp::List GDADAMAX(const arma::vec &pars,
                  const arma::cube &stats,
                  const arma::field<arma::uvec> &actor1,
                  const arma::field<arma::uvec> &actor2,
                  const arma::field<arma::uvec> &dyad,
                  const Rcpp::List &omit_dyad,
                  const arma::vec &interevent_time,
                  std::string model,
                  bool ordinal = false,
                  bool senderRate = true,
                  bool gradient = true,
                  bool hessian = false,
                  Rcpp::Nullable<int> N = R_NilValue,
                  int ncores = 1,
                  int epochs = 1e03,
                  double learning_rate = 0.002,
                  double beta1 = 0.9,
                  double beta2 = 0.999,
                  double epsilon = 0.01){

  arma::uword P = pars.n_elem;
  double loglik,loglik_prev;
  int i = 0;
  arma::vec moment(P,arma::fill::zeros);
  arma::vec inf_norm(P,arma::fill::zeros);
  arma::vec iterations_loglik(epochs,arma::fill::zeros);
  arma::vec step_size(1,arma::fill::zeros);
  arma::vec pars_loc = pars;
  arma::vec grad(P,arma::fill::zeros);
  std::string reason;
  arma::mat iterations_pars(P,epochs,arma::fill::zeros);


  while(i <= (epochs-1)){
    Rcpp::List derv = remDerivatives(pars_loc,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,gradient,hessian,senderRate,N);

    // calculate loglik and gradient with the previous iteration
    loglik = derv["value"];
    grad = Rcpp::as<arma::vec>(derv["gradient"]);

    //Rcpp::Rcout << std::abs(loglik_prev - loglik) << "\n";
    if((i > 1)  & (std::abs(loglik_prev - loglik) < epsilon)){
        reason = "epsilon condition reached";
        break;
    }

    loglik_prev = loglik;

    // store values
    iterations_pars.col(i) = pars_loc;
    iterations_loglik(i) = loglik;

    // updating moment vector
    arma::vec moment_prev = moment;
    moment = (beta1*moment_prev) + (1-beta1)*grad;

    // updating exponentially weighted infinity norm
    arma::vec inf_norm_prev = inf_norm;
    arma::vec abs_grad = arma::abs(grad);
    for(arma::uword p=0; p<P; p++){
      arma::vec inf_norm_both = {beta2*inf_norm_prev(p),abs_grad(p)};
      inf_norm(p) = arma::max(inf_norm_both);
    }

    // updating parameter value
    step_size(0) = learning_rate/(1-std::pow(beta1,i+1));

    // updating gradient
    arma::vec delta = moment / inf_norm;

    // updating parameter value
    pars_loc = pars_loc - (step_size(0) * delta);

    i += 1;
  }

  if(reason.size() == 0){
    reason = "max epochs condition reached";
  }

  return Rcpp::List::create(Rcpp::Named("value") = loglik, 
                            Rcpp::Named("argument") = pars_loc, 
                            Rcpp::Named("gradient") = grad,Rcpp::Named("iterations_loglik") = iterations_loglik(arma::span(0,i-1)), 
                            Rcpp::Named("iterations_pars") = iterations_pars(arma::span::all,arma::span(0,i-1)),
                            Rcpp::Named("iterations") = i, Rcpp::Named("converged") = reason);
}

// /////////////////////////////////////////////////////////////////////////////////
// ///////////(END)               Gradient Descent                (END)/////////////
// /////////////////////////////////////////////////////////////////////////////////



// /////////////////////////////////////////////////////////////////////////////////
// ///////////(START)          Hamiltonian Monte Carlo          (START)/////////////
// /////////////////////////////////////////////////////////////////////////////////


// logPostHMC (for actor-oriented and tie-oriented modeling)
//
// This function calculates the value of the log-posterior density given the loglikelihood and the log-prior density.
//
// @param meanPrior is a vector of prior means with the same dimension as the vector of parameters (pars).
// @param sigmaPrior is a variance and covariance matrix used as a prior.
// @param pars vector of model parameters.
// @param stats array of statistics.
// @param actor1 list of actor1's observed per each time point (attr(reh,"actor1")-1).
// @param actor2 list of actor2's observed per each time point (attr(reh,"actor2")-1).
// @param dyad list of dyads observed per each time point (from the attribute attr(remify,"dyad")-1).
// @param omit_dyad is a list of two objects: vector "time" and matrix "riskset" (or "risksetSender" for the sender model). Two object for handling changing risksets. NULL if no change is defined.
// @param interevent_time vector of interevent times (inside the reh object).
// @param model either "actor" or "tie" oriented model.
// @param ordinal boolean that indicate whether to use the ordinal or interval timing likelihood (default is false)
// @param ncores number of threads to use for the parallelization (default is 1).
// @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false).
// @param N number of actors in the network. This argument is used only in the ReceiverChoice likelihood (model = "actor"). Default is NULL.
//
// @return value of log-posterior density.
double logPostHMC(const arma::vec &meanPrior,
                  const arma::mat &sigmaPrior,
                  const arma::vec &pars,
                  const arma::cube &stats,
                  const arma::field<arma::uvec> &actor1,
                  const arma::field<arma::uvec> &actor2,
                  const arma::field<arma::uvec> &dyad,
                  const Rcpp::List &omit_dyad,
                  const arma::vec &interevent_time,
                  std::string model,
                  bool ordinal = false,
                  int ncores = 1,
                  bool senderRate = true,
                  Rcpp::Nullable<int> N = R_NilValue){

  Rcpp::List derv = remDerivatives(pars,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,false,false,senderRate,N);
  double derv_0 = Rcpp::as<double>(derv[0]);                
  double prior =  arma::accu(0.5 * (pars.t() - meanPrior.t()) * inv(sigmaPrior) * (pars - meanPrior));
  return (prior + derv_0); 
}



// logPostGradientHMC (for actor-oriented and tie-oriented modeling)
//
// This function calculates the value of the gradient of the log-posterior density.
//
// @param meanPrior is a vector of prior means with the same dimension as the vector of parameters (pars).
// @param sigmaPrior is a variance and covariance matrix used as a prior.
// @param pars vector of model parameters.
// @param stats array of statistics.
// @param actor1 list of actor1's observed per each time point (attr(reh,"actor1")-1).
// @param actor2 list of actor2's observed per each time point (attr(reh,"actor2")-1).
// @param dyad list of dyads observed per each time point (from the attribute attr(remify,"dyad")-1).
// @param omit_dyad is a list of two objects: vector "time" and matrix "riskset" (or "risksetSender" for the sender model). Two object for handling changing risksets. NULL if no change is defined.
// @param interevent_time vector of interevent times (inside the reh object).
// @param model either "actor" or "tie" oriented model.
// @param ordinal boolean that indicate whether to use the ordinal or interval timing likelihood (default is false)
// @param ncores number of threads to use for the parallelization (default is 1).
// @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false).
// @param N number of actors in the network. This argument is used only in the ReceiverChoice likelihood (model = "actor"). Default is NULL.
//
// @return value of log-posterior gradient.
arma::vec logPostGradientHMC(const arma::vec &meanPrior,
                              const arma::mat &sigmaPrior,
                              const arma::vec &pars,
                              const arma::cube &stats,
                              const arma::field<arma::uvec> &actor1,
                              const arma::field<arma::uvec> &actor2,
                              const arma::field<arma::uvec> &dyad,
                              const Rcpp::List &omit_dyad,
                              const arma::vec &interevent_time,
                              std::string model,
                              bool ordinal = false,
                              int ncores = 1,
                              bool senderRate = true,
                              Rcpp::Nullable<int> N = R_NilValue){
  Rcpp::List derv = remDerivatives(pars,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,true,false,senderRate,N);
  arma::vec gprior = inv(sigmaPrior)*(pars - meanPrior); // the sign is already changed here
  arma::vec glp = Rcpp::as<arma::vec>(derv[1]);
  return (glp + gprior); 
}

                                  
// iterHMC (for actor-oriented and tie-oriented modeling)
//
// This function does one iteration of the Hamiltonian Monte carlo.
//
// @param L number of leapfrogs. Default (and recommended) value is 100.
// @param epsilon size of the leapfrog (default value is 1e-02).
// @param meanPrior is a vector of prior means with the same dimension as the vector of parameters (pars).
// @param sigmaPrior is a variance and covariance matrix used as a prior.
// @param pars vector of model parameters.
// @param stats array of statistics.
// @param actor1 list of actor1's observed per each time point (attr(reh,"actor1")-1).
// @param actor2 list of actor2's observed per each time point (attr(reh,"actor2")-1).
// @param dyad list of dyads observed per each time point (from the attribute attr(remify,"dyad")-1).
// @param omit_dyad is a list of two objects: vector "time" and matrix "riskset" (or "risksetSender" for the sender model). Two object for handling changing risksets. NULL if no change is defined.
// @param interevent_time vector of interevent times (inside the reh object).
// @param model either "actor" or "tie" oriented model.
// @param ordinal boolean that indicate whether to use the ordinal or interval timing likelihood (default is false)
// @param ncores number of threads to use for the parallelization (default is 1).
// @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false).
// @param N number of actors in the network. This argument is used only in the ReceiverChoice likelihood (model = "actor"). Default is NULL.
//
// @return output of one iteration of the Hamiltonian Monte Carlo
arma::field<arma::vec> iterHMC(arma::uword L,
                  double epsilon,
                  const arma::vec &meanPrior,
                  const arma::mat &sigmaPrior,
                  const arma::vec &pars,
                  const arma::cube &stats,
                  const arma::field<arma::uvec> &actor1,
                  const arma::field<arma::uvec> &actor2,
                  const arma::field<arma::uvec> &dyad,
                  const Rcpp::List &omit_dyad,
                  const arma::vec &interevent_time,
                  std::string model,
                  bool ordinal = false,
                  int ncores = 1,
                  bool senderRate = true,
                  Rcpp::Nullable<int> N = R_NilValue){

  arma::uword P = pars.size(); //number of parameters
  arma::field<arma::vec> out(2); // output object

  // proposing a new value for the momentum variable
  Rcpp::NumericVector draw = Rcpp::rnorm(P);
  arma::vec rP = Rcpp::as<arma::vec>(draw);  
  // current momentum variable
  arma::vec rC = rP;
  
  arma::vec betaC = pars; // current beta
  arma::vec betaP = pars; // new proposed beta (starting value is current beta)

  //leapfrog algorithm, updates via Hamiltonian equations
  rP -= (0.5 * epsilon * logPostGradientHMC(meanPrior,sigmaPrior,betaP,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N));
  for(arma::uword i = 1; i <= L; i++){
    betaP += (epsilon * rP);
    if(i != L) rP -= (epsilon * logPostGradientHMC(meanPrior,sigmaPrior,betaP,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N));
  }
  rP -= (0.5 * epsilon * logPostGradientHMC(meanPrior,sigmaPrior,betaP,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N));
  rP *= (-1);

  //computes final quantities for the acceptance rate
  double U = logPostHMC(meanPrior,sigmaPrior,betaC,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N);
  double propU = logPostHMC(meanPrior,sigmaPrior,betaP,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N);
  double K = 0.5 * arma::accu(arma::pow(rC,2));
  double propK = 0.5 * arma::accu(arma::pow(rP,2));

  //Accepting or rejecting the proposal
  Rcpp::NumericVector randomUnif = Rcpp::runif(1);
  if(randomUnif[0] < exp(U - propU + propK - K)){   // before it was: exp( (-propU - propK)/(-U - K) ) .. exp(propH)/exp(H)
    out(0) = betaP;
    out(1) = propU;
  } else {
    out(0) = betaC;
    out(1) = U;
  }
  return out;
}


// burninHMC
//
// This function performs the burn-in and the thinning at the end of the HMC.
//
// @param samples cube with final draws from HMC algorithm.
// @param loglik matrix of values of the posterior loglikelihood at the different draws.
// @param burnin the number of draws to discard after running the chains.
// @param thin the number of draws to be skipped. For instance, if thin = 10, draws will be selected every 10 generated draws: 1, 11, 21, 31, ...
//
// @return list of two objects: draws and loglik after burnin and thinning step.
Rcpp::List burninHMC(const arma::cube& samples, const arma::mat& loglik, arma::uword burnin, arma::uword thin = 1){

  arma::uword rows = round((samples.n_rows - burnin)/thin); //number of rows of output
  arma::uword nchains = samples.n_slices; // same dimension as loglik.n_cols
  arma::mat out_draws(rows*nchains, samples.n_cols); // output draws, sample.n_cols is equal to the number of parameters
  arma::vec out_loglik(rows*nchains); // output loglik
  arma::uword i,j;
  Rcpp::List out = Rcpp::List();

  for(i = 0; i < nchains; i++){
    arma::uword num = burnin; 
    for(j = 0; j < rows; j++){
      out_draws.row(j+(rows*i)) = samples.slice(i).row(num);
      out_loglik(j+(rows*i)) = loglik(num,i);
      num += thin; // increment index by thinning step
    }
  }
  return Rcpp::List::create(out_draws,out_loglik);
}


// HMC (for actor-oriented and tie-oriented modeling)
//
// This function performs the Hamiltonian Monte Carlo on a tie-oriented or an actor-oriented modeling framework.
//
// @param pars_init is a matrix of dimensions U x nchains where for each column (chain) a random vector of initial values for the parameter is supplied.
// @param nsim is the number of samples from the posterior that have to be generated.
// @param nchains number of chains of length nsim.
// @param burnin is the number of draws to discard after running the chains.
// @param meanPrior is a vector of prior means with the same dimension as the vector of parameters (pars).
// @param sigmaPrior is a variance and covariance matrix used as a prior.
// @param stats array of statistics.
// @param actor1 list of actor1's observed per each time point (attr(reh,"actor1")-1).
// @param actor2 list of actor2's observed per each time point (attr(reh,"actor2")-1).
// @param dyad list of dyads observed per each time point (from the attribute attr(remify,"dyad")-1).
// @param omit_dyad is a list of two objects: vector "time" and matrix "riskset" (or "risksetSender" for the sender model). Two object for handling changing risksets. NULL if no change is defined.
// @param interevent_time vector of interevent times (inside the reh object).
// @param model either "actor" or "tie" oriented model.
// @param ordinal boolean that indicate whether to use the ordinal or interval timing likelihood (default is false).
// @param ncores number of threads to use for the parallelization (default is 1).
// @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false).
// @param N number of actors in the network. This argument is used only in the ReceiverChoice likelihood (model = "actor"). Default is NULL.
// @param thin is the number of draws to be skipped. For instance, if thin = 10, draws will be selected every 10 generated draws: 1, 11, 21, 31, ...
// @param L number of leapfrogs. Default (and recommended) value is 100.
// @param epsilon size of the leapfrog (default value is 1e-02).
//
// @return list of two objects: draws (for all the chains) and vector of values of the log-posterior density.
//
// [[Rcpp::export]]
Rcpp::List HMC(arma::mat pars_init,
                arma::uword nsim,
                arma::uword nchains,
                arma::uword burnin,
                const arma::vec& meanPrior,
                const arma::mat& sigmaPrior,
                const arma::cube &stats,
                const arma::field<arma::uvec> &actor1,
                const arma::field<arma::uvec> &actor2,
                const arma::field<arma::uvec> &dyad,
                const Rcpp::List &omit_dyad,
                const arma::vec &interevent_time,
                std::string model,
                bool ordinal = false,
                int ncores = 1,
                bool senderRate = true,
                Rcpp::Nullable<int> N = R_NilValue,
                arma::uword thin = 1,
                arma::uword L = 100,
                double epsilon = 0.01){ 
  arma::cube array_of_draws(nsim, pars_init.n_rows, nchains); // array of draws from the posterior distribution (chains are by slice)
  arma::mat matrix_of_loglik(nsim,nchains); // matrix of posterior loglikelihood values across chains
  arma::uword j,i;
  Rcpp::List out = Rcpp::List::create(); // output object

  // starting generating chains
  for(j = 0; j < nchains; j++){ //looping through chains
    arma::mat chain_j(pars_init.n_rows,nsim,arma::fill::zeros);
    arma::vec chain_j_loglik(nsim,arma::fill::zeros);
    //[i=0] this step only get the first sample out of the starting value
    arma::field<arma::vec> iter_i_hmc = iterHMC(L,epsilon,meanPrior,sigmaPrior,pars_init.col(j),stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N);
    chain_j.col(0) = iter_i_hmc[0]; // saving draws
    chain_j_loglik(0) = arma::conv_to<double>::from(iter_i_hmc[1]); // saving posterior loglikelihood
    //Rcpp::Rcout << "first value of the chain: " << chain_j.col(0) << "\n";
    for(i = 1; i < nsim; i++){ //looping through iterations of the MCMC
      //Rcpp::Rcout <<i <<"\n"; progress bar on chains progress
      //Then the next step will always be based on the previous one
      arma::field<arma::vec> iter_i_hmc = iterHMC(L,epsilon,meanPrior,sigmaPrior,chain_j.col(i-1),stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N);
      chain_j.col(i) = iter_i_hmc[0]; // saving draws
      chain_j_loglik(i) = arma::conv_to<double>::from(iter_i_hmc[1]); // saving posterior loglikelihood
    }
    array_of_draws.slice(j) = chain_j.t();
    matrix_of_loglik.col(j) = chain_j_loglik;
  }
  //this step performs the burn-in and thinning
  Rcpp::List processed_chains = burninHMC(array_of_draws,matrix_of_loglik,burnin,thin); // it would be ideal to make burninHMC return a matrix
  
  out["draws"] = processed_chains[0];
  out["log_posterior"] = processed_chains[1];

  return out; 
}

// /////////////////////////////////////////////////////////////////////////////////
// ///////////(END)            Hamiltonian Monte Carlo            (END)/////////////
// /////////////////////////////////////////////////////////////////////////////////


// /////////////////////////////////////////////////////////////////////////////////
// ///////////(BEGIN)                     WAIC                  (BEGIN)/////////////
// /////////////////////////////////////////////////////////////////////////////////



// getWAIC 
//
// a function for computing the Watanabe Akaike's Information Criterion after estimating the model with remstimate
//
// @param mu vector of estimates (from MLE or GDADAMAX method). This argument will be an empty vector for method BSIR and HMC.
// @param vcov matrix of covariances estimated from the data (via MLE or GDADAMAX method). This argument will be an empty matrix for method BSIR and HMC.
// @param pars a matrix of posterior draws (from BSIR or HMC method). This argument will be an empty matrix for method MLE and GDADAMAX
// @param stats is cube of M slices. Each slice is a matrix of dimensions D*U (or N*U for the actor-oriented model) with statistics of interest by column and dyads (actors, for the actor-oriented model) by row.
// @param actor1 list of actor1's observed per each time point (attr(reh,"actor1")-1).
// @param actor2 list of actor2's observed per each time point (attr(reh,"actor2")-1).
// @param dyad list of dyads observed per each time point (from the attribute attr(remify,"dyad")-1).
// @param interevent_time the time difference between the current time point and the previous event time.
// @param omit_dyad is a list of two objects: a vector "time" and a matrix "riskset" (or "risksetSender" for the sender model). Two objects for handling changing risksets. The object is NULL if no change of the risk set structure is defined.
// @param model either "actor" or "tie" oriented model.
// @param approach either "Bayesian" (then expecting 'pars' to be a matrix filled with draws from the posterior distribution), "Frequentist" (then expecting 'mu' and 'vcov' to be supplied).
// @param ordinal boolean that indicates whether to use the ordinal or interval timing likelihood (default is false).
// @param ncores number of threads to use for the parallelization (default is 1).
// @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false). Default is true.
// @param nsim number of draws from the posterior distribution. This argument will be used only if mu and vcov are supplied. For remstimate objects based on Bayesian methods, a sample from the posterior distribution will be drawn at R-level and supplied to this function.
//
// @return WAIC calculated on the deviance scale as -2 * elpdWAIC 
//
// [[Rcpp::export]]
double getWAIC(arma::vec mu,
                arma::mat vcov, 
                arma::mat pars,
                const arma::cube &stats,
                const arma::field<arma::uvec> &actor1,
                const arma::field<arma::uvec> &actor2,
                const arma::field<arma::uvec> &dyad,
                const arma::vec &interevent_time,
                const Rcpp::List &omit_dyad,
                std::string model,
                std::string approach,
                bool ordinal = false,
                int ncores = 1,
                bool senderRate = true,
                int nsim = 500){

  // select model
  std::vector<std::string> models = {"tie","actor"};
  std::vector<std::string>::iterator itr_model = std::find(models.begin(), models.end(), model);
  auto which_model = std::distance(models.begin(), itr_model);
  // select approach
  std::vector<std::string> approaches = {"Bayesian","Frequentist"};
  std::vector<std::string>::iterator itr_approach = std::find(approaches.begin(), approaches.end(), approach);
  auto which_approach = std::distance(approaches.begin(), itr_approach);
  if(which_approach == 1){ // generate parameters by approximating their posterior distribution with a Multivariate Normal distribution
    pars = arma::mvnrnd(mu,vcov,nsim);
  }

  // declaring dimensions
  arma::uword M = stats.n_slices; // number of events

  // initializing general indices
  arma::uword m,j;

  // creating object where to store partial output (useful for parallelization)
  arma::vec out(M,arma::fill::zeros);

  // without parallelization ... at the moment 

  // switch between tie-oriented and actor-oriented modeling
  if(which_model == 0){ // tie-oriented modeling (ordinal and interval likelihood supported)
    // omit dyad
    arma::vec riskset_time_vec(M); 
    arma::mat riskset_mat;
    if(omit_dyad.size()>0){
      riskset_time_vec = Rcpp::as<arma::vec>(omit_dyad["time"]);
      riskset_mat = Rcpp::as<arma::mat>(omit_dyad["riskset"]);
    }
    else{
      riskset_time_vec.fill(-1); // to simplify the ifelse in the loop below
    }

    #ifdef _OPENMP
    omp_set_dynamic(0);         // disabling dynamic teams
    omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
    #pragma omp parallel for if(ncores>1) private(m,j) shared(M,riskset_time_vec,riskset_mat,stats,dyad,pars,ordinal,out)
    #endif
    for(m = 0; m < M; m++)
      {
        int riskset_time_m = riskset_time_vec(m);
        arma::vec out_vec_loc_log(pars.n_cols,arma::fill::zeros);
        arma::vec out_vec_loc_p(pars.n_cols,arma::fill::zeros);
        arma::mat stats_m = stats.slice(m).t(); // dimensions : [U*D] we want to access dyads by column
        arma::uvec events_occurred = dyad(m)-1;
        double p_m = 0.0;
        for(j = 0; j < pars.n_cols; j++){
            double lpd_m_j = 0.0;
            // (1)  lpd computation
            arma::vec log_lambda = stats_m.t() * pars.col(j);
            lpd_m_j += arma::accu(log_lambda(events_occurred));
            // dealing with the risk set
            if(riskset_time_m!=(-1)){ // if the 'riskset_time_m' is different from (-1), then a dynamic riskset is observed
              arma::uvec events_at_risk = arma::find(riskset_mat.row(riskset_time_m));
              arma::vec sum_lambda_at_risk = riskset_mat.row(riskset_time_m) * arma::exp(log_lambda);
              if(ordinal){
                lpd_m_j -= log(sum_lambda_at_risk(0));
              }
              else{
                lpd_m_j -= sum_lambda_at_risk(0)*interevent_time(m);
              }
            }
            else{ // loop over all dyad (because for all the time points the riskset is fixed)
              double sum_lambda_at_risk = arma::accu(arma::exp(log_lambda));
              if(ordinal){
                lpd_m_j -= log(sum_lambda_at_risk);
              }
              else{
                lpd_m_j -= sum_lambda_at_risk*interevent_time(m);
              }
            }
            out_vec_loc_log[j] = lpd_m_j;
            out_vec_loc_p[j] = exp(lpd_m_j);
        }
        p_m += mean(out_vec_loc_p);
        out(m) = log(p_m) - var(out_vec_loc_log);
      }
  }
  else if(which_model == 1){ // actor-oriented modeling
    if(senderRate){ // if sender rate model
      // omit_dyad 
      arma::vec riskset_time_vec(M); 
      arma::mat riskset_mat;

      if(omit_dyad.size()>0){
        riskset_time_vec = Rcpp::as<arma::vec>(omit_dyad["time"]);
        riskset_mat = Rcpp::as<arma::mat>(omit_dyad["risksetSender"]);
      }
      else{
        riskset_time_vec.fill(-1); // to simplify the ifelse in the loop below
      }
      
      #ifdef _OPENMP
      omp_set_dynamic(0);         // disabling dynamic teams
      omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
      #pragma omp parallel for if(ncores>1) private(m,j) shared(M,riskset_time_vec,riskset_mat,stats,actor1,pars,ordinal,out)
      #endif
      for(m = 0; m < M; m++)
        {
          int riskset_time_m = riskset_time_vec(m);
          arma::vec out_vec_loc_log(pars.n_cols,arma::fill::zeros);
          arma::vec out_vec_loc_p(pars.n_cols,arma::fill::zeros);
          arma::mat stats_m = stats.slice(m).t(); // dimensions : [U*N] we want to access dyads by column
          arma::uvec sender = actor1(m)-1;
          double p_m = 0.0;
          for(j = 0; j < pars.n_cols; j++){
              double lpd_m_j = 0.0;
              // (1)  lpd computation
              arma::vec log_lambda = stats_m.t() * pars.col(j);
              lpd_m_j += arma::accu(log_lambda(sender));
              // dealing with the risk set
              if(riskset_time_m!=(-1)){ // if the 'riskset_time_m' is different from (-1), then a dynamic riskset is observed
                arma::uvec events_at_risk = arma::find(riskset_mat.row(riskset_time_m));
                arma::vec sum_lambda_at_risk = riskset_mat.row(riskset_time_m) * arma::exp(log_lambda);
                if(ordinal){
                  lpd_m_j -= log(sum_lambda_at_risk(0));
                }
                else{
                  lpd_m_j -= sum_lambda_at_risk(0)*interevent_time(m);
                }
              }
              else{ // loop over all dyad (because for all the time points the riskset is fixed)
                double sum_lambda_at_risk = arma::accu(arma::exp(log_lambda));
                if(ordinal){
                  lpd_m_j -= log(sum_lambda_at_risk);
                }
                else{
                  lpd_m_j -= sum_lambda_at_risk*interevent_time(m);
                }
              }
              out_vec_loc_log[j] = lpd_m_j;
              out_vec_loc_p[j] = exp(lpd_m_j);
          }
          p_m += mean(out_vec_loc_p);
          out(m) = log(p_m) - var(out_vec_loc_log);
        }
    }
    else{ // if receiver choice model
      // omit_dyad 
      arma::vec riskset_time_vec(M); 
      arma::mat riskset_mat;
      if(omit_dyad.size()>0){
        riskset_time_vec = Rcpp::as<arma::vec>(omit_dyad["time"]);
        riskset_mat = Rcpp::as<arma::mat>(omit_dyad["riskset"]);
      }
      else{
        riskset_time_vec.fill(-1); // to simplify the ifelse in the loop below
      }
      arma::uword N = stats.n_rows;

      #ifdef _OPENMP
      omp_set_dynamic(0);         // disabling dynamic teams
      omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
      #pragma omp parallel for if(ncores>1) private(m,j) shared(M,N,riskset_time_vec,riskset_mat,stats,actor1,actor2,pars,out)
      #endif
      for(m = 0; m < M; m++)
        {
          int riskset_time_m = riskset_time_vec(m);
          arma::vec out_vec_loc_log(pars.n_cols,arma::fill::zeros);
          arma::vec out_vec_loc_p(pars.n_cols,arma::fill::zeros);
          arma::mat stats_m = stats.slice(m).t(); // dimensions : [U*N] we want to access dyads by column
          arma::uvec sender = actor1(m)-1;
          arma::uvec receiver = actor2(m)-1;
          double p_m = 0.0;
          for(j = 0; j < pars.n_cols; j++){
            double lpd_m_j = 0.0;
            // (1)  lpd computation
            arma::vec log_lambda = stats_m.t() * pars.col(j);
            lpd_m_j += log_lambda(receiver(0));
            // dealing with the risk set
            arma::vec lambda_at_risk = arma::exp(log_lambda);
            if(riskset_time_m!=(-1)){ // if the 'riskset_time_m' is different from (-1), then a dynamic riskset is observed
              arma::vec riskset_sel = riskset_mat(riskset_time_m,arma::span(sender(0)*N,sender(0)*N+N-1)).t();
              lambda_at_risk %= riskset_sel; // if actors are not at risk to be receiver then there will be zeros in the vector (and this won't affect the sum)
            }
            lambda_at_risk(sender(0)) = 0.0;
            lpd_m_j -= log(arma::accu(lambda_at_risk));
            out_vec_loc_log[j] = lpd_m_j;
            out_vec_loc_p[j] = exp(lpd_m_j);
          }
          p_m += mean(out_vec_loc_p);
          out(m) = log(p_m) - var(out_vec_loc_log);
        }
    }
  }

  // calculating elpdWAIC per time point and WAIC as (-2)*elpdWAIC
  double waic = 0.0;
  waic = arma::accu(out);
  waic *= (-2.0);

  return waic;
}



// /////////////////////////////////////////////////////////////////////////////////
// ///////////(END)                      WAIC                     (END)/////////////
// /////////////////////////////////////////////////////////////////////////////////


// /////////////////////////////////////////////////////////////////////////////////
// ////////////(BEGIN)           Compute diagnostics            (BEGIN)/////////////
// /////////////////////////////////////////////////////////////////////////////////

// computeDiagnostics
//
// function that computes Schoenfeld's residuals and rates given the estimates of the model parameters and the array of statistics.
//
// @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats').
// @param stats is cube of M slices. Each slice is a matrix of dimensions D*U (or N*U for the actor-oriented model) with statistics of interest by column and dyads (actors, for the actor-oriented model) by row. If intercept is present then it will be excluded from the array
// @param actor1 list of actor1's observed per each time point (attr(reh,"actor1")-1).
// @param actor2 list of actor2's observed per each time point (attr(reh,"actor2")-1).
// @param dyad list of dyads observed per each time point (from the attribute attr(remify,"dyad")-1).
// @param omit_dyad is a list of two objects: a vector "time" and a matrix "riskset" (or "risksetSender" for the sender model). Two objects for handling changing risksets. The object is NULL if no change of the risk set structure is defined.
// @param model either "actor" or "tie" oriented model.
// @param N number of actors in the network. This argument is used only in the ReceiverChoice likelihood (model = "actor").
// @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false). Default is true.
// @param ncores number of threads to use for the parallelization (default is 1).
// @param baseline intercept parameter value. Default is 0.
//
// @return list of: residuals, smoothing weights (for estimating a spline on the residuals over time) and rates of actor/dyads over time.
//
// [[Rcpp::export]]
Rcpp::List computeDiagnostics(const arma::vec &pars,
                        const arma::cube &stats,
                        const arma::field<arma::uvec> &actor1,
                        const arma::field<arma::uvec> &actor2,
                        const arma::field<arma::uvec> &dyad,
                        const Rcpp::List &omit_dyad,
                        std::string model,
                        int N,
                        bool senderRate = true,
                        int ncores = 1,
                        double baseline = 0){

  Rcpp::List out = Rcpp::List::create(); // output list object
  std::vector<std::string> which_model = {"tie","actor"};
  arma::uword M = stats.n_slices;
  arma::uword P = stats.n_cols;  
  arma::uword D = stats.n_rows; // is the number of dyads in the tie-oriented model, and the number of actors in the actor-oriented model
  arma::uword m,p,u,z,j;
  int n;
  arma::field<arma::mat> residuals_std(M);
  arma::mat smoothing_weights(M,P,arma::fill::zeros);
  arma::field<arma::vec> rates(M); 

  if(model.compare(which_model[0]) == 0){ // model == "tie"

    arma::vec riskset_time_vec(M); 
    arma::umat riskset_mat;

    if(omit_dyad.size()>0){
      riskset_time_vec = Rcpp::as<arma::vec>(omit_dyad["time"]);
      riskset_mat = Rcpp::as<arma::umat>(omit_dyad["riskset"]);
    }
    else{
      riskset_time_vec.fill(-1);
    }


    #ifdef _OPENMP
    omp_set_dynamic(0);         // disabling dynamic teams
    omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
    #pragma omp parallel for if(ncores>1) private(m,p,u,z,j) shared(M,D,P,dyad,stats,pars,baseline,riskset_time_vec,riskset_mat,residuals_std,rates,smoothing_weights)
    #endif
    for(m = 0; m < M; m++){ // starting from m=1 because at m=0 most of the statistics are zero
      //arma::uword dyad_m = dyad(m);
      arma::uvec dyad_m = dyad(m)-1; // -1 because dyads' IDs must range between 0 and D-1

      arma::mat residuals_mat(dyad_m.n_elem,P,arma::fill::zeros);
      arma::mat residuals_std_m(dyad_m.n_elem,P,arma::fill::zeros);   // residuals_std(m).set_size(dyad_m.n_elem,P);
      arma::vec rates_m(dyad_m.n_elem,arma::fill::zeros);
      int riskset_time_m = riskset_time_vec(m);
      arma::vec expected_stats(P,arma::fill::zeros);
      arma::mat vcov(P,P,arma::fill::zeros); // we skip the baseline (intercept parameter)
      std::vector<int> which_stat_to_keep; //Rcpp::IntegerVector 
 
      arma::mat stats_m = stats.slice(m); // dimensions : [nDyads*nStats]
      arma::vec lambda = arma::exp((stats_m * pars) + baseline); // event rates 
      arma::vec ws = lambda/arma::accu(lambda);

      // calculating expected_stats and residuals_mat
      if(riskset_time_m!=(-1)){
        arma::rowvec riskset_mat_row = arma::conv_to<arma::rowvec>::from(riskset_mat.row(riskset_time_m));
        arma::uvec which_dyads = arma::find(riskset_mat_row); // by default arma::find returns the indices of non-zero elements
        arma::uvec which_p(1);
        rates_m = lambda(which_dyads);
        for(p = 0; p < P; p++){
          which_p(0) = p;
          expected_stats(p) = arma::accu(stats_m(which_dyads,which_p) % ws(which_dyads)); // stats_m.elem(which_dyads+p*D).t()
          residuals_mat.col(p) = stats_m(dyad_m,which_p);
          residuals_mat.col(p) -= expected_stats(p);
        }
      } 
     else{
        rates_m = lambda;
        arma::uvec which_p(1);
        for(p = 0; p < P; p++){
          expected_stats(p) = arma::accu(stats_m.col(p).t() * ws);
          which_p(0) = p;
          residuals_mat.col(p) = stats_m(dyad_m,which_p);
          residuals_mat.col(p) -= expected_stats(p);
        }
      }
      // calculating vcov matrix
      for(u = 0; u < P; u++){
        for(z = 0; z < (u+1); z++){
          if(riskset_time_vec(m)!=(-1)){
            arma::uvec which_dyads = arma::find(riskset_mat.row(riskset_time_vec(m))); 
            vcov(u,z) += sum((stats_m.elem(which_dyads+u*D).t() % stats_m.elem(which_dyads+z*D).t()) * ws(which_dyads));
          }
          else{
            vcov(u,z) += sum((stats_m.col(u).t() % stats_m.col(z).t()) * ws);
          }
          vcov(u,z) -= (expected_stats(u)*expected_stats(z));
          if(u==z){
            if(vcov(u,z)!=0.0){
              which_stat_to_keep.push_back(u);
            }
          }
          if(u != z){
            vcov(z,u) += vcov(u,z);
          }
        }
      }
      // claculate the inverse of the vcov matrix
      arma::uvec which_stats_to_keep_uvec = arma::conv_to<arma::uvec>::from(which_stat_to_keep);
      arma::mat vcov_sel = vcov(which_stats_to_keep_uvec,which_stats_to_keep_uvec);
      arma::mat inverse_vcov = arma::inv_sympd(vcov_sel,arma::inv_opts::allow_approx);
      for(j = 0 ; j < vcov_sel.n_rows; j++){
        int l = which_stat_to_keep[j];
        residuals_std_m.col(l) = residuals_mat.col(l) * inverse_vcov(j,j); // standardizing the residual
        smoothing_weights(m,l) = 1/inverse_vcov(j,j); // saving the inverse variance for the smoothing spline in R
      }
      rates(m) = rates_m;
      residuals_std(m) = residuals_std_m;
    }

    // saving output results
    out["standardized_residuals"] = residuals_std;
    out["smoothing_weights"] = smoothing_weights;
    out["rates"] = rates;
  }
  else if(model.compare(which_model[1]) == 0){ // model  == "actor"

    arma::vec riskset_time_vec(M); 
    arma::mat riskset_mat;

    if(omit_dyad.size()>0){
      riskset_time_vec = Rcpp::as<arma::vec>(omit_dyad["time"]);
      if(senderRate){
        riskset_mat = Rcpp::as<arma::mat>(omit_dyad["risksetSender"]);
      }
      else{
        riskset_mat = Rcpp::as<arma::mat>(omit_dyad["riskset"]);
      }
    }
    else{
      riskset_time_vec.fill(-1);
    }

    if(senderRate){ // sender rate model
      #ifdef _OPENMP
      omp_set_dynamic(0);         // disabling dynamic teams
      omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
      #pragma omp parallel for if(ncores>1) private(m,p,u,z,j) shared(M,P,D,stats,actor1,pars,baseline,riskset_time_vec,riskset_mat,residuals_std,smoothing_weights,rates)
      #endif
      for(m = 0; m < M; m++){
        arma::mat stats_m = stats.slice(m); // dimensions : [nActors*nStats]
        arma::vec expected_stats(P,arma::fill::zeros);
        arma::mat vcov(P,P,arma::fill::zeros); // we skip the baseline (intercept parameter)
        std::vector<int> which_stat_to_keep;
        arma::uvec actor1_m = actor1(m)-1; // sender at t[m], -1 because actors' IDs must range between 0 and N-1
        arma::vec lambda = arma::exp((stats_m * pars) + baseline); // event rates
        arma::vec ws = lambda/sum(lambda); // weights
        arma::mat residuals_mat(actor1_m.n_elem,P);
        arma::mat residuals_std_m(actor1_m.n_elem,P);
        int riskset_time_m = riskset_time_vec(m);
        arma::vec rates_m;
        // calculating expected_stats and residuals_mat
        if(riskset_time_m!=(-1)){
          arma::uvec which_actors = arma::find(riskset_mat.row(riskset_time_m)); 
          rates_m = lambda(which_actors);
          for(p = 0; p < P; p++){
            expected_stats(p) = sum(stats_m.elem(which_actors+p*D).t() * ws(which_actors));
            arma::uvec which_p(1);
            which_p(0) = p;
            residuals_mat.col(p) = stats_m(actor1_m,which_p);
            residuals_mat.col(p) -= expected_stats(p);
          }
        }
        else{
          rates_m = lambda;
          for(p = 0; p < P; p++){
            expected_stats(p) = sum(stats_m.col(p).t() * ws);
            arma::uvec which_p(1);
            which_p(0) = p;
            residuals_mat.col(p) = stats_m(actor1_m,which_p);
            residuals_mat.col(p) -= expected_stats(p);
          }
        }
        // calculating vcov matrix
        for(u = 0; u < P; u++){
          for(z = 0; z < (u+1); z++){
            if(riskset_time_m!=(-1)){
              arma::uvec which_dyads = arma::find(riskset_mat.row(riskset_time_m));
              vcov(u,z) += sum((stats_m.elem(which_dyads+u*D).t() % stats_m.elem(which_dyads+z*D).t()) * ws(which_dyads));
            }
            else{
              vcov(u,z) += sum((stats_m.col(u).t() % stats_m.col(z).t()) * ws);
            }
            vcov(u,z) -= (expected_stats(u)*expected_stats(z));
            if(u==z){
              if(vcov(u,z)!=0.0){
                which_stat_to_keep.push_back(u);
              }
            }
            if(u != z){
              vcov(z,u) += vcov(u,z);
            }
          }
        }
        // claculate the inverse of the vcov matrix
        arma::uvec which_stats_to_keep_uvec = arma::conv_to<arma::uvec>::from(which_stat_to_keep);
        arma::mat vcov_sel = vcov(which_stats_to_keep_uvec,which_stats_to_keep_uvec);
        arma::mat inverse_vcov = arma::inv_sympd(vcov_sel,arma::inv_opts::allow_approx);
        for(j = 0 ; j < vcov_sel.n_rows; j++){
          int l = which_stat_to_keep[j];
          residuals_std_m.col(l) = residuals_mat.col(l) * inverse_vcov(j,j); // standardizing the residual
          smoothing_weights(m,l) = 1/inverse_vcov(j,j); // saving the inverse variance for the smoothing spline in R
        }
        rates(m) = rates_m;
        residuals_std(m) = residuals_std_m;
      }
    }
    else{ // receiver choice model
      #ifdef _OPENMP
      omp_set_dynamic(0);         // disabling dynamic teams
      omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
      #pragma omp parallel for if(ncores>1) private(m,p,u,z,n,j) shared(M,P,N,stats,actor1,actor2,pars,baseline,riskset_time_vec,riskset_mat,residuals_std,smoothing_weights,rates)
      #endif
      for(m = 0; m < M; m++){
        arma::mat stats_m = stats.slice(m); // dimensions : [nActors*nStats]
        arma::vec expected_stats(P,arma::fill::zeros);
        arma::mat vcov(P,P,arma::fill::zeros); // we skip the baseline (intercept parameter)
        std::vector<int> which_stat_to_keep;
        // receiver choice model
        arma::uvec actor1_m_vec = actor1(m)-1; // -1 because actors' IDs must range between 0 and N-1
        arma::uvec actor2_m_vec = actor2(m)-1; // -1 because actors' IDs must range between 0 and N-1
        arma::uword actor1_m = actor1_m_vec(0); // this could be improved (?)
        int actor1_m_int = static_cast<int>(actor1_m);
        int riskset_time_m = riskset_time_vec(m);
        arma::uword actor2_m = actor2_m_vec(0); // this could be improved (?)
        arma::rowvec stats_actor2_m = stats_m.row(actor2_m); // perhaps it is necessary to use the transpose ? .t(); // we select it now before we reduce lambda and stats according to varying risk set
        arma::mat residuals_mat(1,P);
        arma::mat residuals_std_m(1,P);
        arma::vec lambda = arma::exp((stats_m * pars) + baseline); // event rates
        std::vector<int> remove_actors;
        arma::vec ws;
        remove_actors.push_back(actor1_m);
        
        if(riskset_time_m!=(-1)){
          for(n = 0; n<N; n++){
            int dyad_n = remify::getDyadIndex(actor1_m,n,0,N,true); // in the actor-oriented model, D is the number of actors
            if(n!=actor1_m_int && riskset_mat(riskset_time_m,dyad_n)==0){
              remove_actors.push_back(n);
            }
          }
        }
        lambda.shed_rows(arma::conv_to<arma::uvec>::from(remove_actors));
        stats_m.shed_rows(arma::conv_to<arma::uvec>::from(remove_actors));
        rates(m) = lambda; // selection of event rates
        ws = lambda/sum(lambda); // weights
        
        for(p = 0; p < P; p++){
          expected_stats(p) = sum(stats_m.col(p).t() * ws);
          residuals_mat.col(p) = stats_actor2_m(p);
          residuals_mat.col(p) -= expected_stats(p);
        }

        // calculating vcov matrix
        for(u = 0; u < P; u++){
          for(z = 0; z < (u+1); z++){
            vcov(u,z) += sum((stats_m.col(u).t() % stats_m.col(z).t()) * ws); // here we do not have to select anymore because we already removed the actors from the risk set
            vcov(u,z) -= (expected_stats(u)*expected_stats(z));
            if(u==z){
              if(vcov(u,z)!=0.0){
                which_stat_to_keep.push_back(u);
              }
            }
            if(u != z){
              vcov(z,u) = vcov(u,z);
            }
          }
        }

        // claculate the inverse of the vcov matrix
        arma::uvec which_stats_to_keep_uvec = arma::conv_to<arma::uvec>::from(which_stat_to_keep);
        arma::mat vcov_sel = vcov(which_stats_to_keep_uvec,which_stats_to_keep_uvec);
        arma::mat inverse_vcov = arma::inv_sympd(vcov_sel,arma::inv_opts::allow_approx);
        for(j = 0 ; j < vcov_sel.n_rows; j++){
          int l = which_stat_to_keep[j];
          residuals_std_m.col(l) = residuals_mat.col(l) * inverse_vcov(j,j); // standardizing the residual
          smoothing_weights(m,l) = 1/inverse_vcov(j,j); // saving the inverse variance for the smoothing spline in R
        }
        residuals_std(m) = residuals_std_m;
      }
    }
    // saving output results
    out["standardized_residuals"] = residuals_std;
    out["smoothing_weights"] = smoothing_weights;
    out["rates"] = rates;
  }

  return out;                           
}

// /////////////////////////////////////////////////////////////////////////////////
// /////////////(END)            Compute diagnostics            (END)///////////////
// /////////////////////////////////////////////////////////////////////////////////
