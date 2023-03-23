#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <omp.h>
#include <RcppArmadilloExtensions/sample.h> // used for the sample function inside RcppArmadillo
#include <typeinfo>
#include <map>
#include <iterator>
#include <string>
#include <remify.h>
//#include <progress.hpp> // progress bar
//#include <progress_bar.hpp> // progress bar

// /////////////////////////////////////////////////////////////////////////////////
// ///////////(BEGIN)              remDerivatives               (BEGIN)/////////////
// /////////////////////////////////////////////////////////////////////////////////


//' remDerivativesStandard
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values
//'
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param dyad is the vector of observed dyads from the 'remify' object (attr(remify,"dyad"))
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param ordinal boolean that indicate whether to use the ordinal or interval timing likelihood
//' @param ncores integer referring to the number of threads for the parallelization
//' @param gradient boolean true/false whether to return gradient value
//' @param hessian boolean true/false whether to return hessian value
//'
//' @return list of values: loglik, gradient, hessian
//'
// [[Rcpp::export]]
Rcpp::List remDerivativesStandard(const arma::vec &pars,
                                          const arma::cube &stats,
                                          const arma::uvec &dyad, 
                                          const Rcpp::List &omit_dyad,
                                          const arma::vec &interevent_time, // change to intereventTime
                                          bool ordinal = false,
                                          int ncores = 1,
                                          bool gradient = true,
                                          bool hessian = true
                                          ){
    arma::uword U = stats.n_cols; // number of parameters (alternative is stats.n_rows)
    arma::uword D = stats.n_rows;
    arma::uword M = stats.n_slices; // number of events

    arma::uword d,m,l,k;
    arma::vec log_lambda(D,arma::fill::zeros) ;

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
      omp_set_dynamic(0);         // disabling dynamic teams
      omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
      #pragma omp parallel for private(m,d,l,k,log_lambda) shared(M,D,U,loglik,hess,grad,stats,riskset_time_vec,riskset_mat,dyad,gradient,hessian,interevent_time)
      for(m = 0; m < M; m++){
          arma::mat stats_m = stats.slice(m).t(); // dimensions : [U*D] we want to access dyads by column
          log_lambda = stats_m.t() * pars;
          int riskset_time_m = riskset_time_vec(m);
          if(riskset_time_m!=(-1)){ // if the riskset_time_vec[m]=riskset_time_m is different from, then a dynamic riskset is observed
            for(d = 0; d < D; d++){
              if(riskset_mat(riskset_time_m,d)){ // ignoring impossible events that are not in risk set // [[try without ==1 ]]
                  if(dyad(m) == d){ // if the event occurred
                      loglik[m] += log_lambda.at(d);
                      grad.col(m) += stats_m.col(d);
                  }
                  double dtelp = exp(log_lambda.at(d))*interevent_time.at(m);  // change `dtelp` name to something else more understandable
                  loglik[m] -= dtelp;
                  if(gradient){
                      grad.col(m) -= stats_m.col(d)*dtelp;
                  }
                  if(hessian){
                    for(k = 0; k < U; k++){
                      for (l = k; l < U; l++){
                          hess(k,l,m) -= stats_m.at(l,d)*stats_m.at(k,d)*dtelp;
                          hess(l,k,m) = hess(k,l,m);
                      }
                    }
                  }

              }
            }
          }
          else{ // loop over all dyad (because for all the time points the riskset is fixed)
            for(d = 0; d < D; d++){
              if(dyad(m) == d){ // if event occured
                  loglik[m] += log_lambda.at(d);
                  grad.col(m) += stats_m.col(d);
              }
              double dtelp = exp(log_lambda.at(d))*interevent_time.at(m);  // change `dtelp` name to something else more understandable
              loglik[m] -= dtelp;
              if(gradient){
                  grad.col(m) -= stats_m.col(d)*dtelp;
              }
              if(hessian){
                for(k = 0; k < U; k++){
                  for (l = k; l < U; l++){
                      hess(k,l,m) -= stats_m.at(l,d)*stats_m.at(k,d)*dtelp;
                      hess(l,k,m) = hess(k,l,m);
                  }
                }
              }
            }
          }
      }
    }
    else{ // ordinal likelihood
      omp_set_dynamic(0);         // disabling dynamic teams
      omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
      #pragma omp parallel for private(m,d,l,k,log_lambda) shared(M,D,U,loglik,hess,grad,stats,riskset_time_vec,riskset_mat,dyad,gradient,hessian,interevent_time)
      for(m = 0; m < M; m++){
        arma::mat stats_m = stats.slice(m).t(); // dimensions : [U*D] we want to access dyads by column
        log_lambda = stats_m.t() * pars;
        double surv = 0.0;
        arma::vec grad_m(U,arma::fill::zeros);
        arma::mat hess_m(U,U,arma::fill::zeros);
        int riskset_time_m = riskset_time_vec[m];
        if(riskset_time_m!=(-1)){ // if the riskset_time_vec[m]=riskset_time_m is different from, then a dynamic riskset is observed
          for(d = 0; d < D; d++){
            if(riskset_mat(riskset_time_m,d)){ // ignoring impossible events that are not in risk set // [[try without ==1 ]]
              if(dyad(m) == d){ // if event occured
                  loglik[m] += log_lambda.at(d);
              }
              double lambda_d = exp(log_lambda.at(d));
              surv += lambda_d;  // change `dtelp` name to something else more understandable
              if(gradient){
                grad.col(m) += stats_m.col(d); // first component of the gradient
                grad_m += (lambda_d * stats_m.col(d)); // building up the second component of the gradient
              }
              if(hessian){
                for(k = 0; k < U; k++){
                  for (l = k; l < U; l++){
                    hess_m(k,l) -= stats_m.at(l,d)*stats_m.at(k,d)*lambda_d;
                    hess_m(l,k) = hess_m(k,l);
                  }
                }
              }
            }    
          }
        }
      else{ // loop over all the dyads
        for(d = 0; d < D; d++){
          if(dyad(m) == d){ // if event occured
              loglik[m] += log_lambda.at(d);
          }
          double lambda_d = exp(log_lambda.at(d));
          surv += lambda_d;  // change `dtelp` name to something else more understandable
          if(gradient){
            grad.col(m) += stats_m.col(d); // first component of the gradient
            grad_m += (lambda_d * stats_m.col(d)); // building up the second component of the gradient
          }
          if(hessian){
            for(k = 0; k < U; k++){
              for (l = k; l < U; l++){
                hess_m(k,l) -= stats_m.at(l,d)*stats_m.at(k,d)*lambda_d;
                hess_m(l,k) = hess_m(k,l);
              }
            }
          }   
        }
      }
      loglik[m] -= log(surv);
      if(gradient){
        grad_m /= surv;
        grad.col(m) -= grad_m;
      }
      if(hessian){
        for(k = 0; k < U; k++){
          for (l = k; l < U; l++){
            hess_m(k,l) /= surv;
            hess_m(k,l) += (grad_m(k) * grad_m(l))/(std::pow(surv,2));
            hess_m(l,k) = hess_m(k,l);
          }
        }
        hess.slice(m) = hess_m;
      }

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


//' remDerivativesSenderRates
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values for estimating the sender rate parameters for the actor oriented model
//'
//' @param pars is a vector of parameters
//' @param stats is cube of M slices. Each slice is a matrix of dimensions N*U with statistics of interest by column and senders by row.
//' @param actor1 vector of actor1's (column reh$edgelist$actor1)
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param C number of event types 
//' @param D number of dyads
//' @param ordinal boolean, true if the likelihood to use is the ordinal one, interval otherwise
//' @param gradient boolean true/false whether to return gradient value
//' @param hessian boolean true/false whether to return hessian value
//'
//'
//' @return list of values: loglik, grad, fisher information
//'
// [[Rcpp::export]]
Rcpp::List remDerivativesSenderRates(
        const arma::vec &pars,
        const arma::cube &stats,
        const arma::uvec &actor1,         
        const Rcpp::List &omit_dyad,
        const arma::vec &interevent_time,
        int C,
        int D,
        bool ordinal = false,
        bool gradient = true,
        bool hessian  = true){

    int P = stats.n_cols; // number of parameters for send stats
    int M = stats.n_slices; // number of events
    int N = stats.n_rows; //Number of actors
    int n, m;

    // variables for internal computation
    arma::vec lambda_s(N);
    arma::vec expected_stat_m(P);
    Rcpp::IntegerVector dyad_m(3);
    arma::mat fisher_m(P,P);
    double sum_lambda = 0.0;

    //output
    double loglik = 0.0;
    arma::mat fisher(P,P,arma::fill::zeros);
    arma::vec grad(P,arma::fill::zeros);

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

    //loop through all events
    for(m = 0; m < M; m++){
      arma::mat stats_m = stats.slice(m).t(); // dimensions: [P * N]

      //this is exp(beta^T X) dot product
      lambda_s = arma::exp(stats_m.t() * pars);

      //actors in remify edgelist are 0 indexed hence no -1
      arma::uword sender = actor1(m);
      // int type = dyad[2]; // when type will be integrated in the function

      //reset internal variables
      if(gradient | hessian){
        fisher_m.zeros();
        expected_stat_m.zeros();
      }

      // riskset_m
      int riskset_time_m = riskset_time_vec(m);

      loglik += std::log(lambda_s(sender));
      if(gradient){
        grad += stats_m.col(sender);
      }

      // changes in the riskset at m-th event
      if(riskset_time_m!=(-1)){
        sum_lambda = sum(riskset_mat.row(riskset_time_m).t() % lambda_s); 
        if(ordinal){
          loglik-= std::log(sum_lambda);
          if(gradient | hessian){
            for(n = 0; n < N; n++){         //loop throughout all actors
              if(riskset_mat(riskset_time_m,n) == 1){
                expected_stat_m += lambda_s(n) * (stats_m.col(n));
                expected_stat_m /= sum_lambda; //exp(params_s * X_i)*X_i / sum_h (exp(params * X_h))
                if(gradient){
                  grad -= expected_stat_m;
                }
                if(hessian){
                  fisher_m += (expected_stat_m * expected_stat_m.t());
                  fisher_m -= (lambda_s(n) / sum_lambda) *( stats_m.col(n) * (stats_m.col(n).t()));
                  fisher += fisher_m;
                }   
              }
            }
          }
        } 
        else{
          loglik -=  sum_lambda * interevent_time(m);
          if(gradient){
            grad -= interevent_time(m) * stats_m * (riskset_mat.row(riskset_time_m).t() % lambda_s); 
          }
          if(hessian){
            for(n = 0; n < N; n++){         //loop throughout all actors
              fisher_m += (riskset_mat(riskset_time_m,n) * lambda_s(n)) * (stats_m.col(n) * (stats_m.col(n).t()) );
            }
            fisher -= interevent_time(m) * fisher_m;
          }
        }
      }
      else{  //no dynamic riskset
        sum_lambda = sum(lambda_s);
        if(ordinal){
          loglik-= std::log(sum_lambda);
          if(gradient | hessian){
            for(n = 0; n < N; n++){         //loop throughout all actors
              expected_stat_m += lambda_s(n) * (stats_m.col(n));
              expected_stat_m /= sum_lambda; //exp(params_s * X_i)*X_i / sum_h (exp(params * X_h))
              if(gradient){
                grad -= expected_stat_m;
              }
              if(hessian){
                fisher_m += (expected_stat_m * expected_stat_m.t());
                fisher_m -= (lambda_s(n) / sum_lambda) * (stats_m.col(n) * (stats_m.col(n).t()));
                fisher += fisher_m;
              }
            }
          }
        } 
        else{
          loglik -=  sum_lambda * interevent_time(m);
          if(gradient){
            grad -= interevent_time(m) * stats_m * lambda_s; //(6)
          }
          if(hessian){
            for(n = 0; n < N; n++){         //loop throughout all actors
              fisher_m += lambda_s(n)*(stats_m.col(n) * (stats_m.col(n).t()));
            }
            fisher -= interevent_time(m) * fisher_m;
          }
        }
      }
    }

    if(gradient && !hessian){
      return Rcpp::List::create(Rcpp::Named("value") = -loglik, Rcpp::Named("gradient") = -grad);
    }else if(!gradient && !hessian){
      return Rcpp::List::create(Rcpp::Named("value") = -loglik);
    }else{
      return Rcpp::List::create(Rcpp::Named("value") = -loglik, Rcpp::Named("gradient") = -grad, Rcpp::Named("hessian") = -fisher);
    }
}


//' remDerivativesReceiverChoice
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values for estimating the receiver choice parameters for the actor oriented model
//'
//' @param pars is a vector of parameters
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param actor1 vector of actor1's (column reh$edgelist$actor1)
//' @param actor2 vector of actor2's (column reh$edgelist$actor2)
//' @param omit_dyad, list object that takes care of the dynamic rikset (if defined)
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param N the number of actors
//' @param C number of event types 
//' @param D number of dyads
//' @param gradient boolean true/false whether to return gradient value
//' @param hessian boolean true/false whether to return hessian value
//'
//'
//' @return list of values: loglik, grad, fisher
//'
//[[Rcpp::export]]
Rcpp::List remDerivativesReceiverChoice(
        const arma::vec &pars,
        const arma::cube &stats,
        const arma::uvec &actor1,
        const arma::uvec &actor2,
        const Rcpp::List &omit_dyad,
        const arma::vec &interevent_time,
        int N,
        int C,
        int D,
        bool gradient = true,
        bool hessian  = true){

    arma::uword P_d = stats.n_cols; // number of parameters for dyad stats
    arma::uword M = actor2.n_elem; // number of events
    int n;

    // variables for internal computation
    arma::vec lambda_d(N);
    arma::vec expected_stat_m(P_d);
    Rcpp::IntegerVector dyad_m(3);
    arma::mat fisher_m(P_d,P_d);
    arma::uword dyad = 0;
    double denom = 0;

    //output
    double loglik = 0.0;
    arma::mat fisher(P_d,P_d,arma::fill::zeros);
    arma::vec grad(P_d,arma::fill::zeros);

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


    //loop through all events
    for(arma::uword m = 0; m < M; m++){
      arma::mat stats_m = stats.slice(m).t(); // dimensions : [P*N] we want to access dyads by column

      //this is exp(beta^T X) dot product
      lambda_d = arma::exp(stats_m.t() * pars);

      //actors
      arma::uword sender = actor1(m); 
      arma::uword receiver = actor2(m);
      //int type = dyad_m[2]; // when type will be integrated in the function

      //reset internal variables
      if(gradient | hessian){
        fisher_m.zeros();
        expected_stat_m.zeros();
      }
      denom = 0;

      // riskset_m
      int riskset_time_m = riskset_time_vec(m);

      // changes in the riskset at m-th event
      if(riskset_time_m!=(-1)){ 
        for(n = 0; n<N; n++){
            dyad = remify::getDyadIndex(sender,n,0,N,true);
            if(n!=sender && riskset_mat(riskset_time_m,dyad) == 1){ // dynamic riskset
                //loglik
                denom += lambda_d(n); // exp(param_d * X_sender_i) (4)
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
          if(n!=sender){ // 
              //loglik
              denom += lambda_d(n); // exp(param_d * X_sender_i) (4)
              if(hessian){
                fisher_m -= lambda_d(n)*(stats_m.col(n) * stats_m.col(n).t());
              }
              if(gradient | hessian){
              expected_stat_m += lambda_d(n) * stats_m.col(n);
              }
          }
        }
      }

      loglik +=  std::log(lambda_d(receiver)); //(3) params_d^T * X_sender_receiver
      loglik -= std::log(denom);

      if(gradient | hessian){
        expected_stat_m /= denom; //(8)
      }
      if(gradient){
        grad += stats_m.col(receiver); //(7)
        grad -= expected_stat_m;
      }
      if(hessian){
        fisher_m /= denom;
        fisher_m += (expected_stat_m * expected_stat_m.t());
        fisher += fisher_m;
      }     
    }

    if(gradient && !hessian){
      return Rcpp::List::create(Rcpp::Named("value") = -loglik, Rcpp::Named("gradient") = -grad);
    }else if(!gradient && !hessian){
      return Rcpp::List::create(Rcpp::Named("value") = -loglik);
    }else{
      return Rcpp::List::create(Rcpp::Named("value") = -loglik, Rcpp::Named("gradient") = -grad, Rcpp::Named("hessian") = -fisher);
    }
}


//' remDerivatives 
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values
//'
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param actor1 vector of actor1's (column reh$edgelist$actor1)
//' @param actor2 vector of actor2's (column reh$edgelist$actor2)
//' @param dyad vector of dyad (from the attribute attr(remify,"dyad"))
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param model either "actor" or "tie" model
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param ncores number of threads to use for the parallelization
//' @param gradient boolean true/false whether to return gradient value
//' @param hessian boolean true/false whether to return hessian value
//' @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false)
//' @param gradient boolean true/false whether to return gradient value
//' @param hessian boolean true/false whether to return hessian value
//' @param N number of actors. This argument is used only in the ReceiverChoice likelihood (model = "actor")
//' @param C number of event types 
//' @param D number of dyads
//'
//' @return list of values: loglik, gradient, hessian
//'
// [[Rcpp::export]]
Rcpp::List remDerivatives(const arma::vec &pars,
                                  const arma::cube &stats,
                                  const arma::uvec &actor1,
                                  const arma::uvec &actor2,
                                  const arma::uvec &dyad,
                                  const Rcpp::List &omit_dyad,
                                  const arma::vec &interevent_time,
                                  std::string model,
                                  bool ordinal = false,
                                  int ncores = 1,
                                  bool gradient = true,
                                  bool hessian = true,
                                  bool senderRate = true,
                                  Rcpp::Nullable<int> N = R_NilValue,
                                  Rcpp::Nullable<int> C = R_NilValue,
                                  Rcpp::Nullable<int> D = R_NilValue){
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
                  out = remDerivativesReceiverChoice(pars,stats,actor1,actor2,omit_dyad,interevent_time,Rcpp::as<int>(N),Rcpp::as<int>(C),Rcpp::as<int>(D),gradient,hessian);
                  break;
                  }
        case 1 : {
                  out = remDerivativesSenderRates(pars,stats,actor1,omit_dyad,interevent_time,Rcpp::as<int>(C),Rcpp::as<int>(D),ordinal,gradient,hessian);
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


//' GDADAMAX
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values
//'
//' @param pars parameters
//' @param stats array of statistics
//' @param actor1 vector of actor1's (column reh$edgelist$actor1)
//' @param actor2 vector of actor2's (column reh$edgelist$actor2)
//' @param dyad vector of dyad (from the attribute attr(remify,"dyad"))
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined
//' @param interevent_time vector of interevent times (inside the reh object)
//' @param model either "actor" or "tie" model
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false)
//' @param gradient boolean true/false whether to return gradient value
//' @param hessian boolean true/false whether to return hessian value
//' @param N number of actors. This argument is used only in the ReceiverChoice likelihood (model = "actor")
//' @param C number of event types 
//' @param D number of dyads
//' @param ncores number of threads to use for the parallelization
//' @param epochs number of epochs
//' @param learning_rate learning rate
//' @param beta1 hyperparameter beta1
//' @param beta2 hyperparameter beta2
//' @param epsilon hyperparameter eta
//'
//' @return optimization with GDADAM
//'
//' @export
// [[Rcpp::export]]
Rcpp::List GDADAMAX(const arma::vec &pars,
                  const arma::cube &stats,
                  const arma::uvec &actor1,
                  const arma::uvec &actor2,
                  const arma::uvec &dyad,
                  const Rcpp::List &omit_dyad,
                  const arma::vec &interevent_time,
                  std::string model,
                  bool ordinal = false,
                  bool senderRate = true,
                  bool gradient = true,
                  bool hessian = false,
                  Rcpp::Nullable<int> N = R_NilValue,
                  Rcpp::Nullable<int> C = R_NilValue,
                  Rcpp::Nullable<int> D = R_NilValue,
                  int ncores = 1,
                  int epochs = 1e03,
                  double learning_rate = 0.002,
                  double beta1 = 0.9,
                  double beta2 = 0.999,
                  double epsilon = 0.01){

  // link ref: https://machinelearningmastery.com/gradient-descent-optimization-with-adamax-from-scratch/   (cite properly original papers)               
  arma::uword P = pars.n_elem;
  double  alpha = 0.002; 
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
    Rcpp::List derv = remDerivatives(pars_loc,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,gradient,hessian,senderRate,N,C,D);

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
    for(int p=0; p<P; p++){
      arma::vec inf_norm_both = {beta2*inf_norm_prev(p),abs_grad(p)};
      inf_norm(p) = arma::max(inf_norm_both);
    }

    // updating parameter value
    step_size(0) = alpha/(1-std::pow(beta1,i+1));

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


//' logPostHMC
//'
//' This function calculates the value of the log-posterior density given the loglikelihood and the log-prior density
//'
//' @param meanPrior is a vector of prior means with the same dimension as the vector of parameters
//' @param sigmaPrior is a matrix, I have been using a diagonal matrix here with the same dimension as the vector os parameters
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param actor1 vector of actor1's (column reh$edgelist$actor1)
//' @param actor2 vector of actor2's (column reh$edgelist$actor2)
//' @param dyad vector of dyad (from the attribute attr(remify,"dyad"))
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param model either "actor" or "tie" model
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param ncores number of threads to use for the parallelization
//' @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false)
//' @param N number of actors. This argument is used only in the ReceiverChoice likelihood (model = "actor")
//' @param C number of event types 
//' @param D number of dyads
//'
//' @return value of log-posterior density
//'
// [[Rcpp::export]]
double logPostHMC(const arma::vec &meanPrior,
                  const arma::mat &sigmaPrior,
                  const arma::vec &pars,
                  const arma::cube &stats,
                  const arma::uvec &actor1,
                  const arma::uvec &actor2,
                  const arma::uvec &dyad,
                  const Rcpp::List &omit_dyad,
                  const arma::vec &interevent_time,
                  std::string model,
                  bool ordinal = false,
                  int ncores = 1,
                  bool senderRate = true,
                  Rcpp::Nullable<int> N = R_NilValue,
                  Rcpp::Nullable<int> C = R_NilValue,
                  Rcpp::Nullable<int> D = R_NilValue){

  Rcpp::List derv = remDerivatives(pars,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,false,false,senderRate,N,C,D);
  double derv_0 = Rcpp::as<double>(derv[0]);                
  double prior =  arma::accu(0.5 * (pars.t() - meanPrior.t()) * inv(sigmaPrior) * (pars - meanPrior));
  return (prior + derv_0); 
}


//' logPostGradientHMC
//'
//' This function calculates the value of the gradient of the log-posterior density
//'
//' @param meanPrior is a vector of prior means with the same dimension as the vector of parameters
//' @param sigmaPrior is a matrix, I have been using a diagonal matrix here with the same dimension as the vector os parameters
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param actor1 vector of actor1's (column reh$edgelist$actor1)
//' @param actor2 vector of actor2's (column reh$edgelist$actor2)
//' @param dyad vector of dyad (from the attribute attr(remify,"dyad"))
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param model either "actor" or "tie" model
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param ncores number of threads to use for the parallelization
//' @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false)
//' @param N number of actors. This argument is used only in the ReceiverChoice likelihood (model = "actor")
//' @param C number of event types 
//' @param D number of dyads
//'
//' @return value of log-posterior gradient
//'
//' @export
// [[Rcpp::export]]
arma::vec logPostGradientHMC(const arma::vec &meanPrior,
                              const arma::mat &sigmaPrior,
                              const arma::vec &pars,
                              const arma::cube &stats,
                              const arma::uvec &actor1,
                              const arma::uvec &actor2,
                              const arma::uvec &dyad,
                              const Rcpp::List &omit_dyad,
                              const arma::vec &interevent_time,
                              std::string model,
                              bool ordinal = false,
                              int ncores = 1,
                              bool senderRate = true,
                              Rcpp::Nullable<int> N = R_NilValue,
                              Rcpp::Nullable<int> C = R_NilValue,
                              Rcpp::Nullable<int> D = R_NilValue){
  Rcpp::List derv = remDerivatives(pars,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,true,false,senderRate,N,C,D);
  arma::vec gprior = inv(sigmaPrior)*(pars - meanPrior); // the sign is already changed here
  arma::vec glp = Rcpp::as<arma::vec>(derv[1]);
  return (glp + gprior); 
}

                                  
//' iterHMC
//'
//' This function does one iteration of the Hamiltonian Monte carlo
//'
//' @param L number of leapfrogs. Default (and recommended) value is 100.
//' @param epsilon size of the leapfrog. Default value is 1e-02.
//' @param meanPrior is a vector of prior means with the same dimension as the vector of parameters
//' @param sigmaPrior is a matrix, I have been using a diagonal matrix here with the same dimension as the vector os parameters
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param actor1 vector of actor1's (column reh$edgelist$actor1)
//' @param actor2 vector of actor2's (column reh$edgelist$actor2)
//' @param dyad vector of dyad (from the attribute attr(remify,"dyad"))
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param model either "actor" or "tie" model
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param ncores number of threads to use for the parallelization
//' @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false)
//' @param N number of actors. This argument is used only in the ReceiverChoice likelihood (model = "actor")
//' @param C number of event types 
//' @param D number of dyads
//'
// [[Rcpp::export]]
arma::field<arma::vec> iterHMC(arma::uword L,
                  double epsilon,
                  const arma::vec &meanPrior,
                  const arma::mat &sigmaPrior,
                  const arma::vec &pars,
                  const arma::cube &stats,
                  const arma::uvec &actor1,
                  const arma::uvec &actor2,
                  const arma::uvec &dyad,
                  const Rcpp::List &omit_dyad,
                  const arma::vec &interevent_time,
                  std::string model,
                  bool ordinal = false,
                  int ncores = 1,
                  bool senderRate = true,
                  Rcpp::Nullable<int> N = R_NilValue,
                  Rcpp::Nullable<int> C = R_NilValue,
                  Rcpp::Nullable<int> D = R_NilValue){

  arma::uword P = pars.size(); //number of parameters
  arma::field<arma::vec> out(2); // output object

  // proposing a new value for the momentum variable
  arma::vec rP = arma::randn(P); 
  // current momentum variable
  arma::vec rC = rP;
  //Rcpp::Rcout << "rnorm" << rP.t() << "\n";
  
  arma::vec betaC = pars; // current beta
  arma::vec betaP = pars; // new proposed beta (starting value is current beta)

  //leapfrog algorithm, updates via Hamiltonian equations
  rP -= (0.5 * epsilon * logPostGradientHMC(meanPrior,sigmaPrior,betaP,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D));
  for(arma::uword i = 1; i <= L; i++){
    betaP += (epsilon * rP);
    if(i != L) rP -= (epsilon * logPostGradientHMC(meanPrior,sigmaPrior,betaP,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D));
  }
  rP -= (0.5 * epsilon * logPostGradientHMC(meanPrior,sigmaPrior,betaP,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D));
  rP *= (-1);

  //computes final quantities for the acceptance rate
  double U = logPostHMC(meanPrior,sigmaPrior,betaC,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D);
  double propU = logPostHMC(meanPrior,sigmaPrior,betaP,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D);
  double K = 0.5 * arma::accu(arma::pow(rC,2));
  double propK = 0.5 * arma::accu(arma::pow(rP,2));

  //Accepting or rejecting the proposal
  arma::vec randomUnif = arma::randu(1);
  //Rcpp::Rcout << "randomunif: " << randomUnif.t() << "\n";
  if(randomUnif(0) < exp(U - propU + propK - K)){   // before it was: exp( (-propU - propK)/(-U - K) ) .. exp(propH)/exp(H)
    out(0) = betaP;
    out(1) = propU;
    // propU;
  } else {
    out(0) = betaC;
    out(1) = U;
    // U;
  }

  Rcpp::Rcout << "accepted value: " << out(0) << "\n";

  return out;
}


//' burninHMC (to check whether this function experiences issues with the definition of int rows and the following codings)
//'
//' This function performs the burn-in and the thinning at the end of the HMC
//'
//' @param samples cube with final draws
//' @param loglik matrix of values of the posterior loglikelihood at the different draws
//' @param burnin is the number of draws to discard after running the chains
//' @param thin is the number of draws to be skipped. For instance, if thin = 10, draws will be selected every 10 generated draws: 1, 11, 21, 31, ...
//'
//' @return list of two objects: draws and loglik after burnin and thinning step
//'
// [[Rcpp::export]]
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


//' HMC
//'
//' This function performs the Hamiltonian Monte Carlo
//'
//' @param pars_init is a matrix of dimensions U x nchains where for each column (chain) a random vector of initial values for the parameter is supplied.
//' @param nsim is the number of samples from the posterior that have to be generated.
//' @param nchains number of chains of length nsim
//' @param burnin is the number of draws to discard after running the chains
//' @param meanPrior is a vector of prior means with the same dimension as the vector of parameters
//' @param sigmaPrior is a matrix, I have been using a diagonal matrix here with the same dimension as the vector os parameters
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param actor1 vector of actor1's (column reh$edgelist$actor1)
//' @param actor2 vector of actor2's (column reh$edgelist$actor2)
//' @param dyad vector of dyad (from the attribute attr(remify,"dyad"))
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined//' @param interevent_time the time difference between the current time point and the previous event time.//' @param 
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param model either "actor" or "tie" model
//' @param ordinal logic TRUE/FALSE
//' @param ncores number of threads to use for the parallelization
//' @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false)
//' @param N number of actors. This argument is used only in the ReceiverChoice likelihood (model = "actor")
//' @param C number of event types 
//' @param D number of dyads
//' @param thin is the number of draws to be skipped. For instance, if thin = 10, draws will be selected every 10 generated draws: 1, 11, 21, 31, ...
//' @param L number of leapfrogs. Default (and recommended) value is 100.
//' @param epsilon size of the leapfrog. Default value is 1e-02.
//' @param ncores number of threads for parallel computing (default = 1)
//'
//' @return posterior draws
//'
// [[Rcpp::export]]
Rcpp::List HMC(arma::mat pars_init,
                arma::uword nsim,
                arma::uword nchains,
                arma::uword burnin,
                const arma::vec& meanPrior,
                const arma::mat& sigmaPrior,
                const arma::cube &stats,
                const arma::uvec &actor1,
                const arma::uvec &actor2,
                const arma::uvec &dyad,
                const Rcpp::List &omit_dyad,
                const arma::vec &interevent_time,
                std::string model,
                bool ordinal = false,
                int ncores = 1,
                bool senderRate = true,
                Rcpp::Nullable<int> N = R_NilValue,
                Rcpp::Nullable<int> C = R_NilValue,
                Rcpp::Nullable<int> D = R_NilValue,
                arma::uword thin = 1,
                arma::uword L = 100,
                double epsilon = 0.01){ 
  arma::cube array_of_draws(nsim, pars_init.n_rows, nchains); // array of draws from the posterior distribution (chains are by slice)
  arma::mat matrix_of_loglik(nsim,nchains); // matrix of posterior loglikelihood values across chains
  arma::uword j,i;
  Rcpp::List out = Rcpp::List::create(); // output object
  for(j = 0; j < nchains; j++){ //looping through chains
    arma::mat chain_j(pars_init.n_rows,nsim,arma::fill::zeros);
    arma::vec chain_j_loglik(nsim,arma::fill::zeros);
    //[i=0] this step only get the first sample out of the starting value
    arma::field<arma::vec> iter_i_hmc = iterHMC(L,epsilon,meanPrior,sigmaPrior,pars_init.col(j),stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D);
    chain_j.col(0) = iter_i_hmc[0]; // saving draws
    chain_j_loglik(0) = arma::conv_to<double>::from(iter_i_hmc[1]); // saving posterior loglikelihood
    for(i = 1; i < nsim; i++){ //looping through iterations of the MCMC
        //Then the next step will always be based on the previous one
        arma::field<arma::vec> iter_i_hmc = iterHMC(L,epsilon,meanPrior,sigmaPrior,chain_j.col(i-1),stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D);
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



//' emp_dist_longest_batch
//'
//' This function tunes the parameter L for the leapfrog step in the HMC method
//'
//' @param L number of leapfrogs. Default (and recommended) value is 100.
//' @param epsilon size of the leapfrog. Default value is 1e-02.
//' @param meanPrior is a vector of prior means with the same dimension as the vector of parameters
//' @param sigmaPrior is a matrix, I have been using a diagonal matrix here with the same dimension as the vector os parameters
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param actor1 vector of actor1's (column reh$edgelist$actor1)
//' @param actor2 vector of actor2's (column reh$edgelist$actor2)
//' @param dyad vector of dyad (from the attribute attr(remify,"dyad"))
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined//' @param interevent_time the time difference between the current time point and the previous event time.//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param model either "actor" or "tie" model
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param ncores number of threads to use for the parallelization
//' @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false)
//' @param N number of actors. This argument is used only in the ReceiverChoice likelihood (model = "actor")
//' @param C number of event types 
//' @param D number of dyads
//'
//' @export
// [[Rcpp::export]]
arma::vec emp_dist_longest_batch(arma::uword L,
                  double epsilon,
                  const arma::vec &meanPrior,
                  const arma::mat &sigmaPrior,
                  const arma::vec &pars,
                  const arma::cube &stats,
                  const arma::uvec &actor1,
                  const arma::uvec &actor2,
                  const arma::uvec &dyad,
                  const Rcpp::List &omit_dyad,
                  const arma::vec &interevent_time,
                  std::string model,
                  bool ordinal = false,
                  int ncores = 1,
                  bool senderRate = true,
                  Rcpp::Nullable<int> N = R_NilValue,
                  Rcpp::Nullable<int> C = R_NilValue,
                  Rcpp::Nullable<int> D = R_NilValue){
  arma::uword P = pars.size(); //number of parameters
  arma::vec beta_k = pars; // starting point (theta_0)
  arma::vec L_vec(L,arma::fill::zeros);
  arma::uword i,k;

  for(k = 0; k < L; k++){
    // (1) generate momentum
    arma::vec r_k = arma::randn(P);
    arma::vec r_star = r_k;
    arma::vec beta_star = beta_k;
    // (2) compute longest batch
    //arma::mat inv_sigmaPrior = inv(sigmaPrior);
    arma::vec beta0 = beta_k;
    arma::vec r0 = r_k;
    arma::vec beta_diff = (beta0-beta_k);
    arma::uword l = 0;
    double condition_while = arma::accu(beta_diff.t()*r0);
    while(condition_while>=0){
      arma::vec dd= logPostGradientHMC(meanPrior,sigmaPrior,beta0,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D);
      Rcpp::Rcout << "postgr: " << dd.t() << "\n";
      //Rcpp::Rcout << "beta_diff: " << beta_diff.t() << "\n r0:" << r0 << "\n";
      l+=1;
      // START compute leapfrog L=1
      // starting value for r0
      r0 -= 0.5 * epsilon * logPostGradientHMC(meanPrior,sigmaPrior,beta0,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D);
      // iteration L=1
      beta0 += epsilon * r0;
      r0 -= epsilon * logPostGradientHMC(meanPrior,sigmaPrior,beta0,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D);
      // updating r0 for the output
      r0 += 0.5 * epsilon * logPostGradientHMC(meanPrior,sigmaPrior,beta0,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D);
      // END compute leapfrog
      if(l == L){
        beta_star = beta0;
        r_star = r0;
        break;
      }
      arma::vec beta_diff = (beta0-beta_k);
      condition_while = arma::accu(beta_diff.t()*r0);
    }
    //(betaP,r,l)
    arma::uword L_k = l;
    L_vec(k) = L_k;
    // (3) if longest batch L_k is less than L_0, then compute leapfrog on L_0 - L_k
    if(L_k < L){ // compute leapfrog
      // START compute leapfrog
      r_star -= 0.5 * epsilon * logPostGradientHMC(meanPrior,sigmaPrior,beta_star,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D);
      for(i = 1; i <= (L-L_k); i++){
        beta_star += epsilon * r_star;
        r_star -= epsilon * logPostGradientHMC(meanPrior,sigmaPrior,beta_star,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D);
      }
      r_star += 0.5 * epsilon * logPostGradientHMC(meanPrior,sigmaPrior,beta_star,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D);
      // END compute leapfrog
    }
    r_star = - r_star;
    // (4) set next iteration values with probability p, see randomUnif
    //computes final quantities for the acceptance rate
    double U = logPostHMC(meanPrior,sigmaPrior,beta_k,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D);
    double propU = logPostHMC(meanPrior,sigmaPrior,beta_star,stats,actor1,actor2,dyad,omit_dyad,interevent_time,model,ordinal,ncores,senderRate,N,C,D);
    double K = 0.5 * arma::accu(r_k.t() * r_k);
    double propK = 0.5 * arma::accu(r_star.t() * r_star); 

    //Accepting or rejection the proposal
    arma::vec randomUnif = arma::randu(1);
    if(randomUnif[0] < exp(U - propU + propK - K)){
      beta_k = beta_star;
      // propU;
    } else {
      beta_k = beta_k;
      // U;
    }
  }
  return L_vec;
}