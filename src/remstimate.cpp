#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <omp.h>
#include <RcppArmadilloExtensions/sample.h> // used for the sample function inside RcppArmadillo
#include <typeinfo>
#include <map>
#include <iterator>
#include <string>



// /////////////////////////////////////////////////////////////////////////////////
// //////(START)              getDyad-Index/Composition               (START)///////
// /////////////////////////////////////////////////////////////////////////////////


// getDyadIndex - function that in the future version of the pkg will be imported from remify.h
//
// @param actor1 id of actor1 from 0 to N-1
// @param actor2 id of actor2 from 0 to N-1
// @param type id of event type from 0 to C-1
// @param N number of actors
// @param directed bool FALSE/TRUE if the networks is directed (TRUE) or not (FALSE)
//
// @return dyad index according to the combination of id's of actor1/actor2/type
int getDyadIndex(double actor1, double actor2, double type, int N, bool directed) {

    int dyad = -999; // returning impossible index if the dyad is a self-edge (i.e., sender and receiver are the same actor)
    if(actor1 != actor2){
        if(!directed){ // when directed == FALSE we sort actor1 and actor2
            if(actor1 < actor2){
                int dyad_loc = (N*(N-1)/2)*type+(N-1)*actor1+actor2-actor1-1-(actor1*actor1)/2;
                if(actor1>0){
                    dyad_loc += actor1/2;
                }
                dyad = dyad_loc;
            }
            else{
                int dyad_loc = (N*(N-1)/2)*type+(N-1)*actor2+actor1-actor2-1-(actor2*actor2)/2;
                if(actor2>0){
                    dyad_loc += actor2/2;
                }
                dyad = dyad_loc;
            }
        }
        else{ 
            // when directed == TRUE (we do not sort) (actor1 = sender, actor2 = receiver)
            int dyad_loc = N*(N-1)*type+(N-1)*actor1+actor2;
            if(actor2>actor1){
                dyad_loc -= 1;
            }
            dyad = dyad_loc;
        }
    }
    return dyad;
}


// getDyadComposition (only for directed for now) - function that in the future version of the pkg will be imported from remify.h
//
// @param d id of the dyad
// @param C number of event types
// @param N number of actors
// @param D number of dyads
//
// @return dyad index according to the combination of id's of actor1/actor2/type
Rcpp::IntegerVector getDyadComposition(int d, int C, int N, int D) {
  Rcpp::IntegerVector composition(3);
  // Note :
  // (1) this function assumes that all the possible dyads are in the stats object
  // (2) this function is not coded to account for reduced (that omits dyads) arrays of stats
  // (3) this function works only for directed netwroks [[will be updated in the future to the undirected case]]
  double r = d; // this will be finally the receiver
  r += 1;
  int sender,receiver,type = -999;
  double c = 1, s = 1;
  while(c<=C){
    if((r/D)<=(c/C)){
      type = (c-1);
      break;
    }
    c += 1;
  }

  //if(type == (-999)){
  //  Rcpp::Rcout << "error \n"; //errorMessage(0); //adjust error message
  //}

  r -= N*(N-1)*type;

  while(s<=N){
    if((r/(N*(N-1)))<=(s/N)){
      sender = (s-1);
      break;
    }
    s += 1;
  }

  //if(sender == (-999)){
  //  Rcpp::Rcout << "error \n"; //errorMessage(0); //adjust error message
  //}

  arma::mat receiver_vec(N,1);
  receiver_vec.col(0) = arma::linspace(0,N-1,N);
  receiver_vec.shed_row(sender);
  r -= (N-1)*sender;
  receiver = receiver_vec[r-1]; // if either type or sender are not found, the function will stop earlier
  composition = {sender,receiver,type};
  return composition;
}


// /////////////////////////////////////////////////////////////////////////////////
// ////////(END)              getDyad-Index/Composition               (END)/////////
// /////////////////////////////////////////////////////////////////////////////////



// /////////////////////////////////////////////////////////////////////////////////
// ///////////(BEGIN)              remDerivatives               (BEGIN)/////////////
// /////////////////////////////////////////////////////////////////////////////////


//' remDerivativesStandard
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values
//'
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param edgelist is a matrix [M*3] of [time/dyad/weight]
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param ordinal boolean that indicate whether to use the ordinal or interval timing likelihood
//' @param ncores integer referring to the number of threads for the parallelization
//' @param gradient boolean true/false whether to return gradient value
//' @param hessian boolean true/false whether to return hessian value
//'
//' @return list of values: loglik, gradient, hessian
//'
//' @export
// [[Rcpp::export]]
Rcpp::List remDerivativesStandard(const arma::vec &pars,
                                          const arma::cube &stats,
                                          const arma::mat &edgelist, 
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
      #pragma omp parallel for private(m,d,l,k,log_lambda) shared(M,D,U,loglik,hess,grad,stats,riskset_time_vec,riskset_mat,edgelist,gradient,hessian,interevent_time)
      for(m = 0; m < M; m++){
          arma::mat stats_m = stats.slice(m).t(); // dimensions : [U*D] we want to access dyads by column
          log_lambda = stats_m.t() * pars;
          int riskset_time_m = riskset_time_vec(m);
          if(riskset_time_m!=(-1)){ // if the riskset_time_vec[m]=riskset_time_m is different from, then a dynamic riskset is observed
            for(d = 0; d < D; d++){
              if(riskset_mat(riskset_time_m,d)){ // ignoring impossible events that are not in risk set // [[try without ==1 ]]
                  if(edgelist(m,1) == d){ // if the event occurred
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
              if(edgelist(m,1) == d){ // if event occured
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
      #pragma omp parallel for private(m,d,l,k,log_lambda) shared(M,D,U,loglik,hess,grad,stats,riskset_time_vec,riskset_mat,edgelist,gradient,hessian,interevent_time)
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
              if(edgelist(m,1) == d){ // if event occured
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
          if(edgelist(m,1) == d){ // if event occured
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
//' @param edgelist is a matrix [M*3] of [time/dyad/weight]
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param C number of event types 
//' @param D number of dyads
//' @param ordinal boolean, true if the likelihood to use is the ordinal one, interval otherwise
//'
//'
//' @return list of values: loglik, grad, fisher information
//'
// [[Rcpp::export]]
Rcpp::List remDerivativesSenderRates(
        const arma::vec &pars,
        const arma::cube &stats,
        const arma::mat &edgelist, 
        const Rcpp::List &omit_dyad,
        const arma::vec &interevent_time,
        int C,
        int D,
        bool ordinal = false){

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
      dyad_m = getDyadComposition(edgelist(m,1),C,N,D);
      int sender = dyad_m[0];
      // int type = dyad[2]; // when type will be integrated in the function

      //reset internal variables
      fisher_m.zeros();
      expected_stat_m.zeros();

      // riskset_m
      int riskset_time_m = riskset_time_vec(m);

      loglik += std::log(lambda_s(sender));
      grad += stats_m.col(sender);

      // changes in the riskset at m-th event
      if(riskset_time_m!=(-1)){
        sum_lambda = sum(riskset_mat.row(riskset_time_m).t() % lambda_s); 
        if(ordinal){
          loglik-= std::log(sum_lambda);
          for(n = 0; n < N; n++){         //loop throughout all actors
            if(riskset_mat(riskset_time_m,n) == 1){
              expected_stat_m += lambda_s(n) * (stats_m.col(n));
              expected_stat_m /= sum_lambda; //exp(params_s * X_i)*X_i / sum_h (exp(params * X_h))
              grad -= expected_stat_m;
              fisher_m += (expected_stat_m * expected_stat_m.t());
              fisher_m -= (lambda_s(n) / sum_lambda) *( stats_m.col(n) * (stats_m.col(n).t()));
              fisher += fisher_m;
            }
          }
        } 
        else{
          loglik -=  sum_lambda * interevent_time(m);
          grad -= interevent_time(m) * stats_m * (riskset_mat.row(riskset_time_m).t() % lambda_s); 
          for(n = 0; n < N; n++){         //loop throughout all actors
            fisher_m += (riskset_mat(riskset_time_m,n) * lambda_s(n)) * (stats_m.col(n) * (stats_m.col(n).t()) );
          }
          fisher -= interevent_time(m) * fisher_m;
        }
      }
      else{  //no dynamic riskset
        sum_lambda = sum(lambda_s);
        if(ordinal){
          loglik-= std::log(sum_lambda);
          for(n = 0; n < N; n++){         //loop throughout all actors
            expected_stat_m += lambda_s(n) * (stats_m.col(n));
            expected_stat_m /= sum_lambda; //exp(params_s * X_i)*X_i / sum_h (exp(params * X_h))
            grad -= expected_stat_m;
            fisher_m += (expected_stat_m * expected_stat_m.t());
            fisher_m -= (lambda_s(n) / sum_lambda) * (stats_m.col(n) * (stats_m.col(n).t()));
            fisher += fisher_m;
          }
        } 
        else{
          loglik -=  sum_lambda * interevent_time(m);
          grad -= interevent_time(m) * stats_m * lambda_s; //(6)
          for(n = 0; n < N; n++){         //loop throughout all actors
            fisher_m += lambda_s(n)*(stats_m.col(n) * (stats_m.col(n).t()));
          }
          fisher -= interevent_time(m) * fisher_m;
        }
      }
    }
    return Rcpp::List::create(Rcpp::Named("value") = -loglik, Rcpp::Named("gradient") = -grad, Rcpp::Named("hessian") = -fisher);
}


//' remDerivativesReceiverChoice
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values for estimating the receiver choice parameters for the actor oriented model
//'
//' @param pars is a vector of parameters
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param edgelist, output from remify, (note: indices of the actors must start from 0)
//' @param omit_dyad, list object that takes care of the dynamic rikset (if defined)
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param directed, boolean TRUE/FALSE if the network is directed or not
//' @param N the number of actors
//' @param C number of event types 
//' @param D number of dyads
//'
//'
//' @return list of values: loglik, grad, fisher
//'
//[[Rcpp::export]]
Rcpp::List remDerivativesReceiverChoice(
        const arma::vec &pars,
        const arma::cube &stats,
        const arma::mat &edgelist,
        const Rcpp::List &omit_dyad,
        const arma::vec &interevent_time,
        bool directed,
        int N,
        int C,
        int D){

    arma::uword P_d = stats.n_cols; // number of parameters for dyad stats
    arma::uword M = edgelist.n_rows; // number of events
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

      //actors in remify edgelist are 0 indexed hence no -1
      dyad_m = getDyadComposition(edgelist(m,1),C,N,D);
      int sender = dyad_m[0]; 
      int receiver = dyad_m[1];
      //int type = dyad_m[2]; // when type will be integrated in the function

      //reset internal variables
      fisher_m.zeros();
      expected_stat_m.zeros();
      denom = 0;

      // riskset_m
      int riskset_time_m = riskset_time_vec(m);

      // changes in the riskset at m-th event
      if(riskset_time_m!=(-1)){ 
        for(n = 0; n<N; n++){
            dyad = getDyadIndex(sender,n,0,N,directed);
            if(n!=sender && riskset_mat(riskset_time_m,dyad) == 1){ // dynamic riskset
                //loglik
                denom += lambda_d(n); // exp(param_d * X_sender_i) (4)
                fisher_m -= lambda_d(n)*(stats_m.col(n) * stats_m.col(n).t());
                expected_stat_m += lambda_d(n) * stats_m.col(n);
            }
        }
      }
      else{ // no changes in the riskset at m-th event
        for(n = 0;n<N; n++){ //loop throught all actors
          if(n!=sender){ // 
              //loglik
              denom += lambda_d(n); // exp(param_d * X_sender_i) (4)
              fisher_m -= lambda_d(n)*(stats_m.col(n) * stats_m.col(n).t());
              expected_stat_m += lambda_d(n) * stats_m.col(n);
          }
        }
      }

      loglik +=  std::log(lambda_d(receiver)); //(3) params_d^T * X_sender_receiver
      loglik -= std::log(denom);

      grad += stats_m.col(receiver); //(7)
      expected_stat_m /= denom; //(8)
      grad -= expected_stat_m;

      fisher_m /= denom;
      fisher_m += (expected_stat_m * expected_stat_m.t());
      fisher += fisher_m;      
    }

    return Rcpp::List::create(Rcpp::Named("value") = -loglik, Rcpp::Named("gradient") = -grad, Rcpp::Named("hessian") = -fisher);
}


//' remDerivatives 
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values
//'
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param edgelist is a matrix [M*3] of [time/dyad/weight]
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param model either "actor" or "tie" model
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param ncores number of threads to use for the parallelization
//' @param gradient boolean true/false whether to return gradient value
//' @param hessian boolean true/false whether to return hessian value
//' @param senderRate boolean true/false (it is used only when model = "actor") indicates if to estimate the senderRate model (true) or the ReceiverChoice model (false)
//' @param N number of actors. This argument is used only in the ReceiverChoice likelihood (model = "actor")
//' @param C number of event types 
//' @param D number of dyads
//'
//' @return list of values: loglik, gradient, hessian
//'
//' @export
// [[Rcpp::export]]
Rcpp::List remDerivatives(const arma::vec &pars,
                                  const arma::cube &stats,
                                  const arma::mat &edgelist,
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
  // each routine should include a parameter 'fast' if the fast version can be useb
  std::vector<std::string> models = {"tie","actor"};

  std::vector<std::string>::iterator itr = std::find(models.begin(), models.end(), model);
  auto which_model = std::distance(models.begin(), itr);

  switch (which_model)
  {
  case 0: { out = remDerivativesStandard(pars,stats,edgelist,omit_dyad,interevent_time,ordinal,ncores,gradient,hessian);
    break;}

  case 1: { 
      switch (senderRate){
        case 0 : {
                  out = remDerivativesReceiverChoice(pars,stats,edgelist,omit_dyad,interevent_time,true,Rcpp::as<int>(N),Rcpp::as<int>(C),Rcpp::as<int>(D));  // add gradient and hessian switch
                  break;
                  }
        case 1 : {
                  out = remDerivativesSenderRates(pars,stats,edgelist,omit_dyad,interevent_time,Rcpp::as<int>(C),Rcpp::as<int>(D),ordinal); // add gradient and hessian switch
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


//' GD
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values
//'
//' @param pars parameters
//' @param stats array of statistics
//' @param edgelist is a matrix [M*3] of [time/dyad/weight]
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined
//' @param interevent_time vector of interevent times (inside the reh object)
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param model either "actor" or "tie" model
//' @param ncores number of threads to use for the parallelization
//' @param fast TRUE/FALSE whether to perform the fast approach or not
//' @param epochs number of epochs
//' @param learning_rate learning rate
//'
//' @return optimization with Gradient Descent algorithm
//'
//' @export
// [[Rcpp::export]]
Rcpp::List GD(const arma::vec &pars,
              const arma::cube &stats,
              const arma::mat &edgelist,
              const Rcpp::List &omit_dyad,
              const arma::vec &interevent_time,
              std::string model,
              bool ordinal = false,
              int ncores = 1,
              bool fast = false,
              int epochs = 200,
              double learning_rate = 0.001){
    //if loss function oscillates in later epochs, then reduce learning rate for better convergence
    // suggestions: (1) tuning of learning rate, (2) standardization of data
    // further implementation: parallelizing with multiple starting points
    // stopping rule (gradient value/loglik difference/number of epochs)
    arma::vec pars_prev = pars;
    arma::vec pars_next;
    arma::uword P = pars.n_elem;

    //for output
    arma::mat beta(P,epochs,arma::fill::zeros);
    arma::vec loss(epochs,arma::fill::zeros);

    for(int i =0; i<epochs;i++){
        Rcpp::List derv = remDerivatives(pars_prev,stats,edgelist,omit_dyad,interevent_time,model,ordinal,ncores,true,false);
        double loglik = derv["value"];
        arma::vec grad = derv["gradient"];
        pars_next = pars_prev - (learning_rate * grad);
        beta.col(i) = pars_next;
        loss(i) = loglik;
        pars_prev = pars_next;
    }
    return Rcpp::List::create(Rcpp::Named("loss") = loss, Rcpp::Named("coef") = pars_next, Rcpp::Named("betas") = beta);
}


//' GDADAM
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values
//'
//' @param pars parameters
//' @param stats array of statistics
//' @param edgelist is a matrix [M*3] of [time/dyad/weight]
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined
//' @param interevent_time vector of interevent times (inside the reh object)
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param model either "actor" or "tie" model
//' @param ncores number of threads to use for the parallelization
//' @param fast TRUE/FALSE whether to perform the fast approach or not
//' @param epochs number of epochs
//' @param learning_rate learning rate
//' @param beta1 hyperparameter beta1
//' @param beta2 hyperparameter beta2
//' @param eta hyperparameter eta
//'
//' @return optimization with GDADAM
//'
//' @export
// [[Rcpp::export]]
Rcpp::List GDADAM(const arma::vec &pars,
                   const arma::cube &stats,
                   const arma::mat &edgelist,
                   const Rcpp::List &omit_dyad,
                   const arma::vec &interevent_time,
                   std::string model,
                   bool ordinal = false,
                   int ncores = 1,
                   int epochs = 200,
                   double learning_rate = 0.02,
                   double beta1 = 0.9,
                   double beta2 = 0.999,
                   double eta = 0.00000001){

// default values for hyper-parameters recommended in literature (https://arxiv.org/abs/1412.6980)
// β1 =  0.9 for β2 =  0.999
// eta = 10^-8
    arma::uword P = pars.n_elem;
    arma::vec pars_prev = pars;
    arma::vec pars_next;


    arma::vec moment1(P,arma::fill::zeros);
    arma::vec moment2(P,arma::fill::zeros);

    arma::vec moment1_est(P,arma::fill::zeros);
    arma::vec moment2_est(P,arma::fill::zeros);

    //for output
    arma::mat beta(P,epochs,arma::fill::zeros);
    arma::vec loss(epochs,arma::fill::zeros);

    for(int i =0; i<epochs;i++){
        Rcpp::List derv = remDerivatives(pars_prev,stats,edgelist,omit_dyad,interevent_time,model,ordinal,ncores,true,false);
        double loglik = derv["value"];

        arma::vec grad = derv["gradient"];
        moment1 = (beta1*moment1) + (1-beta1)*grad;
        moment2 = (beta2 * moment2) + (1-beta2) * (grad%grad);
        // bias-corrected first and second moment estimates:
        moment1_est = moment1 / (1-pow(beta1,i+1));
        moment2_est = moment2/  (1-pow(beta2,i+1));
        pars_next = pars_prev - (learning_rate * moment1_est)/(sqrt(moment2_est)+eta);

        beta.col(i) = pars_next;
        loss(i) = loglik;

        pars_prev = pars_next;
    }
    return Rcpp::List::create(Rcpp::Named("loss") = loss, Rcpp::Named("coef") = pars_next, Rcpp::Named("betas") = beta);
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
//' @param edgelist is a matrix [M*3] of [time/dyad/weight]
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param model either "actor" or "tie" model
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param ncores number of threads to use for the parallelization
//' @param fast boolean true/false whether to run the fast approach or not
//'
//' @return value of log-posterior density
//'
// [[Rcpp::export]]
double logPostHMC(const arma::vec &meanPrior,
                  const arma::mat &sigmaPrior,
                  const arma::vec &pars,
                  const arma::cube &stats,
                  const arma::mat &edgelist,
                  const Rcpp::List &omit_dyad,
                  const arma::vec &interevent_time,
                  std::string model,
                  bool ordinal = false,
                  int ncores = 1,
                  bool fast = false){

  Rcpp::List derv = remDerivatives(pars,stats,edgelist,omit_dyad,interevent_time,model,ordinal,ncores,fast,false,false);

  double prior = - sum(0.5 * (pars.t() - meanPrior.t()) * inv(sigmaPrior) * (pars - meanPrior));

  return -(prior + derv[0]);
}


//' logPostGradientHMC
//'
//' This function calculates the value of the gradient of the log-posterior density
//'
//' @param meanPrior is a vector of prior means with the same dimension as the vector of parameters
//' @param sigmaPrior is a matrix, I have been using a diagonal matrix here with the same dimension as the vector os parameters
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param edgelist is a matrix [M*3] of [time/dyad/weight]
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param model either "actor" or "tie" model
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param ncores number of threads to use for the parallelization
//' @param fast boolean true/false whether to run the fast approach or not
//'
//' @return value of log-posterior gradient
//'
// [[Rcpp::export]]
arma::vec logPostGradientHMC(const arma::vec &meanPrior,
                              const arma::mat &sigmaPrior,
                              const arma::vec &pars,
                              const arma::cube &stats,
                              const arma::mat &edgelist,
                              const Rcpp::List &omit_dyad,
                              const arma::vec &interevent_time,
                              std::string model,
                              bool ordinal = false,
                              int ncores = 1,
                              bool fast = false){

  Rcpp::List derv = remDerivatives(pars,stats,edgelist,omit_dyad,interevent_time,model,ordinal,ncores,fast,true,false);
  arma::vec gprior = - 0.5 * inv(sigmaPrior) * (pars - meanPrior);
  arma::vec glp = derv[1];

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
//' @param edgelist is a matrix [M*3] of [time/dyad/weight]
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined//' @param interevent_time the time difference between the current time point and the previous event time.//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param model either "actor" or "tie" model
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param ncores number of threads to use for the parallelization
//' @param fast boolean true/false whether to run the fast approach or not
//'
// [[Rcpp::export]]
arma::vec iterHMC(arma::uword L,
                  double epsilon,
                  const arma::vec &meanPrior,
                  const arma::mat &sigmaPrior,
                  const arma::vec &pars,
                  const arma::cube &stats,
                  const arma::mat &edgelist,
                  const Rcpp::List &omit_dyad,
                  const arma::vec &interevent_time,
                  std::string model,
                  bool ordinal = false,
                  int ncores = 1,
                  bool fast = false){

  arma::vec accept; //vector to store sample
  arma::uword N = pars.size(); //number of parameters

  arma::vec r = arma::randn(N); //Rcpp::rnorm(N, 0.0, 1.0); rv's to use in the hamiltonian equations
  arma::vec betaC = pars;
  arma::vec betaP = pars;
  arma::vec rC = r;

  //leapfrog algorithm, updates via Hamiltonian equations
  r = r - 0.5 * epsilon * logPostGradientHMC(meanPrior,sigmaPrior,betaP,stats,edgelist,omit_dyad,interevent_time,model,ordinal,ncores,fast);
  for(arma::uword i = 1; i <= L; i++){
    betaP = betaP + r * epsilon;
    if(i != L) r = r - epsilon * logPostGradientHMC(meanPrior,sigmaPrior,betaP,stats,edgelist,omit_dyad,interevent_time,model,ordinal,ncores,fast);
  }
  r = r - 0.5 * epsilon * logPostGradientHMC(meanPrior,sigmaPrior,betaP,stats,edgelist,omit_dyad,interevent_time,model,ordinal,ncores,fast);
  r = -r;

  //computes final quantities for the acceptance rate
  double U = logPostHMC(meanPrior,sigmaPrior,betaC,stats,edgelist,omit_dyad,interevent_time,model,ordinal,ncores,fast);
  double propU = logPostHMC(meanPrior,sigmaPrior,betaP,stats,edgelist,omit_dyad,interevent_time,model,ordinal,ncores,fast);
  double K = 0.5 * sum(rC.t() * rC);
  double propK = 0.5 * sum(r.t() * r);

  //Accepting or rejection the proposal
  arma::vec randomUnif = arma::randu(1);
  if(randomUnif[0] < exp((-propU - propK)/(-U - K))){
    accept = betaP;
    // propU;
  } else {
    accept = betaC;
    // U;
  }
  return accept;
}


//' burninHMC (to check whether this function experiences issues with the definition of int rows and the following codings)
//'
//' This function performs the burn-in and the thinning at the end of the HMC
//'
//' @param samples cube with final draws
//' @param burnin is the number of draws to discard after running the chains
//' @param thin is the number of draws to be skipped. For instance, if thin = 10, draws will be selected every 10 generated draws: 1, 11, 21, 31, ...
//'
//' @return cube with selected draws
//'
// [[Rcpp::export]]
arma::cube burninHMC(const arma::cube& samples, arma::uword burnin, arma::uword thin = 1){

  arma::uword rows = round((samples.n_rows - burnin)/thin); //number of rows of output
  arma::cube out_cube(rows, samples.n_cols, samples.n_slices); //output

  for(arma::uword i = 0; i < out_cube.n_slices; i++){

    arma::uword num = burnin;

    for(arma::uword j = 0; j < out_cube.n_rows; j++){

      out_cube.slice(i).row(j) = samples.slice(i).row(num);

      num += thin;

    }

  }
  return out_cube;
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
//' @param edgelist is a matrix [M*3] of [time/dyad/weight]
//' @param omit_dyad is a list of two objects: vector "time" and matrix "riskset". Two object for handling changing risksets. NULL if no change is defined//' @param interevent_time the time difference between the current time point and the previous event time.//' @param 
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param ordinal whether to use(TRUE) the ordinal likelihood or not (FALSE) then using the interval likelihood
//' @param model either "actor" or "tie" model
//' @param ordinal logic TRUE/FALSE
//' @param ncores number of threads to use for the parallelization
//' @param fast boolean TRUE/FALSE whether to run the fast approach or not (default = FALSE)
//' @param thin is the number of draws to be skipped. For instance, if thin = 10, draws will be selected every 10 generated draws: 1, 11, 21, 31, ...
//' @param L number of leapfrogs. Default (and recommended) value is 100.
//' @param epsilon size of the leapfrog. Default value is 1e-02.
//' @param ncores number of threads for parallel computing (default = 1)
//'
//' @return posterior draws
//'
// [[Rcpp::export]]
arma::cube HMC(arma::mat pars_init,
                arma::uword nsim,
                arma::uword nchains,
                arma::uword burnin,
                const arma::vec& meanPrior,
                const arma::mat& sigmaPrior,
                const arma::cube &stats,
                const arma::mat &edgelist,
                const Rcpp::List &omit_dyad,
                const arma::vec &interevent_time,
                std::string model,
                bool ordinal = false,
                int ncores = 1,
                bool fast = false,
                arma::uword thin = 1,
                arma::uword L = 100,
                double epsilon = 0.01){

  arma::cube store(nsim, pars_init.n_rows, nchains); //output
  arma::uword j,i;

  for(j = 0; j < nchains; j++){ //looping through chains

    arma::mat aux(pars_init.n_rows, 1,arma::fill::zeros);
    arma::mat chain_j(nsim,pars_init.n_rows,arma::fill::zeros);

    for(i = 0; i < nsim; i++){ //looping through iterations of the MCMC
      Rcpp::Rcout << i << "\n";
      if(i == 0){

        //this step only get the first sample out of the starting value
        aux.col(0) = iterHMC(L,epsilon,meanPrior,sigmaPrior,pars_init.col(j),stats,edgelist,omit_dyad,interevent_time,model,ordinal,ncores,fast);

        chain_j.row(i) = aux.col(0).t();

        continue;

      } else {

        //Then the next step will always be based on the previous one
        aux.col(0) = iterHMC(L,epsilon,meanPrior,sigmaPrior,aux.col(0),stats,edgelist,omit_dyad,interevent_time,model,ordinal,ncores,fast);

        chain_j.row(i) = aux.col(0).t();

      }
    }
    store.slice(j) = chain_j;

  }

  //this does the burn-in and thinning
  arma::cube draws_cube = burninHMC(store,burnin,thin); // it would be ideal to make burninHMC return a matrix
 // arma::mat out_mat = cube2matrix(draws_cube);

  return draws_cube; //out_mat;
}


// /////////////////////////////////////////////////////////////////////////////////
// ///////////(END)            Hamiltonian Monte Carlo            (END)/////////////
// /////////////////////////////////////////////////////////////////////////////////



// /////////////////////////////////////////////////////////////////////////////////
// ///////////(START)           Experimental function           (START)/////////////
// /////////////////////////////////////////////////////////////////////////////////


//' experimental_function (where to try out specific operations at C++ level)
//'
//' the experimental function has no description
//'
//' @param x integer value
//'
//' @return matrix
//'
//' @export
// [[Rcpp::export]]
arma::mat experimental_function(const arma::uword &x){
  arma::mat out(x,x,arma::fill::zeros);
  return out;
}


// /////////////////////////////////////////////////////////////////////////////////
// ///////////(END)           Experimental function               (END)/////////////
// /////////////////////////////////////////////////////////////////////////////////
