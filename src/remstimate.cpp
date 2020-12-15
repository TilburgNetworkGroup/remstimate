#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <omp.h>
#include <RcppArmadilloExtensions/sample.h> // maybe this should be removed
#include <typeinfo>
#include <map>
#include <iterator>
#include <string>

// /////////////////////////////////////////////////////////////////////////////////
// ///////////(BEGIN)     remstimateFAST routine functions      (BEGIN)///////////// 
// /////////// to calculate the objects needed for the fast approach   /////////////
// /////////////////////////////////////////////////////////////////////////////////

//' cube2matrix
//'
//' A function to rearrange the cube of statistics into a matrix.
//'
//' @param stats cube structure of dimensions [D*U*M] filled with statistics values. 
//'
//' @return matrix of dimensions [(M*D)*U]
arma::mat cube2matrix(arma::cube stats){
  arma::uword m,d;
  arma::uword new_row = 0;
  arma::mat out(stats.n_rows*stats.n_slices, stats.n_cols, arma::fill::zeros); 
  for(m = 0; m < stats.n_slices; m++){
    for(d = 0; d < stats.n_rows; d++){
      out.row(new_row) = stats.slice(m).row(d);
      new_row += 1;
    }
  }
  return out;
}

//' getUniqueVectors
//'
//' A function to retrieve only the unique vectors of statistics observed throught times points and dyads. This function is based on the result shown by the Appendix C in the paper 'Hierarchical models for relational event sequences', DuBois et al. 2013 (pp. 308-309).
//'
//' @param stats cube structure of dimensions [D*U*M] filled with statistics values.
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
//' @param edgelist is the preprocessed edgelist dataframe with information about [time,actor1,actor2,type,weight] by row.
//' @param risksetCube array of index position fo dyads, with dimensions [N*N*C]
//' @param M number of observed relational events.
//' @param unique_vectors_stats matrix of unique vectors of statistics (output of getUniqueVectors()).
//' @param stats array of statistics with dimensons [D*U*M]
//'
//' @return vector of q's
//'
//' @export
// [[Rcpp::export]]
arma::vec computeOccurrencies(const Rcpp::DataFrame& edgelist, const arma::ucube& risksetCube, const arma::uword& M, const arma::mat& unique_vectors_stats, const arma::cube& stats){

  arma::uword r,m,d;
  arma::uword R = unique_vectors_stats.n_rows;
  arma::uword U = unique_vectors_stats.n_cols;
  Rcpp::IntegerVector actor1 = edgelist["actor1"];
  Rcpp::IntegerVector actor2 = edgelist["actor2"];
  Rcpp::IntegerVector type = edgelist["type"];
  arma::rowvec stats_event_m(U,arma::fill::zeros);
  arma::vec out(R);

  for(m = 0; m < M; m++){
    d = risksetCube(actor1(m),actor2(m),type(m));
    stats_event_m = stats.slice(m).row(d);
    for(r = 0; r < R; r++){
      if(arma::approx_equal(stats_event_m, unique_vectors_stats.row(r), "absdiff", 0.000001)){
        out(r) += 1;
      }
    }
  }
  return out;
}

// /////////////////////////////////////////////////////////////////////////////////
// /////////////(END)     remstimateFAST routine functions      (END)/////////////// 
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
//' @param event_binary is a matrix [M*D] of 1/0/-1 : 1 indicating the observed dyad and 0 (-1) the non observed dyads that could have (have not) occurred.
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param gradient boolean true/false whether to return gradient value
//' @param hessian boolean true/false whether to return hessian value
//'
//' @return list of values: loglik, gradient, hessian
//'
// [[Rcpp::export]]
Rcpp::List remDerivativesStandard(const arma::vec &pars, const arma::cube &stats, const arma::mat &event_binary, const arma::vec &interevent_time,bool gradient = true,bool hessian = true){
    arma::uword U = pars.n_elem; // number of parameters
    arma::uword D = event_binary.n_cols;
    arma::uword M = event_binary.n_rows; // number of events

    arma::uword d,m,l,k;
    arma::vec log_lambda(D,arma::fill::zeros) ;

    double loglik = 0.0;
    arma::mat hess(U,U,arma::fill::zeros);
    arma::vec grad(U,arma::fill::zeros);
   
    for(m = 0; m < M; m++){
        arma::mat stats_m = stats.slice(m).t(); // dimensions : [U*D] we want to access dyads by column
        log_lambda = stats_m.t() * pars;
        for(d = 0; d < D; d++){
            if(event_binary(m,d)!=-1){ // ignoring impossible events that are not in risk set                
                if(event_binary(m,d) == 1){ // if event occured
                    loglik += log_lambda.at(d);
                    grad += stats_m.col(d);
                }               
                double dtelp = exp(log_lambda.at(d))*interevent_time.at(m);  // change `dtelp` name to something else more understandable
                loglik -= dtelp; 
                if(gradient){
                    grad -= stats_m.col(d)*dtelp;
                }             
                if(hessian){
                  for(k = 0; k < U; k++){
                    for (l = k; l < U; l++){
                        hess(k,l) -= stats_m.at(l,d)*stats_m.at(k,d)*dtelp;
                        hess(l,k) = hess(k,l);
                    }
                  }    
                }      
                
            }            
        }
    
    }

    if(gradient && !hessian){
      return Rcpp::List::create(Rcpp::Named("value") = -loglik, Rcpp::Named("gradient") = -grad);
    }else if(!gradient && !hessian){
      return Rcpp::List::create(Rcpp::Named("value") = -loglik);
    }else{
      return Rcpp::List::create(Rcpp::Named("value") = -loglik, Rcpp::Named("gradient") = -grad, Rcpp::Named("hessian") = -hess);
    }
}



//' remDerivativesFast
//'
//' description of the function here
//'
//' @param pars vector of parameters 
//' @param times_r  former m
//' @param occurrencies_r former q
//' @param unique_vectors_stats former U
//' @param gradient boolean
//' @param hessian boolean
//'
//' @return list of value/gradient/hessian in pars
//'
// [[Rcpp::export]]
Rcpp::List remDerivativesFast(const arma::vec& pars, const arma::vec& times_r, const arma::vec& occurrencies_r, const arma::mat& unique_vectors_stats, bool gradient, bool hessian){
  
  arma::uword U = pars.n_elem;
  arma::uword u,k,l;

  double loglik = sum(occurrencies_r.t()*unique_vectors_stats*pars - times_r.t() * exp(unique_vectors_stats*pars));
  
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


  if(gradient && !hessian){
    return Rcpp::List::create(Rcpp::Named("value") = -loglik, Rcpp::Named("gradient") = -grad);
  }else if(!gradient && !hessian){
    return Rcpp::List::create(Rcpp::Named("value") = -loglik);
  }else{
    return Rcpp::List::create(Rcpp::Named("value") = -loglik, Rcpp::Named("gradient") = -grad, Rcpp::Named("hessian") = -hess);
  }

}


//' remDerivatives
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values
//' 
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param event_binary is a matrix [M*D] of 1/0/-1 : 1 indicating the observed dyad and 0 (-1) the non observed dyads that could have (have not) occurred.
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param times_r used in the fast approach
//' @param occurrencies_r used in the fast approach
//' @param unique_vectors_stats used in the fast approach
//' @param fast boolean true/false whether to run the fast approach or not                               
//' @param gradient boolean true/false whether to return gradient value
//' @param hessian boolean true/false whether to return hessian value
//'
//' @return list of values: loglik, gradient, hessian
//'
//' @export
// [[Rcpp::export]]
Rcpp::List remDerivatives(const arma::vec &pars, 
                                  const arma::cube &stats, 
                                  const arma::mat &event_binary, 
                                  const arma::vec &interevent_time,
                                  const arma::vec &times_r, 
                                  const arma::vec &occurrencies_r, 
                                  const arma::mat &unique_vectors_stats,
                                  bool fast = false,
                                  bool gradient = true,
                                  bool hessian = true){
  Rcpp::List out;

  switch (fast)
  {
  case 1: { out = remDerivativesFast(pars,times_r,occurrencies_r,unique_vectors_stats,gradient,hessian);
    break;}

  case 0: { out = remDerivativesStandard(pars,stats,event_binary,interevent_time,gradient,hessian);
    break;}
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
//' @param event_binary rehBinary (inside the reh object)
//' @param interevent_time vector of interevent times (inside the reh object)
//' @param times_r used in the fast approach
//' @param occurrencies_r used in the fast approach
//' @param unique_vectors_stats used in the fast approach
//' @param fast TRUE/FALSE whether to perform the fast approach or not
//' @param epochs number of epochs
//' @param learning_rate learning rate
//' 
//' @return
//'
//' @export
// [[Rcpp::export]]
Rcpp::List GD(const arma::vec &pars,
              const arma::cube &stats, 
              const arma::mat &event_binary,
              const arma::vec &interevent_time,
              const arma::vec &times_r, 
              const arma::vec &occurrencies_r, 
              const arma::mat &unique_vectors_stats,
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
        Rcpp::List derv = remDerivatives(pars_prev,stats,event_binary,interevent_time,times_r,occurrencies_r,unique_vectors_stats,fast,true,false);
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
//' @param event_binary rehBinary (inside the reh object)
//' @param interevent_time vector of interevent times (inside the reh object)
//' @param times_r used in the fast approach
//' @param occurrencies_r used in the fast approach
//' @param unique_vectors_stats used in the fast approach
//' @param fast TRUE/FALSE whether to perform the fast approach or not
//' @param epochs number of epochs
//' @param learning_rate learning rate
//' @param beta1 hyperparameter beta1
//' @param beta2 hyperparameter beta2
//' @param eta hyperparameter eta
//' 
//' @return
//'
//' @export
// [[Rcpp::export]]
Rcpp::List GDADAM(const arma::vec &pars,
                   const arma::cube &stats, 
                   const arma::mat &event_binary,
                   const arma::vec &interevent_time,
                   const arma::vec &times_r, 
                   const arma::vec &occurrencies_r, 
                   const arma::mat &unique_vectors_stats,
                   bool fast = false,
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
        Rcpp::List derv = remDerivatives(pars_prev,stats,event_binary,interevent_time,times_r,occurrencies_r,unique_vectors_stats,fast,true,false);
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
//' @param event_binary is a matrix [M*D] of 1/0/-1 : 1 indicating the observed dyad and 0 (-1) the non observed dyads that could have (have not) occurred.
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param times_r used in the fast approach
//' @param occurrencies_r used in the fast approach
//' @param unique_vectors_stats used in the fast approach
//' @param fast boolean true/false whether to run the fast approach or not                               
//'
//' @return value of log-posterior density
//'
// [[Rcpp::export]]
double logPostHMC(const arma::vec &meanPrior, 
                  const arma::mat &sigmaPrior,
                  const arma::vec &pars, 
                  const arma::cube &stats, 
                  const arma::mat &event_binary, 
                  const arma::vec &interevent_time,
                  const arma::vec &times_r, 
                  const arma::vec &occurrencies_r, 
                  const arma::mat &unique_vectors_stats,
                  bool fast){
  
  Rcpp::List derv = remDerivatives(pars,stats,event_binary,interevent_time,times_r,occurrencies_r,unique_vectors_stats,fast,false,false);

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
//' @param event_binary is a matrix [M*D] of 1/0/-1 : 1 indicating the observed dyad and 0 (-1) the non observed dyads that could have (have not) occurred.
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param times_r used in the fast approach
//' @param occurrencies_r used in the fast approach
//' @param unique_vectors_stats used in the fast approach
//' @param fast boolean true/false whether to run the fast approach or not                               
//'
//' @return value of log-posterior gradient
//'
// [[Rcpp::export]]
arma::vec logPostGradientHMC(const arma::vec &meanPrior,
                              const arma::mat &sigmaPrior, 
                              const arma::vec &pars, 
                              const arma::cube &stats, 
                              const arma::mat &event_binary, 
                              const arma::vec &interevent_time,
                              const arma::vec &times_r, 
                              const arma::vec &occurrencies_r, 
                              const arma::mat &unique_vectors_stats,
                              bool fast){
  
  Rcpp::List derv = remDerivatives(pars,stats,event_binary,interevent_time,times_r,occurrencies_r,unique_vectors_stats,fast,true,false);
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
//' @param event_binary is a matrix [M*D] of 1/0/-1 : 1 indicating the observed dyad and 0 (-1) the non observed dyads that could have (have not) occurred.
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param times_r used in the fast approach
//' @param occurrencies_r used in the fast approach
//' @param unique_vectors_stats used in the fast approach
//' @param fast boolean true/false whether to run the fast approach or not                               
//'
// [[Rcpp::export]]
arma::vec iterHMC(arma::uword L, 
                  double epsilon, 
                  const arma::vec &meanPrior,
                  const arma::mat &sigmaPrior, 
                  const arma::vec &pars, 
                  const arma::cube &stats, 
                  const arma::mat &event_binary, 
                  const arma::vec &interevent_time,
                  const arma::vec &times_r, 
                  const arma::vec &occurrencies_r, 
                  const arma::mat &unique_vectors_stats,
                  bool fast){
  
  arma::vec accept; //vector to store sample
  arma::uword N = pars.size(); //number of parameters
  
  arma::vec r = arma::randn(N); //Rcpp::rnorm(N, 0.0, 1.0); rv's to use in the hamiltonian equations
  arma::vec betaC = pars;
  arma::vec betaP = pars;
  arma::vec rC = r;
  
  //leapfrog algorithm, updates via Hamiltonian equations
  r = r - 0.5 * epsilon * logPostGradientHMC(meanPrior,sigmaPrior,betaP,stats,event_binary,interevent_time,times_r,occurrencies_r,unique_vectors_stats,fast);
  for(arma::uword i = 1; i <= L; i++){
    betaP = betaP + r * epsilon;
    if(i != L) r = r - epsilon * logPostGradientHMC(meanPrior,sigmaPrior,betaP,stats,event_binary,interevent_time,times_r,occurrencies_r,unique_vectors_stats,fast);
  }
  r = r - 0.5 * epsilon * logPostGradientHMC(meanPrior,sigmaPrior,betaP,stats,event_binary,interevent_time,times_r,occurrencies_r,unique_vectors_stats,fast);
  r = -r;
  
  //computes final quantities for the acceptance rate
  double U = logPostHMC(meanPrior,sigmaPrior,betaC,stats,event_binary,interevent_time,times_r,occurrencies_r,unique_vectors_stats,fast);
  double propU = logPostHMC(meanPrior,sigmaPrior,betaP,stats,event_binary,interevent_time,times_r,occurrencies_r,unique_vectors_stats,fast);
  double K = 0.5 * sum(rC.t() * rC);
  double propK = 0.5 * sum(r.t() * r);
  
  //Accepting or rejection the proposal
  arma::vec randomUnif = arma::randu(1);
  if(randomUnif[0] < exp((-propU - propK)/(-U - K))){ //R::runif(0,1)
    accept = betaP;
  } else {
    accept = betaC;
  }
  return accept;
}       


//' burninHMC (to check whether this function experiences issues with the definition of int rows and the following codings)
//'
//' This function performs the burn-in and the thinning at the end of the HMC
//'
//' @param samples cube with final draws
//' @param n_burnin is the number of draws to discard after running the chains
//' @param n_thin is the number of draws to be skipped. For instance, if n_thin = 10, draws will be selected every 10 generated draws: 1, 11, 21, 31, ...
//'                          
//' @return cube with selected draws
//'
// [[Rcpp::export]]
arma::cube burninHMC(const arma::cube& samples, arma::uword n_burnin, arma::uword n_thin = 1){
  
  arma::uword rows = round((samples.n_rows - n_burnin)/n_thin); //number of rows of output
  arma::cube final(rows, samples.n_cols, samples.n_slices); //output
  
  for(arma::uword i = 0; i < final.n_slices; i++){
    
    arma::uword num = n_burnin;
    
    for(arma::uword j = 0; j < final.n_rows; j++){
      
      final.slice(i).row(j) = samples.slice(i).row(num);
      
      num += n_thin;
      
    }
    
  }

  return final;
}



//' HMC 
//'
//' This function performs the Hamiltonian Monte Carlo
//'
//' @param pars_init is a matrix of dimensions U x n_chains where for each column (chain) a random vector of initial values for the parameter is supplied.
//' @param n_iters is the number of samples from the posterior that have to be generated.
//' @param n_chains number of chains of length n_iters
//' @param n_burnin is the number of draws to discard after running the chains
//' @param meanPrior is a vector of prior means with the same dimension as the vector of parameters
//' @param sigmaPrior is a matrix, I have been using a diagonal matrix here with the same dimension as the vector os parameters
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param event_binary is a matrix [M*D] of 1/0/-1 : 1 indicating the observed dyad and 0 (-1) the non observed dyads that could have (have not) occurred.
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param times_r used in the fast approach
//' @param occurrencies_r used in the fast approach
//' @param unique_vectors_stats used in the fast approach
//' @param fast boolean TRUE/FALSE whether to run the fast approach or not (default = FALSE) 
//' @param n_thin is the number of draws to be skipped. For instance, if n_thin = 10, draws will be selected every 10 generated draws: 1, 11, 21, 31, ...
//' @param L number of leapfrogs. Default (and recommended) value is 100.
//' @param epsilon size of the leapfrog. Default value is 1e-02.
//' @param n_threads number of threads for parallel computing (default = 1)
//'                          
//' @return posterior draws
//'
// [[Rcpp::export]]
arma::cube HMC(arma::mat pars_init, 
                arma::uword n_iters, 
                arma::uword n_chains, 
                arma::uword n_burnin, 
                const arma::vec& meanPrior,
                const arma::mat& sigmaPrior,
                const arma::cube &stats, 
                const arma::mat &event_binary, 
                const arma::vec &interevent_time,
                const arma::vec &times_r, 
                const arma::vec &occurrencies_r, 
                const arma::mat &unique_vectors_stats,
                bool fast = false,
                arma::uword n_thin = 1,
                arma::uword L = 100, 
                double epsilon = 0.01,
                int n_threads = 1){
  
  arma::cube store(n_iters, pars_init.n_rows, n_chains); //output
  arma::uword j,i;

  //omp_set_dynamic(0);         // disabling dynamic teams
  //omp_set_num_threads(n_threads); // number of threads for all consecutive parallel regions
  //#pragma omp parallel for private(j,i) shared(n_iters,n_chains,meanPrior,sigmaPrior,fast,L,epsilon,store,pars_init,stats,event_binary,interevent_time,times_r,occurrencies_r,unique_vectors_stats)
  for(j = 0; j < n_chains; j++){ //looping through chains
    
    arma::mat aux(pars_init.n_rows, 1,arma::fill::zeros);
    arma::mat chain_j(n_iters,pars_init.n_rows,arma::fill::zeros);

    for(i = 0; i < n_iters; i++){ //looping through iterations of the MCMC
      
      if(i == 0){
        
        //this step only get the first sample out of the starting value
        aux.col(0) = iterHMC(L,epsilon,meanPrior,sigmaPrior,pars_init.col(j),stats,event_binary,interevent_time,times_r,occurrencies_r,unique_vectors_stats,fast);
        
        chain_j.row(i) = aux.col(0).t();
        
        continue;
        
      } else {
        
        //Then the next step will always be based on the previous one
        aux.col(0) = iterHMC(L,epsilon,meanPrior,sigmaPrior,aux.col(0),stats,event_binary,interevent_time,times_r,occurrencies_r,unique_vectors_stats,fast);
        
        chain_j.row(i) = aux.col(0).t();
        
      }
    }
    store.slice(j) = chain_j;
    
  }
  
  //this does the burn-in and thinning
  arma::cube out_cube = burninHMC(store,n_burnin,n_thin);
  //arma::mat out_mat = cube2matrix(out_cube); 
  
  return out_cube;
}


// /////////////////////////////////////////////////////////////////////////////////
// ///////////(END)            Hamiltonian Monte Carlo            (END)/////////////
// /////////////////////////////////////////////////////////////////////////////////



//' remDerivativesStandardParallel
//'
//' function that returns a list as an output with loglikelihood/gradient/hessian values at specific parameters' values
//' 
//' @param pars is a vector of parameters (note: the order must be aligned with the column order in 'stats')
//' @param stats is cube of M slices. Each slice is a matrix of dimensions D*U with statistics of interest by column and dyads by row.
//' @param event_binary is a matrix [M*D] of 1/0/-1 : 1 indicating the observed dyad and 0 (-1) the non observed dyads that could have (have not) occurred.
//' @param interevent_time the time difference between the current time point and the previous event time.
//' @param gradient boolean true/false whether to return gradient value
//' @param hessian boolean true/false whether to return hessian value
//' @param n_threads integer
//'
//' @return list of values: loglik, gradient, hessian
//'
//' @export
// [[Rcpp::export]]
Rcpp::List remDerivativesStandardParallel(const arma::vec &pars, const arma::cube &stats, const arma::mat &event_binary, const arma::vec &interevent_time,bool gradient = true,bool hessian = true, int n_threads = 1){
    arma::uword U = pars.n_elem; // number of parameters
    arma::uword D = event_binary.n_cols;
    arma::uword M = event_binary.n_rows; // number of events

    arma::uword d,m,l,k;
    arma::vec log_lambda(D,arma::fill::zeros) ;

    arma::vec loglik(M,arma::fill::zeros);
    arma::cube hess(U,U,M,arma::fill::zeros);
    arma::mat grad(U,M,arma::fill::zeros);
   
    omp_set_dynamic(0);         // disabling dynamic teams
    omp_set_num_threads(n_threads); // number of threads for all consecutive parallel regions
    #pragma omp parallel for private(m,d,l,k,log_lambda) shared(M,D,U,loglik,hess,grad,stats,event_binary,gradient,hessian,interevent_time)
    for(m = 0; m < M; m++){
        arma::mat stats_m = stats.slice(m).t(); // dimensions : [U*D] we want to access dyads by column
        log_lambda = stats_m.t() * pars;
        for(d = 0; d < D; d++){
            if(event_binary(m,d)!=-1){ // ignoring impossible events that are not in risk set                
                if(event_binary(m,d) == 1){ // if event occured
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

    if(gradient && !hessian){
      return Rcpp::List::create(Rcpp::Named("value") = -sum(loglik), Rcpp::Named("gradient") = -sum(grad,1));
    }else if(!gradient && !hessian){
      return Rcpp::List::create(Rcpp::Named("value") = -sum(loglik));
    }else{
      arma::cube H = -sum(hess,2);
      return Rcpp::List::create(Rcpp::Named("value") = -sum(loglik), Rcpp::Named("gradient") = -sum(grad,1), Rcpp::Named("hessian") = H.slice(0));
    }
}




//' generate_mvt
//'
//' @param input
//'
//' @return matrix
//'
//' @export
// [[Rcpp::export]]
arma::mat tryFunction(arma::cube &input){
  
  // method 1: 
  arma::uword rows_0 = input.n_rows*input.n_cols;
input.reshape(rows_0, input.n_slices, 1);
arma::mat C = input.slice(0);
// method 2:
//arma::mat C = arma::reshape( arma::mat(input.memptr(), input.n_elem, 1, false), 5*4, 3);
//
return C;
}