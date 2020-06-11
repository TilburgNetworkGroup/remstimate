#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <omp.h>
#include <RcppArmadilloExtensions/sample.h>
#include <typeinfo>
#include <map>
#include <iterator>
#include <string>


//' useful function that acts like print() in R
#define LOG(x) std::cout << x << "\n"

//' reh (a function to pre-process the REH input)
//'
//' description of remstimate::reh here (most likeliy the function won't be exported as a function in R but only intended to be used internally, since the
//' conversion to integer id's and all that follows can be a kind of pre-processing difficult to digest and not so intuitive)
//'
//' @param edgelist a matrix of columns : time, sender, receiver (if timing is interval, then time column is a sorted continuous vector. If timing is ordinal, then time column is still sorted but is a vector of integers)
//' @param riskset a list where for each time point a matrix of those dyads NOT included in the risk set is provided (as in [sender,receiver]).
//' @param covariates a list of covariates (actor-level, dyadic-level): we should know whether they are actor-level or dyadic-level (see how 'remstats' handle this in the remstats::remstats() function).
//'
//' @return list of objects useful for the analysis of the REH (read more here)
//'
//' @export
// [[Rcpp::export]]  
Rcpp::List reh(arma::mat edgelist, 
                Rcpp::List riskset, 
                Rcpp::List covariates){
    // (1) create a map of actors strings = integer id (Rcpp::Named("actors"))

    // (2) create the new riskset matrix [dyads*2] and [actors*actors] (the latter with column position as cell value) (Rcpp::Named("risksetMatrix") only the squared matrix)

    // (3) convert 'input' edgelist to 'output' edgelist, still with [[time],[sender],[receiver]] (Rcpp::Named("edgelist"))

    // (4) create event binary matrix from the riskset and the edgelist (Rcpp::Named("rehBinary"))

    // (5) save additional vectors/integers/information about the timinig (for example) with a proper Rcpp::Named("") 

    // (6) prepare the output list
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("value_1") = 3); // this is just an example of List with one object 'value_1' that is an integer (=3)

    return out; 
}
