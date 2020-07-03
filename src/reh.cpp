#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <omp.h>
#include <typeinfo>
#include <map>
#include <iterator>
#include <string>

//' getRiskset (obtain permutations of actors' ids and event types).
//'
//' @param actorsID vector of actors' id's.
//' @param typesID vector of types' id's.
//' @param N number of actors in the dataset.
//' @param C number of event types
//'
//' @return matrix of possible dyadic events.
//'
// [[Rcpp::export]]
arma::mat getRiskset(arma::uvec actorsID, arma::uvec typesID, arma::uword N, arma::uword C){
    arma::uword i,j,c;
    arma::mat riskset(N*N,3);
    for(c = 0; c < C; c++){
        for(i = 0; i < N; i++){
            for(j = 0; j < N ; j++){
                    if(j != i){
                riskset(j+i*N+c*(N*N),0) = actorsID(i);
                riskset(j+i*N+c*(N*N),1) = actorsID(j);
                riskset(j+i*N+c*(N*N),2) = typesID(c);
                }
                else {
                    riskset(j+i*N+c*(N*N),0) = NAN;
                    riskset(j+i*N+c*(N*N),1) = NAN;
                    riskset(j+i*N+c*(N*N),2) = NAN;
                }
            }
        }
    }
    return riskset; // in R na.omit() function will be applied to this ouptut so as to filter out the rows (NaN,NaN)
}


//' getRisksetCube
//'
//' @param risksetMatrix output of getRiskset() function
//' @param N number of actors in the dataset.
//' @param C number of event types
//'
//' @return cube of possible combination [sender,receiver,type]: the cell value is the column index in the rehBinary matrix
//'
// [[Rcpp::export]]
arma::ucube getRisksetCube(arma::umat risksetMatrix, arma::uword N, arma::uword C) {
    arma::uword d;
    arma::ucube risksetCube(N,N,C);
    risksetCube.fill(N*N*C); // this is just a number to fill the cube (selfedges at each event type will remain with this value)
    for(d = 0; d < risksetMatrix.n_rows; d++){
            risksetCube(risksetMatrix(d,0),risksetMatrix(d,1),risksetMatrix(d,2)) = d;
        }
    return risksetCube; 
}


//' convertEdgelist
//'
//' @param edgelist is the input data frame with information about [time,sender,receiver,type,weight] by row.
//' @param actorsDictionary dictionary of actors names (input string name = integer id)
//' @param typesDicitonary dictionary of event types (input string name = integer id)
//'
//' @return cube of possible combination [sender,receiver,type]: the cell value is the column index in the rehBinary matrix
//'
// [[Rcpp::export]]
Rcpp::DataFrame convertEdgelist(Rcpp::DataFrame edgelist, Rcpp::List actorsDictionary, Rcpp::List typesDictionary, arma::uword M) {
    arma::uword m;
    Rcpp::NumericVector outputReceiver(M),outputType(M);
    std::vector<int> outputSender(M);
    // edgelist input to be recoded
    std::vector<std::string> stringSender = Rcpp::as<std::vector<std::string>>(edgelist["sender"]);
    std::vector<std::string> stringReceiver = Rcpp::as<std::vector<std::string>>(edgelist["receiver"]);
    std::vector<std::string> stringType = Rcpp::as<std::vector<std::string>>(edgelist["type"]);
    // strings in the dictionaries
    std::vector<std::string> actor = Rcpp::as<std::vector<std::string>>(actorsDictionary["actor"]);
    std::vector<int> actorID = actorsDictionary["actorID"];
    std::vector<std::string> type = Rcpp::as<std::vector<std::string>>(typesDictionary["type"]);
    Rcpp::NumericVector typeID = typesDictionary["typeID"];

    for(m = 0; m < M; m++){
        // find sender
        std::vector<std::string>::iterator i = std::find(actor.begin(), actor.end(), stringSender[m]);
        outputSender[m] = actorID.at(distance(actor.begin(), i));
        // find receiver
        std::vector<std::string>::iterator j = std::find(actor.begin(), actor.end(), stringReceiver[m]);
        outputReceiver[m] = actorID.at(distance(actor.begin(), j));
        // find type 
        std::vector<std::string>::iterator c = std::find(type.begin(), type.end(), stringType[m]);
        outputType[m] = typeID.at(distance(type.begin(), c));
    }
    Rcpp::DataFrame outEdgelist = Rcpp::DataFrame::create(Rcpp::Named("time") = edgelist["time"],
                                                          Rcpp::Named("sender") = outputSender,
                                                          Rcpp::Named("receiver") = outputReceiver,
                                                          Rcpp::Named("type") = outputType,
                                                          Rcpp::Named("weight") = edgelist["weight"]);
    return outEdgelist;
}
