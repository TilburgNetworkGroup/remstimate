#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <omp.h>
#include <typeinfo>
#include <map>
#include <iterator>
#include <string>
#include "messages.h"


#define LOG(x) std::cout << x << "\n"

//' getRisksetMatrix (obtain permutations of actors' ids and event types).
//'
//' @param actorID vector of actors' id's.
//' @param typeID vector of types' id's.
//' @param N number of actors in the dataset.
//' @param C number of event types
//'
//' @return matrix of possible dyadic events.
//'
// [[Rcpp::export]]
arma::mat getRisksetMatrix(arma::uvec actorID, arma::uvec typeID, arma::uword N, arma::uword C){
    arma::uword i,j,c;
    arma::mat riskset(N*N*C,3);
    arma::uvec indices_to_shed(N*C); // this is the vector where to store the indices of selfedges to remove at the end of the function
    indices_to_shed.fill(N*N*C);
    for(c = 0; c < C; c++){
        for(i = 0; i < N; i++){
            for(j = 0; j < N ; j++){
                    if(j != i){
                    riskset(j+i*N+c*(N*N),0) = actorID(i);
                    riskset(j+i*N+c*(N*N),1) = actorID(j);
                    riskset(j+i*N+c*(N*N),2) = typeID(c);
                }
                else {
                    indices_to_shed(j+c*N) = (j+i*N+c*(N*N));
                }
            }
        }
    }
    riskset.shed_rows(indices_to_shed);
    return riskset; 
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


//' convertInputREH
//'
//' @param edgelist is the input data frame with information about [time,sender,receiver,type,weight] by row.
//' @param riskset riskset list with old actors sitring names.
//' @param actorsDictionary dictionary of actors names (input string name = integer id)
//' @param typesDicitonary dictionary of event types (input string name = integer id)
//' @param M number of observed relational events
//'
//' @return cube of possible combination [sender,receiver,type]: the cell value is the column index in the rehBinary matrix
//'
// [[Rcpp::export]]
Rcpp::List convertInputREH(Rcpp::DataFrame edgelist, Rcpp::List riskset, Rcpp::List actorsDictionary, Rcpp::List typesDictionary, arma::uword M) {

    arma::uword m,d;
    Rcpp::NumericVector outputSender(M),outputReceiver(M),outputType(M);
    Rcpp::List outRiskset = Rcpp::List::create();

    // Creating output list object
    Rcpp::List out = Rcpp::List::create();

    // edgelist input to be recoded
    std::vector<std::string> stringSender = Rcpp::as<std::vector<std::string>>(edgelist["sender"]);
    std::vector<std::string> stringReceiver = Rcpp::as<std::vector<std::string>>(edgelist["receiver"]);
    std::vector<std::string> stringType = Rcpp::as<std::vector<std::string>>(edgelist["type"]);
    // strings in the dictionaries
    std::vector<std::string> actor = Rcpp::as<std::vector<std::string>>(actorsDictionary["actor"]);
    std::vector<int> actorID = actorsDictionary["actorID"];
    std::vector<std::string> type = Rcpp::as<std::vector<std::string>>(typesDictionary["type"]);
    std::vector<int> typeID = typesDictionary["typeID"];

    for(m = 0; m < M; m++){
        // (1) converting m-th event in the edgelist input
        // find sender
        std::vector<std::string>::iterator i = std::find(actor.begin(), actor.end(), stringSender[m]);
        outputSender[m] = actorID.at(std::distance(actor.begin(), i));
        // find receiver
        std::vector<std::string>::iterator j = std::find(actor.begin(), actor.end(), stringReceiver[m]);
        outputReceiver[m] = actorID.at(std::distance(actor.begin(), j));
        // find type 
        std::vector<std::string>::iterator c = std::find(type.begin(), type.end(), stringType[m]);
        outputType[m] = typeID.at(std::distance(type.begin(), c));

        // (2) converting m-th object in the riskset input list
        if(!R_IsNaN(riskset[m])){ //
            // we expect a Rcpp::DataFrame (sender,receiver,type)
            Rcpp::StringMatrix omitDyads = riskset[m];

            // number of dyads to omit from the riskset at the m-th time point
            arma::uword D_loc = omitDyads.nrow();

            // converting input senders, receivers and types to std::string vectors
            Rcpp::StringVector omitDyadsSender = omitDyads.column(0);
            Rcpp::StringVector omitDyadsReceiver = omitDyads.column(1);
            Rcpp::StringVector omitDyadsType = omitDyads.column(2);
            std::vector<std::string> stringOmitSender = Rcpp::as<std::vector<std::string>>(omitDyadsSender); // sender
            std::vector<std::string> stringOmitReceiver = Rcpp::as<std::vector<std::string>>(omitDyadsReceiver); // receiver
            std::vector<std::string> stringOmitType = Rcpp::as<std::vector<std::string>>(omitDyadsType); // type

            // allocating space for the converted senders, receivers and types
            Rcpp::NumericVector omitSender(D_loc),omitReceiver(D_loc),omitType(D_loc);

            // find sender
            for(d = 0; d < D_loc; d++){
                std::vector<std::string>::iterator i = std::find(actor.begin(), actor.end(), stringOmitSender[d]);
                omitSender[d] = actorID.at(std::distance(actor.begin(), i));
                // find receiver
                std::vector<std::string>::iterator j = std::find(actor.begin(), actor.end(), stringOmitReceiver[d]);
                omitReceiver[d] = actorID.at(std::distance(actor.begin(), j));
                // find type 
                std::vector<std::string>::iterator c = std::find(type.begin(), type.end(), stringOmitType[d]);
                omitType[d] = typeID.at(std::distance(type.begin(), c));
            }
            Rcpp::DataFrame converted_loc = Rcpp::DataFrame::create(Rcpp::Named("sender") = omitSender,
                                                                    Rcpp::Named("receiver") = omitReceiver,
                                                                    Rcpp::Named("type") = omitType);
            outRiskset.push_back(converted_loc);
        }
        else{
            outRiskset.push_back(R_NaN);
        }

    }
    Rcpp::DataFrame outEdgelist = Rcpp::DataFrame::create(Rcpp::Named("time") = edgelist["time"],
                                                          Rcpp::Named("sender") = outputSender,
                                                          Rcpp::Named("receiver") = outputReceiver,
                                                          Rcpp::Named("type") = outputType,
                                                          Rcpp::Named("weight") = edgelist["weight"]);
    out["edgelist"] = outEdgelist; 
    out["riskset"] = outRiskset;
                                                       
    return out;
}

//' getBinaryREH (a function that returns a utility matrix used in optimization algorithms)
//'
//' @param edgelist edgelist converted according to actorID and typeID
//' @param riskset riskset list converted according to actorID and typeID
//' @param risksetCube arma::cube object [N*N*C] where the cell value returns the column index to use in the outBinaryREH
//' @param M number of observed relational events
//' @param D number of possible dyads (accounting for event types as well)
//'
//' @return utility matrix per row 0 if the event could happen but didn't, 1 if the event happend, -1 if the event couldn't occur
//' 
// [[Rcpp::export]]
arma::mat getBinaryREH(Rcpp::DataFrame edgelist, Rcpp::List riskset, arma::ucube risksetCube, arma::uword M, arma::uword D) {
    arma::uword m,d;
    arma::mat outBinaryREH(M,D,arma::fill::zeros); // by setting the initial values to zero we already handle those
                                                    // relational events that could have occurred but didn't
    Rcpp::NumericVector sender = edgelist["sender"];
    Rcpp::NumericVector receiver = edgelist["receiver"];
    Rcpp::NumericVector type = edgelist["type"];
    for(m = 0; m < M; m++){
        // relational event that occurred
        arma::uword event_m = risksetCube(sender[m],receiver[m],type[m]);
        outBinaryREH(m,event_m) = 1;
        // relational events that couldn't occur
        if(!R_IsNaN(riskset[m])){
            // we expect a Rcpp::DataFrame (sender,receiver,type)
            arma::mat omitDyads = riskset[m];

            // number of dyads to omit from the riskset at the m-th time point
            arma::uword D_loc = omitDyads.n_rows;
            // getting sender, receiver and type combinations to omit from the riskset at the m-th time point
            arma::vec omitSender = omitDyads.col(0); // sender
            arma::vec omitReceiver = omitDyads.col(1); // receiver
            arma::vec omitType = omitDyads.col(2); // type
            for(d = 0; d < D_loc; d++){
                arma::uword event_d = risksetCube(omitSender(d),omitReceiver(d),omitType(d));
                outBinaryREH(m,event_d) = -1;
            }
        }
    }

    return outBinaryREH;
}


//' reh (a function for preprocessing data)
//'
//' @param edgelist is a dataframe of relational events sorted by time: [time,sender,receiver,type,weight]
//' @param riskset is a list of length the number of events, each object a matrix with unobserved dyads (using actors string names)
//' @param covariates list of covariates to be provided according to the input structure working with 'remstats'
//'
//' @return list of objects
//' 
//' @export 
// [[Rcpp::export]]
Rcpp::List reh(Rcpp::DataFrame edgelist, Rcpp::List riskset, Rcpp::List covariates) {
    // [warnings/error messages here] Check here whether inputs arguments have been correctly provided

    // some circumstances we could get into
    // what if type is just 1? what if types are missing?
    // what if weights are not supported? what if some weights are missing?

    // Create output empty list
    Rcpp::List out = Rcpp::List::create();

    // START of the processing

    // intereventTime variable (we should handle in such a way that if Date class variable is provided we can process that data type as well)
    // if(!is.null(start_time)){
    //  if(is.Date(start_time) | is.numeric(start_time)){
        // process the time variable here (getSeconds(), getMinutes(), getHours(), getDays(), getYears(), user should be able to specify the time scale)
        // calculate the intereventTime variable here
        // the starting point (t_0) is important to define: it can either be a day/hour/minute before by convention or a prespeficied by the user
    //  }
    //  else{stop(errorMessage(cond = 1))} // errorMessage() is a function coming from the header file messages.h (loaded on line 9 in this file)
    // }
    // else{ ...  process here when there is no start_time argument define, as it is by default ... }
    Rcpp::NumericVector time_loc = edgelist["time"];
    double t0 = 0.0; //it can be an argument of the function, a check for t0 being the actual minimum value must be provided in the future (just for now it is internally set)
    time_loc.push_front(t0);
    out["intereventTime"] = Rcpp::diff(time_loc); 

    // StringVector of senders
    Rcpp::StringVector senders = edgelist["sender"];

    // StringVector of receivers
    Rcpp::StringVector receivers = edgelist["receiver"];

    //StringVector of senders and receivers (this procedure of preallocation the whole StrigVector works faster than a for loop with of push_front's)
    Rcpp::StringVector senders_and_receivers(senders.length()+receivers.length());
    senders_and_receivers[Rcpp::Range(0,(senders.length()-1))] = senders;
    senders_and_receivers[Rcpp::Range(senders.length(),(senders_and_receivers.length()-1))] = receivers;

    // Finding unique strings in senders_and_receivers 
    Rcpp::StringVector actor = Rcpp::unique(senders_and_receivers);

    // Finding unique strings in event types  
    Rcpp::StringVector vector_of_types = edgelist["type"];
    Rcpp::StringVector type = Rcpp::unique(vector_of_types);

    // Storing some useful dimensions
    out["N"] = actor.length();
    out["C"] = type.length();
    out["D"] = actor.length()*(actor.length()-1)*type.length();
    out["M"] = edgelist.nrows();

    // Creating a dictionary for actors and event types, that is like: 'string_name' = integer (IDentifier)
    Rcpp::List actorsDictionary = Rcpp::List::create(Rcpp::Named("actor") = actor, Rcpp::Named("actorID") = Rcpp::Range(0,actor.length()-1)); 
    out["actorsDictionary"] = actorsDictionary;
    Rcpp::List typesDictionary = Rcpp::List::create(Rcpp::Named("type") = type, Rcpp::Named("typeID") = Rcpp::Range(0,type.length()-1)); 

    // Creating riskset objects (it is not the rehBinary but it just includes all the possible combination of [sender,receiver,type]) ...

    // ... arranged in a matrix [D*3]
    out["risksetMatrix"] = getRisksetMatrix(actorsDictionary["actorID"],typesDictionary["typeID"],out["N"],out["C"]);


    // ... arranged in a cube [N*N*C]
    out["risksetCube"] = getRisksetCube(out["risksetMatrix"],out["N"],out["C"]);

    // Converting input edgelist and riskset list according to the new id's for both actors and event types
    Rcpp::List convertedInput = convertInputREH(edgelist,riskset,actorsDictionary,typesDictionary,out["M"]);
    out["edgelist"] = convertedInput["edgelist"];
    out["riskset"] = convertedInput["riskset"];

    // Create event binary matrix from the riskset and the edgelist, that is rehBinary
    out["rehBinary"] = getBinaryREH(out["edgelist"],out["riskset"],out["risksetCube"],out["M"],out["D"]);
                                    
    // Preprocess covariates here (we want to make 'remstats' understand our input)
    // ...

    // END of the processing and returning output

    return out;
}
