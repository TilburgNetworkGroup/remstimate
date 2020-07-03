#' reh 
#'
#' A function that returns maximum likelihood estimates of REM (interval timing only) by using the optim function
#'
#' @param edgelist is a dataframe of relational events sorted by time: [time,sender,receiver,type,weight]
#' @param riskset is a list of length the number of events, each object a matrix with unobserved dyads (using actors string names)
#' @param covariates list of covariates to be provided according to the input structure working with 'remstats'
#'
#' @return  object list (the function saves also the output of optim)
#' @export
reh <- function(edgelist, riskset, covariates){
    # [warnings/error messages here] Check here whether inputs arguments have been correctly provided

    # some circumstances we could get into
    # what if type is just 1? what if types are missing?
    # what if weights are not supported? what if some weights are missing?

    # Create output empty list
    out <- list()

    # Initializing useful variables
    out$M <- dim(edgelist)[1]
    out$intereventTime <- diff(c(0,edgelist$time))

    # Create a dictionary for actors and event types, that is like: 'string_name' = integer (IDentifier)
    actors <- unique(c(edgelist$sender,edgelist$receiver))
    types <- unique(edgelist$type)
    out$N <- length(actors)
    out$C <- length(types)
    out$D <- out$N*(out$N-1)*out$C
    out$actorsDictionary <- list(actor = actors, actorID = 0:(out$N-1))
    out$typesDictionary <- list(type = types, typeID = 0:(out$C-1) )

    # riskset objects (it is not the rehBinary but it just includes all the possible combination of [sender,receiver,type]) ...

    # ... arranged in a matrix [D*3]
    out$risksetMatrix <- matrix(na.omit(getRiskset(actorID = out$actorsDictionary$actorID, typeID = out$typesDictionary$typeID, N = out$N, C = out$C)),ncol=3)

    # ... arranged in a cube [N*N*C]
    out$risksetCube <- getRisksetCube(risksetMatrix = out$risksetMatrix, N = out$N, C = out$C)

    # Recoded edgelist according to the new id's for both actors and event types
    out$edgelist <- convertEdgelist(edgelist = edgelist, actorsDictionary = out$actorsDictionary, typesDictionary = out$typesDictionary, M = out$M) 

    # Create event binary matrix from the riskset and the edgelist, that is rehBinary
    # ...
                                    
    # Preprocess covariates here (we want to make 'remstats' understand our input)
    # ...
    
}
