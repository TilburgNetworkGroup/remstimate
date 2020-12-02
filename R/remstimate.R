#' remstimate  
#'
#' A function that returns maximum likelihood estimates of REM (interval timing only) by using the optim function
#'
#' @param reh a reh object (list of objects useful for the optimization)
#' @param stats a remstats object (array of statistics)
#' @param method optimization method to be used: trust, optim, GD, BSIR, HMC, BPM(?)
#' @param prior prior distribution when using Bayesian methods
#' @param proposal proposal distribution when using BSIR
#' @param n.chains number of chains to generate
#' @param n.iter number of iterations in each chain
#' @param burnin number of initial iterations to be added as burnin
#' @param parallel parallelize (especially with bayesian methods)
#' @param nthreads number of threads for the parallelization
#' @param seed seed for reproducibility
#' @param silent TRUE/FALSE if FALSE, progress of optimization status will be printed out
#'
#' @return  remstimate S3 object
#' @export
remstimate <- function(reh = NULL,
                       stats = NULL, 
                       method = c("MLE","GD","GDADAM","BSIR","BPM","HMC"), # approach = c("Frequentist","Bayesian"),
                       prior = NULL,
                       proposal = NULL,
                       n.chains = 2,
                       n.iter = 2e03,
                       burnin = 1e03,
                       parallel = FALSE,
                       nthreads = 1,
                       seed = sample(1:1e04,n.chains),
                       silent = TRUE,
                       ...){


    # ... processing input:

    # ... reh

    # ... stats 
    stats <- aperm(stats, perm = c(2,3,1)) # stats reshaped in [D*U*M] (this might change in the future)
    
    # ... method : 
    if(!(method %in% c("MLE","GD","GDADAM","BSIR","BPM","HMC"))){stop("The `method` specified is not available or it is mistyped")}

    # ... prior
    
    # ... proposal

    # ... seed


    # ... creating  an empty lists
    remstimateList <- list()
    out <- list()
    remstimateList$fast <- FALSE

    # ... checking whether the optimization can be speeded up (according to the method proposed by DuBois et al. 2013)
    # ... evaluate if the algorithm is actually going to be faster than by using the standard computation
    # getUniqueVectors() finds the R unique vector of statistics across time points and dyads
    remstimateList$unique_vectors_stats <- getUniqueVectors(stats = stats) 
    
    # number of unique vectors in the array of statistics
    R_obs <- dim(remstimateList$unique_vectors_stats)[1]
    # R_star is the number of unique vectors needed in order to reduce of 25% the number of operations 
    R_star <- round(((reh$M*(2*dim(stats)[2]+reh$D*(2*dim(stats)[2]+1)+2))/(4*(dim(stats)[2]+1)))*0.75) 

    if(R_obs <= R_star) 
    {
        # ...[message]
        remstimateList$fast <- askYesNo(paste("The estimation can be improved and become faster. Do you want to continue?"), default = FALSE, prompts = getOption("askYesNo", gettext(c("Yes", "No"))))
        
        if(is.na(remstimateList$fast)){
                tcltk::tkmessageBox(title = "Alert message",
                message = "Neither 'Yes' nor 'No' was chosen. Therefore, the function stopped.", icon = "info", type = "ok")
        }
        
        # ... IF fast == TRUE
        if(remstimateList$fast){
            # ... unique_vectors_stats is alread created and saved inside remstimateList
            # ... compute times_r : given a unique vector of statistics is the sum for each time point (m=1,...,M) of the corresponding interevent (survival) time (t[m]-t[m-1]) multiplied by the number of dyads that have the same vector of statistics at t[m].
            remstimateList$times_r <- computeTimes(unique_vectors_stats = remstimateList$unique_vectors_stats,
                                                    M = reh$M,
                                                    stats = stats, 
                                                    intereventTime = reh$intereventTime)

            # ... compute occurrencies_r : that is a vector of counts indicating how many times each of the unique vector (indexed by r) of statistics occurred in the network (as in contributing to the hazard in the likelihood).
            remstimateList$occurrencies_r <- computeOccurrencies(edgelist = reh$edgelist, 
                                                                    risksetCube = reh$risksetCube, 
                                                                    M = reh$M, 
                                                                    unique_vectors_stats = remstimateList$unique_vectors_stats, 
                                                                    stats = stats)
                                                                   
        }

        # ... IF fast == FALSE
        else{
            remstimateList$unique_vectors_stats <- matrix(0,2,2) # freeing memory too
            remstimateList$times_r <- 0
            remstimateList$occurrencies_r <- 0
        }
    }

    # ... [1] with trust::trust()
    if(method == "MLE"){
        remstimateList$optimum <- trust::trust(objfun = remDerivatives, 
                                                parinit = rep(0,dim(stats)[2]), 
                                                rinit = 1, 
                                                rmax = 100, 
                                                stats = stats,  
                                                event_binary = reh$rehBinary, 
                                                interevent_time = reh$intereventTime,
                                                times_r = remstimateList$times_r,
                                                occurrencies_r = remstimateList$occurrencies_r,
                                                unique_vectors_stats = remstimateList$unique_vectors_stats,
                                                fast = remstimateList$fast,
                                                gradient = TRUE,
                                                hessian = TRUE) 
        return(remstimateList$optimum)         
        # remstimateList$out$coef <- remstimateList$optimum$pars #(?)
        # remstimateList$out$hessian <- -remstimateList$optimum$hessian # hessian matrix relative
        # remstimateList$out$llik <- -remstimateList$optimum$value # loglikelihood value at MLE values (we take the '-value' because we minimized the '-loglik')
        # remstimateList$out$vcovMatrix <- solve(remstimateList$optimum$hessian)
        # remstimateList$out$AIC <- 2*length(remstimateList$out$coef) - 2*remstimateList$out$llik # AIC
        # remstimateList$out$BIC <- length(remstimateList$out$coef)*log(data$M) - 2*remstimateList$out$llik # BIC
    }

    # ... [2] with Gradient Descent (GD)
    if(method == "GD"){
        remstimateList$optimum <- GD(pars = rep(0,dim(stats)[2]),
                                    stats = stats,  
                                    event_binary = reh$rehBinary, 
                                    interevent_time = reh$intereventTime,
                                    times_r = remstimateList$times_r,
                                    occurrencies_r = remstimateList$occurrencies_r,
                                    unique_vectors_stats = remstimateList$unique_vectors_stats,
                                    fast = remstimateList$fast)
        return(remstimateList$optimum)
    }

    # ... [3] with Gradient Descent ADAM (GDADAM)
    if(method == "GDADAM"){
        remstimateList$optimum <- GDADAM(pars = rep(0,dim(stats)[2]),
                                    stats = stats,  
                                    event_binary = reh$rehBinary, 
                                    interevent_time = reh$intereventTime,
                                    times_r = remstimateList$times_r,
                                    occurrencies_r = remstimateList$occurrencies_r,
                                    unique_vectors_stats = remstimateList$unique_vectors_stats,
                                    fast = remstimateList$fast)
        return(remstimateList$optimum)
    }
   
    # ... [4] with Bayesian Sampling Importance Resampling (BSIR)
    #    if(method == "BSIR"){
            
    #    }

    # ... [5] with Bayesian Posterior Mode (BPM)
    #    if(method == "BPM"){
    #        
    #    }

    # ... [6] Hamiltonian Monte Carlo (HMC)
    #    if(method == "HMC"){
           
    #    }

    # reshape output here
    # define attributes here
    

    return(out)
}
