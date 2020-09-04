#' remstimate  
#'
#' A function that returns maximum likelihood estimates of REM (interval timing only) by using the optim function
#'
#' @param formula formula object (see 'remstats' for now)
#' @param edgelist a a dataframe of relational events sorted by time: [time,sender,receiver,type,weight]
#' @param riskset is a list of length equal to the number of events, each object a matrix with unobserved dyads (using actors string names)
#' @param covariates list of covariates to be provided according to the input structure working with 'remstats'
#' @param stats an array of dimensions n_dyads*statistics*M (if stats is an array for relevent::rem then use aperm(a = stats, perm = c(2,3,1)) )
#' @param approach estimation approach: MLE (Maximum Likelihood Estimates) or Bayesian 
#' @param optimization optimization method to be used: trust/optim will use the corresponding functions, GD, BSIR, BPM, HMC
#' @param fast FALSE/TRUE whether to use the fast calculation of loglikelihood/gradient/hessian values inside the optimization algorithm
#' @param threshold percentage value indicating the minimum value for an actual improvment (in order to decide whether to use the fast approach for the likelihood or not).
#'
#' @return  object list (the function saves also the output of optim)
#' @export
remstimate <- function(formula = NULL,
                       edgelist,
                       riskset = replicate(NROW(edgelist),NaN,simplify=FALSE),
                       covariates = list(a = NULL),
                       stats = NULL,  
                       approach = c("MLE","Bayesian"),
                       optimization = c("trust","optim","GD","BSIR","BPM","HMC"),
                       fast = FALSE,
                       threshold = 0.5,
                       ...){

    # ... processing input _edgelist_ 
    preprocessedREH <- reh(edgelist = edgelist, riskset = riskset, covariates = covariates)

    # ... processing the input _formula_  and getting statistics (_stats_) 
    if(!is.null(formula)){
        stats <- getStats(formula = formula, edgelist = edgelist, dat = preprocessedREH, directed = TRUE, with_type = TRUE, ...)
        stats <- aperm(stats, perm = c(3,2,1))
    }
    else{
        if(is.null(stats)){
            errorMessage(0) # the array of statistics stats is not specified. Please use specify either 'formula' or 'stats'
        }
    }
    # ... processing input _method_ and _optimization_ :  adding controls over the text input by using errorMessage(cond), and threshold parameter (must be in [0,1])
    # ...

    # ... throw warningMessage(cond) according to input _threshold_ value
    # ...

    # ... creating  a local environment and an empty output list object
    remstimateList <- list()
    remstimateList$fast <- fast
    remstimateList$out <- list()

    # ... checking whether the optimization can be speeded up (according to the method proposed by DuBois et al. 2013)
    if(remstimateList$fast){
        # ... evaluate if the algorithm is actually going to be faster than by using the standard computation
        # getUniqueVector() finds the R unique vector of statistics (e.g., statistic for (i,j) and (k,j) are both U=[1,0.5,3], thus (1,0.5,3) will be present only once in the output. The function applies over dyads and time, in turn returning a matrix as output.
        remstimateList$unique_vectors_stats <- getUniqueVectors(stats = stats)
        
        # evalFast(...) output must be a percentage
        actual_improvement <-  0.0 # this is a temporary value set to 0.0
        #evalFast(M =  preprocessedREH$M,
        #                                N = preprocessedREH$N,
        #                                U = dim(stats)[3],
        #                                R = dim(remstimateList$unique_vectors_stats)[1])

        if(actual_improvement < threshold) 
        {
            # ...[message] The improvement is negligible (below the threshold), do you want to use the fast computation though?
            if(interactive()) remstimateList$fast <- askYesNo(paste("The improvement is below the threshold ", threshold,", do you want to use the fast computation though?"), default = FALSE, prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel"))))
            
            if(remstimateList$fast == NA){
                    tkmessageBox(title = "Alert message",
                    message = "Neither 'Yes' nor 'No' was chosen. Thus, the default computation will be used in the optimization", icon = "info", type = "ok")
                    remstimateList$fast <- FALSE
                }
            
            # ... if (remstimateList$fast == TRUE) 
            if(remstimateList$fast){
                # ... unique_vectors_stats is alread created and saved inside remstimateList
                # ... compute times_r : given a unique vector of statistics is the sum for each time point (m=1,...,M) of the corresponding interevent (survival) time (t[m]-t[m-1]) multiplied by the number of dyads that have the same vector of statistics at t[m].
                remstimateList$times_r <- computeTimes(unique_vectors_stats = remstimateList$unique_vectors_stats,
                                                       stats = stats, 
                                                       intereventTime = preprocessedREH$intereventTime)

                # ... compute occurrencies_r : that is a vector of counts indicating how many times each of the unique vector (indexed by r) of statistics occurred in the network (as in contributing to the hazard in the likelihood).
                remstimateList$occurrencies_r <- computeOccurrencies(edgelist = preprocessedREH$edgelist, 
                                                                     risksetMatrix = preprocessedREH$risksetMatrix, 
                                                                     M = preprocessedREH$M, 
                                                                     unique_vectors_stats = remstimateList$unique_vectors_stats, 
                                                                     stats = stats)
            }

            # ... if (remstimateList$fast == FALSE) 
            else{
                remstimateList$unique_vectors_stats <- NULL # free memory
            }
        }
    }

    # ... [1] Maximum Likelihood Estimates (MLE) or Bayesian
    if(approach == "MLE"){
        # ... [0] check that argument optimization contains only [trust,optim,DG]
        # ... [1.1] with trust::trust()
        if(optimization == "trust"){
            if(!remstimateList$fast){
                remstimateList$optimum <- trust::trust(objfun = remDerivatives, 
                                                       parinit = rep(0,dim(stats)[1]), 
                                                       rinit = 1, 
                                                       rmax = 100, 
                                                       stats = stats, #.GlobalEnv$stats, 
                                                       event_binary = preprocessedREH$rehBinary, #.GlobalEnv$binaryREH, 
                                                       interevent_time = preprocessedREH$intereventTime) #.GlobalEnv$intereventTime)  
                                                       return(remstimateList$optimum)         
            }
            else{
                remstimateList$optimum <- trust::trust(objfun = remDerivativesFast, 
                                    parinit = 0, 
                                    rinit = 1, 
                                    rmax = 100, 
                                    times_r = remstimateList$times_r,
                                    occurrencies_r = remstimateList$occurrencies_r,
                                    unique_vectors_stats = remstimateList$unique_vectors_stats)
            }
         #  remstimateList$out$coef <- remstimateList$optimum$pars #(?)
         #  remstimateList$out$hessian <- -remstimateList$optimum$hessian # hessian matrix relative
         #  remstimateList$out$llik <- -remstimateList$optimum$value # loglikelihood value at MLE values (we take the '-value' because we minimized the '-loglik')
         #  remstimateList$out$vcovMatrix <- solve(remstimateList$optimum$hessian)
         #  remstimateList$out$AIC <- 2*length(remstimateList$out$coef) - 2*remstimateList$out$llik # AIC
         #  remstimateList$out$BIC <- length(remstimateList$out$coef)*log(preprocessedREH$M) - 2*remstimateList$out$llik # BIC
        }

        # ... [1.2] with stats::optim()
        if(optimization == "optim"){        
            #ott <- optim(par = pars_init,
            #             fn = nllik,
            #             stats = stats,
            #             event_binary = event_binary,
            #             interevent_time = interevent_time,
            #             hessian = TRUE,
            #             threads = threads,
            #             ...)
            #out$coef <- ott$par # MLE estimates
            #out$hessian <- ott$hessian # hessian matrix
            #out$llik <- -ott$value # loglikelihood value at MLE values (we take the '-value' because we minimized the '-loglik')
            #out$AIC <- 2*length(out$coef) - 2*out$llik # AIC
            #out$BIC <- length(out$coef)*log(length(interevent_time)) - 2*out$llik # BIC
            #out$optim_output <- ott # optim output (if more information are needed)

        }

        # ... [1.3] with  Gradient Descent (GD)
        if(optimization == "DG"){

        }
    }

    if(approach == "Bayesian"){
        # ... [0] check for optimization argument that has to be among these: BSIR, BPM, HMC

        # ... [1.1] Bayesian Sampling Importance Resampling
        if(optimization == "BSIR"){
            # ...
        }

        # ... [1.2] Bayesian Posterior Mode
        if(optimization == "BPM"){
            # ...
        }

        # ... [1.1] Hamiltonian Monte Carlo
        if(optimization == "HMC"){
            # ...
        }
    }

    return(remstimateList$out)
}
