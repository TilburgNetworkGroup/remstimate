#' remstimate  
#'
#' A function that returns maximum likelihood estimates of REM (interval timing only) by using the optim function
#'
#' @param reh a reh object (list of objects useful for the optimization)
#' @param stats a remstats object (array of statistics)
#' @param method optimization method to be used: trust, optim, GD, BSIR, HMC, BPM(?)
#' @param prior prior distribution when using Bayesian methods
#' @param proposal proposal distribution when using BSIR
#' @param n_chains number of chains to generate
#' @param n_iters number of iterations in each chain
#' @param n_burnin number of initial iterations to be added as burnin
#' @param n_thin number of steps to skip in the posterior draws of the HMC
#' @param parallel TRUE/FALSE use parallelization
#' @param seed seed for reproducibility (yet to be integrated in the code)
#' @param silent TRUE/FALSE if FALSE, progress of optimization status will be printed out
#'
#' @return  remstimate S3 object
#' @export
remstimate <- function(reh = NULL,
                       stats = NULL, 
                       method = c("MLE","GD","GDADAM","BSIR","BPM","HMC"), # approach = c("Frequentist","Bayesian"),
                       prior = NULL,
                       proposal = NULL,
                       n_chains = 4,
                       n_iters = 2e03,
                       n_burnin = 1e03,
                       n_thin = 1,
                       parallel = TRUE,
                       seed = sample(1:1e04,n.chains),
                       silent = TRUE,
                       ...){


    # ... processing input:

    # ... reh

    # ... stats 
    stats <- aperm(stats, perm = c(2,3,1)) # stats reshaped in [D*U*M] (this might change in the future)
    
    # ... method : 
    if(!(method %in% c("MLE","GD","GDADAM","BSIR","HMC","BPM"))){stop("The `method` specified is not available or it is mistyped.")}

    # ... prior
    
    # ... proposal

    # ... parallel
    n_threads <- 1
    if(parallel){
        n_cpu_cores <- parallel::detectCores(all.tests = TRUE)
        if(!is.na(n_cpu_cores)) n_threads <- n_cpu_cores - 2
    }


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
    if(method == "BSIR"){
        # (0) create list of objects usefule to perform the bsir
        bsir <- list()
        bsir$log_posterior <- rep(0,n_iters*n_chains)
        bsir$log_proposal <- rep(0,n_iters*n_chains) # proposal can be rmvnorm or rmvt
        bsir$draws <- NULL

        # (1) generate from the proposal (importance distribution)

        # First find location and scale parameters (mle and variances and covariances matrix)
        mle_optimum <- trust::trust(objfun = remDerivatives, 
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
        bsir$mle <- mle_optimum$argument
        bsir$vcov <- solve(mle_optimum$hessian)
        
        if(is.null(proposal)){
            bsir$draws <- mvnfast::rmvt(n = n_iters*n_chains, 
                                        mu = bsir$mle, 
                                        sigma = bsir$vcov,
                                        df = 4,
                                        ncores = n_threads) # df = 4 by default
                                        
            bsir$log_proposal <- mvnfast::dmvt(X = bsir$draws,  
                                                mu = bsir$mle, 
                                                sigma = bsir$vcov,
                                                df = 4,
                                                log = TRUE,
                                                ncores = n_threads)
        }
        else{ 
            # write here when handling a proposal specified by the user
        }

        # (2) calculate log-posterior density give by log_prior + loglik
        ## (2.1) summing log_prior (if NULL it will be non-informative)
        if(!is.null(prior)){#bsir$log_posterior <- as.function(prior)(x=bsir$draws,mean=bsir$mle,vcov=bsir$vcov,...) #write here the handling of the prior specified by the user
        }
        ## (2.2) summing loglik
        loglik <- apply(bsir$draws,1,function(x) (-1)*remDerivatives(pars = x,
                                stats = stats, 
                                event_binary = reh$rehBinary, 
                                interevent_time = reh$intereventTime,
                                times_r = remstimateList$times_r,
                                occurrencies_r = remstimateList$occurrencies_r,
                                unique_vectors_stats = remstimateList$unique_vectors_stats,
                                fast = remstimateList$fast,
                                gradient = FALSE,
                                hessian = FALSE)$value)
        #loglik <- numeric(dim(bsir$draws)[1])   
        #for(iter_draw in 1:dim(bsir$draws)[1])  {loglik[iter_draw] <- remDerivativesStandardParallel(pars=bsir$draws[iter_draw,],
        #                                                                                            stats = stats,
        #                                                                                            event_binary = reh$rehBinary,
        #                                                                                            interevent_time = reh$intereventTime,
        #                                                                                            gradient = FALSE,
        #                                                                                            hessian = FALSE,
        #                                                                                            nthreads = nthreads)$value}     
        bsir$log_posterior <- bsir$log_posterior + loglik

        # (3) calculate Importance resampling log-weights (irlw)
        irlw <- (bsir$log_posterior - max(bsir$log_posterior) + min(bsir$log_proposal)) - bsir$log_proposal
        s_irlw <- log(sum(exp(irlw)) - exp(irlw)) # ISIR (Skare et al., 2003 Scandinavian Journal of Statistics)
        irlw <- irlw - s_irlw
        bsir$irlw <- irlw
        #this is a check from relevent::rem() the function after this message simply draws according to the weights without any particular calculation if the warningis printed out
        #if(any(is.na(iw)|is.nan(iw)|is.na(exp(iw)))||all(exp(iw)==0)){
        #warning("Importance weights seem to have some issues.  Carrying on as best we can.\n")
        #}
        size_post_draws <- round(ifelse(n_chains>1,n_iters,round(n_iters*0.33)))
        post_draws <-sample(x=1:dim(bsir$draws)[1],size=size_post_draws,replace=TRUE,prob=exp(irlw)) # we always draw either n_iters if n_chains>1 or 1/3 of n_iters when n_chains ==1
        bsir$draws <- bsir$draws[post_draws,]
        remstimateList$bsir <- bsir
        
       
        return(remstimateList$bsir)

        #fit$draws<-draws[sel,]
        #fit$lp.draws<-lp[sel]
        #colnames(fit$draws)<-parnam
        #fit$iw<-iw
        #fit$coef.mode<-fit$coef
        #fit$coef<-colMeans(fit$draws)
        #names(fit$coef)<-parnam

    }

    # ... [5] Hamiltonian Monte Carlo (HMC)
    if(method == "HMC"){
        mle_optimum <- trust::trust(objfun = remDerivatives, 
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
        pars_mle <- matrix(rep(mle_optimum$argument,n_chains),nrow=length(mle_optimum$argument),ncol=n_chains)
        pars_init <- pars_mle + matrix(runif(length(mle_optimum$argument)*n_chains,-0.1,0.1),nrow=length(mle_optimum$argument),ncol=n_chains)
        hmc <- HMC(pars_init = pars_init, 
                    n_iters = n_iters, 
                    n_chains = n_chains, 
                    n_burnin = n_burnin, 
                    meanPrior = rep(0,dim(stats)[2]),
                    sigmaPrior = diag(dim(stats)[2]),
                    stats = stats,  
                    event_binary = reh$rehBinary, 
                    interevent_time = reh$intereventTime,
                    times_r = remstimateList$times_r,
                    occurrencies_r = remstimateList$occurrencies_r,
                    unique_vectors_stats = remstimateList$unique_vectors_stats,
                    fast = remstimateList$fast,
                    n_thin = n_thin,
                    L = 100, 
                    epsilon = 0.01,
                    n_threads = n_threads)

        remstimateList$hmc <- hmc
        return(remstimateList$hmc)
    }

    # ... [6] with Bayesian Posterior Mode (BPM)
    #    if(method == "BPM"){
    #        
    #    }

    # reshape output here
    # define attributes here
    

    return(out)
}
