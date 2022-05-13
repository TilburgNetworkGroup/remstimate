#' remstimate  
#'
#' A function for the optimization of tie-oriented (Relational Event Model) or actor-oriented (DyNAM) likelihoods
#'
#' @param reh an \code{reh} object (output object of the function \code{remify::reh})
#' @param stats a \code{remstats} object: when \code{model="tie"}, \code{stats} is an array of statictis with dimensions \code{[M x D x P]}: where \code{M} is the number of events, \code{D} is the number of possible dyads (full riskset), \code{P} is the number of statistics; if \code{model="actor"}, \code{stats} is a list of two arrays named \code{rate} and \code{choice} with dimensions \code{[M x N x P]}, where \code{N} are the actors (senders in the array \code{rate}, receivers in the array \code{choice})
#' @param method optimization method to be used: \code{MLE}, \code{GD}, \code{GDADAM}, \code{BSIR}, \code{HMC}
#' @param model either \code{"actor"} or \code{"tie"} oriented. \code{"tie"} is the default value (literature ref. here to DyNAM and to REM?)
#' @param ncores number of threads for the parallelization (\code{default = 1}, no parallelization)
#' @param prior prior distribution when using Bayesian methods
#' @param nsim  when \code{method = "HMC"} it is the number of simulations (iterations) in each chain, when \code{method = "BSIR"} is the number of samples from the proposal distribution
#' @param nchains number of chains to generate
#' @param burnin number of initial iterations to be added as burnin
#' @param thin number of steps to skip in the posterior draws of the HMC
#' @param init matrix of initial values for \code{"HMC"}
#' @param epochs 200 by defaut. It is the number of iteration used in the methods \code{"GD"} and \code{"GDADAM"}
#' @param learning_rate 0.001 by default. It is the learning rate used in the methods \code{"GD"} and \code{"GDADAM"}
#' @param seed seed for reproducibility (yet to be integrated in the code)
#' @param silent \code{TRUE/FALSE} if \code{FALSE}, progress of optimization status will be printed out
#' @param seed seed for reproducibility (yet to be integrated in the code)
#' @param ... additional parameters. They can be parameters of other functions defined as input in some argument
#'
#' @return  remstimate S3 object
#' @export
remstimate <- function(reh = NULL,
                       stats = NULL, 
                       method = c("MLE","GD","GDADAM","BSIR","HMC"), 
                       model = c("actor","tie"),
                       ncores = 1,
                       prior = NULL,
                       nsim = 1000,
                       nchains = 2,
                       burnin = 500,
                       thin = 10,  
                       init = NULL,                 
                       epochs = 200,
                       learning_rate = 0.001,
                       silent = TRUE,
                       seed = sample(1:1e04,nchains),
                       ...){

    # ... processing input:

    # ... reh
    if(is.null(reh)) stop("missing 'reh' argument.")
    else{
        if(class(reh)!="reh"){
            if(is.data.frame(reh) | is.matrix(reh)){
                reh <- remify::reh(edgelist = reh, ...)
            }
            else{
                stop("Input 'reh' must be either a 'reh' object (from 'remify' package) or the event sequence (matrix or data.frame).")
            }
        }
    }

    # ... model :
    if(is.null(model)) {
        model <- "tie"
        warning("argument 'model' set to 'tie' by default")
    }
    if(!is.null(model) & !(model %in% c("tie","actor"))) stop("argument 'model' must be set to either 'tie' or 'actor'")
    if(!is.null(model) & all(model==c("tie","actor"))) model <- "tie"
    
    # ... method : 
    if(!(method %in% c("MLE","GD","GDADAM","BSIR","HMC"))){stop("The `method` specified is not available or it is mistyped.")}

    # ... type of likelihood
    ordinal <- attr(reh,"ordinal")

    # ... stats 
    # check if 'stats' is an object of class 'remstats'
    model_formula <- variable_names <- NULL
    if(model == "tie")
    {
        if((class(stats)[2]=="remstats") & (class(stats)[1]=="tomstats")){
            variable_names <- as.vector(sapply(dimnames(stats$statistics)[[3]],function(x) sub(pattern = ".x.", replacement = ":", x = x)))
            model_formula <- stats::as.formula(paste("~ ",paste(variable_names,collapse=" + ")))
            stats <- aperm(stats$statistics, perm = c(2,3,1)) # stats reshaped in [D*U*M] (this might change in the future?)
        }
        else if((class(stats)[2]=="remstats") & (class(stats)[1]=="aomstats")){
            stop("'remstats' object supplied cannot work for tie-oriented modeling")
        }
        else if(is.array(stats)){
            M_stats <- dim(stats)[1]
            D_stats <- dim(stats)[2]
            if((M_stats != reh$M) | (D_stats != reh$D)){
                stop("numbers of rows and columns of the array 'stats' must be equal to number of 'time points' and 'dyads' in the network (see dim(reh))")
            }
            else{
                stats <- aperm(stats, perm =c(2,3,1))
            }
            model_formula <- ""
        }
        else{
            stop("tie-oriented modeling: 'stats' must be either a 'tomstats' 'remstats' object or an array of statistics")
        }
    }
    if(model == "actor") 
    {
        if((class(stats)[2]=="remstats") & (class(stats)[1]=="aomstats")){
            variables_rate <- as.vector(sapply(dimnames(stats$statistics$sender_stats)[[3]],function(x) sub(pattern = ".x.", replacement = ":", x = x)))
            variables_choice <- as.vector(sapply(dimnames(stats$statistics$receiver_stats)[[3]],function(x) sub(pattern = ".x.", replacement = ":", x = x)))
            variable_names <- c(variables_rate,variables_choice)
            model_formula <- list(rate_model_formula = stats::as.formula(paste("~ ",paste(variables_rate,collapse=" + "))), choice_model_formula = stats::as.formula(paste("~ ",paste(variables_choice,collapse=" + "))))
            stats <- list(rate = aperm(stats$statistics$sender_stats, perm = c(2,3,1)), choice = aperm(stats$statistics$receiver_stats, perm = c(2,3,1))) # stats reshaped in [D*U*M] (this might change in the future?)
        }
        else if((class(stats)[2]=="remstats") & (class(stats)[1]=="tomstats")){
            stop("'remstats' object supplied cannot work for tie-oriented modeling")
        }
        else if(is.list(stats) & all(c("rate","choice") %in% names(stats))){
             M_rate <- dim(stats$rate)[1]
             N_rate <- dim(stats$rate)[2]
             M_choice <- dim(stats$choice)[1]
             N_choice <- dim(stats$choice)[2]
             if((M_rate != reh$M) & (M_choice != reh$M) & (N_rate != reh$N) & (N_choice != reh$N)){
                 stop("numbers of rows and columns in both arrays ('rate' and 'choice') inside 'stats' must be equal to number of 'time points' and 'actors' in the network (see dim(reh))")
             }
             else{
                 stats <- list(rate = aperm(stats$rate, perm = c(2,3,1)), choice = aperm(stats$choice, perm = c(2,3,1)))
             }
            model_formula <- stats::as.formula(".~.")
        }
        else{
            stop("actor-oriented modeling: 'stats' must be either a 'aomstats' 'remstats' object or a list of two arrays named 'rate' and 'choice'")
        }
    }

    # ... prior
    log <- TRUE
    if(!is.null(prior)){
        additional_input_args <- names(list(...)) #names(as.list(match.call()))[-1]
        args_prior <- match.arg(arg=additional_input_args,choices = methods::formalArgs(prior),several.ok = TRUE)
    }


    # ... ncores
    if(is.null(ncores)) ncores <- 1
    else{
      if(ncores > (floor(parallel::detectCores()) - 2)) stop(cat("'ncores' is recommended to be set at most to ", (floor(parallel::detectCores()) - 2),"."))
    }

    # ... seed
    if(is.null(seed) & !is.null(nchains)) seed <- sample(1:1e04,nchains)


    # ... creating  an empty lists
    remstimateList <- list()
        
    # ... [1] Maximum Likelihood, Gradient Descent (GD), Gradient Descent ADAM (GDADAM)
    if(method %in% c("MLE","GD","GDADAM")){
        if(model == "tie"){ # Relational Event Model (REM)
            if(method == "MLE"){
                optimum_obj <- trust::trust(objfun = remDerivatives, 
                                            parinit = rep(0,dim(stats)[2]), 
                                            rinit = 1, 
                                            rmax = 100, 
                                            stats = stats,  
                                            edgelist = data.matrix(reh$edgelist),
                                            omit_dyad = reh$omit_dyad,
                                            interevent_time = reh$intereventTime,
                                            model = model,
                                            ordinal = ordinal,
                                            ncores = ncores)
            }  
            if(method == "GD"){
                if(is.null(init)){
                    optimum_obj <- GD(pars = stats::rnorm(dim(stats)[2]), # rep(0,dim(stats)[2])
                                            stats = stats,  
                                            edgelist = reh$edgelist,
                                            omit_dyad = reh$omit_dyad,
                                            interevent_time = reh$intereventTime,
                                            model = model,
                                            epochs = epochs,
                                            learning_rate = learning_rate)
                }else{
                    optimum_obj <- GD(pars = init,
                                            stats = stats,  
                                            edgelist = reh$edgelist,
                                            omit_dyad = reh$omit_dyad,
                                            interevent_time = reh$intereventTime,
                                            model = model,
                                            epochs = epochs,
                                            learning_rate = learning_rate)
                }
                optimum_obj$hessian <- remstimate::remDerivativesStandard(par = optimum_obj$argument,
                                                                    stats = stats,  
                                                                    edgelist = data.matrix(reh$edgelist),
                                                                    omit_dyad = reh$omit_dyad,
                                                                    interevent_time = reh$intereventTime,
                                                                    ordinal = ordinal,
                                                                    ncores = ncores)
                optimum_obj$hessian <- optimum_obj$hessian$hessian                                                  
            }   
            if(method == "GDADAM"){
                if(is.null(init)){
                    optimum_obj <- GDADAM(pars = stats::rnorm(dim(stats)[2]), # rep(0,dim(stats)[2])
                                            stats = stats,
                                            edgelist = reh$edgelist,
                                            omit_dyad = reh$omit_dyad,
                                            interevent_time = reh$intereventTime,
                                            model = model,
                                            epochs = epochs,
                                            learning_rate = learning_rate)
                }else{
                    optimum_obj <- GDADAM(pars = init,
                                            stats = stats,
                                            edgelist = reh$edgelist,
                                            omit_dyad = reh$omit_dyad,
                                            interevent_time = reh$intereventTime,
                                            model = model,
                                            epochs = epochs,
                                            learning_rate = learning_rate)
                }
                optimum_obj$hessian <- remstimate::remDerivativesStandard(par = optimum_obj$argument,
                                                                    stats = stats,  
                                                                    edgelist = data.matrix(reh$edgelist),
                                                                    omit_dyad = reh$omit_dyad,
                                                                    interevent_time = reh$intereventTime,
                                                                    ordinal = ordinal,
                                                                    ncores = ncores)
                optimum_obj$hessian <- optimum_obj$hessian$hessian
            }                                  
            remstimateList$coefficients <- optimum_obj$argument 
            remstimateList$loglik <- -optimum_obj$value # loglikelihood value at MLE values (we take the '-value' because we minimized the '-loglik')        
            remstimateList$gradient <-optimum_obj$gradient
            remstimateList$hessian <- optimum_obj$hessian # hessian matrix relative
            remstimateList$vcov <- qr.solve(optimum_obj$hessian) # matrix of variances and covariances
            remstimateList$se <- diag(remstimateList$vcov)**0.5 # standard errors
            names(remstimateList$coefficients) <- names(remstimateList$se) <- rownames(remstimateList$vcov) <- colnames(remstimateList$vcov) <- dimnames(stats)[[2]]
            remstimateList$residual.deviance <- -2*remstimateList$loglik
            remstimateList$null.deviance <- 2*(trust::trust(objfun = remDerivatives, 
                                        parinit = c(0), 
                                        rinit = 1, 
                                        rmax = 100, 
                                        stats = array(1,dim=c(reh$D,1,reh$M)),  
                                        edgelist = reh$edgelist,
                                        omit_dyad = reh$omit_dyad,
                                        interevent_time = reh$intereventTime,
                                        model = model,
                                        ordinal = ordinal,
                                        ncores = ncores)$value)
            remstimateList$model.deviance <- remstimateList$null.deviance -  remstimateList$residual.deviance                     
            remstimateList$df.null <- reh$M
            remstimateList$df.model <- dim(stats)[2]

            remstimateList$AIC <- 2*length(remstimateList$coef) - 2*remstimateList$loglik # AIC
            remstimateList$AICC <- remstimateList$AIC + 2*length(remstimateList$coef)*(length(remstimateList$coef)+1)/(reh$M-length(remstimateList$coef)-1)
            remstimateList$BIC <- length(remstimateList$coef)*log(reh$M) - 2*remstimateList$loglik # BIC
            remstimateList$converged <- optimum_obj$converged
            remstimateList$iterations <- optimum_obj$iteration                                    
        }
        if(model == "actor" ){ # Actor Oriented Model 
            if(method == "MLE"){
                optimum_sender_rate <- trust::trust(objfun = remDerivatives, 
                                    parinit = rep(0,dim(stats$rate)[2]), 
                                    rinit = 1, 
                                    rmax = 100, 
                                    stats = stats$rate, 
                                    edgelist = reh$edgelist,
                                    omit_dyad = reh$omit_dyad,
                                    interevent_time = reh$intereventTime,
                                    model = model,
                                    C = reh$C,
                                    D = reh$D,
                                    ordinal = ordinal,
                                    senderRate = TRUE)
                if(!silent) {print(optimum_sender_rate)}   # this print will be removed in the public version of the package         
                optimum_receiver_choice <- trust::trust(objfun = remDerivatives, 
                                    parinit = rep(0,dim(stats$choice)[2]), 
                                    rinit = 1, 
                                    rmax = 100, 
                                    stats = stats$choice, 
                                    edgelist = reh$edgelist,
                                    omit_dyad = reh$omit_dyad, 
                                    interevent_time = reh$intereventTime,
                                    model = model,
                                    senderRate = FALSE,
                                    N = reh$N,
                                    C = reh$C,
                                    D = reh$D)
            }
            if(method == "GD"){
                stop("Gradient Descent (GD) is not yet available for actor-oriented modeling")
            }
            if(method == "GDADAM"){
                stop("Gradient Descent with ADAM optimization (GDADAM) is not yet available for actor-oriented modeling")
            }

            if(!silent) {print(optimum_receiver_choice)}   # this print will be removed in the public version of the package                   
            remstimateList$coefficients <- c(optimum_sender_rate$argument, optimum_receiver_choice$argument)
            remstimateList$loglik <- -(optimum_sender_rate$value+optimum_receiver_choice$value) # log(L_sender) + log(L_choice)       
            remstimateList$gradient <- rbind(optimum_sender_rate$gradient, optimum_receiver_choice$gradient)
            remstimateList$hessian <- matrix(0,nrow=(dim(stats$rate)[2]+dim(stats$choice)[2]),ncol=(dim(stats$rate)[2]+dim(stats$choice)[2]))
            remstimateList$hessian[1:dim(stats$rate)[2],1:dim(stats$rate)[2]] <- optimum_sender_rate$hessian # hessian matrix for sender rate model
            remstimateList$hessian[-c(1:dim(stats$rate)[2]),-c(1:dim(stats$rate)[2])] <- optimum_receiver_choice$hessian # hessian matrix for choice model
            remstimateList$vcov <- qr.solve(remstimateList$hessian) # matrix of variances and covariances
            remstimateList$se <- diag(remstimateList$vcov)**0.5 # standard errors
            names(remstimateList$coefficients) <- names(remstimateList$se) <- rownames(remstimateList$vcov) <- colnames(remstimateList$vcov) <- colnames(remstimateList$hessian) <- rownames(remstimateList$hessian) <- c(dimnames(stats$rate)[[2]],dimnames(stats$choice)[[2]])
            remstimateList$residual.deviance <- -2*remstimateList$loglik
            remstimateList$null.deviance.sender <- 2*(trust::trust(objfun = remDerivatives, 
                                                                    parinit = c(0), 
                                                                    rinit = 1, 
                                                                    rmax = 100, 
                                                                    stats = array(1,dim=c(dim(stats$rate)[1],1,reh$M)),  
                                                                    edgelist = reh$edgelist,
                                                                    omit_dyad = reh$omit_dyad,
                                                                    interevent_time = reh$intereventTime,
                                                                    model = model,
                                                                    C = reh$C,
                                                                    D = reh$D,
                                                                    ordinal = ordinal,
                                                                    senderRate = TRUE)$value)
                                                     
            remstimateList$null.deviance.choice <- ifelse(length(reh$omit_dyad)>0,2*(trust::trust(objfun = remDerivatives, 
                                                                parinit = c(0), 
                                                                rinit = 1, 
                                                                rmax = 100, 
                                                                stats = array(1,dim=c(dim(stats$choice)[1],1,reh$M)), 
                                                                edgelist = reh$edgelist,
                                                                omit_dyad = reh$omit_dyad, 
                                                                interevent_time = reh$intereventTime,
                                                                model = model,
                                                                senderRate = FALSE,
                                                                N = reh$N,
                                                                C = reh$C,
                                                                D = reh$D)$value),-2*log(1/(reh$N-1)))
            # null.deviance.choice is -2*log(1/(reh$N-1)) when the riskset is always (all the time points) the same
            remstimateList$null.deviance <- remstimateList$null.deviance.choice + remstimateList$null.deviance.choice                       
            remstimateList$model.deviance <- remstimateList$null.deviance -  remstimateList$residual.deviance                     
            remstimateList$df.null <- reh$M
            remstimateList$df.model <- dim(stats$rate)[2] + dim(stats$choice)[2]

            remstimateList$AIC <- 2*length(remstimateList$coefficients) - 2*remstimateList$loglik # AIC
            remstimateList$AICC <- remstimateList$AIC + 2*length(remstimateList$coefficients)*(length(remstimateList$coefficients)+1)/(reh$M-length(remstimateList$coefficients)-1)
            remstimateList$BIC <- length(remstimateList$coefficients)*log(reh$M) - 2*remstimateList$loglik # BIC
            remstimateList$converged <- c("rate" = optimum_sender_rate$converged, "choice" = optimum_receiver_choice$converged)
            remstimateList$iterations <- c("rate" = optimum_sender_rate$iteration, "choice" = optimum_receiver_choice$iteration)
        }            
        # ...
        # in future : if(WAIC) and if(ELPD){if(PSIS)} for predictive power of models
    }
   
    # ... [4] with Bayesian Sampling Importance Resampling (BSIR)
    if(method == "BSIR"){
        # (0) create list of objects usefule to perform the bsir
        bsir <- list()
        if(is.null(nsim)){
            nsim <- 5e03
            warning("'nsim' is set to 5e03 by default")
        }
        bsir$log_posterior <- rep(0,nsim)  # 
        bsir$draws <- list()

        # (1) generate from the proposal (importance distribution)

        # First find location and scale parameters (mles and variances and covariances matrix)
        mle_optimum <- trust::trust(objfun = remDerivatives, 
                                parinit = rep(0,dim(stats)[2]), 
                                rinit = 1, 
                                rmax = 100, 
                                stats = stats,  
                                event_binary = reh$rehBinary, 
                                interevent_time = reh$intereventTime,
                                model = model,
                                ordinal = ordinal,
                                ncores = ncores)
        bsir$mle <- mle_optimum$argument
        bsir$vcov <- qr.solve(mle_optimum$hessian)
        
        # proposal distribution (default proposal is a multivariate Student t): drawing nsim*3 samples
        bsir$draws[[1]] <- mvnfast::rmvt(n = nsim*3, 
                                        mu = bsir$mle, 
                                        sigma = bsir$vcov,
                                        df = 4,
                                        ncores = ncores) # df = 4 by default
                                    
        bsir$log_proposal <- mvnfast::dmvt(X = bsir$draws[[1]],  
                                            mu = bsir$mle, 
                                            sigma = bsir$vcov,
                                            df = 4,
                                            log = TRUE,
                                            ncores = ncores)

        # (2) calculate log-posterior density give by log_prior + loglik
        ## (2.1) summing log_prior (if NULL it will be non-informative) : prior by default can be a multivariate Student t and the user must specify its parameters
        if(!is.null(prior)){
            bsir$log_prior <- prior(bsir$draws[[1]],...) 
            bsir$log_posterior <- bsir$log_prior
        }
        ## (2.2) summing loglik
        loglik <- apply(bsir$draws[[1]],1,function(x) (-1)*remDerivatives(pars = x,
                                stats = stats, 
                                event_binary = reh$rehBinary, 
                                interevent_time = reh$intereventTime,
                                model = model,
                                ordinal = ordinal,
                                ncores = ncores,
                                gradient = FALSE,
                                hessian = FALSE)$value) # this is the most time consuming step to be optimized, avoiding the apply
   
        bsir$log_posterior <- bsir$log_posterior + loglik

        # (3) calculate Importance resampling log-weights (irlw)
        irlw <- (bsir$log_posterior - max(bsir$log_posterior) + min(bsir$log_proposal)) - bsir$log_proposal
        s_irlw <- log(sum(exp(irlw)) - exp(irlw)) # ISIR (Skare et al., 2003 Scandinavian Journal of Statistics)
        irlw <- irlw - s_irlw
        bsir$irw <- exp(irlw) # importance resampling weights (irw)
        #this is a check from relevent::rem() the function after this message simply draws according to the weights without any particular calculation if the warningis printed out
        #if(any(is.na(iw)|is.nan(iw)|is.na(exp(iw)))||all(exp(iw)==0)){
        #warning("Importance weights seem to have some issues.  Carrying on as best we can.\n")
        #}
        post_draws <-sample(x=1:dim(bsir$draws[[1]])[1],size=nsim,replace=TRUE,prob=exp(irlw)) # we always draw nsim (having generated 3*nsim samples from the proposal)
        bsir$draws[[1]] <- bsir$draws[[1]][post_draws,]
        bsir$log_posterior <- bsir$log_posterior[post_draws]
        bsir$post.mode <- bsir$draws[[1]][which.max(bsir$log_posterior),]
        bsir$post.mean <- colMeans(bsir$draws[[1]])
        bsir$vcov <- stats::cov(bsir$draws[[1]])
        bsir$sd <- diag(bsir$vcov)**0.5
        names(bsir$post.mode) <- names(bsir$post.mean) <- rownames(bsir$vcov) <- colnames(bsir$vcov) <- names(bsir$sd) <- dimnames(stats)[[2]]
        
        remstimateList <- bsir
    }

    # ... [5] Hamiltonian Monte Carlo (HMC)
    if(method == "HMC"){
        hmc <- list()
        # remove mles as initial values 

        if(is.null(init)){
            init <- matrix(stats::runif(dim(stats)[2]*nchains,-1,1),nrow=dim(stats)[2],ncol=nchains) # runif was in (-0.1,0.1)
        }
        else{
          #  if()
        }
        mle_optimum <- trust::trust(objfun = remDerivatives, 
                                parinit = rep(0,dim(stats)[2]), 
                                rinit = 1, 
                                rmax = 100, 
                                stats = stats,  
                                event_binary = reh$rehBinary, 
                                interevent_time = reh$intereventTime,
                                model = model,
                                ordinal = ordinal,
                                ncores = ncores)
        #pars_mle <- matrix(rep(mle_optimum$argument,nchains),nrow=length(mle_optimum$argument),ncol=nchains)
        #init <- pars_mle + matrix(runif(length(mle_optimum$argument)*nchains,-0.1,0.1),nrow=length(mle_optimum$argument),ncol=nchains)

        hmc_out <- HMC(pars_init = init, 
                    nsim = (burnin+nsim), 
                    nchains = nchains, 
                    burnin = burnin, 
                    meanPrior = rep(0,dim(stats)[2]),
                    sigmaPrior = diag(dim(stats)[2]),
                    stats = stats,  
                    event_binary = reh$rehBinary, 
                    interevent_time = reh$intereventTime,
                    model = model,
                    thin = thin,
                    L = 100, 
                    epsilon = 0.01,
                    ncores = ncores)
                    # needs to have a model argument (tie or actor string)
        return(hmc_out)
        #remstimateList$hmc <- hmc # output such that one can run gelman plots and statistics // define methods like remstimate.traceplot() remstimate.
                                    # A. Gelman, Carlin, et al. (2013, 267)
                                    # Stan Development Team (2016 Ch 28.) for how Stan calculates Hat, autocorrelations, and ESS.
                                    # Gelman and Rubin (1992) introduce the R-hat statistic

        hmc$post.mode <- hmc_out$draws[which.max(hmc_out$log_posterior),]
        hmc$post.mean <- colMeans(hmc_out$draws)
        hmc$vcov <- stats::cov(hmc_out$draws)
        hmc$sd <- diag(hmc$vcov)**0.5
        names(hmc$post.mode) <- names(hmc$post.mean) <- rownames(hmc$vcov) <- colnames(hmc$vcov) <- names(hmc$sd) <- dimnames(stats)[[2]]
        remstimateList <- c(hmc_out,hmc)
    }

    str_out <- NULL
    # defining structure and attributes based on the 'method'
    if(method %in% c("MLE","GD","GDADAM")){ # frequentist approach
        # defining structure of class 'remstimate'
        str_out <- structure(remstimateList, class = "remstimate")
        # defining attributes
        attr(str_out,"formula") <- model_formula
        attr(str_out,"model") <- model
        attr(str_out,"ordinal") <- ordinal
        attr(str_out, "method") <- method
        attr(str_out, "approach") <- "Frequentist"
        attr(str_out, "statistics") <- variable_names
        if(method %in% c("GD","GDADAM")){                 
            attr(str_out, "epochs") <- epochs
            attr(str_out, "learning rate") <- learning_rate
            attr(str_out, "iterations") <- optimum_obj$iterations
        }
    }
    else{ # bayesian approach
        # defining structure of class 'remstimate'
        str_out <- structure(remstimateList, class = "remstimate")
        # defining attributes
        attr(str_out,"formula") <- model_formula
        attr(str_out,"model") <- model
        attr(str_out,"ordinal") <- ordinal
        attr(str_out, "method") <- method
        attr(str_out, "approach") <- "Bayesian"
        attr(str_out, "statistics") <- variable_names
        attr(str_out, "prior") <- prior # i would like to save also the arguments inside (...), you can use match.args but you to save also the values
        attr(str_out, "nsim") <- nsim # from 'nsim' until 'thin' they will probably defined inside remstimateList
        attr(str_out, "seed") <- seed
        if(method == "HMC"){
            attr(str_out, "nchains") <- nchains
            attr(str_out, "burnin") <- burnin
            attr(str_out, "thin") <- thin 
        }
    }
    return(str_out)
}

# some notes: 
# when creating a method for instance remstimate.traceplot() if(attr(remstimate,"approach")!="Bayesian") stop("method traceplot() is not available for the frequentist approach.")
# if the attr(...,"approach") is changed to "Bayesian" but other bayesian attributes are not found, then print out stop("method traceplot() is available only for Bayeisna approach.")

# short version of summary (without pvalues,just estimates and null deviance residual deviance AIC AICC BIC and WAIC(A), ELPD(A), ELPD (A,PSIS))
# if A_SAP is 1 (default) then just write WAIC, ELPD and ELPD(PSIS) , otherwise write ELPD(A) and ELPD(A,PSIS)

#######################################################################################
#######################################################################################


# print.remstimate
#' @title print.remstimate
#' @rdname print.remstimate
#' @description A function that prints out the estimates returned by a 'remstimate' object.
#' @param x is a \code{remstimate} object 
#' @param ... further arguments to be passed to the print method
#' @method print remstimate
#' @export
print.remstimate<-function(x, ...){
    if (is.null(attr(x,"formula"))) 
        stop("invalid 'remstimate' object:  no 'formula' attribute")
    if (!inherits(x, "remstimate")) 
        warning("calling summary.lm(<fake-remstimate-object>) ...") # check this warning "summary.lm" or "print.lm"
    cat("Relational Event Model",paste("(",attr(x,"model")," oriented)",sep=""),"\n")
    if(attr(x,"approach") == "Frequentist"){
        cat("\nCoefficients:\n\n")
        print(x$coefficients)
    }
    else{ # Bayesian
        cat("\nPosterior Modes:\n\n")
        print(x$post.mode)
    }
    
    cat("\nNull deviance:",x$null.deviance,"\nResidual deviance:",x$residual.deviance,"\n")
    cat("AIC:",x$AIC,"AICC:",x$AICC,"BIC:",x$BIC,"\n\n")
}


#######################################################################################
#######################################################################################


# summary.remstimate
#' @title summary.remstimate
#' @rdname summary.remstimate
#' @description A function that arranges a summary of a 'remstimate' object
#' @param object is a \code{remstimate} object 
#' @param ... further arguments to be passed to the 'summary' method
#' @method summary remstimate
#' @export
summary.remstimate<-function (object, ...) #print.summary.remstimate
{
    # (codings are based on the structure of summary.lm() and summary.glm() from package 'stats')
    if (!inherits(object, "remstimate")) 
        warning("calling summary.remstimate(<fake-remstimate-object>) ...")

    summary_out <- list()
    if(attr(object, "approach") == "Frequentist"){ # Frequentist
        coefsTab <- cbind("Estimate" = object$coefficients,
        "Std. Err" = object$se,
        "z value" = object$coefficients/object$se,
        "Pr(>|z|)" = 2*(1-stats::pnorm(abs(object$coefficients/object$se)))
        )
        rownames(coefsTab) <- attr(object, "statistics")
        summary_out$coefsTab <- coefsTab

    }
    if(attr(object, "approach") == "Bayesian"){ # Bayesian
        coefsTab <- cbind("Post.Mode" = object$post.mode,
        "Post.SD" = object$sd,
        "Q2.5%" = apply(object$draws[[1]],2,stats::quantile,0.025),
        "Q50%" = apply(object$draws[[1]],2,stats::quantile,0.5),
        "Q97.5%" = apply(object$draws[[1]],2,stats::quantile,0.975)
        )
        rownames(coefsTab) <- attr(object, "statistics")
        summary_out$coefsTab <- coefsTab
    }

    keep <- match(c("formula","aic",
		      "contrasts", "df.residual","null.deviance","df.null",
                      "iter", "na.action"), names(object), 0L)

    keep <- list(formula = attr(object,"formula"),
                model = attr(object,"model"),
                ordinal = attr(object,"ordinal"),
                method = attr(object, "method"),
                approach = attr(object, "approach"),
                residual.deviance = object$residual.deviance,
                null.deviance = object$null.deviance,
                model.deviance = object$model.deviance,
                df.residual = (object$df.null - object$df.model),
                df.null = object$df.null,
                df.model = object$df.model,
                AIC = object$AIC,
                AICC = object$AICC,
                BIC = object$BIC,
                learning.rate = attr(object, "learning rate"))                
    summary_out <- do.call(c, list(keep, summary_out))

    if(length(summary_out)==0) stop("invalid 'remstimate' object")

    class(summary_out) <- "summary.remstimate"
    return(summary_out)
}


#######################################################################################
#######################################################################################


# print.summary.remstimate (the Bayesian part is missing in this method)
#' @title print.summary.remstimate
#' @rdname print.summary.remstimate
#' @description A function that prints out a summary of a 'remstimate' object.
#' @param x is a \code{summary.remstimate} object 
#' @param ... further arguments to be passed to the 'print.summary' method
#' @method print summary.remstimate
#' @export
print.summary.remstimate <- function(x, ...)
{
    cat("Relational Event Model",paste("(",x$model," oriented)",sep=""),"\n\n")
    if(x$model == "tie"){
        cat("Call:\n",deparse(x$formula),"\n\n",sep="")
    }
    else if(x$model == "actor"){
       cat("Call rate model:\n",deparse(x$formula$rate_model_formula),"\n\n",sep="")
       cat("Call choice model:\n",deparse(x$formula$choice_model_formula),"\n\n",sep="") 
    }
    second_line <- paste("(",x$method," with ",sep="")
    if(x$ordinal) second_line <- paste(second_line,"ordinal likelihood):",sep="")
    else{
        second_line <- paste(second_line,"interval likelihood):\n\n",sep="")
    }
    if(x$approach == "Frequentist"){
        cat("\nCoefficients",second_line)
    }
    else{ # Bayesian
        cat("\nPosterior Modes",second_line)
    }
    stats::printCoefmat(x$coefsTab, P.values = ifelse(x$approach == "Frequentist",TRUE,FALSE), ...)
    if(x$approach == "Frequentist"){
        cat("Null deviance:", x$null.deviance, "on", x$df.null, "degrees of freedom\n")
        cat("Residual deviance:", x$residual.deviance, "on", x$df.residual, "degrees of freedom\n")
        cat("Chi-square:", x$model.deviance, "on", x$df.model, 
            "degrees of freedom, asymptotic p-value", 1 - stats::pchisq(x$model.deviance, 
                x$df.model), "\n")
        cat("AIC:", x$AIC, "AICC:", x$AICC, "BIC:", x$BIC, "\n")
    }
    if(x$approach == "Bayesian"){
      cat("Log posterior:",x$logpost,"\n")
      cat("Prior parameters:",paste(names(x$prior.param),unlist(x$prior.param),sep="="),"\n")
    }
    #if(x$estimator%in%c("BMCMC","BSIR")){

    #}
}

#######################################################################################
#######################################################################################


# predict.remstimate
#' @title predict.remstimate
#' @rdname predict.remstimate
#' @description A function that returns predictions given a 'remstimate' object.
#' @param object is a \code{remstimate} object 
#' @param ... further arguments to be passed to the 'predict' method depending on the estimator used
#' @method predict remstimate
#' @export
predict.remstimate <- function(object, ...)
{
 # write a predict method here
 stop('this function at the moment does nothing')
}


#######################################################################################
#######################################################################################


# plot.remstimate
#' @title plot.remstimate
#' @rdname plot.remstimate
#' @description A function that returns a plot of diagnostics given a 'remstimate' object and depending on the 'approach' attribute.
#' @param x is a \code{remstimate} object 
#' @param ... further arguments to be passed to the 'predict' method depending on the estimator used
#' @method plot remstimate
#' @export
plot.remstimate <- function(x, ...)
{
 # write a plot method here
 stop('this function at the moment does nothing')
}


#######################################################################################
#######################################################################################

# aic
#' @title aic
#' @description A function that returns the AIC (Akaike's Information Criterion) value in a 'remstimate' object
#' @param object is a \code{remstimate} object 
#' @param ... further arguments to be passed to the 'aic' method
#' @export
aic <- function(object,...){
  UseMethod("aic")
}


#' @describeIn aic AIC (Akaike's Information Criterion) value of a 'remstimate' object
#' @method aic remstimate
#' @export
aic.remstimate <- function(object,...) {if(attr(object, "approach") == "Frequentist"){object$AIC}else{stop("'approach' must be 'Frequentist'")}}


#######################################################################################
#######################################################################################


# aicc
#' @title aicc 
#' @description A function that returns the AICC (Akaike's Information Corrected Criterion) value in a 'remstimate' object
#' @param object is a \code{remstimate} object 
#' @param ... further arguments to be passed to the 'aicc' method
#' @export
aicc <- function(object,...){
  UseMethod("aicc")
}


#' @describeIn aicc AICC (Akaike's Information Corrected Criterion) value of a 'remstimate' object
#' @method aicc remstimate
#' @export
aicc.remstimate <- function(object,...) {if(attr(object, "approach") == "Frequentist"){object$AICC}else{stop("'approach' must be 'Frequentist'")}}


#######################################################################################
#######################################################################################


# bic
#' @title bic
#' @description A function that returns the BIC (Bayesian Information Criterion) value in a 'remstimate' object
#' @param object is a \code{remstimate} object 
#' @param ... further arguments to be passed to the 'bic' method
#' @export
bic <- function(object,...){
  UseMethod("bic")
}


#' @describeIn bic BIC (Bayesian Information Criterion) value of a 'remstimate' object
#' @method bic remstimate
#' @export
bic.remstimate <- function(object,...) {if(attr(object, "approach") == "Frequentist"){object$BIC}else{stop("'approach' must be 'Frequentist'")}}


#######################################################################################
#######################################################################################


# waic
#' @title waic
#' @description A function that returns the WAIC (Watanabe-Akaike's Information Criterion) value in a 'remstimate' object
#' @param object is a \code{remstimate} object 
#' @param ... further arguments to be passed to the 'waic' method
#' @export
waic <- function(object,...){
  UseMethod("waic")
}


#' @describeIn waic WAIC (Watanabe-Akaike's Information Criterion) value of a 'remstimate' object
#' @method waic remstimate
#' @export
waic.remstimate <- function(object,...) {if(attr(object, "approach") == "Frequentist"){if(!is.null(object$WAIC)){object$WAIC}else{stop("'WAIC' value not found.")}}else{stop("'approach' must be 'Frequentist'")}}


#######################################################################################
#######################################################################################
