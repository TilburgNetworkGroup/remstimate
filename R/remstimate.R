#' remstimate  
#'
#' A function for the optimization of tie-oriented (Relational Event Model) or actor-oriented (DyNAM) likelihoods
#'
#' @param reh a processed relational event history (REH). It can be a \code{remify} object (output object of the function \code{remify::remify()}) or the event sequence (as \code{data.frame})
#' @param stats a \code{remstats} object: when \code{model="tie"}, \code{stats} is an array of statistics with dimensions \code{[M x D x P]}: where \code{M} is the number of events, \code{D} is the number of possible dyads (full riskset), \code{P} is the number of statistics; if \code{model="actor"}, \code{stats} is a list of two arrays named \code{rate} and \code{choice} with dimensions \code{[M x N x P]}, where \code{N} are the actors (senders in the array \code{rate}, receivers in the array \code{choice})
#' @param method optimization method to be used. Methods available are: Maximum Likelihood Estimation (\code{"MLE"}), Adaptive Gradient Descent (\code{"GDADAMAX"}), Bayesian Samplin Importance Resampling (\code{"BSIR"}), Hamiltonian Monte Carlo (\code{"HMC"})
#' @param ncores number of threads for the parallelization (default value is \code{ncores = 1}, that means no parallelization)
#' @param prior prior distribution when using the \code{"BSIR"} method
#' @param nsim  when \code{method = "HMC"} it is the number of simulations (iterations) in each chain, when \code{method = "BSIR"} is the number of samples from the proposal distribution
#' @param nchains number of chains to generate in the case of \code{method = "HMC"}
#' @param burnin number of initial iterations to be added as burnin (for \code{method = "HMC"})
#' @param thin number of steps to skip in the posterior draws (for \code{method = "HMC"})
#' @param init vector of initial values if tie-oriented model, or a named list of two vectors ('rate' and 'choice') if actor-oriented model. This argument is used for the methods: \code{"GDADAMAX"} and \code{"HMC"} 
#' @param epochs 1e03 by defaut. It is the number of iteration used in the method \code{"GDADAMAX"}
#' @param epsilon 0.001 by default. It is the inter-iteration difference of the loss function used in the method \code{"GDADAMAX"} and it is used as stop-rule within the algorithm.
#' @param seed seed for reproducibility [[yet to be integrated in the code]]
#' @param silent a \code{TRUE/FALSE} value. If \code{FALSE}, progress of optimization status will be printed out [[yet to be integrated in the code]]
#' @param ... additional parameters. They can be parameters of other functions defined as input in some of the arguments above
#'
#' @return  'remstimate' S3 object
#' @export
remstimate <- function(reh,
                       stats, 
                       method = c("MLE","GDADAMAX","BSIR","HMC"),
                       ncores = 1L,
                       prior = NULL,
                       nsim = 1e03L,
                       nchains = 1L,
                       burnin = 500L,
                       thin = 10L,  
                       init = NULL,                 
                       epochs = 1e03L,
                       epsilon = 0.001,
                       seed = sample(1:1e04,nchains),
                       silent = TRUE,
                       ...){

    # ... processing input:

    # ... reh
    if(!inherits(reh,"remify")){
        if(is.data.frame(reh)){
            reh <- remify::remify(edgelist = reh, ...)
        }
        else{
            stop("'reh' must be either a 'remify' object (from the 'remify' package) or the event sequence (as data.frame).")
        }
    }

    # ... model :
    model <- attr(reh, "model") # attribute from reh object

    # ... method : 
    if(!(method %in% c("MLE","GDADAMAX","BSIR","HMC"))){stop("The `method` specified is not available or it is mistyped.")}

    # ... type of likelihood
    ordinal <- attr(reh,"ordinal")

    # ... directed / undirected network
    if((model == "actor") & !attr(reh,"directed")){stop("actor-oriented modeling can't operate on undirected networks")}

    # ... stats 
    # check if 'stats' is an object of class 'remstats'
    model_formula <- variable_names <- NULL
    if(model == "tie")
    {
        if(inherits(stats,c("remstats","tomstats"))){
            variable_names <- as.vector(sapply(dimnames(stats$statistics)[[3]],function(x) sub(pattern = ".x.", replacement = ":", x = x)))
            model_formula <- stats::as.formula(paste("~ ",paste(variable_names,collapse=" + ")))
            stats <- aperm(stats$statistics, perm = c(2,3,1)) # stats reshaped in [D*U*M] (this might change in the future?)
        }
        else if(inherits(stats,c("remstats","aomstats"))){
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

        if((!is.null(init)) & (method == "GDADAMAX")){
            if(!is.vector(init)){ # there should be more stop()'s as to the structure of the init input based on the modeling framework (tie or actor oriented)
                stop("'init' must be a vector with the starting values of the effects of the statistics")
            }
        }
    }
    if(model == "actor") 
    {
        if(inherits(stats,c("remstats","aomstats"))){
            variables_rate <- as.vector(sapply(dimnames(stats$statistics$sender_stats)[[3]],function(x) sub(pattern = ".x.", replacement = ":", x = x)))
            variables_choice <- as.vector(sapply(dimnames(stats$statistics$receiver_stats)[[3]],function(x) sub(pattern = ".x.", replacement = ":", x = x)))
            variable_names <- c(variables_rate,variables_choice)
            model_formula <- list(rate_model_formula = stats::as.formula(paste("~ ",paste(variables_rate,collapse=" + "))), choice_model_formula = stats::as.formula(paste("~ ",paste(variables_choice,collapse=" + "))))
            stats <- list(sender_rate = aperm(stats$statistics$sender_stats, perm = c(2,3,1)), receiver_choice = aperm(stats$statistics$receiver_stats, perm = c(2,3,1))) # stats reshaped in [D*U*M] (this might change in the future?)
        }
        else if(inherits(stats,c("remstats","tomstats"))){
            stop("'remstats' object supplied cannot work for tie-oriented modeling")
        }
        else if(is.list(stats) & all(c("sender_rate","receiver_choice") %in% names(stats))){
             M_rate <- dim(stats$sender_rate)[1]
             N_rate <- dim(stats$sender_rate)[2]
             M_choice <- dim(stats$receiver_choice)[1]
             N_choice <- dim(stats$receiver_choice)[2]
             if((M_rate != reh$M) & (M_choice != reh$M) & (N_rate != reh$N) & (N_choice != reh$N)){
                 stop("numbers of rows and columns in both arrays ('sender_rate' and 'receiver_choice') inside 'stats' must be equal to number of 'time points' and 'actors' in the network (see dim(reh))")
             }
             else{
                 stats <- list(sender_rate = aperm(stats$sender_rate, perm = c(2,3,1)), choice = aperm(stats$receiver_choice, perm = c(2,3,1)))
             }
            model_formula <- stats::as.formula(".~.")
        }
        else{
            stop("actor-oriented modeling: 'stats' must be either a 'aomstats' 'remstats' object or a list of two arrays named 'rate' and 'choice'")
        }

        if((!is.null(init)) & (method == "GDADAMAX")){
            if(!is.list(init)){
                    stop("'init' must be a list of two vectors named 'sender_rate' and 'receiver_choice', each one with the starting values of the effects of the statistics according to the two different models (rate model and choice model)")
                }
            else{
                if(!all(c("sender_rate","receiver_choice") %in% names(init))){
                    stop("'init' must be a list of two vectors named 'sender_rate' and 'receiver_choice', each one with the starting values of the effects of the statistics according to the two different models (rate model and choice model)")
                }
            }
        }
    }

    # ... epochs and epsilon (parameters for GDADAMAX)
    if(method == "GDADAMAX"){

        # ... epochs
        if(is.null(epochs)){
            epochs <- 1e03L
        }
        else if(is.numeric(epochs) | is.integer(epochs)){
            epochs <- as.integer(epochs) # this converts, for instance 3.2 to 3, 3.7 to 3
            if(epochs<0){
                stop("'epoch' must be a positive number'")
            }
        }
        else{
            stop("'epoch' must be a positive number")
        }


        # ... epsilon
        if(is.null(epsilon)){
            epsilon <- 0.01
        }
        else if((epsilon <= 0)){
            stop("'epsilon' must be a positive number")
        }
    }

    # ... prior
    log <- TRUE
    if(!is.null(prior)){ # [[ yet to be implemented in the code]]
        additional_input_args <- names(list(...)) #names(as.list(match.call()))[-1]
        args_prior <- match.arg(arg=additional_input_args,choices = methods::formalArgs(prior),several.ok = TRUE)
    }


    # ... ncores
    if(is.null(ncores)) ncores <- 1L
    else{
        if((parallel::detectCores() == 2L) & (ncores > 1L))
            stop("'ncores' is recommended to be set at most to 1.")
        else if((parallel::detectCores() > 2L) & (ncores > floor(parallel::detectCores()-2L)))
            stop("'ncores' is recommended to be set at most to: floor(parallel::detectCores()-2L)")
    }

    # ... seed
    if(method == "HMC"){
        if(is.null(seed)){
            if(!is.null(nchains)){
                seed <- sample(1:1e06,nchains)
            }
            else{
                nchains <- 2L
                seed <- sample(1:1e06,nchains)
            }
        }
        else{
            if(is.null(nchains)){
                nchains <- length(seed)
            }
            else{
                if(length(seed) != nchains)
                    stop("the number of chains (`nchains`) must be equal to the number of seeds (`seed`) supplied")
            }
        }
        
    }
    if(method == "BSIR"){
        if(is.null(seed)){
            seed <- sample(1:1e06,nchains)
        }
        else if(length(seed) > 1){
            seed <- seed[1]
            warning("`seed` length is greater than 1. Considering only the first element")
        }
    }

    # ... nsim
    if(is.null(nsim) & (method %in% c("BSIR","HMC"))){
        nsim <- 1e03L
    }

    # ... burnin
    if(is.null(burnin) & (method == "HMC")){
        burnin <- 500L
    }
    else if(method == "HMC"){
        if(burnin > nsim)
            stop("`burnin` value must be lower than the number of simulations (`nsim`)")
        
    }

    # ... thin
    if(is.null(thin) & (method == "HMC")){
        if((burnin == 500L) & (nsim == 1e03L)){ # if the burnin and the nsim value are equal to the default ones then we set thin to its default value
            thin <- 10L
        }
        else{
            # we set a value of thin that is works fine given the nsim and the burnin value
            if((nsim-burnin)<500L){
                thin <- 5L # if the sequences without burnin are shorter than 500, we reduce the thinning, this is arbitrary only when thin argument is NULL
            }
            else{
                thin <- 10L
            }
        }
    }
    else if(!is.null(thin) & (method == "HMC")){
        if(thin <= 0){
            stop("`thin` value must be positive. Set thin = 1 if no thinning of chains is required")
        }
    }


    # ... creating  an empty lists
    remstimateList <- list()
        
    # ... [1] Maximum Likelihood, Gradient Descent ADAMAX (GDADAMAX)
    if(method %in% c("MLE","GDADAMAX")){
        if(model == "tie"){ # Relational Event Model (REM)
            if(method == "MLE"){
                optimum_obj <- trust::trust(objfun = remDerivatives, 
                                            parinit = rep(0,dim(stats)[2]), 
                                            rinit = 1, 
                                            rmax = 100, 
                                            stats = stats,  
                                            actor1 = c(0),
                                            actor2 = c(0),
                                            dyad = attr(reh,"dyad")-1,
                                            omit_dyad = reh$omit_dyad,
                                            interevent_time = reh$intereventTime,
                                            model = model,
                                            ordinal = ordinal,
                                            ncores = ncores)
            }  
            if(method == "GDADAMAX"){
                if(is.null(init)){
                    init <- rep(0,dim(stats)[2])
                }
                optimum_obj <- GDADAMAX(pars = init,
                                    stats = stats,  
                                    actor1 = c(0),
                                    actor2 = c(0),
                                    dyad = attr(reh,"dyad")-1,
                                    omit_dyad = reh$omit_dyad,
                                    interevent_time = reh$intereventTime,
                                    model = model,
                                    ordinal = ordinal,
                                    ncores = ncores,
                                    epochs = epochs,
                                    epsilon = epsilon)
                optimum_obj$hessian <- remstimate::remDerivativesStandard(par = optimum_obj$argument,
                                                                    stats = stats,  
                                                                    actor1 = c(0),
                                                                    actor2 = c(0),
                                                                    dyad = attr(reh,"dyad")-1,
                                                                    omit_dyad = reh$omit_dyad,
                                                                    interevent_time = reh$intereventTime,
                                                                    ordinal = ordinal,
                                                                    ncores = ncores)
                optimum_obj$hessian <- optimum_obj$hessian$hessian                                                  
            }                                     
            remstimateList$coefficients <- optimum_obj$argument 
            remstimateList$loglik <- -optimum_obj$value # loglikelihood value at MLE values (we take the '-value' because we minimized the '-loglik')        
            remstimateList$gradient <- -optimum_obj$gradient
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
                                        actor1 = c(0),
                                        actor2 = c(0),
                                        dyad = attr(reh,"dyad")-1,
                                        omit_dyad = reh$omit_dyad,
                                        interevent_time = reh$intereventTime,
                                        model = model,
                                        ordinal = ordinal,
                                        ncores = ncores)$value)
            remstimateList$model.deviance <- remstimateList$null.deviance -  remstimateList$residual.deviance                     
            remstimateList$df.null <- reh$M
            remstimateList$df.model <- dim(stats)[2]

            remstimateList$AIC <- 2*length(remstimateList$coefficients) - 2*remstimateList$loglik # AIC
            remstimateList$AICC <- remstimateList$AIC + 2*length(remstimateList$coefficients)*(length(remstimateList$coefficients)+1)/(reh$M-length(remstimateList$coefficients)-1)
            remstimateList$BIC <- length(remstimateList$coefficients)*log(reh$M) - 2*remstimateList$loglik # BIC
            remstimateList$converged <- optimum_obj$converged
            remstimateList$iterations <- optimum_obj$iterations                                    
        }
        if(model == "actor" ){ # Actor Oriented Model 
            if(method == "MLE"){
                optimum_sender_rate <- trust::trust(objfun = remDerivatives, 
                                    parinit = rep(0,dim(stats$sender_rate)[2]), 
                                    rinit = 1, 
                                    rmax = 100, 
                                    stats = stats$sender_rate, 
                                    actor1 = reh$edgelist$actor1-1,
                                    actor2 = c(0),
                                    dyad = c(0),
                                    omit_dyad = reh$omit_dyad,
                                    interevent_time = reh$intereventTime,
                                    model = model,
                                    ordinal = ordinal,
                                    senderRate = TRUE,
                                    C = reh$C,
                                    D = reh$D)      
                optimum_receiver_choice <- trust::trust(objfun = remDerivatives, 
                                    parinit = rep(0,dim(stats$receiver_choice)[2]), 
                                    rinit = 1, 
                                    rmax = 100, 
                                    stats = stats$receiver_choice, 
                                    actor1 = reh$edgelist$actor1-1,
                                    actor2 = reh$edgelist$actor2-1,
                                    dyad = c(0),
                                    omit_dyad = reh$omit_dyad, 
                                    interevent_time = reh$intereventTime,
                                    model = model,
                                    senderRate = FALSE,
                                    N = reh$N,
                                    C = reh$C,
                                    D = reh$D)
            }
            if(method == "GDADAMAX"){
                if(is.null(init)){
                    init <- list("sender_rate" = rep(0.0,dim(stats$sender_rate)[2]),"receiver_choice" = rep(0.0,dim(stats$receiver_choice)[2]))
                }
                optimum_sender_rate <- GDADAMAX(pars = init$sender_rate, 
                                        stats = stats$sender_rate,  
                                        actor1 = reh$edgelist$actor1-1,
                                        actor2 = c(0),
                                        dyad = c(0),
                                        omit_dyad = reh$omit_dyad,
                                        interevent_time = reh$intereventTime,
                                        model = model,
                                        ordinal = ordinal,
                                        senderRate = TRUE,
                                        C = reh$C,
                                        D = reh$D,
                                        ncores = ncores,
                                        epochs = epochs,
                                        epsilon = epsilon)                    
                       
                optimum_receiver_choice <- GDADAMAX(pars = init$receiver_choice,
                                        stats = stats$receiver_choice,  
                                        actor1 = reh$edgelist$actor1-1,
                                        actor2 = reh$edgelist$actor2-1,
                                        dyad = c(0),
                                        omit_dyad = reh$omit_dyad,
                                        interevent_time = reh$intereventTime,
                                        model = model,
                                        senderRate = FALSE,
                                        N = reh$N,
                                        C = reh$C,
                                        D = reh$D,
                                        ncores = ncores,
                                        epochs = epochs,
                                        epsilon = epsilon)
                # calculating hessian for the sender rate model
                optimum_sender_rate$hessian <- remstimate::remDerivativesSenderRates(par = optimum_sender_rate$argument,
                                                                    stats = stats$sender_rate,  
                                                                    actor1 = reh$edgelist$actor1-1,
                                                                    actor2 = c(0),
                                                                    dyad = c(0),
                                                                    omit_dyad = reh$omit_dyad,
                                                                    interevent_time = reh$intereventTime,
                                                                    C = reh$C,
                                                                    D = reh$D,
                                                                    ordinal = ordinal,
                                                                    hessian = TRUE)
                optimum_sender_rate$hessian <- optimum_sender_rate$hessian$hessian    

                # calculating hessian for the receiver choice model
                optimum_receiver_choice$hessian <- remstimate::remDerivativesReceiverChoice(par = optimum_sender_rate$argument,
                                                                    stats = stats$receiver_choice,  
                                                                    actor1 = reh$edgelist$actor1-1,
                                                                    actor2 = reh$edgelist$actor2-1,
                                                                    dyad = c(0),
                                                                    omit_dyad = reh$omit_dyad,
                                                                    interevent_time = reh$intereventTime,
                                                                    N = reh$N,
                                                                    C = reh$C,
                                                                    D = reh$D,
                                                                    hessian = TRUE)
                optimum_receiver_choice$hessian <- optimum_receiver_choice$hessian$hessian                         
            }

            remstimateList$sender_rate <- remstimateList$receiver_choice <- list() 

            # coefficients
            remstimateList$sender_rate$coefficients <- optimum_sender_rate$argument
            remstimateList$receiver_choice$coefficients <- optimum_receiver_choice$argument
            # loglik
            remstimateList$sender_rate$loglik <- -optimum_sender_rate$value # log(L_sender)
            remstimateList$receiver_choice$loglik <- -optimum_receiver_choice$value # log(L_choice)  
            # gradient     
            remstimateList$sender_rate$gradient <- optimum_sender_rate$gradient
            remstimateList$receiver_choice$gradient <- optimum_receiver_choice$gradient
            # hessian
            remstimateList$sender_rate$hessian <- optimum_sender_rate$hessian # hessian matrix for sender rate model
            remstimateList$receiver_choice$hessian <- optimum_receiver_choice$hessian # hessian matrix for choice model
            # vcov
            remstimateList$sender_rate$vcov <- qr.solve(remstimateList$sender_rate$hessian) # matrix of variances and covariances for the sender rate model
            remstimateList$receiver_choice$vcov <- qr.solve(remstimateList$receiver_choice$hessian) # matrix of variances and covariances for the receiver choice model
            # standard errors
            remstimateList$sender_rate$se <- diag(remstimateList$sender_rate$vcov)**0.5 # standard errors
            remstimateList$receiver_choice$se <- diag(remstimateList$receiver_choice$vcov)**0.5 # standard errors           
            #re-naming
            names(remstimateList$sender_rate$coefficients) <- names(remstimateList$sender_rate$se) <- rownames(remstimateList$sender_rate$vcov) <- colnames(remstimateList$sender_rate$vcov) <- colnames(remstimateList$sender_rate$hessian) <- rownames(remstimateList$sender_rate$hessian) <- dimnames(stats$sender_rate)[[2]]
            names(remstimateList$receiver_choice$coefficients) <- names(remstimateList$receiver_choice$se) <- rownames(remstimateList$receiver_choice$vcov) <- colnames(remstimateList$receiver_choice$vcov) <- colnames(remstimateList$receiver_choice$hessian) <- rownames(remstimateList$receiver_choice$hessian) <- dimnames(stats$receiver_choice)[[2]]
            # residual deviance
            remstimateList$sender_rate$residual.deviance <- -2*remstimateList$sender_rate$loglik
            remstimateList$sender_rate$residual.deviance <- -2*remstimateList$sender_rate$loglik
            # null deviance
            remstimateList$sender_rate$null.deviance <- 2*(trust::trust(objfun = remDerivatives, 
                                                                    parinit = c(0), 
                                                                    rinit = 1, 
                                                                    rmax = 100, 
                                                                    stats = array(1,dim=c(dim(stats$sender_rate)[1],1,reh$M)),  
                                                                    actor1 = reh$edgelist$actor1-1,
                                                                    actor2 = c(0),
                                                                    dyad = c(0),
                                                                    omit_dyad = reh$omit_dyad,
                                                                    interevent_time = reh$intereventTime,
                                                                    model = model,
                                                                    C = reh$C,
                                                                    D = reh$D,
                                                                    ordinal = ordinal,
                                                                    senderRate = TRUE)$value)
                                                     
            remstimateList$receiver_choice$null.deviance <- ifelse(length(reh$omit_dyad)>0,2*(trust::trust(objfun = remDerivatives, 
                                                                parinit = c(0), 
                                                                rinit = 1, 
                                                                rmax = 100, 
                                                                stats = array(1,dim=c(dim(stats$receiver_choice)[1],1,reh$M)), 
                                                                actor1 = reh$edgelist$actor1-1,
                                                                actor2 = reh$edgelist$actor2-1,
                                                                dyad = c(0),
                                                                omit_dyad = reh$omit_dyad, 
                                                                interevent_time = reh$intereventTime,
                                                                model = model,
                                                                senderRate = FALSE,
                                                                N = reh$N,
                                                                C = reh$C,
                                                                D = reh$D)$value),-2*log(1/(reh$N-1)))
            # null.deviance.choice is -2*log(1/(reh$N-1)) when the riskset is always (all the time points) the same
            # model deviance
            remstimateList$sender_rate$model.deviance <- remstimateList$sender_rate$null.deviance -  remstimateList$sender_rate$residual.deviance 
            remstimateList$receiver_choice$model.deviance <- remstimateList$receiver_choice$null.deviance - remstimateList$receiver_choice$residual.deviance
            # df null
            remstimateList$sender_rate$df.null <- reh$M
            remstimateList$receiver_choice$df.null <- reh$M
            # df model
            remstimateList$sender_rate$df.model <- dim(stats$sender_rate)[2]
            remstimateList$receiver_choice$df.model <- dim(stats$receiver_choice)[2]
            # df residual
            remstimateList$sender_rate$df.residual <- remstimateList$sender_rate$df.null - remstimateList$sender_rate$df.model
            remstimateList$receiver_choice$df.residual <- remstimateList$receiver_choice$df.null - remstimateList$receiver_choice$df.model       
            # AIC
            remstimateList$sender_rate$AIC <- 2*length(remstimateList$sender_rate$coefficients) - 2*remstimateList$sender_rate$loglik 
            remstimateList$receiver_choice$AIC <- 2*length(remstimateList$receiver_choice$coefficients) - 2*remstimateList$receiver_choice$loglik 
            # AICC
            remstimateList$sender_rate$AICC <- remstimateList$sender_rate$AIC + 2*length(remstimateList$sender_rate$coefficients)*(length(remstimateList$sender_rate$coefficients)+1)/(reh$M-length(remstimateList$sender_rate$coefficients)-1)
            remstimateList$receiver_choice$AICC <- remstimateList$receiver_choice$AIC + 2*length(remstimateList$receiver_choice$coefficients)*(length(remstimateList$receiver_choice$coefficients)+1)/(reh$M-length(remstimateList$receiver_choice$coefficients)-1)
            # BIC
            remstimateList$sender_rate$BIC <- length(remstimateList$sender_rate$coefficients)*log(reh$M) - 2*remstimateList$sender_rate$loglik
            remstimateList$receiver_choice$BIC <- length(remstimateList$receiver_choice$coefficients)*log(reh$M) - 2*remstimateList$receiver_choice$loglik
            # converged
            remstimateList$sender_rate$converged <- optimum_sender_rate$converged
            remstimateList$receiver_choice$converged <- optimum_receiver_choice$converged
            # iterations
            remstimateList$sender_rate$iterations <- optimum_sender_rate$iterations
            remstimateList$receiver_choice$iterations <- optimum_receiver_choice$iterations
        }            
        # ...
        # in future : if(WAIC) and if(ELPD){if(PSIS)} for predictive power of models
    }
   
    # ... [4] with Bayesian Sampling Importance Resampling (BSIR)
    if(method == "BSIR"){

        if(model == "tie"){ # Relational Event Model (REM)

            bsir <- list()
            bsir$log_posterior <- rep(0,nsim)  # 
            bsir$draws <- list()  

            # (1) generate from the proposal (importance distribution)

            # First find location and scale parameters (mles and variances and covariances matrix)
            mle_optimum <- trust::trust(objfun = remDerivatives, 
                                        parinit = rep(0,dim(stats)[2]), 
                                        rinit = 1, 
                                        rmax = 100, 
                                        stats = stats,  
                                        actor1 = c(0),
                                        actor2 = c(0),
                                        dyad = attr(reh,"dyad")-1,
                                        omit_dyad = reh$omit_dyad,
                                        interevent_time = reh$intereventTime,
                                        model = model,
                                        ordinal = ordinal,
                                        ncores = ncores)
            bsir$mle <- mle_optimum$argument
            bsir$vcov <- qr.solve(mle_optimum$hessian)
            
            # proposal distribution (default proposal is a multivariate Student t): drawing nsim*3 samples
            bsir$draws <- mvnfast::rmvt(n = nsim*3, 
                                            mu = bsir$mle, 
                                            sigma = bsir$vcov,
                                            df = 4,
                                            ncores = ncores) # df = 4 by default
                                        
            bsir$log_proposal <- mvnfast::dmvt(X = bsir$draws,  
                                                mu = bsir$mle, 
                                                sigma = bsir$vcov,
                                                df = 4,
                                                log = TRUE,
                                                ncores = ncores)

            # (2) calculate log-posterior density give by log_prior + loglik
            ## (2.1) summing log_prior (if NULL it will be non-informative) : prior by default can be a multivariate Student t and the user must specify its parameters
            if(!is.null(prior)){
                bsir$log_prior <- prior(bsir$draws,...) 
                bsir$log_posterior <- bsir$log_prior
            }
            ## (2.2) summing loglik                   
            loglik <- apply(bsir$draws,1,function(x) (-1)*remDerivativesStandard(pars = x,
                                    stats = stats, 
                                    actor1 = c(0),
                                    actor2 = c(0),
                                    dyad = attr(reh,"dyad")-1,
                                    omit_dyad = reh$omit_dyad,
                                    interevent_time = reh$intereventTime,
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

            post_draws <-sample(x=1:dim(bsir$draws)[1],size=nsim,replace=TRUE,prob=exp(irlw)) # we always draw nsim (having generated 3*nsim samples from the proposal)
            bsir$draws <- bsir$draws[post_draws,]
            bsir$log_posterior <- bsir$log_posterior[post_draws]
            bsir$coefficients <- bsir$draws[which.max(bsir$log_posterior),]
            bsir$loglik <- max(bsir$log_posterior)
            bsir$post.mean <- colMeans(bsir$draws)
            bsir$vcov <- stats::cov(bsir$draws)
            bsir$sd <- diag(bsir$vcov)**0.5
            bsir$df.null <- reh$M
            names(bsir$coefficients) <- names(bsir$post.mean) <- rownames(bsir$vcov) <- colnames(bsir$vcov) <- names(bsir$sd) <- dimnames(stats)[[2]]
            
            remstimateList <- bsir
        }
        if(model == "actor"){ # Actor Oriented Model

            bsir <- list()
            senderRate <- c(TRUE,FALSE) # meaning that the first BSIR will be run on the "sender rate" model and the second on the "receiver choice"
            which_model <- c("sender_rate","receiver_choice")
            
            for(i in 1:2){
                bsir_i <- list()
                bsir_i$log_posterior <- rep(0,nsim)
                bsir_i$draws <- list()  

                # (1) generate from the proposal (importance distribution)

                # First find location and scale parameters (mles and variances and covariances matrix)
                mle_optimum <- trust::trust(objfun = remDerivatives, 
                                                    parinit = rep(0,dim(stats[[which_model[i]]])[2]), 
                                                    rinit = 1, 
                                                    rmax = 100, 
                                                    stats = stats[[which_model[i]]], 
                                                    actor1 = reh$edgelist$actor1-1,
                                                    actor2 = reh$edgelist$actor2-1,
                                                    dyad = c(0),
                                                    omit_dyad = reh$omit_dyad,
                                                    interevent_time = reh$intereventTime,
                                                    model = model,
                                                    ordinal = ordinal,
                                                    ncores = ncores,
                                                    senderRate = senderRate[i],
                                                    N = reh$N,
                                                    C = reh$C,
                                                    D = reh$D)                          
                bsir_i$mle <- mle_optimum$argument
                bsir_i$vcov <- qr.solve(mle_optimum$hessian)
                
                # proposal distribution (default proposal is a multivariate Student t): drawing nsim*3 samples
                bsir_i$draws <- mvnfast::rmvt(n = nsim*3, 
                                                mu = bsir_i$mle, 
                                                sigma = bsir_i$vcov,
                                                df = 4,
                                                ncores = ncores) # df = 4 by default
                                            
                bsir_i$log_proposal <- mvnfast::dmvt(X = bsir_i$draws,  
                                                    mu = bsir_i$mle, 
                                                    sigma = bsir_i$vcov,
                                                    df = 4,
                                                    log = TRUE,
                                                    ncores = ncores)

                # (2) calculate log-posterior density give by log_prior + loglik
                ## (2.1) summing log_prior (if NULL it will be non-informative) : prior by default can be a multivariate Student t and the user must specify its parameters
                if(!is.null(prior)){
                    bsir_i$log_prior <- prior[[which_model[i]]](bsir_i$draws,...) 
                    bsir_i$log_posterior <- bsir_i$log_prior
                }
                ## (2.2) summing loglik                   
                loglik <- apply(bsir_i$draws,1,function(x) (-1)*remDerivatives(pars = x, 
                                stats = stats[[which_model[i]]], 
                                actor1 = reh$edgelist$actor1-1,
                                actor2 = reh$edgelist$actor2-1,
                                dyad = c(0),
                                omit_dyad = reh$omit_dyad,
                                interevent_time = reh$intereventTime,
                                model = model,
                                ordinal = ordinal,
                                senderRate = senderRate[i],
                                ncores = ncores,
                                gradient = FALSE,
                                hessian = FALSE,
                                N = reh$N,
                                C = reh$C,
                                D = reh$D)$value) # this is the most time consuming step to be optimized, avoiding the apply        
                bsir_i$log_posterior <- bsir_i$log_posterior + loglik

                # (3) calculate Importance resampling log-weights (irlw)
                irlw <- (bsir_i$log_posterior - max(bsir_i$log_posterior) + min(bsir_i$log_proposal)) - bsir_i$log_proposal
                s_irlw <- log(sum(exp(irlw)) - exp(irlw)) # ISIR (Skare et al., 2003 Scandinavian Journal of Statistics)
                irlw <- irlw - s_irlw
                bsir_i$irw <- exp(irlw) # importance resampling weights (irw)
                #this is a check from relevent::rem() the function after this message simply draws according to the weights without any particular calculation if the warning is printed out
                #if(any(is.na(iw)|is.nan(iw)|is.na(exp(iw)))||all(exp(iw)==0)){
                #warning("Importance weights seem to have some issues.  Carrying on as best we can.\n")
                #}
                post_draws <-sample(x=1:dim(bsir_i$draws)[1],size=nsim,replace=TRUE,prob=exp(irlw)) # we always draw nsim (having generated 3*nsim samples from the proposal)
                bsir_i$draws <- bsir_i$draws[post_draws,]
                bsir_i$log_posterior <- bsir_i$log_posterior[post_draws]
                bsir_i$coefficients <- bsir_i$draws[which.max(bsir_i$log_posterior),]
                bsir_i$loglik <- max(bsir_i$log_posterior)
                bsir_i$post.mean <- colMeans(bsir_i$draws)
                bsir_i$vcov <- stats::cov(bsir_i$draws)
                bsir_i$sd <- diag(bsir_i$vcov)**0.5
                bsir_i$df.null <- reh$M
                names(bsir_i$coefficients) <- names(bsir_i$post.mean) <- rownames(bsir_i$vcov) <- colnames(bsir_i$vcov) <- names(bsir_i$sd) <- dimnames(stats[[which_model[i]]])[[2]]
                
                bsir[[which_model[i]]] <- bsir_i

                rm(bsir_i,mle_optimum,loglik,irlw,s_irlw,post_draws)
            }
            remstimateList <- bsir
        }
    }

    # ... [5] Hamiltonian Monte Carlo (HMC)
    if(method == "HMC"){

        # create an empty list where to store the output
        hmc <- list()

        if(model == "tie"){ # Relational Event Model (REM)
            if(is.null(init)){
                beta_0 <- log(reh$M) - log(sum(reh$intereventTime)*reh$D) # MLE intercept under null model
                init <- matrix(stats::runif((dim(stats)[2])*nchains,-0.1,0.1),nrow=dim(stats)[2],ncol=nchains) # runif was in (-0.1,0.1)
                if(.GlobalEnv$with_intercept_hmc){
                    init[1,] <- init[1,] + beta_0
                }
                else{
                    init <- init + beta_0
                }
                #init <- matrix(rep(.GlobalEnv$mle_hmc,nchains),nrow=dim(stats)[2],ncol=nchains,byrow=FALSE) + matrix(runif(dim(stats)[2]*nchains,-1,1),nrow=dim(stats)[2],ncol=nchains,byrow=FALSE)
                #init[,1] <- .GlobalEnv$mle_hmc
            }
            .GlobalEnv$init_hmc <- init # temporary

            hmc_out <- HMC(pars_init = init, 
                        nsim = (burnin+nsim), 
                        nchains = nchains, 
                        burnin = burnin, 
                        meanPrior = rep(0,dim(stats)[2]),
                        sigmaPrior = diag(dim(stats)[2])*1e05,
                        stats = stats,  
                        actor1 = c(0),
                        actor2 = c(0),
                        dyad = attr(reh,"dyad")-1, 
                        omit_dyad = reh$omit_dyad,
                        interevent_time = reh$intereventTime,
                        model = model,
                        ordinal = ordinal,
                        ncores = ncores,
                        thin = thin,
                        L = .GlobalEnv$L_hmc, # temporary
                        epsilon = .GlobalEnv$epsilon_hmc, # temporary,
                        N = reh$N,
                        C = reh$C,
                        D = reh$D)          
            # [enhancement idea]
            #remstimateList$hmc <- hmc # output such that one can run gelman plots and statistics // define methods like remstimate.traceplot() remstimate.
                                        # A. Gelman, Carlin, et al. (2013, 267)
                                        # Stan Development Team (2016 Ch 28.) for how Stan calculates Hat, autocorrelations, and ESS.
                                        # Gelman and Rubin (1992) introduce the R-hat statistic

            hmc$coefficients <- hmc_out$draws[which.min(hmc_out$log_posterior),] # log_posterior is a vector of posterior nloglik
            hmc$post.mean <- colMeans(hmc_out$draws)
            hmc$vcov <- stats::cov(hmc_out$draws)
            hmc$sd <- diag(hmc$vcov)**0.5
            hmc$loglik <- min(hmc_out$log_posterior)
            hmc$df.null <- reh$M
            names(hmc$coefficients) <- names(hmc$post.mean) <- rownames(hmc$vcov) <- colnames(hmc$vcov) <- names(hmc$sd) <- dimnames(stats)[[2]]
            remstimateList <- c(hmc_out,hmc)
        }
        if(model == "actor"){ # Actor Oriented Model
            hmc <- list()
            senderRate <- c(TRUE,FALSE) # meaning that the first BSIR will be run on the "sender rate" model and the second on the "receiver choice"
            which_model <- c("sender_rate","receiver_choice")
            epsilon_loc <- c(0.01,0.01)
            L_loc <- c(100,100)
            for(i in 1:2){
                hmc_i <- list()  

                # setting startng values for the HMC    
                init_i <- matrix(NA,nrow=dim(stats[[which_model[i]]])[2],ncol=nchains)
                noise_factor <- 1
                for(j in 1:nchains){
                    init_i[,j] <- stats::runif((dim(stats[[which_model[i]]])[2]),-5,5) # runif was in (-0.1,0.1)
                    #mvnfast::rmvn(n=1,mu=rep(0,length(dim(stats[[which_model[i]]])[2]),sigma=noise_factor*diag(length(dim(stats[[which_model[i]]])[2])))
                }

          
            
                hmc_out_i <- HMC(pars_init = init_i, 
                            nsim = (burnin+nsim), 
                            nchains = nchains, 
                            burnin = burnin, 
                            meanPrior = rep(0,dim(stats[[which_model[i]]])[2]),
                            sigmaPrior = diag(dim(stats[[which_model[i]]])[2]),
                            stats = stats[[which_model[i]]],  
                            actor1 = reh$edgelist$actor1-1,
                            actor2 = reh$edgelist$actor2-1,
                            dyad = c(0), 
                            omit_dyad = reh$omit_dyad,
                            interevent_time = reh$intereventTime,
                            model = model,
                            senderRate = senderRate[i],
                            N = reh$N,
                            C = reh$C,
                            D = reh$D,
                            ordinal = ordinal,
                            ncores = ncores,
                            thin = thin,
                            L = L_loc[i],
                            epsilon = epsilon_loc[i])   
                # [enhancement idea]
                #remstimateList$hmc <- hmc # output such that one can run gelman plots and statistics // define methods like remstimate.traceplot() remstimate.
                                            # A. Gelman, Carlin, et al. (2013, 267)
                                            # Stan Development Team (2016 Ch 28.) for how Stan calculates Hat, autocorrelations, and ESS.
                                            # Gelman and Rubin (1992) introduce the R-hat statistic

                hmc_i$coefficients <- hmc_out_i$draws[which.max(hmc_out_i$log_posterior),]
                hmc_i$post.mean <- colMeans(hmc_out_i$draws)
                hmc_i$vcov <- stats::cov(hmc_out_i$draws)
                hmc_i$sd <- diag(hmc_i$vcov)**0.5
                hmc_i$loglik <- max(hmc_out_i$log_posterior)
                hmc_i$df.null <- reh$M
                names(hmc_i$coefficients) <- names(hmc_i$post.mean) <- rownames(hmc_i$vcov) <- colnames(hmc_i$vcov) <- names(hmc_i$sd) <- dimnames(stats[[which_model[i]]])[[2]]
                hmc[[which_model[i]]] <- c(hmc_out_i,hmc_i)
                rm(hmc_out_i,hmc_i)
            }
            remstimateList <-  hmc
        }
    }

    str_out <- NULL
    # defining structure and attributes based on the 'method'
    if(method %in% c("MLE","GDADAMAX")){ # frequentist approach
        # defining structure of class 'remstimate'
        str_out <- structure(remstimateList, class = "remstimate")
        # defining attributes
        attr(str_out,"formula") <- model_formula
        attr(str_out,"model") <- model
        attr(str_out,"ordinal") <- ordinal
        attr(str_out, "method") <- method
        attr(str_out, "approach") <- "Frequentist"
        attr(str_out, "statistics") <- variable_names
        if(method == "GDADAMAX"){                 
            attr(str_out, "epochs") <- epochs
            attr(str_out, "epsilon") <- epsilon
            attr(str_out, "iterations") <- ifelse(model=="tie",remstimateList$iterations,c("sender rate"=remstimateList$sender_rate$iterations, "receiver_choice" = remstimateList$receiver_choice$iterations))
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
        if(attr(x,"model") == "tie"){
            cat("\nCoefficients:\n\n")
            print(x$coefficients)
            cat("\nNull deviance:",x$null.deviance,"\nResidual deviance:",x$residual.deviance,"\n")
            cat("AIC:",x$AIC,"AICC:",x$AICC,"BIC:",x$BIC,"\n\n")
        }
        else if(attr(x,"model") == "actor"){
            # printing sender model 
            cat("\nCoefficients rate model:\n\n")
            print(x$sender_rate$coefficients)
            cat("\nNull deviance:",x$sender_rate$null.deviance,"\nResidual deviance:",x$sender_rate$residual.deviance,"\n")
            cat("AIC:",x$sender_rate$AIC,"AICC:",x$sender_rate$AICC,"BIC:",x$sender_rate$BIC,"\n")
            #
            cat(paste0(rep("-", getOption("width")),collapse = ""))
            #
            # printing receiver model
            cat("\n\nCoefficients choice model:\n\n")
            print(x$receiver_choice$coefficients)
            cat("\nNull deviance:",x$receiver_choice$null.deviance,"\nResidual deviance:",x$receiver_choice$residual.deviance,"\n")
            cat("AIC:",x$receiver_choice$AIC,"AICC:",x$receiver_choice$AICC,"BIC:",x$receiver_choice$BIC,"\n\n")
            #     
        }

    }
    else{ # Bayesian
        if(attr(x,"model") == "tie"){
            cat("\nPosterior Modes:\n\n")
            print(x$coefficients)
        }
        else if(attr(x,"model") == "actor"){
            # printing sender model 
            cat("\nPosterior Modes rate model:\n\n")
            print(x$sender_rate$coefficients)
            cat("\n")
            # write some more info here
            #
            cat(paste0(rep("-", getOption("width")),collapse = ""))
            #
            # printing receiver model
            cat("\n\nPosterior Modes choice model:\n\n")
            print(x$receiver_choice$coefficients)
            # write some more info here
            #
        }
    }
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
        if(attr(object,"model") == "tie"){
            coefsTab <- cbind(object$coefficients,
            object$se,
            object$coefficients/object$se,
            2*(1-stats::pnorm(abs(object$coefficients/object$se))),
            (stats::dnorm(0,mean=object$coefficients,sd=object$se) / stats::dnorm(0,mean=0,sd=object$se*sqrt(object$df.null)))/((stats::dnorm(0,mean=object$coefficients,sd=object$se) / stats::dnorm(0,mean=0,sd=object$se*sqrt(object$df.null)))+1)
            )
            colnames(coefsTab) <- c("Estimate","Std. Err", "z value", "Pr(>|z|)", "Pr(=0)")
            rownames(coefsTab) <- attr(object, "statistics")
            summary_out$coefsTab <- coefsTab
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
                        epsilon = attr(object, "epsilon"))                
            summary_out <- do.call(c, list(keep, summary_out))
        }
        else if(attr(object,"model") == "actor"){
            coefsTab <- list()
            # summary sender model
            coefsTab$sender_rate <- cbind(object$sender_rate$coefficients,
            object$sender_rate$se,
            object$sender_rate$coefficients/object$sender_rate$se,
            2*(1-stats::pnorm(abs(object$sender_rate$coefficients/object$sender_rate$se))),
            (stats::dnorm(0,mean=object$sender_rate$coefficients,sd=object$sender_rate$se) / stats::dnorm(0,mean=0,sd=object$sender_rate$se*sqrt(object$sender_rate$df.null)))/((stats::dnorm(0,mean=object$sender_rate$coefficients,sd=object$sender_rate$se) / stats::dnorm(0,mean=0,sd=object$sender_rate$se*sqrt(object$sender_rate$df.null)))+1)
            )
            colnames(coefsTab$sender_rate) <- c("Estimate","Std. Err", "z value", "Pr(>|z|)", "Pr(=0)")
            rownames(coefsTab$sender_rate) <- names(object$sender_rate$coefficients)
            # summary receiver model
            coefsTab$receiver_choice <- cbind(object$receiver_choice$coefficients,
            object$receiver_choice$se,
            object$receiver_choice$coefficients/object$receiver_choice$se,
            2*(1-stats::pnorm(abs(object$receiver_choice$coefficients/object$receiver_choice$se))),
            (stats::dnorm(0,mean=object$receiver_choice$coefficients,sd=object$receiver_choice$se) / stats::dnorm(0,mean=0,sd=object$receiver_choice$se*sqrt(object$receiver_choice$df.null)))/((stats::dnorm(0,mean=object$receiver_choice$coefficients,sd=object$receiver_choice$se) / stats::dnorm(0,mean=0,sd=object$receiver_choice$se*sqrt(object$receiver_choice$df.null)))+1)
            )
            colnames(coefsTab$receiver_choice) <- c("Estimate","Std. Err", "z value", "Pr(>|z|)", "Pr(=0)")
            rownames(coefsTab$receiver_choice) <- names(object$receiver_choice$coefficients)
            #
            summary_out$coefsTab <- coefsTab
            keep <- match(c("formula",
		      "contrasts", "sender_rate", "receiver_choice", "na.action"), names(object), 0L)

            keep <- list(formula = attr(object,"formula"),
                        model = attr(object,"model"),
                        ordinal = attr(object,"ordinal"),
                        method = attr(object, "method"),
                        approach = attr(object, "approach"),
                        sender_rate = list(
                            residual.deviance = object$sender_rate$residual.deviance, 
                            null.deviance = object$sender_rate$null.deviance,
                            model.deviance = object$sender_rate$model.deviance,
                            df.residual = object$sender_rate$df.residual,
                            df.null = object$sender_rate$df.null,
                            df.model = object$receiver_choice$df.model,
                            AIC = object$sender_rate$AIC,
                            AICC = object$sender_rate$AICC,
                            BIC = object$sender_rate$BIC
                        ),
                        receiver_choice = list(
                            residual.deviance = object$receiver_choice$residual.deviance, 
                            null.deviance = object$receiver_choice$null.deviance,
                            model.deviance = object$receiver_choice$model.deviance,
                            df.residual = object$receiver_choice$df.residual,
                            df.null = object$receiver_choice$df.null,
                            df.model = object$receiver_choice$df.model,
                            AIC = object$receiver_choice$AIC,
                            AICC = object$receiver_choice$AICC,
                            BIC = object$receiver_choice$BIC
                        ),
                        epsilon = attr(object, "epsilon"))                      
            summary_out <- do.call(c, list(keep, summary_out))
        }

    }
    if(attr(object, "approach") == "Bayesian"){ # Bayesian
        if(attr(object,"model") == "tie"){
            coefsTab <- cbind(object$coefficients,
            object$sd,
            apply(object$draws,2,stats::quantile,0.025),
            apply(object$draws,2,stats::quantile,0.5),
            apply(object$draws,2,stats::quantile,0.975),
            (stats::dnorm(0,mean=object$coefficients,sd=object$sd) / stats::dnorm(0,mean=0,sd=object$sd*sqrt(object$df.null)))/((stats::dnorm(0,mean=object$coefficients,sd=object$sd) / stats::dnorm(0,mean=0,sd=object$sd*sqrt(object$df.null)))+1)
            )
            colnames(coefsTab) <- c("Post.Mode","Post.SD","Q2.5%","Q50%","Q97.5%","Pr(=0|x)")
            rownames(coefsTab) <- attr(object, "statistics")
            summary_out$coefsTab <- coefsTab

            keep <- list(formula = attr(object,"formula"),
                        model = attr(object,"model"),
                        ordinal = attr(object,"ordinal"),
                        method = attr(object, "method"),
                        approach = attr(object, "approach"),
                        statistics = attr(object, "statistics"),
                        prior = attr(object, "prior"),
                        nsim = attr(object, "nsim"),
                        seed = attr(object, "seed"),
                        loglik = object$loglik
                        ) 
            keep_HMC <- list()
            if(attr(object, "method") == "HMC"){         
                keep_HMC <- list(
                            nchains = attr(object, "nchains"),
                            burnin = attr(object, "burnin"),
                            thin = attr(object, "thin")         
                            )    
            }                             
            summary_out <- do.call(c, list(keep, keep_HMC, summary_out))
        }
        if(attr(object,"model") == "actor"){
            coefsTab <- list()
            # summary sender model
            coefsTab$sender_rate <- cbind(object$sender_rate$coefficients,
            object$sender_rate$sd,
            apply(object$sender_rate$draws,2,stats::quantile,0.025),
            apply(object$sender_rate$draws,2,stats::quantile,0.5),
            apply(object$sender_rate$draws,2,stats::quantile,0.975),
            (stats::dnorm(0,mean=object$sender_rate$coefficients,sd=object$sender_rate$sd) / stats::dnorm(0,mean=0,sd=object$sender_rate$sd*sqrt(object$sender_rate$df.null)))/((stats::dnorm(0,mean=object$sender_rate$coefficients,sd=object$sender_rate$sd) / stats::dnorm(0,mean=0,sd=object$sender_rate$sd*sqrt(object$sender_rate$df.null)))+1)
            )
            colnames(coefsTab$sender_rate) <- c("Post.Mode","Post.SD","Q2.5%","Q50%","Q97.5%","Pr(=0|x)")
            rownames(coefsTab$sender_rate) <- names(object$sender_rate$coefficients)
            # summary receiver model
            coefsTab$receiver_choice <- cbind(object$receiver_choice$coefficients,
            object$receiver_choice$sd,
            apply(object$receiver_choice$draws,2,stats::quantile,0.025),
            apply(object$receiver_choice$draws,2,stats::quantile,0.5),
            apply(object$receiver_choice$draws,2,stats::quantile,0.975),
            (stats::dnorm(0,mean=object$receiver_choice$coefficients,sd=object$receiver_choice$sd) / stats::dnorm(0,mean=0,sd=object$receiver_choice$sd*sqrt(object$receiver_choice$df.null)))/((stats::dnorm(0,mean=object$receiver_choice$coefficients,sd=object$receiver_choice$sd) / stats::dnorm(0,mean=0,sd=object$receiver_choice$sd*sqrt(object$receiver_choice$df.null)))+1)
            )
            colnames(coefsTab$receiver_choice) <- c("Post.Mode","Post.SD","Q2.5%","Q50%","Q97.5%","Pr(=0|x)")
            rownames(coefsTab$receiver_choice) <- names(object$receiver_choice$coefficients)
            #
            summary_out$coefsTab <- coefsTab
            keep <- match(c("formula",
		      "contrasts", "sender_rate", "receiver_choice", "na.action"), names(object), 0L)

            keep <- list(formula = attr(object,"formula"),
                        model = attr(object,"model"),
                        ordinal = attr(object,"ordinal"),
                        method = attr(object, "method"),
                        approach = attr(object, "approach"),
                        sender_rate = list(
                            #residual.deviance = object$sender_rate$residual.deviance, 
                            #null.deviance = object$sender_rate$null.deviance,
                            #model.deviance = object$sender_rate$model.deviance,
                            #df.residual = object$sender_rate$df.residual,
                            df.null = object$sender_rate$df.null, #,
                            #df.model = object$sender_rate$df.model,
                            loglik = object$sender_rate$loglik
                        ),
                        receiver_choice = list(
                            #residual.deviance = object$receiver_choice$residual.deviance, 
                            #null.deviance = object$receiver_choice$null.deviance,
                            #model.deviance = object$receiver_choice$model.deviance,
                            #df.residual = object$receiver_choice$df.residual,
                            df.null = object$receiver_choice$df.null, #,
                            #df.model = object$receiver_choice$df.model,
                            loglik = object$receiver_choice$loglik
                        ),
                        epsilon = attr(object, "epsilon")) # add other attributes from the objects?                     
            summary_out <- do.call(c, list(keep, summary_out))
        }
    }

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
        stats::printCoefmat(x$coefsTab, P.values = TRUE, signif.stars = FALSE, ...) 
        if(x$approach == "Frequentist"){
            cat("Null deviance:", x$null.deviance, "on", x$df.null, "degrees of freedom\n")
            cat("Residual deviance:", x$residual.deviance, "on", x$df.residual, "degrees of freedom\n")
            cat("Chi-square:", x$model.deviance, "on", x$df.model, 
                "degrees of freedom, asymptotic p-value", 1 - stats::pchisq(x$model.deviance, 
                    x$df.model), "\n")
            cat("AIC:", x$AIC, "AICC:", x$AICC, "BIC:", x$BIC, "\n")
        }
        if(x$approach == "Bayesian"){
        cat("Log posterior:",x$loglik,"\n")
        # cat("Prior parameters:",paste(names(x$prior.param),unlist(x$prior.param),sep="="),"\n")
        }
    }
    else if(x$model == "actor"){
        second_line <- paste("(",x$method," with ",sep="")
        if(x$ordinal){
            second_line <- paste(second_line,"ordinal likelihood):",sep="")
        }
        else{
            second_line <- paste(second_line,"interval likelihood):\n\n",sep="")
        }
        # sender rate summary
        cat("Call rate model:\n\n\t",deparse(x$formula$rate_model_formula),"\n\n",sep="")
        if(x$approach == "Frequentist"){
            cat("\nCoefficients rate model",second_line)
        }
        else{ # Bayesian
            cat("\nPosterior Modes rate model",second_line)
        }
        stats::printCoefmat(x$coefsTab$sender_rate, P.values = ifelse(x$approach == "Frequentist",TRUE,FALSE), signif.stars = FALSE)
        if(x$approach == "Frequentist"){
            cat("Null deviance:", x$sender_rate$null.deviance, "on", x$sender_rate$df.null, "degrees of freedom\n")
            cat("Residual deviance:", x$sender_rate$residual.deviance, "on", x$sender_rate$df.residual, "degrees of freedom\n")
            cat("Chi-square:", x$sender_rate$model.deviance, "on", x$sender_rate$df.model, 
                "degrees of freedom, asymptotic p-value", 1 - stats::pchisq(x$sender_rate$model.deviance, 
                    x$sender_rate$df.model), "\n")
            cat("AIC:", x$sender_rate$AIC, "AICC:", x$sender_rate$AICC, "BIC:", x$sender_rate$BIC, "\n")
        }
        if(x$approach == "Bayesian"){
        cat("Log posterior:",x$sender_rate$loglik,"\n")
        # cat("Prior parameters:",paste(names(x$sender_rate$prior.param),unlist(x$sender_rate$prior.param),sep="="),"\n")
        }
        #
        cat(paste0(rep("-", getOption("width")),collapse = ""),"\n\n")
        #
        # receiver choice summary 
        cat("Call choice model:\n\n\t",deparse(x$formula$choice_model_formula),"\n\n",sep="") 
                if(x$approach == "Frequentist"){
            cat("\nCoefficients choice model",second_line)
        }
        else{ # Bayesian
            cat("\nPosterior Modes choice model",second_line)
        }
        stats::printCoefmat(x$coefsTab$receiver_choice, P.values = ifelse(x$approach == "Frequentist",TRUE,FALSE), signif.stars = FALSE)
        if(x$approach == "Frequentist"){
            cat("Null deviance:", x$receiver_choice$null.deviance, "on", x$receiver_choice$df.null, "degrees of freedom\n")
            cat("Residual deviance:", x$receiver_choice$residual.deviance, "on", x$receiver_choice$df.residual, "degrees of freedom\n")
            cat("Chi-square:", x$receiver_choice$model.deviance, "on", x$receiver_choice$df.model, 
                "degrees of freedom, asymptotic p-value", 1 - stats::pchisq(x$receiver_choice$model.deviance, 
                    x$receiver_choice$df.model), "\n")
            cat("AIC:", x$receiver_choice$AIC, "AICC:", x$receiver_choice$AICC, "BIC:", x$receiver_choice$BIC, "\n")
        }
        if(x$approach == "Bayesian"){
        cat("Log posterior:",x$receiver_choice$loglik,"\n")
        # cat("Prior parameters:",paste(names(x$receiver_choice$prior.param),unlist(x$receiver_choice$prior.param),sep="="),"\n")
        }
    }
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
 return(paste('this function at the moment does nothing'))
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
 return(paste('this function at the moment does nothing'))
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
aic.remstimate <- function(object,...) {
    if(attr(object, "approach") == "Frequentist"){
        return(object$AIC)
    }
    else{
        stop("'approach' must be 'Frequentist'")
    }
}


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
aicc.remstimate <- function(object,...) {
    if(attr(object, "approach") == "Frequentist"){
        return(object$AICC)
    }
    else{
        stop("'approach' must be 'Frequentist'")
    }
}


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
bic.remstimate <- function(object,...) {
    if(attr(object, "approach") == "Frequentist"){
        return(object$BIC)
    }
    else{
        stop("'approach' must be 'Frequentist'")
    }
}


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
waic.remstimate <- function(object,...) {
    if(attr(object, "approach") == "Frequentist"){
        if(!is.null(object$WAIC)){
            return(object$WAIC)
        }
        else{
            return(paste("'WAIC' value not found."))
        }
    }
        else{stop("'approach' must be 'Frequentist'")}   
}


#######################################################################################
#######################################################################################
