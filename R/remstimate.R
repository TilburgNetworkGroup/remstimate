#' remstimate  
#'
#' A function for the optimization of tie-oriented and actor-oriented likelihood. There are four optimization algorithms: two Frequentists, Maximum Likelihood Estimation (\code{MLE}) and Adaptive Gradient Descent (\code{GDADAMAX}), and two Bayesian, Bayesian Sampling Importance Resampling (\code{BSIR}) and Hamiltonian Monte Carlo (\code{HMC}).
#'
#' @param reh a \code{remify} object of the processed relational event history. Output object of the function \code{remify::remify()}.
#' @param stats a \code{remstats} object: when `attr(reh,"model")` is `"tie"`, \code{stats} is an array of statistics with dimensions \code{[M x D x P]}: where \code{M} is the number of events, \code{D} is the number of possible dyads (full riskset), \code{P} is the number of statistics; if `attr(reh,"model")` is `"actor"`, \code{stats} is a list that can contain up to two arrays named \code{"sender_stats"} and \code{"receiver_stats"} with dimensions \code{[M x N x P]}, where \code{N} are the actors (senders in the array \code{"sender_stats"}, receivers in the array \code{"receiver_stats"}). Furthermore, it is possible to only estimate the sender rate model or only the receiver choice model, by using the correct naming of the arrays.
#' @param method the optimization method to estimate model parameters. Methods available are: Maximum Likelihood Estimation (\code{"MLE"}), Adaptive Gradient Descent (\code{"GDADAMAX"}), Bayesian Sampling Importance Resampling (\code{"BSIR"}), Hamiltonian Monte Carlo (\code{"HMC"}). (default method is \code{"MLE"}).
#' @param ncores [\emph{optional}] number of threads for the parallelization. (default value is \code{1}, which means no parallelization)
#' @param prior [\emph{optional}] prior distribution when \code{method} is \code{"BSIR"}. Default value is \code{NULL}, which means that no prior is assumed. For the tie-oriented modeling, the argument \code{prior} is the name of the function in the format \code{name_package::name_density_function}. The parameters of the prior distribution can be supplied as inputs to the remstimate function (e.g., \code{remstimate::remstimate(reh=reh,stats=stats,method="BSIR",ncores=5,prior=mvnfast::dmvn,mu=rep(0,3),sigma=diag(3)*2,log=TRUE)} ). For actor-oriented modeling the argument \code{prior} is a named list of two objects \code{"sender_model"}, which calls the prior function for the sender rate model, and, \code{"receiver_model"}, which calls the prior function for the receiver choice model. For the specification of the prior parameters, the user must define an optional argument called \code{prior_args}, which is also a named list (with names \code{"sender_model"} and \code{"receiver_model"}): each list is a list of objects named after the prior arguments and with value of the prior argument (e.g., \code{prior_args$sender_model = list(mu = rep(1.5,3), sigma = diag(3)*0.5, log = TRUE)}). Finally, both in tie-oriented and actor-oriented modeling prior functions must have an argument that returns the value of the density on a logarithmic scale (i.e., \code{log=TRUE}). \code{log=TRUE} is already set up internally by \code{remstimate()}.
#' @param nsim  [\emph{optional}] when \code{method} is \code{"HMC"}, \code{nsim} is the number of simulations (iterations) in each chain, when \code{method} is \code{"BSIR"}, then \code{nsim} is the number of samples from the proposal distribution. Default value is \code{1000}.
#' @param nchains [\emph{optional}] number of chains to generate in the case of \code{method = "HMC"}. Default value is \code{1}.
#' @param burnin [\emph{optional}] number of initial iterations to be added as burnin for \code{method = "HMC"}. Default value is \code{500}.
#' @param thin [\emph{optional}] number of steps to skip in the posterior draws for \code{method = "HMC"}. Default value is \code{10}. If \code{nsim<100}, thin is set to \code{1}
#' @param init [\emph{optional}] vector of initial values if tie-oriented model, or a named list of two vectors ('sender_model' and 'receiver_model') if both models of the actor-oriented framework are specified. \code{init} can also be a list of only one vector (named 'sender_model' or 'receiver_model'), if the interest is to estimate one specific model of the actor-oriented framework. \code{init} is used for the methods \code{"GDADAMAX"} and \code{"HMC"}. If \code{init} is \code{NULL}, then it will be assigned internally.
#' @param epochs [\emph{optional}] It is the number of iteration used in the method \code{"GDADAMAX"}. Default value is \code{1000}.
#' @param L [\emph{optional}] number of leap-frog steps to use in the method \code{"HMC"}. Default value is \code{50}.
#' @param epsilon [\emph{optional}] It is a parameter used in two methods: if \code{method} is \code{"GDADAMAX"}, it represents the inter-iteration difference of the loss function and it is used as stop-rule within the algorithm (default value is \code{0.001}), if \code{method} is \code{"HMC"} (default value is \code{0.002}), it is a parameter used in the leap-frog algorithm and it is proportional to the step size.
#' @param seed [\emph{optional}] seed value for reproducibility. If \code{NULL}, seed will be assigned by the machine and saved in the output object.
#' @param silent [\emph{optional}-not-yet-implemented] a \code{TRUE/FALSE} value. If \code{FALSE}, progress of optimization status will be printed out. 
#' @param ... additional parameters. They can be parameters of other functions defined as input in some of the arguments above. (e.g., arguments of the \code{prior} distribution)
#'
#' @return  'remstimate' S3 object
#' 
#' @examples 
#' 
#' # ------------------------------------ #
#' #       tie-oriented model: "MLE"      #
#' # ------------------------------------ #
#' 
#' # loading data
#' data(tie_reh)
#'   
#' # specifying linear predictor
#' tie_model <- ~ 1 + 
#'                remstats::indegreeSender()+
#'                remstats::inertia()+
#'                remstats::reciprocity() 
#' 
#' # calculating statistics
#' tie_reh_stats <- remstats::remstats(reh = tie_reh, 
#'                                     tie_effects = tie_model)
#' 
#' # running estimation
#' tie_mle <- remstimate::remstimate(reh = tie_reh,
#'                                   stats = tie_reh_stats,
#'                                   method = "MLE",
#'                                   ncores = 1)
#' # summary
#' summary(tie_mle)
#' 
#' # ------------------------------------ #
#' #      actor-oriented model: "MLE"     #
#' # ------------------------------------ #
#' 
#' # loading data
#' data(ao_reh)
#'   
#' # specifying linear predictor (for sender rate and receiver choice model)
#' rate_model <- ~ 1 + remstats::indegreeSender()
#' choice_model <- ~ remstats::inertia() + remstats::reciprocity()
#' 
#' # calculating statistics
#' ao_reh_stats <- remstats::remstats(reh = ao_reh, 
#'                                    sender_effects = rate_model, 
#'                                    receiver_effects = choice_model)
#' 
#' # running estimation
#' ao_mle <- remstimate::remstimate(reh = ao_reh,
#'                                  stats = ao_reh_stats,
#'                                  method = "MLE",
#'                                  ncores = 1)
#' # summary
#' summary(ao_mle)
#' 
#' # ------------------------------------ #
#' #   for more examples check vignettes  #
#' # ------------------------------------ #
#' 
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
                       L = 50L,
                       epsilon = ifelse(method=="GDADAMAX",0.001,0.002),
                       seed = NULL,
                       silent = TRUE,
                       ...){

    # ... processing input arguments

    # ... additional arguments supplied via the three dots ellipsis
    additional_input_args <- list(...) # we include 'ncores' in the case prior distribution routines allow for parallelization (e.g., the 'mvnfast' package)

    if(!any(names(additional_input_args) %in% "log")){
        additional_input_args$log <- TRUE
    }
    if(!any(names(additional_input_args) %in% "ncores")){
        additional_input_args$ncores <- ncores
    }

    # ... remify object ('reh' input argument)
    if(!inherits(reh,"remify")){
        stop("'reh' must be a 'remify' object (see ?remify::remify).")
    }

    # ... model
    model <- attr(reh, "model") # attribute from reh object
    if(is.null(model) | !(model %in% c("actor","tie"))){
        stop("attribute 'model' of input 'reh' must be either 'actor' or 'tie'")
    }

    # ... active riskset? then overwrite two objects (this prevents from coding many ifelse() to switch between active-oriented riskset objects and full-oriented riskset objects)
    if((attr(reh,"riskset") == "active") & (model == "tie")){
        reh$D <- reh$activeD
        attr(reh,"dyad") <- attr(reh,"dyadIDactive")
    }

    # ... method
    if(!(method %in% c("MLE","GDADAMAX","BSIR","HMC"))){stop("The `method` specified is not available or it is mistyped.")}
    else{
        method  <- match.arg(arg = method, choices = c("MLE","GDADAMAX","BSIR","HMC"), several.ok = FALSE) # default is "MLE"
    }

    # ... type of likelihood
    ordinal <- attr(reh,"ordinal")

    # ... directed / undirected network
    if((model == "actor") & !attr(reh,"directed")){stop("actor-oriented modeling can't operate on undirected networks")}

    # ... ncores
    if(is.null(ncores)) ncores <- 1L
    else if(((parallel::detectCores() > 2L) & (ncores > floor(parallel::detectCores()-2L))) | ((parallel::detectCores() == 2L) & (ncores > 1L))){
            stop("'ncores' is recommended to be set at most to: floor(parallel::detectCores()-2L)")
    }

    # ... stats 
    model_formula <- variable_names <- where_is_baseline <- NULL
    if(model == "tie")
    {
        if(all(inherits(stats,c("remstats","tomstats"),TRUE))){
            variable_names <- dimnames(stats)[[3]] #as.vector(sapply(dimnames(stats)[[3]],function(x) sub(pattern = ".x.", replacement = ":", x = x)))
            model_formula <- stats::as.formula(paste("~ ",paste(variable_names,collapse=" + ")))
            stats <- aperm(stats, perm = c(2,3,1)) # stats reshaped in [D*U*M]
            # is there a baseline term?
            if(any(variable_names %in% c("baseline"))){
                where_is_baseline <- which(variable_names == "baseline")
            }
        }
        else if(all(inherits(stats,c("remstats","aomstats"),TRUE))){
            stop("'remstats' object supplied cannot work for tie-oriented modeling")
        }
        else{
            stop("tie-oriented modeling: 'stats' must be a 'tomstats' 'remstats' object")
        }
    }
    if(model == "actor") 
    {
        model_formula <- list() # becomes a list
        if(all(inherits(stats,c("remstats","aomstats"),TRUE))){
            variables_rate <- variables_choice <- NULL
            if(!is.null(stats$sender_stats)){ # sender model is specified
                variables_rate <- dimnames(stats$sender_stats)[[3]] #as.vector(sapply(dimnames(stats$sender_stats)[[3]],function(x) sub(pattern = ".x.", replacement = ":", x = x)))
                model_formula[["rate_model_formula"]] <- stats::as.formula(paste("~ ",paste(variables_rate,collapse=" + ")))
                stats$sender_stats <- aperm(stats$sender_stats, perm = c(2,3,1)) # stats reshaped in [N*U*M]
                # is there a baseline term?
                where_is_baseline <- which(variables_rate == "baseline")
            }
            if(!is.null(stats$receiver_stats)){ # receiver model is specified
                variables_choice <- dimnames(stats$receiver_stats)[[3]] #as.vector(sapply(dimnames(stats$receiver_stats)[[3]],function(x) sub(pattern = ".x.", replacement = ":", x = x)))
                model_formula[["choice_model_formula"]] <- stats::as.formula(paste("~ ",paste(variables_choice,collapse=" + ")))
                stats$receiver_stats <- aperm(stats$receiver_stats, perm = c(2,3,1)) # stats reshaped in [N*U*M]
            }
            # vector of variable names and list of model formulas
            variable_names <- c(variables_rate,variables_choice)
        }
        else if(all(inherits(stats,c("remstats","tomstats"),TRUE))){
            stop("'remstats' object supplied cannot work for tie-oriented modeling")
        }
        else{
            stop("actor-oriented modeling: 'stats' must be either a 'aomstats' 'remstats' object or a list of two arrays named 'sender_stats' and 'receiver_stats'")
        }
    }

    # ... GDADAMAX: optional arguments
    if(method == "GDADAMAX"){
        # ... epsilon
        if(is.null(epsilon)){
            epsilon <- 0.001
        }
        else if(length(epsilon)>1){
            epsilon <- 0.001
            warning("'epsilon' is set to its default value: 0.001") 
        }
        else if((epsilon <= 0) | !is.numeric(epsilon)){
            epsilon <- 0.001
            warning("'epsilon' is set to its default value: 0.001")
        }
        # ... epochs
        if(is.null(epochs)){
            epochs <- 1e03L
        }
        else if(is.numeric(epochs) | is.integer(epochs)){
            epochs <- as.integer(epochs) # this converts, for instance, 3.2 to 3, 3.7 to 3
            if(epochs<0){
                stop("'epoch' must be a positive number")
            }
        }
        else{
            stop("'epoch' must be a positive number")
        }
    }

    # ... BSIR: optional arguments
    if(method == "BSIR"){
        # ... seed
        if(is.null(seed)){
            seed <- .Random.seed
        }
         else if(length(seed) > 1){
            seed <- seed[1]
            warning("`seed` length is greater than 1. Considering only the first element")
        }

        # ... nsim
        if(is.null(nsim)){
            nsim <- 1e03L # setting default value of 'nsim'
        }

        # ... prior distribution
        list_args <- args_prior <- args_prior_sender <- args_prior_receiver <- NULL
        if(!is.null(prior)){
            if(model == "tie"){                      
                args_prior <- match.arg(arg=names(additional_input_args),choices = methods::formalArgs(prior),several.ok = TRUE)
                list_args <- lapply(args_prior, function(x) additional_input_args[[x]])
                names(list_args) <- args_prior
            }
            else if(model == "actor"){
                list_args <- list()

                # sender model
                args_prior_sender <- match.arg(arg=names(additional_input_args$prior_args$sender_model),choices = methods::formalArgs(prior[["sender_model"]]),several.ok = TRUE)
                list_args[["sender_model"]] <- lapply(args_prior_sender, function(x) additional_input_args$prior_args$sender_model[[x]])
                names(list_args[["sender_model"]]) <- args_prior_sender

                # receiver model
                args_prior_receiver <- match.arg(arg=names(additional_input_args$prior_args$receiver_model),choices = methods::formalArgs(prior[["receiver_model"]]),several.ok = TRUE)
                list_args[["receiver_model"]] <- lapply(args_prior_receiver, function(x) additional_input_args$prior_args$receiver_model[[x]])
                names(list_args[["receiver_model"]]) <- args_prior_receiver
            }
        }
    }        
    
   

    # ... HMC: optional arguments
    if(method == "HMC"){
        # ... seed
        if(is.null(seed)){
                seed <- .Random.seed
        }
        else if(length(seed) > 1){
            seed <- seed[1]
            warning("`seed` length is greater than 1. Considering only the first element")
        }

        # ... nchains
        if(is.null(nchains)){
            nchains <- 1L # setting default value of 'nchains'
        }

        # ... nsim 
        if(is.null(nsim)){
            nsim <- 1e03L # setting default value of 'nsim'
        }

        # ... burnin
        if(is.null(burnin)){
            burnin <- 500L # setting default value of 'burnin'
        }

        # ... thin
        if(is.null(thin)){
            if(nsim >= 100){
                thin <- 10L
                warning("'thin' parameter undefined. Using thin = 10")
            }
            else{
                thin <- 1L # no thinning is applied
                warning("'nsim' is less than 100. No thinning is applied")
            }
        }
        else if(is.numeric(thin)){
            if(thin <= 0){
                stop("`thin` value must be positive. Set thin = 1 if no thinning of chains is required")
            }
            if(thin == 1){
                warning("`thin` value is 1. No thinning applied to chains")
            }
        }

        # ... L (number of leap-frog steps)
        if(is.null(L)){
            L <- 50L # setting default value of 'L'
        }
        else if(length(L)>1){
            L <- 50L
            warning("input 'L' must be a positive (integer) number . 'L' is set to its default value: 50")            
        }
        else if((L<=1) | !is.numeric(L)){
            L <- 50L
            warning("input 'L' must be a positive (integer) number . 'L' is set to its default value: 50")
        }
        # ... epsilon (sie of time step in the leap-frog algorithm)
        if(is.null(epsilon)){
            epsilon <- 0.1/L # 0.1 is the full time, reached by steps of 0.1/L
        }
        else if(length(epsilon)>1){
            epsilon <- 0.1/L
            warning("input 'epsilon' must be a positive number. 'epsilon' is set to its default value: 0.1/L")            
        }
        else if((epsilon <= 0) | !is.numeric(epsilon)){
            epsilon <- 0.1/L
            warning("input 'epsilon' must be a positive number. 'epsilon' is set to its default value: 0.1/L")
        }
        # 0.1 is a [temporary parameter] and it will change in the future releases of 'remstimate'
    }

    # ... init (initial values) for GDADAMAX and HMC
    if(method %in% c("GDADAMAX","HMC")){
        if(!is.null(init)){
            if(model == "actor"){
                if(!is.null(init)){ # if 'init' is not NULL, we check it. otherwise we will define it internally
                    if(!is.list(init)){
                            stop("'init' must be a list of two vectors named 'sender_model' and 'receiver_model', each one with the starting values of the effects of the statistics according to the two different models (rate model and choice model)")
                        }
                    else if(!all(c("sender_model","receiver_model") %in% names(init))){
                            stop("'init' must be a list of two vectors named 'sender_model' and 'receiver_model', each one with the starting values of the effects of the statistics according to the two different models (rate model and choice model)")
                        }
                    else if(length(init$sender_model)!=dim(stats$sender_stats)[2] | length(init$receiver_model)!=dim(stats$receiver_stats)[2]){
                        stop("each element of list 'init' must be equal to number of statistics according to the arrays supplied in 'stats'")
                    }
                }
            }
            else if(model == "tie"){
                if(!is.null(init)){ # if 'init' is not NULL, we check it. otherwise we will define it internally
                    if(!is.vector(init)){
                        stop("'init' must be a vector with the starting values of the effects of the statistics")
                    }
                    else if(length(init)!=dim(stats)[2]){
                        stop("length of vector 'init' must be equal to the number of statistics in 'stats'")   
                    }
                }
            }
        }
        else{
            if(model == "actor"){
                init <- list()
            }
        }
    }

    # ... creating  an empty lists
    remstimateList <- list()
        
    # ... estimating with Maximum Likelihood (MLE) or Gradient Descent ADAMAX (GDADAMAX)
    if(method %in% c("MLE","GDADAMAX")){
        if(model == "tie"){ # Tie-oriented Model  (Relational Event Model)
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
                    beta_0 <- log(reh$M) - log(sum(reh$intereventTime)*reh$D) # MLE of the intercept under only intercept model
                    init <- stats::runif(dim(stats)[2],-0.1,0.1) # previously rep(0,dim(stats)[2])
                    if(any(dimnames(stats)[[2]] == "baseline")){
                        init[which(dimnames(stats)[[2]] == "baseline")] <- init[which(dimnames(stats)[[2]] == "baseline")] + beta_0
                    }
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
                optimum_obj$hessian <- remDerivativesStandard(pars = optimum_obj$argument,
                                                                    stats = stats,
                                                                    dyad = attr(reh,"dyad")-1,
                                                                    omit_dyad = reh$omit_dyad,
                                                                    interevent_time = reh$intereventTime,
                                                                    ordinal = ordinal,
                                                                    ncores = ncores)
                optimum_obj$hessian <- optimum_obj$hessian$hessian 
                                                                
            }            
            # coefficients                         
            remstimateList$coefficients <- as.vector(optimum_obj$argument)
            names(remstimateList$coefficients) <- variable_names
            # loglik
            remstimateList$loglik <- -optimum_obj$value # loglikelihood value at MLE values (we take the '-value' because we minimized the '-loglik')       
            # gradient 
            remstimateList$gradient <- -optimum_obj$gradient
            # hessian
            remstimateList$hessian <- optimum_obj$hessian # hessian matrix relative
            # vcov
            remstimateList$vcov <- qr.solve(optimum_obj$hessian) # matrix of variances and covariances
            # standard errors
            remstimateList$se <- diag(remstimateList$vcov)**0.5 # standard errors
            # re-naming
            names(remstimateList$coefficients) <- names(remstimateList$se) <- rownames(remstimateList$vcov) <- colnames(remstimateList$vcov) <- dimnames(stats)[[2]]
            # residuals
            select_vars <- if(is.null(where_is_baseline)) 1:length(variable_names) else c(1:length(variable_names))[-where_is_baseline]
            baseline_value <- if(is.null(where_is_baseline)) 0 else as.vector(optimum_obj$argument)[where_is_baseline]
            stats_loc <- if(length(select_vars)==1) array(stats[,select_vars,],dim=c(dim(stats)[1],1,dim(stats)[3])) else stats[,select_vars,]
            remstimateList$diagnostics <- list()
            remstimateList$diagnostics$residuals <- computeDiagnostics(pars = as.vector(optimum_obj$argument)[select_vars],
                                                    stats = stats_loc,
                                                    actor1 = c(0),
                                                    actor2 = c(0),
                                                    dyad = attr(reh,"dyad")-1,
                                                    omit_dyad = reh$omit_dyad,
                                                    model = attr(reh,"model"),
                                                    ncores = ncores,
                                                    baseline = baseline_value,
                                                    N = reh$N)
            colnames(remstimateList$diagnostics$residuals$standardized_residuals) <- variable_names[select_vars]     
            colnames(remstimateList$diagnostics$residuals$smoothing_weights) <- colnames(remstimateList$diagnostics$residuals$standardized_residuals)
            # lambdas (event rates)
            remstimateList$diagnostics$rates <- remstimateList$diagnostics$residuals$rates   
            remstimateList$diagnostics$residuals$rates <- NULL                             
            # residual.deviance
            remstimateList$residual.deviance <- -2*remstimateList$loglik        
            # null.deviance
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
                                        ncores = ncores)$value) # [[NOTE: the computation of the null deviance can be simplified. However, if omit_dyad is specified, the simplification should also account for the dynamic riskset]]
            # model.deviance                            
            remstimateList$model.deviance <- remstimateList$null.deviance -  remstimateList$residual.deviance  
            # df.null                   
            remstimateList$df.null <- reh$M
            # df.model
            remstimateList$df.model <- dim(stats)[2]
            # df.residual
            remstimateList$df.residual <- remstimateList$df.null - remstimateList$df.model
            # AIC
            remstimateList$AIC <- 2*length(remstimateList$coefficients) - 2*remstimateList$loglik # AIC
            # AICC
            remstimateList$AICC <- remstimateList$AIC + 2*length(remstimateList$coefficients)*(length(remstimateList$coefficients)+1)/(reh$M-length(remstimateList$coefficients)-1)
            # BIC
            remstimateList$BIC <- length(remstimateList$coefficients)*log(reh$M) - 2*remstimateList$loglik # BIC
            # converged
            remstimateList$converged <- optimum_obj$converged
            # iterations
            remstimateList$iterations <- optimum_obj$iterations                                    
        }
        if(model == "actor" ){ # Actor-oriented Model 
            optimum_model <- list()
            senderRate <- c(TRUE,FALSE) # meaning that the first BSIR will be run on the "sender rate" model and the second on the "receiver choice"
            which_model <- c("sender_model","receiver_model")
            which_stats <- c("sender_stats","receiver_stats")
            for(i in 1:2){
                if(!is.null(stats[[which_stats[i]]])){  # run estimation only if the model is specified
                    if(method == "MLE"){ # Maximum Likelihood Estimates
                        # code here
                        optimum_model[[which_model[i]]] <- trust::trust(objfun = remDerivatives, 
                                    parinit = rep(0,dim(stats[[which_stats[i]]])[2]), 
                                    rinit = 1, 
                                    rmax = 100, 
                                    stats = stats[[which_stats[i]]], 
                                    actor1 = reh$edgelist$actor1-1,
                                    actor2 = reh$edgelist$actor2-1,
                                    dyad = c(0),
                                    omit_dyad = reh$omit_dyad, 
                                    interevent_time = reh$intereventTime,
                                    model = model,
                                    senderRate = senderRate[i],
                                    N = reh$N,
                                    C = ifelse(is.null(reh$C),1,reh$C),
                                    D = reh$D)
                    }
                    else if(method == "GDADAMAX"){ # Gradient Descent
                        if(is.null(init[[which_model[i]]])){
                            init[[which_model[i]]] <- stats::runif(dim(stats[[which_stats[i]]])[2],-0.1,0.1) # rep(0.0,dim(stats[[which_stats[i]]])[2]) 
                        }
                        optimum_model[[which_model[i]]] <- GDADAMAX(pars = init[[which_model[i]]], 
                                                stats = stats[[which_stats[i]]],  
                                                actor1 = reh$edgelist$actor1-1,
                                                actor2 = reh$edgelist$actor2-1,
                                                dyad = c(0),
                                                omit_dyad = reh$omit_dyad,
                                                interevent_time = reh$intereventTime,
                                                model = model,
                                                ordinal = ordinal,
                                                senderRate = senderRate[i],
                                                N = reh$N,
                                                C = ifelse(is.null(reh$C),1,reh$C),
                                                D = reh$D,
                                                ncores = ncores,
                                                epochs = epochs,
                                                epsilon = epsilon)
                        optimum_model[[which_model[i]]]$hessian <- remDerivatives(pars = optimum_model[[which_model[i]]]$argument,
                                                                        stats = stats[[which_stats[i]]],  
                                                                        actor1 = reh$edgelist$actor1-1,
                                                                        actor2 = reh$edgelist$actor2-1,
                                                                        dyad = c(0),
                                                                        omit_dyad = reh$omit_dyad,
                                                                        interevent_time = reh$intereventTime,
                                                                        model = model,
                                                                        ordinal = ordinal,
                                                                        senderRate = senderRate[i],
                                                                        N = reh$N,
                                                                        C = ifelse(is.null(reh$C),1,reh$C),
                                                                        D = reh$D,
                                                                        ncores = ncores,
                                                                        hessian = TRUE) 
                        optimum_model[[which_model[i]]]$hessian <- optimum_model[[which_model[i]]]$hessian$hessian   
                    }

                    # output for the receiver model
                    remstimateList[[which_model[i]]] <- list() 
                    # coefficients
                    remstimateList[[which_model[i]]]$coefficients <- as.vector(optimum_model[[which_model[i]]]$argument)
                    # loglik
                    remstimateList[[which_model[i]]]$loglik <- -optimum_model[[which_model[i]]]$value # log(L_choice)  
                    # gradient     
                    remstimateList[[which_model[i]]]$gradient <- optimum_model[[which_model[i]]]$gradient
                    # hessian
                    remstimateList[[which_model[i]]]$hessian <- optimum_model[[which_model[i]]]$hessian
                    # vcov
                    remstimateList[[which_model[i]]]$vcov <- qr.solve(remstimateList[[which_model[i]]]$hessian) # matrix of variances and covariances for the receiver choice model
                    # standard errors
                    remstimateList[[which_model[i]]]$se <- diag(remstimateList[[which_model[i]]]$vcov)**0.5 # standard errors         
                    # re-naming
                    names(remstimateList[[which_model[i]]]$coefficients) <- names(remstimateList[[which_model[i]]]$se) <- rownames(remstimateList[[which_model[i]]]$vcov) <- colnames(remstimateList[[which_model[i]]]$vcov) <- colnames(remstimateList[[which_model[i]]]$hessian) <- rownames(remstimateList[[which_model[i]]]$hessian) <- if(i == 1) variables_rate else variables_choice
                    # residuals
                    baseline_value <- 0
                    select_vars <- c(1:dim(stats[[which_stats[i]]])[2])
                    if(senderRate[i]){ # only for sender model
                        baseline_value <- if(is.null(where_is_baseline)) 0 else as.vector(optimum_model[[which_model[i]]]$argument)[where_is_baseline]
                        select_vars <- if(is.null(where_is_baseline)) c(1:dim(stats[[which_stats[i]]])[2]) else c(1:dim(stats[[which_stats[i]]])[2])[-where_is_baseline]
                    }
                    stats_loc <- if(length(select_vars)==1) array(stats[[which_stats[i]]][,select_vars,],dim=c(dim(stats[[which_stats[i]]])[1],1,dim(stats[[which_stats[i]]])[3])) else stats[[which_stats[i]]][,select_vars,]
                    remstimateList[[which_model[i]]]$diagnostics <- list()
                    remstimateList[[which_model[i]]]$diagnostics$residuals <- computeDiagnostics(pars = as.vector(optimum_model[[which_model[i]]]$argument)[select_vars],
                                                            stats = stats_loc,
                                                            actor1 = reh$edgelist$actor1-1,
                                                            actor2 = reh$edgelist$actor2-1,
                                                            dyad = c(0),
                                                            omit_dyad = reh$omit_dyad,
                                                            model = attr(reh,"model"),
                                                            N = reh$N,
                                                            senderRate = senderRate[i],
                                                            ncores = ncores,
                                                            baseline = baseline_value)                                      
                    colnames(remstimateList[[which_model[i]]]$diagnostics$residuals$standardized_residuals) <- if(senderRate[i]) variables_rate[select_vars] else variables_choice[select_vars]
                    colnames(remstimateList[[which_model[i]]]$diagnostics$residuals$smoothing_weights) <- colnames(remstimateList[[which_model[i]]]$diagnostics$residuals$standardized_residuals)
                    # lambdas (event rates)
                    remstimateList[[which_model[i]]]$diagnostics$rates <- remstimateList[[which_model[i]]]$diagnostics$residuals$rates   
                    remstimateList[[which_model[i]]]$diagnostics$residuals$rates <- NULL         
                    # residual deviance
                    remstimateList[[which_model[i]]]$residual.deviance <- -2*remstimateList[[which_model[i]]]$loglik
                    # null deviance                                    
                    remstimateList[[which_model[i]]]$null.deviance <- NULL

                    if(senderRate[i]){ # only for sender model
                        remstimateList$sender_model$null.deviance <- 2*(trust::trust(objfun = remDerivatives, 
                                                                                    parinit = c(0), 
                                                                                    rinit = 1, 
                                                                                    rmax = 100, 
                                                                                    stats = array(1,dim=c(dim(stats$sender_stats)[1],1,reh$M)),  
                                                                                    actor1 = reh$edgelist$actor1-1,
                                                                                    actor2 = c(0),
                                                                                    dyad = c(0),
                                                                                    omit_dyad = reh$omit_dyad,
                                                                                    interevent_time = reh$intereventTime,
                                                                                    model = model,
                                                                                    C = ifelse(is.null(reh$C),1,reh$C),
                                                                                    D = reh$D,
                                                                                    ordinal = ordinal,
                                                                                    senderRate = TRUE)$value)
                    }
                    else{ # for receiver model
                        remstimateList$receiver_model$null.deviance <- ifelse(length(reh$omit_dyad)>0,2*(trust::trust(objfun = remDerivatives, 
                                                                                                                    parinit = c(0), 
                                                                                                                    rinit = 1, 
                                                                                                                    rmax = 100, 
                                                                                                                    stats = array(1,dim=c(dim(stats$receiver_stats)[1],1,reh$M)), 
                                                                                                                    actor1 = reh$edgelist$actor1-1,
                                                                                                                    actor2 = reh$edgelist$actor2-1,
                                                                                                                    dyad = c(0),
                                                                                                                    omit_dyad = reh$omit_dyad, 
                                                                                                                    interevent_time = reh$intereventTime,
                                                                                                                    model = model,
                                                                                                                    senderRate = FALSE,
                                                                                                                    N = reh$N,
                                                                                                                    C = ifelse(is.null(reh$C),1,reh$C),
                                                                                                                    D = reh$D)$value),-2*log(1/(reh$N-1)))
                    }
                    
                    # [[NOTE: null.deviance.choice is for the receivder model -2*log(1/(reh$N-1)) when the riskset is always (all the time points) the same.
                    # The computation of the null.deviance can be simplified and the simplification should also account for dynamic riskset if omit_dyad is specified]]
                    
                    # model deviance
                    remstimateList[[which_model[i]]]$model.deviance <- remstimateList[[which_model[i]]]$null.deviance - remstimateList[[which_model[i]]]$residual.deviance
                    # df null
                    remstimateList[[which_model[i]]]$df.null <- reh$M
                    # df model
                    remstimateList[[which_model[i]]]$df.model <- if(i==1) length(variables_rate) else length(variables_choice)
                    # df residual
                    remstimateList[[which_model[i]]]$df.residual <- remstimateList[[which_model[i]]]$df.null - remstimateList[[which_model[i]]]$df.model      
                    # AIC
                    remstimateList[[which_model[i]]]$AIC <- 2*length(remstimateList[[which_model[i]]]$coefficients) - 2*remstimateList[[which_model[i]]]$loglik 
                    # AICC
                    remstimateList[[which_model[i]]]$AICC <- remstimateList[[which_model[i]]]$AIC + 2*length(remstimateList[[which_model[i]]]$coefficients)*(length(remstimateList[[which_model[i]]]$coefficients)+1)/(reh$M-length(remstimateList[[which_model[i]]]$coefficients)-1)
                    # BIC
                    remstimateList[[which_model[i]]]$BIC <- length(remstimateList[[which_model[i]]]$coefficients)*log(reh$M) - 2*remstimateList[[which_model[i]]]$loglik
                    # converged
                    remstimateList[[which_model[i]]]$converged <- optimum_model[[which_model[i]]]$converged
                    # iterations
                    remstimateList[[which_model[i]]]$iterations <- optimum_model[[which_model[i]]]$iterations 
                }
            }                 
        }            
        # ...
        # in future : if(WAIC) and if(ELPD){if(PSIS)} for predictive power of models
    }
   
    # ... estimating with Bayesian Sampling Importance Resampling (BSIR)
    if(method == "BSIR"){
        if(model == "tie"){ # Tie-oriented Model (Relational Event Model)
            bsir <- list()
            bsir$log_posterior <- rep(0,nsim)  # 
            bsir$draws <- list()  
            # (1) generate from the proposal (importance distribution)
            # First find location and scale parameters (mles and variances and covariances matrix) - MLE step
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

            # proposal distribution (default proposal is a multivariate Student t): drawing nsim*3 samples
            bsir$draws <- mvnfast::rmvt(n = nsim*3, 
                                            mu = mle_optimum$argument, 
                                            sigma = qr.solve(mle_optimum$hessian),
                                            df = 4,
                                            ncores = ncores) # df = 4 by default
                                        
            bsir$log_proposal <- mvnfast::dmvt(X = bsir$draws,  
                                                mu = mle_optimum$argument, 
                                                sigma = qr.solve(mle_optimum$hessian),
                                                df = 4,
                                                log = TRUE,
                                                ncores = ncores)
            # (2) calculate log-posterior density give by log_prior + loglik
            ## (2.1) summing log_prior (if NULL it will be non-informative) : prior by default can be a multivariate Student t and the user must specify its parameters
            if(!is.null(prior)){
                bsir$log_prior <- do.call(prior,c(list(bsir$draws),list_args)) 
                bsir$log_posterior <- bsir$log_prior
            }
            ## (2.2) summing loglik 
            loglik <- apply(bsir$draws,1,function(x) (-1)*remDerivativesStandard(pars = x,
                                    stats = stats,
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
            bsir$df.model <- dim(stats)[2]
            bsir$df.residual <- bsir$df.null - bsir$df.model
            names(bsir$coefficients) <- names(bsir$post.mean) <- rownames(bsir$vcov) <- colnames(bsir$vcov) <- names(bsir$sd) <- dimnames(stats)[[2]]
            # residuals
            select_vars <- if(is.null(where_is_baseline)) 1:length(variable_names) else c(1:length(variable_names))[-where_is_baseline]
            stats_loc <- if(length(select_vars)==1) array(stats[,select_vars,],dim=c(dim(stats)[1],1,dim(stats)[3])) else stats[,select_vars,]
            baseline_value <- if(is.null(where_is_baseline)) 0 else bsir$coefficients[where_is_baseline]
            bsir$diagnostics <- list()
            bsir$diagnostics$residuals <- computeDiagnostics(pars = bsir$coefficients[select_vars],
                                                    stats = stats_loc,
                                                    actor1 = c(0),
                                                    actor2 = c(0),
                                                    dyad = attr(reh,"dyad")-1,
                                                    omit_dyad = reh$omit_dyad,
                                                    model = attr(reh,"model"),
                                                    ncores = ncores,
                                                    baseline = baseline_value,
                                                    N = reh$N)
            colnames(bsir$diagnostics$residuals$standardized_residuals) <- variable_names[select_vars]  
            colnames(bsir$diagnostics$residuals$smoothing_weights) <- colnames(bsir$diagnostics$residuals$standardized_residuals)  
            # lambdas (event rates)
            bsir$diagnostics$rates <- bsir$diagnostics$residuals$rates   
            bsir$diagnostics$residuals$rates <- NULL  

            remstimateList <- bsir
        }
        if(model == "actor"){ # Actor-oriented Model
            bsir <- list()
            senderRate <- c(TRUE,FALSE) # meaning that the first BSIR will be run on the "sender rate" model and the second on the "receiver choice"
            which_model <- c("sender_model","receiver_model")
            which_stats <- c("sender_stats","receiver_stats")
            for(i in 1:2){
                if(!is.null(stats[[which_stats[i]]])){  # run BSIR only if the model is specified
                    bsir_i <- list()
                    bsir_i$log_posterior <- rep(0,nsim)
                    bsir_i$draws <- list()  

                    # (1) generate from the proposal (importance distribution)

                    # First find location and scale parameters (mles and variances and covariances matrix)
                    mle_optimum <- trust::trust(objfun = remDerivatives, 
                                                        parinit = rep(0,dim(stats[[which_stats[i]]])[2]), 
                                                        rinit = 1, 
                                                        rmax = 100, 
                                                        stats = stats[[which_stats[i]]], 
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
                                                        C = ifelse(is.null(reh$C),1,reh$C),
                                                        D = reh$D)                          
                    
                    # proposal distribution (default proposal is a multivariate Student t): drawing nsim*3 samples
                    bsir_i$draws <- mvnfast::rmvt(n = nsim*3, 
                                                    mu = mle_optimum$argument, 
                                                    sigma = qr.solve(mle_optimum$hessian),
                                                    df = 4,
                                                    ncores = ncores) # df = 4 by default
                                                
                    bsir_i$log_proposal <- mvnfast::dmvt(X = bsir_i$draws,  
                                                        mu = mle_optimum$argument, 
                                                        sigma = qr.solve(mle_optimum$hessian),
                                                        df = 4,
                                                        log = TRUE,
                                                        ncores = ncores)

                    # (2) calculate log-posterior density give by log_prior + loglik
                    ## (2.1) summing log_prior (if NULL it will be non-informative) : prior by default can be a multivariate Student t and the user must specify its parameters
                    if(!is.null(prior)){
                        bsir_i$log_prior <- do.call(prior[[which_model[i]]],c(list(bsir_i$draws),list_args[[which_model[i]]]))
                        bsir_i$log_posterior <- bsir_i$log_prior
                    }
                    ## (2.2) summing loglik                   
                    loglik <- apply(bsir_i$draws,1,function(x) (-1)*remDerivatives(pars = x, 
                                    stats = stats[[which_stats[i]]], 
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
                                    C = ifelse(is.null(reh$C),1,reh$C),
                                    D = reh$D)$value) # this is the most time consuming step to be optimized, avoiding the apply        
                    bsir_i$log_posterior <- bsir_i$log_posterior + loglik

                    # (3) calculate Importance resampling log-weights (irlw)
                    irlw <- (bsir_i$log_posterior - max(bsir_i$log_posterior) + min(bsir_i$log_proposal)) - bsir_i$log_proposal
                    s_irlw <- log(sum(exp(irlw)) - exp(irlw)) # ISIR (Skare et al., 2003 Scandinavian Journal of Statistics)
                    irlw <- irlw - s_irlw
                    bsir_i$irw <- exp(irlw) # importance resampling weights (irw)
                    #if(any(is.na(irw)|is.nan(irw)|is.na(exp(irw)))||all(exp(irw)==0)) --> Importance weights might present some issues.
                    post_draws <-sample(x=1:dim(bsir_i$draws)[1],size=nsim,replace=TRUE,prob=exp(irlw)) # we always draw nsim (having generated 3*nsim samples from the proposal)
                    bsir_i$draws <- bsir_i$draws[post_draws,]
                    bsir_i$log_posterior <- bsir_i$log_posterior[post_draws]
                    bsir_i$coefficients <- bsir_i$draws[which.max(bsir_i$log_posterior),]
                    bsir_i$loglik <- max(bsir_i$log_posterior)
                    bsir_i$post.mean <- colMeans(bsir_i$draws)
                    bsir_i$vcov <- stats::cov(bsir_i$draws)
                    bsir_i$sd <- diag(bsir_i$vcov)**0.5
                    bsir_i$df.null <- reh$M
                    bsir_i$df.model <- dim(stats[[which_stats[i]]])[2]
                    bsir_i$df.residual <- bsir_i$df.null - bsir_i$df.model
                    names(bsir_i$coefficients) <- names(bsir_i$post.mean) <- rownames(bsir_i$vcov) <- colnames(bsir_i$vcov) <- names(bsir_i$sd) <- if(i == 1) variables_rate else variables_choice
                    # residuals
                    baseline_value <- 0
                    select_vars <- c(1:dim(stats[[which_stats[i]]])[2])
                    if(senderRate[i]){
                        baseline_value <- if(is.null(where_is_baseline)) 0 else bsir_i$coefficients[where_is_baseline]
                        select_vars <- if(is.null(where_is_baseline)) c(1:dim(stats[[which_stats[i]]])[2]) else c(1:dim(stats[[which_stats[i]]])[2])[-where_is_baseline]
                    }
                    stats_loc <- if(length(select_vars)==1) array(stats[[which_stats[i]]][,select_vars,],dim=c(dim(stats[[which_stats[i]]])[1],1,dim(stats[[which_stats[i]]])[3])) else stats[[which_stats[i]]][,select_vars,]
                    bsir_i$diagnostics <- list()
                    bsir_i$diagnostics$residuals <- computeDiagnostics(pars = bsir_i$coefficients[select_vars],
                                                            stats = stats_loc,
                                                            actor1 = reh$edgelist$actor1-1,
                                                            actor2 = reh$edgelist$actor2-1,
                                                            dyad = c(0),
                                                            omit_dyad = reh$omit_dyad,
                                                            model = attr(reh,"model"),
                                                            N = reh$N,
                                                            senderRate = senderRate[i],
                                                            ncores = ncores,
                                                            baseline = baseline_value)
                    colnames(bsir_i$diagnostics$residuals$standardized_residuals) <- if(senderRate[i]) variables_rate[select_vars] else variables_choice[select_vars]
                    colnames(bsir_i$diagnostics$residuals$smoothing_weights) <- colnames(bsir_i$diagnostics$residuals$standardized_residuals)
                    # lambdas (event rates)
                    bsir_i$diagnostics$rates <- bsir_i$diagnostics$residuals$rates   
                    bsir_i$diagnostics$residuals$rates <- NULL         
                    bsir[[which_model[i]]] <- bsir_i
                    # freeing some memory
                    rm(bsir_i,mle_optimum,loglik,irlw,s_irlw,post_draws)
                }
            }
            remstimateList <- bsir
        }
    }

    # ... estimating with Hamiltonian Monte Carlo (HMC)
    if(method == "HMC"){
        hmc <- list() # create an empty list where to store the output
        if(model == "tie"){ # Tie-oriented Model (Relational Event Model)
            if(is.null(init)){
                beta_0 <- log(reh$M) - log(sum(reh$intereventTime)*reh$D) # MLE intercept under null model
                init <- matrix(stats::runif((dim(stats)[2])*nchains,-0.1,0.1),nrow=dim(stats)[2],ncol=nchains)
                if(any(dimnames(stats)[[2]] == "baseline")){
                    init[which(dimnames(stats)[[2]] == "baseline"),] <- init[which(dimnames(stats)[[2]] == "baseline"),] + beta_0
                }
            }
            set.seed(seed)
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
                        L = L,
                        epsilon = epsilon, 
                        N = reh$N,
                        C = ifelse(is.null(reh$C),1,reh$C),
                        D = reh$D)          
            hmc$coefficients <- hmc_out$draws[which.min(hmc_out$log_posterior),] # log_posterior is a vector of posterior nloglik, therefore we take the minimum
            hmc$post.mean <- colMeans(hmc_out$draws)
            hmc$vcov <- stats::cov(hmc_out$draws)
            hmc$sd <- diag(hmc$vcov)**0.5
            hmc$loglik <- -min(hmc_out$log_posterior) # max loglik is -min(log_posterior)
            hmc$df.null <- reh$M
            hmc$df.model <- dim(stats)[2]
            hmc$df.residual <- hmc$df.null - hmc$df.model
            names(hmc$coefficients) <- names(hmc$post.mean) <- rownames(hmc$vcov) <- colnames(hmc$vcov) <- names(hmc$sd) <- dimnames(stats)[[2]]
            # residuals
            select_vars <- if(is.null(where_is_baseline)) 1:length(variable_names) else c(1:length(variable_names))[-where_is_baseline]
            stats_loc <- if(length(select_vars)==1) array(stats[,select_vars,],dim=c(dim(stats)[1],1,dim(stats)[3])) else stats[,select_vars,]
            baseline_value <- if(is.null(where_is_baseline)) 0 else hmc$coefficients[where_is_baseline]
            hmc$diagnostics <- list()
            hmc$diagnostics$residuals <- computeDiagnostics(pars = hmc$coefficients[select_vars],
                                                    stats = stats_loc,
                                                    actor1 = c(0),
                                                    actor2 = c(0),
                                                    dyad = attr(reh,"dyad")-1,
                                                    omit_dyad = reh$omit_dyad,
                                                    model = attr(reh,"model"),
                                                    ncores = ncores,
                                                    baseline = baseline_value,
                                                    N = reh$N)
            colnames(hmc$diagnostics$residuals$standardized_residuals) <- variable_names[select_vars]     
            colnames(hmc$diagnostics$residuals$smoothing_weights) <- colnames(hmc$diagnostics$residuals$standardized_residuals)
            # lambdas (event rates)
            hmc$diagnostics$rates <- hmc$diagnostics$residuals$rates   
            hmc$diagnostics$residuals$rates <- NULL  

            remstimateList <- c(hmc_out,hmc)
        }
        if(model == "actor"){ # Actor-oriented Model
            senderRate <- c(TRUE,FALSE) # meaning that the first BSIR will be run on the "sender rate" model and the second on the "receiver choice"
            which_model <- c("sender_model","receiver_model")
            which_stats <- c("sender_stats","receiver_stats")
            for(i in 1:2){
                if(!is.null(stats[[which_stats[i]]])){
                    hmc_i <- list()  
                    if(is.null(init[[which_model[i]]])){
                        init[[which_model[i]]] <- matrix(stats::runif((dim(stats[[which_stats[i]]])[2])*nchains,-0.1,0.1),nrow=dim(stats[[which_stats[i]]])[2],ncol=nchains)
                    }
                    set.seed(seed)
                    hmc_out_i <- HMC(pars_init = init[[which_model[i]]], 
                                nsim = (burnin+nsim), 
                                nchains = nchains, 
                                burnin = burnin, 
                                meanPrior = rep(0,dim(stats[[which_stats[i]]])[2]),
                                sigmaPrior = diag(dim(stats[[which_stats[i]]])[2]),
                                stats = stats[[which_stats[i]]],  
                                actor1 = reh$edgelist$actor1-1,
                                actor2 = reh$edgelist$actor2-1,
                                dyad = c(0), 
                                omit_dyad = reh$omit_dyad,
                                interevent_time = reh$intereventTime,
                                model = model,
                                senderRate = senderRate[i],
                                N = reh$N,
                                C = ifelse(is.null(reh$C),1,reh$C),
                                D = reh$D,
                                ordinal = ordinal,
                                ncores = ncores,
                                thin = thin,
                                L = L,
                                epsilon = epsilon)   
                    hmc_i$coefficients <- hmc_out_i$draws[which.min(hmc_out_i$log_posterior),] # log_posterior is a vector of posterior nloglik, therefore we take the minimum
                    hmc_i$post.mean <- colMeans(hmc_out_i$draws)
                    hmc_i$vcov <- stats::cov(hmc_out_i$draws)
                    hmc_i$sd <- diag(hmc_i$vcov)**0.5
                    hmc_i$loglik <- -min(hmc_out_i$log_posterior) # max loglik is -min(log_posterior)
                    hmc_i$df.null <- reh$M
                    hmc_i$df.model <- dim(stats[[which_stats[i]]])[2]
                    hmc_i$df.residual <- hmc_i$df.null - hmc_i$df.model
                    names(hmc_i$coefficients) <- names(hmc_i$post.mean) <- rownames(hmc_i$vcov) <- colnames(hmc_i$vcov) <- names(hmc_i$sd) <- if(i == 1) variables_rate else variables_choice
                    # residuals
                    baseline_value <- 0
                    select_vars <- c(1:dim(stats[[which_stats[i]]])[2])
                    if(senderRate[i]){
                        baseline_value <- if(is.null(where_is_baseline)) 0 else hmc_i$coefficients[where_is_baseline]
                        select_vars <- if(is.null(where_is_baseline)) c(1:dim(stats[[which_stats[i]]])[2]) else c(1:dim(stats[[which_stats[i]]])[2])[-where_is_baseline]
                    }
                    stats_loc <- if(length(select_vars)==1) array(stats[[which_stats[i]]][,select_vars,],dim=c(dim(stats[[which_stats[i]]])[1],1,dim(stats[[which_stats[i]]])[3])) else stats[[which_stats[i]]][,select_vars,]
                    hmc_i$diagnostics <- list()
                    hmc_i$diagnostics$residuals <- computeDiagnostics(pars = hmc_i$coefficients[select_vars],
                                                            stats = stats_loc,
                                                            actor1 = reh$edgelist$actor1-1,
                                                            actor2 = reh$edgelist$actor2-1,
                                                            dyad = c(0),
                                                            omit_dyad = reh$omit_dyad,
                                                            model = attr(reh,"model"),
                                                            N = reh$N,
                                                            senderRate = senderRate[i],
                                                            ncores = ncores,
                                                            baseline = baseline_value)
                    colnames(hmc_i$diagnostics$residuals$standardized_residuals) <- if(senderRate[i]) variables_rate[select_vars] else variables_choice[select_vars]
                    colnames(hmc_i$diagnostics$residuals$smoothing_weights) <- colnames(hmc_i$diagnostics$residuals$standardized_residuals)
                    # lambdas (event rates)
                    hmc_i$diagnostics$rates <- hmc_i$diagnostics$residuals$rates   
                    hmc_i$diagnostics$residuals$rates <- NULL  

                    hmc[[which_model[i]]] <- c(hmc_out_i,hmc_i)
                    rm(hmc_out_i,hmc_i)
                }
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
            attr(str_out, "init") <- init
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
            attr(str_out, "init") <- init 
        }
    }
    return(str_out)
}

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
#' 
#' @examples
#' 
#' # ------------------------------------ #
#' #       method 'print' for the         #
#' #       tie-oriented model: "BSIR"     #
#' # ------------------------------------ #
#' 
#' # loading data
#' data(tie_reh)
#'   
#' # specifying linear predictor
#' tie_model <- ~ 1 + 
#'                remstats::indegreeSender()+
#'                remstats::inertia()+
#'                remstats::reciprocity() 
#' 
#' # calculating statistics
#' tie_reh_stats <- remstats::remstats(reh = tie_reh, 
#'                                     tie_effects = tie_model)
#' 
#' # running estimation
#' tie_mle <- remstimate::remstimate(reh = tie_reh,
#'                                   stats = tie_reh_stats,
#'                                   method = "BSIR",
#'                                   nsim = 100,
#'                                   ncores = 1)
#' 
#' # print
#' tie_mle
#' 
#' # ------------------------------------ #
#' #      method 'print' for the          #
#' #      actor-oriented model: "BSIR"    #
#' # ------------------------------------ #
#' 
#' # loading data
#' data(ao_reh)
#'   
#' # specifying linear predictor (for sender rate and receiver choice model)
#' rate_model <- ~ 1 + remstats::indegreeSender()
#' choice_model <- ~ remstats::inertia() + remstats::reciprocity()
#' 
#' # calculating statistics
#' ao_reh_stats <- remstats::remstats(reh = ao_reh, 
#'                                    sender_effects = rate_model, 
#'                                    receiver_effects = choice_model)
#' 
#' # running estimation
#' ao_mle <- remstimate::remstimate(reh = ao_reh,
#'                                  stats = ao_reh_stats,
#'                                  method = "BSIR",
#'                                  nsim = 100,
#'                                  ncores = 1)
#' 
#' # print
#' ao_mle
#' 
#' # ------------------------------------ #
#' #   for more examples check vignettes  #
#' # ------------------------------------ #
#' 
print.remstimate<-function(x, ...){
    if (is.null(attr(x,"formula"))) 
        stop("invalid 'remstimate' object:  no 'formula' attribute")
    #if (!inherits(x, "remstimate")) 
    #    warning("calling print.remstimate(<fake-remstimate-object>) ...") 
    cat("Relational Event Model",paste("(",attr(x,"model")," oriented)",sep=""),"\n")
    if(attr(x,"approach") == "Frequentist"){
        if(attr(x,"model") == "tie"){
            cat("\nCoefficients:\n\n")
            print(x$coefficients)
            cat("\nNull deviance:",x$null.deviance,"\nResidual deviance:",x$residual.deviance,"\n")
            cat("AIC:",x$AIC,"AICC:",x$AICC,"BIC:",x$BIC,"\n\n")
        }
        else if(attr(x,"model") == "actor"){
            if(!is.null(x$sender_model)){  # printing sender model 
                cat("\nCoefficients rate model **for sender**:\n\n")
                print(x$sender_model$coefficients)
                cat("\nNull deviance:",x$sender_model$null.deviance,"\nResidual deviance:",x$sender_model$residual.deviance,"\n")
                cat("AIC:",x$sender_model$AIC,"AICC:",x$sender_model$AICC,"BIC:",x$sender_model$BIC,"\n")
            }
            if(!is.null(x$sender_model) & !is.null(x$receiver_model)){ # if both models are estimated, separate the two print by a "-"
                cat(paste0(rep("-", getOption("width")),collapse = ""))
            }
            if(!is.null(x$receiver_model)){ # printing receiver model
                cat("\n\nCoefficients choice model **for receiver**:\n\n")
                print(x$receiver_model$coefficients)
                cat("\nNull deviance:",x$receiver_model$null.deviance,"\nResidual deviance:",x$receiver_model$residual.deviance,"\n")
                cat("AIC:",x$receiver_model$AIC,"AICC:",x$receiver_model$AICC,"BIC:",x$receiver_model$BIC,"\n\n")
            }   
        }

    }
    else{ # Bayesian
        if(attr(x,"model") == "tie"){
            cat("\nPosterior Modes:\n\n")
            print(x$coefficients)
        }
        else if(attr(x,"model") == "actor"){
            if(!is.null(x$sender_model)){  # printing sender model 
                cat("\nPosterior Modes rate model **for sender**:\n\n")
                print(x$sender_model$coefficients)
                cat("\n")
                # write some more info here
            }
            if(!is.null(x$sender_model) & !is.null(x$receiver_model)){ # if both models are estimated, separate the two print by a "-"
                cat(paste0(rep("-", getOption("width")),collapse = ""))
            }
            if(!is.null(x$receiver_model)){ # printing receiver model
                cat("\n\nPosterior Modes choice model **for sender**:\n\n")
                print(x$receiver_model$coefficients)
                # write some more info here
            }
        }
    }
}


#######################################################################################
#######################################################################################


# summary.remstimate
#' @title Generate the summary of a remstimate object
#' @rdname summary.remstimate
#' @description A function that prints out the summary of a remstimate object
#' @param object is a \code{remstimate} object 
#' @param ... further arguments to be passed to the 'summary' method
#' @method summary remstimate
#' @export
#' 
#' @examples
#' 
#' # ------------------------------------ #
#' #       method 'summary' for the       #
#' #       tie-oriented model: "BSIR"     #
#' # ------------------------------------ #
#' 
#' # loading data
#' data(tie_reh)
#'   
#' # specifying linear predictor
#' tie_model <- ~ 1 + 
#'                remstats::indegreeSender()+
#'                remstats::inertia()+
#'                remstats::reciprocity() 
#' 
#' # calculating statistics
#' tie_reh_stats <- remstats::remstats(reh = tie_reh, 
#'                                     tie_effects = tie_model)
#' 
#' # running estimation
#' tie_mle <- remstimate::remstimate(reh = tie_reh,
#'                                   stats = tie_reh_stats,
#'                                   method = "BSIR",
#'                                   nsim = 100,
#'                                   ncores = 1)
#' 
#' # summary
#' summary(tie_mle)
#' 
#' # ------------------------------------ #
#' #      method 'summary' for the        #
#' #      actor-oriented model: "BSIR"    #
#' # ------------------------------------ #
#' 
#' # loading data
#' data(ao_reh)
#'   
#' # specifying linear predictor (for sender rate and receiver choice model)
#' rate_model <- ~ 1 + remstats::indegreeSender()
#' choice_model <- ~ remstats::inertia() + remstats::reciprocity()
#' 
#' # calculating statistics
#' ao_reh_stats <- remstats::remstats(reh = ao_reh, 
#'                                    sender_effects = rate_model, 
#'                                    receiver_effects = choice_model)
#' 
#' # running estimation
#' ao_mle <- remstimate::remstimate(reh = ao_reh,
#'                                  stats = ao_reh_stats,
#'                                  method = "BSIR",
#'                                  nsim = 100,
#'                                  ncores = 1)
#' 
#' # summary
#' summary(ao_mle)
#' 
#' # ------------------------------------ #
#' #   for more examples check vignettes  #
#' # ------------------------------------ #
#' 
summary.remstimate<-function (object, ...)
{
    # (codings are based on the structure of summary.lm() and summary.glm() from package 'stats')
    #if (!inherits(object, "remstimate")) 
    #    warning("calling summary.remstimate(<fake-remstimate-object>) ...")
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
            if(!is.null(object$sender_model)){ # summary sender model
                coefsTab$sender_model <- cbind(object$sender_model$coefficients,
                object$sender_model$se,
                object$sender_model$coefficients/object$sender_model$se,
                2*(1-stats::pnorm(abs(object$sender_model$coefficients/object$sender_model$se))),
                (stats::dnorm(0,mean=object$sender_model$coefficients,sd=object$sender_model$se) / stats::dnorm(0,mean=0,sd=object$sender_model$se*sqrt(object$sender_model$df.null)))/((stats::dnorm(0,mean=object$sender_model$coefficients,sd=object$sender_model$se) / stats::dnorm(0,mean=0,sd=object$sender_model$se*sqrt(object$sender_model$df.null)))+1)
                )
                colnames(coefsTab$sender_model) <- c("Estimate","Std. Err", "z value", "Pr(>|z|)", "Pr(=0)")
                rownames(coefsTab$sender_model) <- names(object$sender_model$coefficients)
            }
            if(!is.null(object$receiver_model)){ # summary receiver model
                coefsTab$receiver_model <- cbind(object$receiver_model$coefficients,
                object$receiver_model$se,
                object$receiver_model$coefficients/object$receiver_model$se,
                2*(1-stats::pnorm(abs(object$receiver_model$coefficients/object$receiver_model$se))),
                (stats::dnorm(0,mean=object$receiver_model$coefficients,sd=object$receiver_model$se) / stats::dnorm(0,mean=0,sd=object$receiver_model$se*sqrt(object$receiver_model$df.null)))/((stats::dnorm(0,mean=object$receiver_model$coefficients,sd=object$receiver_model$se) / stats::dnorm(0,mean=0,sd=object$receiver_model$se*sqrt(object$receiver_model$df.null)))+1)
                )
                colnames(coefsTab$receiver_model) <- c("Estimate","Std. Err", "z value", "Pr(>|z|)", "Pr(=0)")
                rownames(coefsTab$receiver_model) <- names(object$receiver_model$coefficients)
            }
            summary_out$coefsTab <- coefsTab
            #keep <- match(c("formula",
		    #  "contrasts", "sender_model", "receiver_model", "na.action"), names(object), 0L) # check if it works without this line
            keep <- list(formula = attr(object,"formula"),
                        model = attr(object,"model"),
                        ordinal = attr(object,"ordinal"),
                        method = attr(object, "method"),
                        approach = attr(object, "approach"),
                        epsilon = attr(object, "epsilon"))        
            if(!is.null(object$sender_model)){
                keep$sender_model <- list(residual.deviance = object$sender_model$residual.deviance, 
                                            null.deviance = object$sender_model$null.deviance,
                                            model.deviance = object$sender_model$model.deviance,
                                            df.residual = object$sender_model$df.residual,
                                            df.null = object$sender_model$df.null,
                                            df.model = object$sender_model$df.model,
                                            AIC = object$sender_model$AIC,
                                            AICC = object$sender_model$AICC,
                                            BIC = object$sender_model$BIC
                                        )
            }  
            if(!is.null(object$receiver_model)){
                keep$receiver_model <- list(residual.deviance = object$receiver_model$residual.deviance, 
                                            null.deviance = object$receiver_model$null.deviance,
                                            model.deviance = object$receiver_model$model.deviance,
                                            df.residual = object$receiver_model$df.residual,
                                            df.null = object$receiver_model$df.null,
                                            df.model = object$receiver_model$df.model,
                                            AIC = object$receiver_model$AIC,
                                            AICC = object$receiver_model$AICC,
                                            BIC = object$receiver_model$BIC
                                        )
            }           
            summary_out <- do.call(c, list(keep, summary_out))
        }
    }
    else if(attr(object, "approach") == "Bayesian"){ # Bayesian
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
            if(!is.null(object$sender_model)){ # summary sender model
                coefsTab$sender_model <- cbind(object$sender_model$coefficients,
                object$sender_model$sd,
                apply(object$sender_model$draws,2,stats::quantile,0.025),
                apply(object$sender_model$draws,2,stats::quantile,0.5),
                apply(object$sender_model$draws,2,stats::quantile,0.975),
                (stats::dnorm(0,mean=object$sender_model$coefficients,sd=object$sender_model$sd) / stats::dnorm(0,mean=0,sd=object$sender_model$sd*sqrt(object$sender_model$df.null)))/((stats::dnorm(0,mean=object$sender_model$coefficients,sd=object$sender_model$sd) / stats::dnorm(0,mean=0,sd=object$sender_model$sd*sqrt(object$sender_model$df.null)))+1)
                )
                colnames(coefsTab$sender_model) <- c("Post.Mode","Post.SD","Q2.5%","Q50%","Q97.5%","Pr(=0|x)")
                rownames(coefsTab$sender_model) <- names(object$sender_model$coefficients)
            }
            if(!is.null(object$receiver_model)){ # summary receiver model
                coefsTab$receiver_model <- cbind(object$receiver_model$coefficients,
                object$receiver_model$sd,
                apply(object$receiver_model$draws,2,stats::quantile,0.025),
                apply(object$receiver_model$draws,2,stats::quantile,0.5),
                apply(object$receiver_model$draws,2,stats::quantile,0.975),
                (stats::dnorm(0,mean=object$receiver_model$coefficients,sd=object$receiver_model$sd) / stats::dnorm(0,mean=0,sd=object$receiver_model$sd*sqrt(object$receiver_model$df.null)))/((stats::dnorm(0,mean=object$receiver_model$coefficients,sd=object$receiver_model$sd) / stats::dnorm(0,mean=0,sd=object$receiver_model$sd*sqrt(object$receiver_model$df.null)))+1)
                )
                colnames(coefsTab$receiver_model) <- c("Post.Mode","Post.SD","Q2.5%","Q50%","Q97.5%","Pr(=0|x)")
                rownames(coefsTab$receiver_model) <- names(object$receiver_model$coefficients)
            }
            summary_out$coefsTab <- coefsTab
            #keep <- match(c("formula",
		    #  "contrasts", "sender_model", "receiver_model", "na.action"), names(object), 0L)
            keep <- list(formula = attr(object,"formula"),
                        model = attr(object,"model"),
                        ordinal = attr(object,"ordinal"),
                        method = attr(object, "method"),
                        approach = attr(object, "approach"),
                        epsilon = attr(object, "epsilon")) # add other attributes from the objects?     
            if(!is.null(object$sender_model)){
                keep$sender_model <- list(
                            #residual.deviance = object$sender_model$residual.deviance, 
                            #null.deviance = object$sender_model$null.deviance,
                            #model.deviance = object$sender_model$model.deviance,
                            #df.residual = object$sender_model$df.residual,
                            df.null = object$sender_model$df.null,
                            #df.model = object$sender_model$df.model,
                            loglik = object$sender_model$loglik
                        )
            }
            if(!is.null(object$receiver_model)){
                keep$receiver_model <- list(
                            #residual.deviance = object$receiver_model$residual.deviance, 
                            #null.deviance = object$receiver_model$null.deviance,
                            #model.deviance = object$receiver_model$model.deviance,
                            #df.residual = object$receiver_model$df.residual,
                            df.null = object$receiver_model$df.null,
                            #df.model = object$receiver_model$df.model,
                            loglik = object$receiver_model$loglik
                        )  
            }                
            summary_out <- do.call(c, list(keep, summary_out))
        }
    }

    if(length(summary_out)==0) stop("invalid 'remstimate' object")

    cat("Relational Event Model",paste("(",summary_out$model," oriented)",sep=""),"\n\n")
    if(summary_out$model == "tie"){
        cat("Call:\n",deparse(summary_out$formula),"\n\n",sep="")
        second_line <- paste("(",summary_out$method," with ",sep="")
        if(summary_out$ordinal) second_line <- paste(second_line,"ordinal likelihood):",sep="")
        else{
            second_line <- paste(second_line,"interval likelihood):\n\n",sep="")
        }
        if(summary_out$approach == "Frequentist"){
            cat("\nCoefficients",second_line)
        }
        else{ # Bayesian
            cat("\nPosterior Modes",second_line)
        }
        stats::printCoefmat(summary_out$coefsTab, P.values = TRUE, signif.stars = FALSE, ...) 
        if(summary_out$approach == "Frequentist"){
            cat("Null deviance:", summary_out$null.deviance, "on", summary_out$df.null, "degrees of freedom\n")
            cat("Residual deviance:", summary_out$residual.deviance, "on", summary_out$df.residual, "degrees of freedom\n")
            cat("Chi-square:", summary_out$model.deviance, "on", summary_out$df.model, 
                "degrees of freedom, asymptotic p-value", 1 - stats::pchisq(summary_out$model.deviance, 
                    summary_out$df.model), "\n")
            cat("AIC:", summary_out$AIC, "AICC:", summary_out$AICC, "BIC:", summary_out$BIC, "\n")
        }
        if(summary_out$approach == "Bayesian"){
        cat("Log posterior:",summary_out$loglik,"\n")
        # cat("Prior parameters:",paste(names(summary_out$prior.param),unlist(summary_out$prior.param),sep="="),"\n")
        }
    }
    else if(summary_out$model == "actor"){
        second_line <- paste("(",summary_out$method," with ",sep="")
        if(summary_out$ordinal){
            second_line <- paste(second_line,"ordinal likelihood):",sep="")
        }
        else{
            second_line <- paste(second_line,"interval likelihood):\n\n",sep="")
        }
        if(!is.null(summary_out$sender_model)){
            # sender rate summary
            cat("Call rate model **for sender**:\n\n\t",deparse(summary_out$formula$rate_model_formula),"\n\n",sep="")
            if(summary_out$approach == "Frequentist"){
                cat("\nCoefficients rate model",second_line)
            }
            else{ # Bayesian
                cat("\nPosterior Modes rate model",second_line)
            }
            stats::printCoefmat(summary_out$coefsTab$sender_model, P.values = TRUE, signif.stars = FALSE)
            if(summary_out$approach == "Frequentist"){
                cat("Null deviance:", summary_out$sender_model$null.deviance, "on", summary_out$sender_model$df.null, "degrees of freedom\n")
                cat("Residual deviance:", summary_out$sender_model$residual.deviance, "on", summary_out$sender_model$df.residual, "degrees of freedom\n")
                cat("Chi-square:", summary_out$sender_model$model.deviance, "on", summary_out$sender_model$df.model, 
                    "degrees of freedom, asymptotic p-value", 1 - stats::pchisq(summary_out$sender_model$model.deviance, 
                        summary_out$sender_model$df.model), "\n")
                cat("AIC:", summary_out$sender_model$AIC, "AICC:", summary_out$sender_model$AICC, "BIC:", summary_out$sender_model$BIC, "\n")
            }
            if(summary_out$approach == "Bayesian"){
            cat("Log posterior:",summary_out$sender_model$loglik,"\n")
            # cat("Prior parameters:",paste(names(summary_out$sender_model$prior.param),unlist(summary_out$sender_model$prior.param),sep="="),"\n")
            }
        }
        if(!is.null(summary_out$sender_model) & !is.null(summary_out$receiver_model)){ # if both models are estimated, separate the two print by a "-"
            cat(paste0(rep("-", getOption("width")),collapse = ""),"\n\n")
        }
        if(!is.null(summary_out$receiver_model)){
            # receiver choice summary 
            cat("Call choice model **for receiver**:\n\n\t",deparse(summary_out$formula$choice_model_formula),"\n\n",sep="") 
                    if(summary_out$approach == "Frequentist"){
                cat("\nCoefficients choice model",second_line)
            }
            else{ # Bayesian
                cat("\nPosterior Modes choice model",second_line)
            }
            stats::printCoefmat(summary_out$coefsTab$receiver_model, P.values = TRUE, signif.stars = FALSE)
            if(summary_out$approach == "Frequentist"){
                cat("Null deviance:", summary_out$receiver_model$null.deviance, "on", summary_out$receiver_model$df.null, "degrees of freedom\n")
                cat("Residual deviance:", summary_out$receiver_model$residual.deviance, "on", summary_out$receiver_model$df.residual, "degrees of freedom\n")
                cat("Chi-square:", summary_out$receiver_model$model.deviance, "on", summary_out$receiver_model$df.model, 
                    "degrees of freedom, asymptotic p-value", 1 - stats::pchisq(summary_out$receiver_model$model.deviance, 
                        summary_out$receiver_model$df.model), "\n")
                cat("AIC:", summary_out$receiver_model$AIC, "AICC:", summary_out$receiver_model$AICC, "BIC:", summary_out$receiver_model$BIC, "\n")
            }
            if(summary_out$approach == "Bayesian"){
            cat("Log posterior:",summary_out$receiver_model$loglik,"\n")
            # cat("Prior parameters:",paste(names(summary_out$receiver_model$prior.param),unlist(summary_out$receiver_model$prior.param),sep="="),"\n")
            }
        }
    }

}

#######################################################################################
#######################################################################################


# residuals.remstimate
#' @title residuals.remstimate
#' @rdname residuals.remstimate
#' @description A function that extracts model residuals (Schoenfeld's residuals - add citation) from a 'remstimate' object.
#' @param object is a \code{remstimate} object 
#' @param ... further arguments to be passed to the 'residuals' method depending on the estimator used
#' @method residuals remstimate
#' @export
#' 
#' @examples 
#' 
#' # No examples available at the moment
#' 
residuals.remstimate <- function(object, ...)
{
 # return residuals
 if(attr(object,"model") == "tie"){
    return(object$diagnostics$residuals)
 }
 else{
    return(list(sender_model = object$sender_model$diagnostics$residuals, receiver_model = object$receiver_model$diagnostics$residuals))
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
#' 
#' @examples 
#' 
#' # No examples available at the moment
#' 
predict.remstimate <- function(object, ...)
{
    # write a predict method here
    # - in-sample predictive performance 
    # - out-of-sample predictions
    return(paste('this function at the moment does nothing'))
}

#######################################################################################
#######################################################################################


# plot.remstimate
#' @title plot.remstimate
#' @rdname plot.remstimate
#' @description A function that returns a plot of diagnostics given a 'remstimate' object and depending on the 'approach' attribute.
#' @param x is a \code{remstimate} object
#' @param reh 'remify' object 
#' @param ... further arguments to be passed to the 'predict' method depending on the estimator used
#' @method plot remstimate
#' @export
#' 
#' @examples 
#' 
#' # No examples available at the moment
#' 
plot.remstimate <- function(x, reh, ...)
{
    if(attr(x,"model") != attr(reh,"model")){
        stop("'x' and 'reh' have different attribute 'model'")
    }
    if(attr(x,"model") == "tie"){ # tie-oriented modeling

        # (1) waiting times vs. theoretical distribution
        if(!attr(reh,"ordinal")){
            sum_rates <- lapply(x$diagnostics$rates,sum)
            observed <- sort(reh$intereventTime*unlist(sum_rates))
            theoretical <- stats::qexp(p = c(1:reh$M)/reh$M, rate = 1)
            plot(theoretical,observed, xlab = "Theoretical Quantiles", 
            ylab = "Observed Quantiles") # use bquote() / Observed quatiles : ~(t[m]-t[m-1])*Sigma[e](lambda[e](t[m])) / Theoretical Quantiles ~Exp(1)
            mtext(text = "Q-Q waiting times", side = 3, line = 2,cex=1.5)
            abline(a=0,b=1,lty=2,lwd=1.5)
        }

        # (2) standardized Schoenfeld's residuals
        resids <- residuals(x)
        P <- dim(resids$standardized_residuals)[2] # number of statistics
        #m_exclude <- (reh$M-dim(resids$standardized_residuals)[1])
        for(p in 1:P){
            t_p <- reh$edgelist$time #[-c(1:m_exclude)]
            y_p <- resids$standardized_residuals[,p]
            qrt_p <- quantile(y_p, probs=c(0.25,0.75))
            lb_p <- qrt_p[1] - diff(qrt_p)*3
            ub_p <- qrt_p[2] + diff(qrt_p)*3
            ylim_p <- c(min(y_p),max(y_p))
            if(length(which(y_p<ub_p & y_p>lb_p)) >= floor(0.95*length(y_p))){ # reduce the ylim only if the bound lb_p and ub_p include more than the 90% of the observations
                ylim_p <- c(lb_p,ub_p)
            }
            w_p <- resids$smoothing_weights[,p]
            w_p[w_p<0] <- w_p[w_p>1e04] <- 0.0
            par(mfrow=c(1,1))
            plot(t_p,y_p,xlab = "time", ylab = "scaled Schoenfeld's residuals",ylim=ylim_p) # plotting standardized Schoenfeld's residuals
            lines(smooth.spline(x = t_p, y = y_p,w=w_p,cv=FALSE),lwd=1.5,col=2) # plotting smoothed weighted regression of the residuals
            abline(h=0,col="gray40",lwd=1.5,lty=2)
            mtext(text = colnames(resids$standardized_residuals)[p], side = 3, line = 1,cex=1.5)
        }

        # (3) observed event rates with top % event rate shaded regions
        #observed <- sapply(1:reh$M, function(y) x$diagnostics$rates[[y]][attr(reh,"dyad")[y]]) 
        #thres <- sapply(1:reh$M, function(y) quantile(x$diagnostics$rates[[y]],probs=c(0.75,0.90,0.95))) # M x nPercentiles 
        #discrete_scale_colors <- c("#64f75739", "#e7f75739", "#f7a25739", "#f7575739")
        # plot observed event rate
        #plot(reh$edgelist$time,observed,ylim=c(min(observed),max(observed)),xlab="time",ylab=bquote(lambda(e[m])),cex=0.8,col="#000000ad")
        # plot colored polygons
        # <5%
        #polygon(c(reh$edgelist$time, rev(reh$edgelist$time)), c(thres[3,], rep(max(observed,thres),reh$M)), col = discrete_scale_colors[1],lty=0)
        # 5-10%
        #polygon(c(reh$edgelist$time, rev(reh$edgelist$time)), c(thres[2,], rev(thres[3,])), col = discrete_scale_colors[2],lty=0)
        # 10-25%
        #polygon(c(reh$edgelist$time, rev(reh$edgelist$time)), c(thres[1,], rev(thres[2,])), col = discrete_scale_colors[3],lty=0)
        # >25%
        #polygon(c(reh$edgelist$time, rev(reh$edgelist$time)), c(rep(min(observed,thres),reh$M), rev(thres[1,])), col = discrete_scale_colors[4],lty=0)
        #mtext(text = "Event rates of observed events", side = 3, line = 2,cex=1.5)
        #mtext(text = "(colored regions based on thresholds of rates)", side = 3, line = 1,cex=1)
        # legend
        #legend("bottom", legend=c(">25%","10%-25%","5%-10%","<5%"),col=discrete_scale_colors, lty=1, lwd =rep(20,4), cex=1,horiz=TRUE)
        #dev.flush()

        # (3) top% classification vs. time
        # the claculation of obs_percentile maybe different for 'manual' risk sets
        #obs_percentile <- sapply(1:reh$M, function(y) {
        #ordered_decr <- c(1:reh$D)[order(x$diagnostics$rates[[y]],decreasing=TRUE)]
        #which(ordered_decr==attr(reh,"dyad")[y])/reh$D
        #})
        #classification <- cut(x = obs_percentile, breaks = c(0,0.05,0.10,0.25,1), labels = c("<0.05","0.05-0.1","0.1-0.25",">0.25"), right = TRUE, ordered_result = TRUE)
        #plot(reh$edgelist$time,classification)
        #dev.flush()
        
    }

    else if(attr(x,"model") == "actor"){ # actor-oriented modeling
        # code here for actor-oriented modeling
        which_model <- c("sender_model","receiver_model")
        title_model <- c("Rate model (sender)","Choice model (receiver)")
        resids <- residuals(x)
        for(i in 1:2){
            if(!is.null(x[[which_model[i]]])){
                # (1) waiting times vs. theoretical distribution
                if(!attr(reh,"ordinal") & i==1){
                    sum_rates <- lapply(x[[which_model[i]]]$diagnostics$rates,sum)
                    observed <- sort(reh$intereventTime*unlist(sum_rates))
                    theoretical <- stats::qexp(p = c(1:reh$M)/reh$M, rate = 1)
                    plot(theoretical,observed, xlab = "Theoretical Quantiles", 
                    ylab = "Observed Quantiles") # use bquote() / Observed quatiles : ~(t[m]-t[m-1])*Sigma[e](lambda[e](t[m])) / Theoretical Quantiles ~Exp(1)
                    mtext(text = "Q-Q waiting times", side = 3, line = 2,cex=1.5)
                    mtext(text = title_model[i], side = 3, line = 1,cex=1)
                    abline(a=0,b=1,lty=2,lwd=1.5)
                }

                # (2) standardized Schoenfeld's residuals
                #m_exclude <- (reh$M-dim(resids[[which_model[i]]]$standardized_residuals)[1])
                P <- dim(resids[[which_model[i]]]$standardized_residuals)[2] # number of statistics
                for(p in 1:P){
                    t_p <- reh$edgelist$time #[-c(1:m_exclude)]
                    y_p <- resids[[which_model[i]]]$standardized_residuals[,p]
                    qrt_p <- quantile(y_p, probs=c(0.25,0.75))
                    lb_p <- qrt_p[1] - diff(qrt_p)*3
                    ub_p <- qrt_p[2] + diff(qrt_p)*3
                    ylim_p <- c(min(y_p),max(y_p))
                    if(length(which(y_p<ub_p & y_p>lb_p)) >= floor(0.95*length(y_p))){ # reduce the ylim only if the bound lb_p and ub_p include more than the 90% of the observations
                        ylim_p <- c(lb_p,ub_p)
                    }
                    w_p <- resids[[which_model[i]]]$smoothing_weights[,p]
                    w_p[w_p<0] <- w_p[w_p>1e04] <- 0.0
                    par(mfrow=c(1,1))
                    plot(t_p,y_p,xlab = "time", ylab = "scaled Schoenfeld's residuals",ylim=ylim_p) # plotting standardized Schoenfeld's residuals
                    lines(smooth.spline(x = t_p, y = y_p,w=w_p,cv=FALSE),lwd=1.5,col=2) # plotting smoothed weighted regression of the residuals
                    abline(h=0,col="gray40",lwd=1.5,lty=2)
                    mtext(text = colnames(resids[[which_model[i]]]$standardized_residuals)[p], side = 3, line = 2,cex=1.5)
                    mtext(text = title_model[i], side = 3, line = 1, cex = 1)
                }

                # (3) observed event rates with top % event rate shaded regions
                which_actor <- NULL
                if(i==1){ # sender model
                    which_actor <- reh$edgelist$actor1_ID 
                }
                else{ # receiver model
                    which_actor <- unlist(sapply(1:reh$M, function(y) {
                        s <- rep(0,reh$N)
                        s[reh$edgelist$actor2_ID[y]] <- 1
                        which(s[-reh$edgelist$actor1_ID[y]] == 1)
                        }
                    ))
                }
                #observed <- sapply(1:reh$M, function(y) x[[which_model[i]]]$diagnostics$rates[[y]][which_actor[y]]) 
                #thres <- sapply(1:reh$M, function(y) quantile(x[[which_model[i]]]$diagnostics$rates[[y]],probs=c(0.75,0.90,0.95))) # M x nPercentiles 
                #discrete_scale_colors <- c("#64f75739", "#e7f75739", "#f7a25739", "#f7575739")
                # plot observed event rate
                #plot(reh$edgelist$time,observed,ylim=c(min(observed),max(observed)),xlab="time",ylab="event rate",cex=0.8,col="#000000ad")
                # plot colored polygons
                # <5%
                #polygon(c(reh$edgelist$time, rev(reh$edgelist$time)), c(thres[3,], rep(max(observed,thres),reh$M)), col = discrete_scale_colors[1],lty=0)
                # 5-10%
                #polygon(c(reh$edgelist$time, rev(reh$edgelist$time)), c(thres[2,], rev(thres[3,])), col = discrete_scale_colors[2],lty=0)
                # 10-25%
                #polygon(c(reh$edgelist$time, rev(reh$edgelist$time)), c(thres[1,], rev(thres[2,])), col = discrete_scale_colors[3],lty=0)
                # >25%
                #polygon(c(reh$edgelist$time, rev(reh$edgelist$time)), c(rep(min(observed,thres),reh$M), rev(thres[1,])), col = discrete_scale_colors[4],lty=0)
                #mtext(text = "Observed event rate", side = 3, line = 2,cex=1.5)
                #mtext(text = "(colored regions based on thresholds of rates)", side = 3, line = 1,cex=1)
                # legend 
                #legend("bottom", legend=c(">25%","10%-25%","5%-10%","<5%"),col=discrete_scale_colors, lty=1, lwd =rep(20,4), cex=1,horiz=TRUE)
                #dev.flush()

                # (3) top% classification vs. time
                # the claculation of obs_percentile maybe different for 'manual' risk sets
                #obs_percentile <- sapply(1:reh$M, function(y) {
                #ordered_decr <- c(1:reh$D)[order(x[[which_model[i]]]$diagnostics$rates[[y]],decreasing=TRUE)]
                #which(ordered_decr==which_actor[y])/reh$D
                #})
                #classification <- cut(x = obs_percentile, breaks = c(0,0.05,0.10,0.25,1), labels = c("<0.05","0.05-0.1","0.1-0.25",">0.25"), right = TRUE, ordered_result = TRUE)
                #plot(reh$edgelist$time,classification)
                #dev.flush()
            }
            
        }

    }

}


#######################################################################################
#######################################################################################

# aic
#' @title aic
#' @description A function that returns the AIC (Akaike's Information Criterion) value in a 'remstimate' object
#' @param object is a \code{remstimate} object 
#' @param ... further arguments to be passed to the 'aic' method
#' @export
#' 
#' @examples
#' 
#' # ------------------------------------ #
#' #       tie-oriented model: "MLE"      #
#' # ------------------------------------ #
#' 
#' # loading data
#' data(tie_reh)
#'   
#' # specifying linear predictor
#' tie_model <- ~ 1 + 
#'                remstats::indegreeSender()+
#'                remstats::inertia()+
#'                remstats::reciprocity() 
#' 
#' # calculating statistics
#' tie_reh_stats <- remstats::remstats(reh = tie_reh, 
#'                                     tie_effects = tie_model)
#' 
#' # running estimation
#' tie_mle <- remstimate::remstimate(reh = tie_reh,
#'                                   stats = tie_reh_stats,
#'                                   method = "MLE",
#'                                   ncores = 1)
#' 
#' # AIC
#' aic(tie_mle) #
#' 
aic <- function(object,...){
  UseMethod("aic")
}


#' @describeIn aic AIC (Akaike's Information Criterion) value of a 'remstimate' object
#' @method aic remstimate
#' @export
aic.remstimate <- function(object,...) {
    if(attr(object, "approach") == "Frequentist"){
        if(attr(object,"model") == "tie"){
            return(object$AIC)
        }
        else if(attr(object,"model") == "actor"){
            AIC <- NULL
            if(!is.null(object$sender_model)){
                AIC <- c("sender model" = object$sender_model$AIC)
            }
            if(!is.null(object$receiver_model)){
                AIC <- c(AIC, "receiver model" = object$receiver_model$AIC)
            }
            return(AIC)
        }
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
#' 
#' @examples
#' 
#' # ------------------------------------ #
#' #       tie-oriented model: "MLE"      #
#' # ------------------------------------ #
#' 
#' # loading data
#' data(tie_reh)
#'   
#' # specifying linear predictor
#' tie_model <- ~ 1 + 
#'                remstats::indegreeSender()+
#'                remstats::inertia()+
#'                remstats::reciprocity() 
#' 
#' # calculating statistics
#' tie_reh_stats <- remstats::remstats(reh = tie_reh, 
#'                                     tie_effects = tie_model)
#' 
#' # running estimation
#' tie_mle <- remstimate::remstimate(reh = tie_reh,
#'                                   stats = tie_reh_stats,
#'                                   method = "MLE",
#'                                   ncores = 1)
#' 
#' # AICC
#' aicc(tie_mle) 
#' 
aicc <- function(object,...){
  UseMethod("aicc")
}


#' @describeIn aicc AICC (Akaike's Information Corrected Criterion) value of a 'remstimate' object
#' @method aicc remstimate
#' @export
aicc.remstimate <- function(object,...) {
    if(attr(object, "approach") == "Frequentist"){
        if(attr(object,"model") == "tie"){
            return(object$AICC)
        }
        else if(attr(object,"model") == "actor"){
            AICC <- NULL
            if(!is.null(object$sender_model)){
                AICC <- c("sender model" = object$sender_model$AICC)
            }
            if(!is.null(object$receiver_model)){
                AICC <- c(AICC, "receiver model" = object$receiver_model$AICC)
            }
            return(AICC)
        }
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
#' 
#' @examples
#' 
#' # ------------------------------------ #
#' #       tie-oriented model: "MLE"      #
#' # ------------------------------------ #
#' 
#' # loading data
#' data(tie_reh)
#'   
#' # specifying linear predictor
#' tie_model <- ~ 1 + 
#'                remstats::indegreeSender()+
#'                remstats::inertia()+
#'                remstats::reciprocity() 
#' 
#' # calculating statistics
#' tie_reh_stats <- remstats::remstats(reh = tie_reh, 
#'                                     tie_effects = tie_model)
#' 
#' # running estimation
#' tie_mle <- remstimate::remstimate(reh = tie_reh,
#'                                   stats = tie_reh_stats,
#'                                   method = "MLE",
#'                                   ncores = 1)
#' 
#' # BIC
#' bic(tie_mle) 
#' 
bic <- function(object,...){
  UseMethod("bic")
}


#' @describeIn bic BIC (Bayesian Information Criterion) value of a 'remstimate' object
#' @method bic remstimate
#' @export
bic.remstimate <- function(object,...) {
    if(attr(object, "approach") == "Frequentist"){
        if(attr(object,"model") == "tie"){
            return(object$BIC)
        }
        else if(attr(object,"model") == "actor"){
            BIC <- NULL
            if(!is.null(object$sender_model)){
                BIC <- c("sender model" = object$sender_model$BIC)
            }
            if(!is.null(object$receiver_model)){
                BIC <- c(BIC, "receiver model" = object$receiver_model$BIC)
            }
            return(BIC)
        }
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
#' 
#' @examples 
#' 
#' # No examples available at the moment
#' 
waic <- function(object,...){
  UseMethod("waic")
}


#' @describeIn waic WAIC (Watanabe-Akaike's Information Criterion) value of a 'remstimate' object
#' @method waic remstimate
#' @export
waic.remstimate <- function(object,...) {
    if(attr(object, "approach") == "Frequentist"){
        if(attr(object,"model") == "tie"){
            return(object$WAIC)
        }
        else if(attr(object,"model") == "actor"){
            WAIC <- NULL
            if(!is.null(object$sender_model)){
                WAIC <- c("sender model" = object$sender_model$WAIC)
            }
            if(!is.null(object$receiver_model)){
                WAIC <- c(WAIC, "receiver model" = object$receiver_model$WAIC)
            }
            return(WAIC)
        }
    }
    else{stop("'approach' must be 'Frequentist'")}   
}


#######################################################################################
#######################################################################################




