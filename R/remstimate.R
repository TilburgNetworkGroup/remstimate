#' remstimate  
#'
#' A function for the optimization of tie-oriented and actor-oriented likelihood. There are four optimization algorithms: two Frequentists, Maximum Likelihood Estimation (\code{MLE}) and Adaptive Gradient Descent (\code{GDADAMAX}), and two Bayesian, Bayesian Sampling Importance Resampling (\code{BSIR}) and Hamiltonian Monte Carlo (\code{HMC}).
#'
#' @param reh a \code{remify} object of the processed relational event history. Output object of the function \code{remify::remify()}.
#' @param stats a \code{remstats} object: when `attr(reh,"model")` is `"tie"`, \code{stats} is an array of statistics with dimensions \code{[M x D x P]}: where \code{M} is the number of events, \code{D} is the number of possible dyads (full riskset), \code{P} is the number of statistics; if `attr(reh,"model")` is `"actor"`, \code{stats} is a list that can contain up to two arrays named \code{"sender_stats"} and \code{"receiver_stats"} with dimensions \code{[M x N x P]}, where \code{N} are the actors (senders in the array \code{"sender_stats"}, receivers in the array \code{"receiver_stats"}). Furthermore, it is possible to only estimate the sender rate model or only the receiver choice model, by using the correct naming of the arrays.
#' @param method the optimization method to estimate model parameters. Methods available are: Maximum Likelihood Estimation (\code{"MLE"}, and also the default method), Adaptive Gradient Descent (\code{"GDADAMAX"}), Bayesian Sampling Importance Resampling (\code{"BSIR"}), Hamiltonian Monte Carlo (\code{"HMC"}). (default method is \code{"MLE"}).
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
#' data(tie_data)
#' 
#' # processing event sequence with remify
#' tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
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
#' data(ao_data)
#' 
#' # processing event sequence with remify
#' ao_reh <- remify::remify(edgelist = ao_data$edgelist, model = "actor")
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
                       ncores = attr(reh,"ncores"),
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
    if(!is.null(attr(stats,"model"))){
        if(attr(stats,"model")!=attr(reh,"model")){
            stop("attribute 'model' of input 'reh' and input 'stats' must be the same")
        }
    }

    # ... active riskset? then overwrite two objects (this prevents from coding many ifelse() to switch between active-oriented riskset objects and full-oriented riskset objects)
    if((attr(reh,"riskset") == "active")){
        reh$D <- reh$activeD
        if(model == "tie"){
            attr(reh,"dyadID") <- attr(reh,"dyadIDactive")
            reh$omit_dyad <- list() # because "reh$omit_dyad$time" and "reh$omit_dyad$riskset" for riskset="active" are obsolete (will be removed from remify output in the future 3.x.x version)
            # check line 113 remify.cpp
        }
    }

    # ... method
    if(is.null(method)){
        method <- "MLE" # method "MLE" as default
    }
    if(!any(method %in% c("MLE","GDADAMAX","BSIR","HMC"))){stop("The `method` specified is not available or it is mistyped.")}
    else{
        method  <- match.arg(arg = method, choices = c("MLE","GDADAMAX","BSIR","HMC"), several.ok = FALSE) # default is "MLE"
    }

    # ... type of likelihood
    ordinal <- attr(reh,"ordinal")

    # ... directed / undirected network
    if((model == "actor") & !attr(reh,"directed")){stop("actor-oriented modeling can't operate on undirected networks")}

    # ... ncores
    if(is.null(ncores)) ncores <- attr(reh,"ncores")
    if(((parallel::detectCores() > 2L) & (ncores > floor(parallel::detectCores()-2L))) | ((parallel::detectCores() == 2L) & (ncores > 1L))){
            stop("'ncores' is recommended to be set at most to: floor(parallel::detectCores()-2L)")
    }

    # ... omit dyad
    if(is.null(reh$omit_dyad)){
        reh$omit_dyad <- list()
    }

    # ... processing start and stop values and method from ("subset" and "method" attribute of remstats object)
    # we do this now because later the stats object will change dimensions and won't be a remstats object anymore
    if(all(inherits(stats,c("remstats","tomstats"),TRUE)) | all(inherits(stats,c("remstats","aomstats"),TRUE))){
        stats_attr_method <- attr(stats,"method")
        if(is.null(stats_attr_method)){
            stop("attribute 'method' not found inside object 'remstats'. Input argument 'stats' must be an object of class 'remstats' from the package 'remstats' (>=3.2.0)")
        }
        omit_dyad_receiver <- NULL
        if(stats_attr_method == "pe"){
            if(!is.null(attr(reh,"evenly_spaced_interevent_time"))){
                reh$intereventTime <- attr(reh,"evenly_spaced_interevent_time")
            }
            if(!is.null(reh$E)){ # reh$E is NULL only when there are no simultaneous events
                reh$M <- reh$E # overwriting dimension (we can do it because remstimate works only with reh$M so if the method is "pt", reh$M will remain so. For method "pe" we assign reh$E to reh$M
            }
        }
        if(!is.null(attr(stats,"subset"))){
            start_stop <- as.numeric(unlist(attr(stats,"subset")))
        }
        else{
            start_stop <- c(1,reh$M)
        }
    }
    else{
        stop("'stats' must be a 'remstats' object from the package 'remstats' (>= 3.2.0), suitable for tie-oriented modeling ('tomstats') or actor-oriented modeling ('aomstats')")
    }

    # ... stats 
    model_formula <- variables_names <- where_is_baseline <- NULL
    if(model == "tie")
    {   
        if(all(inherits(stats,c("remstats","tomstats"),TRUE))){
            if(!is.null(dimnames(stats)[[3]])){
                variables_names <- dimnames(stats)[[3]]
            }
            if(is.null(attr(stats,"formula"))){
                model_formula <- stats::as.formula(paste("~ ",paste(variables_names,collapse=" + ")))
            }
            else{
                model_formula <- attr(stats,"formula")
            }
            # is there a baseline term?
            if(any(tolower(variables_names) %in% c("baseline"))){
                where_is_baseline <- which(variables_names == "baseline")
            }
            stats <- aperm(stats, perm = c(2,3,1)) # stats reshaped in [D*U*M]
            # start and stop for tie-oriented model
            if(stats_attr_method == "pt"){
                attr(reh,"dyadID") <- attr(reh,"dyadID")[start_stop[1]:start_stop[2]] # this is already working for dyadIDactive because we reassing attribute dyadID in line 123
                attr(reh,"actor1ID") <- attr(reh,"actor1ID")[start_stop[1]:start_stop[2]]
                attr(reh,"actor2ID") <- attr(reh,"actor2ID")[start_stop[1]:start_stop[2]] 
                if((length(reh$omit_dyad)>0) & !is.null(attr(reh,"indices_simultaneous_events"))){
                    reh$omit_dyad$time <-  reh$omit_dyad$time[-attr(reh,"indices_simultaneous_events")] 
                }
            }
            else if(stats_attr_method =="pe"){
                attr(reh,"dyadID") <- unlist(attr(reh,"dyadID"))[start_stop[1]:start_stop[2]] # this is already working for dyadIDactive because we reassing attribute dyadID in line 123
                attr(reh,"actor1ID") <- unlist(attr(reh,"actor1ID"))[start_stop[1]:start_stop[2]]
                attr(reh,"actor2ID") <- unlist(attr(reh,"actor2ID"))[start_stop[1]:start_stop[2]] 
            } 
            reh$intereventTime <- reh$intereventTime[start_stop[1]:start_stop[2]] # in line 168 we already re-assigned the intereventTime variable in case of method="pe", so line 224 is a valid processing for "pt" and "pe"
            reh$M <- diff(start_stop)+1
            if(length(reh$omit_dyad)>0){
                reh$omit_dyad$time <- reh$omit_dyad$time[start_stop[1]:start_stop[2]]
            }
        }
        else if(all(inherits(stats,c("remstats","aomstats"),TRUE))){
            stop("'remstats' object supplied cannot work for tie-oriented modeling")
        }

        # .. check on dimensions
        if(length(attr(reh,"dyadID")) != dim(stats)[3]){ # if the dimension of the processed intereventTime are different from the dimensions of the input stats object, then throw error
            stop("the number of time points (or number of events) doesn't match the (row) dimension of the 'remstats' object")
        }
    }
    if(model == "actor") 
    {
        model_formula <- list() # becomes a list
        if(all(inherits(stats,c("remstats","aomstats"),TRUE))){
            variables_rate <- variables_choice <- NULL
            if(!is.null(stats$sender_stats)){ # sender model is specified
                variables_rate <- dimnames(stats$sender_stats)[[3]]
                if(!is.null(attr(stats,"formula")$rate)){
                    model_formula[["rate_model_formula"]] <- attr(stats,"formula")$rate
                }
                else{
                    model_formula[["rate_model_formula"]] <- stats::as.formula(paste("~ ",paste(variables_rate,collapse=" + ")))
                }
                # is there a baseline term?
                if(any(tolower(variables_rate) %in% c("baseline"))){
                    where_is_baseline <- which(variables_rate == "baseline")
                }
                stats$sender_stats <- aperm(stats$sender_stats, perm = c(2,3,1)) # stats reshaped in [N*U*M]
            }
            if(!is.null(stats$receiver_stats)){ # receiver model is specified
                variables_choice <- dimnames(stats$receiver_stats)[[3]] 
                if(!is.null(attr(stats,"formula")$choice)){
                    model_formula[["choice_model_formula"]] <- attr(stats,"formula")$choice
                }
                else{
                    model_formula[["choice_model_formula"]] <- stats::as.formula(paste("~ ",paste(variables_choice,collapse=" + ")))
                }
                stats$receiver_stats <- aperm(stats$receiver_stats, perm = c(2,3,1)) # stats reshaped in [N*U*M]
            }

            # vector of variable names and list of model formulas
            variables_names <- list(sender_model = variables_rate, receiver_model = variables_choice)

            # start and stop for actor-oriented model
            if(stats_attr_method == "pt"){
                attr(reh,"actor1ID") <- attr(reh,"actor1ID")[start_stop[1]:start_stop[2]]
                attr(reh,"actor2ID") <- unlist(attr(reh,"actor2ID")[start_stop[1]:start_stop[2]]) # unlist here because of receiver-choice model
                if(is.null(reh$E)){ # this scenario can still happen
                    reh$E <- reh$M
                }
                if(length(reh$omit_dyad)>0){

                    # for the receiver model 
                    if(!is.null(stats$receiver_stats)){
                        start_stop_time <- unique(reh$edgelist$time)[start_stop]
                        lb_time <- min(which(reh$edgelist$time>=start_stop_time[1]))
                        ub_time <- max(which(reh$edgelist$time<=start_stop_time[2]))
                        omit_dyad_receiver <- list(time = reh$omit_dyad$time[lb_time:ub_time], riskset = reh$omit_dyad$riskset) # for the receiver model
                    }

                    # for the sender model (we process now the sender model because this will modify the reh$omit_dyad$time object)
                    if(!is.null(stats$sender_stats)){
                        if(!is.null(attr(reh,"indices_simultaneous_events"))){
                            reh$omit_dyad$time <-  reh$omit_dyad$time[-attr(reh,"indices_simultaneous_events")][start_stop[1]:start_stop[2]] # for the sender model
                        }
                        else{
                            reh$omit_dyad$time <-  reh$omit_dyad$time[start_stop[1]:start_stop[2]] # for the sender model
                        }
                    }
                }
            }
            else if(stats_attr_method == "pe"){
                attr(reh,"actor1ID") <- unlist(attr(reh,"actor1ID"))[start_stop[1]:start_stop[2]]
                attr(reh,"actor2ID") <- unlist(attr(reh,"actor2ID"))[start_stop[1]:start_stop[2]]
                if(length(reh$omit_dyad)>0){
                    if(!is.null(stats$receiver_stats)){
                        omit_dyad_receiver <- list(time = reh$omit_dyad$time[start_stop[1]:start_stop[2]]  , riskset = reh$omit_dyad$riskset)
                    } 
                    reh$omit_dyad$time <- reh$omit_dyad$time[start_stop[1]:start_stop[2]]        
                }
            } 
            if(!ordinal){
                reh$intereventTime <- reh$intereventTime[start_stop[1]:start_stop[2]]
            }
            reh$M <- diff(start_stop)+1

            # .. check on dimensions
            no_correct_dimensions <- FALSE
            if(!is.null(stats$sender_stats)){
                if((length(attr(reh,"actor1ID")) != dim(stats$sender_stats)[3])){ # if the dimension of the processed intereventTime is different from the dimensions of the input stats object, then throw error
                    no_correct_dimensions <- TRUE
                }
            }
            if(!is.null(stats$receiver_stats)){
                if((length(attr(reh,"actor2ID")) != dim(stats$receiver_stats)[3])){ # if the dimension of the edgelist is different from the dimensions of the input stats object, then throw error
                    no_correct_dimensions <- TRUE
                }
            }
            if(no_correct_dimensions){
                stop("the number of time points (or number of events) doesn't match the (row) dimension of the 'remstats' object")
            }

        }
        else if(all(inherits(stats,c("remstats","tomstats"),TRUE))){
            stop("'remstats' object supplied cannot work for actor-oriented modeling")
        }
    }

    # ... adjusting the intereventTime
    if(ordinal){
        reh$intereventTime <- c(1) # we can assign a vector of length 1 because the intereventTime will not be used from any ordinal likelihood
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
    
    # ... setting up objects for actor-oriented model
    if(model == "actor"){
        senderRate <- c(TRUE,FALSE) # meaning that the first BSIR will be run on the "sender rate" model and the second on the "receiver choice"
        which_model <- c("sender_model","receiver_model")
        which_stats <- c("sender_stats","receiver_stats")
        actor1ID_condition <- c(TRUE,FALSE)
        if(stats_attr_method == "pe"){
            actor1ID_condition <- as.logical(actor1ID_condition*FALSE) # actor1ID is needed unlist both for sender and receiver model
        }
    }

    # ... estimating with Maximum Likelihood (MLE) or Gradient Descent ADAMAX (GDADAMAX)
    if(method %in% c("MLE","GDADAMAX")){
        if(model == "tie"){ # Tie-oriented Model  (Relational Event Model)
            if(method == "MLE"){
                optimum_obj <- trust::trust(objfun = remDerivatives, 
                                            parinit = rep(0,dim(stats)[2]), 
                                            rinit = 1, 
                                            rmax = 100, 
                                            stats = stats,  
                                            actor1 = list(),
                                            actor2 = list(),
                                            dyad = attr(reh,"dyadID"),
                                            omit_dyad = reh$omit_dyad,
                                            interevent_time = reh$intereventTime,
                                            model = model,
                                            ordinal = ordinal,
                                            ncores = ncores)
            }  
            if(method == "GDADAMAX"){
                if(is.null(init)){
                    init <- stats::runif(dim(stats)[2],-0.1,0.1) # previously rep(0,dim(stats)[2])
                    if(any(variables_names == "baseline")){
                        init[which(variables_names == "baseline")] <- init[which(variables_names == "baseline")] + (log(reh$M) - log(sum(reh$intereventTime)*reh$D))
                    }
                }
                optimum_obj <- GDADAMAX(pars = init,
                                    stats = stats,  
                                    actor1 = list(),
                                    actor2 = list(),
                                    dyad = attr(reh,"dyadID"),
                                    omit_dyad = reh$omit_dyad,
                                    interevent_time = reh$intereventTime,
                                    model = model,
                                    ordinal = ordinal,
                                    ncores = ncores,
                                    epochs = epochs,
                                    epsilon = epsilon)
                optimum_obj$hessian <- remDerivativesStandard(pars = optimum_obj$argument,
                                                                    stats = stats,
                                                                    dyad = attr(reh,"dyadID"),
                                                                    omit_dyad = reh$omit_dyad,
                                                                    interevent_time = reh$intereventTime,
                                                                    ordinal = ordinal,
                                                                    ncores = ncores)
                optimum_obj$hessian <- optimum_obj$hessian$hessian 
                                                                
            }            
            # coefficients                         
            remstimateList$coefficients <- as.vector(optimum_obj$argument)
            names(remstimateList$coefficients) <- variables_names
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
            names(remstimateList$coefficients) <- names(remstimateList$se) <- rownames(remstimateList$vcov) <- colnames(remstimateList$vcov) <- variables_names                      
            # residual.deviance
            remstimateList$residual.deviance <- -2*remstimateList$loglik        
            # null.deviance
            remstimateList$null.deviance <- 2*(trust::trust(objfun = remDerivatives, 
                                        parinit = c(0), 
                                        rinit = 1, 
                                        rmax = 100, 
                                        stats = array(1,dim=c(reh$D,1,reh$M)),  
                                        actor1 = list(),
                                        actor2 = list(),
                                        dyad = attr(reh,"dyadID"),
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
        if(model == "actor"){ # Actor-oriented Model 
            optimum_model <- list()
            for(i in 1:2){
                if(!is.null(stats[[which_stats[i]]])){  # run estimation only if the model is specified
                    actor1ID_ls <- if(actor1ID_condition[i]) attr(reh,"actor1ID") else unlist(attr(reh,"actor1ID"))
                    omit_dyad_actor <- if(senderRate[i]) reh$omit_dyad else omit_dyad_receiver
                    if(method == "MLE"){ # Maximum Likelihood Estimates
                        optimum_model[[which_model[i]]] <- trust::trust(objfun = remDerivatives, 
                                    parinit = rep(0,dim(stats[[which_stats[i]]])[2]), 
                                    rinit = 1, 
                                    rmax = 100, 
                                    stats = stats[[which_stats[i]]], 
                                    actor1 = actor1ID_ls,
                                    actor2 = attr(reh,"actor2ID"),
                                    dyad = list(),
                                    omit_dyad = omit_dyad_actor, 
                                    interevent_time = reh$intereventTime,
                                    model = model,
                                    ordinal = ordinal,
                                    senderRate = senderRate[i],
                                    N = reh$N,
                                    ncores = ncores)
                    }
                    else if(method == "GDADAMAX"){ # Gradient Descent
                        if(is.null(init[[which_model[i]]])){
                            init[[which_model[i]]] <- stats::runif(dim(stats[[which_stats[i]]])[2],-0.1,0.1) # rep(0.0,dim(stats[[which_stats[i]]])[2]) 
                        }
                        optimum_model[[which_model[i]]] <- GDADAMAX(pars = init[[which_model[i]]], 
                                                stats = stats[[which_stats[i]]],  
                                                actor1 = actor1ID_ls,
                                                actor2 = attr(reh,"actor2ID"),
                                                dyad = list(),
                                                omit_dyad = omit_dyad_actor,
                                                interevent_time = reh$intereventTime,
                                                model = model,
                                                ordinal = ordinal,
                                                senderRate = senderRate[i],
                                                N = reh$N,
                                                ncores = ncores,
                                                epochs = epochs,
                                                epsilon = epsilon)
                        optimum_model[[which_model[i]]]$hessian <- remDerivatives(pars = optimum_model[[which_model[i]]]$argument,
                                                                        stats = stats[[which_stats[i]]],  
                                                                        actor1 = actor1ID_ls,
                                                                        actor2 = attr(reh,"actor2ID"),
                                                                        dyad = list(),
                                                                        omit_dyad = omit_dyad_actor,
                                                                        interevent_time = reh$intereventTime,
                                                                        model = model,
                                                                        ordinal = ordinal,
                                                                        senderRate = senderRate[i],
                                                                        N = reh$N,
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
                    # residual deviance
                    remstimateList[[which_model[i]]]$residual.deviance <- -2*remstimateList[[which_model[i]]]$loglik
                    # null deviance                                    
                    remstimateList[[which_model[i]]]$null.deviance <- NULL
                    
                    if(senderRate[i]){ # only for sender model
                        remstimateList$sender_model$null.deviance <- 2*(trust::trust(objfun = remDerivatives, 
                                                                                    parinit = c(0), 
                                                                                    rinit = 1, 
                                                                                    rmax = 100, 
                                                                                    stats = array(1,dim=c(dim(stats$sender_stats)[1],1,dim(stats$sender_stats)[3])),  
                                                                                    actor1 = actor1ID_ls,
                                                                                    actor2 = list(),
                                                                                    dyad = list(),
                                                                                    omit_dyad = omit_dyad_actor,
                                                                                    interevent_time = reh$intereventTime,
                                                                                    model = model,
                                                                                    ordinal = ordinal,
                                                                                    ncores = ncores,
                                                                                    senderRate = TRUE)$value)
                    }
                    else{ # for receiver model
                        remstimateList$receiver_model$null.deviance <- ifelse(length(reh$omit_dyad)>0,2*(trust::trust(objfun = remDerivatives, 
                                                                                                                    parinit = c(0), 
                                                                                                                    rinit = 1, 
                                                                                                                    rmax = 100, 
                                                                                                                    stats = array(1,dim=c(dim(stats$receiver_stats)[1],1,dim(stats$receiver_stats)[3])), 
                                                                                                                    actor1 = actor1ID_ls,
                                                                                                                    actor2 = attr(reh,"actor2ID"),
                                                                                                                    dyad = list(),
                                                                                                                    omit_dyad = omit_dyad_actor, 
                                                                                                                    interevent_time = reh$intereventTime,
                                                                                                                    model = model,
                                                                                                                    ncores = ncores,
                                                                                                                    senderRate = FALSE,
                                                                                                                    N = reh$N)$value),-2*log(1/(reh$N-1)))
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
                                        actor1 = list(),
                                        actor2 = list(),
                                        dyad = attr(reh,"dyadID"),
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
                                    dyad = attr(reh,"dyadID"),
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
            names(bsir$coefficients) <- names(bsir$post.mean) <- rownames(bsir$vcov) <- colnames(bsir$vcov) <- names(bsir$sd) <- variables_names
            remstimateList <- bsir
        }
        if(model == "actor"){ # Actor-oriented Model
            bsir <- list()
            for(i in 1:2){
                if(!is.null(stats[[which_stats[i]]])){  # run BSIR only if the model is specified
                    bsir_i <- list()
                    bsir_i$log_posterior <- rep(0,nsim)
                    bsir_i$draws <- list()  
                    actor1ID_ls <- if(actor1ID_condition[i]) attr(reh,"actor1ID") else unlist(attr(reh,"actor1ID"))
                    omit_dyad_actor <- if(senderRate[i]) reh$omit_dyad else omit_dyad_receiver
                    # (1) generate from the proposal (importance distribution)

                    # First find location and scale parameters (mles and variances and covariances matrix)
                    mle_optimum <- trust::trust(objfun = remDerivatives, 
                                                        parinit = rep(0,dim(stats[[which_stats[i]]])[2]), 
                                                        rinit = 1, 
                                                        rmax = 100, 
                                                        stats = stats[[which_stats[i]]], 
                                                        actor1 = actor1ID_ls,
                                                        actor2 = attr(reh,"actor2ID"),
                                                        dyad = list(),
                                                        omit_dyad = omit_dyad_actor,
                                                        interevent_time = reh$intereventTime,
                                                        model = model,
                                                        ordinal = ordinal,
                                                        ncores = ncores,
                                                        senderRate = senderRate[i],
                                                        N = reh$N)                          
                    
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
                                    actor1 = actor1ID_ls,
                                    actor2 = attr(reh,"actor2ID"),
                                    dyad = list(),
                                    omit_dyad = omit_dyad_actor,
                                    interevent_time = reh$intereventTime,
                                    model = model,
                                    ordinal = ordinal,
                                    senderRate = senderRate[i],
                                    ncores = ncores,
                                    gradient = FALSE,
                                    hessian = FALSE,
                                    N = reh$N)$value) # this is the most time consuming step to be optimized, avoiding the apply        
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
                init <- matrix(stats::runif((dim(stats)[2])*nchains,-0.1,0.1),nrow=dim(stats)[2],ncol=nchains)
                if(any(variables_names == "baseline")){
                    init[which(variables_names == "baseline"),] <- init[which(variables_names == "baseline"),] + (log(reh$M) - log(sum(reh$intereventTime)*reh$D))
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
                        actor1 = list(),
                        actor2 = list(),
                        dyad = attr(reh,"dyadID"), 
                        omit_dyad = reh$omit_dyad,
                        interevent_time = reh$intereventTime,
                        model = model,
                        ordinal = ordinal,
                        ncores = ncores,
                        thin = thin,
                        L = L,
                        epsilon = epsilon, 
                        N = reh$N)          
            hmc$coefficients <- hmc_out$draws[which.min(hmc_out$log_posterior),] # log_posterior is a vector of posterior nloglik, therefore we take the minimum
            hmc$post.mean <- colMeans(hmc_out$draws)
            hmc$vcov <- stats::cov(hmc_out$draws)
            hmc$sd <- diag(hmc$vcov)**0.5
            hmc$loglik <- -min(hmc_out$log_posterior) # max loglik is -min(log_posterior)
            hmc$df.null <- reh$M
            hmc$df.model <- dim(stats)[2]
            hmc$df.residual <- hmc$df.null - hmc$df.model
            names(hmc$coefficients) <- names(hmc$post.mean) <- rownames(hmc$vcov) <- colnames(hmc$vcov) <- names(hmc$sd) <- variables_names

            remstimateList <- c(hmc_out,hmc)
        }
        if(model == "actor"){ # Actor-oriented Model
            for(i in 1:2){
                if(!is.null(stats[[which_stats[i]]])){
                    actor1ID_ls <- if(actor1ID_condition[i]) attr(reh,"actor1ID") else unlist(attr(reh,"actor1ID"))
                    omit_dyad_actor <- if(senderRate[i]) reh$omit_dyad else omit_dyad_receiver
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
                                actor1 = actor1ID_ls,
                                actor2 = attr(reh,"actor2ID"),
                                dyad = list(), 
                                omit_dyad = omit_dyad_actor,
                                interevent_time = reh$intereventTime,
                                model = model,
                                senderRate = senderRate[i],
                                N = reh$N,
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
        attr(str_out, "statistics") <- variables_names
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
        attr(str_out, "statistics") <- variables_names
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
    # additional attributes useful for remstimate methods
    attr(str_out,"where_is_baseline") <- where_is_baseline
    attr(str_out,"ncores") <- ncores
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
#' data(tie_data)
#' 
#' # processing event sequence with remify
#' tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
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
#' data(ao_data)
#' 
#' # processing event sequence with remify
#' ao_reh <- remify::remify(edgelist = ao_data$edgelist, model = "actor")
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
#' @description A function that returns the summary of a remstimate object
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
#' data(tie_data)
#' 
#' # processing event sequence with remify
#' tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
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
#' data(ao_data)
#' 
#' # processing event sequence with remify
#' ao_reh <- remify::remify(edgelist = ao_data$edgelist, model = "actor")
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
                        epsilon = attr(object, "epsilon"),
                        chiP = 1 - stats::pchisq(object$model.deviance, object$df.model))                
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
                                            BIC = object$sender_model$BIC,
                                            chiP = 1 - stats::pchisq(object$sender_model$model.deviance, object$sender_model$df.model)
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
                                            BIC = object$receiver_model$BIC,
                                            chiP = 1 - stats::pchisq(object$receiver_model$model.deviance, object$receiver_model$df.model)
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
    else{
        if(length(summary_out)==0) stop("invalid 'remstimate' object")
    }
    class(summary_out) <- "summary.remstimate"

    cat("Relational Event Model",paste("(",summary_out$model," oriented)",sep=""),"\n\n")
    if(summary_out$model == "tie"){
        cat("Call:\n",deparse(summary_out$formula),"\n\n",sep="")
        second_line <- paste("(",summary_out$method," with ",sep="")
        if(summary_out$ordinal) second_line <- paste(second_line,"ordinal likelihood):\n\n",sep="")
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

    invisible(summary_out)
}



# diagnostics
#' @title diagnostics
#' @description A function that returns the diagnostics of a \code{remstimate} object. The output object of the method \code{diagnostics} contains the residuals of the model estimated in the \code{remstimate} object, and the event rates estimated from the model at each tiem point. For tie-oriented modeling frameworks the object contains: a list \code{residuals} with two objects, \code{standardized_residuals} containing standardized Schoenfeld's residuals (Schoenfeld, D., 1982, <doi:10.2307/2335876>; Grambsch, P. M., & Therneau, T. M., 1994, <doi:10.2307/2337123>; Winnett, A., & Sasieni, P., 2001, <jstor.org/stable/2673500>), and \code{smoothing_weights} (a matrix of weights used for the red smooth splines in the plot of the residuals), an array structure \code{rates} with the event rates estimated from the optimized model parameters, and \code{.reh.processed} which is a pseudo-hidden object containing a further processed \code{remify} object that helps speed up the plotting function \code{plot.remstimate} and that the user is not supposed to modify. As to the actor-oriented modeling frameworks, in the diagnostics output there are two main list objects named after \code{sender_model} and \code{receiver_model}. After selecting the model, the structure of diagnostics is the same as for the tie-oriented model. Each model's diagnostics (sender or receiver) is available only if the corresponding model is found in the \code{remstimate} object.
#' @param object is a \code{remstimate} object 
#' @param reh is a \code{remify} object, the same used for the 'remstimate' object
#' @param stats is a \code{remstats} object, the same used for the 'remstimate' object
#' @param ... further arguments to be passed to the 'diagnostics' method
#' @export
#' 
#' @return A object of class \code{"remstimate","diagnostics"} with standardized Schoenfeld's residuals and estimated event rates given the optimized model parameters. 
#'
#' 
#' @examples
#' 
#' # ------------------------------------ #
#' #       tie-oriented model: "MLE"      #
#' # ------------------------------------ #
#' 
#' # loading data
#' data(tie_data)
#' 
#' # processing event sequence with remify
#' tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
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
#' # diagnostics
#' tie_diagnostics <- diagnostics(object = tie_mle, reh = tie_reh, stats = tie_reh_stats)
#' names(tie_diagnostics)
#' 
diagnostics <- function(object,reh,stats,...){
  UseMethod("diagnostics")
}


#' @describeIn diagnostics diagnostics of a 'remstimate' object
#' @method diagnostics remstimate
#' @export
diagnostics.remstimate <- function(object,reh,stats,...) {

    # processing dyadID actor1ID and actor2ID depending on pt/pe and on start and stop values

    # ... processing input arguments

    # ... remify object ('reh' input argument)
    if(!inherits(reh,"remify")){
        stop("'reh' must be a 'remify' object (see ?remify::remify).")
    }

    # ... model
    model <- attr(reh, "model") # attribute from reh object
    if(is.null(model) | !(model %in% c("actor","tie"))){
        stop("attribute 'model' of input 'reh' must be either 'actor' or 'tie'")
    }
    if(!is.null(attr(stats,"model"))){
        if(attr(stats,"model")!=attr(reh,"model")){
            stop("attribute 'model' of input 'reh' and input 'stats' must be the same")
        }
    }

    # ... active riskset? then overwrite two objects (this prevents from coding many ifelse() to switch between active-oriented riskset objects and full-oriented riskset objects)
    if((attr(reh,"riskset") == "active")){
        reh$D <- reh$activeD
        if(model == "tie"){
            attr(reh,"dyadID") <- attr(reh,"dyadIDactive")
            reh$omit_dyad <- list() # because "reh$omit_dyad$time" and "reh$omit_dyad$riskset" for riskset="active" are obsolete (will be removed from remify output in the future 3.x.x version)
            # check line 113 remify.cpp
        }
    }

    # ... type of likelihood
    ordinal <- attr(reh,"ordinal")

    # ... ncores
    if(!is.null(attr(remstimate,"ncores"))){
        ncores <- attr(remstimate,"ncores")
    }

    # ... omit dyad
    if(is.null(reh$omit_dyad)){
        reh$omit_dyad <- list()
    }

    # ... processing start and stop values and method from ("subset" and "method" attribute of remstats object)
    # we do this now because later the stats object will change dimensions and won't be a remstats object anymore
    if(all(inherits(stats,c("remstats","tomstats"),TRUE)) | all(inherits(stats,c("remstats","aomstats"),TRUE))){
        stats_attr_method <- attr(stats,"method")
        if(is.null(stats_attr_method)){
            stop("attribute 'method' not found inside object 'remstats'. Input argument 'stats' must be an object of class 'remstats' from the package 'remstats' (>=3.2.0)")
        }
        omit_dyad_receiver <- NULL
        if(stats_attr_method == "pe"){
            if(!is.null(attr(reh,"evenly_spaced_interevent_time"))){
                reh$intereventTime <- attr(reh,"evenly_spaced_interevent_time")
            }
            if(!is.null(reh$E)){ # reh$E is NULL only when there are no simultaneous events
                reh$M <- reh$E # overwriting dimension (we can do it because remstimate works only with reh$M so if the method is "pt", reh$M will remain so. For method "pe" we assign reh$E to reh$M
            }
        }
        if(!is.null(attr(stats,"subset"))){
            start_stop <- as.numeric(unlist(attr(stats,"subset")))
        }
        else{
            start_stop <- c(1,reh$M)
        }
    }
    else{
        stop("'stats' must be a 'remstats' object from the package 'remstats' (>= 3.2.0), suitable for tie-oriented modeling ('tomstats') or actor-oriented modeling ('aomstats')")
    }

    # ... stats 
    model_formula <- variables_names <- where_is_baseline <- NULL
    if(model == "tie")
    {   
        if(all(inherits(stats,c("remstats","tomstats"),TRUE))){
            if(!is.null(dimnames(stats)[[3]])){
                variables_names <- dimnames(stats)[[3]]
            }
            if(is.null(attr(stats,"formula"))){
                model_formula <- stats::as.formula(paste("~ ",paste(variables_names,collapse=" + ")))
            }
            else{
                model_formula <- attr(stats,"formula")
            }
            # is there a baseline term?
            if(any(tolower(variables_names) %in% c("baseline"))){
                where_is_baseline <- which(variables_names == "baseline")
            }
            stats <- aperm(stats, perm = c(2,3,1)) # stats reshaped in [D*U*M]
            # start and stop for tie-oriented model
            if(stats_attr_method == "pt"){
                attr(reh,"dyadID") <- attr(reh,"dyadID")[start_stop[1]:start_stop[2]] # this is already working for dyadIDactive because we reassing attribute dyadID in line 123
                attr(reh,"actor1ID") <- attr(reh,"actor1ID")[start_stop[1]:start_stop[2]]
                attr(reh,"actor2ID") <- attr(reh,"actor2ID")[start_stop[1]:start_stop[2]] 
                if((length(reh$omit_dyad)>0) & !is.null(attr(reh,"indices_simultaneous_events"))){
                    reh$omit_dyad$time <-  reh$omit_dyad$time[-attr(reh,"indices_simultaneous_events")] 
                }
            }
            else if(stats_attr_method =="pe"){
                attr(reh,"dyadID") <- unlist(attr(reh,"dyadID"))[start_stop[1]:start_stop[2]] # this is already working for dyadIDactive because we reassing attribute dyadID in line 123
                attr(reh,"actor1ID") <- unlist(attr(reh,"actor1ID"))[start_stop[1]:start_stop[2]]
                attr(reh,"actor2ID") <- unlist(attr(reh,"actor2ID"))[start_stop[1]:start_stop[2]] 
            } 
            reh$intereventTime <- reh$intereventTime[start_stop[1]:start_stop[2]] # in line 168 we already re-assigned the intereventTime variable in case of method="pe", so line 224 is a valid processing for "pt" and "pe"
            reh$M <- diff(start_stop)+1
            if(length(reh$omit_dyad)>0){
                reh$omit_dyad$time <- reh$omit_dyad$time[start_stop[1]:start_stop[2]]
            }
        }
        else if(all(inherits(stats,c("remstats","aomstats"),TRUE))){
            stop("'remstats' object supplied cannot work for tie-oriented modeling")
        }

        # .. check on dimensions
        if(length(attr(reh,"dyadID")) != dim(stats)[3]){ # if the dimension of the processed intereventTime are different from the dimensions of the input stats object, then throw error
            stop("the number of time points (or number of events) doesn't match the (row) dimension of the 'remstats' object")
        }
    }
    if(model == "actor") 
    {
        model_formula <- list() # becomes a list
        if(all(inherits(stats,c("remstats","aomstats"),TRUE))){
            variables_rate <- variables_choice <- NULL
            if(!is.null(stats$sender_stats)){ # sender model is specified
                variables_rate <- dimnames(stats$sender_stats)[[3]]
                if(!is.null(attr(stats,"formula")$rate)){
                    model_formula[["rate_model_formula"]] <- attr(stats,"formula")$rate
                }
                else{
                    model_formula[["rate_model_formula"]] <- stats::as.formula(paste("~ ",paste(variables_rate,collapse=" + ")))
                }
                # is there a baseline term?
                if(any(tolower(variables_rate) %in% c("baseline"))){
                    where_is_baseline <- which(variables_rate == "baseline")
                }
                stats$sender_stats <- aperm(stats$sender_stats, perm = c(2,3,1)) # stats reshaped in [N*U*M]
            }
            if(!is.null(stats$receiver_stats)){ # receiver model is specified
                variables_choice <- dimnames(stats$receiver_stats)[[3]] 
                if(!is.null(attr(stats,"formula")$choice)){
                    model_formula[["choice_model_formula"]] <- attr(stats,"formula")$choice
                }
                else{
                    model_formula[["choice_model_formula"]] <- stats::as.formula(paste("~ ",paste(variables_choice,collapse=" + ")))
                }
                stats$receiver_stats <- aperm(stats$receiver_stats, perm = c(2,3,1)) # stats reshaped in [N*U*M]
            }

            # vector of variable names and list of model formulas
            variables_names <- list(sender_model = variables_rate, receiver_model = variables_choice)

            # start and stop for actor-oriented model
            if(stats_attr_method == "pt"){
                attr(reh,"actor1ID") <- attr(reh,"actor1ID")[start_stop[1]:start_stop[2]]
                attr(reh,"actor2ID") <- unlist(attr(reh,"actor2ID")[start_stop[1]:start_stop[2]]) # unlist here because of receiver-choice model
                if(is.null(reh$E)){ # this scenario can still happen
                    reh$E <- reh$M
                }
                if(length(reh$omit_dyad)>0){

                    # for the receiver model 
                    if(!is.null(stats$receiver_stats)){
                        start_stop_time <- unique(reh$edgelist$time)[start_stop]
                        lb_time <- min(which(reh$edgelist$time>=start_stop_time[1]))
                        ub_time <- max(which(reh$edgelist$time<=start_stop_time[2]))
                        omit_dyad_receiver <- list(time = reh$omit_dyad$time[lb_time:ub_time], riskset = reh$omit_dyad$riskset) # for the receiver model
                    }

                    # for the sender model (we process now the sender model because this will modify the reh$omit_dyad$time object)
                    if(!is.null(stats$sender_stats)){
                        if(!is.null(attr(reh,"indices_simultaneous_events"))){
                            reh$omit_dyad$time <-  reh$omit_dyad$time[-attr(reh,"indices_simultaneous_events")][start_stop[1]:start_stop[2]] # for the sender model
                        }
                        else{
                            reh$omit_dyad$time <-  reh$omit_dyad$time[start_stop[1]:start_stop[2]] # for the sender model
                        }
                    }
                }
            }
            else if(stats_attr_method == "pe"){
                attr(reh,"actor1ID") <- unlist(attr(reh,"actor1ID"))[start_stop[1]:start_stop[2]]
                attr(reh,"actor2ID") <- unlist(attr(reh,"actor2ID"))[start_stop[1]:start_stop[2]]
                if(length(reh$omit_dyad)>0){
                    if(!is.null(stats$receiver_stats)){
                        omit_dyad_receiver <- list(time = reh$omit_dyad$time[start_stop[1]:start_stop[2]]  , riskset = reh$omit_dyad$riskset)
                    } 
                    reh$omit_dyad$time <- reh$omit_dyad$time[start_stop[1]:start_stop[2]]        
                }
            } 
            if(!ordinal){
                reh$intereventTime <- reh$intereventTime[start_stop[1]:start_stop[2]]
            }
            reh$M <- diff(start_stop)+1

            # .. check on dimensions
            no_correct_dimensions <- FALSE
            if(!is.null(stats$sender_stats)){
                if((length(attr(reh,"actor1ID")) != dim(stats$sender_stats)[3])){ # if the dimension of the processed intereventTime is different from the dimensions of the input stats object, then throw error
                    no_correct_dimensions <- TRUE
                }
            }
            if(!is.null(stats$receiver_stats)){
                if((length(attr(reh,"actor2ID")) != dim(stats$receiver_stats)[3])){ # if the dimension of the edgelist is different from the dimensions of the input stats object, then throw error
                    no_correct_dimensions <- TRUE
                }
            }
            if(no_correct_dimensions){
                stop("the number of time points (or number of events) doesn't match the (row) dimension of the 'remstats' object")
            }

        }
        else if(all(inherits(stats,c("remstats","tomstats"),TRUE))){
            stop("'remstats' object supplied cannot work for actor-oriented modeling")
        }
    }

    # ... adjusting the intereventTime
    if(ordinal){
        reh$intereventTime <- c(1) # we can assign a vector of length 1 because the intereventTime will not be used from any ordinal likelihood
    }

    ## diagnostics per model type (tie-oriented and actor-oriented)

    if(attr(object,"model") == "tie"){ # tie-oriented modeling
        # [[CHECK]] on variables names from object and stats (throw an error if they do not match)
        # [[CHECK]] on reh and object characteristics (throw an error if they do not match)
        variables_names <- attr(object, "statistics")
        where_is_baseline <- attr(object,"where_is_baseline")
        select_vars <- if(is.null(where_is_baseline)) 1:length(variables_names) else c(1:length(variables_names))[-where_is_baseline]
        baseline_value <- if(is.null(where_is_baseline)) 0 else as.vector(object$coefficients)[where_is_baseline]
        stats <- if(length(select_vars)==1) array(stats[,select_vars,],dim=c(dim(stats)[1],1,dim(stats)[3])) else stats[,select_vars,]
        diagnostics <- list()
        diagnostics$residuals <- computeDiagnostics(pars = as.vector(object$coefficients)[select_vars],
                                                stats = stats,
                                                actor1 = list(),
                                                actor2 = list(),
                                                dyad = attr(reh,"dyadID"),
                                                omit_dyad = reh$omit_dyad,
                                                model = attr(reh,"model"),
                                                ncores = attr(object,"ncores"),
                                                baseline = baseline_value,
                                                N = reh$N
                                                )
        colnames(diagnostics$residuals$smoothing_weights) <- variables_names[select_vars]
        # lambdas (event rates)
        diagnostics$rates <- diagnostics$residuals$rates   
        diagnostics$residuals$rates <- NULL  
    }
    else if(attr(object,"model") == "actor"){ # actor-oriented modeling
        # [[CHECK]] on variables names from object and stats (throw an error if they do not match)
        # [[CHECK]] on reh and object characteristics (throw an error if they do not match)
        variables_names <- attr(object, "statistics")
        where_is_baseline <- attr(object,"where_is_baseline")
        senderRate <- c(TRUE,FALSE)
        which_model <- c("sender_model","receiver_model")
        which_stats <- c("sender_stats","receiver_stats") 
        actor1ID_condition <- c(TRUE,FALSE)
        if(stats_attr_method == "pe"){
            actor1ID_condition <- as.logical(actor1ID_condition*FALSE) # actor1ID is needed unlist both for sender and receiver model
        }
        diagnostics <- list()
        # residuals
        for(i in 1:2){
            if(!is.null(stats[[which_stats[i]]])){
                actor1ID_ls <- if(actor1ID_condition[i]) attr(reh,"actor1ID") else unlist(attr(reh,"actor1ID"))
                omit_dyad_actor <- if(senderRate[i]) reh$omit_dyad else omit_dyad_receiver
                diagnostics[[which_model[i]]] <- list()
                baseline_value <- 0
                select_vars <- c(1:dim(stats[[which_stats[i]]])[2])
                if(senderRate[i]){ # only for sender model
                    baseline_value <- if(is.null(where_is_baseline)) 0 else as.vector(object[[which_model[i]]]$coefficients)[where_is_baseline]
                    select_vars <- if(is.null(where_is_baseline)) c(1:dim(stats[[which_stats[i]]])[2]) else c(1:dim(stats[[which_stats[i]]])[2])[-where_is_baseline]
                }
                stats[[which_stats[i]]] <- if(length(select_vars)==1) array(stats[[which_stats[i]]][,select_vars,],dim=c(dim(stats[[which_stats[i]]])[1],1,dim(stats[[which_stats[i]]])[3])) else stats[[which_stats[i]]][,select_vars,]
                diagnostics[[which_model[i]]] <- list()
                diagnostics[[which_model[i]]]$residuals <- computeDiagnostics(pars = as.vector(object[[which_model[i]]]$coefficients)[select_vars],
                                                        stats = stats[[which_stats[i]]],
                                                        actor1 = actor1ID_ls,
                                                        actor2 = attr(reh,"actor2ID"),
                                                        dyad = list(),
                                                        omit_dyad = omit_dyad_actor,
                                                        model = attr(reh,"model"),
                                                        N = reh$N,
                                                        senderRate = senderRate[i],
                                                        ncores = attr(object,"ncores"),
                                                        baseline = baseline_value)                                   
                colnames(diagnostics[[which_model[i]]]$residuals$smoothing_weights) <- if(senderRate[i]) variables_names[["sender_model"]][select_vars] else variables_names[["receiver_model"]][select_vars]
                # lambdas (event rates)
                diagnostics[[which_model[i]]]$rates <- diagnostics[[which_model[i]]]$residuals$rates   
                diagnostics[[which_model[i]]]$residuals$rates <- NULL 
            }
        }
    }
    diagnostics$.reh.processed <- reh   
    diagnostics$.reh.processed$stats.method <- stats_attr_method
    class(diagnostics) <- c("remstimate","diagnostics")
    return(diagnostics) 
}


#######################################################################################
#######################################################################################


# predict.remstimate
#' @title predict.remstimate
#' @rdname predict.remstimate
#' @description A function that returns in-sample A-steps ahead predictions given a 'remstimate' object.
#' @param object is a \code{remstimate} object 
#' @param A is the number of steps ahead to predict
#' @param diagnostics a object of class \code{"remstimate","diagnostics"}
#' @param ... further arguments to be passed to the 'predict' method depending on the estimator used
#' @method predict remstimate
#' @export
#' 
#' @examples 
#' 
#' # No examples available at the moment
#' 
predict.remstimate <- function(object, A, diagnostics, ...)
{
    # - In-sample predictions

    
    # A-steps ahead A = 1, ... 
    # if in-sample, comparison with actual observed dyads
    #if out-of-sample, find comparison measures

    # Out-of-sample predictions 
    #  code here an out-of-sample approach for predictions

    return(paste('this function at the moment does nothing'))
}

#######################################################################################
#######################################################################################


# plot.remstimate
#' @title plot.remstimate
#' @rdname plot.remstimate
#' @description A function that returns a plot of diagnostics given a 'remstimate' object and depending on the 'approach' attribute.
#' @param x is a \code{remstimate} object
#' @param which one or more numbers between 1 and 2. Plots described in order: (1) two plots: a Q-Q plot of the waiting times where theoretical quantiles (Exponential distribution with rate 1) are plotted against observed quantiles (these are calculated as the multiplication at each time point between the sum of the event rates and the corresponding waiting time, which should be distributed as an exponential with rate 1). Next to the q-q plot, a density plot of the rescaled waiting times (in red) vs. the theoretical distribution (exponential distribution with rate 1, in black). The observed density is truncated at the 99th percentile of the waiting times, (2) standardized Schoenfeld's residuals (per each variable in the model, excluding the baseline) with smoothed weighted spline (line in red). The Schoenfeld's residuals help understand the potential presence of time dependence of the effects of statistics specified in the model.
#' @param reh 'remify' object, the same used for the 'remstimate' object
#' @param diagnostics is a \code{"remstimate" "diagnostics"} object
#' @param ... further arguments to be passed to the 'plot' method depending: for instance the remstats object with statistics ('stats')
#' @method plot remstimate
#' @export
#' 
#' @examples 
#' 
#' # ------------------------------------ #
#' #       tie-oriented model: "MLE"      #
#' # ------------------------------------ #
#' 
#' # loading data
#' data(tie_data)
#' 
#' # processing event sequence with remify
#' tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
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
#' # diagnostics
#' tie_diagnostics <- diagnostics(object = tie_mle, reh = tie_reh, stats = tie_reh_stats)
#' 
#' # plot
#' plot(x = tie_mle, reh  = tie_reh, diagnostics = tie_diagnostics)
#' 
plot.remstimate <- function(x, 
                            reh, 
                            which = c(1:2), 
                            diagnostics = NULL, ...)
{
    # which plot
    selected <- which
    which <- rep(FALSE,2)
    which[selected] <- TRUE

    if(attr(x,"model") != attr(reh,"model")){
        stop("'x' and 'reh' have different attribute 'model'")
    }
    # if diagnostics is NULL, then look for object 'stats' in '...' argument and compute diagnostics
    if(is.null(diagnostics)){
        additional_input_args <- list(...) # check for stats argument
        if(!any(names(additional_input_args) %in% "stats")){
            stop("'stats' must be provided if argument 'diagnostics' is NULL")
        }
        else{
            diagnostics <- diagnostics(object = x, reh = reh, stats = additional_input_args$stats)
        }
    }
    else if(!all(inherits(diagnostics,c("remstimate","diagnostics"),TRUE))){
            stop("'diagnostics' must be an object of class 'remstimate' 'diagnostics'")
    }

    # overwriting reh with one coming from diagnostics(), because there pt/pe methods and start/stop from remstats are already processed in there
    reh <- diagnostics$.reh.processed

    # saving current graphic parameters
    op <- par(no.readonly = TRUE)
    on.exit(expr = par(op))

    # tie-oriented modeling
    if(attr(x,"model") == "tie"){ 

        # (1) waiting times vs. theoretical distribution
        if(which[1L]){
            if(!attr(reh,"ordinal")){
                sum_rates <- lapply(diagnostics$rates,sum)
                observed <- sort(reh$intereventTime*unlist(sum_rates))
                theoretical <- stats::qexp(p = c(1:reh$M)/reh$M, rate = 1)
                par(mfrow=c(1,2))
                plot(theoretical,observed, xlab = "Theoretical Quantiles", 
                ylab = "Observed Quantiles",cex=0.8) # use bquote() / Observed quatiles : ~(t[m]-t[m-1])*Sigma[e](lambda[e](t[m])) / Theoretical Quantiles ~Exp(1)
                mtext(text = "Q-Q waiting times", side = 3, line = 2,cex=1.5)
                abline(a=0,b=1,lty=2,lwd=1.5)
                density_observed <- density(observed)
                density_observed$y <- density_observed$y/max(density_observed$y) # rescaling max density to 1 (to compare with dexp)
                curve(dexp,from=min(observed),to=as.numeric(quantile(observed,probs=c(0.99))),col=1,lwd=1.5,xlab="waiting times")
                lines(density_observed,col=2,lwd=1.5)
                mtext(text = "Density plot of waiting times", side = 3, line = 2,cex=1.5)
                legend("topright", legend = c("Theoretical density","Observed density"), lwd=c(1.5,1.5), lty = c(1,1), col = c(1,2))
                par(op)
            }
        }

        # (2) standardized Schoenfeld's residuals
        if(which[2L]){
            P <- dim(diagnostics$residuals$standardized_residuals[[1]])[2] # number of statistics
            n_repeats_per_time_point <- rep(1,dim(diagnostics$residuals$standardized_residuals)[1]) # for attr(residuals, "stats.method") = "pe"
            if(reh$stats.method == "pt"){
                n_repeats_per_time_point <- sapply(1:dim(diagnostics$residuals$standardized_residuals)[1], function(v) dim(diagnostics$residuals$standardized_residuals[[v]])[1])
            }
            if(!attr(reh,"ordinal")){
                t_p <- rep(cumsum(reh$intereventTime),n_repeats_per_time_point)
            }
            else{
                t_p <- rep(cumsum(rep(1,reh$M)),n_repeats_per_time_point)
            }
            for(p in 1:P){
                y_p <-  unlist(sapply(1:dim(diagnostics$residuals$standardized_residuals)[1], function(v) diagnostics$residuals$standardized_residuals[[v]][,p]))
                qrt_p <- quantile(y_p, probs=c(0.25,0.75))
                lb_p <- qrt_p[1] - diff(qrt_p)*1.5
                ub_p <- qrt_p[2] + diff(qrt_p)*1.5
                ylim_p <- c(min(y_p),max(y_p))
                #if(length(which(y_p<ub_p & y_p>lb_p)) >= floor(0.95*length(y_p))){ # reduce the ylim only if the bound lb_p and ub_p include more than the 90% of the observations
                    ylim_p <- c(lb_p,ub_p)
                #}
                w_p <- rep(diagnostics$residuals$smoothing_weights[,p],n_repeats_per_time_point)
                w_p[w_p<0] <- w_p[w_p>1e04] <- 0.0
                if(all(w_p <= 0.0)){
                    w_p <- rep(1/length(w_p),length(w_p)) # assigning equal weights
                }
                par(mfrow=c(1,1))
                plot(t_p,y_p,xlab = "time", ylab = "scaled Schoenfeld's residuals",ylim=ylim_p,col=grDevices::rgb(128,128,128,200,maxColorValue = 255)) # standardized Schoenfeld's residuals
                lines(smooth.spline(x = t_p, y = y_p,w=w_p,cv=NA),lwd=3.5,col=2) # smoothed weighted spline of the residuals
                abline(h=0,col="black",lwd=3,lty=2)
                mtext(text = colnames(diagnostics$residuals$smoothing_weights)[p], side = 3, line = 1,cex=1.5)
                par(op)
            }
        }
    }

    # actor-oriented modeling
    else if(attr(x,"model") == "actor"){ 
        which_model <- c("sender_model","receiver_model")
        title_model <- c("Rate model (sender)","Choice model (receiver)")
        for(i in 1:2){
            if(!is.null(x[[which_model[i]]])){

                # (1) waiting times vs. theoretical distribution
                if(which[1L]){
                    if(!attr(reh,"ordinal") & i==1){
                        sum_rates <- lapply(diagnostics[[which_model[i]]]$rates,sum)
                        observed <- sort(reh$intereventTime*unlist(sum_rates))
                        theoretical <- stats::qexp(p = c(1:reh$M)/reh$M, rate = 1)
                        par(mfrow=c(1,2))
                        plot(theoretical,observed, xlab = "Theoretical Quantiles", 
                        ylab = "Observed Quantiles",cex=0.8) # use bquote() / Observed quatiles : ~(t[m]-t[m-1])*Sigma[e](lambda[e](t[m])) / Theoretical Quantiles ~Exp(1)
                        mtext(text = "Q-Q waiting times", side = 3, line = 2,cex=1.5)
                        mtext(text = title_model[i], side = 3, line = 1,cex=1)
                        abline(a=0,b=1,lty=2,lwd=1.5)
                        density_observed <- density(observed)
                        density_observed$y <- density_observed$y/max(density_observed$y) # rescaling max density to 1 (to compare with dexp)
                        curve(dexp,from=min(observed),to=as.numeric(quantile(observed,probs=c(0.99))),col=1,lwd=1.5,xlab="waiting times",ylab="Density")
                        lines(density_observed,col=2,lwd=1.5)
                        mtext(text = "Density plot of waiting times", side = 3, line = 2,cex=1.5)
                        mtext(text = title_model[i], side = 3, line = 1,cex=1)
                        legend("topright", legend = c("Theoretical density","Observed density"), lwd=c(1.5,1.5), lty = c(1,1), col = c(1,2))
                        par(op)
                    }
                }

                # (2) standardized Schoenfeld's residuals
                if(which[2L]){
                    P <- dim(diagnostics[[which_model[i]]]$residuals$standardized_residuals[[1]])[2] # number of statistics
                    n_repeats_per_time_point <- rep(1,dim(diagnostics[[which_model[i]]]$residuals$standardized_residuals)[1]) # for attr(residuals, "stats.method") = "pe"
                    if(i == 1){
                        if(reh$stats.method == "pt"){
                            n_repeats_per_time_point <- sapply(1:dim(diagnostics[[which_model[i]]]$residuals$standardized_residuals)[1], function(v) dim(diagnostics[[which_model[i]]]$residuals$standardized_residuals[[v]])[1])
                        }
                        if(!attr(reh,"ordinal")){
                            t_p <- rep(cumsum(reh$intereventTime),n_repeats_per_time_point)
                        }
                        else{
                            t_p <- rep(cumsum(rep(1,reh$M)),n_repeats_per_time_point)
                        }
                    }
                    else if(i == 2){
                        t_p <- cumsum(n_repeats_per_time_point)
                    }
                    for(p in 1:P){
                        y_p <- unlist(sapply(1:dim(diagnostics[[which_model[i]]]$residuals$standardized_residuals)[1], function(v) diagnostics[[which_model[i]]]$residuals$standardized_residuals[[v]][,p]))
                        qrt_p <- quantile(y_p, probs=c(0.25,0.75))
                        lb_p <- qrt_p[1] - diff(qrt_p)*1.5
                        ub_p <- qrt_p[2] + diff(qrt_p)*1.5
                        #ylim_p <- c(min(y_p),max(y_p))
                        #if(length(which(y_p<ub_p & y_p>lb_p)) >= floor(0.95*length(y_p))){ # reduce the ylim only if the bound lb_p and ub_p include more than the 95% of the observations
                            ylim_p <- c(lb_p,ub_p)
                        #}
                        w_p <- rep(diagnostics[[which_model[i]]]$residuals$smoothing_weights[,p],n_repeats_per_time_point)
                        w_p[w_p<0] <- w_p[w_p>1e04] <- 0.0
                        if(all(w_p <= 0.0)){
                            w_p <- rep(1/length(w_p),length(w_p)) # assigning equal weights
                        }
                        par(mfrow=c(1,1))
                        plot(t_p,y_p,xlab = "time", ylab = "scaled Schoenfeld's residuals",ylim=ylim_p,col=grDevices::rgb(128,128,128,200,maxColorValue = 255)) # standardized Schoenfeld's residuals
                        lines(smooth.spline(x = t_p, y = y_p,w=w_p,cv=NA),lwd=3.5,col=2) # smoothed weighted spline of the residuals
                        abline(h=0,col="black",lwd=3,lty=2)
                        mtext(text = colnames(diagnostics[[which_model[i]]]$residuals$smoothing_weights)[p], side = 3, line = 2,cex=1.5)
                        mtext(text = title_model[i], side = 3, line = 1, cex = 1)
                        par(op)
                    }
                }
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
#' data(tie_data)
#' 
#' # processing event sequence with remify
#' tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
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
#' data(tie_data)
#' 
#' # processing event sequence with remify
#' tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
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
#' data(tie_data)
#' 
#' # processing event sequence with remify
#' tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
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




