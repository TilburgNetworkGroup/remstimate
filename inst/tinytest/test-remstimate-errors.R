## testing errors of remstimate ##

# loading data for tie-oriented modeling, defining linear predictor and computing statistics
data(tie_data)

# processing data
tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")

# specifying linear predictor
tie_model <- ~ 1 + remstats::inertia()

# calculating statistics
tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)

## input `reh` is a data.frame 
mle_loc <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "MLE",
                        ncores = 1L)

# loading data for actor-oriented modeling, defining linear predictor and computing statistics                     
data(ao_data)

# processing data with simultaneous events (two events observed at the same time point)
ao_reh <- remify::remify(edgelist = ao_data$edgelist, model = "actor")

# specifying linear predictor (for rate and choice model)
rate_model <- ~ 1 + remstats::indegreeSender()
choice_model <- ~ remstats::inertia() + remstats::reciprocity()

# calculating statistics
ao_reh_stats <- remstats::remstats(reh = ao_reh, sender_effects = rate_model, receiver_effects = choice_model)

# ao_mle_loc 
ao_mle_loc <- remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        method = "MLE",
                        ncores = 1L)

# testing errors

## input `reh` is not a `remify` object
expect_error(remstimate::remstimate(reh = data.matrix(ao_data$edgelist),
                        stats = tie_reh_stats,
                        method = "MLE",
                        ncores = 1L),
"'reh' must be a 'remify' object (see ?remify::remify).",
fixed = TRUE
)         

## input `method` is not available or it is mistyped
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "mle",
                        ncores = 1L),
"The `method` specified is not available or it is mistyped.",
fixed = TRUE
)

## directed = FALSE and model = "actor"
attr(ao_reh,"directed") <- FALSE
expect_error(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        method = "MLE",
                        ncores = 1L),
"actor-oriented modeling can't operate on undirected networks",                      
fixed = TRUE
)
attr(ao_reh,"directed") <- TRUE


## tests on the stats object

# actor-oriented modeling and tomstats object
attr(tie_reh,"model") <- "actor" 
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "MLE",
                        ncores = 1L),
"attribute 'model' of input 'reh' and input 'stats' must be the same",
fixed = TRUE
) 

# actor-oriented modeling and stats object that is not a remstats object
class(tie_reh_stats) <- "none"
attr(tie_reh,"model") <- "tie"
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "MLE",
                        ncores = 1L),
"'stats' must be a 'remstats' object from the package 'remstats' (>= 3.2.0), suitable for tie-oriented modeling ('tomstats') or actor-oriented modeling ('aomstats')",
fixed = TRUE
) 

# tie-oriented modeling with a aomstats object
class(tie_reh_stats) <- c("aomstats","remstats")
attr(tie_reh,"model") <- "tie" 
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "MLE",
                        ncores = 1L),
"'remstats' object supplied cannot work for tie-oriented modeling",
fixed = TRUE
) 

# tie-oriented modeling with a stats object that is not a remstats object
class(tie_reh_stats) <- ""
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "MLE",
                        ncores = 1L),
"'stats' must be a 'remstats' object from the package 'remstats' (>= 3.2.0), suitable for tie-oriented modeling ('tomstats') or actor-oriented modeling ('aomstats')",
fixed = TRUE
) 
class(tie_reh_stats) <- c("tomstats","remstats")

# model == "tie" AND method %in% c("GDADAMAX","HMC") AND !is.null(init) AND !is.vector(init)
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "GDADAMAX",
                        ncores = 1L,
                        init = matrix(1:16,nrow=4,ncol=4)),
"'init' must be a vector with the starting values of the effects of the statistics",
fixed = TRUE
) 

# model == "tie" AND method %in% c("GDADAMAX","HMC") AND !is.null(init) AND length(init)!=dim(stats)[2]
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "GDADAMAX",
                        ncores = 1L,
                        init = c(1)),
"length of vector 'init' must be equal to the number of statistics in 'stats'",
fixed = TRUE
) 

#  model == "tie" - plot function of not available effects 
expect_error(plot(x = mle_loc, reh = tie_reh, which = 2, effects = "reciprocity", stats = tie_reh_stats),
"effects not found in object 'remstimate'",
fixed = TRUE
)

# model == "actor" - plot function of not available effects
expect_silent(plot(1:2))
expect_error(plot(x = ao_mle_loc, reh = ao_reh, which = 2, sender_effects = "reciprocity", receiver_effects = "indegreeSender", stats = ao_reh_stats),
"effects not found in object 'remstimate'",
fixed = TRUE
) 


# model == "actor" AND method %in% c("GDADAMAX","HMC") AND !is.null(init) AND !is.list(init)
expect_error(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        method = "GDADAMAX",
                        ncores = 1L,
                        init = matrix(1:16,nrow=4,ncol=4)),
"'init' must be a list of two vectors named 'sender_model' and 'receiver_model', each one with the starting values of the effects of the statistics according to the two different models (rate model and choice model)",
fixed = TRUE
) 

# model == "actor" AND method %in% c("GDADAMAX","HMC") AND !is.null(init) AND !all(c("sender_model","receiver_model") %in% names(init))
expect_error(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        method = "GDADAMAX",
                        ncores = 1L,
                        init = list(s = c(1,1),r = c(1,1))),
"'init' must be a list of two vectors named 'sender_model' and 'receiver_model', each one with the starting values of the effects of the statistics according to the two different models (rate model and choice model)",
fixed = TRUE
) 

# model == "actor" AND method %in% c("GDADAMAX","HMC") AND !is.null(init) AND length(init$sender_model)!=dim(stats$sender_model)[2] | length(init$receiver_model)!=dim(stats$receiver_model)[2])
expect_error(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        method = "GDADAMAX",
                        ncores = 1L,
                        init = list(sender_model = c(1),receiver_model = c(1))),
"each element of list 'init' must be equal to number of statistics according to the arrays supplied in 'stats'",
fixed = TRUE
) 


## input `epochs` is negative
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "GDADAMAX",
                        ncores = 1L,
                        epochs = -29),
"'epoch' must be a positive number",
fixed = TRUE
)  

## input `epochs` is a character
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "GDADAMAX",
                        ncores = 1L,
                        epochs = "ss"),
"'epoch' must be a positive number",
fixed = TRUE
)

# estimate with epochs = NULL and epsilon = NULL
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "GDADAMAX",
                        ncores = 1L,
                        epochs = NULL,
                        epsilon = NULL))


## when parallel::detectCores() == 2 and input `ncores` is greater than 1
ncores_loc <- if(parallel::detectCores() <= 2L) 10L else 2e03L
# run the test only if there are 2 cores in the machine
expect_error(remstimate::remstimate(reh = tie_reh,
                      stats = tie_reh_stats,
                      method = "MLE",
                      ncores = ncores_loc),
"'ncores' is recommended to be set at most to: floor(parallel::detectCores()-2L)",
fixed = TRUE
)


# estimate with NULL ncores 
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "MLE",
                        ncores = NULL))

# tests on seed, nchains, nsim, burnin, thin, L and espilon parameters

## seed == NULL, nchains = 3L and method = "HMC"
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        seed = NULL,
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L))
## seed == NULL, nchains = NULL and method = "HMC"                      
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        seed = NULL,
                        nchains = NULL,
                        nsim = 10L,
                        burnin = 5L))         

## seed = 1234, nchains = NULL, method = "HMC"   
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        seed = 1234,
                        nchains = NULL,
                        nsim = 10L,
                        burnin = 5L))

## seed = NULL, method = "BSIR" 
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "BSIR",
                        seed = NULL,
                        nsim = 10L))

## nsim = NULL, method  = "BSIR" (or "HMC")                    
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "BSIR",
                        nsim = NULL,
                        seed = 1234)) 
## burnin = NULL, method = "HMC" 
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        nsim = 10L,
                        thin = 5L,
                        burnin = NULL,
                        seed = 1234))          

## thin <= 0, method == "HMC"                               
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        nsim = 10L,
                        burnin = 5L,
                        thin = 0),
"`thin` value must be positive. Set thin = 1 if no thinning of chains is required",
fixed = TRUE
) 

# model of input reh is not actor or tie
## thin <= 0, method == "HMC"           
attr(tie_reh,"model") <- ""                    
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        nsim = 10L,
                        burnin = 5L,
                        thin = 0),
"attribute 'model' of input 'reh' must be either 'actor' or 'tie'",
fixed = TRUE
) 
attr(tie_reh,"model") <- "tie" 

## attribute "formula" is NULL
tie_mle <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "MLE",
                        ncores = 1L)
attr(tie_mle,"formula") <- NULL
expect_error(print(tie_mle),
"invalid 'remstimate' object:  no 'formula' attribute",
fixed = TRUE
) 

# invalid 'remstimate' object
attr(tie_mle, "approach") <- ""
expect_error(summary(tie_mle),
"invalid 'remstimate' object",
fixed = TRUE)


# invalid 'remstimate' object
attr(tie_mle, "model") <- "actor"
expect_error(plot(tie_mle,tie_reh),"'x' and 'reh' have different attribute 'model'",fixed = TRUE)

# stats
attr(tie_mle, "model") <- "tie"
expect_error(plot(x = tie_mle, reh = tie_reh,diagnostics = NULL),"'stats' must be provided if argument 'diagnostics' is NULL",fixed = TRUE)

# diagnostics object
expect_error(plot(x = tie_mle, reh = tie_reh,diagnostics = list(a = "foo")),"'diagnostics' must be an object of class 'remstimate' 'diagnostics'",fixed = TRUE)
 
# number of time points doesn't mach dimensions of remstats object
## tie-oriented modeling
tie_reh_stats <- tie_reh_stats[1:10,,]
class(tie_reh_stats) <- c("tomstats", "remstats")
attr(tie_reh_stats,"method") <- "pt"
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "MLE",
                        ncores = 1L),
"the number of time points (or number of events) doesn't match the (row) dimension of the 'remstats' object",
fixed = TRUE)
## actor-oriented modeling
data(ao_data)
# processing data with simultaneous events (two events observed at the same time point)
ao_reh <- remify::remify(edgelist = ao_data$edgelist, model = "actor")
rate_model <- ~ 1 + remstats::indegreeSender()
choice_model <- ~ remstats::inertia() + remstats::reciprocity()
ao_reh_stats <- remstats::remstats(reh = ao_reh, sender_effects = rate_model, receiver_effects = choice_model)
ao_reh_stats$sender_stats <- ao_reh_stats$sender_stats[1:10,,]
ao_reh_stats$receiver_stats <- ao_reh_stats$receiver_stats[1:10,,]
class(ao_reh_stats) <- c("aomstats","remstats")
attr(ao_reh_stats,"method") <- "pt"
expect_error(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        method = "MLE",
                        ncores = 1L),
"the number of time points (or number of events) doesn't match the (row) dimension of the 'remstats' object",
fixed = TRUE)

# actor-oriented modeling with a tomstats object
ao_reh_stats <- remstats::remstats(reh = ao_reh, sender_effects = rate_model, receiver_effects = choice_model)
class(ao_reh_stats) <- c("tomstats","remstats")
expect_error(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        method = "MLE",
                        ncores = 1L),
"'remstats' object supplied cannot work for actor-oriented modeling",
fixed = TRUE
) 

# method of remstats is not available
attr(tie_reh_stats,"method") <- NULL
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "MLE",
                        ncores = 1L),
"attribute 'method' not found inside object 'remstats'. Input argument 'stats' must be an object of class 'remstats' from the package 'remstats' (>=3.2.0)",
fixed = TRUE)