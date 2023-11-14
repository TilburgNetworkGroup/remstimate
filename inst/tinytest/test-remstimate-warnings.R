## testing warnings of remstimate ##
  
# loading data
data(tie_reh)

# specifying linear predictor
tie_model <- ~ 1 + remstats::inertia()

# calculating statistics
tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)


## seed = c(1234,4321), nchains = 2L, method = "HMC"
expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        seed = c(1234,4321),
                        nchains = 2L,
                        nsim = 10L,
                        burnin = 5L),
"`seed` length is greater than 1. Considering only the first element",
fixed = TRUE)  

## seed = c(1234,4321), method = "BSIR"                      
expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "BSIR",
                        nsim = 10L,
                        seed = c(1234,4321)),
"`seed` length is greater than 1. Considering only the first element",
fixed = TRUE
)   

## method = "HMC", thin = NULL, nsim >= 100
expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        nchains = 2L,
                        nsim = 100L,
                        burnin = 5L,
                        thin = NULL),
"'thin' parameter undefined. Using thin = 10",
fixed = TRUE)  

## method = "HMC", thin = NULL, nsim < 100
expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        nchains = 2L,
                        nsim = 10L,
                        burnin = 5L,
                        thin = NULL),
"'nsim' is less than 100. No thinning is applied",
fixed = TRUE)  

## method = "HMC", thin = 1, nsim >= 100
expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        nchains = 2L,
                        nsim = 100,
                        burnin = 5L,
                        thin = 1L ),
"`thin` value is 1. No thinning applied to chains",
fixed = TRUE)  

## method = "HMC", length(L)>1 OR L <= 1 OR !is.numeric(L)
expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        nchains = 2L,
                        nsim = 10L,
                        burnin = 5L,
                        thin = 1L,
                        L = c(1,1)),
"input 'L' must be a positive (integer) number . 'L' is set to its default value: 50",
fixed = TRUE)  
expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        nchains = 2L,
                        nsim = 10L,
                        burnin = 5L,
                        thin = 1L,
                        L = 0),
"input 'L' must be a positive (integer) number . 'L' is set to its default value: 50",
fixed = TRUE) 
expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        nchains = 2L,
                        nsim = 10L,
                        burnin = 5L,
                        thin = 1L,
                        L = FALSE),
"input 'L' must be a positive (integer) number . 'L' is set to its default value: 50",
fixed = TRUE) 

## method = "HMC", length(epsilon)>1 OR epsilon <= 0 OR !is.numeric(epsilon)
expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        nchains = 2L,
                        nsim = 10L,
                        burnin = 5L,
                        thin = 1L,
                        epsilon = c(0.001,0.001)),
"input 'epsilon' must be a positive number. 'epsilon' is set to its default value: 0.1/L",
fixed = TRUE) 

expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        nchains = 2L,
                        nsim = 10L,
                        burnin = 5L,
                        thin = 1L,
                        epsilon = -1),
"input 'epsilon' must be a positive number. 'epsilon' is set to its default value: 0.1/L",
fixed = TRUE) 

expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "HMC",
                        nchains = 2L,
                        nsim = 10L,
                        burnin = 5L,
                        thin = 1L,
                        epsilon = FALSE),
"input 'epsilon' must be a positive number. 'epsilon' is set to its default value: 0.1/L",
fixed = TRUE) 


## method "GDADAMAX" and epsilon <= 0 OR length(epsilon) > 1 OR !is.numeric(epsilon)
expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "GDADAMAX",
                        epsilon = -1),
"'epsilon' is set to its default value: 0.001",
fixed = TRUE)  
expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "GDADAMAX",
                        epsilon = c(1,1)),
"'epsilon' is set to its default value: 0.001",
fixed = TRUE) 
expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        method = "GDADAMAX",
                        epsilon = TRUE),
"'epsilon' is set to its default value: 0.001",
fixed = TRUE) 