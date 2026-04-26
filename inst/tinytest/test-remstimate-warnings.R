## testing warnings of remstimate ##

# loading data
data(tie_data)

# processing data
tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")

# specifying linear predictor
tie_model <- ~ 1 + remstats::inertia()

# calculating statistics
tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)


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

