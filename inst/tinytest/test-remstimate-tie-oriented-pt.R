## testing tie-oriented modeling (methods) ##

library(tinytest)

# loading data
data(tie_data)

set_seed <- 23929

# processing data with simultaneous events (two or more events observed at the same time point)
tie_data$edgelist$time <- floor(tie_data$edgelist$time)
tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")

# specifying linear predictor
tie_model <- ~ 1  + remstats::indegreeSender()+remstats::inertia()+remstats::reciprocity()

# calculating statistics
tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)

# (1) method = "MLE"
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE"))
# testing method = NULL (default is "MLE")
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = NULL))

tie_mle <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE")
expect_silent(tie_mle)
expect_inherits(tie_mle,"remstimate")
expect_length(tie_mle,18)
expect_identical(names(tie_mle),c("coefficients","loglik","gradient","hessian","vcov","se","residual.deviance","null.deviance","model.deviance","df.null","df.model","df.residual","AIC","AICC","BIC","converged","iterations","sampled"))
expect_length(attributes(tie_mle),11)
expect_identical(names(attributes(tie_mle)),c("names","class","formula","model","ordinal","method","approach","statistics","where_is_baseline","ncores","sampled"))
expect_identical(attr(tie_mle,"approach"),"Frequentist")
expect_silent(print(tie_mle))
expect_silent(summary(tie_mle))
#expect_silent(diagnostics(object = tie_mle, reh = tie_reh, stats = tie_reh_stats))
tie_reh_diagnostics <- diagnostics(object = tie_mle, reh = tie_reh, stats = tie_reh_stats)
expect_silent(plot(x = tie_mle,reh = tie_reh, diagnostics = tie_reh_diagnostics))
expect_silent(plot(x = tie_mle,reh = tie_reh, diagnostics = tie_reh_diagnostics, which = 2, effects = "inertia")) # plotting a specific effect
expect_silent(plot(x = tie_mle,reh = tie_reh, stats = tie_reh_stats)) # without supplying diagnostics but supplying the array of stats
expect_silent(AIC(tie_mle))
expect_silent(AICC(tie_mle))
expect_silent(BIC(tie_mle))
expect_error(WAIC(tie_mle))

# error when 'diagnostics' has at least one statistic's name wrong
colnames(tie_reh_diagnostics$residuals$smoothing_weights)[1] <- "INDEGREEsender"
expect_error(plot(x = tie_mle,reh = tie_reh, diagnostics = tie_reh_diagnostics),
"one or more effects not found inside the object 'diagnostics'.",
fixed=TRUE)

# test diagnostics with only baseline
tie_stats_baseline <- remstats::remstats(reh = tie_reh, tie_effects = ~1)
tie_mle_only_baseline <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_stats_baseline,
                        ncores = 1L,
                        method = "MLE")
expect_silent(diagnostics(object = tie_mle_only_baseline, reh = tie_reh, stats = tie_stats_baseline))


# tests on "WAIC" for "MLE"

# WAIC = TRUE
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE))
# WAIC = NULL
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        WAIC = NULL))
# WAIC not a logical scalar
expect_error(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        WAIC = 20))

# nsimWAIC is not a numeric scalar
expect_warning(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = "text"))
# nsimWAIC is supplied
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100))
tie_mle_with_waic <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100)
expect_silent(print(tie_mle_with_waic))
expect_silent(summary(tie_mle_with_waic))
expect_silent(WAIC(tie_mle_with_waic))


# ordinal likelihood: WAIC for ordinal likelihod + testing omit_dyad routine
#omit_dyad <- NULL
#omit_dyad[[1]] <- list(time = c(120,158), dyad = data.frame(actor1=c(NA,"4"),actor2=c("4",NA),type=c(NA,NA))) # excluding actor 4 from the risk set between (observed) time points 120 and 158
tie_reh_ordinal <- remify::remify(edgelist = tie_data$edgelist, model = "tie", ordinal = TRUE)
tie_reh_stats_ordinal <- remstats::remstats(reh = tie_reh_ordinal, tie_effects = tie_model)
expect_silent(remstimate::remstimate(reh = tie_reh_ordinal,
                        stats = tie_reh_stats_ordinal,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100))
# interval likelihood: WAIC for ordinal likelihod + testing omit_dyad routine in estimation (remstimate) and diagnostics method
tie_reh_omit_test <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
expect_silent(remstimate::remstimate(reh = tie_reh_omit_test,
                        stats = tie_reh_stats_ordinal,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100))
tie_mle_omit_test <- remstimate::remstimate(reh = tie_reh_omit_test,
                        stats = tie_reh_stats_ordinal,
                        ncores = 1L)
expect_silent(diagnostics(object = tie_mle_omit_test, reh = tie_reh_omit_test, stats = tie_reh_stats_ordinal))


# (2) method  = "HMC"
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L,
                        seed = set_seed))
tie_hmc <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 3L,
                        nsim = 10L,
                        burnin = 5L,
                        seed = set_seed)
expect_silent(tie_hmc)
expect_inherits(tie_hmc,"remstimate")
expect_length(tie_hmc,9)
expect_identical(names(tie_hmc),c("draws","log_posterior","coefficients","post.mean","vcov","sd","loglik","sampled","df.null"))
expect_length(attributes(tie_hmc),18)
expect_identical(names(attributes(tie_hmc)),c("names","class","formula","model","ordinal","method","approach","statistics","where_is_baseline","ncores","sampled","nsim","nchains","burnin","thin","L","epsilon","seed"))
expect_identical(attr(tie_hmc,"approach"),"Bayesian")
expect_silent(print(tie_hmc))
expect_silent(summary(tie_hmc))
expect_silent(diagnostics(object = tie_hmc, reh = tie_reh, stats = tie_reh_stats))
tie_reh_diagnostics <- diagnostics(object = tie_hmc, reh = tie_reh, stats = tie_reh_stats)
expect_silent(plot(x = tie_hmc,reh = tie_reh, diagnostics = tie_reh_diagnostics))
expect_error(AIC(tie_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(AICC(tie_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(BIC(tie_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(WAIC(tie_hmc))

# tests on "WAIC" for "HMC"
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 100L,
                        burnin = 5L,
                        seed = set_seed,
                        WAIC = TRUE,
                        nsimWAIC = 10))



# testing with omit_dyad
tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie") # removing actor "4" from time=120 to time=148
# calculating statistics
tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)
# (1) method = "MLE"
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE"))

# Risk set "active"

# testing estimation methods with active riskset
tie_reh <- remify::remify(edgelist = tie_data$edgelist,
                            model = "tie",
                            riskset="active")
# calculating statistics
tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)
attr(tie_reh_stats,"formula") <- NULL

# (1) method = "MLE"
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE"))
tie_mle <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE")
expect_silent(diagnostics(object = tie_mle, reh = tie_reh, stats = tie_reh_stats))

# (2) method  = "HMC"
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L,
                        seed = set_seed))



