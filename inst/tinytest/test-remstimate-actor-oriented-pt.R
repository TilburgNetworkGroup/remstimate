## testing actor-oriented modeling ##

# loading data
data(ao_data)

set_seed <- 23929

# processing data with simultaneous events (two events observed at the same time point)
ao_data$edgelist$time <- floor(ao_data$edgelist$time)
ao_data$edgelist$time[seq(5,95,by=5)] <- ao_data$edgelist$time[seq(5,95,by=5)-1]
ao_reh <- remify(edgelist = ao_data$edgelist, model = "actor")

# specifying linear predictor (for rate and choice model)
rate_model <- ~ 1 + indegreeSender()
choice_model <- ~ inertia() + reciprocity()

# calculating statistics
ao_reh_stats <- remstats(reh = ao_reh,
                                    sender_effects = rate_model,
                                    receiver_effects = choice_model)

# (1) method = "MLE"
expect_silent(remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "MLE"))
ao_mle <- remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "MLE")
expect_silent(ao_mle)
expect_inherits(ao_mle,"remstimate")
expect_length(ao_mle,2)
expect_identical(names(ao_mle),c("sender_model","receiver_model"))
expect_length(ao_mle$sender_model,17)
expect_length(ao_mle$receiver_model,17)
expect_identical(names(ao_mle$sender_model),c("coefficients","loglik","gradient","hessian","vcov","se","residual.deviance","AIC","AICC","BIC","converged","iterations","df.null","df.model","df.residual","null.deviance","model.deviance"))
expect_identical(names(ao_mle$receiver_model),c("coefficients","loglik","gradient","hessian","vcov","se","residual.deviance","AIC","AICC","BIC","converged","iterations","df.null","df.model","df.residual","null.deviance","model.deviance"))
expect_length(attributes(ao_mle),9)
expect_identical(names(attributes(ao_mle)),c("names","class","model","ordinal","method","approach","formula","statistics","ncores"))
expect_identical(attr(ao_mle,"approach"),"Frequentist")
expect_silent(print(ao_mle))
expect_silent(summary(ao_mle))
#expect_silent(diagnostics(object = ao_mle, reh = ao_reh, stats = ao_reh_stats))
#ao_reh_diagnostics <- diagnostics(object = ao_mle, reh = ao_reh, stats = ao_reh_stats)
#expect_silent(plot(x = ao_mle,reh = ao_reh, diagnostics = ao_reh_diagnostics))
#expect_silent(plot(x = ao_mle,reh = ao_reh, diagnostics = ao_reh_diagnostics, which = 2, sender_effects = "indegreeSender", receiver_effects = "inertia")) # plotting specific effects from both models
#expect_silent(plot(x = ao_mle,reh = ao_reh, diagnostics = ao_reh_diagnostics, which = 2, sender_effects = NA, receiver_effects = "inertia")) # plotting specific effect from one model only

expect_silent(AIC(ao_mle))
expect_silent(AICC(ao_mle))
expect_silent(BIC(ao_mle))
expect_error(WAIC(ao_mle))

# # error when 'diagnostics' has at least one statistic's name wrong
# colnames(ao_reh_diagnostics$sender_model$residuals$smoothing_weights)[1] <- "INDEGREEsender"
# expect_error(plot(x = ao_mle,reh = ao_reh, diagnostics = ao_reh_diagnostics),
# "one or more effects not found inside the object 'diagnostics'.",
# fixed=TRUE)

# test diagnostics with only baseline
# ao_stats_baseline <- remstats(reh = ao_reh, sender_effects = ~1 , method = "pt")
# ao_mle_only_baseline <- remstimate(reh = ao_reh,
#                         stats = ao_stats_baseline,
#                         ncores = 1L,
#                         method = "MLE")
# expect_silent(diagnostics(object = ao_mle_only_baseline, reh = ao_reh, stats = ao_stats_baseline))



# tests on "WAIC" for "MLE"

# WAIC = TRUE
expect_silent(remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE))
# nsimWAIC is not a numeric scalar
expect_warning(remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = "text"))
# nsimWAIC is supplied
expect_silent(remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100))
ao_mle_with_waic <- remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100)
expect_silent(print(ao_mle_with_waic))
expect_silent(summary(ao_mle_with_waic))
expect_silent(WAIC(ao_mle_with_waic))

# WAIC for ordinal likelihod + testing omit_dyad routine for ordinal likelihood and tie-oriented model
# omit_dyad <- list()
# omit_dyad[[1]] <- list(time = c(330,585), dyad = data.frame(actor1=c(NA,"4"),actor2=c("4",NA),type=c(NA,NA))) # excluding actor 4 from the risk set between (observed) time points 330 and 585
ao_reh_ordinal <- remify(edgelist = ao_data$edgelist, model = "actor", ordinal = TRUE)
ao_reh_stats_ordinal <- remstats(reh = ao_reh_ordinal, sender_effects = rate_model, receiver_effects = choice_model)
expect_silent(remstimate(reh = ao_reh_ordinal,
                        stats = ao_reh_stats_ordinal,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100))
#testing only ordinal likelihood for actor-oriented model
ao_mle_ordinal_omit_dyad <- remstimate(reh = ao_reh_ordinal,
                        stats = ao_reh_stats_ordinal,
                        ncores = 1L)
# WAIC for interval likelihood with omit_dyad
ao_reh_omit <- remify(edgelist = ao_data$edgelist, model = "actor")
ao_reh_stats_omit <- remstats(reh = ao_reh_omit, sender_effects = rate_model, receiver_effects = choice_model)
expect_silent(remstimate(reh = ao_reh_omit,
                        stats = ao_reh_stats_omit,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100))
expect_silent(diagnostics(object = ao_mle_ordinal_omit_dyad,reh = ao_reh_ordinal,stats = ao_reh_stats_ordinal))

# (4) method  = "HMC"
expect_silent(remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 100L,
                        burnin = 5L,
                        seed = set_seed))
ao_hmc <- remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L,
                        L = 50,
                        seed = set_seed)
ao_hmc_two_chains <- remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 2L,
                        nsim = 10L,
                        burnin = 5L,
                        seed = set_seed)
expect_silent(ao_hmc)
expect_inherits(ao_hmc,"remstimate")
expect_length(ao_hmc,2)
expect_identical(names(ao_hmc),c("sender_model","receiver_model"))
expect_length(ao_hmc$sender_model,9)
expect_length(ao_hmc$receiver_model,9)
expect_identical(names(ao_hmc$sender_model),c("coefficients","post.mean","vcov","sd","loglik","draws","df.null","df.model","df.residual"))
expect_identical(names(ao_hmc$receiver_model),c("coefficients","post.mean","vcov","sd","loglik","draws","df.null","df.model","df.residual"))
expect_length(attributes(ao_hmc),16)
expect_identical(names(attributes(ao_hmc)),c("names","class","model","ordinal","method","approach","formula","statistics","ncores","nsim","nchains","burnin","thin","L","epsilon","seed"))
expect_identical(attr(ao_hmc,"approach"),"Bayesian")
expect_silent(print(ao_hmc))
expect_silent(summary(ao_hmc))
expect_silent(diagnostics(object = ao_hmc, reh = ao_reh, stats = ao_reh_stats))
ao_reh_diagnostics <- diagnostics(object = ao_hmc, reh = ao_reh, stats = ao_reh_stats)
ao_reh_diagnostics <- diagnostics(object = ao_hmc_two_chains, reh = ao_reh, stats = ao_reh_stats) # two chains
expect_silent(plot(x = ao_hmc_two_chains,reh = ao_reh, diagnostics = ao_reh_diagnostics)) # two chains
expect_error(BIC(ao_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(AICC(ao_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(BIC(ao_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(WAIC(ao_hmc))

# tests on "WAIC" for "HMC"
expect_silent(remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 100L,
                        burnin = 5L,
                        seed = set_seed,
                        WAIC = TRUE,
                        nsimWAIC = 10))

# ordinal likelihood (actor-oriented modeling)
ao_reh <- remify(edgelist = ao_data$edgelist, model = "actor", ordinal = TRUE)
ao_reh_stats <- remstats(reh = ao_reh, sender_effects = rate_model, receiver_effects = choice_model)
expect_silent(remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "MLE"))


# Risk set "active"

# testing estimation methods with active riskset
ao_reh <- remify(edgelist = ao_data$edgelist, model = "actor", riskset = "active")

# calculating statistics
ao_reh_stats <- remstats(reh = ao_reh, sender_effects = rate_model, receiver_effects = choice_model)
attr(ao_reh_stats,"formula") <- NULL

# (1) method = "MLE"
rate_model_two_effects <- ~ 1 + indegreeSender() + outdegreeSender()
ao_reh_stats_2 <- remstats(reh = ao_reh, sender_effects = rate_model_two_effects, receiver_effects = choice_model)
expect_silent(remstimate(reh = ao_reh,
                        stats = ao_reh_stats_2,
                        ncores = 1L,
                        method = "MLE"))

ao_mle <- remstimate(reh = ao_reh,
                        stats = ao_reh_stats_2,
                        ncores = 1L,
                        method = "MLE")
expect_silent(diagnostics(object = ao_mle, reh = ao_reh, stats = ao_reh_stats_2))

# (4) method  = "HMC"
expect_silent(remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L,
                        seed = set_seed))



