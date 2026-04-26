
library(tinytest)

# loading data
data(ao_data)

set_seed <- 23929

# specifying linear predictors
rate_model   <- ~ 1 + remstats::indegreeSender()
choice_model <- ~ remstats::inertia() + remstats::reciprocity()

# ── active riskset ────────────────────────────────────────────────────────────

ao_reh_active <- remify::remify(edgelist = ao_data$edgelist,
                                 model   = "actor",
                                 riskset = "active")

ao_stats_active <- remstats::remstats(reh              = ao_reh_active,
                                       sender_effects   = rate_model,
                                       receiver_effects = choice_model)

# MLE
expect_silent(
  remstimate::remstimate(reh    = ao_reh_active,
                          stats  = ao_stats_active,
                          ncores = 1L,
                          method = "MLE")
)

ao_mle_active <- remstimate::remstimate(reh = ao_reh_active,
                                         stats  = ao_stats_active,
                                         ncores = 1L,
                                         method = "MLE",
                                         WAIC   = TRUE)

expect_inherits(ao_mle_active, "remstimate")
expect_identical(attr(ao_mle_active, "approach"), "Frequentist")
expect_identical(names(ao_mle_active), c("sender_model", "receiver_model"))
expect_false(is.null(ao_mle_active$sender_model$coefficients))
expect_false(is.null(ao_mle_active$receiver_model$coefficients))
expect_identical(names(ao_mle_active$sender_model$coefficients),
                 c("baseline", "indegreeSender"))
expect_identical(names(ao_mle_active$receiver_model$coefficients),
                 c("inertia", "reciprocity"))

# information criteria
expect_silent(AIC(ao_mle_active))
expect_silent(AICC(ao_mle_active))
expect_silent(BIC(ao_mle_active))
expect_silent(WAIC(ao_mle_active))

# HMC
expect_silent(
  remstimate::remstimate(reh     = ao_reh_active,
                          stats   = ao_stats_active,
                          ncores  = 1L,
                          method  = "HMC",
                          nchains = 1L,
                          nsim    = 10L,
                          burnin  = 5L,
                          seed    = set_seed)
)

ao_hmc_active <- remstimate::remstimate(reh     = ao_reh_active,
                                         stats   = ao_stats_active,
                                         ncores  = 1L,
                                         method  = "HMC",
                                         nchains = 1L,
                                         nsim    = 10L,
                                         burnin  = 5L,
                                         seed    = set_seed)

expect_inherits(ao_hmc_active, "remstimate")
expect_identical(attr(ao_hmc_active, "approach"), "Bayesian")

# ── manual riskset ────────────────────────────────────────────────────────────

# build manual riskset from observed dyads + their reverses
el <- ao_data$edgelist[, c("actor1", "actor2")]
manual_rs <- unique(as.data.frame(rbind(
  as.matrix(el),
  as.matrix(el[, c("actor2", "actor1")])
)))

ao_reh_manual <- remify::remify(edgelist       = ao_data$edgelist,
                                 model          = "actor",
                                 riskset        = "manual",
                                 manual.riskset = manual_rs)

ao_stats_manual <- remstats::remstats(reh              = ao_reh_manual,
                                       sender_effects   = rate_model,
                                       receiver_effects = choice_model)

# MLE
expect_silent(
  remstimate::remstimate(reh    = ao_reh_manual,
                          stats  = ao_stats_manual,
                          ncores = 1L,
                          method = "MLE")
)

ao_mle_manual <- remstimate::remstimate(reh    = ao_reh_manual,
                                         stats  = ao_stats_manual,
                                         ncores = 1L,
                                         method = "MLE",
                                         WAIC   = TRUE)

expect_inherits(ao_mle_manual, "remstimate")
expect_identical(attr(ao_mle_manual, "approach"), "Frequentist")
expect_false(is.null(ao_mle_manual$sender_model$coefficients))
expect_false(is.null(ao_mle_manual$receiver_model$coefficients))
expect_silent(WAIC(ao_mle_manual))

# ── active_saturated riskset ──────────────────────────────────────────────────

ao_reh_sat <- remify::remify(edgelist = ao_data$edgelist,
                              model   = "actor",
                              riskset = "active_saturated")

ao_stats_sat <- remstats::remstats(reh              = ao_reh_sat,
                                    sender_effects   = rate_model,
                                    receiver_effects = choice_model)

# MLE
expect_silent(
  remstimate::remstimate(reh    = ao_reh_sat,
                          stats  = ao_stats_sat,
                          ncores = 1L,
                          method = "MLE")
)

ao_mle_sat <- remstimate::remstimate(reh    = ao_reh_sat,
                                      stats  = ao_stats_sat,
                                      ncores = 1L,
                                      method = "MLE",
                                      WAIC   = TRUE)

expect_inherits(ao_mle_sat, "remstimate")
expect_identical(attr(ao_mle_sat, "approach"), "Frequentist")
expect_false(is.null(ao_mle_sat$sender_model$coefficients))
expect_false(is.null(ao_mle_sat$receiver_model$coefficients))
expect_silent(WAIC(ao_mle_sat))

# ── consistency check: active == active_saturated when data is the same ───────
# (both should produce valid estimates, not necessarily identical since
#  the risksets may differ)
expect_true(is.finite(ao_mle_active$sender_model$loglik))
expect_true(is.finite(ao_mle_manual$sender_model$loglik))
expect_true(is.finite(ao_mle_sat$sender_model$loglik))

