
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

ao_mle_active <- remstimate::remstimate(reh    = ao_reh_active,
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

# ── consistency check: all riskset variants produce finite log-likelihoods ────
expect_true(is.finite(ao_mle_active$sender_model$loglik))
expect_true(is.finite(ao_mle_manual$sender_model$loglik))
expect_true(is.finite(ao_mle_sat$sender_model$loglik))

# ── validation: compare active riskset MLE with clogit / Poisson GLM ─────────
#
# Receiver choice model  → conditional logistic regression (survival::clogit)
#   one stratum per event; case = observed receiver, controls = other valid receivers
#
# Sender rate model      → Poisson GLM with offset(log(interevent_time))
#   one row per (actor, event) for actors in sender_riskset

if (requireNamespace("survival", quietly = TRUE)) {
  library(survival)

  M   <- nrow(ao_data$edgelist)
  iet <- ao_reh_active$intereventTime            # interevent times
  a1  <- ao_data$edgelist$actor1                 # 1-indexed observed senders
  a2  <- ao_data$edgelist$actor2                 # 1-indexed observed receivers

  sender_rs   <- ao_reh_active$sender_riskset    # 1-indexed vector of valid senders
  receiver_rs <- ao_reh_active$receiver_riskset  # list[N]: valid receivers per sender (1-indexed)

  # Note: with ao_data all 5 actors appear as both senders and receivers, so the
  # active riskset is effectively full here. The validation still confirms that
  # the index-based code path produces estimates identical to the reference methods.

  # sender stats:   array[event, actor, param],  params = (baseline, indegreeSender)
  # receiver stats: array[event, actor, param],  params = (inertia, reciprocity)
  s_stats <- ao_stats_active$sender_stats
  r_stats <- ao_stats_active$receiver_stats
  start_stop1 <- as.numeric(attr(ao_stats_active,"subset"))
  length1 <- dim(r_stats)[1]

  # ── receiver choice: clogit ───────────────────────────────────────────────

  recv_rows <- vector("list", length1)
  for (m in seq_len(length1)) {
    s   <- a1[m+start_stop1[1]-1]
    rec <- a2[m+start_stop1[1]-1]
    rs  <- receiver_rs[[s]]          # valid receivers for this sender (1-indexed)
    recv_rows[[m]] <- data.frame(
      stratum     = m,
      y           = as.integer(rs == rec),
      inertia     = r_stats[m, rs, 1],
      reciprocity = r_stats[m, rs, 2]
    )
  }
  df_recv <- do.call(rbind, recv_rows)

  fit_clogit <- survival::clogit(
    y ~ inertia + reciprocity + survival::strata(stratum),
    data = df_recv
  )

  expect_equal(
    unname(coef(ao_mle_active$receiver_model)[c("inertia", "reciprocity")]),
    unname(coef(fit_clogit)),
    tolerance = 1e-4
  )

  # ── sender rate: Poisson GLM ──────────────────────────────────────────────

  send_rows <- vector("list", length1)
  for (m in seq_len(length1)) {
    s <- a1[m+start_stop1[1]-1]
    send_rows[[m]] <- data.frame(
      y              = as.integer(sender_rs == s),
      log_iet        = log(iet[m+start_stop1[1]-1]),
      indegreeSender = s_stats[m, sender_rs, 2]   # col 1 = baseline (all 1s)
    )
  }
  df_send <- do.call(rbind, send_rows)

  fit_poisson <- glm(
    y ~ indegreeSender + offset(log_iet),
    data   = df_send,
    family = poisson
  )

  # (Intercept) in the Poisson model corresponds to the baseline coefficient
  expect_equal(
    unname(coef(ao_mle_active$sender_model)[c("baseline", "indegreeSender")]),
    unname(coef(fit_poisson)[c("(Intercept)", "indegreeSender")]),
    tolerance = 1e-4
  )

}
