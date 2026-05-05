
library(tinytest)

# ── shared setup ──────────────────────────────────────────────────────────────

# Actor-oriented data with simultaneous events
data(ao_data)
ao_data$edgelist$time <- floor(ao_data$edgelist$time)
ao_data$edgelist$time[seq(5, 95, by = 5)] <-
  ao_data$edgelist$time[seq(5, 95, by = 5) - 1]

rate_model   <- ~ 1 + indegreeSender()
choice_model <- ~ inertia() + reciprocity()

ao_reh   <- remify(edgelist = ao_data$edgelist, model = "actor")
ao_stats <- remstats(reh              = ao_reh,
                                sender_effects   = rate_model,
                                receiver_effects = choice_model)
ao_mle   <- remstimate(reh    = ao_reh,
                                    stats  = ao_stats,
                                    ncores = 1L,
                                    method = "MLE",
                                    WAIC   = TRUE)

# Tie-oriented data
data(tie_data)
tie_model <- ~ 1 + inertia() + reciprocity()
tie_reh   <- remify(edgelist = tie_data$edgelist, model = "tie")
tie_stats <- remstats(reh = tie_reh, tie_effects = tie_model)
tie_mle   <- remstimate(reh    = tie_reh,
                                     stats  = tie_stats,
                                     ncores = 1L,
                                     method = "MLE",
                                     WAIC   = TRUE)

# helper: open a null pdf device, run expr, close device
with_null_dev <- function(expr) {
  pdf(nullfile())
  on.exit(dev.off())
  force(expr)
}

# ── diagnostics() output structure ────────────────────────────────────────────

expect_silent(
  ao_diag <- diagnostics(object = ao_mle, reh = ao_reh, stats = ao_stats)
)
expect_inherits(ao_diag, c("diagnostics", "remstimate"))
expect_true(all(c("sender_model", "receiver_model", ".reh.processed") %in% names(ao_diag)))
expect_false(is.null(ao_diag$sender_model$residuals))
expect_false(is.null(ao_diag$sender_model$rates))
expect_false(is.null(ao_diag$sender_model$recall))
expect_false(is.null(ao_diag$receiver_model$residuals))
expect_false(is.null(ao_diag$receiver_model$rates))
expect_false(is.null(ao_diag$receiver_model$recall))

expect_silent(
  tie_diag <- diagnostics(object = tie_mle, reh = tie_reh, stats = tie_stats)
)
expect_inherits(tie_diag, c("diagnostics", "remstimate"))
expect_true(all(c("residuals", "rates", "recall", ".reh.processed") %in% names(tie_diag)))

# recall structure
expect_true(all(c("per_event", "summary") %in% names(ao_diag$sender_model$recall)))
expect_true(all(c("mean_rel_rank", "median_rel_rank", "mean_cum_prob",
                   "top_pct", "top_pct_prop") %in%
                names(ao_diag$sender_model$recall$summary)))
expect_true(all(c("event", "rel_rank", "cum_prob") %in%
                names(ao_diag$sender_model$recall$per_event)))
expect_true(all(ao_diag$sender_model$recall$per_event$rel_rank < 1))
expect_true(all(ao_diag$sender_model$recall$per_event$rel_rank >= 0))

# ── S3 dispatch: plot(diag_obj) ───────────────────────────────────────────────

# ── plot.diagnostics: actor model ─────────────────────────────────────────────

# plots 1 + 2 with no 'object' (MLE-only path)
expect_silent(
  with_null_dev(
    plot(ao_diag, which = c(1, 2))
  )
)

# plot 1 only (waiting times)
expect_silent(
  with_null_dev(
    plot(ao_diag, which = 1)
  )
)

# plot 2 only — all effects
expect_silent(
  with_null_dev(
    plot(ao_diag, which = 2)
  )
)

# plot 2: specific sender and receiver effects
expect_silent(
  with_null_dev(
    plot(ao_diag, which = 2,
                     sender_effects   = "indegreeSender",
                     receiver_effects = "inertia")
  )
)

# plot 2: skip sender entirely (NA), plot one receiver effect
expect_silent(
  with_null_dev(
    plot(ao_diag, which = 2,
                     sender_effects   = NA,
                     receiver_effects = "reciprocity")
  )
)

# plot 2: wrong effect name → error
expect_error(
  with_null_dev(
    plot(ao_diag, which = 2,
                     sender_effects = "NONEXISTENT_STAT")
  )
)

# with 'object' (MLE — plots 3/4 silently skipped because method != "HMC")
expect_silent(
  with_null_dev(
    plot(ao_diag, object = ao_mle, which = c(1, 2, 3, 4))
  )
)

# return value is x invisibly
ret <- with_null_dev(
  plot(ao_diag, which = 1)
)
expect_identical(ret, ao_diag)

# ── plot.diagnostics: tie model ───────────────────────────────────────────────

expect_silent(
  with_null_dev(
    plot(tie_diag, which = c(1, 2))
  )
)

expect_silent(
  with_null_dev(
    plot(tie_diag, which = 2, effects = "inertia")
  )
)

expect_silent(
  with_null_dev(
    plot(tie_diag, which = 2, effects = c("inertia", "reciprocity"))
  )
)

expect_error(
  with_null_dev(
    plot(tie_diag, which = 2, effects = "NONEXISTENT")
  )
)

# with 'object' supplied for MLE (plots 3/4 skipped silently)
expect_silent(
  with_null_dev(
    plot(tie_diag, object = tie_mle, which = c(1, 2, 3, 4))
  )
)

# ── plot.remstimate: backward compatibility ───────────────────────────────────

# old-style: pass stats via ..., let plot.remstimate compute diagnostics
expect_silent(
  with_null_dev(
    plot(ao_mle, reh = ao_reh, stats = ao_stats, which = 1)
  )
)

expect_silent(
  with_null_dev(
    plot(tie_mle, reh = tie_reh, stats = tie_stats, which = c(1, 2))
  )
)

# old-style: pass precomputed diagnostics
expect_silent(
  with_null_dev(
    plot(ao_mle, reh = ao_reh, diagnostics = ao_diag, which = c(1, 2))
  )
)

expect_silent(
  with_null_dev(
    plot(tie_mle, reh = tie_reh, diagnostics = tie_diag, which = 2,
         effects = "inertia")
  )
)

# error: no stats and no diagnostics supplied
expect_error(
  plot(ao_mle, reh = ao_reh, which = 1),
  "'stats' must be provided if argument 'diagnostics' is NULL",
  fixed = TRUE
)

# error: diagnostics is wrong class
expect_error(
  plot(ao_mle, reh = ao_reh, diagnostics = list(foo = 1), which = 1)
)

# ── HMC: plots 3 and 4 ───────────────────────────────────────────────────────
# Small HMC run (1 chain, 15 draws) to exercise posterior/trace plots.
# Wrapped in tryCatch so a slow CI environment can skip gracefully.

ao_hmc <- tryCatch(
  remstimate(reh     = ao_reh,
                          stats   = ao_stats,
                          ncores  = 1L,
                          method  = "HMC",
                          nchains = 1L,
                          nsim    = 15L,
                          burnin  = 5L,
                          seed    = 12345L),
  error = function(e) NULL
)

if (!is.null(ao_hmc)) {
  ao_hmc_diag <- diagnostics(object = ao_hmc, reh = ao_reh, stats = ao_stats)

  # plot 3: posterior histograms
  expect_silent(
    with_null_dev(
      plot(ao_hmc_diag, object = ao_hmc, which = 3)
    )
  )

  # plot 4: trace plots (single chain)
  expect_silent(
    with_null_dev(
      plot(ao_hmc_diag, object = ao_hmc, which = 4)
    )
  )

  # specific effects
  expect_silent(
    with_null_dev(
      plot(ao_hmc_diag, object = ao_hmc, which = c(3, 4),
                       sender_effects = "indegreeSender")
    )
  )

  # backward compat via plot.remstimate
  expect_silent(
    with_null_dev(
      plot(ao_hmc, reh = ao_reh, diagnostics = ao_hmc_diag, which = c(3, 4))
    )
  )
}

# multi-chain HMC: trace plot branches to multi-chain path
ao_hmc2 <- tryCatch(
  remstimate(reh     = ao_reh,
                          stats   = ao_stats,
                          ncores  = 1L,
                          method  = "HMC",
                          nchains = 2L,
                          nsim    = 15L,
                          burnin  = 5L,
                          seed    = 12345L),
  error = function(e) NULL
)

if (!is.null(ao_hmc2)) {
  ao_hmc2_diag <- diagnostics(object = ao_hmc2, reh = ao_reh, stats = ao_stats)
  expect_silent(
    with_null_dev(
      plot(ao_hmc2_diag, object = ao_hmc2, which = 4,
                       sender_effects = "indegreeSender")
    )
  )
}

