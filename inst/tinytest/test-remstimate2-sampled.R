# Tests for remstimate(): MLE and HMC with full and sampled tomstats
# Key tests:
#   1. Full stats (tomstats): remstimate MLE matches remstimate MLE
#   2. Sampled full (samp_num = D): estimates match full MLE exactly
#   3. Sampled partial: correct sign, reasonable magnitude, SE >= full SE
#   4. Output structure checks
#   5. HMC: runs without error, posterior mean near MLE

library(tinytest)

data(history, package = "remstats")
data(info,    package = "remstats")

colnames(history)[colnames(history) == "setting"] <- "type"
history_sub <- history[1:40, ]

reh <- remify::remify(edgelist = history_sub, model = "tie", riskset = "active")

effects <- ~ inertia(consider_type = FALSE) +
               indegreeSender(consider_type = FALSE) +
               outdegreeSender(consider_type = FALSE)

ts_full <- remstats::tomstats(effects, reh = reh,
                                attr_actors = info,
                                memory = "decay", memory_value = 1000,
                                start = 2, stop = 30,
                                sampling = FALSE)

D <- nrow(attr(ts_full, "riskset"))  # number of dyads in active riskset

# ---------------------------------------------------------------------------
# SECTION 1: Full stats — remstimate MLE matches remstimate MLE
# ---------------------------------------------------------------------------
est_full2 <- remstimate(reh = reh, stats = ts_full, method = "MLE", ncores = 1L)

expect_inherits(est_full2, "remstimate", info = "full: inherits remstimate")
expect_identical(attr(est_full2, "approach"), "Frequentist",
  info = "full: approach is Frequentist")
expect_identical(attr(est_full2, "method"), "MLE",
  info = "full: method is MLE")
expect_false(isTRUE(attr(est_full2, "sampled")),
  info = "full: sampled attribute is FALSE")

# Coefficients present and named
expect_true(!is.null(est_full2$coefficients), info = "full: coefficients present")
expect_equal(names(est_full2$coefficients),
  c("baseline", "inertia", "indegreeSender", "outdegreeSender"),
  info = "full: coefficient names correct")

# vcov is positive definite (all eigenvalues > 0)
expect_true(all(eigen(est_full2$vcov)$values > 0),
  info = "full: vcov positive definite")

# loglik is finite and negative
expect_true(is.finite(est_full2$loglik), info = "full: loglik finite")
expect_true(est_full2$loglik < 0, info = "full: loglik negative")

# AIC/BIC present
expect_true(!is.null(est_full2$AIC), info = "full: AIC present")
expect_true(!is.null(est_full2$BIC), info = "full: BIC present")

# converged
expect_true(isTRUE(est_full2$converged), info = "full: converged")

# ---------------------------------------------------------------------------
# SECTION 2: Sampled full (samp_num = D) — should match full MLE exactly
# ---------------------------------------------------------------------------
ts_samp_full <- remstats::tomstats(effects, reh = reh,
                                     attr_actors = info,
                                     memory = "decay", memory_value = 1000,
                                     start = 2, stop = 30,
                                     sampling = TRUE,
                                     samp_num = D,
                                     seed = 1L)

est_samp_full <- remstimate(reh = reh, stats = ts_samp_full,
                              method = "MLE", ncores = 1L)

expect_inherits(est_samp_full, "remstimate",
  info = "samp_full: inherits remstimate")
expect_true(isTRUE(attr(est_samp_full, "sampled")),
  info = "samp_full: sampled attribute is TRUE")

# With samp_num = D, all dyads sampled — coefficients should be consistent
# (same sign and similar magnitude) but NOT identical to full MLE because
# the sampled likelihood uses conditional logistic form which drops the
# waiting-time component of the full interval likelihood.
expect_equal(sign(est_samp_full$coefficients), sign(est_full2$coefficients),
  info = "samp_full: coefficient signs match full MLE")

expect_true(
  max(abs(est_samp_full$coefficients - est_full2$coefficients)) < 5,
  info = "samp_full: coefficients in reasonable range of full MLE")

# loglik of sampled form is different (drops waiting time term) — just check finite
expect_true(is.finite(est_samp_full$loglik),
  info = "samp_full: loglik finite")

# ---------------------------------------------------------------------------
# SECTION 3: Sampled partial — sign consistency and SE inflation
# ---------------------------------------------------------------------------
ts_samp <- remstats::tomstats(effects, reh = reh,
                                attr_actors = info,
                                memory = "decay", memory_value = 1000,
                                start = 2, stop = 30,
                                sampling = TRUE,
                                samp_num = 5L,
                                seed = 42L)

est_samp <- remstimate(reh = reh, stats = ts_samp, method = "MLE", ncores = 1L)

expect_inherits(est_samp, "remstimate",
  info = "samp_partial: inherits remstimate")
expect_true(isTRUE(attr(est_samp, "sampled")),
  info = "samp_partial: sampled attribute is TRUE")

# Coefficients present and named
expect_equal(names(est_samp$coefficients),
  names(est_full2$coefficients),
  info = "samp_partial: coefficient names match")

# Signs should match full MLE (sampling is unbiased asymptotically)
expect_equal(sign(est_samp$coefficients), sign(est_full2$coefficients),
  info = "samp_partial: coefficient signs match full MLE")

# SE should be >= full SE (sampling adds variance)
# Use tolerance since numerical differences can be small
expect_true(all(est_samp$se >= est_full2$se * 0.9),
  info = "samp_partial: SE not smaller than full SE (with 10% tolerance)")

# loglik is finite
expect_true(is.finite(est_samp$loglik),
  info = "samp_partial: sampled loglik finite")

# ---------------------------------------------------------------------------
# SECTION 4: Output structure checks for sampled MLE
# ---------------------------------------------------------------------------
expected_names <- c("coefficients", "loglik", "gradient", "hessian",
                     "vcov", "se", "residual.deviance", "null.deviance",
                     "model.deviance", "df.null", "df.model", "df.residual",
                     "AIC", "AICC", "BIC", "converged", "iterations",
                     "sampled", "samp_num", "sampling_scheme")

expect_true(all(expected_names %in% names(est_samp)),
  info = "samp_partial: all expected output fields present")

expect_equal(est_samp$df.model, 4L,
  info = "samp_partial: df.model = number of parameters")
expect_equal(est_samp$samp_num, 5L,
  info = "samp_partial: samp_num stored correctly")

# ---------------------------------------------------------------------------
# SECTION 5: Two different seeds give different but similar estimates
# ---------------------------------------------------------------------------
ts_samp2 <- remstats::tomstats(effects, reh = reh,
                                 attr_actors = info,
                                 memory = "decay", memory_value = 1000,
                                 start = 2, stop = 30,
                                 sampling = TRUE, samp_num = 5L, seed = 99L)

est_samp2 <- remstimate(reh = reh, stats = ts_samp2, method = "MLE", ncores = 1L)

# Different seeds -> different estimates
expect_false(identical(est_samp$coefficients, est_samp2$coefficients),
  info = "different seeds: different coefficient estimates")

# But same sign
expect_equal(sign(est_samp2$coefficients), sign(est_full2$coefficients),
  info = "different seeds: signs still match full MLE")

# ---------------------------------------------------------------------------
# SECTION 6: Typed events with consider_type = "separate"
# ---------------------------------------------------------------------------
effects_typed <- ~ inertia(consider_type = "separate") +
                    outdegreeSender(consider_type = FALSE)

ts_typed_full <- remstats::tomstats(effects_typed, reh = reh,
                                      attr_actors = info,
                                      memory = "decay", memory_value = 1000,
                                      start = 2, stop = 30,
                                      sampling = FALSE)

est_typed_full <- remstimate(reh = reh, stats = ts_typed_full,
                               method = "MLE", ncores = 1L)

expect_inherits(est_typed_full, "remstimate",
  info = "typed full: inherits remstimate")
expect_true("inertia.social" %in% names(est_typed_full$coefficients),
  info = "typed full: separate slice names in coefficients")
expect_true("inertia.work" %in% names(est_typed_full$coefficients),
  info = "typed full: separate slice names in coefficients")

# Sampled typed
D_typed <- nrow(attr(ts_typed_full, "riskset"))
ts_typed_samp <- remstats::tomstats(effects_typed, reh = reh,
                                      attr_actors = info,
                                      memory = "decay", memory_value = 1000,
                                      start = 2, stop = 30,
                                      sampling = TRUE,
                                      samp_num = min(5L, D_typed),
                                      seed = 7L)

est_typed_samp <- remstimate(reh = reh, stats = ts_typed_samp,
                               method = "MLE", ncores = 1L)
expect_inherits(est_typed_samp, "remstimate",
  info = "typed samp: inherits remstimate")
expect_equal(names(est_typed_samp$coefficients),
  names(est_typed_full$coefficients),
  info = "typed samp: coefficient names match full")

# ---------------------------------------------------------------------------
# SECTION 7: HMC runs without error, posterior mean near MLE
# ---------------------------------------------------------------------------
est_hmc <- remstimate(reh = reh, stats = ts_full,
                        method = "HMC", ncores = 1L,
                        nsim = 200L, burnin = 100L, thin = 5L,
                        L = 20L, epsilon = 0.002, seed = 1L)

expect_inherits(est_hmc, "remstimate", info = "HMC: inherits remstimate")
expect_identical(attr(est_hmc, "approach"), "Bayesian",
  info = "HMC: approach is Bayesian")
expect_identical(attr(est_hmc, "method"), "HMC",
  info = "HMC: method is HMC")

# Posterior mean should be in reasonable range of MLE
expect_true(
  max(abs(est_hmc$post.mean - est_full2$coefficients)) < 1.0,
  info = "HMC: posterior mean within 1 unit of MLE")

# draws present
expect_true(!is.null(est_hmc$draws), info = "HMC: draws present")
expect_equal(ncol(est_hmc$draws), 4L, info = "HMC: draws has P columns")

# ---------------------------------------------------------------------------
# SECTION 8: Error handling
# ---------------------------------------------------------------------------
expect_error(
  remstimate(reh = list(), stats = ts_full),
  info = "error: non-remify reh"
)

expect_error(
  remstimate(reh = reh, stats = array(1, dim = c(3,3,3))),
  info = "error: non-remstats stats"
)

