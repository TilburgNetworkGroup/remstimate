# Validation tests: remstimate2 MLE vs glm/clogit on stacked data
#
# For interval timing:
#   glm(obs ~ -1 + baseline + inertia + ... + offset(log_interevent),
#       family = poisson) should match remstimate2 MLE
#
# For ordinal timing:
#   survival::clogit(obs ~ -1 + inertia + ... + strata(event))
#   should match remstimate2 MLE (no baseline in ordinal)

library(tinytest)

data(history, package = "remstats")
data(info,    package = "remstats")
colnames(history)[colnames(history) == "setting"] <- "type"
history_sub <- history[1:40, ]

tol <- 1e-4

# ---------------------------------------------------------------------------
# SECTION 1: Interval likelihood — remstimate2 MLE vs Poisson GLM
# ---------------------------------------------------------------------------
reh_int <- remify::remify2(edgelist = history_sub, model = "tie",
                            riskset = "active")

effects <- ~ inertia(consider_type = FALSE) +
               indegreeSender(consider_type = FALSE) +
               outdegreeSender(consider_type = FALSE)

ts_int <- remstats::tomstats2(effects, reh = reh_int,
                               attr_actors = info,
                               memory = "decay", memory_value = 1000,
                               start = 2, stop = 30,
                               sampling = FALSE)

est_int <- remstimate2(reh = reh_int, stats = ts_int,
                        method = "MLE", ncores = 1L)

stacked_int <- remstats::stack_stats(ts_int, reh_int)
df_int <- stacked_int$remstats_stack

# Poisson GLM with log inter-event time offset
fit_glm <- glm(obs ~ -1 + baseline + inertia + indegreeSender + outdegreeSender +
                 offset(log_interevent),
               family = poisson, data = df_int)

# Coefficients should match
expect_equal(
  unname(est_int$coefficients),
  unname(coef(fit_glm)),
  tolerance = tol,
  info = "interval: remstimate2 MLE matches Poisson GLM coefficients"
)
# Note: loglik differs by Poisson constant log(y!) — coefficients are the valid check

# ---------------------------------------------------------------------------
# SECTION 2: Ordinal likelihood — remstimate2 MLE vs survival::clogit
# ---------------------------------------------------------------------------
reh_ord <- remify::remify2(edgelist = history_sub, model = "tie",
                            riskset = "active", ordinal = TRUE)

effects_no_baseline <- ~ inertia(consider_type = FALSE) +
                           indegreeSender(consider_type = FALSE) +
                           outdegreeSender(consider_type = FALSE)

ts_ord <- remstats::tomstats2(effects_no_baseline, reh = reh_ord,
                               attr_actors = info,
                               memory = "decay", memory_value = 1000,
                               start = 2, stop = 30,
                               sampling = FALSE)

est_ord <- remstimate2(reh = reh_ord, stats = ts_ord,
                        method = "MLE", ncores = 1L)

stacked_ord <- remstats::stack_stats(ts_ord, reh_ord)
df_ord <- stacked_ord$remstats_stack

# Conditional logistic regression stratified by event
library(survival)
fit_clogit <- survival::clogit(
  obs ~ -1 + inertia + indegreeSender + outdegreeSender + survival::strata(event),
  data = df_ord
)

# Coefficients should match
expect_equal(
  unname(est_ord$coefficients),
  unname(coef(fit_clogit)),
  tolerance = tol,
  info = "ordinal: remstimate2 MLE matches clogit coefficients"
)

# ---------------------------------------------------------------------------
# SECTION 3: Interval — full riskset
# ---------------------------------------------------------------------------
reh_full <- remify::remify2(edgelist = history_sub, model = "tie",
                             riskset = "full")

ts_full_rs <- remstats::tomstats2(effects, reh = reh_full,
                                   attr_actors = info,
                                   memory = "decay", memory_value = 1000,
                                   start = 2, stop = 30,
                                   sampling = FALSE)

est_full_rs <- remstimate2(reh = reh_full, stats = ts_full_rs,
                            method = "MLE", ncores = 1L)

stacked_full_rs <- remstats::stack_stats(ts_full_rs, reh_full)
df_full_rs <- stacked_full_rs$remstats_stack

fit_glm_full <- glm(obs ~ -1 + baseline + inertia + indegreeSender +
                      outdegreeSender + offset(log_interevent),
                    family = poisson, data = df_full_rs)

expect_equal(
  unname(est_full_rs$coefficients),
  unname(coef(fit_glm_full)),
  tolerance = tol,
  info = "interval full riskset: remstimate2 MLE matches Poisson GLM"
)

# ---------------------------------------------------------------------------
# SECTION 4: Typed events — consider_type = "separate", interval
# ---------------------------------------------------------------------------
effects_typed <- ~ inertia(consider_type = "separate") +
                    outdegreeSender(consider_type = FALSE)

ts_typed <- remstats::tomstats2(effects_typed, reh = reh_int,
                                 attr_actors = info,
                                 memory = "decay", memory_value = 1000,
                                 start = 2, stop = 30,
                                 sampling = FALSE)

est_typed <- remstimate2(reh = reh_int, stats = ts_typed,
                          method = "MLE", ncores = 1L)

stacked_typed <- remstats::stack_stats(ts_typed, reh_int)
df_typed <- stacked_typed$remstats_stack

# Build formula dynamically from coefficient names
stat_cols <- names(est_typed$coefficients)
fmla_typed <- as.formula(paste(
  "obs ~ -1 +", paste(stat_cols, collapse = " + "), "+ offset(log_interevent)"
))

fit_glm_typed <- glm(fmla_typed, family = poisson, data = df_typed)

expect_equal(
  unname(est_typed$coefficients),
  unname(coef(fit_glm_typed)),
  tolerance = tol,
  info = "typed interval: remstimate2 MLE matches Poisson GLM (consider_type = 'separate')"
)

# Also check typed ordinal
ts_typed_ord <- remstats::tomstats2(effects_typed, reh = reh_ord,
                                     attr_actors = info,
                                     memory = "decay", memory_value = 1000,
                                     start = 2, stop = 30,
                                     sampling = FALSE)

est_typed_ord <- remstimate2(reh = reh_ord, stats = ts_typed_ord,
                              method = "MLE", ncores = 1L)

stacked_typed_ord <- remstats::stack_stats(ts_typed_ord, reh_ord)
df_typed_ord <- stacked_typed_ord$remstats_stack

stat_cols_ord <- setdiff(names(est_typed_ord$coefficients), "baseline")
fmla_typed_ord <- as.formula(paste(
  "obs ~ -1 +", paste(stat_cols_ord, collapse = " + "),
  "+ survival::strata(event)"
))

fit_clogit_typed <- survival::clogit(fmla_typed_ord, data = df_typed_ord)

expect_equal(
  unname(est_typed_ord$coefficients),
  unname(coef(fit_clogit_typed)),
  tolerance = tol,
  info = "typed ordinal: remstimate2 MLE matches clogit (consider_type = 'separate')"
)

# ---------------------------------------------------------------------------
# SECTION 5: Sampled stats — remstimate2 vs weighted Poisson GLM
# ---------------------------------------------------------------------------
samp_num <- 10L

ts_samp_int <- remstats::tomstats2(effects, reh = reh_int,
                                    attr_actors = info,
                                    memory = "decay", memory_value = 1000,
                                    start = 2, stop = 30,
                                    sampling = TRUE, samp_num = samp_num, seed = 1L)

est_samp_int <- remstimate2(reh = reh_int, stats = ts_samp_int,
                             method = "MLE", ncores = 1L)

stacked_samp_int <- remstats::stack_stats(ts_samp_int, reh_int)
df_samp_int <- stacked_samp_int$remstats_stack

fit_glm_samp <- glm(obs ~ -1 + baseline + inertia + indegreeSender +
                      outdegreeSender + offset(log_interevent),
                    weights = weight,
                    family = poisson, data = df_samp_int)

expect_equal(
  unname(est_samp_int$coefficients),
  unname(coef(fit_glm_samp)),
  tolerance = tol,
  info = "sampled interval: remstimate2 MLE matches weighted Poisson GLM"
)

# ---------------------------------------------------------------------------
# SECTION 6: Sampled stats — remstimate2 vs weighted clogit (ordinal)
# ---------------------------------------------------------------------------
ts_samp_ord <- remstats::tomstats2(effects_no_baseline, reh = reh_ord,
                                    attr_actors = info,
                                    memory = "decay", memory_value = 1000,
                                    start = 2, stop = 30,
                                    sampling = TRUE, samp_num = samp_num, seed = 1L)

est_samp_ord <- remstimate2(reh = reh_ord, stats = ts_samp_ord,
                             method = "MLE", ncores = 1L)

stacked_samp_ord <- remstats::stack_stats(ts_samp_ord, reh_ord)
df_samp_ord <- stacked_samp_ord$remstats_stack

fit_clogit_samp <- survival::clogit(
  obs ~ -1 + inertia + indegreeSender + outdegreeSender +
    survival::strata(event),
  weights = weight,
  method = "approximate",
  data = df_samp_ord
)

expect_equal(
  unname(est_samp_ord$coefficients),
  unname(coef(fit_clogit_samp)),
  tolerance = tol,
  info = "sampled ordinal: remstimate2 MLE matches weighted clogit"
)

