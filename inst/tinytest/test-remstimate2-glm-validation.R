# Validation tests: remstimate MLE vs glm/clogit on stacked data
#
# For interval timing:
#   glm(obs ~ -1 + baseline + inertia + ... + offset(log_interevent),
#       family = poisson) should match remstimate MLE
#
# For ordinal timing:
#   survival::clogit(obs ~ -1 + inertia + ... + strata(event))
#   should match remstimate MLE (no baseline in ordinal)

library(tinytest)

data(history, package = "remstats")
data(info,    package = "remstats")
colnames(history)[colnames(history) == "setting"] <- "type"
history_sub <- history[1:40, ]

tol <- 1e-4

# ---------------------------------------------------------------------------
# SECTION 1: Interval likelihood — remstimate MLE vs Poisson GLM
# ---------------------------------------------------------------------------
reh_int <- remify(edgelist = history_sub, model = "tie",
                            riskset = "active")

effects <- ~ inertia(consider_type = FALSE) +
               indegreeSender(consider_type = FALSE) +
               outdegreeSender(consider_type = FALSE)

ts_int <- tomstats(effects, reh = reh_int,
                               attr_actors = info,
                               memory = "decay", memory_value = 1000,
                               first = 2, last = 30,
                               sampling = FALSE)

est_int <- remstimate(reh = reh_int, stats = ts_int,
                        method = "MLE", ncores = 1L)

stacked_int <- stack_stats(ts_int, reh_int)
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
  info = "interval: remstimate MLE matches Poisson GLM coefficients"
)
# Note: loglik differs by Poisson constant log(y!) — coefficients are the valid check

# ---------------------------------------------------------------------------
# SECTION 2: Ordinal likelihood — remstimate MLE vs survival::clogit
# ---------------------------------------------------------------------------
reh_ord <- remify(edgelist = history_sub, model = "tie",
                            riskset = "active", ordinal = TRUE)

effects_no_baseline <- ~ inertia(consider_type = FALSE) +
                           indegreeSender(consider_type = FALSE) +
                           outdegreeSender(consider_type = FALSE)

ts_ord <- tomstats(effects_no_baseline, reh = reh_ord,
                               attr_actors = info,
                               memory = "decay", memory_value = 1000,
                               first = 2, last = 30,
                               sampling = FALSE)

est_ord <- remstimate(reh = reh_ord, stats = ts_ord,
                        method = "MLE", ncores = 1L)

stacked_ord <- stack_stats(ts_ord, reh_ord)
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
  info = "ordinal: remstimate MLE matches clogit coefficients"
)

# ---------------------------------------------------------------------------
# SECTION 3: Interval — full riskset
# ---------------------------------------------------------------------------
reh_full <- remify(edgelist = history_sub, model = "tie",
                             riskset = "full")

ts_full_rs <- tomstats(effects, reh = reh_full,
                                   attr_actors = info,
                                   memory = "decay", memory_value = 1000,
                                   first = 2, last = 30,
                                   sampling = FALSE)

est_full_rs <- remstimate(reh = reh_full, stats = ts_full_rs,
                            method = "MLE", ncores = 1L)

stacked_full_rs <- stack_stats(ts_full_rs, reh_full)
df_full_rs <- stacked_full_rs$remstats_stack

fit_glm_full <- glm(obs ~ -1 + baseline + inertia + indegreeSender +
                      outdegreeSender + offset(log_interevent),
                    family = poisson, data = df_full_rs)

expect_equal(
  unname(est_full_rs$coefficients),
  unname(coef(fit_glm_full)),
  tolerance = tol,
  info = "interval full riskset: remstimate MLE matches Poisson GLM"
)

# ---------------------------------------------------------------------------
# SECTION 4: Typed events — consider_type = "separate", interval
# ---------------------------------------------------------------------------
effects_typed <- ~ inertia(consider_type = "separate") +
                    outdegreeSender(consider_type = FALSE)

ts_typed <- tomstats(effects_typed, reh = reh_int,
                                 attr_actors = info,
                                 memory = "decay", memory_value = 1000,
                                 first = 2, last = 30,
                                 sampling = FALSE)

est_typed <- remstimate(reh = reh_int, stats = ts_typed,
                          method = "MLE", ncores = 1L)

stacked_typed <- stack_stats(ts_typed, reh_int)
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
  info = "typed interval: remstimate MLE matches Poisson GLM (consider_type = 'separate')"
)

# Also check typed ordinal
ts_typed_ord <- tomstats(effects_typed, reh = reh_ord,
                                     attr_actors = info,
                                     memory = "decay", memory_value = 1000,
                                     first = 2, last = 30,
                                     sampling = FALSE)

est_typed_ord <- remstimate(reh = reh_ord, stats = ts_typed_ord,
                              method = "MLE", ncores = 1L)

stacked_typed_ord <- stack_stats(ts_typed_ord, reh_ord)
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
  info = "typed ordinal: remstimate MLE matches clogit (consider_type = 'separate')"
)

# ---------------------------------------------------------------------------
# SECTION 5: Sampled stats — remstimate vs weighted Poisson GLM
# ---------------------------------------------------------------------------
samp_num <- 10L

ts_samp_int <- tomstats(effects, reh = reh_int,
                                    attr_actors = info,
                                    memory = "decay", memory_value = 1000,
                                    first = 2, last = 30,
                                    sampling = TRUE, samp_num = samp_num, seed = 1L)

est_samp_int <- remstimate(reh = reh_int, stats = ts_samp_int,
                             method = "MLE", ncores = 1L)

stacked_samp_int <- stack_stats(ts_samp_int, reh_int)
df_samp_int <- stacked_samp_int$remstats_stack

fit_glm_samp <- glm(obs ~ -1 + baseline + inertia + indegreeSender +
                      outdegreeSender + offset(log_interevent),
                    weights = weight,
                    family = poisson, data = df_samp_int)

expect_equal(
  unname(est_samp_int$coefficients),
  unname(coef(fit_glm_samp)),
  tolerance = tol,
  info = "sampled interval: remstimate MLE matches weighted Poisson GLM"
)

# ---------------------------------------------------------------------------
# SECTION 6: Sampled stats — remstimate vs weighted clogit (ordinal)
# ---------------------------------------------------------------------------
ts_samp_ord <- tomstats(effects_no_baseline, reh = reh_ord,
                                    attr_actors = info,
                                    memory = "decay", memory_value = 1000,
                                    first = 2, last = 30,
                                    sampling = TRUE, samp_num = samp_num, seed = 1L)

est_samp_ord <- remstimate(reh = reh_ord, stats = ts_samp_ord,
                             method = "MLE", ncores = 1L)

stacked_samp_ord <- stack_stats(ts_samp_ord, reh_ord)
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
  info = "sampled ordinal: remstimate MLE matches weighted clogit"
)



# ---------------------------------------------------------------------------
# SECTION 7: AOM — interval, sender + receiver stacks vs Poisson GLM / clogit
# ---------------------------------------------------------------------------
reh_actor_test <- remify(
  edgelist = history_sub,
  model    = "actor",
  directed = TRUE
)

sender_effects   <- ~ indegreeSender()
receiver_effects <- ~ inertia(consider_type = "ignore") +
  indegreeReceiver(consider_type = "ignore")

ts_aom_int <- aomstats(
  reh              = reh_actor_test,
  sender_effects   = sender_effects,
  receiver_effects = receiver_effects,
  memory           = "decay",
  memory_value     = 1000,
  start            = 2,
  stop             = 30
)

est_aom_int <- remstimate(
  reh    = reh_actor_test,
  stats  = ts_aom_int,
  method = "MLE",
  ncores = 1L
)

stacked_aom_int <- stack_stats(stats = ts_aom_int, reh = reh_actor_test)

# -- Sender: Poisson GLM with log_interevent offset --------------------------
df_sender_int <- stacked_aom_int$sender_stack

fit_glm_sender <- glm(
  obs ~ -1 + baseline + indegreeSender + offset(log_interevent),
  family = poisson,
  data   = df_sender_int
)

expect_equal(
  unname(est_aom_int$sender_model$coefficients),
  unname(coef(fit_glm_sender)),
  tolerance = tol,
  info = "AOM interval: sender remstimate MLE matches Poisson GLM"
)

# -- Receiver: clogit stratified by event ------------------------------------
# No log_interevent: receiver model is a pure choice model
df_receiver_int <- stacked_aom_int$receiver_stack

fit_clogit_receiver <- survival::clogit(
  obs ~ -1 + inertia + indegreeReceiver + survival::strata(event),
  data = df_receiver_int
)

expect_equal(
  unname(est_aom_int$receiver_model$coefficients),
  unname(coef(fit_clogit_receiver)),
  tolerance = tol,
  info = "AOM interval: receiver remstimate MLE matches clogit"
)

# ---------------------------------------------------------------------------
# SECTION 8: AOM — ordinal, sender clogit + receiver clogit
# ---------------------------------------------------------------------------
reh_actor_ord <- remify(
  edgelist = history_sub,
  model    = "actor",
  directed = TRUE,
  ordinal  = TRUE
)

ts_aom_ord <- aomstats(
  reh              = reh_actor_ord,
  sender_effects   = ~ indegreeSender(),
  receiver_effects = ~ inertia(consider_type = "ignore") +
    outdegreeReceiver(consider_type = "ignore"),
  memory           = "decay",
  memory_value     = 1000,
  start            = 2,
  stop             = 30
)

est_aom_ord <- remstimate(
  reh    = reh_actor_ord,
  stats  = ts_aom_ord,
  method = "MLE",
  ncores = 1L
)

stacked_aom_ord <- stack_stats(ts_aom_ord, reh_actor_ord)

# -- Sender: clogit (ordinal => no log_interevent, no baseline) --------------
df_sender_ord <- stacked_aom_ord$sender_stack

# baseline column is present but dropped for ordinal (absorbed by strata)
fit_clogit_sender <- survival::clogit(
  obs ~ -1 + indegreeSender + survival::strata(event),
  data = df_sender_ord
)

expect_equal(
  unname(est_aom_ord$sender_model$coefficients),
  unname(coef(fit_clogit_sender)),
  tolerance = tol,
  info = "AOM ordinal: sender remstimate MLE matches clogit"
)

# -- Receiver: clogit (same as interval case, receiver is always ordinal) ----
df_receiver_ord <- stacked_aom_ord$receiver_stack

fit_clogit_receiver_ord <- survival::clogit(
  obs ~ -1 + inertia + outdegreeReceiver + survival::strata(event),
  data = df_receiver_ord
)

expect_equal(
  unname(est_aom_ord$receiver_model$coefficients),
  unname(coef(fit_clogit_receiver_ord)),
  tolerance = tol,
  info = "AOM ordinal: receiver remstimate MLE matches clogit"
)

# ---------------------------------------------------------------------------
# SECTION 9: AOM — NULL sender or receiver stats
# ---------------------------------------------------------------------------
ts_aom_rec_only <- aomstats(
  reh              = reh_actor_test,
  receiver_effects = ~ inertia(consider_type = FALSE),
  memory           = "decay",
  memory_value     = 1000,
  start            = 2,
  stop             = 30
)

stacked_rec_only <- stack_stats(ts_aom_rec_only, reh_actor_test)

expect_null(
  stacked_rec_only$sender_stack,
  info = "AOM: sender_stack is NULL when no sender_effects specified"
)
expect_true(
  is.data.frame(stacked_rec_only$receiver_stack),
  info = "AOM: receiver_stack is a data frame when only receiver_effects specified"
)

ts_aom_send_only <- aomstats(
  reh            = reh_actor_test,
  sender_effects = ~ indegreeSender(),
  memory         = "decay",
  memory_value   = 1000,
  start          = 2,
  stop           = 30
)

stacked_send_only <- stack_stats(ts_aom_send_only, reh_actor_test)

expect_true(
  is.data.frame(stacked_send_only$sender_stack),
  info = "AOM: sender_stack is a data frame when only sender_effects specified"
)
expect_null(
  stacked_send_only$receiver_stack,
  info = "AOM: receiver_stack is NULL when no receiver_effects specified"
)

