# GAM backend (mgcv) - interval (Poisson) and ordinal (stratified cox.ph)
#
# The ordinal GAM maximises the stratified Cox partial likelihood, which with a
# constant dummy time and one stratum per risk set IS the conditional-logit
# likelihood. So an ordinal GAM with no smooth terms must reproduce
# survival::clogit exactly - that equivalence is what the first section pins.

library(tinytest)

if (!requireNamespace("mgcv", quietly = TRUE) ||
    !requireNamespace("survival", quietly = TRUE)) exit_file("mgcv/survival")

library(remify)
library(remstats)
library(survival)   # clogit() calls coxph() unqualified

data(history, package = "remstats")
colnames(history)[colnames(history) == "setting"] <- "type"
history_sub <- history[1:40, ]

tol <- 1e-4

effects_bl <- ~ 1 + inertia(consider_type = FALSE) +
                  indegreeSender(consider_type = FALSE) +
                  outdegreeSender(consider_type = FALSE)
effects_no_bl <- ~ inertia(consider_type = FALSE) +
                     indegreeSender(consider_type = FALSE) +
                     outdegreeSender(consider_type = FALSE)

# ---------------------------------------------------------------------------
# SECTION 1: ordinal GAM with no smooths == survival::clogit
# ---------------------------------------------------------------------------
reh_ord <- remify(edgelist = history_sub, model = "tie",
                  riskset = "active", ordinal = TRUE)

ts_ord <- tomstats(effects_no_bl, reh = reh_ord,
                   memory = "decay", memory_value = 1000,
                   first = 2, last = 30, sampling = FALSE)

fit_ord <- remstimate(reh = reh_ord, stats = ts_ord, gam = list(smooths = NULL))

df_ord <- stack_stats(ts_ord, reh_ord)$remstats_stack
fit_clogit <- survival::clogit(
  obs ~ -1 + inertia + indegreeSender + outdegreeSender +
    survival::strata(time_index),
  data = df_ord
)

expect_equal(unname(coef(fit_ord)), unname(coef(fit_clogit)), tolerance = tol,
             info = "ordinal GAM (no smooths) matches clogit coefficients")

expect_equal(as.numeric(fit_ord$loglik), as.numeric(fit_clogit$loglik[2L]),
             tolerance = tol,
             info = "ordinal GAM log-likelihood is the conditional-logit loglik")

# null = beta 0 conditional logit = -sum_m log|R_m|; for an unsampled design
# with a constant risk-set size that is -M*log(D)
M_ord <- fit_ord$df.null
D_ord <- nrow(df_ord) / M_ord
expect_equal(fit_ord$null.deviance, -2 * as.numeric(fit_clogit$loglik[1L]),
             tolerance = tol,
             info = "ordinal GAM null deviance matches clogit's null")
expect_equal(fit_ord$null.deviance, 2 * M_ord * log(D_ord), tolerance = tol,
             info = "ordinal GAM null deviance equals 2*M*log(D)")

# the ordinal MLE targets the same likelihood
est_ord <- remstimate(reh = reh_ord, stats = ts_ord, method = "MLE", ncores = 1L)
expect_equal(unname(coef(fit_ord)), unname(est_ord$coefficients),
             tolerance = tol,
             info = "ordinal GAM (no smooths) matches the C++ MLE")

expect_true(isTRUE(attr(fit_ord, "ordinal")),
            info = "ordinal attribute is set on an ordinal GAM")
expect_equal(attr(fit_ord, "engine"), "gam",
             info = "ordinal GAM always uses gam()")
expect_equal(attr(fit_ord, "family"), "cox.ph",
             info = "ordinal GAM records the cox.ph family")

# ---------------------------------------------------------------------------
# SECTION 2: sampled ordinal GAM == weighted clogit
# ---------------------------------------------------------------------------
ts_samp_ord <- tomstats(effects_no_bl, reh = reh_ord,
                        memory = "decay", memory_value = 1000,
                        first = 2, last = 30,
                        sampling = TRUE, samp_num = 10L, seed = 1L)

fit_samp_ord <- remstimate(reh = reh_ord, stats = ts_samp_ord,
                           gam = list(smooths = NULL))

df_samp_ord <- stack_stats(ts_samp_ord, reh_ord)$remstats_stack
fit_clogit_samp <- survival::clogit(
  obs ~ -1 + inertia + indegreeSender + outdegreeSender +
    survival::strata(time_index),
  weights = weight, method = "approximate", data = df_samp_ord
)

# the case-control correction enters as offset(-log pi_d) on the linear
# predictor, which is algebraically the weighted risk-set sum clogit forms
expect_equal(unname(coef(fit_samp_ord)), unname(coef(fit_clogit_samp)),
             tolerance = tol,
             info = "sampled ordinal GAM matches weighted clogit")

# remstats samples the alternatives uniformly, so the observed dyad is absent
# from many sampled risk sets. A risk set with no event contributes nothing to a
# partial likelihood, so it must not be counted in the null either.
n_informative <- sum(tapply(df_samp_ord$obs, df_samp_ord$time_index, sum) > 0)
expect_true(n_informative < fit_samp_ord$df.null,
            info = "the sampled design really does have empty risk sets")
expect_equal(fit_samp_ord$null.deviance,
             -2 * as.numeric(fit_clogit_samp$loglik[1L]), tolerance = tol,
             info = "sampled ordinal null deviance matches weighted clogit")

est_samp_ord <- remstimate(reh = reh_ord, stats = ts_samp_ord,
                           method = "MLE", ncores = 1L)
expect_equal(fit_samp_ord$null.deviance, est_samp_ord$null.deviance,
             tolerance = tol,
             info = "sampled ordinal null deviance matches the C++ MLE")

# ---------------------------------------------------------------------------
# SECTION 3: ordinal GAM with a smooth term
# ---------------------------------------------------------------------------
fit_sm <- remstimate(reh = reh_ord, stats = ts_ord,
                     gam = list(smooths = "inertia", k = 5))

expect_true(inherits(fit_sm, "remstimate_gam"),
            info = "smoothed ordinal fit keeps the remstimate_gam class")
expect_false(is.null(fit_sm$smooth_map),
             info = "smoothed ordinal fit records a smooth map")
expect_equal(fit_sm$smooth_map$statistic, "inertia",
             info = "the smooth is built on the requested statistic")
expect_true(all(is.finite(c(fit_sm$AIC, fit_sm$BIC, fit_sm$AICC, fit_sm$loglik))),
            info = "information criteria are finite for a smoothed ordinal fit")
expect_true(fit_sm$df.model > 0 && fit_sm$df.model < length(coef(fit_sm)) + 1,
            info = "effective df lies between 0 and the number of coefficients")
expect_false(isTRUE(attr(fit_sm, "coef_contract_1to1")),
             info = "a smoothed fit advertises the broken 1:1 coef contract")

# the smoothed statistic must not also appear as a parametric coefficient
expect_true(!"inertia" %in% names(coef(fit_sm, type = "parametric")),
            info = "a smoothed statistic has no parametric coefficient")

# S3 surface
expect_silent(invisible(capture.output(print(fit_sm))))
invisible(capture.output(s <- summary(fit_sm)))
expect_true(!is.null(s$smoothTab),
            info = "summary reports the smooth table for an ordinal GAM")
expect_true(isTRUE(s$ordinal), info = "summary carries the ordinal flag")
expect_equal(as.numeric(logLik(fit_sm)), fit_sm$loglik,
             info = "logLik method returns the stored partial log-likelihood")

dg <- diagnostics(fit_sm, reh_ord)
expect_true(inherits(dg, "diagnostics_gam"),
            info = "ordinal GAM diagnostics return a diagnostics_gam object")
expect_false(is.null(dg$recall),
             info = "recall is computed for an ordinal GAM")

# ---------------------------------------------------------------------------
# SECTION 4: ordinal guard rails
# ---------------------------------------------------------------------------
expect_error(
  remstimate(reh = reh_ord, stats = ts_ord,
             gam = list(smooths = "inertia", k = 5, engine = "bam")),
  pattern = "bam",
  info = "bam is refused for ordinal models (no general families)"
)

expect_error(
  remstimate(reh = reh_ord, stats = ts_ord,
             gam = list(smooths = "inertia", k = 5, method = "fREML")),
  pattern = "fREML",
  info = "fREML is refused for ordinal models"
)

# a statistic constant within every risk set is absorbed by that risk set's
# baseline hazard, so it is dropped rather than handed to mgcv (which fails
# with a bare 'non-conformable arguments')
df_alias <- df_ord
df_alias$baseline <- 1                              # constant everywhere
df_alias$eventcov <- as.numeric(df_alias$time_index) # constant within a risk set
stacked_alias <- structure(
  list(remstats_stack = df_alias,
       stat_names = c("inertia", "indegreeSender", "outdegreeSender",
                      "baseline", "eventcov"),
       ordinal = TRUE, model = "tie"),
  class = "remstats_stacked"
)

expect_warning(
  fit_alias <- remstimate(reh = reh_ord, stats = stacked_alias,
                          gam = list(smooths = NULL)),
  pattern = "constant within every risk set",
  info = "risk-set-constant statistics are dropped with a warning"
)
expect_equal(attr(fit_alias, "statistics"),
             c("inertia", "indegreeSender", "outdegreeSender"),
             info = "only the identifiable statistics survive")

# but an aliased statistic the caller explicitly asked to smooth is an error,
# not a silent drop: 'bs' / 'k' / 'pc' are positional, so removing an entry
# from 'smooths' would shift the remaining settings onto other statistics
expect_error(
  remstimate(reh = reh_ord, stats = stacked_alias,
             gam = list(smooths = c("eventcov", "inertia"), k = c(5, 8))),
  pattern = "constant within every risk set",
  info = "smoothing a risk-set-constant statistic is refused"
)

# more than one event in a risk set puts cox.ph on Peto's approximation
df_ties <- df_ord
first_rs <- df_ties$time_index == df_ties$time_index[1L]
df_ties$obs[which(first_rs & df_ties$obs == 0L)[1L]] <- 1L
stacked_ties <- structure(
  list(remstats_stack = df_ties,
       stat_names = c("inertia", "indegreeSender", "outdegreeSender"),
       ordinal = TRUE, model = "tie"),
  class = "remstats_stacked"
)
expect_warning(
  remstimate(reh = reh_ord, stats = stacked_ties, gam = list(smooths = NULL)),
  pattern = "Peto",
  info = "ties in a risk set warn that the partial likelihood is approximate"
)

# ---------------------------------------------------------------------------
# SECTION 5: interval GAM
# ---------------------------------------------------------------------------
reh_int <- remify(edgelist = history_sub, model = "tie", riskset = "active")
ts_int  <- tomstats(effects_bl, reh = reh_int,
                    memory = "decay", memory_value = 1000,
                    first = 2, last = 30, sampling = FALSE)

fit_int <- remstimate(reh = reh_int, stats = ts_int, gam = list(smooths = NULL))
df_int  <- stack_stats(ts_int, reh_int)$remstats_stack
fit_pois <- glm(obs ~ -1 + baseline + inertia + indegreeSender +
                  outdegreeSender + offset(log_interevent),
                family = poisson, data = df_int)

expect_equal(unname(coef(fit_int)), unname(coef(fit_pois)), tolerance = tol,
             info = "interval GAM (no smooths) still matches the Poisson GLM")
expect_false(isTRUE(attr(fit_int, "ordinal")),
             info = "interval fit is not flagged ordinal")
expect_equal(attr(fit_int, "family"), "poisson",
             info = "interval fit records the poisson family")
expect_true(isTRUE(fit_int$converged),
            info = "interval fit reports convergence")

# interval still requires a baseline column; ordinal must not. Note remstats
# adds 'baseline' for interval timing whenever the effects formula keeps its
# (implicit) intercept, so suppressing it takes an explicit -1.
ts_int_no_bl <- tomstats(~ -1 + inertia(consider_type = FALSE) +
                             indegreeSender(consider_type = FALSE) +
                             outdegreeSender(consider_type = FALSE),
                         reh = reh_int, memory = "decay", memory_value = 1000,
                         first = 2, last = 30, sampling = FALSE)

expect_error(
  remstimate(reh = reh_int, stats = ts_int_no_bl, gam = list(smooths = NULL)),
  pattern = "baseline",
  info = "interval GAM still demands an intercept column"
)

# the same statistics under ordinal timing need no baseline at all
expect_silent(
  remstimate(reh = reh_ord, stats = ts_ord, gam = list(smooths = NULL))
)

# ---------------------------------------------------------------------------
# SECTION 6: shape constraints (mgcv::scasm)
# ---------------------------------------------------------------------------
# scasm() is a recent addition to mgcv; skip the section rather than fail the
# whole file on an older installation.
if (!exists("scasm", envir = asNamespace("mgcv"), inherits = FALSE))
  exit_file("mgcv without scasm()")

# Partial effect of a smooth on a grid of the statistic, holding the other
# statistics at the values of one design row. This is the quantity a shape
partial_effect <- function(fit, stat, n = 25L) {
  df   <- fit$stacked_data
  grid <- df[rep(1L, n), , drop = FALSE]
  grid[[stat]] <- seq(min(df[[stat]], na.rm = TRUE),
                      max(df[[stat]], na.rm = TRUE), length.out = n)
  as.numeric(stats::predict(fit$backend_fit, newdata = grid,
                            type = "terms")[, paste0("s(", stat, ")")])
}

# In this design the inertia effect decreases, so "m-" is the constraint the
# data can meet; "m+" is used below to test the infeasible branch.
fit_mono <- remstimate(reh = reh_int, stats = ts_int,
                       gam = list(smooths = "inertia", k = 5,
                                  constraints = "m-"))

expect_equal(attr(fit_mono, "engine"), "scasm",
             info = "a shape constraint routes the fit to scasm()")
expect_equal(attr(fit_mono, "constrained"), "inertia",
             info = "the constrained statistic is recorded on the fit")
expect_equal(attr(fit_mono, "constraints")$inertia, "m-",
             info = "the mgcv constraint code is recorded on the fit")
expect_equal(fit_mono$smooth_map$shape, "m-",
             info = "the smooth map reports the constraint")
expect_true(isTRUE(fit_mono$converged),
            info = "a constrained fit reports Fellner-Schall convergence")

f_mono <- partial_effect(fit_mono, "inertia")
expect_true(all(diff(f_mono) <= 1e-8),
            info = "an 'm-' smooth really is monotone decreasing")

# and the constraint is what makes the difference: unconstrained, the same
# smooth is free to move either way and the engine stays gam()
fit_free <- remstimate(reh = reh_int, stats = ts_int,
                       gam = list(smooths = "inertia", k = 5))
expect_equal(attr(fit_free, "engine"), "gam",
             info = "without constraints the engine is unchanged")
expect_equal(fit_free$smooth_map$shape, "",
             info = "an unconstrained smooth records no shape")

# ordinal timing: cox.ph under scasm
fit_mono_ord <- remstimate(reh = reh_ord, stats = ts_ord,
                           gam = list(smooths = "inertia", k = 5,
                                      constraints = "m-"))
expect_equal(attr(fit_mono_ord, "engine"), "scasm",
             info = "ordinal shape constraints also use scasm()")
expect_equal(attr(fit_mono_ord, "family"), "cox.ph",
             info = "a constrained ordinal fit keeps the cox.ph family")
expect_true(all(diff(partial_effect(fit_mono_ord, "inertia")) <= 1e-8),
            info = "an ordinal 'm-' smooth is monotone decreasing too")
expect_true(all(is.finite(c(fit_mono_ord$AIC, fit_mono_ord$BIC,
                            fit_mono_ord$loglik))),
            info = "information criteria are finite for a constrained fit")

# S3 surface holds up
expect_silent(invisible(capture.output(print(fit_mono))))
invisible(capture.output(s_con <- summary(fit_mono)))
expect_true(!is.null(s_con$smoothTab),
            info = "summary reports the smooth table for a constrained GAM")
expect_true(inherits(diagnostics(fit_mono, reh_int), "diagnostics_gam"),
            info = "a constrained GAM still yields diagnostics")

# named form: the constraint reaches only the smooth it names
fit_named <- remstimate(reh = reh_int, stats = ts_int,
                        gam = list(smooths = c("inertia", "indegreeSender"),
                                   k = 5, constraints = c(inertia = "m-")))
expect_equal(fit_named$smooth_map$shape, c("m-", ""),
             info = "a named constraint applies only to the smooth it names")

# codes combine, both in one string and as a list element
fit_combo <- remstimate(reh = reh_int, stats = ts_int,
                        gam = list(smooths = "inertia", k = 5,
                                   constraints = "m- c+"))
expect_equal(fit_combo$smooth_map$shape, "m- c+",
             info = "several constraint codes can be combined for one smooth")
fit_list <- remstimate(reh = reh_int, stats = ts_int,
                       gam = list(smooths = c("inertia", "indegreeSender"),
                                  k = 5,
                                  constraints = list(c("m-", "c+"), "none")))
expect_equal(fit_list$smooth_map$shape, c("m- c+", ""),
             info = "the list form takes a character vector per smooth")

f_combo <- partial_effect(fit_combo, "inertia")
expect_true(all(diff(f_combo) <= 1e-8),
            info = "the combined smooth is still monotone decreasing")
expect_true(all(diff(f_combo, differences = 2) >= -1e-8),
            info = "the combined smooth is also convex")

# ---------------------------------------------------------------------------
# SECTION 7: shape-constraint guard rails
# ---------------------------------------------------------------------------
# gam()/bam() build the 'sc' basis but ignore its constraint matrix, so they
# would return a silently unconstrained fit
expect_error(
  remstimate(reh = reh_int, stats = ts_int,
             gam = list(smooths = "inertia", k = 5, constraints = "m-",
                        engine = "gam")),
  pattern = "cannot impose shape constraints",
  info = "engine = 'gam' is refused when constraints are given"
)
expect_error(
  remstimate(reh = reh_int, stats = ts_int,
             gam = list(smooths = "inertia", k = 5, constraints = "m-",
                        engine = "bam")),
  pattern = "cannot impose shape constraints",
  info = "engine = 'bam' is refused when constraints are given"
)

# scasm() has no smoothness criterion to choose - it is Fellner-Schall only -
# so a 'method' would be swallowed by '...' and silently ignored
expect_error(
  remstimate(reh = reh_int, stats = ts_int,
             gam = list(smooths = "inertia", k = 5, constraints = "m-",
                        method = "REML")),
  pattern = "does not apply to engine = 'scasm'",
  info = "'method' is refused for a constrained fit"
)

# only bs = "sc" carries a constraint
expect_error(
  remstimate(reh = reh_int, stats = ts_int,
             gam = list(smooths = "inertia", k = 5, constraints = "m-",
                        bs = "cr")),
  pattern = "cannot carry the constraint",
  info = "an incompatible basis is refused for a constrained smooth"
)
expect_equal(
  attr(remstimate(reh = reh_int, stats = ts_int,
                  gam = list(smooths = "inertia", k = 5, constraints = "m-",
                             bs = "sc")), "engine"),
  "scasm",
  info = "bs = 'sc' may be given explicitly"
)

expect_error(
  remstimate(reh = reh_int, stats = ts_int,
             gam = list(smooths = "inertia", k = 5,
                        constraints = "increasing")),
  pattern = "unknown shape constraint",
  info = "only mgcv's own codes are accepted"
)
expect_error(
  remstimate(reh = reh_int, stats = ts_int,
             gam = list(smooths = "inertia", k = 5, constraints = "m+ m-")),
  pattern = "contradictory shape constraints",
  info = "a smooth cannot be both increasing and decreasing"
)
expect_error(
  remstimate(reh = reh_int, stats = ts_int,
             gam = list(smooths = "inertia", k = 5,
                        constraints = c(baseline = "m-"))),
  pattern = "not among 'smooths'",
  info = "a constraint must name a smoothed statistic"
)
expect_error(
  remstimate(reh = reh_int, stats = ts_int,
             gam = list(smooths = c("inertia", "indegreeSender"), k = 5,
                        constraints = c("m-", "m-", "m-"))),
  pattern = "must be length 1",
  info = "a positional constraint vector must match 'smooths'"
)
expect_error(
  remstimate(reh = reh_int, stats = ts_int,
             gam = list(smooths = NULL, constraints = "m-")),
  pattern = "'smooths' is empty",
  info = "there is nothing to constrain without a smooth"
)

# a constraint the data run flatly against pins the fit to the edge of the
# feasible set, where mgcv fails with a bare "subscript out of bounds"
expect_error(
  remstimate(reh = reh_int, stats = ts_int,
             gam = list(smooths = "inertia", k = 5, constraints = "m+")),
  pattern = "binds over the whole range",
  info = "an infeasible constraint gives an actionable error"
)

# f(pc) = 0 sits on the boundary of f >= 0, so a positivity constraint drops the
# point constraint instead of handing mgcv a singular problem
expect_true(
  !grepl("pc", remstimate:::.gam_smooth_terms("inertia", "tp", 5, 0,
                                              list("+")), fixed = TRUE),
  info = "a '+' constraint drops the point constraint"
)
expect_true(
  grepl("pc = 0", remstimate:::.gam_smooth_terms("inertia", "tp", 5, 0,
                                                 list("m-")), fixed = TRUE),
  info = "a monotonicity constraint keeps the point constraint"
)
