# GAM backend (mgcv) - time-varying effects, s(t, by = STATISTIC)
#
# A time-varying effect replaces the constant coefficient of a statistic by a
# smooth function of the event time. Two mgcv properties make or break it, and
# both are pinned below:
#
#   (1) a smooth with a numeric, varying 'by' variable is NOT centred, so the
#       smooth carries the whole coefficient function - that is what makes the
#       fitted curve readable as beta(t);
#   (2) because of (1) the statistic must not ALSO enter linearly. It would be
#       exactly collinear with the constant component of the by-smooth, and mgcv
#       resolves that silently: the linear coefficient comes back as 0 with
#       standard error 0 and a NaN p-value, with no warning at all.
#
# Section 1 is the regression test for both.

library(tinytest)

if (!requireNamespace("mgcv", quietly = TRUE)) exit_file("mgcv")

library(remify)
library(remstats)

data(history, package = "remstats")
colnames(history)[colnames(history) == "setting"] <- "type"
history_sub <- history[1:60, ]

tol <- 1e-6

effects_bl <- ~ 1 + inertia(consider_type = FALSE) +
                  indegreeSender(consider_type = FALSE)

reh <- remify(edgelist = history_sub, model = "tie", riskset = "active")
ts  <- tomstats(effects_bl, reh = reh, first = 2, last = 40)

fit <- remstimate(reh = reh, stats = ts, gam = list(tve = "inertia"))

# ---------------------------------------------------------------------------
# SECTION 1: the two identification properties
# ---------------------------------------------------------------------------
sm_by <- Filter(function(s) !identical(s$by, "NA"), fit$backend_fit$smooth)
sm_bl <- Filter(function(s)  identical(s$by, "NA"), fit$backend_fit$smooth)

expect_equal(length(sm_by), 1L,
             info = "one s(t, by = ...) term per 'tve' statistic")
expect_equal(sm_by[[1]]$by, "inertia",
             info = "the by variable is the statistic named in 'tve'")

# (1) uncentred: mgcv absorbs no constraint, so the whole basis is kept
expect_equal(attr(sm_by[[1]], "nCons"), 0L,
             info = "s(t, by = X) is NOT centred - it carries the level of beta(t)")
expect_equal(sm_by[[1]]$df, sm_by[[1]]$bs.dim,
             info = "no constraint absorbed, so df == k for a time-varying effect")
expect_true(sm_by[[1]]$null.space.dim >= 2L,
            info = "the by-smooth's null space keeps the constant in time")

# the time-varying baseline, by contrast, IS centred - that is what keeps
# 'baseline' interpretable as the row-average log-rate
expect_equal(length(sm_bl), 1L, info = "baseline = TRUE adds one s(t) term")
expect_equal(attr(sm_bl[[1]], "nCons"), 1L,
             info = "s(t) (no by variable) is centred by mgcv")

# (2) the statistic must not also enter linearly
expect_false("inertia" %in% names(coef(fit, type = "parametric")),
             info = "a 'tve' statistic is removed from the linear terms")
expect_true("indegreeSender" %in% names(coef(fit, type = "parametric")),
            info = "statistics not named in 'tve' still enter linearly")
expect_true("baseline" %in% names(coef(fit, type = "parametric")),
            info = "baseline still carries the overall log-rate")
expect_false(any(!is.finite(fit$se)),
             info = "no aliased coefficient (an aliased one has se exactly 0/NaN)")

# ---------------------------------------------------------------------------
# SECTION 2: wiring - the fit is the mgcv fit of the intended formula
# ---------------------------------------------------------------------------
df_ref <- stack_stats(ts, reh, add_actors = FALSE)$remstats_stack
df_ref$samp_offset <- 0
df_ref$t <- as.numeric(reh$edgelist$time)[df_ref$time_index]

ref <- mgcv::gam(
  obs ~ -1 + offset(log_interevent + samp_offset) +
    s(t, bs = "tp") + s(t, bs = "tp", by = inertia) +
    baseline + indegreeSender,
  family = stats::poisson(), data = df_ref, method = "REML")

expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = tol,
             info = "coefficients match the direct mgcv fit")
expect_equal(sum(fit$backend_fit$edf), sum(ref$edf), tolerance = tol,
             info = "effective degrees of freedom match the direct mgcv fit")

# the REM log-likelihood is the Poisson one minus the sum(log dt) over events
off_obs <- sum(df_ref$log_interevent[df_ref$obs == 1L])
expect_equal(as.numeric(logLik(fit)),
             as.numeric(stats::logLik(ref)) - off_obs, tolerance = tol,
             info = "loglik drops the offset constant, as in the nonlinear path")

# the time axis is the event time, matched through time_index
expect_equal(unname(attr(fit, "time_range")),
             unname(range(df_ref$t)), tolerance = tol,
             info = "t is reh$edgelist$time indexed by time_index")

# ---------------------------------------------------------------------------
# SECTION 3: guards
# ---------------------------------------------------------------------------
expect_error(remstimate(reh, ts, gam = list(tve = "inertia",
                                            smooths = "indegreeSender")),
             pattern = "cannot both be given",
             info = "a fit is either a nonlinear-effects or a time-varying one")

reh_ord <- remify(edgelist = history_sub, model = "tie", riskset = "active",
                  ordinal = TRUE)
ts_ord  <- tomstats(~ inertia(consider_type = FALSE) +
                      indegreeSender(consider_type = FALSE),
                    reh = reh_ord, first = 2, last = 40)
expect_error(remstimate(reh_ord, ts_ord, gam = list(tve = "inertia")),
             pattern = "not available for ordinal models",
             info = "ordinal remify replaces the event times by their ranks")

expect_error(remstimate(reh, ts, gam = list(tve = "baseline")),
             pattern = "cannot name 'baseline'",
             info = "the time-varying baseline is requested with baseline =")
expect_error(remstimate(reh, ts, gam = list(tve = "inertia",
                                            constraints = "m+")),
             pattern = "does not apply to 'tve'",
             info = "shape constraints are a nonlinear-effects feature")
expect_error(remstimate(reh, ts, gam = list(tve = "inertia", pc = 0)),
             pattern = "'pc' does not apply",
             info = "a point constraint would pin beta(t) to zero at t = pc")
expect_error(remstimate(reh, ts, gam = list(tve = "inertia", bs = "re")),
             pattern = "not available for a time-varying effect")
expect_error(remstimate(reh, ts, gam = list(tve = "inertia", bs = "sc")),
             pattern = "not available for a time-varying effect")
expect_error(remstimate(reh, ts, gam = list(tve = "notAStatistic")),
             pattern = "tve name\\(s\\) not found")
expect_error(remstimate(reh, ts, gam = list(tve = c("inertia", "inertia"))),
             pattern = "more than once")
expect_error(remstimate(reh, ts, gam = list(tve = "inertia", k = 500)),
             pattern = "distinct event times",
             info = "k cannot exceed the number of risk sets")
expect_error(remstimate(reh, ts, gam = list(tve = "inertia", baseline = "yes")),
             pattern = "must be TRUE")
expect_error(remstimate(reh, ts, gam = list(smooths = "inertia",
                                            baseline = FALSE)),
             pattern = "'baseline' was given but 'tve' is empty",
             info = "'baseline' only configures a time-varying-effects model")

# ---------------------------------------------------------------------------
# SECTION 4: the baseline switch
# ---------------------------------------------------------------------------
fit_nobl <- remstimate(reh, ts, gam = list(tve = "inertia", baseline = FALSE))
expect_equal(length(fit_nobl$backend_fit$smooth), 1L,
             info = "baseline = FALSE leaves only the s(t, by = ...) term")
expect_equal(fit_nobl$backend_fit$smooth[[1]]$by, "inertia")

fit_k <- remstimate(reh, ts, gam = list(tve = "inertia", baseline = 12))
expect_equal(fit_k$backend_fit$smooth[[1]]$bs.dim, 12L,
             info = "a numeric 'baseline' sets the basis dimension of s(t)")
expect_true(identical(fit_k$backend_fit$smooth[[1]]$by, "NA"),
            info = "the baseline smooth is listed first")

# ---------------------------------------------------------------------------
# SECTION 5: object contract and methods
# ---------------------------------------------------------------------------
expect_equal(attr(fit, "gam_type"), "tve")
expect_equal(attr(fit, "tve"), "inertia")
expect_equal(length(attr(fit, "smooths")), 0L,
             info = "'smooths' is empty in a time-varying fit")
expect_false(attr(fit, "coef_contract_1to1"),
             info = "the coefficient vector holds one entry per basis function")
expect_true(all(c("statistic", "by", "term") %in% names(fit$smooth_map)))
expect_equal(fit$smooth_map$statistic[fit$smooth_map$term == "s(t):inertia"],
             "inertia",
             info = "smooth_map resolves a time-varying term to its by variable")
expect_true(is.na(fit$smooth_map$by[fit$smooth_map$term == "s(t)"]),
            info = "the baseline smooth has no by variable")
expect_equal(sort(fit$linear), sort(c("baseline", "indegreeSender")))

# print() and summary() write straight to stdout, so both are captured
printed <- capture.output(print(fit))
expect_true(any(grepl("Time-varying effects", printed, fixed = TRUE)),
            info = "print() labels the terms as time-varying effects")
expect_true(any(grepl("s(t):inertia", printed, fixed = TRUE)))
expect_true(length(coef(fit)) > length(coef(fit, type = "parametric")))

sm_tab <- NULL
summarised <- capture.output(sm_tab <- summary(fit)$smoothTab)
expect_true(any(grepl("time-varying effects", summarised, fixed = TRUE)),
            info = "summary() heads the term table as time-varying effects")
expect_true(all(c("s(t)", "s(t):inertia") %in% rownames(sm_tab)),
            info = "both time smooths are tested as whole terms in summary()")

dg <- diagnostics(fit, reh)
expect_true(inherits(dg, "diagnostics_gam"))
expect_true(!is.null(dg$recall))
expect_true(!is.null(dg$concurvity),
            info = "several smooths of t: concurvity is reported")

pdf(NULL)
expect_silent(plot(fit, reh, which = 9))
dev.off()

# ---------------------------------------------------------------------------
# SECTION 6: recovery of a known beta(t)
# ---------------------------------------------------------------------------
# Statistics and risk sets are taken from the real design; only the response is
# simulated, from a known time-varying coefficient. The design is replicated so
# the check is powered rather than testing a smoother against noise.
if (at_home()) {

  set.seed(20260806)
  df_sim <- df_ref
  Tmax   <- max(df_sim$t)
  beta_true <- function(tt) 0.6 * sin(2 * pi * tt / Tmax)
  b0 <- -3.5

  R   <- 40
  big <- do.call(rbind, replicate(R, df_sim, simplify = FALSE))
  big$time_index <- big$time_index +
    rep((seq_len(R) - 1L) * 10000L, each = nrow(df_sim))
  eta <- b0 + beta_true(big$t) * big$inertia + 0.2 * big$indegreeSender
  big$obs <- rpois(nrow(big), exp(eta + big$log_interevent))

  sim <- mgcv::gam(
    obs ~ -1 + offset(log_interevent) + baseline + indegreeSender +
      s(t, bs = "tp", by = inertia, k = 12),
    family = stats::poisson(), data = big, method = "REML")

  grid <- data.frame(t = seq(min(big$t), Tmax, length.out = 9),
                     inertia = 1, indegreeSender = 0, baseline = 1,
                     log_interevent = 0)
  pr  <- stats::predict(sim, newdata = grid, type = "terms", se.fit = TRUE)
  est <- pr$fit[, "s(t):inertia"]
  se  <- pr$se.fit[, "s(t):inertia"]
  truth <- beta_true(grid$t)

  expect_true(max(abs(est - truth)) < 0.12,
              info = "the fitted curve recovers the true beta(t)")
  expect_true(mean(abs(est - truth) < 1.96 * se) >= 0.75,
              info = "pointwise intervals cover the truth at most grid points")
  expect_true(abs(stats::coef(sim)[["baseline"]] - b0) < 0.1,
              info = "the intercept is not absorbed by the time-varying effect")
}

# ---------------------------------------------------------------------------
# SECTION 7: a constant effect is recovered as a constant
# ---------------------------------------------------------------------------
# With select = TRUE the whole term can shrink away, so a statistic whose effect
# does not move over time should end up with an effective df near a straight
# line rather than a wiggly curve.
fit_sel <- remstimate(reh, ts, gam = list(tve = "inertia", select = TRUE))
edf_by <- with(fit_sel$smooth_map, edf[!is.na(by)])
expect_true(edf_by < 5,
            info = "select = TRUE lets a near-constant effect shrink")

# ---------------------------------------------------------------------------
# SECTION 8: case-control sampled statistics
# ---------------------------------------------------------------------------
ts_samp <- tryCatch(
  tomstats(effects_bl, reh = reh, first = 2, last = 40,
           sampling = TRUE, samp_num = 5),
  error = function(e) NULL)

if (!is.null(ts_samp) && inherits(ts_samp, "tomstats_sampled")) {
  fit_samp <- remstimate(reh, ts_samp, gam = list(tve = "inertia"))
  expect_true(inherits(fit_samp, "remstimate_gam"))
  expect_equal(attr(fit_samp, "gam_type"), "tve")
  expect_true(grepl("samp_offset", deparse(attr(fit_samp, "formula"))[1]) ||
              any(grepl("samp_offset", deparse(attr(fit_samp, "formula")))),
              info = "the case-control correction stays in the offset")
  expect_false("inertia" %in% names(coef(fit_samp, type = "parametric")))
}