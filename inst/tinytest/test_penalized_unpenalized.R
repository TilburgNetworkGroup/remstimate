library(tinytest)
library(remify)
library(remstats)
library(remstimate)

# ─────────────────────────────────────────────────────────────────────────────
# Penalised backends (GLMNET / SHRINKEM): which statistics are left unpenalised.
#
# Default rule (structural): a statistic is intercept-like - and therefore
# exempt from the penalty - when it is a 0/1 indicator. This covers 'baseline',
# the durem 'baseline.start' / 'baseline.end' process intercepts and FEtype_*
# type dummies. Two additive controls adjust the default:
#   * unpenalized : ADD names to the exemption
#   * penalized   : REMOVE names from it (force-penalise a 0/1 effect such as a
#                   p-shift dummy). 'penalized' wins when a name is in both.
# ─────────────────────────────────────────────────────────────────────────────

# ── 1. .intercept_like_stats(): structural 0/1 detection ──────────────────────

df_mix <- data.frame(
  baseline       = rep(1, 8),                 # all ones  -> indicator
  baseline.start = rep(c(0, 1), 4),           # 0/1       -> indicator
  psABAB.end     = c(0, 0, 1, 0, 1, 0, 0, 1), # 0/1 dummy -> indicator (a real effect!)
  inertia        = c(0, 2, 7, 1, 0, 3, 5, 2), # counts    -> NOT indicator
  reciprocity    = rnorm(8)                   # continuous-> NOT indicator
)
sn <- colnames(df_mix)

expect_equal(
  sort(remstimate:::.intercept_like_stats(df_mix, sn)),
  sort(c("baseline", "baseline.start", "psABAB.end")),
  info = "indicator detection flags all 0/1 columns (incl. a p-shift dummy)")

# a statistic not present in the design is silently ignored
expect_equal(
  sort(remstimate:::.intercept_like_stats(df_mix, c(sn, "does_not_exist"))),
  sort(c("baseline", "baseline.start", "psABAB.end")),
  info = "indicator detection ignores absent statistic names")

# ── 1c. .recall_ranks(): average ranks under ties (no stack-order inflation) ───

# distinct probabilities: rank follows the ordering; cum is cumulative-through
rr1 <- remstimate:::.recall_ranks(c(0.5, 0.3, 0.2), pos = 1L)
expect_equal(rr1$rank, 1,   info = "recall_ranks: highest prob -> rank 1")
expect_equal(rr1$cum,  0.5, info = "recall_ranks: cum prob through the top")
rr2 <- remstimate:::.recall_ranks(c(0.5, 0.3, 0.2), pos = 3L)
expect_equal(rr2$rank, 3,   info = "recall_ranks: lowest prob -> rank 3")
expect_equal(rr2$cum,  1.0, info = "recall_ranks: cum prob through the bottom")

# FULL TIE: every position gets the SAME average rank (D_t+1)/2, independent of
# where the observed outcome sits in the risk set. This is the durem fix - a
# model that cannot discriminate scores ~chance instead of being handed rank 1.
probs_tie <- rep(0.25, 4)
ranks_tie <- vapply(1:4, function(p) remstimate:::.recall_ranks(probs_tie, p)$rank, numeric(1))
expect_equal(unique(ranks_tie), 2.5,
             info = "recall_ranks: full tie -> average rank, position-independent")
expect_equal(1 - ranks_tie[1] / 4, 0.375,
             info = "recall_ranks: tied outcome scores ~chance rel_rank, not ~1")

# partial tie: tied outcomes share the mean of their ranks
rr3 <- remstimate:::.recall_ranks(c(0.4, 0.2, 0.2, 0.2), pos = 2L)
expect_equal(rr3$rank, 3, info = "recall_ranks: tied block gets the mean rank")

# ── 2. .resolve_unpenalized(): additive unpenalized / subtractive penalized ────

# default: exactly the intercept-like set
expect_equal(
  sort(remstimate:::.resolve_unpenalized(df_mix, sn)),
  sort(c("baseline", "baseline.start", "psABAB.end")),
  info = "resolver default == intercept-like set")

# penalized REMOVES a 0/1 effect from the exemption (the p-shift use case)
expect_false(
  "psABAB.end" %in% remstimate:::.resolve_unpenalized(df_mix, sn, penalized = "psABAB.end"),
  info = "penalized forces a 0/1 dummy back into the penalty")
expect_true(
  all(c("baseline", "baseline.start") %in%
        remstimate:::.resolve_unpenalized(df_mix, sn, penalized = "psABAB.end")),
  info = "penalized leaves the other intercepts exempt")

# unpenalized ADDS a continuous effect to the exemption
expect_true(
  "reciprocity" %in% remstimate:::.resolve_unpenalized(df_mix, sn, unpenalized = "reciprocity"),
  info = "unpenalized adds a continuous effect to the exemption")
expect_true(
  all(c("baseline", "baseline.start", "psABAB.end") %in%
        remstimate:::.resolve_unpenalized(df_mix, sn, unpenalized = "reciprocity")),
  info = "unpenalized is additive: the default set is retained")

# penalized wins over unpenalized when a name is given in both
expect_false(
  "psABAB.end" %in% remstimate:::.resolve_unpenalized(
    df_mix, sn, unpenalized = "psABAB.end", penalized = "psABAB.end"),
  info = "penalized takes precedence over unpenalized")

# matching is EXACT: a base name that is not itself a statistic does nothing
expect_true(
  "psABAB.end" %in% remstimate:::.resolve_unpenalized(df_mix, sn, penalized = "psABAB"),
  info = "exact matching: 'psABAB' does not touch 'psABAB.end'")
expect_false(
  "psABAB.end" %in% remstimate:::.resolve_unpenalized(df_mix, sn, penalized = "psABAB.end"),
  info = "exact matching: the full name 'psABAB.end' works")

# ── 2b. .check_penalty_names(): warn on names that match no statistic ──────────

expect_warning(
  remstimate:::.check_penalty_names(penalized = "psABAB", valid = sn),
  "not found among the model statistics",
  info = "unknown penalty name (typo / missing suffix) warns")
expect_silent(
  remstimate:::.check_penalty_names(penalized = "psABAB.end",
                                    unpenalized = "reciprocity", valid = sn))

# ── 3. GLMNET integration: baseline is the (unpenalised) intercept ────────────

runs_ok <- function(expr) tryCatch({ force(expr); TRUE }, error = function(e) FALSE)

if (requireNamespace("glmnet", quietly = TRUE)) {

  el <- data.frame(
    time   = 1:16,
    actor1 = c(1, 2, 1, 3, 2, 1, 3, 2, 4, 1, 2, 3, 4, 1, 2, 5),
    actor2 = c(2, 3, 4, 1, 1, 3, 2, 4, 1, 5, 4, 1, 2, 3, 5, 1)
  )
  reh   <- remify(el, model = "tie", ordinal = FALSE)
  # several varying effects so the design comfortably keeps >= 2 penalisable
  # statistics even if one happens to be degenerate and is dropped
  stats <- remstats(reh, tie_effects = ~ inertia() + reciprocity() +
                                          indegreeSender() + outdegreeReceiver())

  set.seed(1)
  fit_l <- remstimate(reh, stats, penalty = list(alpha = 1))

  expect_inherits(fit_l, "remstimate_glmnet",   info = "glmnet: class")
  expect_true("baseline" %in% names(fit_l$coefficients),
              info = "glmnet: baseline retained in coefficients")
  expect_true("baseline" %in% fit_l$unpenalized,
              info = "glmnet: baseline reported as unpenalised")
  # the baseline is the intercept (~ log base rate); it must NOT be shrunk to 0
  expect_true(abs(fit_l$coefficients[["baseline"]]) > 1e-6,
              info = "glmnet: baseline is a real level, not shrunk to 0")

  # override plumbing is reflected in the reported unpenalised set
  set.seed(1)
  fit_add <- remstimate(reh, stats, penalty = list(alpha = 1, unpenalized = "inertia"))
  expect_true("inertia" %in% fit_add$unpenalized,
              info = "glmnet: unpenalized ADDS 'inertia' to the exemption")

  set.seed(1)
  fit_sub <- remstimate(reh, stats, penalty = list(alpha = 1, penalized = "inertia"))
  expect_false("inertia" %in% fit_sub$unpenalized,
              info = "glmnet: penalized keeps 'inertia' in the penalty")

  # diagnostics route through the shared MLE machinery and carry an 'event' index
  diag_l <- diagnostics(fit_l, reh, stats = stats)
  expect_inherits(diag_l, "diagnostics",              info = "glmnet: diagnostics class")
  expect_false(is.null(diag_l$recall),                info = "glmnet: recall present")
  expect_true("event" %in% names(diag_l$recall$per_event),
              info = "glmnet: recall carries an 'event' column for plot.diagnostics")

  # plot.diagnostics interface must run without error (render to a throwaway device)
  pf <- tempfile(fileext = ".pdf"); grDevices::pdf(pf)
  expect_true(runs_ok(plot(fit_l, reh, stats = stats, which = 3)),
              info = "glmnet: recall plot (which = 3) runs")
  expect_true(runs_ok(plot(fit_l, reh, stats = stats, which = 9)),
              info = "glmnet: regularisation-path plot (which = 9) runs")
  grDevices::dev.off(); unlink(pf)
}

# ── 4. GLMNET duration diagnostics route through the durem recall machinery ───
#    Regression guard: a durem fit must NOT hit diagnostics.remstimate()'s
#    tomstats/aomstats check; it goes through .point_durem_diagnostics instead.

if (at_home() && requireNamespace("glmnet", quietly = TRUE)) {

  el_d <- data.frame(
    time   = c(1, 2, 5, 8, 11, 14),
    actor1 = c("A", "B", "A", "C", "B", "A"),
    actor2 = c("B", "C", "C", "A", "A", "B"),
    end    = c(3, 4, 7, 10, 13, 16)
  )
  ok_d <- tryCatch({
    suppressWarnings({
      reh_d   <- remify(el_d, duration = TRUE)
      stats_d <- remstats(reh_d,
                          start_effects = ~ 1 + inertia() + reciprocity(scaling = "std"),
                          end_effects   = ~ 1 + inertia(),
                          psi_start = 1, psi_end = 1)
      set.seed(1)
      fit_d   <- remstimate(reh_d, stats_d, penalty = list(alpha = 1))
    })
    TRUE
  }, error = function(e) FALSE)

  if (ok_d && inherits(fit_d, "remstimate_glmnet")) {
    diag_d <- diagnostics(fit_d, reh_d, stats_d)
    expect_inherits(diag_d, "diagnostics_durem",
                    info = "glmnet durem: diagnostics route to durem machinery, not tie")
    expect_false(is.null(diag_d$recall_joint),
                 info = "glmnet durem: joint start/end recall present")

    pf <- tempfile(fileext = ".pdf"); grDevices::pdf(pf)
    expect_true(runs_ok(plot(fit_d, reh_d, stats = stats_d, which = 1)),
                info = "glmnet durem: recall plot routes to plot.diagnostics_durem")
    grDevices::dev.off(); unlink(pf)
  }
}

# ── 5. MIXREM diagnostics route through plot.diagnostics (slow: at_home only) ──

if (at_home() && requireNamespace("flexmix", quietly = TRUE)) {

  el2 <- data.frame(
    time   = 1:14,
    actor1 = c(1, 2, 1, 3, 2, 1, 3, 2, 4, 1, 2, 3, 4, 1),
    actor2 = c(2, 3, 4, 1, 1, 3, 2, 4, 1, 5, 4, 1, 2, 3)
  )
  reh2   <- remify(el2, model = "tie", ordinal = FALSE)
  stats2 <- remstats(reh2, tie_effects = ~ inertia() + reciprocity())

  fit_mx <- tryCatch(
    remstimate(reh2, stats2, mixture = list(k = 2, random = ~ (1 + inertia | dyad))),
    error = function(e) NULL)

  if (!is.null(fit_mx) && inherits(fit_mx, "remstimate_mixrem")) {
    diag_mx <- diagnostics(fit_mx, reh2, stats = stats2)
    expect_inherits(diag_mx, "diagnostics_mixrem",  info = "mixrem: diagnostics class")
    expect_inherits(diag_mx, "diagnostics",         info = "mixrem: routes through plot.diagnostics")
    expect_false(is.null(diag_mx$recall),           info = "mixrem: posterior-weighted recall present")
    expect_true("event" %in% names(diag_mx$recall$per_event),
                info = "mixrem: recall carries an 'event' column")
    expect_false(is.null(diag_mx$recall_by_component),
                info = "mixrem: per-component recall present")
    expect_false(is.null(diag_mx$.reh.processed),
                info = "mixrem: .reh.processed stored for plot.diagnostics")

    pf <- tempfile(fileext = ".pdf"); grDevices::pdf(pf)
  }
}
