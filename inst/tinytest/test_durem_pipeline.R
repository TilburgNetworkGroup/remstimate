# ── Smoke tests: durem pipeline (remstats + remstimate) ──────────────────────
# Checks that the full pipeline runs without errors for key configurations.
# Not a correctness check — just a "does it blow up" guard.

library(tinytest)

# ── Test data ────────────────────────────────────────────────────────────────

el_simple <- data.frame(
  time     = c(1, 2, 3, 4, 5, 6),
  actor1   = c("A", "B", "A", "B", "A", "B"),
  actor2   = c("B", "C", "B", "C", "B", "C"),
  duration = c(2, 2, 2, 3, 3, 3)
)

el_typed <- data.frame(
  time     = c(1, 2, 3, 4, 6, 6),
  actor1   = c("A", "B", "A", "B", "A", "B"),
  actor2   = c("B", "C", "B", "C", "B", "C"),
  type     = c("X", "X", "Y", "Y", "X", "Y"),
  duration = c(2, 2, 2, 3, 3, 3)
)


# ── 1. Basic durem, interval timing ─────────────────────────────────────────

reh <- remify(el_simple, duration = TRUE, model = "tie")
stats <- remstats(reh, start_effects = ~ inertia(), psi_start = 1)
stacked <- stack_stats(stats, reh)

expect_true(inherits(stacked, "remstats_stacked"),
            info = "stack_stats returns remstats_stacked")
expect_true("time_index" %in% colnames(stacked$remstats_stack),
            info = "stacked output has 'time' column")
expect_true(!"event" %in% colnames(stacked$remstats_stack),
            info = "stacked output no longer has 'event' column")
expect_true(nrow(stacked$remstats_stack) > 0L,
            info = "stacked output has rows")

fit <- remstimate(reh, stats, method = "MLE")
expect_true(inherits(fit, "remstimate_durem"),
            info = "MLE fit returns remstimate_durem")
expect_true(!is.null(fit$coefficients),
            info = "MLE fit has coefficients")


# ── 2. Basic durem, ordinal timing ──────────────────────────────────────────

reh_ord <- remify(el_simple, duration = TRUE, model = "tie", ordinal = TRUE)
stats_ord <- remstats(reh_ord, start_effects = ~ inertia(), psi_start = 1)
stacked_ord <- stack_stats(stats_ord, reh_ord)

expect_true(!"log_interevent" %in% colnames(stacked_ord$remstats_stack),
            info = "ordinal stacked output has no log_interevent")

fit_ord <- remstimate(reh_ord, stats_ord, method = "MLE")
expect_true(inherits(fit_ord, "remstimate_durem"),
            info = "ordinal MLE fit returns remstimate_durem")
expect_equal(attr(fit_ord, "engine"), "clogit",
             info = "ordinal uses clogit engine")


# ── 3. Typed durem, ext=TRUE, interval ──────────────────────────────────────

reh_te <- remify(el_typed, duration = TRUE, extend_riskset_by_type = TRUE,
                 riskset = "active", model = "tie")
stats_te <- remstats(reh_te,
                     start_effects = ~ inertia(consider_type = "interact"),
                     psi_start = 1)
stacked_te <- stack_stats(stats_te, reh_te)

expect_true("type" %in% colnames(stacked_te$remstats_stack),
            info = "ext=TRUE stacked output has 'type' column")
expect_true(nrow(stacked_te$remstats_stack) > 0L,
            info = "ext=TRUE stacked output has rows")

fit_te <- remstimate(reh_te, stats_te, method = "MLE")
expect_true(!is.null(fit_te$coefficients),
            info = "ext=TRUE MLE fit has coefficients")


# ── 4. type_exclusive = TRUE vs FALSE ────────────────────────────────────────

reh_excl <- remify(el_typed, duration = TRUE, extend_riskset_by_type = TRUE,
                   riskset = "active", model = "tie", type_exclusive = TRUE)
stats_excl <- remstats(reh_excl,
                       start_effects = ~ inertia(consider_type = "interact"),
                       psi_start = 1)
stacked_excl <- stack_stats(stats_excl, reh_excl)

reh_indep <- remify(el_typed, duration = TRUE, extend_riskset_by_type = TRUE,
                    riskset = "active", model = "tie", type_exclusive = FALSE)
stats_indep <- remstats(reh_indep,
                        start_effects = ~ inertia(consider_type = "interact"),
                        psi_start = 1)
stacked_indep <- stack_stats(stats_indep, reh_indep)

expect_true(nrow(stacked_excl$remstats_stack) <
              nrow(stacked_indep$remstats_stack),
            info = "type_exclusive=TRUE produces fewer rows than FALSE")


# ── 5. Typed durem, ext=TRUE, ordinal ────────────────────────────────────────

reh_te_ord <- remify(el_typed, duration = TRUE, extend_riskset_by_type = TRUE,
                     riskset = "active", model = "tie", ordinal = TRUE)
stats_te_ord <- remstats(reh_te_ord,
                         start_effects = ~ inertia(consider_type = "interact"),
                         psi_start = 1)
stacked_te_ord <- stack_stats(stats_te_ord, reh_te_ord)

expect_true(nrow(stacked_te_ord$remstats_stack) > 0L,
            info = "ordinal + ext=TRUE stacked output has rows")

fit_te_ord <- remstimate(reh_te_ord, stats_te_ord, method = "MLE")
expect_true(!is.null(fit_te_ord$coefficients),
            info = "ordinal + ext=TRUE MLE fit has coefficients")


# ── 6. Diagnostics: basic durem ─────────────────────────────────────────────

diag <- diagnostics(fit, reh, stats)
expect_true(inherits(diag, "diagnostics_durem"),
            info = "diagnostics returns diagnostics_durem")
expect_true(!is.null(diag$recall_joint),
            info = "diagnostics has joint recall")


# ── 7. Diagnostics: typed durem (per-type recall) ───────────────────────────

diag_te <- diagnostics(fit_te, reh_te, stats_te)
expect_true(!is.null(diag_te$recall_joint),
            info = "ext=TRUE diagnostics has joint recall")
expect_true(!is.null(diag_te$recall_by_type),
            info = "ext=TRUE diagnostics has per-type recall")
expect_true(length(diag_te$recall_by_type) > 0L,
            info = "per-type recall has entries")


# ── 8. No duplicate dyads within same time point ────────────────────────────

df <- stacked_te$remstats_stack
for (tp in unique(df$time)) {
  dyads_at_tp <- df$dyad[df$time == tp]
  expect_true(length(dyads_at_tp) == length(unique(dyads_at_tp)),
              info = paste("no duplicate dyads at time", tp))
}


# ── 9. Time column matches subset sequence ──────────────────────────────────

sub <- attr(stats, "subset")
first <- sub[["first"]]
last  <- sub[["last"]]
time_vals <- sort(unique(stacked$remstats_stack$time))
expect_true(min(time_vals) >= first,
            info = "time values >= first from subset")
expect_true(max(time_vals) <= last,
            info = "time values <= last from subset")


# ── 10. Print methods run without error ─────────────────────────────────────

expect_silent(capture.output(print(diag)),
              info = "print.diagnostics_durem runs for basic durem")
expect_silent(capture.output(print(diag_te)),
              info = "print.diagnostics_durem runs for ext=TRUE")
expect_silent(capture.output(summary(fit)),
              info = "summary.remstimate_durem runs for interval")
expect_silent(capture.output(summary(fit_ord)),
              info = "summary.remstimate_durem runs for ordinal")


# ── 11. GLMNET + durem diagnostics ──────────────────────────────────────────

if (requireNamespace("glmnet", quietly = TRUE)) {
  fit_glmnet <- remstimate(reh, stats, method = "GLMNET")
  expect_true(inherits(fit_glmnet, "remstimate"),
              info = "GLMNET durem fit works")

  diag_glmnet <- diagnostics(fit_glmnet)
  expect_true(!is.null(diag_glmnet$recall),
              info = "GLMNET durem diagnostics has recall")
  expect_silent(capture.output(print(diag_glmnet)),
                info = "print diagnostics_remstimate runs for GLMNET")
}


# ── 12. GLMM + durem diagnostics ────────────────────────────────────────────

if (requireNamespace("lme4", quietly = TRUE)) {
  fit_glmm <- remstimate(reh, stats, method = "GLMM",
                          random = ~ (1 | actor1))
  expect_true(inherits(fit_glmm, "remstimate"),
              info = "GLMM durem fit works")

  diag_glmm <- diagnostics(fit_glmm)
  expect_true(!is.null(diag_glmm$recall),
              info = "GLMM durem diagnostics has recall")
  expect_true(isTRUE(diag_glmm$use_ranef),
              info = "GLMM diagnostics uses random effects")
}


# ── 13. MIXREM + durem diagnostics ──────────────────────────────────────────

if (requireNamespace("flexmix", quietly = TRUE)) {
  fit_mix <- remstimate(reh, stats, method = "MIXREM",
                         random = ~ (1 | dyad), k = 2)
  expect_true(inherits(fit_mix, "remstimate"),
              info = "MIXREM durem fit works")

  diag_mix <- diagnostics(fit_mix)
  expect_true(!is.null(diag_mix$recall),
              info = "MIXREM durem diagnostics has recall")
  expect_true(!is.null(diag_mix$recall_by_component),
              info = "MIXREM diagnostics has per-component recall")
  expect_silent(capture.output(print(diag_mix)),
                info = "print diagnostics_mixrem runs")
}


# ── 14. GLMNET + typed durem diagnostics ────────────────────────────────────

if (requireNamespace("glmnet", quietly = TRUE)) {
  # GLMNET uses add_actors=FALSE, so no type column — no per-type recall
  fit_glmnet_te <- remstimate(reh_te, stats_te, method = "GLMNET")
  diag_glmnet_te <- diagnostics(fit_glmnet_te)
  expect_true(!is.null(diag_glmnet_te$recall),
              info = "GLMNET typed durem diagnostics has recall")
}

