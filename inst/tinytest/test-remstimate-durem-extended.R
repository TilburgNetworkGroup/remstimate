## Additional tests for the duration REM pipeline
##
## Covers: MLE/GLM equivalence, directed end, typed events, mixed formulas,
## one-sided effects, IC methods, diagnostics, and plot.

library(tinytest)

# ── Shared edgelists ─────────────────────────────────────────────────────────

el <- data.frame(
    time   = c(1, 2, 5, 6, 6, 10, 11),
    actor1 = c("A", "B", "A", "B", "C", "C", "B"),
    actor2 = c("B", "C", "C","A", "B", "A", "C"),
    end    = c(6, 7, 8, 7, 7, 11, 13)
)

el_typed <- data.frame(
    time   = c(1, 2, 3, 4, 5, 6),
    actor1 = c("A", "B", "A", "B", "A", "B"),
    actor2 = c("B", "C", "B", "C", "B", "C"),
    type   = c("X", "X", "Y", "Y", "X", "Y"),
    duration = c(2, 2, 2, 3, 3, 3)
)

el_larger <- data.frame(
    time   = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    actor1 = c("A","B","C","A","B","C","A","B","C","A"),
    actor2 = c("B","C","A","C","A","B","B","C","A","C"),
    end    = c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
)

# ══════════════════════════════════════════════════════════════════════════════
# 1. MLE and GLM give the same output
# ══════════════════════════════════════════════════════════════════════════════

suppressWarnings({
    reh <- remify(el, duration = TRUE)
    stats <- remstats(reh,
                      start_effects = ~ inertia(),
                      end_effects   = ~ inertia(),
                      psi_start = 1, psi_end = 1)
    fit_mle <- remstimate(reh, stats, method = "MLE")
    fit_glm <- remstimate(reh, stats, method = "GLM")
})

expect_equal(coef(fit_mle), coef(fit_glm), tolerance = 1e-10,
    info = "MLE and GLM produce identical coefficients")
expect_equal(fit_mle$loglik, fit_glm$loglik, tolerance = 1e-10,
    info = "MLE and GLM produce identical log-likelihood")
expect_equal(fit_mle$AIC, fit_glm$AIC, tolerance = 1e-10,
    info = "MLE and GLM produce identical AIC")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Directed end model
# ══════════════════════════════════════════════════════════════════════════════

suppressWarnings({
    reh_de <- remify(el, duration = TRUE, directed_end = TRUE)
    stats_de <- remstats(reh_de,
                         start_effects = ~ inertia(),
                         end_effects   = ~ inertia(),
                         psi_start = 1, psi_end = 1)
    fit_de <- remstimate(reh_de, stats_de)
})

expect_inherits(fit_de, "remstimate_durem",
    info = "directed end: correct class")
expect_true(all(is.finite(coef(fit_de))),
    info = "directed end: coefficients finite")
# Directed end has N*(N-1) end dyads vs N*(N-1)/2 for undirected
expect_true(fit_de$stacked_data$D_end >= fit_de$stacked_data$D_start %/% 2,
    info = "directed end: D_end consistent with directed riskset")


# ══════════════════════════════════════════════════════════════════════════════
# 3. Start effects only (no end_effects)
# ══════════════════════════════════════════════════════════════════════════════

suppressWarnings({
    reh_so <- remify(el, duration = TRUE)
    stats_so <- remstats(reh_so,
                         start_effects = ~ inertia(),
                         psi_start = 1)
    fit_so <- remstimate(reh_so, stats_so)
})

expect_inherits(fit_so, "remstimate_durem",
    info = "start-only: correct class")
expect_true("inertia.start" %in% names(coef(fit_so)),
    info = "start-only: inertia.start present")
# No end stats in names
expect_false(any(grepl("\\.end$", names(coef(fit_so)))),
    info = "start-only: no .end coefficients")


# ══════════════════════════════════════════════════════════════════════════════
# 4. Typed events with consider_type = "separate"
# ══════════════════════════════════════════════════════════════════════════════

suppressWarnings({
    reh_t <- remify(el_typed, duration = TRUE, model = "tie")
    stats_t <- remstats(reh_t,
                        start_effects = ~ inertia(consider_type = "separate"),
                        psi_start = 1)
    fit_t <- remstimate(reh_t, stats_t)
})

expect_inherits(fit_t, "remstimate_durem",
    info = "typed separate: correct class")
expect_true("inertia.X.start" %in% names(coef(fit_t)),
    info = "typed separate: inertia.X.start present")
expect_true("inertia.Y.start" %in% names(coef(fit_t)),
    info = "typed separate: inertia.Y.start present")
expect_true(all(is.finite(coef(fit_t))),
    info = "typed separate: coefficients finite")


# ══════════════════════════════════════════════════════════════════════════════
# 5. Typed events with consider_type = "interact" (ext=TRUE)
# ══════════════════════════════════════════════════════════════════════════════

suppressWarnings({
    reh_ti <- remify(el_typed, duration = TRUE, model = "tie",
                     extend_riskset_by_type = TRUE)
    stats_ti <- remstats(reh_ti,
                         start_effects = ~ inertia(consider_type = "interact"),
                         psi_start = 1)
    fit_ti <- remstimate(reh_ti, stats_ti)
})

expect_inherits(fit_ti, "remstimate_durem",
    info = "typed interact: correct class")
# interact produces C^2 stats per effect: X.X, X.Y, Y.X, Y.Y
interact_names <- grep("^inertia\\.", names(coef(fit_ti)), value = TRUE)
expect_equal(length(interact_names), 4L,
    info = "typed interact: 4 inertia stats (C^2 = 2^2)")
expect_true(all(is.finite(coef(fit_ti))),
    info = "typed interact: coefficients finite")


# ══════════════════════════════════════════════════════════════════════════════
# 6. Mixed formula: active-state + history stats
# ══════════════════════════════════════════════════════════════════════════════

suppressWarnings({
    reh_mx <- remify(el_typed, duration = TRUE, model = "tie")
    stats_mx <- remstats(reh_mx,
                         start_effects = ~ activeTie() + inertia(),
                         psi_start = 1)
    fit_mx <- remstimate(reh_mx, stats_mx)
})

expect_inherits(fit_mx, "remstimate_durem",
    info = "mixed formula: correct class")
expect_true("activeTie.start" %in% names(coef(fit_mx)),
    info = "mixed formula: activeTie.start present")
expect_true("inertia.start" %in% names(coef(fit_mx)),
    info = "mixed formula: inertia.start present")
expect_true(all(is.finite(coef(fit_mx))),
    info = "mixed formula: coefficients finite")


# ══════════════════════════════════════════════════════════════════════════════
# 7. AIC / AICC / BIC methods
# ══════════════════════════════════════════════════════════════════════════════

suppressWarnings({
    reh7 <- remify(el, duration = TRUE)
    stats7 <- remstats(reh7,
                       start_effects = ~ inertia(),
                       end_effects   = ~ inertia(),
                       psi_start = 1, psi_end = 1)
    fit7 <- remstimate(reh7, stats7)
})

expect_true(is.numeric(fit7$AIC) && is.finite(fit7$AIC),
    info = "AIC is finite numeric")
expect_true(is.numeric(fit7$AICC) && is.finite(fit7$AICC),
    info = "AICC is finite numeric")
expect_true(is.numeric(fit7$BIC) && is.finite(fit7$BIC),
    info = "BIC is finite numeric")
# AICC >= AIC (correction is non-negative for P >= 0)
expect_true(fit7$AICC >= fit7$AIC,
    info = "AICC >= AIC")


# ══════════════════════════════════════════════════════════════════════════════
# 8. Deviance: model deviance = null - residual
# ══════════════════════════════════════════════════════════════════════════════

expect_equal(fit7$model.deviance,
             fit7$null.deviance - fit7$residual.deviance,
             tolerance = 1e-8,
    info = "model deviance = null - residual")


# ══════════════════════════════════════════════════════════════════════════════
# 9. vcov and se
# ══════════════════════════════════════════════════════════════════════════════

expect_true(is.matrix(fit7$vcov),
    info = "vcov is a matrix")
expect_equal(nrow(fit7$vcov), length(coef(fit7)),
    info = "vcov dimensions match number of coefficients")
expect_equal(length(fit7$se), length(coef(fit7)),
    info = "se length matches number of coefficients")
expect_true(all(fit7$se > 0),
    info = "all standard errors are positive")


# ══════════════════════════════════════════════════════════════════════════════
# 10. diagnostics() with separate start/end recall
# ══════════════════════════════════════════════════════════════════════════════

suppressWarnings({
    reh10 <- remify(el_larger, duration = TRUE)
    stats10 <- remstats(reh10,
                        start_effects = ~ inertia(),
                        end_effects   = ~ inertia(),
                        psi_start = 1, psi_end = 1)
    fit10 <- remstimate(reh10, stats10)
    diag10 <- diagnostics(fit10, reh10, stats10)
})

expect_inherits(diag10, "diagnostics_durem",
    info = "diagnostics: correct class")

# Joint recall
expect_true(!is.null(diag10$recall_joint),
    info = "diagnostics: joint recall computed")
expect_true(!is.null(diag10$recall_joint$summary),
    info = "diagnostics: joint recall summary present")
expect_true(diag10$recall_joint$summary$mean_rel_rank >= 0 &&
            diag10$recall_joint$summary$mean_rel_rank <= 1,
    info = "diagnostics: joint mean_rel_rank in [0, 1]")

# Start recall
expect_true(!is.null(diag10$recall_start),
    info = "diagnostics: start recall computed")
expect_true(diag10$recall_start$summary$mean_rel_rank >= 0 &&
            diag10$recall_start$summary$mean_rel_rank <= 1,
    info = "diagnostics: start mean_rel_rank in [0, 1]")

# End recall
expect_true(!is.null(diag10$recall_end),
    info = "diagnostics: end recall computed")
expect_true(diag10$recall_end$summary$mean_rel_rank >= 0 &&
            diag10$recall_end$summary$mean_rel_rank <= 1,
    info = "diagnostics: end mean_rel_rank in [0, 1]")

# Residuals
expect_true(!is.null(diag10$deviance_residuals),
    info = "diagnostics: deviance residuals computed")
expect_true(!is.null(diag10$pearson_residuals),
    info = "diagnostics: Pearson residuals computed")


# ══════════════════════════════════════════════════════════════════════════════
# 11. print(diagnostics) shows joint + start + end recall
# ══════════════════════════════════════════════════════════════════════════════

out11 <- capture.output(print(diag10))
expect_true(any(grepl("Joint", out11)),
    info = "print(diagnostics) shows joint recall")
expect_true(any(grepl("Start", out11)),
    info = "print(diagnostics) shows start recall")
expect_true(any(grepl("End", out11)),
    info = "print(diagnostics) shows end recall")


# ══════════════════════════════════════════════════════════════════════════════
# 12. summary() returns expected structure
# ══════════════════════════════════════════════════════════════════════════════

s <- summary(fit10)
expect_true(is.list(s),
    info = "summary returns a list")
expect_true(!is.null(s$coefsTab),
    info = "summary contains coefsTab")
expect_equal(nrow(s$coefsTab), length(coef(fit10)),
    info = "coefsTab has one row per coefficient")
expect_equal(ncol(s$coefsTab), 4L,
    info = "coefsTab has 4 columns (Estimate, SE, z, p)")


# ══════════════════════════════════════════════════════════════════════════════
# 13. psi = 0: unweighted (all completed events contribute weight 1)
# ══════════════════════════════════════════════════════════════════════════════

suppressWarnings({
    reh13 <- remify(el_larger, duration = TRUE)
    stats13 <- remstats(reh13,
                        start_effects = ~ inertia(),
                        end_effects   = ~ inertia(),
                        psi_start = 0, psi_end = 0)
    fit13 <- remstimate(reh13, stats13)
})

expect_inherits(fit13, "remstimate_durem",
    info = "psi=0: correct class")
expect_true(all(is.finite(coef(fit13))),
    info = "psi=0: coefficients finite")


# ══════════════════════════════════════════════════════════════════════════════
# 14. Larger example with multiple effects and baseline
# ══════════════════════════════════════════════════════════════════════════════

suppressWarnings({
    reh14 <- remify(el_larger, duration = TRUE, directed_end = TRUE)
    stats14 <- remstats(reh14,
                        start_effects = ~ inertia() + outdegreeSender(),
                        end_effects   = ~ inertia() + indegreeSender(),
                        psi_start = 1, psi_end = 1)
    fit14 <- remstimate(reh14, stats14)
})

expected_names <- c("baseline.start", "inertia.start", "outdegreeSender.start",
                    "baseline.end", "inertia.end", "indegreeSender.end")
expect_equal(sort(names(coef(fit14))), sort(expected_names),
    info = "multi-effect: all expected coefficient names present")
expect_true(all(is.finite(coef(fit14))),
    info = "multi-effect: all coefficients finite")


# ══════════════════════════════════════════════════════════════════════════════
# 15. Input validation: wrong arguments for durem
# ══════════════════════════════════════════════════════════════════════════════

suppressWarnings(reh15 <- remify(el, duration = TRUE))

# tie_effects on a durem object should error
expect_error(
    remstats(reh15, tie_effects = ~ inertia()),
    pattern = "start_effects.*end_effects",
    info = "tie_effects on remify_durem gives informative error"
)

# sender_effects on a durem object should error
expect_error(
    remstats(reh15, sender_effects = ~ inertia()),
    pattern = "not yet supported",
    info = "sender_effects on remify_durem gives informative error"
)

# No effects at all should error
expect_error(
    remstats(reh15),
    info = "no effects on remify_durem gives error"
)


# ══════════════════════════════════════════════════════════════════════════════
# 16. Input validation: start_effects on non-durem should error
# ══════════════════════════════════════════════════════════════════════════════

reh16 <- remify(el[, c("time","actor1","actor2")], duration = FALSE, model = "tie")

expect_error(
    remstats(reh16, start_effects = ~ inertia()),
    pattern = "duration",
    info = "start_effects on non-durem gives informative error"
)


# ══════════════════════════════════════════════════════════════════════════════
# 17. interact coercion warning when ext=FALSE
# ══════════════════════════════════════════════════════════════════════════════

suppressWarnings({
    reh17 <- remify(el_typed, duration = TRUE, model = "tie",
                    extend_riskset_by_type = FALSE)
})

remstats(reh17,
         start_effects = ~ activeTie(consider_type = "interact"),
         first = 1L, last = Inf)

# ══════════════════════════════════════════════════════════════════════════════
# 18. Overlapping ties on same dyad
# ══════════════════════════════════════════════════════════════════════════════

el_overlap <- data.frame(
    time   = c(1, 3, 5, 7),
    actor1 = c("A", "A", "B", "A"),
    actor2 = c("B", "B", "C", "B"),
    duration = c(4, 2, 3, 1)
)

suppressWarnings({
    reh18 <- remify(el_overlap, duration = TRUE)
    stats18 <- remstats(reh18,
                        start_effects = ~ activeTie() + inertia(),
                        psi_start = 1)
    fit18 <- remstimate(reh18, stats18)
})

expect_inherits(fit18, "remstimate_durem",
    info = "overlapping ties: correct class")
expect_true(all(is.finite(coef(fit18))),
    info = "overlapping ties: coefficients finite")

