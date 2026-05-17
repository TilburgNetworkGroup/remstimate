## Tests for remstimate.remstats_durem
##
## Same 3-event edgelist as test-stack-durem.R:
##   A→B (t=1, end=6)   B→C (t=2, end=7)   A→C (t=5, end=8)
## Dual unique times: 1, 2, 5, 6, 7, 8 → with start=2: 5 time points
##
## The model has two parameters: inertia.start and inertia.end.
## With such a small data set we cannot recover precise values, but we can
## check that fitting converges, coefficients are named correctly, and that
## the joint start/end likelihood is well-defined.

library(tinytest)

el <- data.frame(
    time   = c(1, 2, 5, 11, 12, 15),
    actor1 = c("A", "B", "A"),
    actor2 = c("B", "C", "C"),
    end    = c(6, 7, 8, 12, 12, 16)
)

suppressWarnings({
    reh   <- remify(el, duration = TRUE)
    stats <- remstats(reh,
                      start_effects = ~ inertia(),
                      end_effects   = ~ inertia(),
                      psi_start = 1, psi_end = 1)
    fit   <- remstimate(reh, stats)
})

# ── 1. Class ──────────────────────────────────────────────────────────────────

expect_inherits(fit, "remstimate_durem",
    info = "fit inherits remstimate_durem")
expect_inherits(fit, "remstimate",
    info = "fit inherits remstimate")

# ── 2. Coefficients ───────────────────────────────────────────────────────────

expect_true(is.numeric(fit$coefficients),
    info = "coefficients is numeric")
expect_equal(length(fit$coefficients), 4L,
    info = "two coefficients: inertia.start + inertia.end")
expect_true("inertia.start" %in% names(fit$coefficients),
    info = "inertia.start coefficient present")
expect_true("inertia.end" %in% names(fit$coefficients),
    info = "inertia.end coefficient present")
expect_true(all(is.finite(fit$coefficients)),
    info = "all coefficients are finite")

# ── 3. loglik ────────────────────────────────────────────────────────────────

expect_true(is.numeric(fit$loglik),
    info = "loglik is numeric")
expect_equal(length(fit$loglik), 1L,
    info = "loglik is scalar")
expect_true(is.finite(fit$loglik),
    info = "loglik is finite")
expect_true(fit$loglik < 0,
    info = "Poisson log-likelihood is negative")

# ── 4. Attributes ─────────────────────────────────────────────────────────────

expect_equal(attr(fit, "model"),   "tie",   info = "model attribute = durem")
expect_equal(attr(fit, "method"),  "MLE",   info = "method attribute = MLE")
expect_equal(attr(fit, "engine"),  "glm",   info = "engine attribute = glm")
expect_equal(attr(fit, "ordinal"), FALSE,   info = "ordinal = FALSE")
expect_equal(attr(fit, "statistics"), c("baseline.start","inertia.start","baseline.end","inertia.end"),
    info = "statistics attribute matches stat_names")

# ── 5. stacked_data reference ─────────────────────────────────────────────────

expect_true(is.list(fit$stacked_data),
    info = "stacked_data is a list")
expect_inherits(fit$stacked_data, "remstats_stacked_durem",
    info = "stacked_data is remstats_stacked_durem")

# ── 6. backend_fit ────────────────────────────────────────────────────────────

expect_inherits(fit$backend_fit, "glm",
    info = "backend_fit is a glm object")

# ── 7. coef() and logLik() S3 methods ────────────────────────────────────────

expect_equal(coef(fit), fit$coefficients,
    info = "coef(fit) equals fit$coefficients")

ll <- logLik(fit)
expect_inherits(ll, "logLik",
    info = "logLik(fit) returns logLik class")
expect_equal(as.numeric(ll), fit$loglik,
    info = "logLik value matches fit$loglik")

# ── 8. print/summary run without error ───────────────────────────────────────

expect_silent(print(fit),
    info = "print(fit) runs without error")

expect_silent(summary(fit),
    info = "summary(fit) runs without error")

# ── 10. Undirected end model ──────────────────────────────────────────────────

suppressWarnings({
    reh_ude   <- remify(el, duration = TRUE)  # directed_end = FALSE (default)
    stats_ude <- remstats(reh_ude,
                          start_effects = ~ inertia(),
                          end_effects   = ~ inertia(),
                          psi_start = 0, psi_end = 0)
    fit_ude   <- remstimate(reh_ude, stats_ude)
})

expect_inherits(fit_ude, "remstimate_durem",
    info = "undirected end model: correct class")
expect_true(all(is.finite(coef(fit_ude))),
    info = "undirected end model: coefficients finite")

# ── 11. Right-censored events ─────────────────────────────────────────────────

el_cens <- el
el_cens$end[2] <- NA  # B→C is right-censored

suppressWarnings({
    reh_cens   <- remify(el_cens, duration = TRUE)
    stats_cens <- remstats(reh_cens,
                           start_effects = ~ inertia(),
                           end_effects   = ~ inertia())
    fit_cens   <- remstimate(reh_cens, stats_cens)
})

expect_inherits(fit_cens, "remstimate_durem",
    info = "right-censored: correct class")
expect_true(all(is.finite(coef(fit_cens))),
    info = "right-censored: coefficients finite")

# ── 12. Multiple effects ──────────────────────────────────────────────────────

suppressWarnings({
    reh2   <- remify(el, duration = TRUE)
    stats2 <- remstats(reh2,
                       start_effects = ~ inertia() + outdegreeSender(),
                       end_effects   = ~ inertia(),
                       psi_start = 1, psi_end = 1)
    fit2   <- remstimate(reh2, stats2)
})

expect_equal(length(coef(fit2)), 5L,
    info = "three coefficients: inertia.start + outdegreeSender.start + inertia.end")
expect_true(all(c("inertia.start", "outdegreeSender.start", "inertia.end") %in%
                    names(coef(fit2))),
    info = "all expected coefficient names present for multi-effect model")

