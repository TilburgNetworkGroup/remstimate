# ── Duration REM estimation backends ─────────────────────────────────────────
#
# All durem estimation goes through stack_stats() first, then dispatches to
# an R backend.  The C++ MLE/HMC path is not used because the durem riskset
# is time-varying (dyads toggle between start-risk and end-risk).
#
# Backends:
#   .remstimate_durem_glm   — Poisson GLM (interval) or clogit (ordinal)
#   .remstimate_durem_brms  — Bayesian via brms (interval only; HMC on the stacked Poisson)
#
# Called from remstimate() when reh inherits "remify_durem".


# ── Poisson GLM / clogit backend (MLE) ────────────────────────────────────

#' @importFrom survival strata clogit coxph
#' @keywords internal
.remstimate_durem_glm <- function(stacked) {

    df         <- stacked$remstats_stack
    stat_names <- stacked$stat_names
    ordinal    <- isTRUE(stacked$ordinal)

    if (length(stat_names) == 0L)
        stop("No statistics found — check start_effects / end_effects.",
             call. = FALSE)

    if (ordinal) {
        # ── Ordinal: conditional logistic regression ──────────────────────
        if (!requireNamespace("survival", quietly = TRUE))
            stop("Package 'survival' is required for ordinal duration models. ",
                 "Install it with install.packages('survival').", call. = FALSE)

        formula_obj <- stats::as.formula(paste0(
            "obs ~ ",
            paste(stat_names, collapse = " + "),
            " + strata(time_index)"
        ))

        fit <- survival::clogit(formula_obj, data = df)

        coefs      <- stats::coef(fit)
        loglik_val <- as.numeric(fit$loglik[2L])
        null_loglik <- as.numeric(fit$loglik[1L])
        P          <- length(coefs)
        M          <- stacked$E

        vcov_mat <- tryCatch(
            stats::vcov(fit),
            error = function(e) matrix(NA_real_, P, P)
        )
        se <- sqrt(diag(vcov_mat))
        names(se) <- names(coefs)

        aic_val  <- stats::AIC(fit)
        bic_val  <- stats::BIC(fit)
        aicc_val <- aic_val + 2 * P * (P + 1) / max(M - P - 1, 1)

    } else {
        # ── Interval: Poisson GLM ────────────────────────────────────────
        formula_obj <- stats::as.formula(paste0(
            "obs ~ -1 + offset(log_interevent) + ",
            paste(stat_names, collapse = " + ")
        ))

        fit <- stats::glm(formula_obj, family = stats::poisson(), data = df)

        coefs      <- stats::coef(fit)
        loglik_val <- as.numeric(stats::logLik(fit))
        P          <- length(coefs)
        M          <- stacked$E

        vcov_mat <- tryCatch(
            stats::vcov(fit),
            error = function(e) matrix(NA_real_, P, P)
        )
        se <- sqrt(diag(vcov_mat))
        names(se) <- names(coefs)

        null_formula <- stats::as.formula("obs ~ -1 + offset(log_interevent)")
        null_fit     <- stats::glm(null_formula, family = stats::poisson(), data = df)
        null_loglik  <- as.numeric(stats::logLik(null_fit))

        aic_val  <- stats::AIC(fit)
        bic_val  <- stats::BIC(fit)
        aicc_val <- aic_val + 2 * P * (P + 1) / max(M - P - 1, 1)
    }

    where_bl <- which(tolower(stat_names) %in% c("baseline.start", "baseline"))
    where_bl <- if (length(where_bl)) where_bl[1L] else NULL

    res <- list(
        coefficients      = coefs,
        loglik            = loglik_val,
        gradient          = NULL,
        hessian           = NULL,
        vcov              = vcov_mat,
        se                = se,
        residual.deviance = -2 * loglik_val,
        null.deviance     = -2 * null_loglik,
        model.deviance    = -2 * null_loglik - (-2 * loglik_val),
        AIC               = aic_val,
        AICC              = aicc_val,
        BIC               = bic_val,
        df.null           = M,
        df.model          = P,
        df.residual       = M - P,
        converged         = if (ordinal) !is.null(coefs) else fit$converged,
        iterations        = if (ordinal) fit$iter else fit$iter,
        stacked_data      = stacked,
        backend_fit       = fit
    )

    structure(
        res,
        class             = c("remstimate_durem", "remstimate"),
        formula           = formula_obj,
        model             = "tie",
        approach          = "Frequentist",
        method            = "MLE",
        engine            = if (ordinal) "clogit" else "glm",
        ordinal           = ordinal,
        statistics        = stat_names,
        where_is_baseline = where_bl,
        ncores            = 1L
    )
}


# ── brms backend (HMC) ──────────────────────────────────────────────────────

#' @keywords internal
.remstimate_durem_brms <- function(stacked,
                                   prior   = NULL,
                                   nsim    = 500L,
                                   nchains = 1L,
                                   burnin  = 500L,
                                   thin    = 10L,
                                   seed    = NULL,
                                   ...) {

    if (!requireNamespace("brms", quietly = TRUE))
        stop("Package 'brms' is required for HMC estimation of duration models. ",
             "Install it with install.packages('brms').", call. = FALSE)

    df         <- stacked$remstats_stack
    stat_names <- stacked$stat_names

    if (isTRUE(stacked$ordinal))
        stop("Ordinal + HMC is not yet supported for duration models. ",
             "Use method = 'MLE' for ordinal durem (clogit backend).",
             call. = FALSE)

    if (length(stat_names) == 0L)
        stop("No statistics found — check start_effects / end_effects.",
             call. = FALSE)

    formula_obj <- stats::as.formula(paste0(
        "obs ~ -1 + offset(log_interevent) + ",
        paste(stat_names, collapse = " + ")
    ))

    if (is.null(prior))
        prior <- brms::set_prior("normal(0, 10)", class = "b")

    fit <- brms::brm(
        formula  = formula_obj,
        family   = stats::poisson(),
        data     = df,
        prior    = prior,
        iter     = nsim + burnin,
        warmup   = burnin,
        chains   = nchains,
        thin     = thin,
        seed     = seed,
        ...
    )

    draws <- brms::as_draws_matrix(fit)
    b_cols <- grep("^b_", colnames(draws))
    coef_draws <- as.matrix(draws[, b_cols])
    colnames(coef_draws) <- stat_names

    coefs    <- colMeans(coef_draws)
    vcov_mat <- stats::cov(coef_draws)
    sd_vec   <- sqrt(diag(vcov_mat))
    P        <- length(coefs)
    M        <- stacked$E

    names(coefs) <- names(sd_vec) <-
        rownames(vcov_mat) <- colnames(vcov_mat) <- stat_names

    res <- list(
        coefficients = coefs,
        post.mean    = coefs,
        vcov         = vcov_mat,
        sd           = sd_vec,
        loglik       = NULL,
        draws        = coef_draws,
        df.null      = M,
        df.model     = P,
        df.residual  = M - P,
        stacked_data = stacked,
        backend_fit  = fit
    )

    structure(
        res,
        class      = c("remstimate_durem", "remstimate"),
        formula    = formula_obj,
        model      = "tie",
        approach   = "Bayesian",
        method     = "HMC",
        engine     = "brms",
        ordinal    = FALSE,
        statistics = stat_names,
        ncores     = 1L
    )
}


# ── S3 methods ──────────────────────────────────────────────────────────────

#' @export
#' @method print remstimate_durem
print.remstimate_durem <- function(x, ...) {

    approach <- attr(x, "approach")
    engine   <- attr(x, "engine")
    method   <- attr(x, "method")

    cat("Relational Event Model (tie oriented, duration)\n")
    cat("Estimation:", method, "[", engine, "]\n")

    if (approach == "Frequentist") {
        cat("\nCoefficients:\n\n")
        print(x$coefficients)
        cat("\nNull deviance:    ", x$null.deviance, "\n")
        cat("Residual deviance:", x$residual.deviance, "\n")
        cat("AIC:", x$AIC, " AICC:", x$AICC, " BIC:", x$BIC, "\n\n")
    } else {
        cat("\nPosterior means:\n\n")
        print(x$post.mean)
        if (!is.null(x$sd)) {
            cat("\nPosterior SD:\n\n")
            print(x$sd)
        }
        cat("\n")
    }
    invisible(x)
}

#' @export
#' @method summary remstimate_durem
summary.remstimate_durem <- function(object, ...) {

    approach   <- attr(object, "approach")
    engine     <- attr(object, "engine")
    method     <- attr(object, "method")
    stat_names <- attr(object, "statistics")

    summary_out <- list()
    summary_out$approach <- approach
    summary_out$method   <- method
    summary_out$engine   <- engine

    cat("Relational Event Model (tie oriented, duration)\n")
    cat("Estimation:", method, "[", engine, "]\n")
    cat(strrep("-", 50), "\n")

    if (approach == "Frequentist") {
        coefs <- object$coefficients
        se    <- object$se
        z     <- coefs / se
        p     <- 2 * (1 - stats::pnorm(abs(z)))

        coefsTab <- cbind(coefs, se, z, p)
        colnames(coefsTab) <- c("Estimate", "Std. Err", "z value", "Pr(>|z|)")
        rownames(coefsTab) <- stat_names
        summary_out$coefsTab <- coefsTab

        cat("\nCoefficients:\n")
        stats::printCoefmat(coefsTab, signif.stars = TRUE)

        cat("\nNull deviance:    ", object$null.deviance, " on", object$df.null,
            "degrees of freedom\n")
        cat("Residual deviance:", object$residual.deviance, " on", object$df.residual,
            "degrees of freedom\n")
        cat("Chi-square:", object$model.deviance, "on", object$df.model,
            "degrees of freedom, p-value",
            1 - stats::pchisq(object$model.deviance, object$df.model), "\n")
        cat("AIC:", object$AIC, " AICC:", object$AICC, " BIC:", object$BIC, "\n")

        summary_out$AIC               <- object$AIC
        summary_out$AICC              <- object$AICC
        summary_out$BIC               <- object$BIC
        summary_out$null.deviance     <- object$null.deviance
        summary_out$residual.deviance <- object$residual.deviance
        summary_out$model.deviance    <- object$model.deviance
        summary_out$loglik            <- object$loglik

    } else {
        draws <- object$draws
        coefs <- object$post.mean
        sd    <- object$sd
        q025  <- apply(draws, 2, stats::quantile, 0.025)
        q975  <- apply(draws, 2, stats::quantile, 0.975)

        coefsTab <- cbind(coefs, sd, q025, q975)
        colnames(coefsTab) <- c("Post. Mean", "Post. SD", "2.5%", "97.5%")
        rownames(coefsTab) <- stat_names
        summary_out$coefsTab <- coefsTab

        cat("\nCoefficients:\n")
        print(round(coefsTab, 4))
        cat("\n")
    }

    invisible(summary_out)
}

#' @export
#' @method coef remstimate_durem
coef.remstimate_durem <- function(object, ...) object$coefficients

#' @export
#' @method logLik remstimate_durem
logLik.remstimate_durem <- function(object, ...) object$loglik


# ── diagnostics ─────────────────────────────────────────────────────────────

#' @export
#' @method diagnostics remstimate_durem
diagnostics.remstimate_durem <- function(object, reh, stats, top_pct = 0.05, ...) {

    approach <- attr(object, "approach")
    if (approach == "Bayesian")
        warning("Diagnostics for Bayesian duration models use posterior means.",
                call. = FALSE)

    fit <- object$backend_fit
    out <- list()

    # ── Deviance and Pearson residuals ──────────────────────────────────────
    if (inherits(fit, "glm")) {
        out$deviance_residuals <- stats::residuals(fit, type = "deviance")
        out$pearson_residuals  <- stats::residuals(fit, type = "pearson")
        dev_res <- out$deviance_residuals
        out$residual_summary <- data.frame(
            min    = min(dev_res),
            q1     = stats::quantile(dev_res, 0.25),
            median = stats::median(dev_res),
            q3     = stats::quantile(dev_res, 0.75),
            max    = max(dev_res),
            row.names = NULL
        )
    } else if (inherits(fit, "coxph")) {
        out$deviance_residuals <- stats::residuals(fit, type = "deviance")
        dev_res <- out$deviance_residuals
        out$residual_summary <- data.frame(
            min    = min(dev_res),
            q1     = stats::quantile(dev_res, 0.25),
            median = stats::median(dev_res),
            q3     = stats::quantile(dev_res, 0.75),
            max    = max(dev_res),
            row.names = NULL
        )
    }

    # ── Recall ───────────────────────────────────────────────────────────────
    stacked <- object$stacked_data
    if (!is.null(stacked)) {
        df         <- stacked$remstats_stack
        stat_names <- stacked$stat_names
        sn_start   <- stacked$stat_names_start
        sn_end     <- stacked$stat_names_end
        coefs      <- object$coefficients

        X  <- as.matrix(df[, stat_names, drop = FALSE])
        lp <- as.numeric(X %*% coefs[stat_names])

        obs_idx   <- which(df$obs == 1L)
        event_ids <- df$time

        # Identify start vs end rows by nonzero start/end stats.
        # Start rows have baseline.start = 1 (or at least one start stat != 0)
        # and all end stats == 0; vice versa for end rows.
        is_start <- rep(FALSE, nrow(df))
        is_end   <- rep(FALSE, nrow(df))
        if (!is.null(sn_start) && length(sn_start) > 0L)
            is_start <- rowSums(abs(df[, sn_start, drop = FALSE])) > 0
        if (!is.null(sn_end) && length(sn_end) > 0L)
            is_end <- rowSums(abs(df[, sn_end, drop = FALSE])) > 0

        # ── Joint recall: all competing dyads ────────────────────────────────
        out$recall_joint <- .recall_block(lp, obs_idx, event_ids, top_pct)

        # ── Start recall: observed starts ranked among start-risk dyads ──────
        start_obs <- intersect(obs_idx, which(is_start))
        if (length(start_obs) > 0L) {
            # Only keep start-risk rows for the events that have an observed start
            start_events <- unique(event_ids[start_obs])
            start_mask   <- which(is_start & event_ids %in% start_events)
            out$recall_start <- .recall_block(
                lp[start_mask], match(start_obs, start_mask),
                event_ids[start_mask], top_pct
            )
        }

        # ── End recall: observed ends ranked among end-risk dyads ────────────
        end_obs <- intersect(obs_idx, which(is_end))
        if (length(end_obs) > 0L) {
            end_events <- unique(event_ids[end_obs])
            end_mask   <- which(is_end & event_ids %in% end_events)
            out$recall_end <- .recall_block(
                lp[end_mask], match(end_obs, end_mask),
                event_ids[end_mask], top_pct
            )
        }

        # ── Per-type recall (ext=TRUE) ──────────────────────────────────────
        if ("type" %in% names(df) && !all(is.na(df$type))) {
            types <- sort(unique(df$type[!is.na(df$type)]))
            if (length(types) > 1L) {
                out$recall_by_type <- list()
                for (tp in types) {
                    tp_mask <- which(df$type == tp)
                    tp_obs  <- intersect(obs_idx, tp_mask)
                    if (length(tp_obs) == 0L) next
                    out$recall_by_type[[tp]] <- .recall_block(
                        lp[tp_mask], match(tp_obs, tp_mask),
                        event_ids[tp_mask], top_pct
                    )
                }

                # Per-type × start/end
                if (length(start_obs) > 0L) {
                    out$recall_start_by_type <- list()
                    for (tp in types) {
                        tp_start <- which(is_start & df$type == tp)
                        tp_obs_s <- intersect(obs_idx, tp_start)
                        if (length(tp_obs_s) == 0L) next
                        tp_events <- unique(event_ids[tp_obs_s])
                        tp_mask   <- which(is_start & df$type == tp &
                                             event_ids %in% tp_events)
                        out$recall_start_by_type[[tp]] <- .recall_block(
                            lp[tp_mask], match(tp_obs_s, tp_mask),
                            event_ids[tp_mask], top_pct
                        )
                    }
                }

                if (length(end_obs) > 0L) {
                    out$recall_end_by_type <- list()
                    for (tp in types) {
                        tp_end   <- which(is_end & df$type == tp)
                        tp_obs_e <- intersect(obs_idx, tp_end)
                        if (length(tp_obs_e) == 0L) next
                        tp_events <- unique(event_ids[tp_obs_e])
                        tp_mask   <- which(is_end & df$type == tp &
                                             event_ids %in% tp_events)
                        out$recall_end_by_type[[tp]] <- .recall_block(
                            lp[tp_mask], match(tp_obs_e, tp_mask),
                            event_ids[tp_mask], top_pct
                        )
                    }
                }
            }
        }
    }

    out$.reh.processed <- reh
    class(out) <- c("diagnostics_durem", "diagnostics", "remstimate")
    out
}


# ── print.diagnostics_durem ─────────────────────────────────────────────────

#' @export
#' @method print diagnostics_durem
print.diagnostics_durem <- function(x, ...) {
    reh <- x$.reh.processed
    cat("Diagnostics for a Relational Event Model (duration)\n")
    cat(strrep("-", 50), "\n")
    if (!is.null(reh$N)) cat(sprintf("%-12s: %d\n", "Actors", reh$N))

    if (!is.null(x$residual_summary)) {
        cat("\nDeviance residuals:\n")
        print(x$residual_summary, row.names = FALSE)
    }

    if (!is.null(x$recall_joint) || !is.null(x$recall_start) || !is.null(x$recall_end)) {
        cat("\nRecall:\n")
        .print_recall_summary(x$recall_joint, "Joint")
        .print_recall_summary(x$recall_start, "Start")
        .print_recall_summary(x$recall_end,   "End")
    }

    if (!is.null(x$recall_by_type)) {
        cat("\nRecall by type:\n")
        for (tp in names(x$recall_by_type))
            .print_recall_summary(x$recall_by_type[[tp]], paste0("  ", tp))
    }
    if (!is.null(x$recall_start_by_type)) {
        cat("  Start by type:\n")
        for (tp in names(x$recall_start_by_type))
            .print_recall_summary(x$recall_start_by_type[[tp]], paste0("    ", tp))
    }
    if (!is.null(x$recall_end_by_type)) {
        cat("  End by type:\n")
        for (tp in names(x$recall_end_by_type))
            .print_recall_summary(x$recall_end_by_type[[tp]], paste0("    ", tp))
    }

    invisible(x)
}


# ── plot ────────────────────────────────────────────────────────────────────

#' @export
#' @method plot remstimate_durem
plot.remstimate_durem <- function(x, reh = NULL, stats = NULL,
                                  diagnostics_object = NULL,
                                  which = 1L, ...) {

    if (is.null(diagnostics_object) && !is.null(reh) && !is.null(stats))
        diagnostics_object <- diagnostics(x, reh, stats)

    if (is.null(diagnostics_object)) {
        # Fallback: coefficient CI plot
        coefs <- x$coefficients
        se    <- x$se
        if (is.null(se)) {
            graphics::barplot(coefs, main = "Coefficients", las = 2)
            return(invisible(x))
        }
        ci_lo <- coefs - 1.96 * se
        ci_hi <- coefs + 1.96 * se
        ord   <- seq_along(coefs)
        graphics::plot(ord, coefs, ylim = range(c(ci_lo, ci_hi)),
                       xaxt = "n", xlab = "", ylab = "Estimate",
                       main = "Coefficients with 95% CI", pch = 19)
        graphics::segments(ord, ci_lo, ord, ci_hi)
        graphics::abline(h = 0, lty = 2, col = "grey")
        graphics::axis(1, at = ord, labels = names(coefs), las = 2)
        return(invisible(x))
    }

    if (which == 1L) {
        # Joint recall
        .plot_recall_scatter(diagnostics_object$recall_joint, "Recall: joint (all competing dyads)")
    } else if (which == 2L) {
        # Deviance residuals
        if (!is.null(diagnostics_object$deviance_residuals)) {
            graphics::plot(diagnostics_object$deviance_residuals,
                           xlab = "Observation", ylab = "Deviance residual",
                           main = "Deviance residuals", pch = ".",
                           col = grDevices::adjustcolor("black", 0.3), ...)
            graphics::abline(h = 0, lty = 2, col = "grey")
        }
    } else if (which == 3L) {
        # Start recall
        .plot_recall_scatter(diagnostics_object$recall_start, "Recall: start events")
    } else if (which == 4L) {
        # End recall
        .plot_recall_scatter(diagnostics_object$recall_end, "Recall: end events")
    } else if (which == 5L) {
        # All recalls side by side
        has_start <- !is.null(diagnostics_object$recall_start)
        has_end   <- !is.null(diagnostics_object$recall_end)
        n_panels  <- 1L + has_start + has_end
        old_par   <- graphics::par(mfrow = c(1, n_panels))
        on.exit(graphics::par(old_par))
        .plot_recall_scatter(diagnostics_object$recall_joint, "Joint")
        if (has_start) .plot_recall_scatter(diagnostics_object$recall_start, "Start")
        if (has_end)   .plot_recall_scatter(diagnostics_object$recall_end,   "End")
    } else if (which == 6L) {
        # Per-type recall panels
        rbt <- diagnostics_object$recall_by_type
        if (!is.null(rbt) && length(rbt) > 0L) {
            n_types  <- length(rbt)
            old_par  <- graphics::par(mfrow = c(1, n_types))
            on.exit(graphics::par(old_par))
            for (tp in names(rbt))
                .plot_recall_scatter(rbt[[tp]], paste("Type:", tp))
        }
    } else if (which == 7L) {
        # Per-type start recall panels
        rbt <- diagnostics_object$recall_start_by_type
        if (!is.null(rbt) && length(rbt) > 0L) {
            n_types  <- length(rbt)
            old_par  <- graphics::par(mfrow = c(1, n_types))
            on.exit(graphics::par(old_par))
            for (tp in names(rbt))
                .plot_recall_scatter(rbt[[tp]], paste("Start type:", tp))
        }
    } else if (which == 8L) {
        # Per-type end recall panels
        rbt <- diagnostics_object$recall_end_by_type
        if (!is.null(rbt) && length(rbt) > 0L) {
            n_types  <- length(rbt)
            old_par  <- graphics::par(mfrow = c(1, n_types))
            on.exit(graphics::par(old_par))
            for (tp in names(rbt))
                .plot_recall_scatter(rbt[[tp]], paste("End type:", tp))
        }
    } else if (which == 9L) {
        # Probability ratio: joint
        .plot_probratio_scatter(diagnostics_object$recall_joint,
                                "Prob ratio: joint", ...)
    } else if (which == 10L) {
        # Probability ratio: start + end side by side
        has_start <- !is.null(diagnostics_object$recall_start)
        has_end   <- !is.null(diagnostics_object$recall_end)
        n_panels  <- 1L + has_start + has_end
        old_par   <- graphics::par(mfrow = c(1, n_panels))
        on.exit(graphics::par(old_par))
        .plot_probratio_scatter(diagnostics_object$recall_joint, "Joint")
        if (has_start) .plot_probratio_scatter(diagnostics_object$recall_start, "Start")
        if (has_end)   .plot_probratio_scatter(diagnostics_object$recall_end, "End")
    }

    invisible(x)
}
