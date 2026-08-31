# ── Duration REM estimation backends ─────────────────────────────────────────
#
# All durem estimation goes through stack_stats() first, then dispatches to
# an R backend.  The C++ MLE/HMC path is not used because the durem riskset
# is time-varying (dyads toggle between start-risk and end-risk).
#
# Backends:
#   .remstimate_glm         - Poisson GLM (interval) or clogit (ordinal)
#
# Called from remstimate() when reh inherits "remify_durem".


# ── clogit helper ─────────────────────────────────────────────────────────
#
# 'survival' is only a Suggested dependency, so it is never attached to the
# search path. survival::clogit() is implemented by rewriting the problem as a
# bare coxph() call with a Surv(...) response and a strata() term, then
# eval()'ing that call in the caller's frame / the formula's environment. With
# survival unattached, none of those bare names resolve, so clogit fails on
# whichever it reaches first:
#   Error in strata(time_index) : could not find function "strata"
#   Error in coxph(...)          : could not find function "coxph"
#
# Whether it fails is platform-dependent: it only works when some other package
# (coxme, lme4, ...) happened to attach 'survival' first, which is why R-hub
# Windows passed but macOS / clean CRAN checks fail.
#
# Fix: bind the survival functions clogit relies on in this frame, and point the
# formula's environment here, so every bare name resolves during both clogit's
# own eval() and coxph's model.frame() - on every platform, without attaching
# survival. The bare names (not survival::) are kept in the formula so coxph's
# "specials" detection, which matches on the unqualified symbol, still
# recognises strata().
#
#' @keywords internal
.durem_clogit <- function(formula_obj, df) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("Package 'survival' is required for ordinal (clogit) duration models. ",
         "Install it with install.packages('survival').", call. = FALSE)

  coxph  <- survival::coxph
  Surv   <- survival::Surv
  strata <- survival::strata
  clogit <- survival::clogit
  environment(formula_obj) <- environment()

  clogit(formula_obj, data = df)
}


# ── Poisson GLM / clogit backend (MLE) ────────────────────────────────────

#' @keywords internal
.remstimate_glm <- function(stacked) {

  df         <- stacked$remstats_stack
  stat_names <- stacked$stat_names
  ordinal    <- isTRUE(stacked$ordinal)
  # This backend also serves plain pre-stacked tie models, which are not
  # duration models; the flag drives the print/summary header.
  is_durem   <- inherits(stacked, "remstats_stacked_durem") ||
                identical(stacked$model, "durem")

  # case-control sampling correction: offset += log(weight) = -log(pi_d)
  # cases have weight 1 -> 0; controls weight 1/pi -> -log(pi). Absent => 0.
  df$.samp_off <- if ("weight" %in% names(df)) log(df$weight) else 0

    if (length(stat_names) == 0L)
        stop("No statistics found - check start_effects / end_effects.",
             call. = FALSE)

    if (ordinal) {
        # ── Ordinal: conditional logistic regression ──────────────────────
        formula_obj <- stats::as.formula(paste0(
          "obs ~ ",
          paste(stat_names, collapse = " + "),
          " + strata(time_index) + offset(.samp_off)"
        ))

        fit <- .durem_clogit(formula_obj, df)

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
        "obs ~ -1 + offset(log_interevent + .samp_off) + ",
        paste(stat_names, collapse = " + ")
      ))
        fit <- stats::glm(formula_obj, family = stats::poisson(), data = df)

        coefs <- stats::coef(fit)
        P     <- length(coefs)
        M     <- stacked$E

        # Poisson loglik -> REM loglik: drop the offset constant Σ log(dt) over events
        off_obs    <- sum(df$log_interevent[df$obs == 1L])
        loglik_val <- as.numeric(stats::logLik(fit)) - off_obs

        vcov_mat <- tryCatch(stats::vcov(fit),
                             error = function(e) matrix(NA_real_, P, P))
        se <- sqrt(diag(vcov_mat)); names(se) <- names(coefs)

        # baseline-intercept null (matches C++), same offset correction
        null_fit    <- stats::glm(obs ~ offset(log_interevent + .samp_off),
                                  family = stats::poisson(), data = df)
        null_loglik <- as.numeric(stats::logLik(null_fit)) - off_obs

        aic_val  <- 2 * P - 2 * loglik_val
        bic_val  <- P * log(M) - 2 * loglik_val            # n = E, not nrow
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
        durem             = is_durem,
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

.stacked_block_out <- function(fit, ll, ll0, P, E, eng) {
  cf <- stats::coef(fit)
  vc <- tryCatch(stats::vcov(fit), error = function(e) matrix(NA_real_, P, P))
  se <- sqrt(diag(vc)); names(se) <- names(cf)
  list(coefficients = cf, loglik = ll, se = se, vcov = vc,
       residual.deviance = -2*ll, null.deviance = -2*ll0,
       model.deviance = (-2*ll0) - (-2*ll),
       AIC = 2*P - 2*ll, BIC = P*log(E) - 2*ll,
       AICC = (2*P - 2*ll) + 2*P*(P+1)/max(E-P-1, 1),
       df.null = E, df.model = P, df.residual = E - P,
       converged = TRUE, backend_fit = fit)
}

.fit_rate_block <- function(df, sn, E, ordinal) {
  P <- length(sn)
  if (ordinal) {
    f <- stats::as.formula(paste0("obs ~ ", paste(sn, collapse=" + "), " + strata(time_index)"))
    fit <- .durem_clogit(f, df)
    .stacked_block_out(fit, fit$loglik[2L], fit$loglik[1L], P, E, "clogit")
  } else {
    f   <- stats::as.formula(paste0("obs ~ -1 + offset(log_interevent) + ", paste(sn, collapse=" + ")))
    fit <- stats::glm(f, family = stats::poisson(), data = df)
    off <- sum(df$log_interevent[df$obs == 1L])
    ll  <- as.numeric(stats::logLik(fit)) - off
    ll0 <- as.numeric(stats::logLik(
      stats::glm(obs ~ offset(log_interevent), family = stats::poisson(), data = df))) - off
    .stacked_block_out(fit, ll, ll0, P, E, "glm")
  }
}

.fit_choice_block <- function(df, sn, E) {            # always clogit
  P <- length(sn)
  f <- stats::as.formula(paste0("obs ~ ", paste(sn, collapse=" + "), " + strata(time_index)"))
  fit <- .durem_clogit(f, df)
  .stacked_block_out(fit, fit$loglik[2L], fit$loglik[1L], P, E, "clogit")
}

.remstimate_stacked_glm_actor <- function(stacked) {
  E <- stacked$E; ordinal <- isTRUE(stacked$ordinal)
  out <- list(
    sender_model   = .fit_rate_block(stacked$sender_stack,   stacked$sender_stat_names,   E, ordinal),
    receiver_model = .fit_choice_block(stacked$receiver_stack, stacked$receiver_stat_names, E)
  )
  structure(out, class = "remstimate",
            model = "actor", ordinal = ordinal, method = "MLE", approach = "Frequentist",
            formula = list(
              rate_model_formula   = stats::reformulate(stacked$sender_stat_names),
              choice_model_formula = stats::reformulate(stacked$receiver_stat_names)),
            statistics = list(sender_model   = stacked$sender_stat_names,
                              receiver_model = stacked$receiver_stat_names),
            ncores = 1L)
}


# ── S3 methods ──────────────────────────────────────────────────────────────

#' @export
#' @method print remstimate_durem
print.remstimate_durem <- function(x, ...) {

    approach <- attr(x, "approach")
    engine   <- attr(x, "engine")
    method   <- attr(x, "method")

    cat(if (isFALSE(attr(x, "durem")))
            "Relational Event Model (tie oriented)\n"
        else "Relational Event Model (tie oriented, duration)\n")
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

    cat(if (isFALSE(attr(object, "durem")))
            "Relational Event Model (tie oriented)\n"
        else "Relational Event Model (tie oriented, duration)\n")
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
logLik.remstimate_durem <- function(object, ...) {
  structure(object$loglik, class = "logLik",
            df = length(object$coefficients),
            nobs = object$df.null)
}

# ── diagnostics ─────────────────────────────────────────────────────────────

# Shared durem recall/surprise computation. Given a linear predictor `lp`
# aligned to the rows of the stacked frame `df`, the start/end statistic names,
# and the remify object, build the joint / start / end recall tables, their
# surprises and offender tables, and the per-type breakdowns. Used by both the
# MLE path (lp = X %*% coef) and the GLMM path (lp = predict(fit), incl. BLUPs),
# so the two never drift apart.
.durem_recall_tables <- function(lp, df, reh, sn_start, sn_end,
                                 top_pct = 0.05, surprise_threshold = 0.2) {
    out <- list()

    obs_idx   <- which(df$obs == 1L)
    event_ids <- df$time_index
    # Per-row dyad identity comes straight from the stacked design (resolved
    # per-process in stack_stats: `incl` for starts, `rs_end` for ends). The old
    # reh$edgelist_dual$.eidx[df$time_index] round-trip indexed by epoch ordinal,
    # which is shared by every observed event in a multi-event epoch (~1/3 of
    # events), so it mislabelled them and inflated per-dyad end/joint counts.
    # .recall_block picks ids[obs] per observed row, so a per-row label vector
    # attributes every event correctly.
    if (is.null(df$actor1) || is.null(df$actor2))
        stop("stacked design lacks actor1/actor2 columns; ",
             "re-stack with add_actors = TRUE.", call. = FALSE)
    df$eidx   <- paste(df$actor1, "->", df$actor2)

    # Identify start vs end rows by nonzero start/end stats. Intersect with the
    # columns actually present: a degenerate start/end statistic may have been
    # dropped from the (GLMM) design by .drop_constant_stats.
    sn_start <- intersect(sn_start, colnames(df))
    sn_end   <- intersect(sn_end,   colnames(df))
    # Start/end membership is exact in the `process` column set by stack_stats;
    # prefer it over inferring the process from nonzero stats (which is fragile
    # if the baseline/start/end intercept is ever dropped from the design).
    if (!is.null(df$process)) {
        is_start <- df$process == "start"
        is_end   <- df$process == "end"
    } else {
        is_start <- if (length(sn_start) > 0L) rowSums(abs(df[, sn_start, drop = FALSE])) > 0 else rep(FALSE, nrow(df))
        is_end   <- if (length(sn_end)   > 0L) rowSums(abs(df[, sn_end,   drop = FALSE])) > 0 else rep(FALSE, nrow(df))
    }

    # ── Joint recall: all competing dyads ────────────────────────────────
    out$recall_joint <- .recall_block(lp, obs_idx, event_ids, top_pct, ids = df$eidx)

    # ── Start recall: observed starts ranked among start-risk dyads ──────
    start_obs <- intersect(obs_idx, which(is_start))
    if (length(start_obs) > 0L) {
        start_events <- unique(event_ids[start_obs])
        start_mask   <- which(is_start & event_ids %in% start_events)
        out$recall_start <- .recall_block(
            lp[start_mask], match(start_obs, start_mask),
            event_ids[start_mask], top_pct, ids = df$eidx[start_mask]
        )
    }

    # ── End recall: observed ends ranked among end-risk dyads ────────────
    # min_D_t = 2: single-candidate end epochs carry no diagnostic signal.
    end_obs <- intersect(obs_idx, which(is_end))
    if (length(end_obs) > 0L) {
      end_events <- unique(event_ids[end_obs])
      end_mask   <- which(is_end & event_ids %in% end_events)
      out$recall_end <- .recall_block(
        lp[end_mask], match(end_obs, end_mask),
        event_ids[end_mask], top_pct, ids = df$eidx[end_mask],
        min_D_t = 2
      )
    }

    # ── Surprises: most poorly-predicted observed events ─────────────────
    out$surprises_joint <- .surprises_from_recall(out$recall_joint, surprise_threshold)
    out$surprises_start <- .surprises_from_recall(out$recall_start, surprise_threshold)
    out$surprises_end   <- .surprises_from_recall(out$recall_end,   surprise_threshold)

    # `eidx` now carries the "actor1 -> actor2" label directly (see above), so
    # feed it straight to .offender_table -- no .durem_dyad_labels round-trip.
    out$surprise_offenders_joint <- .offender_table(
      ids_surprises = out$surprises_joint$eidx,
      ids_all       = out$recall_joint$per_event$eidx
    )
    if (!is.null(out$recall_start))
      out$surprise_offenders_start <- .offender_table(
        ids_surprises = out$surprises_start$eidx,
        ids_all       = out$recall_start$per_event$eidx
      )
    if (!is.null(out$recall_end))
      out$surprise_offenders_end <- .offender_table(
        ids_surprises = out$surprises_end$eidx,
        ids_all       = out$recall_end$per_event$eidx
      )
    out$surprise_threshold <- surprise_threshold

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
                    event_ids[tp_mask], top_pct, ids = df$eidx[tp_mask]
                )
            }

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
                        event_ids[tp_mask], top_pct, ids = df$eidx[tp_mask]
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
                  event_ids[tp_mask], top_pct, ids = df$eidx[tp_mask],
                  min_D_t = 2
                )
              }
            }

            surprise_lists <- list(
              surprises_by_type       = out$recall_by_type,
              surprises_start_by_type = out$recall_start_by_type,
              surprises_end_by_type   = out$recall_end_by_type
            )
            for (nm in names(surprise_lists)) {
              rc_list <- surprise_lists[[nm]]
              if (is.null(rc_list) || length(rc_list) == 0L) next
              out[[nm]] <- lapply(rc_list, .surprises_from_recall,
                                  threshold = surprise_threshold)
            }
        }
    }

    out
}

#' @export
#' @method diagnostics remstimate_durem
diagnostics.remstimate_durem <- function(object, reh, stats, top_pct = 0.05,
                                         surprise_threshold = 0.2, ...) {

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
        coefs      <- object$coefficients

        X  <- as.matrix(df[, stat_names, drop = FALSE])
        lp <- as.numeric(X %*% coefs[stat_names])

        out <- c(out, .durem_recall_tables(
            lp, df, reh,
            sn_start = stacked$stat_names_start,
            sn_end   = stacked$stat_names_end,
            top_pct  = top_pct,
            surprise_threshold = surprise_threshold))
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
      .print_recall_summary(x$recall_joint, "Joint", x$surprises_joint, x$surprise_threshold)
      .print_recall_summary(x$recall_start, "Start", x$surprises_start, x$surprise_threshold)
      .print_recall_summary(x$recall_end,   "End",   x$surprises_end,   x$surprise_threshold)
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


# ── summary.diagnostics_durem ────────────────────────────────────────────────
# print.diagnostics_durem gives the one-line-per-recall overview; summary() adds
# the full per-recall table (mean AND median rank, cum prob, prob ratio, log
# loss, top-pct) for joint/start/end, mirroring summary.diagnostics for the
# tie/actor models. Without this method summary(diag_durem) fell through to
# summary.diagnostics, which looks for $recall / $residuals in the tie layout the
# durem object does not use, and printed almost nothing.

#' @export
#' @method summary diagnostics_durem
summary.diagnostics_durem <- function(object, ...) {
    reh <- object$.reh.processed

    cat("Diagnostics for a Relational Event Model (duration)\n")
    cat(strrep("-", 50), "\n", sep = "")
    if (!is.null(reh$N)) cat(sprintf("%-12s: %d\n", "Actors", reh$N))

    if (!is.null(object$residual_summary)) {
        cat("\nDeviance residuals:\n")
        print(object$residual_summary, row.names = FALSE)
    }

    # one summary row per recall table
    .row <- function(rc, label) {
        if (is.null(rc)) return(NULL)
        rs <- rc$summary
        data.frame(
            recall          = label,
            n_events        = nrow(rc$per_event),
            mean_rel_rank   = round(rs$mean_rel_rank,   4),
            median_rel_rank = round(rs$median_rel_rank, 4),
            mean_cum_prob   = round(rs$mean_cum_prob,   4),
            mean_prob_ratio = round(rs$mean_prob_ratio, 4),
            mean_log_loss   = round(rs$mean_log_loss,   4),
            top_pct         = rs$top_pct,
            top_pct_prop    = round(rs$top_pct_prop,    4),
            row.names = NULL, stringsAsFactors = FALSE
        )
    }

    rec <- do.call(rbind, list(
        .row(object$recall_joint, "joint"),
        .row(object$recall_start, "start"),
        .row(object$recall_end,   "end")
    ))
    if (!is.null(rec)) {
        cat("\nRecall:\n")
        print(rec, row.names = FALSE)
    }

    # per-type (joint) recall, when the design carries multiple event types
    if (!is.null(object$recall_by_type) && length(object$recall_by_type)) {
        bt <- do.call(rbind, lapply(names(object$recall_by_type), function(tp)
            .row(object$recall_by_type[[tp]], paste0("type:", tp))))
        if (!is.null(bt)) {
            cat("\nRecall by type (joint):\n")
            print(bt, row.names = FALSE)
        }
    }

    invisible(object)
}


# ── plot ────────────────────────────────────────────────────────────────────

#' @export
#' @method plot diagnostics_durem
plot.diagnostics_durem <- function(x, which = 5L, ...) {
    # x is the durem diagnostics object (recall_joint / recall_start /
    # recall_end / *_by_type / deviance_residuals). Plotting straight from the
    # diagnostics object mirrors plot.diagnostics for tie/actor, so plot(diag)
    # works uniformly across model types.
    diagnostics_object <- x

    if (which == 1L) {
        # Joint recall
        .plot_recall_scatter(diagnostics_object$recall_joint, "Recall: start and end events combined")
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
    } else if (which == 11L) {
        # Random-effect normality Q-Q (GLMM duration models only)
        re_list <- diagnostics_object$ranef
        if (!is.null(re_list)) {
            panels <- list()
            for (nm in names(re_list)) {
                re <- re_list[[nm]]
                if (is.null(re)) next
                for (g in names(re)) {
                    v <- suppressWarnings(as.numeric(unlist(re[[g]], use.names = FALSE)))
                    v <- v[is.finite(v)]
                    if (length(v) < 2L) next
                    lbl <- if (length(re_list) > 1L) paste0(nm, ": ", g) else g
                    panels[[lbl]] <- v
                }
            }
            if (length(panels)) {
                nc      <- min(length(panels), 3L)
                old_par <- graphics::par(mfrow = c(ceiling(length(panels) / nc), nc))
                on.exit(graphics::par(old_par))
                for (p in seq_along(panels)) {
                    stats::qqnorm(panels[[p]], main = "", pch = 16, cex = 0.8)
                    stats::qqline(panels[[p]], col = 2, lwd = 2)
                    graphics::mtext(paste("Random effects -", names(panels)[p]),
                                    side = 3, line = 1, cex = 1.1)
                }
            }
        }
    }

    invisible(x)
}

#' @export
#' @method plot remstimate_durem
plot.remstimate_durem <- function(x, reh = NULL, stats = NULL,
                                  diagnostics_object = NULL,
                                  which = 5L, ...) {

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

    # Delegate to plot.diagnostics_durem so plot(fit, ...) and plot(diagnostics)
    # share one code path.
    plot.diagnostics_durem(diagnostics_object, which = which, ...)
    invisible(x)
}
