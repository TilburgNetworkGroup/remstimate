# GLMM backend — lme4::glmer or glmmTMB::glmmTMB
#
# remstimate(reh, stats, method = "GLMM", random = ~ (1|actor1), engine = "lme4")

.remstimate_glmm <- function(reh, stats,
                             random  = NULL,
                             engine  = c("lme4", "glmmTMB"),
                             verbose = FALSE,
                             ...) {
  engine <- match.arg(engine)

  s <- .remstimate_make_stack(reh, stats, add_actors = TRUE)

  # Ordinal models need coxme (conditional logit + random effects)
  if (s$ordinal) {
    engine <- "coxme"
  }

  if (!requireNamespace(engine, quietly = TRUE))
    stop("install.packages('", engine, "')")
  if (is.null(random))
    stop("provide a random formula, e.g. ~ (1 | actor1)")

  if (s$model == "tie") {

    fit   <- .glmm_fit_one(s$df, s$stat_names, s$ordinal, random, engine, verbose, ...)
    coefs <- .glmm_fixef(fit, engine)

    .remstimate_wrap(
      coefficients = coefs,
      stat_names   = s$stat_names,
      loglik       = .glmm_loglik(fit),
      stacked_data = s$df,
      backend_fit  = fit,
      model        = "tie",
      method       = "GLMM",
      engine       = engine,
      ordinal      = s$ordinal
    )

  } else {

    fit_s <- fit_r <- NULL
    if (!is.null(s$df$sender))
      fit_s <- .glmm_fit_one(s$df$sender, s$stat_names$sender_model,
                              s$ordinal, random, engine, verbose, ...)
    if (!is.null(s$df$receiver))
      fit_r <- .glmm_fit_one(s$df$receiver, s$stat_names$receiver_model,
                              s$ordinal, random, engine, verbose, ...)

    .remstimate_wrap(
      coefficients = list(sender_model   = .glmm_fixef(fit_s, engine),
                          receiver_model = .glmm_fixef(fit_r, engine)),
      stat_names   = s$stat_names,
      loglik       = c(sender = .glmm_loglik(fit_s), receiver = .glmm_loglik(fit_r)),
      stacked_data = s$df,
      backend_fit  = list(sender_model = fit_s, receiver_model = fit_r),
      model        = "actor",
      method       = "GLMM",
      engine       = engine,
      ordinal      = s$ordinal
    )
  }
}

.glmm_fit_one <- function(df, stat_names, ordinal, random, engine, verbose, ...) {
  if (engine == "coxme") {
    if (!requireNamespace("coxme", quietly = TRUE))
      stop("install.packages('coxme')")
    df$.surv_time <- rep(1, nrow(df))
    formula_obj <- stats::as.formula(paste0(
      "survival::Surv(.surv_time, obs) ~ ",
      paste(stat_names, collapse = " + "),
      " + ", deparse(random[[2]]),
      " + strata(time_index)"
    ))
    coxme::coxme(formula_obj, data = df)
  } else if (engine == "lme4") {
    lme4::glmer(fml, family = familie, data = df, ...)
  } else {
    glmmTMB::glmmTMB(fml, family = familie, data = df, ...)
  }
}

.glmm_fixef <- function(fit, engine) {
  if (is.null(fit)) return(NULL)
  if (engine == "lme4")    return(lme4::fixef(fit))
  if (engine == "glmmTMB") return(glmmTMB::fixef(fit)$cond)
  coef(fit)
}

.glmm_loglik <- function(fit) {
  if (is.null(fit)) return(NULL)
  tryCatch(as.numeric(logLik(fit)), error = function(e) NULL)
}

#' @export
print.remstimate_glmm <- function(x, ...) {
  engine <- attr(x, "engine")
  cat("REM —", attr(x, "model"), "model — GLMM [", engine, "]\n\n")

  if (attr(x, "model") == "tie") {
    cat("Fixed effects:\n"); print(x$coefficients)
    cat("\nRandom effects:\n")
    tryCatch({
      vc <- if (engine == "lme4") lme4::VarCorr(x$backend_fit) else
            glmmTMB::VarCorr(x$backend_fit)
      print(vc, comp = "Variance")
    }, error = function(e) NULL)
    if (!is.null(x$loglik)) cat("\nlogLik:", round(x$loglik, 4), "\n")
    if (engine == "lme4" && isTRUE(tryCatch(lme4::isSingular(x$backend_fit), error = function(e) FALSE)))
      message("singular fit — consider simplifying the random structure")
  } else {
    for (deel in c("sender_model", "receiver_model")) {
      if (is.null(x[[deel]])) next
      cat(deel, "— fixed effects:\n"); print(x[[deel]]$coefficients); cat("\n")
    }
  }
  invisible(x)
}

#' @export
summary.remstimate_glmm <- function(object, ...) {
  if (attr(object, "model") == "tie") summary(object$backend_fit, ...)
  else lapply(list(object$sender_model$backend_fit,
                   object$receiver_model$backend_fit), summary)
}

#' @export
coef.remstimate_glmm   <- function(object, ...) object$coefficients
#' @export
logLik.remstimate_glmm <- function(object, ...) object$loglik

#' @export
#' @method diagnostics remstimate_glmm
diagnostics.remstimate_glmm <- function(object, reh = NULL, stats = NULL,
                                         top_pct = 0.05,
                                         use_ranef = TRUE, ...) {
  model <- attr(object, "model")
  if (model == "actor") {
    warning("Diagnostics for actor-oriented GLMM not yet supported.",
            call. = FALSE)
    return(invisible(NULL))
  }

  df <- object$stacked_data
  if (is.null(df)) stop("No stacked data stored in fit object.", call. = FALSE)

  fit <- object$backend_fit
  if (use_ranef && !is.null(fit)) {
    # predict() includes random effects (BLUPs) for dyad/actor-specific predictions
    lp <- tryCatch(
      stats::predict(fit, type = "link"),
      error = function(e) NULL
    )
  } else {
    lp <- NULL
  }

  # Fallback: fixed effects only
  if (is.null(lp)) {
    stat_names <- attr(object, "statistics")
    coefs      <- object$coefficients
    X  <- as.matrix(df[, stat_names, drop = FALSE])
    lp <- as.numeric(X %*% coefs[stat_names])
  }

  out <- .diagnostics_recall(lp, df, top_pct)
  out$use_ranef <- use_ranef && !is.null(fit)
  class(out) <- c("diagnostics_remstimate", "diagnostics", "remstimate")
  out
}

#' @export
#' @method plot remstimate_glmm
plot.remstimate_glmm <- function(x, reh = NULL, stats = NULL,
                                  diagnostics_object = NULL,
                                  which = 1L, ...) {
  if (is.null(diagnostics_object))
    diagnostics_object <- diagnostics(x, reh, stats)
  if (is.null(diagnostics_object)) return(invisible(x))

  if (which == 1L) {
    lbl <- if (isTRUE(diagnostics_object$use_ranef))
      "Recall: GLMM (incl. random effects)" else "Recall: GLMM (fixed effects only)"
    .plot_recall_scatter(diagnostics_object$recall, lbl, ...)
  } else if (which == 2L && !is.null(diagnostics_object$recall_by_type)) {
    rbt <- diagnostics_object$recall_by_type
    old_par <- graphics::par(mfrow = c(1, length(rbt)))
    on.exit(graphics::par(old_par))
    for (tp in names(rbt))
      .plot_recall_scatter(rbt[[tp]], paste("Type:", tp), ...)
  } else if (which == 3L) {
    lbl <- if (isTRUE(diagnostics_object$use_ranef))
      "Prob ratio: GLMM (incl. random effects)" else "Prob ratio: GLMM (fixed effects only)"
    .plot_probratio_scatter(diagnostics_object$recall, lbl, ...)
  }
  invisible(x)
}
