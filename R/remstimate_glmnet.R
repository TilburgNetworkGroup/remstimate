# GLMNET backend - ridge / lasso / elastic net via glmnet::cv.glmnet
#
# remstimate(reh, stats, method = "GLMNET", alpha = 1, lambda_select = "1se")

.remstimate_glmnet <- function(reh, stats,
                                alpha         = 1,
                                nfolds        = 10L,
                                lambda_select = c("1se", "min"),
                                unpenalized   = NULL,
                                penalized     = NULL,
                                ...) {

  lambda_select <- match.arg(lambda_select)
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("install.packages('glmnet')")

  s <- .remstimate_make_stack(reh, stats, add_actors = FALSE)

  .check_penalty_names(unpenalized, penalized,
    valid = if (s$model == "tie") s$stat_names
            else union(s$stat_names$sender_model, s$stat_names$receiver_model))

  if (s$model == "tie") {

    res <- .glmnet_fit_one(df = s$df, stat_names = s$stat_names, ordinal = s$ordinal,
                           alpha = alpha, nfolds = nfolds, lambda_select = lambda_select,
                           unpenalized = unpenalized, penalized = penalized, ...)

    .remstimate_wrap(
      coefficients = res$coefs,
      stat_names   = s$stat_names,
      stacked_data = s$df,
      backend_fit  = res$cv_fit,
      model        = "tie",
      method       = "GLMNET",
      engine       = "glmnet",
      ordinal      = s$ordinal,
      extra = list(
        lambda_min    = res$lambda_min,
        lambda_1se    = res$lambda_1se,
        lambda_sel    = res$lambda_sel,
        lambda_select = lambda_select,
        alpha         = alpha,
        unpenalized   = res$unpenalized
      )
    )

  } else {

    res_s <- res_r <- NULL
    if (!is.null(s$df$sender))
      res_s <- .glmnet_fit_one(s$df$sender, s$stat_names$sender_model,
                                s$ordinal, alpha, nfolds, lambda_select,
                                unpenalized = unpenalized, penalized = penalized, ...)
    if (!is.null(s$df$receiver))
      res_r <- .glmnet_fit_one(s$df$receiver, s$stat_names$receiver_model,
                                s$ordinal, alpha, nfolds, lambda_select,
                                unpenalized = unpenalized, penalized = penalized, ...)

    .remstimate_wrap(
      coefficients = list(sender_model   = if (!is.null(res_s)) res_s$coefs else NULL,
                          receiver_model = if (!is.null(res_r)) res_r$coefs else NULL),
      stat_names   = s$stat_names,
      stacked_data = s$df,
      backend_fit  = list(sender_model   = if (!is.null(res_s)) res_s$cv_fit else NULL,
                          receiver_model = if (!is.null(res_r)) res_r$cv_fit else NULL),
      model   = "actor",
      method  = "GLMNET",
      engine  = "glmnet",
      ordinal = s$ordinal,
      extra   = list(alpha = alpha, lambda_select = lambda_select,
                     unpenalized = list(
                       sender_model   = if (!is.null(res_s)) res_s$unpenalized else NULL,
                       receiver_model = if (!is.null(res_r)) res_r$unpenalized else NULL))
    )
  }
}

.glmnet_fit_one <- function(df, stat_names, ordinal, alpha, nfolds, lambda_select,
                            unpenalized = NULL, penalized = NULL, ...) {
  mm <- .remstimate_model_matrix(df, stat_names, ordinal)
  X  <- mm$X
  cn <- colnames(X)

  # Which columns form the intercept / baseline structure and must NOT be
  # penalised. Default (shared with SHRINKEM): the intercept-like indicator
  # columns, then the caller's additive controls - penalty = list(unpenalized =,
  # penalized =). See .resolve_unpenalized().
  unpen <- .resolve_unpenalized(X, cn, unpenalized, penalized)

  # A constant column (in practice only the all-ones 'baseline', which
  # .drop_constant_stats deliberately keeps) cannot be a glmnet predictor:
  # glmnet standardises its design matrix, so a zero-variance column collapses to
  # a 0 coefficient while its level leaks into glmnet's intrinsic intercept -
  # exactly the 'baseline shrunk to 0' symptom, which also inflates the other
  # slopes. Route a single all-ones column to that (unpenalised) intercept and
  # map it back afterwards. Non-constant intercept columns (start / end, type
  # dummies) stay in X with penalty.factor 0.
  is_const <- vapply(cn, function(nm) length(unique(X[, nm][!is.na(X[, nm])])) <= 1L,
                     logical(1))
  is_ones  <- vapply(cn, function(nm) { z <- X[, nm][!is.na(X[, nm])]
                                        length(z) > 0L && all(z == 1) }, logical(1))
  const_cols    <- cn[is_const]
  ones_cols     <- cn[is_const & is_ones]
  intercept_col <- if (length(ones_cols)) ones_cols[1L] else NULL

  Xfit <- X[, setdiff(cn, const_cols), drop = FALSE]
  if (ncol(Xfit) < 2L)
    stop("glmnet needs >= 2 varying statistics to penalise; this model has: ",
         paste(colnames(Xfit), collapse = ", "),
         ". For a model with a single varying effect, do not use 'penalty'.",
         call. = FALSE)

  pen_factor <- as.numeric(!(colnames(Xfit) %in% unpen))   # 0 = unpenalised

  cv_fit <- glmnet::cv.glmnet(Xfit, mm$y,
                               family         = if (!ordinal) "poisson" else "binomial",
                               alpha          = alpha,
                               nfolds         = nfolds,
                               offset         = mm$offset,
                               weights        = mm$weights,
                               penalty.factor = pen_factor,
                               intercept      = !is.null(intercept_col),
                               ...)

  lambda_sel <- if (lambda_select == "min") cv_fit$lambda.min else cv_fit$lambda.1se
  ruw        <- glmnet::coef.glmnet(cv_fit, s = lambda_sel)
  raw        <- as.vector(ruw); names(raw) <- rownames(ruw)

  # Rebuild the coefficient vector in the ORIGINAL statistic order so it aligns
  # 1:1 with dimnames(stats) for diagnostics.remstimate. The all-ones intercept
  # column (if any) takes glmnet's fitted intercept; every other statistic takes
  # its (possibly zero) glmnet estimate.
  coefs <- stats::setNames(numeric(length(stat_names)), stat_names)
  pen   <- raw[names(raw) != "(Intercept)"]
  coefs[names(pen)] <- pen
  if (!is.null(intercept_col)) coefs[[intercept_col]] <- raw[["(Intercept)"]]

  list(coefs       = coefs,
       cv_fit      = cv_fit,
       lambda_min  = cv_fit$lambda.min,
       lambda_1se  = cv_fit$lambda.1se,
       lambda_sel  = lambda_sel,
       unpenalized = union(colnames(Xfit)[pen_factor == 0], intercept_col))
}

#' @export
print.remstimate_glmnet <- function(x, ...) {
  methode <- if (x$alpha == 1) "lasso" else if (x$alpha == 0) "ridge" else
    paste0("elastic net (alpha=", x$alpha, ")")
  cat("REM -", attr(x, "model"), "model - GLMNET [", methode, "]\n\n")
  lbl_sel <- if ((x$lambda_select %||% "1se") == "min") "lambda.min" else "lambda.1se"
  if (attr(x, "model") == "tie") {
    cat("Coefficients (", lbl_sel, "):\n", sep = "")
    print(x$coefficients)
    cat("\nlambda.min =", round(x$lambda_min, 6),
        " | lambda.1se =", round(x$lambda_1se, 6), "\n")
    if (length(x$unpenalized))
      cat("Unpenalised (intercept structure):",
          paste(x$unpenalized, collapse = ", "), "\n")
  } else {
    for (deel in c("sender_model", "receiver_model")) {
      if (is.null(x[[deel]])) next
      cat(deel, ":\n"); print(x[[deel]]$coefficients)
      up <- x$unpenalized[[deel]]
      if (length(up))
        cat("  unpenalised:", paste(up, collapse = ", "), "\n")
      cat("\n")
    }
  }
  invisible(x)
}

#' @export
summary.remstimate_glmnet <- function(object, ...) {
  if (attr(object, "model") == "tie") print(object$backend_fit)
  else lapply(object$backend_fit, print)
}

#' @export
coef.remstimate_glmnet <- function(object, ...) object$coefficients

# A penalised fit is, for diagnostic purposes, an ordinary REM fit whose
# coefficients happen to be regularised (with the baseline left unpenalised, see
# .glmnet_fit_one). So diagnostics are obtained exactly as for MLE: drop the
# GLMNET subclass and dispatch diagnostics.remstimate, which builds the
# Schoenfeld residuals, waiting-time rates and recall from object$coefficients /
# attr(object,"where_is_baseline"). This works unchanged for tie and actor
# models. The selected penalty and mixing parameter are attached for reference.
#' @export
#' @method diagnostics remstimate_glmnet
diagnostics.remstimate_glmnet <- function(object, reh, stats, top_pct = 0.05, ...) {
  # Duration models are stacked as a tie-like design but need the start/end
  # recall structure, not the tie Schoenfeld / waiting-time diagnostics (and
  # diagnostics.remstimate() rejects a remstats_durem object). Route them to the
  # shared point-estimate durem machinery.
  if (inherits(reh, "remify_durem") || inherits(stats, "remstats_durem"))
    return(.point_durem_diagnostics(
      object, reh, stats, top_pct = top_pct,
      surprise_threshold = list(...)$surprise_threshold %||% 0.2))

  object2 <- object
  class(object2) <- "remstimate"
  diag_obj <- diagnostics(object2, reh, stats, top_pct = top_pct, ...)
  diag_obj$.alpha      <- object$alpha
  diag_obj$.lambda_sel <- object$lambda_sel
  diag_obj
}

# GLMNET has no dedicated diagnostics-plot method: diagnostics() returns a
# c("diagnostics","remstimate") object, so plot()/print() dispatch to
# plot.diagnostics / print.diagnostics exactly as for MLE. The plot method below
# only adds the glmnet-specific regularisation-path panel (which = 9); every
# other 'which' is delegated to plot.remstimate -> plot.diagnostics.
#' @export
#' @method plot remstimate_glmnet
plot.remstimate_glmnet <- function(x, reh, diagnostics = NULL,
                                    which = 1:3, effects = NULL,
                                    sender_effects = NULL, receiver_effects = NULL,
                                    ...) {
  # Duration models use the scalar-'which' durem recall plotter (1-11), a
  # different numbering from the tie/actor plot.diagnostics interface.
  dots <- list(...)
  if (inherits(reh, "remify_durem") || inherits(dots$stats, "remstats_durem")) {
    d <- diagnostics %||% diagnostics(x, reh, dots$stats)
    for (w in which) plot.diagnostics_durem(d, which = w)
    return(invisible(x))
  }

  base_which <- which[which <= 8L]
  if (length(base_which)) {
    x_base <- x; class(x_base) <- "remstimate"
    plot.remstimate(x_base, reh, diagnostics = diagnostics,
                    which = base_which, effects = effects,
                    sender_effects = sender_effects,
                    receiver_effects = receiver_effects, ...)
  }

  # (9) cross-validation / regularisation path from cv.glmnet, with the selected
  # lambda marked. One panel per sub-model (tie: one; actor: sender + receiver).
  if (9L %in% which) {
    methode_lbl <- if (x$alpha == 1) "Lasso" else if (x$alpha == 0) "Ridge" else
      sprintf("Elastic net (alpha=%g)", x$alpha)

    cv_fits <- if (inherits(x$backend_fit, "cv.glmnet")) list(tie = x$backend_fit)
    else Filter(Negate(is.null), x$backend_fit)

    op <- graphics::par(mfrow = c(1L, length(cv_fits)), no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)

    for (nm in names(cv_fits)) {
      plot(cv_fits[[nm]])
      graphics::title(main = paste(methode_lbl, if (nm == "tie") "" else paste("-", nm)),
                      line = 3)
      graphics::abline(v = log(x$lambda_sel %||% cv_fits[[nm]]$lambda.1se),
                       col = "darkorange", lwd = 2, lty = 2)
    }
  }
  invisible(x)
}
