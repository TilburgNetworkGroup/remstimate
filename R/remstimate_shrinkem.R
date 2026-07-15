# SHRINKEM backend - approximate Bayesian regularization (ridge / lasso /
# horseshoe) via the shrinkem package, applied to the MLE coefficients and
# their error covariance. The intercept / baseline structure is left
# unpenalized by default (the same 0/1-indicator rule as GLMNET, via
# .intercept_like_stats); adjust with penalty = list(unpenalized = , penalized = )
# (add to / remove from the exemption).
#
# remstimate(reh, stats, penalty = list(prior = "horseshoe"), approach = "Bayesian")

# --------------------------------------------------------------------------
# Estimability screen for the Bayesian (shrinkem) backend.
#
# shrinkem regularizes a finite MLE (coef, vcov). When the statistics are
# collinear the information matrix is singular, that MLE does not exist, and the
# trust-region optimiser used by .remstimate_mle_cpp() diverges along the flat
# ridge until its internal eigen() aborts with
#   "infinite or missing values in 'x'".
# Screen the design first and stop with a message naming the offending
# statistics. 'X' is the stacked design for one model block; constant columns
# (e.g. the all-ones baseline / intercept structure) are ignored.
# --------------------------------------------------------------------------
.assert_estimable_stats <- function(X, stat_names, block = NULL) {
  where <- if (is.null(block)) "" else sprintf(" in the %s model", block)

  dup <- unique(stat_names[duplicated(stat_names)])
  if (length(dup))
    stop(sprintf(
      "penalised REM: duplicated statistic%s%s: %s. The design is rank-deficient, so the MLE required by the Bayesian (shrinkem) backend does not exist. Rename or drop the duplicate before estimating.",
      if (length(dup) > 1L) "s" else "", where, paste(dup, collapse = ", ")),
      call. = FALSE)

  X    <- as.matrix(X)
  keep <- vapply(seq_len(ncol(X)), function(j) {
    z <- X[, j]; z <- z[is.finite(z)]
    length(z) > 0L && stats::sd(z) > 0
  }, logical(1))
  Xc <- X[, keep, drop = FALSE]
  Xc <- Xc[stats::complete.cases(Xc), , drop = FALSE]
  if (ncol(Xc) < 2L || nrow(Xc) < 2L) return(invisible(TRUE))

  qrX <- qr(scale(Xc), tol = 1e-7)
  if (qrX$rank < ncol(Xc)) {
    aliased <- colnames(Xc)[qrX$pivot[(qrX$rank + 1L):ncol(Xc)]]
    stop(sprintf(
      "penalised REM: collinear statistic%s%s: %s. The design is rank-deficient, so the MLE required by the Bayesian (shrinkem) backend does not exist. Drop or respecify %s before estimating.",
      if (length(aliased) > 1L) "s" else "", where, paste(aliased, collapse = ", "),
      if (length(aliased) > 1L) "them" else "it"),
      call. = FALSE)
  }
  invisible(TRUE)
}

.remstimate_shrinkem <- function(reh, stats,
                                 type        = "horseshoe",
                                 ncores      = 1L,
                                 seed        = NULL,
                                 unpenalized = NULL,
                                 penalized   = NULL,
                                 ...) {

  if (!requireNamespace("shrinkem", quietly = TRUE))
    stop("install.packages('shrinkem')")
  if (!is.null(seed)) set.seed(seed)

  is_durem <- inherits(stats, "remstats_durem")
  model    <- if (inherits(stats, "aomstats")) "actor" else "tie"  # durem -> "tie"
  ordinal  <- isTRUE(attr(reh, "ordinal"))

  # Stacked design: resolves the unpenalised set (intercept-structure statistics,
  # 0/1 indicators, shared with the GLMNET backend) and is screened for
  # rank-deficiency before the MLE runs.
  s <- .remstimate_make_stack(reh, stats, add_actors = FALSE)

  .check_penalty_names(unpenalized, penalized,
    valid = if (model == "tie") s$stat_names
            else union(s$stat_names$sender_model, s$stat_names$receiver_model))

  # Fail loud on a rank-deficient design *before* the MLE. shrinkem regularizes a
  # finite MLE (coef, vcov); that MLE does not exist when statistics are collinear
  # (e.g. a duplicated statistic), and .remstimate_mle_cpp()'s trust-region
  # optimiser would otherwise diverge along the flat ridge until eigen() aborts
  # with an opaque error. Duration models fit via glm(), which drops aliased
  # terms to NA without crashing, so they are exempt.
  if (!is_durem) {
    if (model == "tie") {
      .assert_estimable_stats(s$df[, s$stat_names, drop = FALSE], s$stat_names)
    } else {
      if (!is.null(s$df$sender))
        .assert_estimable_stats(s$df$sender[, s$stat_names$sender_model, drop = FALSE],
                                s$stat_names$sender_model, "sender")
      if (!is.null(s$df$receiver))
        .assert_estimable_stats(s$df$receiver[, s$stat_names$receiver_model, drop = FALSE],
                                s$stat_names$receiver_model, "receiver")
    }
  }

  # The (coef, vcov) that shrinkem regularizes come from the MLE. The C++ MLE
  # backend only accepts tie/actor stats, so for duration models we source the
  # estimates from the durem GLM instead (a single tie-style coefficient block).
  mle <- if (is_durem) {
    .remstimate_glm(stack_stats(stats, reh, add_actors = FALSE))
  } else {
    .remstimate_mle_cpp(reh, stats, ncores = ncores)
  }

  if (model == "tie") {

    unpen_names <- .resolve_unpenalized(s$df, s$stat_names, unpenalized, penalized)
    res <- .shrinkem_fit_one(mle$coefficients, mle$vcov, type = type,
                             unpenalized = unpen_names, ...)

    .remstimate_wrap(
      coefficients = res$coefs,
      stat_names   = names(mle$coefficients),
      loglik       = mle$loglik,
      backend_fit  = res$fit,
      model        = "tie",
      method       = "SHRINKEM",
      engine       = "shrinkem",
      ordinal      = ordinal,
      extra = list(
        shrinkem_type = type,
        estimates     = res$estimates,
        selected      = res$selected,
        unpenalized   = res$unpenalized
      )
    )

  } else {

    unpen_s <- .resolve_unpenalized(s$df$sender,   s$stat_names$sender_model,   unpenalized, penalized)
    unpen_r <- .resolve_unpenalized(s$df$receiver, s$stat_names$receiver_model, unpenalized, penalized)

    res_s <- .shrinkem_fit_one(mle$sender_model$coefficients,
                               mle$sender_model$vcov, type = type,
                               unpenalized = unpen_s, ...)
    res_r <- .shrinkem_fit_one(mle$receiver_model$coefficients,
                               mle$receiver_model$vcov, type = type,
                               unpenalized = unpen_r, ...)

    .remstimate_wrap(
      coefficients = list(sender_model   = res_s$coefs,
                          receiver_model = res_r$coefs),
      stat_names   = attr(mle, "statistics"),
      backend_fit  = list(sender_model   = res_s$fit,
                          receiver_model = res_r$fit),
      model   = "actor",
      method  = "SHRINKEM",
      engine  = "shrinkem",
      ordinal = ordinal,
      extra = list(
        shrinkem_type = type,
        estimates     = list(sender_model   = res_s$estimates,
                             receiver_model = res_r$estimates),
        selected      = list(sender_model   = res_s$selected,
                             receiver_model = res_r$selected),
        unpenalized   = list(sender_model   = res_s$unpenalized,
                             receiver_model = res_r$unpenalized)
      )
    )
  }
}

# regularize one (coef, vcov) block. 'unpenalized' is a character vector of
# statistic names left unpenalised (the intercept structure); when NULL it falls
# back to an exact 'baseline'. shrinkem itself takes a logical vector.
.shrinkem_fit_one <- function(coef, Sigma, type, unpenalized = NULL, ...) {
  stat_names  <- names(coef)
  unpen_names <- if (!is.null(unpenalized)) unpenalized
                 else stat_names[.remstimate_find_baseline(stat_names)]
  unpen <- stat_names %in% unpen_names

  fit <- shrinkem::shrinkem(x = coef, Sigma = Sigma, type = type,
                            unpenalized = unpen, ...)

  coefs <- fit$estimates$shrunk.mode
  names(coefs) <- stat_names
  list(coefs       = coefs,
       fit         = fit,
       estimates   = fit$estimates,
       selected    = fit$estimates$nonzero,
       unpenalized = stat_names[unpen])
}

# round numeric columns of a shrinkem estimates table for printing
.shrinkem_round_est <- function(est, digits = 3) {
  num <- vapply(est, is.numeric, logical(1))
  est[num] <- lapply(est[num], round, digits)
  est
}

#' @export
coef.remstimate_shrinkem <- function(object, ...) object$coefficients

#' @export
print.remstimate_shrinkem <- function(x, ...) {
  cat("REM -", attr(x, "model"), "model - SHRINKEM [", x$shrinkem_type, "]\n\n")
  if (attr(x, "model") == "tie") {
    cat("Shrunken coefficients (posterior mode):\n")
    print(data.frame(estimate = round(x$coefficients, 3), selected = x$selected))
    if (length(x$unpenalized))
      cat("Unpenalised (intercept structure):",
          paste(x$unpenalized, collapse = ", "), "\n")
  } else {
    for (deel in c("sender_model", "receiver_model")) {
      if (is.null(x[[deel]])) next
      cat(deel, ":\n")
      print(data.frame(estimate = round(x[[deel]]$coefficients, 3),
                       selected = x$selected[[deel]]))
      up <- x$unpenalized[[deel]]
      if (length(up)) cat("  unpenalised:", paste(up, collapse = ", "), "\n")
      cat("\n")
    }
  }
  invisible(x)
}

#' @export
summary.remstimate_shrinkem <- function(object, ...) {
  cat("REM -", attr(object, "model"), "model - SHRINKEM [",
      object$shrinkem_type, "] regularization\n\n")
  if (attr(object, "model") == "tie") {
    print(.shrinkem_round_est(object$estimates))
  } else {
    for (deel in c("sender_model", "receiver_model")) {
      est <- object$estimates[[deel]]
      if (is.null(est)) next
      cat(deel, ":\n"); print(.shrinkem_round_est(est)); cat("\n")
    }
  }
  invisible(object)
}

# diagnostics: delegate to the base MLE recall (same path as glmm); the object
# is MLE-shaped with the shrunken coefficients, so no extra work is needed
#' @export
#' @method diagnostics remstimate_shrinkem
diagnostics.remstimate_shrinkem <- function(object, reh, stats, top_pct = 0.05, ...) {
  # Duration models need the start/end recall structure (see the GLMNET note).
  if (inherits(reh, "remify_durem") || inherits(stats, "remstats_durem")) {
    out <- .point_durem_diagnostics(
      object, reh, stats, top_pct = top_pct,
      surprise_threshold = list(...)$surprise_threshold %||% 0.2)
    out$.engine <- "shrinkem"
    return(out)
  }

  class(object) <- "remstimate"
  diag_obj <- diagnostics(object, reh, stats, top_pct = top_pct, ...)
  diag_obj$.engine <- "shrinkem"
  diag_obj
}

# plot: which <= 5 delegate to base; which == 6 shows the shrinkage map
#' @export
#' @method plot remstimate_shrinkem
plot.remstimate_shrinkem <- function(x, reh, diagnostics = NULL,
                                     which = 1:3, effects = NULL,
                                     sender_effects = NULL, receiver_effects = NULL,
                                     ...) {
  # Duration models use the scalar-'which' durem recall plotter (1-11).
  dots <- list(...)
  if (inherits(reh, "remify_durem") || inherits(dots$stats, "remstats_durem")) {
    d <- diagnostics %||% diagnostics(x, reh, dots$stats)
    for (w in which) plot.diagnostics_durem(d, which = w)
    return(invisible(x))
  }

  basis_which <- which[which <= 5L]
  if (length(basis_which)) {
    x_basis <- x; class(x_basis) <- "remstimate"
    plot.remstimate(x_basis, reh, diagnostics = diagnostics,
                    which = basis_which, effects = effects,
                    sender_effects = sender_effects,
                    receiver_effects = receiver_effects, ...)
  }

  if (6L %in% which) {
    parts <- if (attr(x, "model") == "tie") list(tie = x$estimates)
             else Filter(Negate(is.null), x$estimates)

    op <- graphics::par(mfrow = c(1L, length(parts)), no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)

    for (nm in names(parts)) {
      est <- parts[[nm]]
      graphics::plot(est$input.est, est$shrunk.mode,
                     xlab = "MLE estimate", ylab = "Shrunken estimate",
                     main = paste("Shrinkage", if (nm == "tie") "" else paste("-", nm)),
                     pch = 16,
                     col = ifelse(est$nonzero, "black",
                                  grDevices::adjustcolor("black", 0.35)))
      graphics::abline(0, 1, col = "darkorange", lwd = 2, lty = 2)
      graphics::abline(h = 0, col = "grey80")
    }
  }
  invisible(x)
}
