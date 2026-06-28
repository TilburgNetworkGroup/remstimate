# SHRINKEM backend — approximate Bayesian regularization (ridge / lasso /
# horseshoe) via the shrinkem package, applied to the MLE coefficients and
# their error covariance. The baseline rate (if present) is left unpenalized.
#
# remstimate(reh, stats, penalty = list(prior = "horseshoe"), approach = "Bayesian")

.remstimate_shrinkem <- function(reh, stats,
                                 type   = "horseshoe",
                                 ncores = 1L,
                                 seed   = NULL,
                                 ...) {

  if (!requireNamespace("shrinkem", quietly = TRUE))
    stop("install.packages('shrinkem')")
  if (!is.null(seed)) set.seed(seed)

  model   <- if (inherits(stats, "aomstats")) "actor" else "tie"
  ordinal <- isTRUE(attr(reh, "ordinal"))

  # the MLE supplies the (coef, vcov) that shrinkem regularizes
  mle <- .remstimate_mle_cpp(reh, stats, ncores = ncores)

  if (model == "tie") {

    res <- .shrinkem_fit_one(mle$coefficients, mle$vcov, type = type, ...)

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
        selected      = res$selected
      )
    )

  } else {

    res_s <- .shrinkem_fit_one(mle$sender_model$coefficients,
                               mle$sender_model$vcov, type = type, ...)
    res_r <- .shrinkem_fit_one(mle$receiver_model$coefficients,
                               mle$receiver_model$vcov, type = type, ...)

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
                             receiver_model = res_r$selected)
      )
    )
  }
}

# regularize one (coef, vcov) block; baseline left unpenalized if present
.shrinkem_fit_one <- function(coef, Sigma, type, ...) {
  stat_names <- names(coef)
  unpen <- rep(FALSE, length(coef))
  bl    <- .remstimate_find_baseline(stat_names)
  if (!is.null(bl)) unpen[bl] <- TRUE

  fit <- shrinkem::shrinkem(x = coef, Sigma = Sigma, type = type,
                            unpenalized = unpen, ...)

  coefs <- fit$estimates$shrunk.mode
  names(coefs) <- stat_names
  list(coefs     = coefs,
       fit       = fit,
       estimates = fit$estimates,
       selected  = fit$estimates$nonzero)
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
  cat("REM —", attr(x, "model"), "model — SHRINKEM [", x$shrinkem_type, "]\n\n")
  if (attr(x, "model") == "tie") {
    cat("Shrunken coefficients (posterior mode):\n")
    print(data.frame(estimate = round(x$coefficients, 3), selected = x$selected))
  } else {
    for (deel in c("sender_model", "receiver_model")) {
      if (is.null(x[[deel]])) next
      cat(deel, ":\n")
      print(data.frame(estimate = round(x[[deel]]$coefficients, 3),
                       selected = x$selected[[deel]]))
      cat("\n")
    }
  }
  invisible(x)
}

#' @export
summary.remstimate_shrinkem <- function(object, ...) {
  cat("REM —", attr(object, "model"), "model — SHRINKEM [",
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
