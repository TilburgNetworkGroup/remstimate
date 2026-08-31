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
                                 group       = NULL,
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
    # Column scales from the stacked design, on each statistic's structural
    # support. For duration models the MLE is sourced from a separate
    # stack_stats() call, but both stacks are built from the same reh / stats
    # with add_actors = FALSE, so the row set - and hence each column's SD - is
    # the same.
    scl <- .penalty_scales(s$df, s$stat_names, process = s$df$process)
    res <- .shrinkem_fit_one(mle$coefficients, mle$vcov, type = type,
                             unpenalized = unpen_names,
                             scale = scl, group = group, ...)

    uit <- .remstimate_wrap(
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
        cred_level    = res$cred_level,
        unpenalized   = res$unpenalized
      )
    )
    # plot.diagnostics draws its posterior-density (which = 4) and trace
    # (which = 5) panels from object$draws and attr(object, "nchains"). Exposing
    # the shrinkem chain in the same place makes those panels work for this
    # backend as well, instead of being HMC-only.
    uit$draws <- .align_draws(res$fit$draws$beta, names(uit$coefficients))
    attr(uit, "nchains") <- 1L
    uit

  } else {

    unpen_s <- .resolve_unpenalized(s$df$sender,   s$stat_names$sender_model,   unpenalized, penalized)
    unpen_r <- .resolve_unpenalized(s$df$receiver, s$stat_names$receiver_model, unpenalized, penalized)

    scl_s <- .penalty_scales(s$df$sender,   s$stat_names$sender_model)
    scl_r <- .penalty_scales(s$df$receiver, s$stat_names$receiver_model)

    res_s <- .shrinkem_fit_one(mle$sender_model$coefficients,
                               mle$sender_model$vcov, type = type,
                               unpenalized = unpen_s,
                               scale = scl_s, group = group, ...)
    res_r <- .shrinkem_fit_one(mle$receiver_model$coefficients,
                               mle$receiver_model$vcov, type = type,
                               unpenalized = unpen_r,
                               scale = scl_r, group = group, ...)

    uit <- .remstimate_wrap(
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
        cred_level    = res_s$cred_level,
        unpenalized   = list(sender_model   = res_s$unpenalized,
                             receiver_model = res_r$unpenalized)
      )
    )
    # per-component chains: plot.diagnostics reads object[[model]]$draws
    uit$sender_model$draws <-
      .align_draws(res_s$fit$draws$beta, names(uit$sender_model$coefficients))
    uit$receiver_model$draws <-
      .align_draws(res_r$fit$draws$beta, names(uit$receiver_model$coefficients))
    attr(uit, "nchains") <- 1L
    uit
  }
}

# regularize one (coef, vcov) block. 'unpenalized' is a character vector of
# statistic names left unpenalised (the intercept structure); when NULL it falls
# back to an exact 'baseline'. shrinkem itself takes a logical vector.
.shrinkem_fit_one <- function(coef, Sigma, type, unpenalized = NULL,
                              scale = NULL, group = NULL, ...) {
  stat_names  <- names(coef)
  unpen_names <- if (!is.null(unpenalized)) unpenalized
                 else stat_names[.remstimate_find_baseline(stat_names)]
  unpen <- stat_names %in% unpen_names

  # Put the estimates on a comparable scale before the prior sees them (see
  # .penalty_scales). With D = diag(scale), beta* = D beta and Sigma* = D Sigma D
  # is an exact linear reparameterisation - maximum likelihood is equivariant and
  # the delta method is exact for a linear map - so this is the same Gaussian
  # approximation a fit on a standardised design would have produced, and only
  # the meaning of the prior changes. The shrunken quantities are mapped back
  # with D^-1 below; selection is invariant because D is positive diagonal.
  s <- if (is.null(scale)) rep(1, length(coef)) else as.numeric(scale[stat_names])
  s[!is.finite(s) | s <= 0] <- 1

  # One shrinkage parameter per group (shrinkem estimates lambda2 per group).
  # 'group' is keyed by statistic name; statistics it does not mention share a
  # single default group.
  grp <- 1
  if (!is.null(group)) {
    grp <- if (!is.null(names(group))) as.character(group[stat_names])
           else as.character(group)
    grp[is.na(grp)] <- ".default"
  }

  # the credible level shrinkem used, so that print() and plot() can label the
  # interval with the level that actually produced it rather than assuming 0.95
  cred <- list(...)$cred.level %||% 0.95

  fit <- shrinkem::shrinkem(x = coef * s, Sigma = Sigma * outer(s, s),
                            type = type, unpenalized = unpen, group = grp, ...)

  # Back to the original scale. 'nonzero' is scale-invariant and is left alone.
  back <- c("input.est", "shrunk.mean", "shrunk.median", "shrunk.mode",
            "shrunk.lower", "shrunk.upper")
  for (kolom in intersect(back, names(fit$estimates)))
    fit$estimates[[kolom]] <- fit$estimates[[kolom]] / s
  if (!is.null(fit$draws$beta)) fit$draws$beta <- sweep(fit$draws$beta, 2L, s, "/")
  if (!is.null(fit$input.est))  fit$input.est  <- fit$input.est / s

  coefs <- fit$estimates$shrunk.mode
  names(coefs) <- stat_names
  list(coefs       = coefs,
       fit         = fit,
       estimates   = fit$estimates,
       selected    = fit$estimates$nonzero,
       cred_level  = cred,
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

# plot.diagnostics labels its posterior-density and trace panels from
# names(object$coefficients) but indexes the chain by POSITION
# (object$draws[, which_effects[p]]). shrinkem names the columns of draws$beta
# after the estimate it was given, so the two orders agree by construction -
# reorder anyway, so that a future change to the coefficient vector cannot
# silently plot one effect under another's label.
.align_draws <- function(draws, coef_names) {
  if (is.null(draws) || is.null(colnames(draws))) return(draws)
  if (!all(coef_names %in% colnames(draws))) return(draws)
  draws[, coef_names, drop = FALSE]
}

# Posterior mode with its credible interval. The interval is what the 'nonzero'
# flag is derived from (0 inside or outside it), so printing the interval itself
# is strictly more informative and makes the level explicit; 'nonzero' remains
# available through summary().
.shrinkem_coef_table <- function(coefs, est, cred, digits = 3) {
  out <- data.frame(estimate = round(coefs, digits))
  if (!is.null(est) && all(c("shrunk.lower", "shrunk.upper") %in% names(est))) {
    lo <- 100 * (1 - cred) / 2
    out[[paste0(format(lo, trim = TRUE), "%")]]       <- round(est$shrunk.lower, digits)
    out[[paste0(format(100 - lo, trim = TRUE), "%")]] <- round(est$shrunk.upper, digits)
  }
  out
}

#' @export
print.remstimate_shrinkem <- function(x, digits = 3, ...) {
  cred <- x$cred_level %||% 0.95
  kop  <- sprintf("Shrunken coefficients (posterior mode, %s%% credible interval):",
                  format(100 * cred, trim = TRUE))
  cat("REM -", attr(x, "model"), "model - SHRINKEM [", x$shrinkem_type, "]\n\n")
  if (attr(x, "model") == "tie") {
    cat(kop, "\n")
    print(.shrinkem_coef_table(x$coefficients, x$estimates, cred, digits))
    if (length(x$unpenalized))
      cat("Unpenalised (intercept structure):",
          paste(x$unpenalized, collapse = ", "), "\n")
  } else {
    for (deel in c("sender_model", "receiver_model")) {
      if (is.null(x[[deel]])) next
      cat(deel, ":\n"); cat(kop, "\n")
      print(.shrinkem_coef_table(x[[deel]]$coefficients, x$estimates[[deel]],
                                 cred, digits))
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
plot.remstimate_shrinkem <- function(x, reh = NULL, diagnostics = NULL,
                                     which, effects = NULL,
                                     sender_effects = NULL, receiver_effects = NULL,
                                     ...) {
  # Without 'reh' only the parameter panels are available, so plot(fit) shows the
  # posterior densities and traces of the shrinkem chain; with 'reh' the default
  # is the usual model-fit set.
  if (missing(which)) which <- if (is.null(reh)) 4:6 else c(1:3, 6L)
  # 'which' follows plot.diagnostics throughout, except for 6, which is this
  # backend's own MLE-vs-shrunk panel. Posterior densities (4) and traces (5)
  # are handled by plot.diagnostics from object$draws, which this backend now
  # populates. 6 needs only the fitted object, so it is drawn before 'reh'.
  if (6L %in% which) {
    .plot_shrinkage_panel(x)
    which <- setdiff(which, 6L)
    if (!length(which)) return(invisible(x))
  }
  # 4 and 5 read object$draws only; plot.remstimate short-circuits those without
  # building a diagnostics object, so 'reh' is required only for the rest.
  if (is.null(reh) && !all(which %in% c(4L, 5L)))
    stop("'reh' is required for the model-fit plots; which = 4:5 (posterior ",
         "density and traces) can be drawn without it.", call. = FALSE)

  # Duration models use the scalar-'which' durem recall plotter (1-11).
  dots <- list(...)
  if (inherits(reh, "remify_durem") || inherits(dots$stats, "remstats_durem")) {
    d <- diagnostics %||% diagnostics(x, reh, dots$stats)
    for (w in which) plot.diagnostics_durem(d, which = w)
    return(invisible(x))
  }

  basis_which <- setdiff(which, 6L)
  if (length(basis_which)) {
    x_basis <- x; class(x_basis) <- "remstimate"
    plot.remstimate(x_basis, reh, diagnostics = diagnostics,
                    which = basis_which, effects = effects,
                    sender_effects = sender_effects,
                    receiver_effects = receiver_effects, ...)
  }

  invisible(x)
}

# MLE vs shrunk estimate, one panel per model component.
.plot_shrinkage_panel <- function(x) {
  parts <- if (attr(x, "model") == "tie") list(tie = x$estimates)
           else Filter(Negate(is.null), x$estimates)
  unpen <- if (attr(x, "model") == "tie") list(tie = x$unpenalized) else x$unpenalized

  op <- graphics::par(mfrow = c(1L, length(parts)), no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  for (nm in names(parts)) {
    est  <- parts[[nm]]
    naam <- rownames(est)

    # Unpenalised effects (the baseline above all) are never shrunk and can sit
    # far from zero, which flattens the scale for the effects that were. Drop
    # them: the panel is about what the prior did.
    vrij <- !(naam %in% (unpen[[nm]] %||% character(0)))
    if (any(vrij) && !all(vrij)) {
      est <- est[vrij, , drop = FALSE]; naam <- naam[vrij]
    }
    weggelaten <- sum(!vrij)

    # One padded range on both axes, so the y = x reference really is the
    # 45-degree line and the points are not pinned to the panel edges.
    rng <- range(c(est$input.est, est$shrunk.mode, 0), finite = TRUE)
    pad <- 0.15 * diff(rng)
    if (!is.finite(pad) || pad <= 0) pad <- 0.5
    lim <- rng + c(-pad, pad)

    graphics::plot(est$input.est, est$shrunk.mode,
                   xlim = lim, ylim = lim,
                   xlab = "MLE estimate", ylab = "Shrunken estimate",
                   main = paste("Shrinkage", if (nm == "tie") "" else paste("-", nm)),
                   sub = if (weggelaten)
                     sprintf("%d unpenalised effect%s omitted", weggelaten,
                             if (weggelaten > 1L) "s" else "") else "",
                   pch = 16,
                   col = ifelse(est$nonzero, "black",
                                grDevices::adjustcolor("black", 0.35)))
    graphics::abline(0, 1, col = "darkorange", lwd = 2, lty = 2)
    graphics::abline(h = 0, col = "grey80")

    # label each point, pushing labels inward near the right-hand edge
    if (length(naam)) {
      # Stagger labels above / below in x order. Deterministic rather than
      # randomised, so the same fit always yields the same figure.
      r <- rank(est$input.est, ties.method = "first")
      kant <- ifelse(r %% 2L == 1L, 3L, 1L)
      graphics::text(est$input.est, est$shrunk.mode, labels = naam,
                     pos = kant, cex = 0.65, xpd = NA,
                     col = ifelse(est$nonzero, "black", "grey45"))
    }
  }
  invisible(NULL)
}
