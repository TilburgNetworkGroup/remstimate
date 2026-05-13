# MIXREM backend — eindige mengsels via flexmix::flexmix
#
# remstimate(reh, stats, method = "MIXREM",
#            random = ~ (1 + inertia | dyad), k = 2)
#
# random specificeert welke effecten per component variëren en de groeperingsvariabele
# k als vector geeft meerdere fits; gebruik bic_table() om te vergelijken

.remstimate_mixrem <- function(reh, stats,
                                random      = NULL,
                                k           = 2L,
                                concomitant = NULL,
                                nrep        = 3L,
                                ...) {

  if (!requireNamespace("flexmix", quietly = TRUE))
    stop("install.packages('flexmix')")
  if (is.null(random))
    stop("specify clustering, e.g. random = ~ (1 + inertia | dyad)")

  s       <- .remstimate_make_stack(reh, stats, add_actors = TRUE)
  geparsed <- .parse_mixrem_formula(random)

  if (length(k) == 1L) {
    if (s$model == "tie")
      return(.mixrem_fit_one(s$df, s$stat_names, s$ordinal,
                              geparsed, k, concomitant, nrep, model = "tie", ...))

    fit_s <- fit_r <- NULL
    if (!is.null(s$df$sender))
      fit_s <- .mixrem_fit_one(s$df$sender, s$stat_names$sender_model, s$ordinal,
                                geparsed, k, concomitant, nrep, model = "actor_sender", ...)
    if (!is.null(s$df$receiver))
      fit_r <- .mixrem_fit_one(s$df$receiver, s$stat_names$receiver_model, s$ordinal,
                                geparsed, k, concomitant, nrep, model = "actor_receiver", ...)
    return(list(sender_model = fit_s, receiver_model = fit_r))
  }

  # meerdere k-waarden
  fits <- lapply(k, function(ki) {
    if (s$model == "tie")
      .mixrem_fit_one(s$df, s$stat_names, s$ordinal,
                      geparsed, ki, concomitant, nrep, model = "tie", ...)
    else
      list(
        sender_model   = if (!is.null(s$df$sender))
          .mixrem_fit_one(s$df$sender, s$stat_names$sender_model, s$ordinal,
                          geparsed, ki, concomitant, nrep, model = "actor_sender", ...) else NULL,
        receiver_model = if (!is.null(s$df$receiver))
          .mixrem_fit_one(s$df$receiver, s$stat_names$receiver_model, s$ordinal,
                          geparsed, ki, concomitant, nrep, model = "actor_receiver", ...) else NULL
      )
  })
  names(fits) <- paste0("k", k)
  structure(fits, class = c("remstimate_mixrem_list", "list"), k = k)
}

.mixrem_fit_one <- function(df, stat_names, ordinal, geparsed, k,
                              concomitant, nrep, model = "tie", ...) {

  vaste_kant <- .remstimate_fixed_rhs(stat_names, ordinal)
  if (ordinal) {
    df$time   <- factor(df$time)
    vaste_kant <- sub("^-1 \\+ ", "-1 + time + ", vaste_kant)
    familie    <- flexmix::FLXMRglm(family = "binomial")
  } else {
    familie <- flexmix::FLXMRglm(family = "poisson")
  }

  fml  <- as.formula(paste("obs ~", vaste_kant, "|", geparsed$group))
  conc <- if (!is.null(concomitant)) flexmix::FLXPmultinom(concomitant) else flexmix::FLXPconstant()

  fit <- flexmix::flexmix(fml, data = df, k = k,
                           model       = familie,
                           concomitant = conc,
                           control     = list(nrep = nrep), ...)

  coef_mat <- flexmix::parameters(fit)
  if (is.null(dim(coef_mat)))
    coef_mat <- matrix(coef_mat, ncol = 1, dimnames = list(names(coef_mat), "Comp.1"))

  kansen <- flexmix::prior(fit)
  volgorde <- order(kansen, decreasing = TRUE)
  coef_mat <- coef_mat[, volgorde, drop = FALSE]
  kansen   <- kansen[volgorde]
  colnames(coef_mat) <- paste0("Component.", seq_len(k))

  .remstimate_wrap(
    coefficients = coef_mat,
    stat_names   = stat_names,
    loglik       = flexmix::logLik(fit),
    formula      = fml,
    stacked_data = df,
    backend_fit  = fit,
    model        = sub("_.*", "", model),
    method       = "MIXREM",
    engine       = "flexmix",
    ordinal      = ordinal,
    extra        = list(k = k, prior_probs = kansen, bic = flexmix::BIC(fit))
  )
}

.parse_mixrem_formula <- function(random) {
  rhs   <- deparse(random[[2]])
  groep_match <- regmatches(rhs,
    regexpr("\\|\\s*([A-Za-z_.][A-Za-z_.0-9]*)", rhs, perl = TRUE))
  if (!length(groep_match))
    stop("random must contain a grouping variable after |, e.g. ~ (1 + inertia | dyad)")

  groep  <- trimws(sub("^\\|\\s*", "", groep_match))
  binnen <- sub("\\|.*", "", sub("^\\(", "", rhs))
  termen <- trimws(strsplit(binnen, "\\+")[[1]])
  termen <- termen[termen != "1" & nchar(termen) > 0]

  list(effects = if (length(termen)) termen else NULL, group = groep)
}

#' @export
bic_table <- function(x, ...) UseMethod("bic_table")

#' @export
bic_table.remstimate_mixrem_list <- function(x, ...) {
  k_vals <- attr(x, "k")
  bics   <- sapply(x, function(fit)
    if (inherits(fit, "remstimate_mixrem")) fit$bic else NA_real_)
  out <- data.frame(k = k_vals, BIC = bics)
  out$delta_BIC <- out$BIC - min(out$BIC, na.rm = TRUE)
  out[order(out$BIC), ]
}

#' @export
print.remstimate_mixrem <- function(x, ...) {
  cat("REM —", attr(x, "model"), "model — MIXREM [flexmix, k =", x$k, "]\n\n")
  cat("Mixing proportions:\n"); print(round(x$prior_probs, 4))
  cat("\nCoefficients per component:\n"); print(round(x$coefficients, 4))
  cat("\nlogLik:", round(as.numeric(x$loglik), 4), " BIC:", round(x$bic, 2), "\n")
  invisible(x)
}

#' @export
print.remstimate_mixrem_list <- function(x, ...) {
  cat("MIXREM fits voor k =", paste(attr(x, "k"), collapse = ", "), "\n")
  print(bic_table(x))
  invisible(x)
}

#' @export
summary.remstimate_mixrem  <- function(object, ...) summary(object$backend_fit, ...)
#' @export
coef.remstimate_mixrem     <- function(object, ...) object$coefficients
#' @export
logLik.remstimate_mixrem   <- function(object, ...) object$loglik

#' @export
#' @method diagnostics remstimate_mixrem
diagnostics.remstimate_mixrem <- function(object, reh = NULL, stats = NULL,
                                           top_pct = 0.05, ...) {
  model <- attr(object, "model")
  if (model == "actor") {
    warning("Diagnostics for actor-oriented MIXREM not yet supported.",
            call. = FALSE)
    return(invisible(NULL))
  }

  df <- object$stacked_data
  if (is.null(df)) stop("No stacked data stored in fit object.", call. = FALSE)

  stat_names <- attr(object, "statistics")
  coef_mat   <- object$coefficients  # [P x K]
  fit        <- object$backend_fit
  K          <- ncol(coef_mat)

  # Posterior-weighted linear predictor
  X <- as.matrix(df[, stat_names, drop = FALSE])

  # Subset coef_mat to stat_names rows (in case of extra rows like time FE)
  coef_rows <- intersect(rownames(coef_mat), stat_names)
  if (length(coef_rows) == 0L) {
    # Fallback: assume row order matches stat_names
    coef_rows <- stat_names[seq_len(min(nrow(coef_mat), length(stat_names)))]
  }

  if (!is.null(fit) && K > 1L) {
    post <- tryCatch(flexmix::posterior(fit), error = function(e) NULL)
    if (!is.null(post)) {
      # LP per component, then posterior-weighted average
      lp_k <- lapply(seq_len(K), function(k) {
        bk <- coef_mat[coef_rows, k]
        as.numeric(X[, coef_rows, drop = FALSE] %*% bk)
      })
      lp <- rowSums(vapply(seq_len(K), function(k) post[, k] * lp_k[[k]],
                           numeric(nrow(X))))
    } else {
      # Fallback: use first (largest) component
      bk <- coef_mat[coef_rows, 1L]
      lp <- as.numeric(X[, coef_rows, drop = FALSE] %*% bk)
    }
  } else {
    # Single component
    bk <- if (is.matrix(coef_mat)) coef_mat[coef_rows, 1L] else coef_mat[coef_rows]
    lp <- as.numeric(X[, coef_rows, drop = FALSE] %*% bk)
  }

  out <- .diagnostics_recall(lp, df, top_pct)

  # Per-component recall
  if (K > 1L) {
    out$recall_by_component <- list()
    for (k in seq_len(K)) {
      bk <- coef_mat[coef_rows, k]
      lp_k <- as.numeric(X[, coef_rows, drop = FALSE] %*% bk)
      out$recall_by_component[[paste0("Component.", k)]] <-
        .recall_block(lp_k, which(df$obs == 1L), df$time, top_pct)
    }
  }

  class(out) <- c("diagnostics_mixrem", "diagnostics_remstimate",
                   "diagnostics", "remstimate")
  out
}

#' @export
#' @method print diagnostics_mixrem
print.diagnostics_mixrem <- function(x, ...) {
  cat("Diagnostics for Relational Event Model (MIXREM)\n")
  cat(strrep("-", 50), "\n")

  if (!is.null(x$recall)) {
    cat("\nRecall (posterior-weighted):\n")
    .print_recall_summary(x$recall, "Joint")
  }
  if (!is.null(x$recall_by_type)) {
    cat("\nRecall by type:\n")
    for (tp in names(x$recall_by_type))
      .print_recall_summary(x$recall_by_type[[tp]], paste0("  ", tp))
  }
  if (!is.null(x$recall_by_component)) {
    cat("\nRecall by component:\n")
    for (k in names(x$recall_by_component))
      .print_recall_summary(x$recall_by_component[[k]], paste0("  ", k))
  }
  invisible(x)
}

#' @export
#' @method plot remstimate_mixrem
plot.remstimate_mixrem <- function(x, reh = NULL, stats = NULL,
                                    diagnostics_object = NULL,
                                    which = 1L, ...) {
  if (is.null(diagnostics_object))
    diagnostics_object <- diagnostics(x, reh, stats)
  if (is.null(diagnostics_object)) return(invisible(x))

  if (which == 1L) {
    .plot_recall_scatter(diagnostics_object$recall,
                         "Recall: MIXREM (posterior-weighted)", ...)
  } else if (which == 2L && !is.null(diagnostics_object$recall_by_component)) {
    rbt <- diagnostics_object$recall_by_component
    old_par <- graphics::par(mfrow = c(1, length(rbt)))
    on.exit(graphics::par(old_par))
    for (k in names(rbt))
      .plot_recall_scatter(rbt[[k]], k, ...)
  } else if (which == 3L && !is.null(diagnostics_object$recall_by_type)) {
    rbt <- diagnostics_object$recall_by_type
    old_par <- graphics::par(mfrow = c(1, length(rbt)))
    on.exit(graphics::par(old_par))
    for (tp in names(rbt))
      .plot_recall_scatter(rbt[[tp]], paste("Type:", tp), ...)
  } else if (which == 4L) {
    .plot_probratio_scatter(diagnostics_object$recall,
                             "Prob ratio: MIXREM (posterior-weighted)", ...)
  }
  invisible(x)
}
