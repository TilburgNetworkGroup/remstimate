# MIXREM backend - eindige mengsels via flexmix::flexmix
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

  vaste_klant <- .remstimate_fixed_rhs(stat_names, ordinal)
  if (ordinal) {
    df$time_index <- factor(df$time_index)
    df$obs_fail <- 1 - df$obs                    # flexmix binomial needs cbind(succ, fail)
    vaste_klant   <- sub("^-1 \\+ ", "-1 + time_index + ", vaste_klant)
    familie     <- flexmix::FLXMRglm(family = "binomial")
    respons     <- "cbind(obs, obs_fail)"
  } else {
    familie <- flexmix::FLXMRglm(family = "poisson")
    respons <- "obs"
  }

  fml  <- as.formula(paste(respons, "~", vaste_klant, "|", geparsed$group))
  conc <- if (!is.null(concomitant)) flexmix::FLXPmultinom(concomitant) else flexmix::FLXPconstant()

  fit <- flexmix::flexmix(fml, data = df, k = k,
                           model       = familie,
                           concomitant = conc,
                           control     = list(nrep = nrep), ...)

  coef_mat <- flexmix::parameters(fit)
  if (is.null(dim(coef_mat)))
    coef_mat <- matrix(coef_mat, ncol = 1, dimnames = list(names(coef_mat), "Comp.1"))
  rownames(coef_mat) <- gsub("^coef\\.", "", rownames(coef_mat))

  kansen <- flexmix::prior(fit)
  volgorde <- order(kansen, decreasing = TRUE)
  coef_mat <- coef_mat[, volgorde, drop = FALSE]
  kansen   <- kansen[volgorde]
  colnames(coef_mat) <- paste0("Component.", seq_len(k))

  # drop ordinal time-stratum dummies - report substantive statistics only
  keep     <- intersect(stat_names, rownames(coef_mat))
  coef_mat <- coef_mat[keep, , drop = FALSE]

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
    extra        = list(k = k, prior_probs = kansen, bic = stats::BIC(fit))
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
  cat("REM -", attr(x, "model"), "model - MIXREM [flexmix, k =", x$k, "]\n\n")
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

# MIXREM diagnostics. A finite mixture is the discrete analogue of the GLMM
# random effects: instead of continuous dyad/actor BLUPs, dyads/actors are
# clustered into components that share the same parameters. So diagnostics are
# built the same way as the GLMM random-effect-aware recall - a within-event
# ranking under a component-aware linear predictor - and routed through
# plot.diagnostics:
#   $recall              posterior-weighted recall (the cluster-aware analogue of
#                        GLMM's $recall_ranef)          -> plot.diagnostics which 3 / 6
#   $recall_by_type      per event type (typed events)  -> which 8
#   $recall_by_component per mixture component           -> which 9
# There is no single global linear predictor for a mixture, so the Schoenfeld
# residuals (which 2) and waiting-time Q-Q (which 1) are not produced; those
# panels are simply absent and plot.diagnostics skips them.
#' @export
#' @method diagnostics remstimate_mixrem
diagnostics.remstimate_mixrem <- function(object, reh, stats = NULL,
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
  rownames(coef_mat) <- gsub("^coef\\.", "", rownames(coef_mat))
  fit        <- object$backend_fit
  K          <- ncol(coef_mat)

  X <- as.matrix(df[, stat_names, drop = FALSE])

  # Subset coef_mat to stat_names rows (in case of extra rows like time FE)
  coef_rows <- intersect(rownames(coef_mat), stat_names)
  if (length(coef_rows) == 0L) {
    # Fallback: assume row order matches stat_names
    coef_rows <- stat_names[seq_len(min(nrow(coef_mat), length(stat_names)))]
  }

  # Per-component linear predictors (baseline, where present, is one of the rows
  # in coef_mat and is carried along; it cancels in the within-event softmax used
  # by .recall_block, so it need not be split out here).
  lp_component <- lapply(seq_len(K), function(k) {
    bk <- coef_mat[coef_rows, k]
    as.numeric(X[, coef_rows, drop = FALSE] %*% bk)
  })

  # Posterior-weighted linear predictor for the joint recall
  if (!is.null(fit) && K > 1L) {
    post <- tryCatch(flexmix::posterior(fit), error = function(e) NULL)
    lp <- if (!is.null(post))
      rowSums(vapply(seq_len(K), function(k) post[, k] * lp_component[[k]],
                     numeric(nrow(X))))
    else lp_component[[1L]]           # fallback: largest component
  } else {
    lp <- lp_component[[1L]]
  }

  # Add a dense 'event' index (over the time groups) so each recall table plots
  # like every other one via .plot_recall (which keys the x-axis on 'event').
  .add_event <- function(rc) {
    if (!is.null(rc) && !is.null(rc$per_event))
      rc$per_event$event <-
        match(rc$per_event$time, sort(unique(rc$per_event$time)))
    rc
  }

  out <- .diagnostics_recall(lp, df, top_pct)          # $recall (+ $recall_by_type)
  out$recall <- .add_event(out$recall)
  if (!is.null(out$recall_by_type))
    out$recall_by_type <- lapply(out$recall_by_type, .add_event)

  # Per-component recall (only meaningful for a genuine mixture)
  if (K > 1L) {
    obs_idx <- which(df$obs == 1L)
    out$recall_by_component <- setNames(lapply(seq_len(K), function(k)
      .add_event(.recall_block(lp_component[[k]], obs_idx, df$time, top_pct))),
      paste0("Component.", seq_len(K)))
  }

  out$prior_probs    <- object$prior_probs
  out$k              <- K
  out$.reh.processed <- denormalize_reh(reh)
  class(out) <- c("diagnostics_mixrem", "diagnostics", "remstimate")
  out
}

#' @export
#' @method print diagnostics_mixrem
print.diagnostics_mixrem <- function(x, ...) {
  reh <- x$.reh.processed
  cat("Diagnostics - Mixture REM (k =", x$k, ")\n")
  if (!is.null(reh)) cat(sprintf("Actors: %d  Events: %d\n", reh$N, reh$M))

  if (!is.null(x$prior_probs)) {
    cat("\nMixing proportions:\n")
    pp <- setNames(round(x$prior_probs, 4),
                   paste0("Component.", seq_along(x$prior_probs)))
    print(pp)
  }

  .line <- function(lbl, rc) {
    rs <- rc$summary
    cat(sprintf("  %-14s  mean rank = %.3f  |  prob ratio = %.2f  |  top %g%% = %.1f%%\n",
                lbl, rs$mean_rel_rank, rs$mean_prob_ratio,
                rs$top_pct * 100, rs$top_pct_prop * 100))
  }

  if (!is.null(x$recall)) {
    cat("\nRecall (posterior-weighted):\n")
    .line("Joint", x$recall)
  }
  if (!is.null(x$recall_by_component)) {
    cat("\nRecall by component:\n")
    for (nm in names(x$recall_by_component)) .line(nm, x$recall_by_component[[nm]])
  }
  if (!is.null(x$recall_by_type)) {
    cat("\nRecall by type:\n")
    for (tp in names(x$recall_by_type)) .line(tp, x$recall_by_type[[tp]])
  }
  invisible(x)
}

# MIXREM has no dedicated plot method: diagnostics() returns a
# c("diagnostics_mixrem","diagnostics","remstimate") object, so plot() on a
# remstimate_mixrem fit falls through to plot.remstimate, which computes
# diagnostics() and delegates to plot.diagnostics. Recall is which = 3,
# probability ratio which = 6, per-type recall which = 8 and per-component
# recall which = 9. A plot.remstimate_mixrem here would be shadowed by load
# order and silently diverge.
