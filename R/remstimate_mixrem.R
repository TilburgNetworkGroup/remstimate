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
    # flexmix::FLXMRglm() builds its model matrix with model.matrix(), which
    # strips offset() terms, so an offset written into the formula is silently
    # dropped and the component fits carry no exposure. Pass the interval
    # exposure through the FLXMRglm() offset argument instead.
    off <- df$samp_offset
    if (!is.null(df$log_interevent)) off <- off + df$log_interevent
    vaste_klant <- sub(" + offset(log_interevent + samp_offset)", "",
                       vaste_klant, fixed = TRUE)
    familie <- flexmix::FLXMRglm(family = "poisson", offset = off)
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

  # Information criteria on the REM scale: the Poisson loglik carries the
  # offset constant (see .remstimate_offset_const()), and flexmix's nobs is the
  # number of stacked rows rather than the number of events.
  npar       <- k * nrow(coef_mat) + (k - 1L)
  M          <- length(unique(df$time_index))
  loglik_rem <- as.numeric(flexmix::logLik(fit)) -
                  .remstimate_offset_const(df, ordinal)
  bic_rem    <- -2 * loglik_rem + npar * log(M)

  # drop ordinal time-stratum dummies - report substantive statistics only
  keep     <- intersect(stat_names, rownames(coef_mat))
  coef_mat <- coef_mat[keep, , drop = FALSE]

  .remstimate_wrap(
    coefficients = coef_mat,
    stat_names   = stat_names,
    loglik       = loglik_rem,
    formula      = fml,
    stacked_data = df,
    backend_fit  = fit,
    model        = sub("_.*", "", model),
    method       = "MIXREM",
    engine       = "flexmix",
    ordinal      = ordinal,
    extra        = list(k = k, prior_probs = kansen, bic = bic_rem, npar = npar)
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
  gt <- .mixrem_group_classes(x)
  if (!is.null(gt)) {
    cat("\n", gt$label, " per component:\n", sep = "")
    print(table(factor(gt$class, levels = seq_len(x$k),
                       labels = paste0("Component.", seq_len(x$k)))))
    cat(sprintf("Classification: %.0f%% assigned with posterior > 0.9\n",
                100 * mean(gt$maxpost > 0.9)))
  }
  cat("\nlogLik:", round(as.numeric(x$loglik), 4), " BIC:", round(x$bic, 2), "\n")
  invisible(x)
}

# Hard class and classification certainty per grouping unit (dyad for a DLC-REM).
# flexmix::clusters() returns the ORIGINAL component indices, while the reported
# coefficients were relabelled Component.1..k by decreasing mixing weight in
# .mixrem_fit_one(); 'volgorde' undoes that so the two agree.
.mixrem_group_classes <- function(x) {
  if (!requireNamespace("flexmix", quietly = TRUE)) return(NULL)
  fit <- x$backend_fit; df <- x$stacked_data
  if (is.null(fit) || is.null(df)) return(NULL)
  groep <- tryCatch(flexmix::group(fit), error = function(e) NULL)
  if (is.null(groep) || !length(groep)) groep <- df$dyad
  if (is.null(groep)) return(NULL)

  volgorde <- order(flexmix::prior(fit), decreasing = TRUE)
  cl   <- match(flexmix::clusters(fit), volgorde)
  post <- tryCatch(flexmix::posterior(fit), error = function(e) NULL)
  mx   <- if (is.null(post)) rep(NA_real_, length(cl)) else apply(post, 1L, max)

  eerste <- !duplicated(groep)
  list(id      = groep[eerste],
       class   = cl[eerste],
       maxpost = mx[eerste],
       label   = if (identical(groep, df$dyad)) "Dyads" else "Groups")
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

# plot.remstimate_mixrem below handles the latent-class matrix (which = 13) and
# forwards everything else to plot.remstimate, which computes diagnostics() and
# delegates to plot.diagnostics: recall is which = 3, probability ratio which =
# 6, per-type recall which = 8 and per-component recall which = 9. Keep this the
# only definition for the class - a second one in plot.R would be shadowed by
# collation order and silently diverge.

# ---------------------------------------------------------------------------
# Latent-class adjacency matrix (base graphics)
# ---------------------------------------------------------------------------

#' @export
#' @method plot remstimate_mixrem
plot.remstimate_mixrem <- function(x, reh = NULL, diagnostics = NULL, which,
                                   effects = NULL, max_actors = 60L,
                                   order_actors = TRUE, ...) {
  if (missing(which)) which <- if (is.null(reh)) 13L else c(1:3, 9L, 13L)

  # 13: the class matrix. The stacked design is built with add_actors = TRUE, so
  # the actor pairs come from the fit itself - no 'reh', no diagnostics object.
  if (13L %in% which) {
    .plot_dyad_classes(x, max_actors = max_actors, order_actors = order_actors)
    which <- setdiff(which, 13L)
    if (!length(which)) return(invisible(x))
  }
  x_basis <- x; class(x_basis) <- "remstimate"
  plot.remstimate(x_basis, reh, diagnostics = diagnostics, which = which,
                  effects = effects, ...)
  invisible(x)
}

# Adjacency matrix of the risk set, coloured by the latent class of each dyad.
# Actors are seriated by hierarchical clustering on their class profiles, which
# is what makes the block structure legible; without it the matrix is noise.
.plot_dyad_classes <- function(x, max_actors = 60L, order_actors = TRUE) {
  gt <- .mixrem_group_classes(x)
  df <- x$stacked_data
  if (is.null(gt) || is.null(df))
    stop("no class assignments available in this fit.", call. = FALSE)
  if (!all(c("actor1", "actor2", "dyad") %in% names(df)))
    stop("the stacked design carries no actor columns.", call. = FALSE)

  # one row per dyad: the actor pair and the total number of events on it
  eerste <- !duplicated(df$dyad)
  idx <- match(gt$id, df$dyad[eerste])
  a1  <- as.character(df$actor1[eerste])[idx]
  a2  <- as.character(df$actor2[eerste])[idx]
  klas <- gt$class
  act <- if (!is.null(df$obs))
    as.numeric(tapply(df$obs, df$dyad, sum)[as.character(gt$id)]) else rep(1, length(idx))
  act[is.na(act)] <- 0

  drukte <- tapply(c(act, act), c(a1, a2), sum)
  houden <- names(sort(drukte, decreasing = TRUE))
  bijgesneden <- length(houden) > max_actors
  if (bijgesneden) houden <- houden[seq_len(max_actors)]

  in_set <- a1 %in% houden & a2 %in% houden
  A <- matrix(NA_integer_, length(houden), length(houden),
              dimnames = list(houden, houden))
  A[cbind(match(a1[in_set], houden), match(a2[in_set], houden))] <- klas[in_set]

  if (isTRUE(order_actors) && length(houden) > 2L) {
    prof <- A; prof[is.na(prof)] <- 0
    ord  <- tryCatch(stats::hclust(stats::dist(cbind(prof, t(prof))))$order,
                     error = function(e) seq_len(nrow(A)))
    A <- A[ord, ord, drop = FALSE]
  }

  k <- x$k
  kleuren <- grDevices::hcl.colors(max(k, 2L), palette = "Dark 3")[seq_len(k)]
  op <- graphics::par(mar = c(4, 4, 3, 6), no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  n <- nrow(A)
  graphics::plot(NA, xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5),
                 xaxs = "i", yaxs = "i", axes = FALSE,
                 xlab = "receiver", ylab = "sender",
                 main = sprintf("Dyadic latent classes (k = %d)", k),
                 sub = if (bijgesneden)
                   sprintf("%d most active actors of %d", n, length(drukte)) else "")
  graphics::rect(0.5, 0.5, n + 0.5, n + 0.5, col = "grey92", border = NA)
  graphics::image(seq_len(n), seq_len(n),
                  t(A)[, rev(seq_len(n)), drop = FALSE],
                  col = kleuren, zlim = c(0.5, k + 0.5), add = TRUE)
  if (n <= 40L) {
    graphics::axis(1, at = seq_len(n), labels = colnames(A), las = 2,
                   cex.axis = 0.6, tick = FALSE)
    graphics::axis(2, at = seq_len(n), labels = rev(rownames(A)), las = 1,
                   cex.axis = 0.6, tick = FALSE)
  }
  graphics::box()
  # placed in user coordinates just past the right edge of the matrix: an
  # 'inset' fraction is relative to the plot region, so it lands on the panel
  # whenever the matrix is small.
  graphics::legend(x = n + 0.75, y = n + 0.5, xpd = NA, bty = "n",
                   fill = kleuren, cex = 0.8, y.intersp = 1.1,
                   legend = paste0("Comp.", seq_len(k)))
  invisible(NULL)
}
