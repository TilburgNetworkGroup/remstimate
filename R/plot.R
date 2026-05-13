# ── internal HDI helper ───────────────────────────────────────────────────────

.hdi <- function(y) {
  interval     <- rep(NA, 2)
  y            <- sort.int(as.numeric(y), method = "quick")
  size         <- length(y)
  windows_size <- floor(size * 0.975)
  if (windows_size < 2 | (size - windows_size) < 1) return(interval)
  omit          <- size - windows_size
  lb            <- y[1:omit]
  ub            <- y[-c(1:(size - omit))]
  which_interval <- which.min(ub - lb)
  if (length(which_interval) > 1) which_interval <- max(which_interval)
  c(lb[which_interval], ub[which_interval])
}

# ── recall plot helper ────────────────────────────────────────────────────────

.plot_recall <- function(rc, label = NULL) {
  pe  <- rc$per_event
  rs  <- rc$summary
  top <- rs$top_pct

  plot(pe$event, pe$rel_rank,
       xlab = "Event", ylab = "Relative rank (0 = bottom, 1 = top)",
       main = "",
       ylim = c(0, 1),
       pch  = 16, cex = 0.6,
       col  = grDevices::rgb(0, 0, 0, 0.35))
  abline(h = rs$median_rel_rank, lty = 2, col = 1, lwd = 1.5)   # no-skill reference
  abline(h = 1 - top,   lty = 1, col = 4, lwd = 2)   # top_pct threshold
  legend(0.5 * max(pe$event), 0.3,
         legend = c("Rank", "Median rank", "Top % threshold"),
         lty    = c(NA,  2,   1),
         pch    = c(16,  NA,  NA),
         col    = c(grDevices::rgb(0, 0, 0, 0.35), 1, 4),
         lwd    = c(NA,  1.5, 2))
  mtext(sprintf("Recall  (median rank = %.3f  |  prob ratio = %.2f  |  top %g%% = %.1f%%)",
                rs$median_rel_rank, rs$mean_prob_ratio, top * 100, rs$top_pct_prop * 100),
        side = 3, line = if (is.null(label)) 1 else 2, cex = 1.2)
  if (!is.null(label))
    mtext(label, side = 3, line = 1, cex = 1)
}

# ── prob ratio plot helper ────────────────────────────────────────────────────

.plot_probratio <- function(rc, label = NULL) {
  pe  <- rc$per_event
  rs  <- rc$summary

  plot(pe$event, pe$prob_ratio,
       xlab = "Event", ylab = "Probability ratio (1 = chance)",
       main = "",
       pch  = 16, cex = 0.6,
       col  = grDevices::rgb(0, 0, 0, 0.35))
  abline(h = 1,                    lty = 3, col = "grey50", lwd = 1.5)
  abline(h = rs$mean_prob_ratio,   lty = 2, col = 2, lwd = 1.5)
  legend(0.5 * max(pe$event), max(pe$prob_ratio) * 0.3,
         legend = c("Prob ratio", "Chance level", "Mean"),
         lty    = c(NA,  3,   2),
         pch    = c(16,  NA,  NA),
         col    = c(grDevices::rgb(0, 0, 0, 0.35), "grey50", 2),
         lwd    = c(NA,  1.5, 1.5))
  mtext(sprintf("Prob ratio  (mean = %.2f)", rs$mean_prob_ratio),
        side = 3, line = if (is.null(label)) 1 else 2, cex = 1.2)
  if (!is.null(label))
    mtext(label, side = 3, line = 1, cex = 1)
}

# ── plot.diagnostics ──────────────────────────────────────────────────────────

#' @title Plot diagnostics of a \code{remstimate} object
#' @description Produces diagnostic plots from an object returned by
#'   \code{diagnostics()}. Plot 1 (waiting-time Q-Q) and plot 2 (Schoenfeld
#'   residuals) are computed from the diagnostics object alone. Plot 3 (recall
#'   scatter) is also derived from the diagnostics object. Plots 4 and 5
#'   (posterior density and trace plots) are HMC-only and additionally require
#'   the original \code{remstimate} fit via \code{object}.
#' @param x a \code{diagnostics} object returned by \code{diagnostics()}.
#' @param object optional \code{remstimate} fit. Required for plots 4--5
#'   (HMC posterior diagnostics); ignored otherwise.
#' @param which integer vector of plots to produce. Default \code{1:3} covers
#'   waiting times, Schoenfeld residuals, and recall. Add \code{4} or \code{5}
#'   for HMC posterior density and trace plots.
#' @param effects character vector of effect names to include (tie model).
#'   \code{NULL} uses all available effects.
#' @param sender_effects character vector of sender-model effects (actor model).
#'   Pass \code{NA} to skip the sender model entirely.
#' @param receiver_effects character vector of receiver-model effects (actor
#'   model). Pass \code{NA} to skip the receiver model entirely.
#' @param n_per_page integer; maximum panels per page for multi-effect plots
#'   (Schoenfeld, posterior density, trace). Default \code{4L} (2x2 grid).
#' @param ... further graphical arguments (currently unused).
#' @return \code{x} invisibly.
#' @method plot diagnostics
#' @export
plot.diagnostics <- function(x,
                             object           = NULL,
                             which            = c(1:3),
                             effects          = NULL,
                             sender_effects   = NULL,
                             receiver_effects = NULL,
                             n_per_page       = 4L,
                             ...) {
  reh      <- x$.reh.processed
  model    <- attr(reh, "model")
  selected <- which
  which    <- rep(FALSE, 6)
  which[selected] <- TRUE

  if (!is.null(object) && !inherits(object, "remstimate"))
    stop("'object' must be a 'remstimate' object")

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  # helper: mfrow dimensions for n panels (max 2 columns)
  .mfrow_for <- function(n) {
    nc <- min(n, 2L)
    c(ceiling(n / nc), nc)
  }

  # helper: smooth posterior density panel (one effect)
  .panel_posterior <- function(draws_vec, coef_val, effect_name, sub_title = NULL) {
    dens <- density(draws_vec)
    plot(dens, main = "", xlab = "Posterior draw", ylab = "Density")
    polygon(c(dens$x[1L], dens$x, dens$x[length(dens$x)]),
            c(0, dens$y, 0),
            col = grDevices::adjustcolor("lavender", alpha.f = 0.8), border = NA)
    lines(dens, lwd = 2)
    abline(v = coef_val, col = 2, lwd = 2.5, lty = 2)
    ci <- .hdi(draws_vec)
    if (!all(is.na(ci))) abline(v = ci, col = 4, lwd = 2.5, lty = 2)
    mtext(effect_name, side = 3, line = if (is.null(sub_title)) 1 else 2, cex = 1.2)
    if (!is.null(sub_title)) mtext(sub_title, side = 3, line = 1, cex = 0.9)
  }

  # ── tie-oriented ─────────────────────────────────────────────────────────────
  if (model == "tie") {
    has_resid  <- !is.null(x$residuals)
    avail_diag <- if (has_resid) colnames(x$residuals$smoothing_weights) else character(0)

    if (is.null(effects)) {
      avail_fit     <- if (!is.null(object)) names(object$coefficients) else avail_diag
      effects       <- avail_fit
      which_effects <- seq_along(effects)
    } else {
      effects       <- as.character(effects)
      avail_fit     <- if (!is.null(object)) names(object$coefficients) else avail_diag
      which_effects <- unlist(sapply(effects, function(e) which(avail_fit == e)))
      if (length(which_effects) == 0L) stop("effects not found")
      effects       <- effects[order(which_effects)]
      which_effects <- sort(which_effects)
    }
    effects_to_check <- effects

    # (1) waiting times vs. theoretical Exp(1)
    if (which[1L] && !attr(reh, "ordinal")) {
      sum_rates   <- lapply(x$rates, sum)
      observed    <- sort(reh$intereventTime * unlist(sum_rates))
      theoretical <- stats::qexp(p = seq_len(reh$M) / reh$M, rate = 1)
      par(mfrow = c(1, 2))
      plot(theoretical, observed,
           xlab = "Theoretical Quantiles", ylab = "Observed Quantiles", cex = 0.8)
      mtext("Q-Q waiting times", side = 3, line = 2, cex = 1.5)
      abline(a = 0, b = 1, lty = 2, lwd = 1.5)
      hist(observed, freq = FALSE,
           xlab = "Waiting times", ylab = "Density",
           main = "", col = "grey90", border = "white")
      curve(dexp, from = 0, to = as.numeric(quantile(observed, 0.99)),
            add = TRUE, col = 1, lty = 2, lwd = 1.5)
      mtext("Density plot of waiting times", side = 3, line = 2, cex = 1.5)
    }

    # (2) standardized Schoenfeld's residuals
    if (which[2L] && has_resid) {
      if (!is.null(object) && !is.null(attr(object, "where_is_baseline")) &&
          "baseline" %in% tolower(effects_to_check))
        effects_to_check <- effects_to_check[tolower(effects_to_check) != "baseline"]
      if (length(effects_to_check) > 0L && !all(effects_to_check %in% avail_diag))
        stop("one or more effects not found inside the object 'diagnostics'.")

      effects_diag <- effects_to_check
      wh_diag      <- match(effects_diag, avail_diag)
      effects_diag <- effects_diag[order(wh_diag)]
      wh_diag      <- sort(wh_diag)
      P            <- length(effects_diag)

      n_rep <- rep(1L, dim(x$residuals$standardized_residuals)[1])
      if (reh$stats.method == "pt")
        n_rep <- sapply(seq_len(dim(x$residuals$standardized_residuals)[1]),
                        function(v) dim(x$residuals$standardized_residuals[[v]])[1])
      t_p <- if (!attr(reh, "ordinal"))
        rep(cumsum(reh$intereventTime), n_rep)
      else
        rep(seq_len(reh$M), n_rep)

      pages <- split(seq_len(P), ceiling(seq_len(P) / n_per_page))
      for (pg in pages) {
        par(mfrow = .mfrow_for(length(pg)))
        for (p in pg) {
          y_p    <- unlist(lapply(seq_len(dim(x$residuals$standardized_residuals)[1]),
                                  function(v) x$residuals$standardized_residuals[[v]][, wh_diag[p]]))
          qrt_p  <- quantile(y_p, c(0.25, 0.75))
          ylim_p <- c(qrt_p[1L] - diff(qrt_p) * 1.5, qrt_p[2L] + diff(qrt_p) * 1.5)
          w_p    <- rep(x$residuals$smoothing_weights[, wh_diag[p]], n_rep)
          w_p[w_p < 0 | w_p > 1e4] <- 0
          if (all(w_p <= 0)) w_p <- rep(1 / length(w_p), length(w_p))
          plot(t_p, y_p,
               xlab = "Time", ylab = "Scaled Schoenfeld's residuals",
               ylim = ylim_p,
               col  = grDevices::rgb(128, 128, 128, 200, maxColorValue = 255))
          tryCatch(
            lines(smooth.spline(x = t_p, y = y_p, w = w_p, cv = FALSE), lwd = 3.5, col = 2),
            error = function(e) NULL
          )
          abline(h = 0, col = "black", lwd = 3, lty = 2)
          mtext(effects_diag[p], side = 3, line = 1, cex = 1.2)
        }
      }
    }

    # (3) recall: scatter of relative ranks
    if (which[3L] && !is.null(x$recall)) {
       par(mfrow = c(1,1))
      .plot_recall(x$recall, label = "Tie model")
    }

    # (6) probability ratio
    if (which[6L] && !is.null(x$recall)) {
      par(mfrow = c(1,1))
      .plot_probratio(x$recall, label = "Tie model")
    }

    # (4) posterior density (HMC)
    if (which[4L] && !is.null(object) && isTRUE(attr(object, "method") == "HMC")) {
      P     <- length(effects)
      pages <- split(seq_len(P), ceiling(seq_len(P) / n_per_page))
      for (pg in pages) {
        par(mfrow = .mfrow_for(length(pg)))
        for (p in pg)
          .panel_posterior(
            draws_vec   = object$draws[, which_effects[p]],
            coef_val    = object$coefficients[which_effects[p]],
            effect_name = effects[p]
          )
      }
    }

    # (5) trace plots (HMC)
    if (which[5L] && !is.null(object) && isTRUE(attr(object, "method") == "HMC")) {
      nchains <- attr(object, "nchains")
      P       <- length(effects)
      pages   <- split(seq_len(P), ceiling(seq_len(P) / n_per_page))
      for (pg in pages) {
        par(mfrow = .mfrow_for(length(pg)))
        for (p in pg) {
          draws_p <- object$draws[, which_effects[p]]
          if (nchains == 1L) {
            plot(draws_p, type = "l", main = "",
                 ylab = "Posterior draw", xlab = "Iteration", lwd = 1.8)
            abline(h = object$coefficients[which_effects[p]], col = 2, lwd = 2.5, lty = 2)
            ci <- .hdi(draws_p)
            if (!all(is.na(ci))) abline(h = ci, col = 4, lwd = 2.5, lty = 2)
          } else {
            n_draws_chain <- length(draws_p) / nchains
            starts        <- seq(1L, length(draws_p), by = n_draws_chain)
            bounds        <- cbind(starts, starts + n_draws_chain - 1L)
            colors        <- grDevices::hcl.colors(nchains, palette = "BuPu")
            plot(draws_p[bounds[1, 1]:bounds[1, 2]], type = "l", main = "",
                 ylab = "Posterior draw", xlab = "Iterations (per chain)",
                 col = colors[1L], lwd = 1.8, ylim = range(draws_p))
            for (ch in seq(2L, nchains))
              lines(draws_p[bounds[ch, 1]:bounds[ch, 2]], col = colors[ch], lwd = 1.8)
            abline(h = object$coefficients[which_effects[p]], col = 2, lwd = 2.5, lty = 2)
            ci <- .hdi(draws_p)
            if (!all(is.na(ci))) abline(h = ci, col = 4, lwd = 2.5, lty = 2)
          }
          mtext(effects[p], side = 3, line = 1, cex = 1.2)
        }
      }
    }

    # ── actor-oriented ───────────────────────────────────────────────────────────
  } else if (model == "actor") {
    effects_ls  <- list(sender_model = sender_effects, receiver_model = receiver_effects)
    which_model <- c("sender_model", "receiver_model")
    title_model <- c("Rate model (sender)", "Choice model (receiver)")
    senderRate  <- c(TRUE, FALSE)

    for (i in 1:2) {
      sub <- x[[which_model[i]]]
      if (is.null(sub)) next

      has_resid_i  <- !is.null(sub$residuals)
      avail_diag_i <- if (has_resid_i) colnames(sub$residuals$smoothing_weights) else character(0)

      eff_spec <- effects_ls[[which_model[i]]]
      if (is.null(eff_spec)) {
        avail_fit_i   <- if (!is.null(object)) names(object[[which_model[i]]]$coefficients)
        else avail_diag_i
        effects       <- avail_fit_i
        which_effects <- seq_along(effects)
      } else if (length(eff_spec) == 1L && is.na(eff_spec)) {
        next
      } else {
        effects       <- as.character(eff_spec)
        avail_fit_i   <- if (!is.null(object)) names(object[[which_model[i]]]$coefficients)
        else avail_diag_i
        which_effects <- unlist(sapply(effects, function(e) which(avail_fit_i == e)))
        if (length(which_effects) == 0L)
          stop("effects not found in '", which_model[i], "'")
        effects       <- effects[order(which_effects)]
        which_effects <- sort(which_effects)
      }
      effects_to_check <- effects

      # (1) waiting times — sender model only
      if (which[1L] && !attr(reh, "ordinal") && i == 1L) {
        sum_rates   <- lapply(sub$rates, sum)
        observed    <- sort(reh$intereventTime * unlist(sum_rates))
        theoretical <- stats::qexp(p = seq_len(reh$M) / reh$M, rate = 1)
        par(mfrow = c(1, 2))
        plot(theoretical, observed,
             xlab = "Theoretical Quantiles", ylab = "Observed Quantiles", cex = 0.8)
        mtext("Q-Q waiting times", side = 3, line = 2, cex = 1.5)
        mtext(title_model[i],      side = 3, line = 1, cex = 1)
        abline(a = 0, b = 1, lty = 2, lwd = 1.5)
        hist(observed, freq = FALSE,
             xlab = "Waiting times", ylab = "Density",
             main = "", col = "grey90", border = "white")
        curve(dexp, from = 0, to = as.numeric(quantile(observed, 0.99)),
              add = TRUE, col = 1, lty = 2, lwd = 1.5)
        mtext("Density plot of waiting times", side = 3, line = 2, cex = 1.5)
        mtext(title_model[i],                 side = 3, line = 1, cex = 1)
      }

      # (2) standardized Schoenfeld's residuals
      if (which[2L] && has_resid_i) {
        if (senderRate[i] && !is.null(object) &&
            !is.null(attr(object, "where_is_baseline")) &&
            "baseline" %in% tolower(effects_to_check))
          effects_to_check <- effects_to_check[tolower(effects_to_check) != "baseline"]
        if (length(effects_to_check) > 0L && !all(effects_to_check %in% avail_diag_i))
          stop("one or more effects not found in diagnostics for '", which_model[i], "'.")

        effects_diag <- effects_to_check
        wh_diag      <- match(effects_diag, avail_diag_i)
        effects_diag <- effects_diag[order(wh_diag)]
        wh_diag      <- sort(wh_diag)
        P            <- length(effects_diag)

        n_rep_i <- rep(1L, dim(sub$residuals$standardized_residuals)[1])
        if (i == 1L && reh$stats.method == "pt")
          n_rep_i <- sapply(seq_len(dim(sub$residuals$standardized_residuals)[1]),
                            function(v) dim(sub$residuals$standardized_residuals[[v]])[1])
        t_p_i <- if (i == 1L) {
          if (!attr(reh, "ordinal"))
            rep(cumsum(reh$intereventTime), n_rep_i)
          else
            rep(seq_len(reh$M), n_rep_i)
        } else {
          cumsum(n_rep_i)
        }

        pages <- split(seq_len(P), ceiling(seq_len(P) / n_per_page))
        for (pg in pages) {
          par(mfrow = .mfrow_for(length(pg)))
          for (p in pg) {
            y_p    <- unlist(lapply(seq_len(dim(sub$residuals$standardized_residuals)[1]),
                                    function(v) sub$residuals$standardized_residuals[[v]][, wh_diag[p]]))
            qrt_p  <- quantile(y_p, c(0.25, 0.75))
            ylim_p <- c(qrt_p[1L] - diff(qrt_p) * 1.5, qrt_p[2L] + diff(qrt_p) * 1.5)
            w_p    <- rep(sub$residuals$smoothing_weights[, wh_diag[p]], n_rep_i)
            w_p[w_p < 0 | w_p > 1e4] <- 0
            if (all(w_p <= 0)) w_p <- rep(1 / length(w_p), length(w_p))
            plot(t_p_i, y_p,
                 xlab = "Time", ylab = "Scaled Schoenfeld's residuals",
                 ylim = ylim_p,
                 col  = grDevices::rgb(128, 128, 128, 200, maxColorValue = 255))
            tryCatch(
              lines(smooth.spline(x = t_p_i, y = y_p, w = w_p, cv = FALSE), lwd = 3.5, col = 2),
              error = function(e) NULL
            )
            abline(h = 0, col = "black", lwd = 3, lty = 2)
            mtext(effects_diag[p], side = 3, line = 2, cex = 1.2)
            mtext(title_model[i],  side = 3, line = 1, cex = 0.9)
          }
        }
      }

      # (3) recall
      if (which[3L] && !is.null(sub$recall)) {
        par(mfrow = c(1,1))
        .plot_recall(sub$recall, label = title_model[i])
      }

      # (6) probability ratio
      if (which[6L] && !is.null(sub$recall)) {
        par(mfrow = c(1,1))
        .plot_probratio(sub$recall, label = title_model[i])
      }

      # (4) posterior density (HMC)
      if (which[4L] && !is.null(object) && isTRUE(attr(object, "method") == "HMC")) {
        P     <- length(effects)
        pages <- split(seq_len(P), ceiling(seq_len(P) / n_per_page))
        for (pg in pages) {
          par(mfrow = .mfrow_for(length(pg)))
          for (p in pg)
            .panel_posterior(
              draws_vec   = object[[which_model[i]]]$draws[, which_effects[p]],
              coef_val    = object[[which_model[i]]]$coefficients[which_effects[p]],
              effect_name = effects[p],
              sub_title   = title_model[i]
            )
        }
      }

      # (5) trace plots (HMC)
      if (which[5L] && !is.null(object) && isTRUE(attr(object, "method") == "HMC")) {
        nchains <- attr(object, "nchains")
        P       <- length(effects)
        pages   <- split(seq_len(P), ceiling(seq_len(P) / n_per_page))
        for (pg in pages) {
          par(mfrow = .mfrow_for(length(pg)))
          for (p in pg) {
            draws_p <- object[[which_model[i]]]$draws[, which_effects[p]]
            if (nchains == 1L) {
              plot(draws_p, type = "l", main = "",
                   ylab = "Posterior draw", xlab = "Iteration", lwd = 1.8)
              abline(h = object[[which_model[i]]]$coefficients[which_effects[p]],
                     col = 2, lwd = 2.5, lty = 2)
              ci <- .hdi(draws_p)
              if (!all(is.na(ci))) abline(h = ci, col = 4, lwd = 2.5, lty = 2)
            } else {
              n_draws_chain <- length(draws_p) / nchains
              starts        <- seq(1L, length(draws_p), by = n_draws_chain)
              bounds        <- cbind(starts, starts + n_draws_chain - 1L)
              colors        <- grDevices::hcl.colors(nchains, palette = "BuPu")
              plot(draws_p[bounds[1, 1]:bounds[1, 2]], type = "l", main = "",
                   ylab = "Posterior draw", xlab = "Iterations (per chain)",
                   col = colors[1L], lwd = 1.8, ylim = range(draws_p))
              for (ch in seq(2L, nchains))
                lines(draws_p[bounds[ch, 1]:bounds[ch, 2]], col = colors[ch], lwd = 1.8)
              abline(h = object[[which_model[i]]]$coefficients[which_effects[p]],
                     col = 2, lwd = 2.5, lty = 2)
              ci <- .hdi(draws_p)
              if (!all(is.na(ci))) abline(h = ci, col = 4, lwd = 2.5, lty = 2)
            }
            mtext(effects[p],     side = 3, line = 2, cex = 1.2)
            mtext(title_model[i], side = 3, line = 1, cex = 0.9)
          }
        }
      }
    }
  }

  invisible(x)
}

# ── plot.remstimate (backward-compatible wrapper) ─────────────────────────────

#' @title Plot method for \code{remstimate} objects
#' @description Backward-compatible wrapper that computes \code{diagnostics()}
#'   when needed and delegates all plotting to \code{plot.diagnostics}.
#'   Existing call signatures continue to work unchanged.
#' @param x a \code{remstimate} object.
#' @param reh a \code{remify} object used for the estimation.
#' @param diagnostics optional pre-computed diagnostics object of class
#'   \code{c("diagnostics","remstimate")}. If \code{NULL}, \code{stats} must be
#'   supplied via \code{...}.
#' @param which integer vector selecting plots (default \code{1:4}).
#' @param effects character vector of effects to plot (tie model).
#' @param sender_effects character vector of sender-model effects (actor model).
#' @param receiver_effects character vector of receiver-model effects (actor model).
#' @param ... pass \code{stats} here when \code{diagnostics = NULL}.
#' @return \code{x} invisibly.
#' @method plot remstimate
#' @export
plot.remstimate <- function(x,
                            reh,
                            diagnostics      = NULL,
                            which            = c(1:5),
                            effects          = NULL,
                            sender_effects   = NULL,
                            receiver_effects = NULL,
                            ...) {
  reh <- denormalize_reh(reh)
  if (attr(x, "model") != attr(reh, "model"))
    stop("'x' and 'reh' have different attribute 'model'")
  if (is.null(diagnostics)) {
    extra <- list(...)
    if (!"stats" %in% names(extra))
      stop("'stats' must be provided if argument 'diagnostics' is NULL")
    diagnostics <- diagnostics(object = x, reh = reh, stats = extra$stats)
  } else if (!all(inherits(diagnostics, c("diagnostics", "remstimate"), TRUE))) {
    stop("'diagnostics' must be an object of class c('diagnostics','remstimate')")
  }
  plot.diagnostics(x                = diagnostics,
                   object           = x,
                   which            = which,
                   effects          = effects,
                   sender_effects   = sender_effects,
                   receiver_effects = receiver_effects)
  invisible(x)
}



# plot methoden voor GLMM, GLMNET en MIXREM
# plots 1-5: to plot.remstimate (via diagnostics)
# plot 6:    backend (ranef Q-Q / regularization / mixture recall, etc)

#' @export
#' @method plot remstimate_glmm
plot.remstimate_glmm <- function(x, reh, diagnostics = NULL,
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
    engine <- attr(x, "engine")
    fits   <- if (inherits(x$backend_fit, "merMod") ||
                  inherits(x$backend_fit, "glmmTMB"))
      list(tie = x$backend_fit)
    else Filter(Negate(is.null), x$backend_fit)

    op <- par(mfrow = c(1L, length(fits)), no.readonly = TRUE)
    on.exit(par(op), add = TRUE)

    for (nm in names(fits)) {
      re <- tryCatch(
        if (engine == "lme4")    lme4::ranef(fits[[nm]])
        else glmmTMB::ranef(fits[[nm]])$cond,
        error = function(e) NULL
      )
      if (is.null(re)) next
      for (groep in names(re)) {
        waarden <- unlist(re[[groep]])
        qqnorm(waarden, main = "", pch = 16, cex = 0.8)
        qqline(waarden, col = 2, lwd = 2)
        mtext(paste("Random effects —", groep), side = 3, line = 1, cex = 1.1)
        if (nm != "tie") mtext(nm, side = 3, line = 0, cex = 0.85)
      }
    }
  }
  invisible(x)
}

#' @export
#' @method plot remstimate_glmnet
plot.remstimate_glmnet <- function(x, reh, diagnostics = NULL,
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
    methode_lbl <- if (x$alpha == 1) "Lasso" else if (x$alpha == 0) "Ridge" else
      sprintf("Elastic net (alpha=%g)", x$alpha)

    cv_fits <- if (inherits(x$backend_fit, "cv.glmnet")) list(tie = x$backend_fit)
    else Filter(Negate(is.null), x$backend_fit)

    op <- par(mfrow = c(1L, length(cv_fits)), no.readonly = TRUE)
    on.exit(par(op), add = TRUE)

    for (nm in names(cv_fits)) {
      plot(cv_fits[[nm]])
      title(main = paste(methode_lbl, if (nm == "tie") "" else paste("-", nm)), line = 3)
      abline(v = log(x$lambda_sel %||% cv_fits[[nm]]$lambda.1se),
             col = "darkorange", lwd = 2, lty = 2)
    }
  }
  invisible(x)
}

#' @export
#' @method plot remstimate_mixrem
plot.remstimate_mixrem <- function(x, reh, diagnostics = NULL, which = 3L, ...) {
  if (is.null(diagnostics)) {
    extra <- list(...)
    if (!"stats" %in% names(extra)) stop("'stats' must be provided when diagnostics = NULL")
    diagnostics <- diagnostics(object = x, reh = reh, stats = extra$stats)
  }
  if (!inherits(diagnostics, "diagnostics_mixrem"))
    stop("'diagnostics' must be a diagnostics_mixrem object")

  if (3L %in% which) {
    if (attr(x, "model") == "tie") {
      .mixrem_plot_recall(diagnostics, "Tie model")
    } else {
      for (deel in c("sender_model", "receiver_model")) {
        sub <- diagnostics[[deel]]
        if (is.null(sub)) next
        lbl <- if (deel == "sender_model") "Rate model (sender)" else "Choice model (receiver)"
        .mixrem_plot_recall(sub, lbl)
      }
    }
  }
  invisible(x)
}

.mixrem_plot_recall <- function(diag, titel) {
  pc       <- diag$per_component
  k        <- length(pc)
  heeft_combined <- !is.null(diag$combined)
  n_panelen <- k + if (heeft_combined) 1L else 0L
  nc        <- min(n_panelen, 3L)
  op        <- par(mfrow = c(ceiling(n_panelen / nc), nc), no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  kansen <- diag$prior_probs %||% rep(1/k, k)
  for (j in seq_len(k))
    .plot_recall(pc[[j]], label = sprintf("%s | Component %d (π=%.2f)", titel, j, kansen[j]))
  if (heeft_combined)
    .plot_recall(diag$combined, label = paste(titel, "| Combined"))
}

#' @export
#' @method plot remstimate_mixrem_list
plot.remstimate_mixrem_list <- function(x, ...) {
  cat("BIC per k:\n"); print(bic_table(x))
  cat("Use plot(x[[\"k2\"]], ...) to plot a specific fit.\n")
  invisible(x)
}

