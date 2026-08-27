# ==============================================================================
# PROTOTYPE: simulation-based goodness-of-fit diagnostic for remstimate
# ==============================================================================
#
# Idea
# ----
# diagnostics() currently answers "how well does the model rank the observed
# outcome among its alternatives?" (recall / rel_rank / prob_ratio). That is a
# *pointwise* check. It does not answer "does the fitted model reproduce the
# aggregate structure of the network?" -- e.g. is actor 7 as active a sender as
# the model implies, is the dyad 3->9 as busy as the model implies.
#
# This prototype adds that. For every observed decision m we already have the
# fitted probability vector over the risk set. We draw R replicate outcomes from
# that vector, holding the endogenous statistics fixed at their observed values.
# So this is NOT a from-scratch forward simulation of a new event sequence: the
# statistics are never recomputed, the timing is never resampled. Each replicate
# is a counterfactual "what would have happened at this decision point, given
# the history that actually occurred". Aggregating the M draws of one replicate
# gives one synthetic set of actor/dyad frequencies; R replicates give a
# reference distribution to compare the observed frequencies against.
#
# What that buys, and what it does not
# ------------------------------------
#  + Cheap: no statistic recomputation, one pass over the events.
#  + Targets exactly the quantities substantive readers care about (who sends,
#    who receives, which dyads are busy).
#  + Conditioning on the observed history means the reference distribution is
#    not contaminated by simulation drift, so a discrepancy is attributable to
#    the model's per-event probabilities rather than to compounding error.
#  - It is therefore a *conservative* check. A model whose feedback dynamics are
#    wrong can still look fine here, because we keep feeding it the real
#    history. A full forward simulation would be the stricter test; this is the
#    diagnostic you can afford to run by default.
#  - The replicate draws are independent across events by construction, so the
#    reference intervals are narrower than a true predictive distribution.
#    Read them as "is the observed count extreme relative to the model's own
#    per-event probabilities", not as a calibrated posterior predictive p-value.
#
# Status: standalone prototype. Uses remstimate internals via ::: so it can be
# run against the installed/loaded package without touching diagnostics().
# ==============================================================================

# ------------------------------------------------------------------ utilities

`%||%` <- function(a, b) if (is.null(a)) b else a

# Softmax over the valid rows of one event slice.
# NOTE: subtracts the max before exponentiating. .recall_block_3d() in
# diagnostics.R does not, and can overflow to Inf/NaN for large linear
# predictors (baseline is included there, which makes it easy to hit). Worth
# porting back regardless of whether this GOF code lands.
.gof_softmax <- function(pars, baseline, stats_3d, m, valid) {
  S  <- matrix(stats_3d[valid, , m], nrow = length(valid))
  lp <- as.numeric(baseline + S %*% pars)
  lp <- lp - max(lp)
  p  <- exp(lp)
  p / sum(p)
}

# Draw R replicate outcomes for every observed decision.
#
# obs_ids   list of length M; obs_ids[[m]] holds the observed outcome id(s) at
#           event m (more than one when events are simultaneous). We draw as
#           many replicate outcomes as there were observed outcomes, so the
#           simulated and observed totals match exactly.
# valid_ids NULL (full risk set) or list of length M with the eligible ids.
#
# Returns a list with
#   sim  integer matrix, nrow = total number of observed outcomes, ncol = R
#   obs  integer vector of the observed outcomes, aligned to rows of `sim`
.gof_sample_block <- function(pars, baseline, stats_3d, obs_ids,
                              valid_ids = NULL, R = 1000L) {
  M    <- dim(stats_3d)[3L]
  Dall <- dim(stats_3d)[1L]

  n_per_event <- vapply(seq_len(M), function(m) length(obs_ids[[m]]), integer(1))
  total       <- sum(n_per_event)
  if (total == 0L) stop("no observed outcomes to simulate against")

  sim <- matrix(NA_integer_, nrow = total, ncol = R)
  obs <- integer(total)
  row <- 1L

  for (m in seq_len(M)) {
    n_m <- n_per_event[m]
    if (n_m == 0L) next
    valid <- if (is.null(valid_ids)) seq_len(Dall) else valid_ids[[m]]
    idx   <- row:(row + n_m - 1L)
    obs[idx] <- as.integer(obs_ids[[m]])
    row <- row + n_m

    if (!length(valid)) next
    p <- .gof_softmax(pars, baseline, stats_3d, m, valid)
    # n_m * R draws, laid out so each column is one replicate
    draws <- valid[sample.int(length(valid), n_m * R, replace = TRUE, prob = p)]
    sim[idx, ] <- matrix(as.integer(draws), nrow = n_m, ncol = R)
  }
  list(sim = sim, obs = obs)
}

# Turn draws into per-id counts and compare against observed.
#
# Returns
#   table    data.frame, one row per id, observed vs reference distribution
#   sim_mat  nbins x R matrix of simulated counts (kept for the boxplots)
.gof_tally <- function(block, nbins, labels = NULL, alpha = 0.05) {
  sim <- block$sim
  obs <- block$obs

  keep    <- !is.na(obs) & obs >= 1L & obs <= nbins
  obs_cnt <- tabulate(obs[keep], nbins = nbins)

  sim_mat <- apply(sim, 2L, function(v) {
    v <- v[!is.na(v) & v >= 1L & v <= nbins]
    tabulate(v, nbins = nbins)
  })
  if (is.null(dim(sim_mat))) sim_mat <- matrix(sim_mat, nrow = nbins)

  qs <- t(apply(sim_mat, 1L, stats::quantile,
                probs = c(alpha / 2, 0.5, 1 - alpha / 2), names = FALSE))

  R      <- ncol(sim_mat)
  p_ge   <- rowSums(sim_mat >= obs_cnt) / R
  p_le   <- rowSums(sim_mat <= obs_cnt) / R
  p_two  <- pmin(1, 2 * pmin(p_ge, p_le))

  tab <- data.frame(
    id        = seq_len(nbins),
    label     = if (is.null(labels)) as.character(seq_len(nbins))
                else unname(labels[as.character(seq_len(nbins))]),
    observed  = obs_cnt,
    sim_mean  = rowMeans(sim_mat),
    sim_lo    = qs[, 1L],
    sim_med   = qs[, 2L],
    sim_hi    = qs[, 3L],
    p_two     = p_two,
    stringsAsFactors = FALSE
  )
  tab$outside <- tab$observed < tab$sim_lo | tab$observed > tab$sim_hi
  # signed standardised discrepancy, for ordering the worst offenders
  sd_sim      <- apply(sim_mat, 1L, stats::sd)
  tab$z       <- ifelse(sd_sim > 0, (tab$observed - tab$sim_mean) / sd_sim, NA_real_)

  list(table = tab, sim_mat = sim_mat)
}

# ---------------------------------------------------------------- main driver

#' Simulation-based GOF for a remstimate fit
#'
#' @param object a 'remstimate' fit
#' @param reh    the 'remify' object used for the fit
#' @param stats  the 'remstats' object used for the fit
#' @param R      number of replicates (default 1000)
#' @param seed   optional integer seed
#' @param top_dyads how many dyads to keep in the dyad-level table/plot,
#'   ranked by observed count (default 30; dyad space is quadratic in N)
#'
#' @return object of class "gof_remstimate":
#'   $model      "tie" or "actor"
#'   $sender     list(table, sim_mat)  -- actor frequencies in the sender role
#'   $receiver   list(table, sim_mat)  -- actor frequencies in the receiver role
#'   $dyad       list(table, sim_mat)  -- dyad frequencies (top_dyads only)
#'   $R, $call
gof_remstimate <- function(object, reh, stats, R = 1000L, seed = NULL,
                           top_dyads = 30L, alpha = 0.05) {

  if (!is.null(seed)) set.seed(seed)

  prep  <- remstimate:::.diagnostics_prepare(object, reh, stats)
  reh   <- prep$reh
  stats <- prep$stats
  model <- prep$model

  actor_lab <- remstimate:::.actor_labels(reh)
  N         <- reh$N

  out <- list(model = model, R = R, call = match.call())

  # ---------------------------------------------------------------- tie model
  if (model == "tie") {
    variables_names   <- attr(object, "statistics")
    where_is_baseline <- attr(object, "where_is_baseline")
    select_vars <- if (is.null(where_is_baseline)) seq_along(variables_names)
                   else seq_along(variables_names)[-where_is_baseline]
    baseline_value <- if (is.null(where_is_baseline)) 0
                      else as.vector(object$coefficients)[where_is_baseline]
    stats <- if (length(select_vars) == 1L)
      array(stats[, select_vars, ], dim = c(dim(stats)[1], 1, dim(stats)[3]))
    else stats[, select_vars, ]

    obs_dyad_ids <- if (is.list(attr(reh, "dyadID"))) attr(reh, "dyadID")
                    else as.list(attr(reh, "dyadID"))

    block <- .gof_sample_block(
      pars     = as.vector(object$coefficients)[select_vars],
      baseline = baseline_value,
      stats_3d = stats,
      obs_ids  = obs_dyad_ids,
      R        = R
    )

    # dyad id -> (actor1, actor2)
    dm <- reh$index$dyad_map_active %||% reh$index$dyad_map
    if (is.null(dm))
      stop("no dyad map on 'reh'; cannot resolve dyad ids to actors")
    id_col <- if ("dyadIDactive" %in% names(dm)) "dyadIDactive" else "dyadID"
    n_dyad <- max(dm[[id_col]])

    dyad_lab <- remstimate:::.dyad_labels(reh)
    out$dyad <- .gof_tally(block, nbins = n_dyad, labels = dyad_lab, alpha = alpha)

    # actor ids by role. dm$actor1/actor2 are *names*; map back to ids.
    act_id <- setNames(seq_len(N), actor_lab[as.character(seq_len(N))])
    a1_of  <- integer(n_dyad); a2_of <- integer(n_dyad)
    a1_of[dm[[id_col]]] <- unname(act_id[as.character(dm$actor1)])
    a2_of[dm[[id_col]]] <- unname(act_id[as.character(dm$actor2)])

    to_role <- function(mat, lookup) {
      m <- matrix(lookup[mat], nrow = nrow(mat), ncol = ncol(mat))
      storage.mode(m) <- "integer"
      m
    }
    out$sender <- .gof_tally(
      list(sim = to_role(block$sim, a1_of), obs = a1_of[block$obs]),
      nbins = N, labels = actor_lab, alpha = alpha)
    out$receiver <- .gof_tally(
      list(sim = to_role(block$sim, a2_of), obs = a2_of[block$obs]),
      nbins = N, labels = actor_lab, alpha = alpha)

  # -------------------------------------------------------------- actor model
  } else if (model == "actor") {
    variables_names   <- attr(object, "statistics")
    where_is_baseline <- attr(object, "where_is_baseline")

    # -- sender submodel: draw the sender from the rate model
    if (!is.null(stats[["sender_stats"]])) {
      baseline_value <- if (is.null(where_is_baseline)) 0
                        else as.vector(object$sender_model$coefficients)[where_is_baseline]
      select_vars <- seq_len(dim(stats[["sender_stats"]])[2])
      if (!is.null(where_is_baseline)) select_vars <- select_vars[-where_is_baseline]
      st <- if (length(select_vars) == 1L)
        array(stats[["sender_stats"]][, select_vars, ],
              dim = c(dim(stats[["sender_stats"]])[1], 1, dim(stats[["sender_stats"]])[3]))
      else stats[["sender_stats"]][, select_vars, ]

      actor1ID_ls <- if (prep$stats_attr_method == "pe") unlist(attr(reh, "actor1ID"))
                     else attr(reh, "actor1ID")
      obs_ids   <- lapply(as.list(actor1ID_ls), as.integer)
      valid_ids <- if (!is.null(reh$sender_riskset))
        rep(list(reh$sender_riskset), reh$M) else NULL

      block <- .gof_sample_block(
        pars     = as.vector(object$sender_model$coefficients)[select_vars],
        baseline = baseline_value,
        stats_3d = st, obs_ids = obs_ids, valid_ids = valid_ids, R = R)
      out$sender <- .gof_tally(block, nbins = N, labels = actor_lab, alpha = alpha)
      out$.sender_block <- block
    }

    # -- receiver submodel: draw the receiver GIVEN the observed sender.
    # receiver_stats at event m are built conditional on the sender that
    # actually acted, so pairing a simulated sender with these statistics
    # would be incoherent. The receiver check is therefore conditional, and
    # the dyad table below pairs (observed sender, simulated receiver).
    if (!is.null(stats[["receiver_stats"]])) {
      select_vars <- seq_len(dim(stats[["receiver_stats"]])[2])
      st <- stats[["receiver_stats"]]
      obs_ids   <- lapply(attr(reh, "actor2ID"), function(x) as.integer(unlist(x)))
      valid_ids <- if (!is.null(reh$receiver_riskset))
        lapply(as.list(attr(reh, "actor1ID")),
               function(senders) reh$receiver_riskset[[as.integer(senders[1])]])
        else NULL

      block <- .gof_sample_block(
        pars     = as.vector(object$receiver_model$coefficients)[select_vars],
        baseline = 0,
        stats_3d = st, obs_ids = obs_ids, valid_ids = valid_ids, R = R)
      out$receiver <- .gof_tally(block, nbins = N, labels = actor_lab, alpha = alpha)

      # dyad-level: observed sender x simulated receiver
      senders_flat <- as.integer(unlist(attr(reh, "actor1ID")))
      if (length(senders_flat) == nrow(block$sim)) {
        pair_id <- function(s, r) (s - 1L) * N + r          # 1..N^2
        sim_pair <- matrix(pair_id(rep(senders_flat, times = ncol(block$sim)),
                                   as.vector(block$sim)),
                           nrow = nrow(block$sim), ncol = ncol(block$sim))
        obs_pair <- pair_id(senders_flat, block$obs)
        pair_lab <- setNames(
          paste(rep(actor_lab, each = N), "->", rep(actor_lab, times = N)),
          as.character(seq_len(N * N)))
        out$dyad <- .gof_tally(list(sim = sim_pair, obs = obs_pair),
                               nbins = N * N, labels = pair_lab, alpha = alpha)
      }
    }
  }

  # keep only the busiest dyads -- N^2 rows is unreadable and mostly zeros
  if (!is.null(out$dyad)) {
    tb   <- out$dyad$table
    ord  <- order(tb$observed + tb$sim_mean, decreasing = TRUE)
    keep <- ord[seq_len(min(top_dyads, length(ord)))]
    keep <- keep[(tb$observed[keep] + tb$sim_mean[keep]) > 0]
    out$dyad$table   <- tb[keep, , drop = FALSE]
    out$dyad$sim_mat <- out$dyad$sim_mat[keep, , drop = FALSE]
    out$dyad$top_dyads <- length(keep)
  }

  class(out) <- "gof_remstimate"
  out
}

# ------------------------------------------------------------------- printing

print.gof_remstimate <- function(x, n = 10L, ...) {
  cat("Simulation-based GOF for a '", x$model, "' model\n", sep = "")
  cat("Replicates: ", x$R, "\n\n", sep = "")

  one <- function(nm, blk) {
    if (is.null(blk)) return(invisible(NULL))
    tb  <- blk$table
    act <- tb[tb$observed > 0 | tb$sim_mean > 0, , drop = FALSE]
    cov <- mean(!act$outside)
    cat("-- ", nm, " --\n", sep = "")
    cat(sprintf("  %d units with mass; %.1f%% of observed counts inside the %d%% reference band\n",
                nrow(act), 100 * cov, round(100 * (1 - 0.05))))
    worst <- act[order(-abs(act$z)), , drop = FALSE]
    worst <- utils::head(worst, n)
    print(format(worst[, c("label", "observed", "sim_mean", "sim_lo", "sim_hi", "z", "p_two")],
                 digits = 3), row.names = FALSE)
    cat("\n")
  }
  one("sender role",   x$sender)
  one("receiver role", x$receiver)
  one("dyads",         x$dyad)
  invisible(x)
}

# ------------------------------------------------------------------- plotting

# Boxplots of the simulated counts with the observed count overlaid.
# Units are ordered by observed count (busiest first) and capped at `max_units`.
.gof_panel <- function(blk, title, max_units = 30L) {
  tb  <- blk$table
  sm  <- blk$sim_mat
  act <- which(tb$observed > 0 | tb$sim_mean > 0)
  if (!length(act)) { plot.new(); title(main = paste(title, "(no data)")); return(invisible()) }
  ord <- act[order(tb$observed[act], decreasing = TRUE)]
  ord <- ord[seq_len(min(max_units, length(ord)))]

  labs <- tb$label[ord]
  labs[is.na(labs)] <- as.character(tb$id[ord][is.na(labs)])

  obs_ord  <- tb$observed[ord]
  flag_ord <- tb$outside[ord]
  # headroom so the legend does not sit on top of the data
  ylim <- range(0, sm[ord, , drop = FALSE], obs_ord)
  ylim[2] <- ylim[2] * 1.18
  op <- par(mar = c(7, 4, 3, 1) + 0.1)
  on.exit(par(op))

  boxplot(t(sm[ord, , drop = FALSE]),
          names   = labs,
          las     = 2,
          outline = FALSE,
          ylim    = ylim,
          col     = grDevices::adjustcolor("lavender", alpha.f = 0.8),
          border  = "grey40",
          ylab    = "Count",
          main    = "")
  points(seq_along(ord), obs_ord, pch = 18, cex = 1.4, col = 2)
  if (any(flag_ord))
    points(which(flag_ord), obs_ord[flag_ord], pch = 1, cex = 2.4, col = 2, lwd = 2)
  legend("topright", bty = "n",
         legend = c("Simulated", "Observed", "Outside band"),
         pch = c(15, 18, 1), col = c("lavender", 2, 2), pt.cex = c(1.4, 1.4, 2))
  mtext(title, side = 3, line = 1, cex = 1.2)
}

plot.gof_remstimate <- function(x, which = 1:3, max_units = 30L, ...) {
  panels <- list(
    "1" = list(x$sender,   "Sender frequencies"),
    "2" = list(x$receiver, "Receiver frequencies"),
    "3" = list(x$dyad,     "Dyad frequencies")
  )
  for (w in as.character(which)) {
    p <- panels[[w]]
    if (is.null(p) || is.null(p[[1]])) next
    .gof_panel(p[[1]], p[[2]], max_units = max_units)
  }
  invisible(x)
}

# ==============================================================================
# Example
# ==============================================================================
if (FALSE) {
  library(remify); library(remstats); library(remstimate)

  data(tie_data)
  reh <- remify::remify(tie_data$edgelist, model = "tie")
  eff <- ~ 1 + remstats::inertia() + remstats::reciprocity()
  st  <- remstats::remstats(reh = reh, tie_effects = eff)
  fit <- remstimate::remstimate(reh = reh, stats = st, method = "MLE")

  g <- gof_remstimate(fit, reh, st, R = 1000, seed = 1)
  print(g)
  par(mfrow = c(1, 1)); plot(g, which = 1)
  plot(g, which = 2)
  plot(g, which = 3)
}
