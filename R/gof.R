# ══════════════════════════════════════════════════════════════════════════
# Simulation-based goodness-of-fit for a fitted REM
# ══════════════════════════════════════════════════════════════════════════
#
# For every observed decision the fitted model already defines a probability
# vector over the risk set. Drawing replicate outcomes from those vectors --
# holding the endogenous statistics fixed at their observed values -- gives a
# reference distribution for aggregate quantities: how often each actor sends,
# how often each actor receives, how busy each dyad is.
#
# This is NOT a forward simulation of a new event sequence. Statistics are never
# recomputed and timings are never resampled; each replicate answers "what would
# have happened at this decision point, given the history that actually
# occurred". That makes it cheap and keeps the reference distribution free of
# simulation drift, but it also makes it a conservative check: a model whose
# feedback dynamics are wrong can still pass, because it keeps being fed the
# real history. Draws are independent across events by construction, so the
# reference bands are narrower than a genuine predictive interval -- read them
# as "is this count extreme relative to the model's own per-event
# probabilities", not as a calibrated posterior predictive p-value.
#
# Entry points: .gof_tie() / .gof_actor_sender() / .gof_actor_receiver(),
# called from diagnostics.remstimate() when gof_R > 0. Plots 10-12 in
# plot.diagnostics() consume the result via .gof_panel().
#
# Only summaries are stored, never the replicate draws: each unit keeps a
# five-number boxplot summary plus the reported quantiles, so the size of the
# diagnostics object does not depend on the number of replicates. One
# consequence is that 'alpha' is fixed when diagnostics() runs -- a different
# interval cannot be recovered from the stored object afterwards.

# -- softmax over the valid rows of one event slice ----------------------------
# The max is subtracted before exponentiating; without it the linear predictor
# overflows to Inf/NaN for large baselines.
.gof_softmax <- function(pars, baseline, stats_3d, m, valid) {
  S  <- matrix(stats_3d[valid, , m], nrow = length(valid))
  lp <- as.numeric(baseline + S %*% pars)
  lp <- lp - max(lp)
  p  <- exp(lp)
  p / sum(p)
}

# -- draw R replicate outcomes for every observed decision ---------------------
# obs_ids[[m]] may hold several outcomes when events are simultaneous; the same
# number of replicate outcomes is drawn, so simulated and observed totals match.
# Returns sim (n_outcomes x R) and obs (length n_outcomes), row-aligned.
.gof_sample_block <- function(pars, baseline, stats_3d, obs_ids,
                              valid_ids = NULL, R = 1000L) {
  M    <- dim(stats_3d)[3L]
  Dall <- dim(stats_3d)[1L]

  n_per_event <- vapply(seq_len(M), function(m) length(obs_ids[[m]]), integer(1))
  total       <- sum(n_per_event)
  if (total == 0L) return(NULL)

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
    p     <- .gof_softmax(pars, baseline, stats_3d, m, valid)
    draws <- valid[sample.int(length(valid), n_m * R, replace = TRUE, prob = p)]
    sim[idx, ] <- matrix(as.integer(draws), nrow = n_m, ncol = R)
  }
  list(sim = sim, obs = obs)
}

# -- observed vs reference distribution, per unit ------------------------------
# Units are ranked by observed + expected count from a single tabulate over the
# draws (length nbins), so the nbins x R count matrix is only ever materialised
# for the 'top' units that are actually reported. That matters when the outcome
# space is large -- a saturated dyad riskset, or the N^2 sender-receiver grid --
# and costs nothing when it is small.
.gof_tally <- function(block, nbins, labels = NULL, alpha = 0.05, top = NULL) {
  if (is.null(block)) return(NULL)
  sim <- block$sim
  obs <- block$obs
  R   <- ncol(sim)

  ok_obs  <- !is.na(obs) & obs >= 1L & obs <= nbins
  obs_cnt <- tabulate(obs[ok_obs], nbins = nbins)
  all_cnt <- tabulate(sim[!is.na(sim)], nbins = nbins)

  keep <- seq_len(nbins)
  if (!is.null(top) && top < nbins)
    keep <- order(obs_cnt + all_cnt / R, decreasing = TRUE)[seq_len(top)]
  keep <- keep[(obs_cnt[keep] + all_cnt[keep]) > 0]
  if (!length(keep)) return(NULL)

  nb  <- length(keep)
  pos <- integer(nbins)
  pos[keep] <- seq_along(keep)

  b     <- pos[sim]
  repid <- rep(seq_len(R), each = nrow(sim))
  ok    <- !is.na(b) & b > 0L
  cnt   <- matrix(tabulate((repid[ok] - 1L) * nb + b[ok], nbins = nb * R),
                  nrow = nb, ncol = R)

  qs   <- t(apply(cnt, 1L, stats::quantile,
                  probs = c(alpha / 2, 0.5, 1 - alpha / 2), names = FALSE))
  mn   <- rowMeans(cnt)
  sdv  <- apply(cnt, 1L, stats::sd)
  p_ge <- rowSums(cnt >= obs_cnt[keep]) / R
  p_le <- rowSums(cnt <= obs_cnt[keep]) / R

  list(
    table = data.frame(
      id       = keep,
      label    = if (is.null(labels)) as.character(keep)
                 else unname(labels[as.character(keep)]),
      observed = obs_cnt[keep],
      sim_mean = mn,
      sim_lo   = qs[, 1L],
      sim_med  = qs[, 2L],
      sim_hi   = qs[, 3L],
      p_two    = pmin(1, 2 * pmin(p_ge, p_le)),
      z        = ifelse(sdv > 0, (obs_cnt[keep] - mn) / sdv, NA_real_),
      stringsAsFactors = FALSE
    ),
    # five-number summary per unit (5 x nb), drawn directly by bxp()
    bxp_stats = vapply(seq_len(nb),
                       function(i) grDevices::boxplot.stats(cnt[i, ])$stats,
                       numeric(5)),
    R     = R,
    alpha = alpha
  )
}

# -- tie model -----------------------------------------------------------------
# stats_3d / pars / baseline are the already-subset objects used for recall.
# The dyad space here is whatever .diagnostics_prepare() left in place: for a
# riskset of "active", "active_saturated" or "manual" that is the reduced
# (dyadIDactive) space, not N(N-1).
.gof_tie <- function(pars, baseline, stats_3d, obs_dyad_ids, reh,
                     R = 1000L, top_dyads = 30L, alpha = 0.05) {
  block <- .gof_sample_block(pars, baseline, stats_3d, obs_dyad_ids, R = R)
  if (is.null(block)) return(NULL)

  dm <- reh$index$dyad_map_active %||% reh$index$dyad_map
  if (is.null(dm)) return(NULL)
  id_col <- if ("dyadIDactive" %in% names(dm)) "dyadIDactive" else "dyadID"
  n_dyad <- max(dm[[id_col]])
  N      <- reh$N

  act_lab <- .actor_labels(reh)
  # dm$actor1 / dm$actor2 are actor *names*; map back to ids via the dictionary
  act_id <- stats::setNames(seq_len(N), act_lab[as.character(seq_len(N))])
  a1_of  <- integer(n_dyad); a2_of <- integer(n_dyad)
  a1_of[dm[[id_col]]] <- unname(act_id[as.character(dm$actor1)])
  a2_of[dm[[id_col]]] <- unname(act_id[as.character(dm$actor2)])

  to_role <- function(mat, lookup) {
    m <- matrix(lookup[mat], nrow = nrow(mat), ncol = ncol(mat))
    storage.mode(m) <- "integer"
    m
  }
  list(
    sender = .gof_tally(list(sim = to_role(block$sim, a1_of),
                             obs = a1_of[block$obs]),
                        nbins = N, labels = act_lab, alpha = alpha),
    receiver = .gof_tally(list(sim = to_role(block$sim, a2_of),
                               obs = a2_of[block$obs]),
                          nbins = N, labels = act_lab, alpha = alpha),
    dyad = .gof_tally(block, nbins = n_dyad, labels = .dyad_labels(reh),
                      alpha = alpha, top = top_dyads),
    R = R
  )
}

# -- actor model, sender submodel ----------------------------------------------
.gof_actor_sender <- function(pars, baseline, stats_3d, obs_ids, valid_ids,
                              reh, R = 1000L, alpha = 0.05) {
  block <- .gof_sample_block(pars, baseline, stats_3d, obs_ids, valid_ids, R = R)
  .gof_tally(block, nbins = reh$N, labels = .actor_labels(reh), alpha = alpha)
}

# -- actor model, receiver submodel --------------------------------------------
# receiver_stats at event m are built conditional on the sender that actually
# acted, so the receiver is drawn given the observed sender. The dyad table
# therefore pairs the observed sender with the simulated receiver.
.gof_actor_receiver <- function(pars, stats_3d, obs_ids, valid_ids, reh,
                                R = 1000L, top_dyads = 30L, alpha = 0.05) {
  block <- .gof_sample_block(pars, 0, stats_3d, obs_ids, valid_ids, R = R)
  if (is.null(block)) return(NULL)
  N       <- reh$N
  act_lab <- .actor_labels(reh)

  out <- list(receiver = .gof_tally(block, nbins = N, labels = act_lab,
                                    alpha = alpha))

  senders <- as.integer(unlist(attr(reh, "actor1ID")))
  if (length(senders) == nrow(block$sim)) {
    pair     <- function(s, r) (s - 1L) * N + r
    sim_pair <- matrix(pair(rep(senders, times = ncol(block$sim)),
                            as.vector(block$sim)),
                       nrow = nrow(block$sim), ncol = ncol(block$sim))
    pair_lab <- stats::setNames(
      paste(rep(act_lab, each = N), "->", rep(act_lab, times = N)),
      as.character(seq_len(N * N)))
    out$dyad <- .gof_tally(list(sim = sim_pair, obs = pair(senders, block$obs)),
                           nbins = N * N, labels = pair_lab,
                           alpha = alpha, top = top_dyads)
  }
  out
}

# ══════════════════════════════════════════════════════════════════════════
# plotting (plots 10-12 of plot.diagnostics)
# ══════════════════════════════════════════════════════════════════════════

.gof_panel <- function(blk, title, max_units = 30L) {
  if (is.null(blk) || !nrow(blk$table)) {
    plot.new()
    mtext(paste(title, "(no data)"), side = 3, line = 1, cex = 1.2)
    return(invisible(NULL))
  }
  tb  <- blk$table
  ord <- order(tb$observed, decreasing = TRUE)
  ord <- ord[seq_len(min(max_units, length(ord)))]

  labs <- tb$label[ord]
  labs[is.na(labs)] <- as.character(tb$id[ord][is.na(labs)])
  obs_ord <- tb$observed[ord]
  st      <- blk$bxp_stats[, ord, drop = FALSE]

  ylim <- range(0, st, obs_ord)
  ylim[2] <- ylim[2] * 1.18   # headroom for the legend

  old_par <- par(mar = c(7, 4, 4, 2) + 0.1)
  on.exit(par(old_par))

  bxp(list(stats = st,
           n     = rep(blk$R, length(ord)),
           conf  = matrix(NA_real_, nrow = 2L, ncol = length(ord)),
           out   = numeric(0),
           group = numeric(0),
           names = labs),
      outline = FALSE, las = 2, ylim = ylim, ylab = "Count", main = "",
      boxfill = grDevices::adjustcolor("lavender", alpha.f = 0.8),
      border  = "grey40")
  points(seq_along(ord), obs_ord, pch = 18, cex = 1.4, col = 2)
  legend("topright", bty = "n",
         legend = c("Simulated", "Observed"),
         pch    = c(15, 18),
         col    = c("lavender", 2),
         pt.cex = c(1.4, 1.4))
  mtext(title, side = 3, line = 1, cex = 1.2)
}

# ══════════════════════════════════════════════════════════════════════════
# print helper (one line inside print.diagnostics, under Recall)
# ══════════════════════════════════════════════════════════════════════════

.gof_print_line <- function(gof) {
  if (is.null(gof)) return(invisible(NULL))
  part <- function(blk, nm) {
    if (is.null(blk) || !nrow(blk$table)) return(NULL)
    tb  <- blk$table
    inb <- sum(!(tb$observed < tb$sim_lo | tb$observed > tb$sim_hi))
    sprintf("%s %d/%d", nm, inb, nrow(tb))
  }
  bits <- c(part(gof$sender, "sender"),
            part(gof$receiver, "receiver"),
            part(gof$dyad, "dyad"))
  if (!length(bits)) return(invisible(NULL))
  ref   <- gof$sender %||% gof$receiver %||% gof$dyad
  alpha <- ref$alpha %||% 0.05
  R     <- gof$R %||% ref$R
  cat(sprintf("  GOF        : within %g%% band -- %s  (%s replicates)\n",
              100 * (1 - alpha), paste(bits, collapse = "  |  "), format(R)))
  invisible(NULL)
}
