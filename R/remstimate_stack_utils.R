# shared helpers for external remstimate backends (GLMM, GLMNET, MIXREM)

# Filter a recall table down to the most poorly-predicted (surprising) events.
# threshold: keep rows with rel_rank <= threshold (low rel_rank = observed event
# ranked near the bottom of the risk set = surprising). Sorted most-surprising first.
.surprises_from_recall <- function(rc, threshold = 0.2) {
  if (is.null(rc)) return(NULL)
  pe <- rc$per_event
  out <- pe[pe$rel_rank <= threshold, , drop = FALSE]
  out[order(out$rel_rank), , drop = FALSE]
}

# Resolve durem eidx values to "actor1 -> actor2" dyad labels via reh$edgelist.
.durem_dyad_labels <- function(reh, eidx) {
  if (is.null(eidx) || length(eidx) == 0L) return(character(0))
  paste(reh$edgelist$actor1[eidx], "->", reh$edgelist$actor2[eidx])
}

# remstats::stack_stats() drops the statistic's dimname when a SINGLE statistic
# is stacked, returning a default 'V2'-style column. The formula built from
# stat_names (which is correct, taken from the array dimnames) then cannot find
# that column and coxme/glm(m) fail with 'invalid type (closure) for variable ...'.
# Re-attach the names to the auto-generated columns; strictly a no-op unless the
# non-id columns are all default 'V<number>' names and their count matches.
.stack_name_stats <- function(df, stat_names) {
  if (is.null(df) || !length(stat_names)) return(df)
  id_cols <- c("time_index", "time", "obs", "dyad", "actor", "actor1", "actor2",
               "type", "weight", "log_interevent", "samp_offset", "eidx",
               "process", ".surv_time")
  stat_cols <- setdiff(colnames(df), id_cols)
  if (length(stat_cols) == length(stat_names) &&
      length(stat_cols) > 0L && all(grepl("^V[0-9]+$", stat_cols))) {
    colnames(df)[match(stat_cols, colnames(df))] <- stat_names
  }
  df
}

# Detect statistics that have no variation across the stacked design (constant
# columns), warn by name, and drop them. A statistic can vary over the full
# riskset yet still be degenerate on the decisions that are actually stacked --
# e.g. reciprocity on a duration end process for strictly one-directional data.
# Such columns are not identifiable: MLE (clogit/glm) silently returns NA, while
# some GLMM/GLMNET backends error ('one of the covariates is a constant').
# Dropping them lets every stacked backend proceed the way MLE effectively does.
# 'baseline' (the all-ones intercept column) is excluded, as it is constant by
# construction. Returns the (possibly reduced) design and stat_names.
.drop_constant_stats <- function(df, stat_names, label = NULL) {
  if (is.null(df) || !length(stat_names))
    return(list(df = df, stat_names = stat_names))
  sn <- stat_names[tolower(stat_names) != "baseline"]
  sn <- intersect(sn, colnames(df))
  const <- sn[vapply(sn, function(v) {
    x <- df[[v]]
    length(unique(x[!is.na(x)])) <= 1L
  }, logical(1))]
  if (length(const)) {
    tag <- if (is.null(label)) "" else paste0(" [", label, "]")
    warning("Statistic(s) with no variation in the stacked design", tag, ": ",
            paste(const, collapse = ", "),
            ". Not identifiable (constant column); dropped from the fit.",
            call. = FALSE)
    df         <- df[, setdiff(colnames(df), const), drop = FALSE]
    stat_names <- setdiff(stat_names, const)
  }
  list(df = df, stat_names = stat_names)
}

.remstimate_make_stack <- function(reh, stats, add_actors = TRUE) {

  # ── Pre-stacked pass-through (idempotent) ──────────────────────────────────
  if (inherits(stats, "remstats_stacked")) {
    if (!isTRUE(stats$model == "actor")) {
      df <- stats$remstats_stack
      if (is.null(df$samp_offset))
        df$samp_offset <- if ("weight" %in% names(df)) log(df$weight) else 0
      return(list(df = df, stat_names = stats$stat_names,
                  ordinal = isTRUE(stats$ordinal), model = "tie"))
    }
    df_s <- stats$sender_stack; df_r <- stats$receiver_stack
    if (!is.null(df_s) && is.null(df_s$samp_offset)) df_s$samp_offset <- 0
    if (!is.null(df_r) && is.null(df_r$samp_offset)) df_r$samp_offset <- 0
    sn_s <- stats$sender_stat_names;   if (is.null(sn_s)) sn_s <- character(0)
    sn_r <- stats$receiver_stat_names; if (is.null(sn_r)) sn_r <- character(0)
    return(list(df = list(sender = df_s, receiver = df_r),
                stat_names = list(sender_model = sn_s, receiver_model = sn_r),
                ordinal = isTRUE(stats$ordinal), model = "actor"))
  }

  model     <- if (inherits(stats, "aomstats")) "actor" else "tie"
  stacking_stats <- remstats::stack_stats(stats, reh, add_actors = add_actors)

  if (model == "tie") {
    df <- stacking_stats$remstats_stack
    # sampling correction: offset = -log(π_d) per row
    # for uniform case-control sampling this upweights each control by (D-1)/K
    if (inherits(stats, "tomstats_sampled")) {
      sp_mat         <- attr(stats, "samp_prob")   # [M × samp_num], case col = 1
      df$samp_offset <- as.vector(t(-log(sp_mat))) # event-major flatten
    } else {
      df$samp_offset <- rep(0, nrow(df))
    }
    if (inherits(stats, "remstats_durem")) {
      stat_names <- stats$stacked$stat_names
    }else{
      stat_names <- dimnames(stats)[[3]]
    }
    df <- .stack_name_stats(df, stat_names)
    cc <- .drop_constant_stats(df, stat_names)
    df <- cc$df; stat_names <- cc$stat_names
    list(df = df, stat_names = stat_names, ordinal = stacking_stats$ordinal, model = model)
  } else {
    df_s <- stacking_stats$sender_stack
    df_r <- stacking_stats$receiver_stack
    if (!is.null(df_s)) df_s$samp_offset <- rep(0, nrow(df_s))
    if (!is.null(df_r)) df_r$samp_offset <- rep(0, nrow(df_r))
    sn_s <- if (!is.null(df_s)) dimnames(stats$sender_stats)[[3]]   else character(0)
    sn_r <- if (!is.null(df_r)) dimnames(stats$receiver_stats)[[3]] else character(0)
    df_s <- .stack_name_stats(df_s, sn_s)
    df_r <- .stack_name_stats(df_r, sn_r)
    cs <- .drop_constant_stats(df_s, sn_s, "sender");   df_s <- cs$df; sn_s <- cs$stat_names
    cr <- .drop_constant_stats(df_r, sn_r, "receiver"); df_r <- cr$df; sn_r <- cr$stat_names
    list(
      df = list(sender = df_s, receiver = df_r),
      stat_names = list(sender_model = sn_s, receiver_model = sn_r),
      ordinal = stacking_stats$ordinal,
      model   = model
    )
  }
}

# rhs string for glmer / glmmTMB / flexmix formula
# interval timing: Poisson with log_interevent + samp_offset
# ordinal timing:  binomial, no time offset but sampling correction still applies
.remstimate_fixed_rhs <- function(stat_names, ordinal) {
  rhs <- paste(c("-1", stat_names), collapse = " + ")
  if (!ordinal)
    paste(rhs, "+ offset(log_interevent + samp_offset)")
  else
    paste(rhs, "+ offset(samp_offset)")
}

# model matrix for glmnet
.remstimate_model_matrix <- function(df, stat_names, ordinal) {
  off <- df$samp_offset
  if (!ordinal && !is.null(df$log_interevent))
    off <- off + df$log_interevent
  list(
    X       = as.matrix(df[, stat_names, drop = FALSE]),
    y       = df$obs,
    offset  = off,
    weights = if ("weight" %in% names(df)) df$weight else NULL
  )
}

.remstimate_find_baseline <- function(stat_names) {
  idx <- which(tolower(stat_names) == "baseline")
  if (length(idx)) idx[[1L]] else NULL
}

# Statistics that form the intercept / baseline structure and are left
# unpenalised by default in the penalised backends (GLMNET, SHRINKEM). The rule
# is structural, not name-based: a statistic is intercept-like when it is an
# indicator (all non-NA values in {0,1}). This captures the overall 'baseline'
# (all ones), the duration start/end process intercepts ('baseline.start' /
# 'baseline.end') and fixed-effect type dummies (e.g. 'FEtype_Y'), while genuine
# effects - counts (inertia) and standardised / continuous statistics - are
# penalised. 'data' may be a data.frame (stacked design) or a numeric matrix
# (model matrix); only the columns named in 'stat_names' are inspected.
.intercept_like_stats <- function(data, stat_names) {
  present <- intersect(stat_names, colnames(data))
  if (!length(present)) return(character(0))
  getcol <- if (is.data.frame(data)) function(nm) data[[nm]] else function(nm) data[, nm]
  present[vapply(present, function(nm) {
    u <- unique(getcol(nm)); u <- u[!is.na(u)]
    length(u) > 0L && all(u %in% c(0, 1))
  }, logical(1))]
}

# Resolve the final set of unpenalised statistic names for a penalised backend.
# Start from the auto-detected intercept structure (.intercept_like_stats), then
# apply the user's two additive controls from the 'penalty' list:
#   * unpenalized : names ADDED to the exemption (e.g. a continuous effect you
#                   want to keep out of the penalty)
#   * penalized   : names REMOVED from the exemption, i.e. forced BACK into the
#                   penalty even though they are 0/1 indicators (e.g. a p-shift
#                   dummy such as 'psABAB.end', or a fixed-effect dummy you do
#                   want to shrink)
# 'penalized' wins over 'unpenalized' when a name appears in both. Only names
# actually present in 'stat_names' are returned.
# Statistic names must be given exactly as they appear in the remstats object
# (e.g. 'psABAB.end', not 'psABAB'). Matching is exact.
.resolve_unpenalized <- function(data, stat_names,
                                  unpenalized = NULL, penalized = NULL) {
  default <- .intercept_like_stats(data, stat_names)
  final   <- setdiff(union(default, unpenalized), penalized)
  intersect(final, stat_names)
}

# Warn (once, at the backend entry) about names passed to penalty's 'unpenalized'
# / 'penalized' that do not match any model statistic - almost always a typo or a
# missing duration / type suffix. Such names are silently ignored otherwise.
.check_penalty_names <- function(unpenalized = NULL, penalized = NULL, valid = character(0)) {
  bad <- setdiff(unique(c(unpenalized, penalized)), valid)
  if (length(bad))
    warning("penalty name(s) not found among the model statistics and ignored: ",
            paste(bad, collapse = ", "), ".\n  Available statistics: ",
            paste(valid, collapse = ", "), ".", call. = FALSE)
  invisible(NULL)
}

# Duration-model diagnostics for the point-estimate penalised backends (GLMNET,
# SHRINKEM). A duration fit is stacked as a tie-like design but needs the
# start/end recall structure rather than the tie Schoenfeld / waiting-time
# diagnostics, and diagnostics.remstimate()'s .diagnostics_prepare() rejects a
# remstats_durem object outright. So we build the durem recall tables directly
# from the (penalised / shrunken) coefficients, mirroring the MLE durem path and
# the fixed-effects fallback of the GLMM durem diagnostics. Returns a
# 'diagnostics_durem' object so plot()/print() dispatch as for MLE durem.
.point_durem_diagnostics <- function(object, reh, stats, top_pct = 0.05,
                                     surprise_threshold = 0.2) {
  df <- object$stacked_data
  # SHRINKEM does not store the stacked design; rebuild it (identical to the one
  # the MLE used). GLMNET stores it as a data.frame and is reused as-is.
  if (is.null(df) || !is.data.frame(df))
    df <- .remstimate_make_stack(reh, stats, add_actors = FALSE)$df
  if (is.null(df) || !is.data.frame(df))
    stop("Duration diagnostics are only available for tie-oriented duration models.",
         call. = FALSE)

  sn    <- attr(object, "statistics")
  coefs <- object$coefficients
  sn    <- intersect(sn, colnames(df))
  X     <- as.matrix(df[, sn, drop = FALSE])
  lp    <- as.numeric(X %*% coefs[sn])

  out <- .durem_recall_tables(
    lp, df, reh,
    sn_start = stats$stacked$stat_names_start,
    sn_end   = stats$stacked$stat_names_end,
    top_pct  = top_pct,
    surprise_threshold = surprise_threshold)
  out$.reh.processed <- reh
  class(out) <- c("diagnostics_durem", "diagnostics", "remstimate")
  out
}

# build the return object; actor model mirrors MLE structure so
# diagnostics.remstimate works via inheritance without changes
.remstimate_wrap <- function(coefficients,
                              stat_names,
                              loglik       = NULL,
                              formula      = NULL,
                              stacked_data = NULL,
                              backend_fit  = NULL,
                              model        = "tie",
                              method       = "GLMM",
                              engine       = NA_character_,
                              ordinal      = FALSE,
                              extra        = list()) {

  if (model == "actor") {
    kern <- list(
      sender_model = list(
        coefficients = coefficients$sender_model,
        backend_fit  = if (is.list(backend_fit)) backend_fit$sender_model else backend_fit
      ),
      receiver_model = list(
        coefficients = coefficients$receiver_model,
        backend_fit  = if (is.list(backend_fit)) backend_fit$receiver_model else NULL
      ),
      loglik       = loglik,
      formula      = formula,
      stacked_data = stacked_data
    )
    waar_is_baseline <- .remstimate_find_baseline(stat_names$sender_model)
  } else {
    kern <- list(
      coefficients = coefficients,
      loglik       = loglik,
      formula      = formula,
      stacked_data = stacked_data,
      backend_fit  = backend_fit
    )
    waar_is_baseline <- .remstimate_find_baseline(stat_names)
  }

  structure(
    c(kern, extra),
    class             = c(paste0("remstimate_", tolower(method)), "remstimate"),
    model             = model,
    approach          = "Frequentist",
    method            = method,
    engine            = engine,
    ordinal           = ordinal,
    statistics        = stat_names,
    where_is_baseline = waar_is_baseline,
    ncores            = 1L
  )
}

# extract grouping vars from e.g. ~ (1|actor1) + (1|actor2)
.parse_random_groups <- function(random) {
  if (is.null(random)) return(character(0))
  rhs  <- deparse(random[[2]])
  hits <- gregexpr("\\|\\s*([A-Za-z_.][A-Za-z_.0-9]*)", rhs, perl = TRUE)
  trimws(sub("^\\|\\s*", "", regmatches(rhs, hits)[[1]]))
}


# ── Shared recall computation ────────────────────────────────────────────────
# Used by diagnostics methods for GLMM, GLMNET, MIXREM, and (optionally) durem.

# Fair ranking of observed outcomes within a risk set. Given the normalised
# within-event probabilities and the position(s) of the observed outcome(s),
# returns the average rank (1 = most likely; tied outcomes share the mean of
# their ranks) and the cumulative probability of outcomes at least as likely as
# the observed. Average ranks are essential: without them a model that cannot
# discriminate (all probabilities equal) is silently handed rank 1 by the
# arbitrary storage order of the risk set, inflating recall well above chance.
# With average ranks such a model scores ~0.5 (chance), as it should.
.recall_ranks <- function(probs, pos) {
  r_all <- rank(-probs, ties.method = "average")   # 1 = highest probability
  p_obs <- probs[pos]
  list(rank = r_all[pos],
       cum  = vapply(p_obs, function(p) sum(probs[probs > p]) + p, numeric(1)))
}

# Rescale an average rank (1 = top of a D_t-sized risk set, D_t = bottom) onto
# [0, 1]. Dividing the beaten-competitor count (0 .. D_t-1) by D_t-1 rather than
# D_t makes the best possible outcome score exactly 1 (instead of 1 - 1/D_t) and
# the worst 0; a non-discriminating model still averages ~0.5. A singleton risk
# set (D_t == 1) has no competitors to rank against, so it carries no diagnostic
# signal and is scored at the neutral 0.5 (avoids a 0/0 and never spuriously
# flags a forced, information-free prediction as a surprise). rank may be a
# vector (simultaneous observed events); D_t is the scalar risk-set size.
.rel_rank <- function(rank, D_t) {
  if (D_t > 1L) 1 - (rank - 1) / (D_t - 1) else rep(0.5, length(rank))
}

# Core: rank observed events within time-point groups
.recall_block <- function(lp, obs_idx, event_ids, top_pct = 0.05, ids = NULL, min_D_t = 1) {
  events <- unique(event_ids)
  rows <- vector("list", length(events))
  for (i in seq_along(events)) {
    ev   <- events[i]
    mask <- which(event_ids == ev)
    obs  <- intersect(obs_idx, mask)
    if (length(obs) == 0L) next
    D_t  <- length(mask)
    if (D_t < min_D_t) next
    obs_pos <- match(obs, mask)
    lp_ev   <- lp[mask]
    probs   <- exp(lp_ev)
    probs   <- probs / sum(probs)
    rr        <- .recall_ranks(probs, obs_pos)
    obs_probs <- probs[obs_pos]
    rows[[i]] <- data.frame(
      time       = ev,
      rel_rank   = .rel_rank(rr$rank, D_t),
      cum_prob   = rr$cum,
      obs_prob   = obs_probs,
      prob_ratio = obs_probs * D_t,
      log_loss   = -log(obs_probs),
      D_t        = D_t,
      eidx       = if (!is.null(ids)) ids[obs] else NA_integer_
    )
  }
  pe <- do.call(rbind, rows)
  if (is.null(pe) || nrow(pe) == 0L) return(NULL)
  list(
    per_event = pe,
    summary   = data.frame(
      mean_rel_rank   = mean(pe$rel_rank),
      median_rel_rank = stats::median(pe$rel_rank),
      mean_cum_prob   = mean(pe$cum_prob),
      mean_prob_ratio = mean(pe$prob_ratio),
      mean_log_loss   = mean(pe$log_loss),
      top_pct         = top_pct,
      top_pct_prop    = mean(pe$rel_rank >= 1 - top_pct)
    )
  )
}

# Full recall: joint + per-type (if type column exists)
.diagnostics_recall <- function(lp, df, top_pct = 0.05) {
  obs_idx   <- which(df$obs == 1L)
  event_ids <- df$time
  out <- list()

  out$recall <- .recall_block(lp, obs_idx, event_ids, top_pct)

  if ("type" %in% names(df) && !all(is.na(df$type))) {
    types <- sort(unique(df$type[!is.na(df$type)]))
    if (length(types) > 1L) {
      out$recall_by_type <- list()
      for (tp in types) {
        tp_mask <- which(df$type == tp)
        tp_obs  <- intersect(obs_idx, tp_mask)
        if (length(tp_obs) == 0L) next
        out$recall_by_type[[tp]] <- .recall_block(
          lp[tp_mask], match(tp_obs, tp_mask),
          event_ids[tp_mask], top_pct
        )
      }
    }
  }

  out
}

# Print helper for recall output
.print_recall_summary <- function(rc, label, surprises = NULL, threshold = NULL) {
  if (is.null(rc)) return()
  rs  <- rc$summary
  pct <- round(rs$top_pct_prop * 100, 1)
  low_str <- ""
  if (!is.null(threshold)) {
    n_pe    <- nrow(rc$per_event)
    n_sur   <- if (is.null(surprises)) 0L else nrow(surprises)
    low_pct <- if (n_pe > 0) round(100 * n_sur / n_pe, 1) else NA
    low_str <- sprintf(" | lowest %g%% = %s%%", threshold * 100, low_pct)
  }
  cat(sprintf("  %-10s: mean rank = %.3f | prob ratio = %.2f | top %g%% = %s%%%s\n",
              label, rs$mean_rel_rank, rs$mean_prob_ratio,
              rs$top_pct * 100, pct, low_str))
}

# Plot helper for recall output (relative rank). Self-contained margins (so it
# does not inherit a crowded top margin from a previous plot) plus a median line
# and a smoothed trend over time, matching the tie/actor recall plot style.
.plot_recall_scatter <- function(rc, main, ...) {
  if (is.null(rc)) return()
  pe <- rc$per_event
  if (is.null(pe) || nrow(pe) == 0L) return()

  op <- graphics::par(mar = c(4.5, 4.5, 4, 1.5))
  on.exit(graphics::par(op), add = TRUE)

  graphics::plot(pe$time, pe$rel_rank,
                 xlab = "Time point", ylab = "Relative rank (0 = bottom, 1 = top)",
                 main = "", ylim = c(0, 1),
                 pch = 16, cex = 0.5, col = grDevices::adjustcolor("black", 0.3), ...)

  med <- stats::median(pe$rel_rank, na.rm = TRUE)
  graphics::abline(h = med, lty = 2, col = "blue", lwd = 1.5)

  # smoothed trend of the relative rank over time (runmed then a spline), to
  # reveal stretches where recall was systematically better/worse
  ord       <- order(pe$time)
  x_ord     <- pe$time[ord]
  y_ord     <- pe$rel_rank[ord]
  has_trend <- FALSE
  if (length(y_ord) >= 4L) {
    k <- max(3L, min(length(y_ord) - 1L, round(length(y_ord) * 0.1)))
    if (k %% 2 == 0) k <- k + 1L   # runmed needs an odd bandwidth
    y_med <- tryCatch(stats::runmed(y_ord, k = k), error = function(e) NULL)
    if (!is.null(y_med)) {
      sm <- tryCatch(suppressWarnings(stats::smooth.spline(x_ord, y_med)),
                     error = function(e) NULL)
      if (!is.null(sm)) graphics::lines(sm$x, sm$y, col = "firebrick", lwd = 3.5)
      else              graphics::lines(x_ord, y_med, col = "firebrick", lwd = 3.5)
      has_trend <- TRUE
    }
  }

  graphics::legend("bottomright", inset = 0.02, bg = "white", cex = 0.8,
                   legend = c("Rank", "Median rank", if (has_trend) "Smoothed trend"),
                   lty    = c(NA, 2, if (has_trend) 1),
                   pch    = c(16, NA, if (has_trend) NA),
                   col    = c(grDevices::adjustcolor("black", 0.3), "blue",
                              if (has_trend) "firebrick"),
                   lwd    = c(NA, 1.5, if (has_trend) 3.5))

  graphics::title(main = main, line = 1.6, cex.main = 1.1)
}

# Plot helper for probability ratio
.plot_probratio_scatter <- function(rc, main, ...) {
  if (is.null(rc)) return()
  pe <- rc$per_event
  if (is.null(pe) || nrow(pe) == 0L) return()

  op <- graphics::par(mar = c(4.5, 4.5, 4, 1.5))
  on.exit(graphics::par(op), add = TRUE)

  graphics::plot(pe$time, pe$prob_ratio,
                 xlab = "Time point", ylab = "Probability ratio",
                 main = "",
                 pch = 16, cex = 0.5, col = grDevices::adjustcolor("black", 0.3), ...)
  graphics::abline(h = 1, col = "grey50", lty = 3)
  graphics::abline(h = mean(pe$prob_ratio), col = "red", lty = 2)
  graphics::title(main = main, line = 1.6, cex.main = 1.1)
}


# Tabulate how often each observed id (dyad or actor) shows up among surprises
# vs. among all observed events, so repeat offenders surface immediately.
.offender_table <- function(ids_surprises, ids_all, labels = NULL) {
  if (length(ids_all) == 0L) return(NULL)
  tab_all <- table(ids_all)
  tab_sur <- table(ids_surprises)
  ids <- names(tab_all)
  n_total     <- as.integer(tab_all[ids])
  n_surprises <- as.integer(tab_sur[ids])
  n_surprises[is.na(n_surprises)] <- 0L
  out <- data.frame(
    id          = ids,
    n_surprises = n_surprises,
    n_total     = n_total,
    prop        = n_surprises / n_total,
    stringsAsFactors = FALSE
  )
  if (!is.null(labels)) out$label <- unname(labels[out$id])
  out[order(-out$n_surprises, -out$prop), ]
}

.dyad_labels <- function(reh) {
  dm <- reh$index$dyad_map_active %||% reh$index$dyad_map
  if (is.null(dm)) return(NULL)
  id_col <- if ("dyadIDactive" %in% names(dm)) "dyadIDactive" else "dyadID"
  setNames(paste(dm$actor1, "->", dm$actor2), dm[[id_col]])
}

.actor_labels <- function(reh) {
  ad <- reh$meta$dictionary$actors
  if (is.null(ad)) return(NULL)
  setNames(ad$actorName, ad$actorID)
}

# Look up id_list[[event]][sub_index] per row — recovers the sender for a
# receiver-model surprise row (receiver choice is conditional on sender).
.lookup_by_event <- function(event, sub_index, id_list) {
  mapply(function(ev, si) {
    v <- id_list[[ev]]
    if (si > length(v)) return(NA_integer_)
    as.integer(v[si])
  }, event, sub_index)
}

# Build "SenderName -> ReceiverName" labels from raw actor IDs.
.actor_pair_labels <- function(reh, sender_ids, receiver_ids) {
  ad <- reh$meta$dictionary$actors
  if (is.null(ad)) return(NULL)
  nm <- setNames(ad$actorName, as.character(ad$actorID))
  paste(unname(nm[as.character(sender_ids)]), "->", unname(nm[as.character(receiver_ids)]))
}
