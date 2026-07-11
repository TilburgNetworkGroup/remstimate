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
    list(df = df, stat_names = stat_names, ordinal = stacking_stats$ordinal, model = model)
  } else {
    df_s <- stacking_stats$sender_stack
    df_r <- stacking_stats$receiver_stack
    if (!is.null(df_s)) df_s$samp_offset <- rep(0, nrow(df_s))
    if (!is.null(df_r)) df_r$samp_offset <- rep(0, nrow(df_r))
    list(
      df = list(sender = df_s, receiver = df_r),
      stat_names = list(
        sender_model   = if (!is.null(df_s)) dimnames(stats$sender_stats)[[3]] else character(0),
        receiver_model = if (!is.null(df_r)) dimnames(stats$receiver_stats)[[3]] else character(0)
      ),
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
    ord     <- order(probs, decreasing = TRUE)
    rnks    <- match(obs_pos, ord)
    rnks    <- rnks[!is.na(rnks)]
    if (length(rnks) == 0L) next
    obs_probs <- probs[obs_pos]
    rows[[i]] <- data.frame(
      time       = ev,
      rel_rank   = 1 - rnks / D_t,
      cum_prob   = cumsum(probs[ord])[rnks],
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

# Plot helper for recall output (relative rank)
.plot_recall_scatter <- function(rc, main, ...) {
  if (is.null(rc)) return()
  pe <- rc$per_event
  graphics::plot(pe$time, pe$rel_rank,
                 xlab = "Time point", ylab = "Relative rank",
                 main = main,
                 pch = ".", col = grDevices::adjustcolor("black", 0.3), ...)
  graphics::abline(h = mean(pe$rel_rank), col = "red", lty = 2)
}

# Plot helper for probability ratio
.plot_probratio_scatter <- function(rc, main, ...) {
  if (is.null(rc)) return()
  pe <- rc$per_event
  graphics::plot(pe$time, pe$prob_ratio,
                 xlab = "Time point", ylab = "Probability ratio",
                 main = main,
                 pch = ".", col = grDevices::adjustcolor("black", 0.3), ...)
  graphics::abline(h = 1, col = "grey50", lty = 3)
  graphics::abline(h = mean(pe$prob_ratio), col = "red", lty = 2)
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
