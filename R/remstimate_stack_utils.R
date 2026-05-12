# shared helpers for external remstimate backends (GLMM, GLMNET, MIXREM)

.remstimate_make_stack <- function(reh, stats, add_actors = TRUE) {

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
      stat_names <- c(dimnames(stats$start_stats)[[3]],dimnames(stats$end_stats)[[3]])
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
  list(
    X       = as.matrix(df[, stat_names, drop = FALSE]),
    y       = df$obs,
    offset  = if (!ordinal) df$log_interevent + df$samp_offset else df$samp_offset,
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
.recall_block <- function(lp, obs_idx, event_ids, top_pct = 0.05) {
  events <- unique(event_ids)
  rows <- vector("list", length(events))
  for (i in seq_along(events)) {
    ev   <- events[i]
    mask <- which(event_ids == ev)
    obs  <- intersect(obs_idx, mask)
    if (length(obs) == 0L) next
    obs_pos <- match(obs, mask)
    lp_ev   <- lp[mask]
    probs   <- exp(lp_ev)
    probs   <- probs / sum(probs)
    ord     <- order(probs, decreasing = TRUE)
    rnks    <- match(obs_pos, ord)
    rnks    <- rnks[!is.na(rnks)]
    if (length(rnks) == 0L) next
    rows[[i]] <- data.frame(
      time     = ev,
      rel_rank = 1 - rnks / length(mask),
      cum_prob = cumsum(probs[ord])[rnks]
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
.print_recall_summary <- function(rc, label) {
  if (is.null(rc)) return()
  rs  <- rc$summary
  pct <- round(rs$top_pct_prop * 100, 1)
  cat(sprintf("  %-10s: mean rank = %.3f | top %g%% = %s%%\n",
              label, rs$mean_rel_rank, rs$top_pct * 100, pct))
}

# Plot helper for recall output
.plot_recall_scatter <- function(rc, main, ...) {
  if (is.null(rc)) return()
  pe <- rc$per_event
  graphics::plot(pe$time, pe$rel_rank,
                 xlab = "Time point", ylab = "Relative rank",
                 main = main,
                 pch = ".", col = grDevices::adjustcolor("black", 0.3), ...)
  graphics::abline(h = mean(pe$rel_rank), col = "red", lty = 2)
}
