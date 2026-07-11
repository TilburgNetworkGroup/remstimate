# -- recall helper -------------------------------------------------------------
.recall_block_3d <- function(pars, baseline, stats_3d, obs_ids,
                             valid_ids = NULL, top_pct = 0.05) {
  M    <- dim(stats_3d)[3L]
  rows <- vector("list", M)
  for (m in seq_len(M)) {
    obs   <- obs_ids[[m]]
    if (!length(obs)) next
    valid <- if (is.null(valid_ids)) seq_len(dim(stats_3d)[1L]) else valid_ids[[m]]
    if (!length(valid)) next
    S     <- matrix(stats_3d[valid, , m], nrow = length(valid))
    probs <- exp(as.numeric(baseline + S %*% pars))
    probs <- probs / sum(probs)
    D_t   <- length(valid)
    ord   <- order(probs, decreasing = TRUE)
    pos_all   <- match(obs, valid)
    keep      <- which(!is.na(pos_all))
    if (!length(keep)) next
    pos       <- pos_all[keep]
    sub_index <- keep
    rnks      <- match(pos, ord)
    keep2     <- which(!is.na(rnks))
    if (!length(keep2)) next
    rnks      <- rnks[keep2]
    pos       <- pos[keep2]
    sub_index <- sub_index[keep2]
    obs_probs <- probs[pos]
    rows[[m]] <- data.frame(
      event      = m,
      sub_index  = sub_index,
      obs_id     = obs[sub_index],
      rel_rank   = 1 - rnks / D_t,
      cum_prob   = cumsum(probs[ord])[rnks],
      obs_prob   = obs_probs,
      prob_ratio = obs_probs * D_t,
      log_loss   = -log(obs_probs)
    )
  }
  pe <- do.call(rbind, rows)
  list(
    per_event = pe,
    summary   = data.frame(
      mean_rel_rank   = mean(pe$rel_rank),
      median_rel_rank = median(pe$rel_rank),
      mean_cum_prob   = mean(pe$cum_prob),
      mean_prob_ratio = mean(pe$prob_ratio),
      mean_log_loss   = mean(pe$log_loss),
      top_pct         = top_pct,
      top_pct_prop    = mean(pe$rel_rank >= 1 - top_pct)
    )
  )
}

# Actor-level offender table for the tie model: how often does an actor
# appear among surprises as actor1 (sender role) / actor2 (receiver role) /
# either, relative to how often it appears in that role among all observed
# events. dyad_map must have columns actor1, actor2 aligned to obs_id (dyad ID).
.tie_actor_offenders <- function(surprises, per_event, dyad_map, id_col) {
  if (is.null(surprises) || nrow(surprises) == 0L) return(NULL)

  a1_sur <- dyad_map$actor1[match(surprises$obs_id, dyad_map[[id_col]])]
  a2_sur <- dyad_map$actor2[match(surprises$obs_id, dyad_map[[id_col]])]
  a1_all <- dyad_map$actor1[match(per_event$obs_id, dyad_map[[id_col]])]
  a2_all <- dyad_map$actor2[match(per_event$obs_id, dyad_map[[id_col]])]

  list(
    sender   = .offender_table(a1_sur, a1_all),
    receiver = .offender_table(a2_sur, a2_all),
    either   = .offender_table(c(a1_sur, a2_sur), c(a1_all, a2_all))
  )
}

# diagnostics
#' @title Compute the diagnostics of a \code{remstimate} object. Diagnostics are based on point estimates, also for Bayesian fit.
#' @param object is a \code{remstimate} object.
#' @param reh is a \code{remify} object, the same used for the 'remstimate' object.
#' @param stats is a \code{remstats} object, the same used for the 'remstimate' object.
#' @param top_pct numeric scalar in (0,1): threshold for the top-percentage recall summary (default 0.05).
#' @param surprise_threshold numeric scalar in (0,1): rows of the recall table with
#'   \code{rel_rank <= surprise_threshold} are collected into \code{$surprises}
#'   (default 0.20). Lower \code{rel_rank} means the observed outcome was ranked
#'   further from the top of the predicted probabilities, i.e. more surprising.
#' @param ... further arguments to be passed to the 'diagnostics' method.
#' @details
#' Every recall table (\code{$recall$per_event} for tie/actor,
#' \code{$recall_joint}/\code{$recall_start}/\code{$recall_end} for durem) and every
#' \code{$surprises} table derived from it shares this column layout:
#' \describe{
#'   \item{event}{Row position within the analyzed subset (tie/actor only).}
#'   \item{edgelist_id_row}{The corresponding row index into \code{reh$edgelist_id}
#'     (\code{event} resolved to the full, untrimmed object), computed as
#'     \code{start_stop[1] - 1 + event}.
#'   \item{sub_index}{Which simultaneous observation at that \code{event} this row
#'     is, when multiple events share a time point.}
#'   \item{obs_id}{The observed outcome as an ID. Tie model: the observed dyad ID.
#'     Actor sender submodel: the observed sender actor ID. Actor receiver
#'     submodel: the observed receiver actor ID (see \code{sender_id} below for
#'     the paired sender).}
#'   \item{sender_id, obs_label, sender_label, dyad_label}{Present on
#'     \code{$surprises} only. Human-readable resolutions of \code{obs_id}
#'     (and, for the actor receiver submodel, the paired sender), added for
#'     readability -- not present on the raw \code{$recall} tables.}
#'   \item{rel_rank}{\code{1 - rank/D_t}: 1 = observed outcome ranked most likely
#'     of all risk-set alternatives; 0 = ranked last (maximally surprising).}
#'   \item{cum_prob}{Cumulative predicted probability of everything ranked
#'     at-or-above the observed outcome.}
#'   \item{obs_prob}{Predicted probability assigned to the observed outcome.}
#'   \item{prob_ratio}{\code{obs_prob * D_t}: ratio vs. uniform/random guessing.
#'     Near 1 means the model had no real signal for this outcome (e.g. a
#'     cold-start actor/dyad); well below 1 means the model actively favored a
#'     different outcome (a genuine misspecification signal).}
#'   \item{log_loss}{\code{-log(obs_prob)}, the per-observation negative
#'     log-likelihood contribution.}
#'   \item{D_t}{Size of the risk set (number of eligible alternatives) at that
#'     decision. For durem end-process recall, epochs with \code{D_t = 1} are
#'     excluded by default, since a single-candidate softmax is always exactly
#'     1 regardless of model fit and carries no diagnostic information.}
#' }
#' }
#' \code{$surprise_offenders} (tie, actor sender/receiver) and
#' \code{$surprise_offenders_dyad} (actor receiver only) tabulate, per dyad or
#' actor, how often it appears among surprises relative to how often it was
#' actually observed (\code{prop = n_surprises / n_total}), sorted worst first.
#'
#' For durem, \code{obs_id}/\code{event}/\code{sub_index} are replaced by
#' \code{eidx}, a row index directly into \code{reh$edgelist} (i.e. already
#' resolved to the original input edgelist row -- no offset arithmetic needed).
#' \code{$surprise_offenders_joint}/\code{_start}/\code{_end} tabulate, per
#' dyad (\code{"actor1 -> actor2"}), how often it appears among surprises
#' relative to how often it was actually observed in that recall table.
#' @export
diagnostics <- function(object, reh, stats, ...) {
  UseMethod("diagnostics")
}

denormalize_reh <- function(reh) {
  if (is.null(reh$meta)) return(reh)
  attr(reh, "model")    <- reh$meta$model
  attr(reh, "directed") <- reh$meta$directed
  attr(reh, "ordinal")  <- reh$meta$ordinal
  attr(reh, "weighted") <- reh$meta$weighted
  attr(reh, "riskset")  <- reh$meta$riskset_source %||% reh$meta$riskset
  attr(reh, "origin")   <- reh$meta$origin
  attr(reh, "dyadID")        <- reh$ids$dyad
  attr(reh, "dyadIDactive")  <- reh$ids$dyad_active
  attr(reh, "actor1ID")      <- reh$ids$actor1
  attr(reh, "actor2ID")      <- reh$ids$actor2
  attr(reh, "indices_simultaneous_events") <- reh$indices_simultaneous_events
  reh
}

# ══════════════════════════════════════════════════════════════════════════
# .diagnostics_prepare() — shared setup extracted from diagnostics.remstimate()
# ══════════════════════════════════════════════════════════════════════════

.diagnostics_prepare <- function(object, reh, stats) {
  reh <- denormalize_reh(reh)
  if (!inherits(reh, "remify")) {
    stop("'reh' must be a 'remify' object (see ?remify::remify).")
  }
  model <- ifelse(is.null(attr(reh, "model")), "", attr(reh, "model"))
  if (!(model %in% c("actor", "tie"))) {
    stop("attribute 'model' of input 'reh' must be either 'actor' or 'tie'")
  }
  model_object <- ifelse(is.null(attr(object, "model")), "", attr(object, "model"))
  model_stats  <- ifelse(is.null(attr(stats,  "model")), "", attr(stats,  "model"))
  if ((model_stats != model) | (model_object != model)) {
    stop("attribute 'model' of input 'object' and 'stats' must be the same as the attribute of input 'reh'")
  }
  # active / manual (reduced) riskset. attr(reh,"riskset") is riskset_source,
  # so it can be "active", "active_saturated" or "manual" - all use the reduced
  # dyad indexing (activeD / dyadIDactive). Matching only "active" here left the
  # saturated/manual cases with a full-riskset dyadID and an uncleared omit_dyad
  # (which lacks a 'riskset' element), causing computeDiagnostics to fail with
  # 'Index out of bounds: [index=riskset]'.
  if (grepl("^active", attr(reh, "riskset")) || identical(attr(reh, "riskset"), "manual")) {
    reh$D <- reh$activeD
    if (model == "tie") {
      attr(reh, "dyadID") <- attr(reh, "dyadIDactive")
      reh$omit_dyad <- list()
    }
  }
  ordinal <- attr(reh, "ordinal")
  ncores  <- NULL   # [new] defensive default; only ever reached if attr(object,"ncores") is NULL,
  #       which does not occur for any remstimate object produced by this package
  if (!is.null(attr(object, "ncores"))) ncores <- attr(object, "ncores")
  if (is.null(reh$omit_dyad)) reh$omit_dyad <- list()

  # start/stop and method
  if (all(inherits(stats, c("remstats","tomstats"), TRUE)) |
      all(inherits(stats, c("remstats","aomstats"), TRUE))) {
    stats_attr_method <- attr(stats, "method")
    if (is.null(stats_attr_method)) {
      stop("attribute 'method' not found inside object 'remstats'. Input argument 'stats' must be an object of class 'remstats' from the package 'remstats' (>=3.2.0)")
    }
    omit_dyad_receiver <- NULL
    if (stats_attr_method == "pe") {
      if (!is.null(attr(reh, "evenly_spaced_interevent_time")))
        reh$intereventTime <- attr(reh, "evenly_spaced_interevent_time")
      if (!is.null(reh$E)) reh$M <- reh$E
    }
    start_stop <- if (!is.null(attr(stats, "subset")))
      as.numeric(unlist(attr(stats, "subset"))) else c(1, reh$M)
  } else {
    stop("'stats' must be a 'remstats' object from the package 'remstats' (>= 3.2.0), suitable for tie-oriented modeling ('tomstats') or actor-oriented modeling ('aomstats')")
  }

  model_formula <- variables_names <- where_is_baseline <- NULL

  # -- tie-oriented -------------------------------------------------------------
  if (model == "tie") {
    if (all(inherits(stats, c("remstats","tomstats"), TRUE))) {
      if (!is.null(dimnames(stats)[[3]])) variables_names <- dimnames(stats)[[3]]
      model_formula <- if (is.null(attr(stats, "formula")))
        stats::as.formula(paste("~ ", paste(variables_names, collapse=" + ")))
      else attr(stats, "formula")
      if (any(tolower(variables_names) %in% c("baseline")))
        where_is_baseline <- which(variables_names == "baseline")
      stats <- aperm(stats, perm = c(2,3,1))
      if (stats_attr_method == "pt") {
        attr(reh,"dyadID")  <- attr(reh,"dyadID")[start_stop[1]:start_stop[2]]
        attr(reh,"actor1ID") <- attr(reh,"actor1ID")[start_stop[1]:start_stop[2]]
        attr(reh,"actor2ID") <- attr(reh,"actor2ID")[start_stop[1]:start_stop[2]]
        if ((length(reh$omit_dyad) > 0) & !is.null(attr(reh,"indices_simultaneous_events")))
          reh$omit_dyad$time <- reh$omit_dyad$time[-attr(reh,"indices_simultaneous_events")]
      } else if (stats_attr_method == "pe") {
        attr(reh,"dyadID")  <- unlist(attr(reh,"dyadID"))[start_stop[1]:start_stop[2]]
        attr(reh,"actor1ID") <- unlist(attr(reh,"actor1ID"))[start_stop[1]:start_stop[2]]
        attr(reh,"actor2ID") <- unlist(attr(reh,"actor2ID"))[start_stop[1]:start_stop[2]]
      }
      reh$intereventTime <- reh$intereventTime[start_stop[1]:start_stop[2]]
      reh$M <- diff(start_stop) + 1
      if (length(reh$omit_dyad) > 0)
        reh$omit_dyad$time <- reh$omit_dyad$time[start_stop[1]:start_stop[2]]
    } else if (all(inherits(stats, c("remstats","aomstats"), TRUE))) {
      stop("'remstats' object supplied cannot work for tie-oriented modeling")
    }
    if (length(attr(reh,"dyadID")) != dim(stats)[3])
      stop("the number of time points (or number of events) doesn't match the (row) dimension of the 'remstats' object")
  }

  # -- actor-oriented -----------------------------------------------------------
  if (model == "actor") {
    model_formula <- list()
    if (all(inherits(stats, c("remstats","aomstats"), TRUE))) {
      variables_rate <- variables_choice <- NULL
      if (!is.null(stats$sender_stats)) {
        variables_rate <- dimnames(stats$sender_stats)[[3]]
        model_formula[["rate_model_formula"]] <- if (!is.null(attr(stats,"formula")$rate))
          attr(stats,"formula")$rate else
            stats::as.formula(paste("~ ", paste(variables_rate, collapse=" + ")))
        if (any(tolower(variables_rate) %in% c("baseline")))
          where_is_baseline <- which(variables_rate == "baseline")
        stats$sender_stats <- aperm(stats$sender_stats, perm = c(2,3,1))
      }
      if (!is.null(stats$receiver_stats)) {
        variables_choice <- dimnames(stats$receiver_stats)[[3]]
        model_formula[["choice_model_formula"]] <- if (!is.null(attr(stats,"formula")$choice))
          attr(stats,"formula")$choice else
            stats::as.formula(paste("~ ", paste(variables_choice, collapse=" + ")))
        stats$receiver_stats <- aperm(stats$receiver_stats, perm = c(2,3,1))
      }
      variables_names <- list(sender_model = variables_rate, receiver_model = variables_choice)
      if (stats_attr_method == "pt") {
        attr(reh,"actor1ID") <- attr(reh,"actor1ID")[start_stop[1]:start_stop[2]]
        attr(reh,"actor2ID") <- attr(reh,"actor2ID")[start_stop[1]:start_stop[2]]  # keep as list
        if (is.null(reh$E)) reh$E <- reh$M
        if (length(reh$omit_dyad) > 0) {
          if (!is.null(stats$receiver_stats)) {
            start_stop_time <- unique(reh$edgelist$time)[start_stop]
            lb_time <- min(which(reh$edgelist$time >= start_stop_time[1]))
            ub_time <- max(which(reh$edgelist$time <= start_stop_time[2]))
            omit_dyad_receiver <- list(time = reh$omit_dyad$time[lb_time:ub_time],
                                       riskset = reh$omit_dyad$riskset)
          }
          if (!is.null(stats$sender_stats)) {
            reh$omit_dyad$time <- if (!is.null(attr(reh,"indices_simultaneous_events")))
              reh$omit_dyad$time[-attr(reh,"indices_simultaneous_events")][start_stop[1]:start_stop[2]]
            else
              reh$omit_dyad$time[start_stop[1]:start_stop[2]]
          }
        }
      } else if (stats_attr_method == "pe") {
        attr(reh,"actor1ID") <- unlist(attr(reh,"actor1ID"))[start_stop[1]:start_stop[2]]
        attr(reh,"actor2ID") <- unlist(attr(reh,"actor2ID"))[start_stop[1]:start_stop[2]]
        if (length(reh$omit_dyad) > 0) {
          if (!is.null(stats$receiver_stats))
            omit_dyad_receiver <- list(time = reh$omit_dyad$time[start_stop[1]:start_stop[2]],
                                       riskset = reh$omit_dyad$riskset)
          reh$omit_dyad$time <- reh$omit_dyad$time[start_stop[1]:start_stop[2]]
        }
      }
      if (!ordinal) reh$intereventTime <- reh$intereventTime[start_stop[1]:start_stop[2]]
      reh$M <- diff(start_stop) + 1
      no_correct_dimensions <- FALSE
      if (!is.null(stats$sender_stats) &&
          length(attr(reh,"actor1ID")) != dim(stats$sender_stats)[3])
        no_correct_dimensions <- TRUE
      if (!is.null(stats$receiver_stats) &&
          length(attr(reh,"actor2ID")) != dim(stats$receiver_stats)[3])
        no_correct_dimensions <- TRUE
      if (no_correct_dimensions)
        stop("the number of time points (or number of events) doesn't match the (row) dimension of the 'remstats' object")
    } else if (all(inherits(stats, c("remstats","tomstats"), TRUE))) {
      stop("'remstats' object supplied cannot work for actor-oriented modeling")
    }
  }

  if (ordinal) reh$intereventTime <- c(1)

  list(
    reh = reh, stats = stats, model = model, ordinal = ordinal, ncores = ncores,
    stats_attr_method = stats_attr_method, start_stop = start_stop,
    model_formula = model_formula, variables_names = variables_names,
    where_is_baseline = where_is_baseline, omit_dyad_receiver = omit_dyad_receiver
  )
}

# ══════════════════════════════════════════════════════════════════════════
# diagnostics.remstimate() — trimmed to call the shared prep
# ══════════════════════════════════════════════════════════════════════════

#' @describeIn diagnostics diagnostics of a 'remstimate' object
#' @method diagnostics remstimate
#' @export
diagnostics.remstimate <- function(object, reh, stats, top_pct = 0.05, surprise_threshold = 0.2, ...) {
  prep <- .diagnostics_prepare(object, reh, stats)
  reh                <- prep$reh
  stats              <- prep$stats
  ordinal            <- prep$ordinal
  ncores             <- prep$ncores
  stats_attr_method  <- prep$stats_attr_method
  start_stop         <- prep$start_stop
  model_formula      <- prep$model_formula
  variables_names    <- prep$variables_names
  where_is_baseline  <- prep$where_is_baseline
  omit_dyad_receiver <- prep$omit_dyad_receiver

  # -- tie model diagnostics ----------------------------------------------------
  if (attr(object,"model") == "tie") {
    length_comparison <- length(object$coefficients) == dim(stats)[2]
    name_comparison <- length_comparison &&
      all(names(object$coefficients) == dimnames(stats)[[2]])
    if (!name_comparison | !length_comparison)
      stop("input 'object' not compatible with input 'stats'")
    variables_names   <- attr(object, "statistics")
    where_is_baseline <- attr(object, "where_is_baseline")
    select_vars    <- if (is.null(where_is_baseline)) seq_along(variables_names) else
      seq_along(variables_names)[-where_is_baseline]
    baseline_value <- if (is.null(where_is_baseline)) 0 else
      as.vector(object$coefficients)[where_is_baseline]
    stats <- if (length(select_vars) == 1)
      array(stats[,select_vars,], dim = c(dim(stats)[1], 1, dim(stats)[3]))
    else stats[,select_vars,]
    diagnostics <- list()
    diagnostics$residuals <- computeDiagnostics(
      pars      = as.vector(object$coefficients)[select_vars],
      stats     = stats,
      actor1    = list(),
      actor2    = list(),
      dyad      = attr(reh,"dyadID"),
      omit_dyad = reh$omit_dyad,
      model     = attr(reh,"model"),
      ncores    = ncores,
      baseline  = baseline_value,
      N         = reh$N
    )
    colnames(diagnostics$residuals$smoothing_weights) <- variables_names[select_vars]
    diagnostics$rates <- diagnostics$residuals$rates
    diagnostics$residuals$rates <- NULL
    if (!is.null(attr(object,"where_is_baseline")) & (length(select_vars) < 1))
      diagnostics$residuals <- NULL
    # recall
    obs_dyad_ids <- if (is.list(attr(reh,"dyadID"))) attr(reh,"dyadID") else
      as.list(attr(reh,"dyadID"))
    diagnostics$recall <- .recall_block_3d(
      pars     = as.vector(object$coefficients)[select_vars],
      baseline = baseline_value,
      stats_3d = stats,
      obs_ids  = obs_dyad_ids,
      top_pct  = top_pct
    )
    if (!is.null(diagnostics$recall$per_event)) {
      if (stats_attr_method == "pt") {
        diagnostics$recall$per_event$edgelist_id_row <-
          start_stop[1] - 1 + diagnostics$recall$per_event$event
      } else {
        diagnostics$recall$per_event$edgelist_id_row <- NA_integer_
        warning("edgelist_id_row cannot be resolved for stats method 'pe'; ",
                "returning NA.", call. = FALSE)
      }
    }
    diagnostics$surprises <- .surprises_from_recall(diagnostics$recall, surprise_threshold)
    diagnostics$surprise_threshold <- surprise_threshold
    diagnostics$surprise_offenders <- .offender_table(
      ids_surprises = diagnostics$surprises$obs_id,
      ids_all       = unlist(obs_dyad_ids),
      labels        = .dyad_labels(reh)
    )
    dm <- reh$index$dyad_map_active %||% reh$index$dyad_map
    if (!is.null(dm)) {
      id_col <- if ("dyadIDactive" %in% names(dm)) "dyadIDactive" else "dyadID"
      diagnostics$surprise_offenders_actor <- .tie_actor_offenders(
        diagnostics$surprises, diagnostics$recall$per_event, dm, id_col
      )
    }
    if (!is.null(diagnostics$surprises) && nrow(diagnostics$surprises) > 0L) {
      dyad_lab <- .dyad_labels(reh)
      diagnostics$surprises$dyad_label <-
        unname(dyad_lab[as.character(diagnostics$surprises$obs_id)])
    }

    if (!is.null(diagnostics$recall$per_event)) {
      if (stats_attr_method == "pt") {
        diagnostics$recall$per_event$edgelist_id_row <-
          start_stop[1] - 1 + diagnostics$recall$per_event$event
        diagnostics$recall$per_event$time <-
          unique(reh$edgelist$time)[diagnostics$recall$per_event$edgelist_id_row]
      } else {
        diagnostics$recall$per_event$edgelist_id_row <- NA_integer_
        diagnostics$recall$per_event$time <- NA_real_
        warning("edgelist_id_row cannot be resolved for stats method 'pe'; ",
                "returning NA.", call. = FALSE)
      }
    }

    # -- actor model diagnostics --------------------------------------------------
  } else if (attr(object,"model") == "actor") {
    compare_input_sender <- compare_input_receiver <- FALSE
    if (!is.null(stats[["sender_stats"]])) {
      lc <- length(object[["sender_model"]]$coefficients) == dim(stats[["sender_stats"]])[2]
      nc <- lc && all(names(object[["sender_model"]]$coefficients) ==
                        dimnames(stats[["sender_stats"]])[[2]])
      if (!nc | !lc) compare_input_sender <- TRUE
    }
    if (!is.null(stats[["receiver_stats"]])) {
      lc <- length(object[["receiver_model"]]$coefficients) == dim(stats[["receiver_stats"]])[2]
      nc <- lc && all(names(object[["receiver_model"]]$coefficients) ==
                        dimnames(stats[["receiver_stats"]])[[2]])
      if (!nc | !lc) compare_input_receiver <- TRUE
    }
    if (compare_input_sender | compare_input_receiver)
      stop("input 'object' not compatible with input 'stats'")
    variables_names   <- attr(object, "statistics")
    where_is_baseline <- attr(object, "where_is_baseline")
    senderRate   <- c(TRUE,  FALSE)
    which_model  <- c("sender_model",   "receiver_model")
    which_stats  <- c("sender_stats",   "receiver_stats")
    actor1ID_condition <- c(TRUE, FALSE)
    if (stats_attr_method == "pe")
      actor1ID_condition <- c(FALSE, FALSE)
    diagnostics <- list()
    for (i in 1:2) {
      if (is.null(stats[[which_stats[i]]])) next
      actor1ID_ls    <- if (actor1ID_condition[i]) attr(reh,"actor1ID") else
        unlist(attr(reh,"actor1ID"))
      omit_dyad_actor <- if (senderRate[i]) reh$omit_dyad else omit_dyad_receiver
      baseline_value <- 0
      select_vars    <- seq_len(dim(stats[[which_stats[i]]])[2])
      if (senderRate[i]) {
        baseline_value <- if (is.null(where_is_baseline)) 0 else
          as.vector(object[[which_model[i]]]$coefficients)[where_is_baseline]
        select_vars    <- if (is.null(where_is_baseline)) select_vars else
          select_vars[-where_is_baseline]
      }
      stats[[which_stats[i]]] <- if (length(select_vars) == 1)
        array(stats[[which_stats[i]]][,select_vars,],
              dim = c(dim(stats[[which_stats[i]]])[1], 1, dim(stats[[which_stats[i]]])[3]))
      else stats[[which_stats[i]]][,select_vars,]
      diagnostics[[which_model[i]]] <- list()
      diagnostics[[which_model[i]]]$residuals <- computeDiagnostics(
        pars       = as.vector(object[[which_model[i]]]$coefficients)[select_vars],
        stats      = stats[[which_stats[i]]],
        actor1     = actor1ID_ls,
        actor2     = attr(reh,"actor2ID"),
        dyad       = list(),
        omit_dyad  = omit_dyad_actor,
        model      = attr(reh,"model"),
        N          = reh$N,
        senderRate = senderRate[i],
        ncores     = ncores,
        baseline   = baseline_value
      )
      colnames(diagnostics[[which_model[i]]]$residuals$smoothing_weights) <-
        if (senderRate[i]) variables_names[["sender_model"]][select_vars] else
          variables_names[["receiver_model"]][select_vars]
      diagnostics[[which_model[i]]]$rates <- diagnostics[[which_model[i]]]$residuals$rates
      diagnostics[[which_model[i]]]$residuals$rates <- NULL
      if (senderRate[i] & !is.null(where_is_baseline) & length(select_vars) < 1)
        diagnostics[[which_model[i]]]$residuals <- NULL
      # recall
      # actor1ID and actor2ID are 1-indexed (C++ subtracts 1 internally)
      if (senderRate[i]) {
        obs_ids_rc   <- lapply(as.list(actor1ID_ls), as.integer)
        valid_ids_rc <- if (!is.null(reh$sender_riskset))
          rep(list(reh$sender_riskset), reh$M) else NULL
      } else {
        obs_ids_rc   <- lapply(attr(reh,"actor2ID"), function(x) as.integer(unlist(x)))
        valid_ids_rc <- if (!is.null(reh$receiver_riskset))
          lapply(as.list(attr(reh,"actor1ID")), function(senders)
            reh$receiver_riskset[[as.integer(senders[1])]])
        else NULL
      }
      diagnostics[[which_model[i]]]$recall <- .recall_block_3d(
        pars      = as.vector(object[[which_model[i]]]$coefficients)[select_vars],
        baseline  = baseline_value,
        stats_3d  = stats[[which_stats[i]]],
        obs_ids   = obs_ids_rc,
        valid_ids = valid_ids_rc,
        top_pct   = top_pct
      )
      if (!is.null(diagnostics[[which_model[i]]]$recall$per_event)) {
        if (stats_attr_method == "pt") {
          diagnostics[[which_model[i]]]$recall$per_event$edgelist_id_row <-
            start_stop[1] - 1 + diagnostics[[which_model[i]]]$recall$per_event$event
        } else {
          diagnostics[[which_model[i]]]$recall$per_event$edgelist_id_row <- NA_integer_
          warning("edgelist_id_row cannot be resolved for stats method 'pe'; ",
                  "returning NA.", call. = FALSE)
        }
      }
      diagnostics[[which_model[i]]]$surprises <- .surprises_from_recall(
        diagnostics[[which_model[i]]]$recall, surprise_threshold)
      diagnostics[[which_model[i]]]$surprise_threshold <- surprise_threshold
      diagnostics[[which_model[i]]]$surprise_offenders <- .offender_table(
        ids_surprises = diagnostics[[which_model[i]]]$surprises$obs_id,
        ids_all       = unlist(obs_ids_rc),
        labels        = .actor_labels(reh)
      )
      if (!senderRate[i]) {
        sender_ids_surprises <- .lookup_by_event(
          diagnostics[[which_model[i]]]$surprises$event,
          diagnostics[[which_model[i]]]$surprises$sub_index,
          attr(reh, "actor1ID")
        )
        diagnostics[[which_model[i]]]$surprises$sender_id <- sender_ids_surprises

        dyad_surprises <- .actor_pair_labels(reh, sender_ids_surprises,
                                             diagnostics[[which_model[i]]]$surprises$obs_id)
        dyad_all <- .actor_pair_labels(reh, unlist(attr(reh, "actor1ID")), unlist(obs_ids_rc))

        diagnostics[[which_model[i]]]$surprise_offenders_dyad <- .offender_table(
          ids_surprises = dyad_surprises,
          ids_all       = dyad_all
        )
        if (!is.null(diagnostics[[which_model[i]]]$surprises) &&
            nrow(diagnostics[[which_model[i]]]$surprises) > 0L) {
          act_lab <- .actor_labels(reh)
          diagnostics[[which_model[i]]]$surprises$obs_label <-
            unname(act_lab[as.character(diagnostics[[which_model[i]]]$surprises$obs_id)])
          if (!senderRate[i]) {
            diagnostics[[which_model[i]]]$surprises$sender_label <-
              unname(act_lab[as.character(diagnostics[[which_model[i]]]$surprises$sender_id)])
            diagnostics[[which_model[i]]]$surprises$dyad_label <- paste(
              diagnostics[[which_model[i]]]$surprises$sender_label, "->",
              diagnostics[[which_model[i]]]$surprises$obs_label)
          }
        }
      }
      if (!is.null(diagnostics[[which_model[i]]]$recall$per_event)) {
        if (stats_attr_method == "pt") {
          diagnostics[[which_model[i]]]$recall$per_event$edgelist_id_row <-
            start_stop[1] - 1 + diagnostics[[which_model[i]]]$recall$per_event$event
          diagnostics[[which_model[i]]]$recall$per_event$time <-
            unique(reh$edgelist$time)[diagnostics[[which_model[i]]]$recall$per_event$edgelist_id_row]
        } else {
          diagnostics[[which_model[i]]]$recall$per_event$edgelist_id_row <- NA_integer_
          diagnostics[[which_model[i]]]$recall$per_event$time <- NA_real_
          warning("edgelist_id_row cannot be resolved for stats method 'pe'; ",
                  "returning NA.", call. = FALSE)
        }
      }
    }
  }
  diagnostics$.reh.processed <- reh
  diagnostics$.reh.processed$stats.method <- stats_attr_method
  class(diagnostics) <- c("diagnostics","remstimate")
  return(diagnostics)
}


.interp_coef_at_time <- function(cf_block, mid_time, query_time) {
  coefmat <- cf_block$coefficients
  sep     <- cf_block$separated
  P <- ncol(coefmat)
  out <- matrix(NA_real_, length(query_time), P, dimnames = list(NULL, colnames(coefmat)))
  for (p in seq_len(P)) {
    y <- coefmat[, p]
    y[sep[, p]] <- NA
    valid <- which(!is.na(y))
    if (!length(valid)) next
    if (length(valid) == 1L) { out[, p] <- y[valid]; next }
    out[, p] <- stats::approx(x = mid_time[valid], y = y[valid],
                              xout = query_time, rule = 2, method = "linear")$y
  }
  out
}

.recall_block_3d_varying <- function(pars_mat, baseline_vec, stats_3d, obs_ids,
                                     valid_ids = NULL, top_pct = 0.05) {
  M    <- dim(stats_3d)[3L]
  rows <- vector("list", M)
  for (m in seq_len(M)) {
    obs <- obs_ids[[m]]
    if (!length(obs)) next
    valid <- if (is.null(valid_ids)) seq_len(dim(stats_3d)[1L]) else valid_ids[[m]]
    if (!length(valid)) next
    S    <- matrix(stats_3d[valid, , m], nrow = length(valid))
    pars <- pars_mat[m, ]
    base <- if (length(baseline_vec) > 1L) baseline_vec[m] else baseline_vec
    probs <- exp(as.numeric(base + S %*% pars))
    probs <- probs / sum(probs)
    D_t   <- length(valid)
    ord   <- order(probs, decreasing = TRUE)
    pos_all <- match(obs, valid); keep <- which(!is.na(pos_all))
    if (!length(keep)) next
    pos <- pos_all[keep]; sub_index <- keep
    rnks <- match(pos, ord); keep2 <- which(!is.na(rnks))
    if (!length(keep2)) next
    rnks <- rnks[keep2]; pos <- pos[keep2]; sub_index <- sub_index[keep2]
    obs_probs <- probs[pos]
    rows[[m]] <- data.frame(
      event = m, sub_index = sub_index, obs_id = obs[sub_index],
      rel_rank = 1 - rnks / D_t, cum_prob = cumsum(probs[ord])[rnks],
      obs_prob = obs_probs, prob_ratio = obs_probs * D_t, log_loss = -log(obs_probs)
    )
  }
  pe <- do.call(rbind, rows)
  list(per_event = pe, summary = data.frame(
    mean_rel_rank = mean(pe$rel_rank), median_rel_rank = median(pe$rel_rank),
    mean_cum_prob = mean(pe$cum_prob), mean_prob_ratio = mean(pe$prob_ratio),
    mean_log_loss = mean(pe$log_loss), top_pct = top_pct,
    top_pct_prop = mean(pe$rel_rank >= 1 - top_pct)
  ))
}

#' @export
diagnostics.remstimate_window <- function(object, reh, stats, top_pct = 0.05, k = 10, ...) {
  is_tie <- object$type == "tie"
  cf <- coef.remstimate_window(object, k = k)
  mid_time <- rowMeans(object$windows[, c("start_time", "end_time")])

  anchor_fit <- Find(function(f) inherits(f, "remstimate"), object$fits)
  if (is.null(anchor_fit)) stop("No successful fits in this remstimate_window object.", call. = FALSE)

  prep <- .diagnostics_prepare(anchor_fit, reh, stats)
  reh <- prep$reh; stats <- prep$stats
  stats_attr_method <- prep$stats_attr_method; start_stop <- prep$start_stop

  time_axis   <- if (stats_attr_method == "pt") unique(reh$edgelist$time) else reh$edgelist$time
  event_times <- time_axis[seq_len(reh$M) + (start_stop[1] - 1L)]

  .varying_block <- function(cf_block, stats_3d, obs_ids, valid_ids = NULL) {
    var_names <- colnames(cf_block$coefficients)
    wb  <- which(tolower(var_names) == "baseline")
    sel <- if (length(wb)) seq_along(var_names)[-wb] else seq_along(var_names)

    pars_mat <- .interp_coef_at_time(
      list(coefficients = cf_block$coefficients[, sel, drop = FALSE],
           separated    = cf_block$separated[,    sel, drop = FALSE]),
      mid_time, event_times)
    baseline_vec <- if (length(wb))
      .interp_coef_at_time(
        list(coefficients = cf_block$coefficients[, wb, drop = FALSE],
             separated    = cf_block$separated[,    wb, drop = FALSE]),
        mid_time, event_times)[, 1]
    else 0

    stats_sel <- if (length(sel) == 1L)
      array(stats_3d[, sel, ], dim = c(dim(stats_3d)[1], 1, dim(stats_3d)[3]))
    else stats_3d[, sel, , drop = FALSE]

    recall <- .recall_block_3d_varying(pars_mat, baseline_vec, stats_sel, obs_ids, valid_ids, top_pct)
    recall$per_event$time <- event_times[recall$per_event$event]
    recall
  }

  if (is_tie) {
    obs_dyad_ids <- if (is.list(attr(reh,"dyadID"))) attr(reh,"dyadID") else as.list(attr(reh,"dyadID"))
    out <- list(recall = .varying_block(cf, stats, obs_dyad_ids), windows = object$windows, type = "tie")
  } else {
    out <- list(windows = object$windows, type = "actor")
    if (!is.null(stats$sender_stats)) {
      obs_ids   <- lapply(as.list(attr(reh,"actor1ID")), as.integer)
      valid_ids <- if (!is.null(reh$sender_riskset)) rep(list(reh$sender_riskset), reh$M) else NULL
      out$sender <- list(recall = .varying_block(cf$sender, stats$sender_stats, obs_ids, valid_ids))
    }
    if (!is.null(stats$receiver_stats)) {
      obs_ids   <- lapply(attr(reh,"actor2ID"), function(x) as.integer(unlist(x)))
      valid_ids <- if (!is.null(reh$receiver_riskset))
        lapply(as.list(attr(reh,"actor1ID")), function(s) reh$receiver_riskset[[as.integer(s[1])]])
      else NULL
      out$receiver <- list(recall = .varying_block(cf$receiver, stats$receiver_stats, obs_ids, valid_ids))
    }
  }
  class(out) <- c("diagnostics.remstimate_window", "list")
  out
}

.plot_recall_window <- function(recall, windows, title_prefix = "") {
  pe <- recall$per_event
  if (is.null(pe) || nrow(pe) == 0L) {
    plot.new(); title(main = paste0(title_prefix, "Recall (no data)")); return(invisible())
  }

  old_par <- par(mar = c(4, 4, 5, 2) + 0.1)
  on.exit(par(old_par))

  plot(pe$event, pe$rel_rank, xaxt = "n",
       xlab = "Event", ylab = "Relative rank (0 = bottom, 1 = top)", main = "",
       ylim = c(0, 1), pch = 16, cex = 0.6, col = grDevices::rgb(0, 0, 0, 0.35))

  # tryCatch(
  #   lines(smooth.spline(x = pe$event, y = pe$rel_rank, cv = FALSE), lwd = 3.5, col = 2),
  #   error = function(e) NULL
  # )

  ticks <- pretty(pe$event)
  ticks <- ticks[ticks >= min(pe$event) & ticks <= max(pe$event)]
  axis(1, at = ticks)
  time_at_ticks <- stats::approx(pe$event, pe$time, xout = ticks, rule = 2)$y
  axis(3, at = ticks, labels = round(time_at_ticks, 1))
  mtext("Time", side = 3, line = 2.2, cex = 0.8)

  med <- stats::median(pe$rel_rank, na.rm = TRUE)
  abline(h = med, lty = 2, col = "blue", lwd = 1.5)

  ord   <- order(pe$event)
  x_ord <- pe$event[ord]
  y_ord <- pe$rel_rank[ord]

  k <- max(3L, min(length(y_ord) - 1L, round(length(y_ord) * 0.1)))
  if (k %% 2 == 0) k <- k + 1L   # runmed requires an odd bandwidth

  y_med <- tryCatch(stats::runmed(y_ord, k = k), error = function(e) NULL)
  if (!is.null(y_med)) {
    sm <- tryCatch(stats::smooth.spline(x_ord, y_med), error = function(e) NULL)
    if (!is.null(sm)) lines(sm$x, sm$y, col = "firebrick", lwd = 4)
    else               lines(x_ord, y_med, col = "firebrick", lwd = 4)  # fallback if spline fails
  }

  if (nrow(windows) > 1)
    abline(v = windows$end_event[-nrow(windows)], lty = 3, col = "grey70")

  legend("bottomright", inset = 0.02, bg = "white", cex = 0.8,
         legend = c("Rank", "Median rank", "Smoothed trend"),
         lty = c(NA, 2, 1), pch = c(16, NA, NA),
         col = c(grDevices::rgb(0, 0, 0, 0.35), "blue", "firebrick"), lwd = c(NA, 1.5, 2))

  mtext(paste0(title_prefix, sprintf("Recall (median rank = %.3f)", med)),
        side = 3, line = 3.3, cex = 1.1)
}

#' @export
plot.diagnostics.remstimate_window <- function(x, ...) {
  if (x$type == "tie") {
    .plot_recall_window(x$recall, x$windows)
  } else {
    if (!is.null(x$sender))   .plot_recall_window(x$sender$recall,   x$windows, "sender: ")
    if (!is.null(x$sender) && !is.null(x$receiver) && interactive())
      invisible(readline("Press <Enter> for receiver model plot..."))
    if (!is.null(x$receiver)) .plot_recall_window(x$receiver$recall, x$windows, "receiver: ")
  }
  invisible(x)
}


# -- print.diagnostics ---------------------------------------------------------

#' @method print diagnostics
#' @export
print.diagnostics <- function(x, ...) {
  reh   <- x$.reh.processed
  model <- reh$meta$model
  M     <- reh$M
  N     <- reh$N

  cat("Diagnostics for a Relational Event Model\n")
  cat(strrep("-", 42), "\n", sep = "")
  cat(sprintf("%-12s: %s\n", "Model",  model))
  cat(sprintf("%-12s: %d\n", "Actors", N))
  cat(sprintf("%-12s: %d\n", "Events", M))

  .print_submodel <- function(sub, label) {
    if (is.null(sub)) return(invisible(NULL))
    cat("\n", label, "\n", sep = "")
    if (!is.null(sub$residuals))
      cat("  Statistics :", paste(colnames(sub$residuals$smoothing_weights),
                                   collapse = ", "), "\n")
    if (!is.null(sub$recall)) {
      rs  <- sub$recall$summary
      pct <- round(rs$top_pct_prop * 100, 1)
      low_str <- ""
      if (!is.null(sub$surprise_threshold)) {
        n_pe    <- nrow(sub$recall$per_event)
        n_sur   <- if (is.null(sub$surprises)) 0L else nrow(sub$surprises)
        low_pct <- if (n_pe > 0) round(100 * n_sur / n_pe, 1) else NA
        low_str <- sprintf("  |  lowest %g%% = %s%%", sub$surprise_threshold * 100, low_pct)
      }
      cat(sprintf("  Recall     : mean rank = %.3f  |  prob ratio = %.2f  |  top %g%% = %s%%%s\n",
                  rs$mean_rel_rank, rs$mean_prob_ratio, rs$top_pct * 100, pct, low_str))
    }
  }

  if (model == "tie") {
    .print_submodel(x, "Tie model")
  } else if (model == "actor") {
    .print_submodel(x$sender_model,   "Sender rate model")
    .print_submodel(x$receiver_model, "Receiver choice model")
  }
  invisible(x)
}

# -- summary.diagnostics -------------------------------------------------------

#' @method summary diagnostics
#' @export
summary.diagnostics <- function(object, ...) {
  reh   <- object$.reh.processed
  model <- reh$meta$model

  # helper: Schoenfeld residuals summary per statistic
  .resid_summary <- function(resid) {
    if (is.null(resid)) return(NULL)
    sr <- resid$standardized_residuals
    # sr is either a list (pt) or a matrix (pe)
    mat <- if (is.list(sr)) do.call(rbind, lapply(sr, as.matrix)) else as.matrix(sr)
    colnames(mat) <- colnames(resid$smoothing_weights)
    data.frame(
      statistic = colnames(mat),
      mean      = round(colMeans(mat, na.rm = TRUE), 4),
      sd        = round(apply(mat, 2, sd, na.rm = TRUE), 4),
      row.names = NULL
    )
  }

  # helper: recall summary formatted
  .recall_summary <- function(rc) {
    if (is.null(rc)) return(NULL)
    rs <- rc$summary
    data.frame(
      mean_rel_rank   = round(rs$mean_rel_rank, 4),
      median_rel_rank = round(rs$median_rel_rank, 4),
      mean_cum_prob   = round(rs$mean_cum_prob, 4),
      mean_prob_ratio = round(rs$mean_prob_ratio, 4),
      mean_log_loss   = round(rs$mean_log_loss, 4),
      top_pct         = rs$top_pct,
      top_pct_prop    = round(rs$top_pct_prop, 4),
      row.names       = NULL
    )
  }

  .print_submodel_summary <- function(sub, label) {
    if (is.null(sub)) return(invisible(NULL))
    cat("\n", label, "\n", strrep("-", nchar(label)), "\n", sep = "")

    rs  <- .resid_summary(sub$residuals)
    if (!is.null(rs)) {
      cat("\nSchoenfeld residuals:\n")
      print(rs, row.names = FALSE)
    }

    rc <- .recall_summary(sub$recall)
    if (!is.null(rc)) {
      cat("\nRecall:\n")
      print(rc, row.names = FALSE)
    }
  }

  cat("Diagnostics for a Relational Event Model\n")
  cat(strrep("-", 42), "\n", sep = "")
  cat(sprintf("%-12s: %s\n", "Model",  model))
  cat(sprintf("%-12s: %d\n", "Actors", reh$N))
  cat(sprintf("%-12s: %d\n", "Events", reh$M))

  if (model == "tie") {
    .print_submodel_summary(object, "Tie model")
  } else if (model == "actor") {
    .print_submodel_summary(object$sender_model,   "Sender rate model")
    .print_submodel_summary(object$receiver_model, "Receiver choice model")
  }
  invisible(object)
}


# diagnostics voor GLMM, GLMNET en MIXREM backends

# GLMM / GLMNET: gooi de subklasse weg en roep diagnostics.remstimate aan;
# alle benodigde attributen (statistics, where_is_baseline, ncores) zijn
# al gezet door .remstimate_wrap()

#' @export
#' @method diagnostics remstimate_glmm
diagnostics.remstimate_glmm <- function(object, reh, stats, top_pct = 0.05, ...) {
  class(object) <- "remstimate"
  diag_obj <- diagnostics(object, reh, stats, top_pct = top_pct, ...)
  diag_obj$.engine <- attr(object, "engine") %||% "lme4"
  diag_obj
}

#' @export
#' @method diagnostics remstimate_glmnet
diagnostics.remstimate_glmnet <- function(object, reh, stats, top_pct = 0.05, ...) {

  object2 <- object
  class(object2) <- "remstimate"
  diag_obj <- diagnostics(object2, reh, stats, top_pct = top_pct, ...)
  diag_obj$.alpha      <- object$alpha
  diag_obj$.lambda_sel <- object$lambda_sel
  diag_obj
}

#' @export
#' @method diagnostics remstimate_mixrem
diagnostics.remstimate_mixrem <- function(object, reh, stats, top_pct = 0.05, ...) {
  if (attr(object, "model") == "tie")
    .mixrem_diag_tie(object, reh, stats, top_pct)
  else
    list(
      sender_model   = .mixrem_diag_actor(object, reh, stats, "sender_model",   top_pct),
      receiver_model = .mixrem_diag_actor(object, reh, stats, "receiver_model", top_pct)
    )
}

.mixrem_diag_tie <- function(object, reh, stats, top_pct) {

  reh         <- denormalize_reh(reh)
  subset_idx  <- as.numeric(unlist(attr(stats, "subset") %||% list(1, reh$M)))
  stats_perm  <- aperm(stats, c(2, 3, 1))

  dyad_ids     <- attr(reh, "dyadIDactive") %||% attr(reh, "dyadID")
  dyad_ids     <- dyad_ids[subset_idx[1]:subset_idx[2]]
  waargenomen  <- if (is.list(dyad_ids)) dyad_ids else as.list(dyad_ids)

  stat_names        <- dimnames(stats)[[3]]
  waar_is_baseline  <- .remstimate_find_baseline(stat_names)
  selectie          <- if (is.null(waar_is_baseline)) seq_along(stat_names) else
    seq_along(stat_names)[-waar_is_baseline]

  stats_sub <- if (length(selectie) == 1L)
    array(stats_perm[, selectie, ], dim = c(dim(stats_perm)[1], 1, dim(stats_perm)[3]))
  else
    stats_perm[, selectie, , drop = FALSE]

  k <- ncol(object$coefficients)

  per_component <- setNames(lapply(seq_len(k), function(j) {
    basis_j <- if (is.null(waar_is_baseline)) 0 else object$coefficients[waar_is_baseline, j]
    .recall_block_3d(
      pars     = as.vector(object$coefficients[selectie, j]),
      baseline = basis_j,
      stats_3d = stats_sub,
      obs_ids  = waargenomen,
      top_pct  = top_pct
    )
  }), paste0("Component.", seq_len(k)))

  # gecombineerde recall via posterior component assignments
  gecombineerd <- tryCatch({
    clusters <- flexmix::clusters(object$backend_fit)
    M    <- dim(stats_sub)[3L]
    rijen <- vector("list", M)
    for (m in seq_len(M)) {
      d <- waargenomen[[m]]
      if (!length(d)) next
      j     <- clusters[d[1L]]
      basis_j <- if (is.null(waar_is_baseline)) 0 else object$coefficients[waar_is_baseline, j]
      pars_j  <- as.vector(object$coefficients[selectie, j])
      geldig  <- seq_len(dim(stats_sub)[1L])
      S       <- matrix(stats_sub[geldig, , m], nrow = length(geldig))
      kansen  <- exp(as.numeric(basis_j + S %*% pars_j))
      kansen  <- kansen / sum(kansen)
      volgorde <- order(kansen, decreasing = TRUE)
      pos     <- match(d, geldig);       pos  <- pos[!is.na(pos)]
      rnks    <- match(pos, volgorde);   rnks <- rnks[!is.na(rnks)]
      if (!length(rnks)) next
      obs_probs_j <- kansen[pos[!is.na(match(pos, seq_along(kansen)))]]
      rijen[[m]] <- data.frame(
        event      = m,
        rel_rank   = 1 - rnks / length(geldig),
        cum_prob   = cumsum(kansen[volgorde])[rnks],
        obs_prob   = obs_probs_j,
        prob_ratio = obs_probs_j * length(geldig),
        log_loss   = -log(obs_probs_j),
        component  = j
      )
    }
    pe <- do.call(rbind, rijen)
    list(per_event = pe,
         summary   = data.frame(
           mean_rel_rank   = mean(pe$rel_rank),
           median_rel_rank = median(pe$rel_rank),
           mean_cum_prob   = mean(pe$cum_prob),
           mean_prob_ratio = mean(pe$prob_ratio),
           mean_log_loss   = mean(pe$log_loss),
           top_pct         = top_pct,
           top_pct_prop    = mean(pe$rel_rank >= 1 - top_pct)
         ))
  }, error = function(e) NULL)

  out <- list(per_component = per_component, combined = gecombineerd,
              k = k, prior_probs = object$prior_probs)
  out$.reh.processed <- reh
  class(out) <- c("diagnostics_mixrem", "diagnostics", "remstimate")
  out
}

.mixrem_diag_actor <- function(object, reh, stats, deel, top_pct) {
  sub_obj <- object[[deel]]
  if (is.null(sub_obj)) return(NULL)

  reh        <- denormalize_reh(reh)
  subset_idx <- as.numeric(unlist(attr(stats, "subset") %||% list(1, reh$M)))
  is_sender  <- (deel == "sender_model")

  stats_arr  <- if (is_sender) aperm(stats$sender_stats, c(2, 3, 1)) else
    aperm(stats$receiver_stats, c(2, 3, 1))
  stat_names <- if (is_sender) dimnames(stats$sender_stats)[[3]] else
    dimnames(stats$receiver_stats)[[3]]

  location_baseline <- if (is_sender) .remstimate_find_baseline(stat_names) else NULL
  selectie         <- if (is.null(location_baseline)) seq_along(stat_names) else
    seq_along(stat_names)[-location_baseline]
  stats_sub <- if (length(selectie) == 1L)
    array(stats_arr[, selectie, ], dim = c(dim(stats_arr)[1], 1, dim(stats_arr)[3]))
  else stats_arr[, selectie, , drop = FALSE]

  bereik <- subset_idx[1]:subset_idx[2]
  if (is_sender) {
    waargenomen <- lapply(as.list(attr(reh, "actor1ID")[bereik]), as.integer)
    geldig_ids  <- if (!is.null(reh$sender_riskset))
      rep(list(reh$sender_riskset), length(bereik)) else NULL
  } else {
    waargenomen <- lapply(attr(reh, "actor2ID")[bereik], function(x) as.integer(unlist(x)))
    geldig_ids  <- if (!is.null(reh$receiver_riskset))
      lapply(as.list(attr(reh, "actor1ID")[bereik]),
             function(s) reh$receiver_riskset[[as.integer(s[1])]])
    else NULL
  }

  k <- ncol(sub_obj$coefficients)
  per_component <- setNames(lapply(seq_len(k), function(j) {
    basis_j <- if (is.null(location_baseline) || !is_sender) 0 else
      sub_obj$coefficients[location_baseline, j]
    .recall_block_3d(
      pars      = as.vector(sub_obj$coefficients[selectie, j]),
      baseline  = basis_j,
      stats_3d  = stats_sub,
      obs_ids   = waargenomen,
      valid_ids = geldig_ids,
      top_pct   = top_pct
    )
  }), paste0("Component.", seq_len(k)))

  list(per_component = per_component, k = k)
}

#' @export
print.diagnostics_mixrem <- function(x, ...) {
  reh <- x$.reh.processed
  cat("Diagnostics — Mixture REM (k =", x$k, ")\n")
  cat(sprintf("Actors: %d  Events: %d\n", reh$N, reh$M))
  cat("\nMixing proportions:\n")
  pp <- setNames(round(x$prior_probs, 4), paste0("Component.", seq_along(x$prior_probs)))
  print(pp)
  cat("\nRecall per component:\n")
  for (nm in names(x$per_component)) {
    rs <- x$per_component[[nm]]$summary
    cat(sprintf("  %-14s  mean rank = %.3f  |  prob ratio = %.2f  |  top 5%% = %.1f%%\n",
                nm, rs$mean_rel_rank, rs$mean_prob_ratio, rs$top_pct_prop * 100))
  }
  if (!is.null(x$combined)) {
    rs <- x$combined$summary
    cat(sprintf("  %-14s  mean rank = %.3f  |  prob ratio = %.2f  |  top 5%% = %.1f%%\n",
                "Combined", rs$mean_rel_rank, rs$mean_prob_ratio, rs$top_pct_prop * 100))
  }
  invisible(x)
}




