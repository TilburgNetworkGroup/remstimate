# -- recall helper -------------------------------------------------------------
.recall_block <- function(pars, baseline, stats_3d, obs_ids,
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
    ord   <- order(probs, decreasing = TRUE)
    pos   <- match(obs, valid); pos <- pos[!is.na(pos)]
    rnks  <- match(pos, ord);   rnks <- rnks[!is.na(rnks)]
    if (!length(rnks)) next
    rows[[m]] <- data.frame(
      event    = m,
      rel_rank = 1 - rnks / length(valid),
      cum_prob = cumsum(probs[ord])[rnks]
    )
  }
  pe <- do.call(rbind, rows)
  list(
    per_event = pe,
    summary   = data.frame(
      mean_rel_rank   = mean(pe$rel_rank),
      median_rel_rank = median(pe$rel_rank),
      mean_cum_prob   = mean(pe$cum_prob),
      top_pct         = top_pct,
      top_pct_prop    = mean(pe$rel_rank >= 1 - top_pct)
    )
  )
}

# diagnostics
#' @title Compute the diagnostics of a \code{remstimate} object
#' @param object is a \code{remstimate} object.
#' @param reh is a \code{remify} object, the same used for the 'remstimate' object.
#' @param stats is a \code{remstats} object, the same used for the 'remstimate' object.
#' @param top_pct numeric scalar in (0,1): threshold for the top-percentage recall summary (default 0.05).
#' @param ... further arguments to be passed to the 'diagnostics' method.
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

#' @describeIn diagnostics diagnostics of a 'remstimate' object
#' @method diagnostics remstimate
#' @export
diagnostics.remstimate <- function(object, reh, stats, top_pct = 0.05, ...) {
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
  # active riskset
  if (attr(reh, "riskset") == "active") {
    reh$D <- reh$activeD
    if (model == "tie") {
      attr(reh, "dyadID") <- attr(reh, "dyadIDactive")
      reh$omit_dyad <- list()
    }
  }
  ordinal <- attr(reh, "ordinal")
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
    diagnostics$recall <- .recall_block(
      pars     = as.vector(object$coefficients)[select_vars],
      baseline = baseline_value,
      stats_3d = stats,
      obs_ids  = obs_dyad_ids,
      top_pct  = top_pct
    )

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
      diagnostics[[which_model[i]]]$recall <- .recall_block(
        pars      = as.vector(object[[which_model[i]]]$coefficients)[select_vars],
        baseline  = baseline_value,
        stats_3d  = stats[[which_stats[i]]],
        obs_ids   = obs_ids_rc,
        valid_ids = valid_ids_rc,
        top_pct   = top_pct
      )
    }
  }
  diagnostics$.reh.processed <- reh
  diagnostics$.reh.processed$stats.method <- stats_attr_method
  class(diagnostics) <- c("diagnostics","remstimate")
  return(diagnostics)
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
      cat(sprintf("  Recall     : mean rank = %.3f  |  top %g%% = %s%%\n",
                  rs$mean_rel_rank, rs$top_pct * 100, pct))
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
