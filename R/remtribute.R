#' @title remtribute
#'
#' @description Model an event attribute (type, weight, or any
#'   event-level attribute) conditional on the observed dyad sequence.
#'   This implements the conditional decomposition of Brandes, Lerner &
#'   Snijders (2009): first model which dyads interact (via
#'   \code{\link{remstimate}}), then model the event attribute given
#'   the observed dyad.
#'
#'   The \code{reh} object can be either a tie-oriented or
#'   an actor-oriented REM regardless of how the dyad sequence
#'   was modeled upstream. Because the event attribute is modeled given the
#'   observed dyad (which can be modeled using a \code{tie} model or \code{actor}
#'   model), the attribute model always uses \emph{tie-oriented} (dyad-level)
#'   statistics as predictor variables.
#'
#'   Statistics can be supplied in two ways: either pass a pre-computed
#'   \code{tomstats} object via \code{stats}, or pass a one-sided
#'   formula via \code{effects} and let \code{remtribute} compute the
#'   statistics internally. When the \code{reh} was created with
#'   \code{model = "actor"}, the \code{effects} route is recommended
#'   because \code{tomstats} requires a tie-oriented \code{reh}.
#'
#' @param reh A \code{remify} object whose \code{$edgelist} contains the
#'   attribute column (carried through via \code{event_type} or
#'   \code{event_attributes} in \code{\link[remify]{remify}}).
#'   Accepted for both \code{model = "tie"} and \code{model = "actor"}.
#' @param stats A \code{tomstats} object (M x D x p), or \code{NULL}.
#'   When provided, statistics are subsetted internally to the observed
#'   dyad at each event. Ignored when \code{effects} is supplied.
#' @param effects A one-sided formula with tie-model effects, e.g.
#'   \code{~ inertia() + reciprocity()}. When provided, \code{tomstats}
#'   is called internally to compute the statistics. This is the
#'   recommended interface when \code{reh} uses an actor-oriented model,
#'   because it avoids the need for a separate tie-oriented \code{reh}.
#' @param attribute Character string: column name in
#'   \code{reh$edgelist} that contains the event attribute to model.
#'   Default \code{"type"}.
#' @param attribute_type The type of the attribute:
#'   \code{"nominal"} (unordered categories, multinomial logit),
#'   \code{"ordinal"} (ordered categories, cumulative link model), or
#'   \code{"numeric"} (continuous, linear regression).
#' @param attr_actors Optional data frame of actor-level covariates,
#'   passed to \code{tomstats} when \code{effects} is used.
#' @param memory Memory type for statistics computation. Passed to
#'   \code{tomstats} when \code{effects} is used. Default \code{"full"}.
#' @param memory_value Memory parameter value. Passed to
#'   \code{tomstats} when \code{effects} is used.
#' @param ... Additional arguments passed to the fitting backend
#'   (\code{nnet::multinom}, \code{MASS::polr}, or \code{stats::glm}).
#'
#' @return An object of class \code{remtribute} containing the fitted
#'   model, coefficients, and metadata.
#'
#' @details
#' The function extracts the event attribute from the remified
#' edgelist and, for each event, subsets the statistics array to the
#' observed dyad. This yields an M x p design matrix of dyad-level
#' covariates. The attribute is then modeled as a function of
#' these covariates using the appropriate backend:
#' \itemize{
#'   \item \code{nominal}: \code{nnet::multinom}
#'   \item \code{ordinal}: \code{MASS::polr}
#'   \item \code{numeric}: \code{stats::glm} (Gaussian)
#' }
#'
#' Note that even when the upstream event model is actor-oriented
#' (sender rate + receiver choice), the attribute model uses
#' \code{tomstats} (tie-oriented) statistics. This is because
#' \code{aomstats} decomposes statistics into sender-level and
#' receiver-level arrays, while the attribute model conditions on the
#' full observed dyad and requires dyad-level covariates.
#'
#' @references
#' Brandes, U., Lerner, J., & Snijders, T. A. B. (2009). Networks
#'   evolving step by step: Statistical analysis of dyadic event data.
#'   In \emph{2009 International Conference on Advances in Social Network
#'   Analysis and Mining} (pp. 200--205). IEEE.
#'
#' Arena, G., Mulder, J., & Leenders, R. Th. A. J. (2023). Weighting
#'   the past: An extended relational event model for negative and
#'   positive events. \emph{Journal of the Royal Statistical Society
#'   Series A}, 189(1), 359--381.
#'
#' @examples
#' \donttest{
#' # --- Example 1: Tie model + pre-computed stats ---
#' if (requireNamespace("remdata", quietly = TRUE)) {
#'   data(hypertext, package = "remdata")
#'   # numeric attribute
#'   reh_text <- remify::remify(hypertext[1:1000,], model = "tie",
#'     directed = FALSE, riskset = "active",
#'     event_attributes = "duration")
#'   remstats_text <- remstats::remstats(reh_text,
#'     tie_effects = ~ inertia(scaling = "std") + totaldegreeDyad(scaling = "std"),
#'     first = 50)
#'   # fit model for event rate
#'   remstimate_tie <- remstimate(reh_text, remstats_text)
#'   # fit model for event attribute
#'   fit_attr <- remtribute(reh_text, stats = remstats_text, attribute = "duration",
#'    attribute_type = "numeric")
#'
#'   # nominal attribute
#'   hypertext$long <- hypertext$duration>20
#'   reh_text <- remify::remify(hypertext[1:1000,], model = "tie", directed = FALSE,
#'     riskset = "active", event_attributes = "long")
#'   remstats_text <- remstats::remstats(reh_text,
#'     tie_effects = ~ inertia(scaling = "std") +
#'     totaldegreeDyad(scaling = "std"),
#'     first = 50)
#'   remstimate_text <- remstimate(reh_text, remstats_text)
#'   summary(remstimate_text)
#'   fit_attr <- remtribute(reh_text, stats = remstats_text, attribute = "long",
#'     attribute_type = "nominal")
#'
#'   # --- Example 2: Actor model + effects formula ---
#'   reh_act <- remify::remify(edgelist = hypertext[1:1000,], model = "actor",
#'                           event_attributes = "duration")
#'   stats_act <- remstats::remstats(reh_act, sender_effects = ~ outdegreeSender(),
#'     receiver_effects = ~ inertia())
#'   fit_act <- remstimate(reh_act, stats_act)
#'
#'   # Use effects formula -- no need for a separate tie-model reh
#'   fit_type <- remtribute(reh_act, effects = ~ inertia() + reciprocity(),
#'     attribute = "duration",
#'     attribute_type = "numeric")
#'   summary(fit_type)
#' }
#' }
#' @export
remtribute <- function(reh,
                        stats = NULL,
                        effects = NULL,
                        attribute = "type",
                        attribute_type = c("nominal", "ordinal", "numeric"),
                        attr_actors = NULL,
                        memory = "full",
                        memory_value = Inf,
                        ...) {

  attribute_type <- match.arg(attribute_type)

  # ── 1. Validate inputs ──────────────────────────────────────────────────
  if (!inherits(reh, "remify"))
    stop("'reh' must be a 'remify' object.", call. = FALSE)

  if (is.null(stats) && is.null(effects))
    stop("Provide either 'stats' (a tomstats object) or 'effects' (a formula).",
         call. = FALSE)

  if (!is.null(stats) && !inherits(stats, c("tomstats", "tomstats_sampled")))
    stop("'stats' must be a 'tomstats' object.", call. = FALSE)

  if (!attribute %in% names(reh$edgelist))
    stop("Column '", attribute, "' not found in reh$edgelist. ",
         "Use event_type or event_attributes in remify() to carry it through.",
         call. = FALSE)

  # ── 1b. Compute stats from effects formula if needed ────────────────────
  if (!is.null(effects)) {
    if (!inherits(effects, "formula"))
      stop("'effects' must be a one-sided formula, e.g. ~ inertia().",
           call. = FALSE)

    model <- if (!is.null(reh$meta)) reh$meta$model else (attr(reh, "model") %||% "tie")

    if (model == "actor") {
      # Re-remify as tie model so tomstats can compute dyad-level stats
      directed <- if (!is.null(reh$meta)) isTRUE(reh$meta$directed) else isTRUE(attr(reh, "directed"))
      riskset_src <- if (!is.null(reh$meta)) (reh$meta$riskset_source %||% reh$meta$riskset) else attr(reh, "riskset")
      reh_tie <- remify::remify(
        edgelist = reh$edgelist[, intersect(names(reh$edgelist), c("time", "actor1", "actor2", "type", "weight", attribute))],
        model = "tie",
        directed = directed,
        ordinal = if (!is.null(reh$meta)) isTRUE(reh$meta$ordinal) else isTRUE(attr(reh, "ordinal")),
        actors = reh$actors %||% NULL,
        riskset = if (grepl("^active", riskset_src) || identical(riskset_src, "manual")) riskset_src else "full",
        event_type = if (!is.null(reh$meta) && isTRUE(reh$meta$with_type)) "type" else NULL,
        event_attributes = if (attribute != "type") attribute else NULL
      )
    } else {
      reh_tie <- reh
    }

    stats <- remstats::tomstats(
      effects, reh = reh_tie,
      attr_actors = attr_actors,
      memory = memory,
      memory_value = memory_value
    )

    # Use the tie-model reh for dyad indexing
    reh_for_dyads <- reh_tie
  } else {
    reh_for_dyads <- reh
  }

  # ── 2. Raw outcome vector (per event, in edgelist order) ────────────────
  y_full <- reh$edgelist[[attribute]]

  # ── 3. Get observed dyad indices (per stats time-point row) ─────────────
  riskset <- if (!is.null(reh_for_dyads$meta)) {
    reh_for_dyads$meta$riskset_source %||% reh_for_dyads$meta$riskset
  } else {
    attr(reh_for_dyads, "riskset")
  }
  reduced_riskset <- grepl("^active", riskset) || identical(riskset, "manual")

  if (!is.null(reh_for_dyads$ids)) {
    dyad_ids <- if (reduced_riskset) reh_for_dyads$ids$dyad_active
                else reh_for_dyads$ids$dyad
  } else {
    dyad_ids <- if (reduced_riskset)
      (attr(reh_for_dyads, "dyadIDactive") %||% attr(reh_for_dyads, "dyadID"))
    else attr(reh_for_dyads, "dyadID")
  }

  # `dyad_ids` holds the observed dyad(s) for each time-point row of `stats`.
  # With method = "pt" and simultaneous events it is a *list* (a vector of
  # dyads per row); with one event per time point it is an atomic vector.
  # Normalise to a list so the per-event expansion is uniform, then flatten
  # to one entry per event together with the stats row each event maps to.
  dyad_list   <- if (is.list(dyad_ids)) lapply(dyad_ids, as.integer)
                 else as.list(as.integer(dyad_ids))
  n_rows_full <- length(dyad_list)                 # rows over the FULL history
  n_events    <- sum(lengths(dyad_list))

  if (length(y_full) != n_events)
    stop("Cannot align attribute '", attribute, "' (", length(y_full),
         " edgelist rows) with the observed dyad sequence (", n_events,
         " events). 'stats' and 'reh' appear to describe different event ",
         "histories - rebuild both from the same 'reh'.", call. = FALSE)

  row_of_event  <- rep.int(seq_len(n_rows_full), lengths(dyad_list))
  dyad_of_event <- unlist(dyad_list, use.names = FALSE)

  # ── 3b. Restrict to the window that `stats` was computed on ─────────────
  # `stats` is already sliced to its window (dim(stats)[1] rows). The subset
  # attribute gives the start/stop *row* indices into the full history, so we
  # keep events whose row falls inside it and re-base the row index to 1.
  subset_attr <- attr(stats, "subset")
  start_stop  <- if (!is.null(subset_attr)) as.integer(unlist(subset_attr))
                 else c(1L, n_rows_full)

  keep          <- row_of_event >= start_stop[1] & row_of_event <= start_stop[2]
  stats_row     <- row_of_event[keep] - start_stop[1] + 1L
  dyad_of_event <- dyad_of_event[keep]
  y             <- y_full[keep]
  M             <- length(y)

  if (M == 0L)
    stop("No events fall inside the stats window (subset = ",
         paste(start_stop, collapse = ":"), ").", call. = FALSE)

  # Guard against stats/riskset dyad-dimension mismatch before indexing, so a
  # bad riskset fails here rather than with a cryptic subscript error.
  if (any(dyad_of_event > dim(stats)[2], na.rm = TRUE))
    stop("Observed dyad id exceeds the number of dyad columns in 'stats' (",
         dim(stats)[2], "). The remify and remstats risksets differ - ",
         "rebuild both from the same 'reh'.", call. = FALSE)

  # ── 3c. Coerce/validate outcome now that it is windowed ─────────────────
  if (attribute_type == "nominal") {
    y <- as.factor(y)
    if (nlevels(y) < 2L)
      stop("'nominal' requires at least 2 categories; found ", nlevels(y),
           ".", call. = FALSE)
  } else if (attribute_type == "ordinal") {
    y <- as.ordered(y)
    if (nlevels(y) < 2L)
      stop("'ordinal' requires at least 2 ordered categories; found ",
           nlevels(y), ".", call. = FALSE)
  } else {
    y <- as.numeric(y)
    if (all(is.na(y)))
      stop("All values in '", attribute, "' are NA.", call. = FALSE)
  }

  # ── 4. Build the design matrix: stats[row, dyad, ] per event ────────────
  stat_names <- dimnames(stats)[[3]]
  if (is.null(stat_names))
    stop("'stats' has no dimname labels on third dimension.", call. = FALSE)
  p <- length(stat_names)

  X <- matrix(NA_real_, nrow = M, ncol = p)
  for (i in seq_len(M)) {
    X[i, ] <- stats[stats_row[i], dyad_of_event[i], ]
  }
  colnames(X) <- stat_names

  # Remove constant columns: they are collinear with the intercept that the
  # backend already estimates. `baseline` is constant *by construction* (it is
  # the REM intercept, equal to 1 for every dyad), so drop it silently; only
  # report other statistics that turned out constant, since that is unexpected
  # and usually signals a degenerate effect at the observed dyads.
  const_cols <- apply(X, 2, function(col) length(unique(col)) == 1L)
  if (any(const_cols)) {
    is_baseline    <- tolower(stat_names) == "baseline"
    report_dropped <- const_cols & !is_baseline
    if (any(report_dropped))
      message("Dropping constant statistic(s): ",
              paste(stat_names[report_dropped], collapse = ", "))
    X <- X[, !const_cols, drop = FALSE]
    stat_names <- stat_names[!const_cols]
    p <- length(stat_names)
  }

  if (p == 0L)
    stop("No non-constant statistics remaining after subsetting to ",
         "observed dyads.", call. = FALSE)

  # ── 5. Fit model ────────────────────────────────────────────────────────
  df <- as.data.frame(X)
  df$.y <- y

  fml <- stats::as.formula(paste(".y ~", paste(stat_names, collapse = " + ")))

  fit <- switch(attribute_type,
    nominal = {
      if (!requireNamespace("nnet", quietly = TRUE))
        stop("Install 'nnet' for nominal outcome models.", call. = FALSE)
      nnet::multinom(fml, data = df, trace = FALSE, ...)
    },
    ordinal = {
      if (!requireNamespace("MASS", quietly = TRUE))
        stop("Install 'MASS' for ordinal outcome models.", call. = FALSE)
      MASS::polr(fml, data = df, Hess = TRUE, ...)
    },
    numeric = {
      stats::glm(fml, data = df, family = gaussian(), ...)
    }
  )

  # ── 6. Extract results ──────────────────────────────────────────────────
  coefs <- stats::coef(fit)
  vcov  <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  loglik <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)

  out <- list(
    coefficients    = coefs,
    vcov            = vcov,
    loglik          = loglik,
    fit             = fit,
    attribute       = attribute,
    attribute_type  = attribute_type,
    n_events        = M,
    stat_names      = stat_names,
    formula         = fml,
    data            = df
  )

  if (attribute_type == "nominal") {
    out$levels <- levels(y)
    out$n_levels <- nlevels(y)
  } else if (attribute_type == "ordinal") {
    out$levels <- levels(y)
    out$n_levels <- nlevels(y)
  }

  # Model fit statistics
  out$AIC <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)
  out$BIC <- tryCatch(stats::BIC(fit), error = function(e) NA_real_)

  structure(out, class = "remtribute")
}


# ══════════════════════════════════════════════════════════════════════════════
# S3 methods
# ══════════════════════════════════════════════════════════════════════════════

#' @export
print.remtribute <- function(x, ...) {
  cat("Relational Event Attribute Model\n")
  cat("\n")
  cat("  Attribute:     ", x$attribute, "\n")
  cat("  Type:          ", x$attribute_type, "\n")
  cat("  Events:        ", x$n_events, "\n")
  if (!is.null(x$n_levels))
    cat("  Levels:        ", x$n_levels,
        paste0("(", paste(x$levels, collapse = ", "), ")"), "\n")
  cat("  Statistics:    ", paste(x$stat_names, collapse = ", "), "\n")
  cat("  Log-likelihood:", round(x$loglik, 3), "\n")
  cat("  AIC:           ", round(x$AIC, 3), "\n")
  cat("  BIC:           ", round(x$BIC, 3), "\n")
  cat("\nCoefficients:\n")
  print(x$coefficients)
  invisible(x)
}

#' @export
summary.remtribute <- function(object, ...) {
  cat("Relational Event Outcome Model\n")
  cat("\n")
  cat("  Attribute:     ", object$attribute, "\n")
  cat("  Type:          ", object$attribute_type, "\n")
  cat("  Events:        ", object$n_events, "\n")
  if (!is.null(object$n_levels))
    cat("  Levels:        ", object$n_levels,
        paste0("(", paste(object$levels, collapse = ", "), ")"), "\n")
  cat("  Statistics:    ", paste(object$stat_names, collapse = ", "), "\n\n")

  cat("--- Backend model summary ---\n\n")
  summary(object$fit, ...)
}

#' @export
coef.remtribute <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.remtribute <- function(object, ...) {
  object$vcov
}

#' @export
logLik.remtribute <- function(object, ...) {
  object$loglik
}

# NULL-coalescing (if not already available in package namespace)
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}
