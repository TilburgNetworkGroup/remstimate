# ─────────────────────────────────────────────────────────────────────────────
# remstimate_durem.R
# Duration Relational Event Model (DuREM) — remstimate S3 method
#
# Dispatched from remstimate(reh, stats) when stats is a remstats_durem object.
# Fits the joint start/end Poisson GLM by:
#   1. Stacking statistics via stack_stats.remstats_durem
#   2. Fitting glm(obs ~ -1 + offset(log_interevent) + <stat_names>, poisson)
#   3. Returning a remstimate_durem / remstimate S3 object
# ─────────────────────────────────────────────────────────────────────────────


#' Fit a Duration Relational Event Model via Poisson GLM
#'
#' @description
#' S3 method for \code{\link{remstimate}} dispatched when \code{stats} is a
#' \code{remstats_durem} object.  Calls
#' \code{\link{stack_stats.remstats_durem}} to build the stacked long-format
#' data frame and fits the joint start / end model as a Poisson GLM with a
#' log inter-event-time offset.
#'
#' @param reh    A \code{remify_durem} object (output of
#'   \code{remify(edgelist, duration = TRUE)}).
#' @param stats  A \code{remstats_durem} object (output of
#'   \code{remstats(reh, start_effects = ..., end_effects = ...)}).
#' @param method Ignored; DuREM always uses Poisson MLE.  Accepted for
#'   compatibility with the \code{remstimate} generic signature.
#' @param ...    Additional arguments (currently unused).
#'
#' @return An S3 object of class \code{c("remstimate_durem", "remstimate")}
#'   with elements:
#'   \describe{
#'     \item{\code{coefficients}}{Named numeric vector of fitted parameters.}
#'     \item{\code{loglik}}{Scalar log-likelihood.}
#'     \item{\code{formula}}{The \code{glm} formula used.}
#'     \item{\code{stacked_data}}{The \code{remstats_stacked_durem} list from
#'       \code{stack_stats}.}
#'     \item{\code{backend_fit}}{The \code{glm} fit object.}
#'   }
#'   Attributes: \code{model = "durem"}, \code{method = "MLE"},
#'   \code{engine = "glm"}, \code{ordinal = FALSE},
#'   \code{statistics = <stat_names>}.
#'
#' @seealso \code{\link{remstimate}}, \code{\link{stack_stats.remstats_durem}},
#'   \code{\link{duremstimate}}
#'
#' @export
#' @method remstimate remstats_durem
remstimate.remstats_durem <- function(reh,
                                      stats,
                                      method = "MLE",
                                      ...) {

  if (!inherits(reh, "remify_durem"))
    stop("'reh' must be a remify_durem object (from remify(edgelist, duration = TRUE)).")

  if (!inherits(stats, "remstats_durem"))
    stop("'stats' must be a remstats_durem object.")

  # ── 1. Stack statistics ────────────────────────────────────────────────────
  stacked    <- stack_stats(stats, reh, add_actors = FALSE)
  df         <- stacked$remstats_stack
  stat_names <- stacked$stat_names

  if (length(stat_names) == 0L)
    stop("No statistics found in stats object — check start_effects / end_effects.")

  # ── 2. Build Poisson GLM formula ──────────────────────────────────────────
  # offset = log_interevent (log inter-event time at each unique time point)
  formula_obj <- as.formula(paste0(
    "obs ~ -1 + offset(log_interevent) + ",
    paste(stat_names, collapse = " + ")
  ))

  # ── 3. Fit Poisson GLM ────────────────────────────────────────────────────
  fit        <- glm(formula_obj, family = poisson(), data = df)
  coefs      <- coef(fit)
  loglik_val <- as.numeric(logLik(fit))

  # ── 4. Return structured remstimate_durem object ──────────────────────────
  structure(
    list(
      coefficients = coefs,
      loglik       = loglik_val,
      formula      = formula_obj,
      stacked_data = stacked,
      backend_fit  = fit
    ),
    class      = c("remstimate_durem", "remstimate"),
    model      = "durem",
    approach   = "Frequentist",
    method     = "MLE",
    engine     = "glm",
    ordinal    = FALSE,
    statistics = stat_names
  )
}


# ── coef / logLik / print methods ─────────────────────────────────────────────

#' @export
#' @method coef remstimate_durem
coef.remstimate_durem <- function(object, ...) object$coefficients

#' @export
#' @method logLik remstimate_durem
logLik.remstimate_durem <- function(object, ...) {
  structure(object$loglik,
            df    = length(object$coefficients),
            class = "logLik")
}

#' @export
#' @method print remstimate_durem
print.remstimate_durem <- function(x, ...) {
  cat("Duration Relational Event Model (DuREM) — Poisson MLE\n")
  cat("------------------------------------------------------\n")
  cat("Coefficients:\n")
  print(x$coefficients)
  cat("\nLog-likelihood:", x$loglik, "\n")
  invisible(x)
}

#' @export
#' @method summary remstimate_durem
summary.remstimate_durem <- function(object, ...) {
  summary(object$backend_fit, ...)
}
