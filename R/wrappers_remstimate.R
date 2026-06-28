.get_durem_or_standard_stat_names <- function(stats) {
  if (inherits(stats, "remstats_durem")) {
    c(if (!is.null(stats$start_stats)) dimnames(stats$start_stats)[[3]],
      if (!is.null(stats$end_stats)) dimnames(stats$end_stats)[[3]])
  } else if (inherits(stats, "aomstats")) {
    c(dimnames(stats$sender_stats)[[3]], dimnames(stats$receiver_stats)[[3]])
  } else {
    dimnames(stats)[[3]]
  }
}

#' Dyadic Latent Class REM
#'
#' Fits a mixture REM where dyads are assigned to K latent classes,
#' each with class-specific coefficients.
#'
#' @param reh A \code{remify} object.
#' @param stats A \code{remstats} object.
#' @param k Number of latent classes (default 2).
#' @param nrep Number of random restarts (default 3).
#' @param ... Additional arguments passed to \code{remstimate}.
#' @return A \code{remstimate_mixrem} object.
#' @references Lakdawala, Leenders, and Mulder (2025). Not all bonds are created
#' equal: Dyadic latent class models for relational event data. \emph{arXiv preprint (2501.04418)}.
#' @examples
#' \donttest{
#' # Tie-oriented frailty
#' # EXAMPLES
#' }
#'
#' @export
dlcrem <- function(reh, stats, k = 2L, nrep = 3L, ...) {
  random <- stats::as.formula(
    paste0("~ (1 | dyad)")
  )
  remstimate(reh, stats, mixture = list(k = k, random = random), nrep = nrep, ...)
}


#' Frailty REM
#'
#' Fits a relational event model with actor-level random intercepts (frailty)
#' to account for unobserved heterogeneity in actor activity and
#' attractiveness. This is a convenience wrapper around
#' \code{\link{remstimate}} with \code{method = "GLMM"} and a default
#' random-effects formula.
#'
#' For \strong{tie-oriented models}, each actor receives a random intercept
#' for sending (\code{actor1}) and receiving (\code{actor2}), i.e.
#' \code{(1 | actor1) + (1 | actor2)}. In directed
#' models these capture sender activity and receiver attractiveness; in
#' undirected models both reflect general sociability.
#'
#' For \strong{actor-oriented models}, the sender rate model receives
#' \code{(1 | actor1)} and the receiver choice model receives
#' \code{(1 | actor2)}.
#'
#' For \strong{duration models} (\code{remify_durem}), random intercepts
#' are crossed with the \code{process} indicator (start vs end), giving
#' each actor separate frailty terms for event initiation and termination.
#'
#' When \code{extend_riskset_by_type = TRUE} (only possible for a tie-oriented
#' model), random intercepts are further
#' crossed with the event \code{type}, allowing actor heterogeneity to vary
#' across event types.
#'
#' For interval timing, estimation uses \code{lme4} or \code{glmmTMB}
#' (Poisson GLMM). For ordinal timing or actor-oriented receiver choice
#' models, \code{coxme} is used automatically (conditional logit with
#' random effects).
#'
#' Users who need non-default random-effects structures (e.g., random
#' slopes, dyad-level intercepts) should call
#' \code{\link{remstimate}(..., method = "GLMM", random = ...)} directly.
#'
#' @param reh A \code{remify} or \code{remify_durem} object.
#' @param stats A \code{remstats} object (\code{tomstats}, \code{aomstats},
#'   or \code{remstats_durem}).
#' @param approach Either \code{frequentist} or \code{Bayesian} estimation,
#' is supported.
#' @param engine For \code{frequentist}, the backend for interval models:
#' \code{"lme4"} (default) or
#'   \code{"glmmTMB"}. Ignored for ordinal models, which always use
#'   \code{coxme}. Bayesian frailty is not offered this version
#' @param ... Additional arguments passed to \code{\link{remstimate}}.
#'
#' @return A \code{remstimate_glmm} object. See \code{\link{remstimate}}
#'   for details on the return structure.
#'
#' @seealso \code{\link{remstimate}} for the general estimation interface,
#'   \code{\link{dlcrem}} for dyadic latent class models.
#'
#' @examples
#' \donttest{
#' # Tie-oriented frailty
#' reh <- remify::remify(tie_data$edgelist, model = "tie")
#' stats <- remstats::remstats(reh, tie_effects = ~ inertia() + reciprocity())
#' fit <- frailty_rem(reh, stats)
#' summary(fit)
#' lme4::ranef(fit$backend_fit)
#'
#' # Actor-oriented frailty
#' reh_ao <- remify::remify(tie_data$edgelist, model = "actor")
#' stats_ao <- remstats::remstats(reh_ao,
#'   sender_effects = ~ 1 + indegreeSender(),
#'   receiver_effects = ~ inertia())
#' fit_ao <- frailty_rem(reh_ao, stats_ao)
#'
#' # Duration model frailty
#' reh_dur <- remify::remify(edgelist, duration = TRUE)
#' stats_dur <- remstats::remstats(reh_dur,
#'   start_effects = ~ inertia(), end_effects = ~ inertia())
#' fit_dur <- frailty_rem(reh_dur, stats_dur)
#' }
#'
#' @export
frailty_rem <- function(reh, stats, approach = c("frequentist","Bayesian"),
                        engine = "lme4", ...) {
  is_durem    <- inherits(reh, "remify_durem")
  is_actor    <- inherits(stats, "aomstats")
  ext_by_type <- isTRUE(reh$meta$with_type_riskset)

  # Split start/end random intercepts only when BOTH processes are in the
  # stack. A start-only or end-only model has a single process level, so
  # `:process` is degenerate (collapses to actorN).
  split_process <- FALSE
  if (is_durem) {
    proc <- stats$stacked$remstats_stack$process
    split_process <- !is.null(proc) && length(unique(proc)) > 1L
  }

  if (is_actor) {
    random <- list(sender = ~ (1 | actor), receiver = ~ (1 | actor))
  } else {
    proc_tag <- if (is_durem && split_process) ":process" else ""
    type_tag <- if (ext_by_type) ":type" else ""
    random <- stats::as.formula(sprintf(
      "~ (1 | actor1%s%s) + (1 | actor2%s%s)",
      proc_tag, type_tag, proc_tag, type_tag
    ))
  }

  remstimate(reh, stats, approach = approach, random = random,
             engine = engine, ...)
}


