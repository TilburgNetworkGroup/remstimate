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

#' Finite-mixture REM (latent classes)
#'
#' Fits a finite-mixture relational event model: events are assigned to
#' \code{k} latent classes, each with class-specific coefficients, via
#' \pkg{flexmix}. This is the general mixture interface; \code{\link{dlcrem}}
#' is the dyadic-latent-class special case.
#'
#' @param reh A \code{remify} object.
#' @param stats A \code{remstats} object.
#' @param random Clustering formula, e.g. \code{~ (1 | dyad)} or
#'   \code{~ (1 + inertia | dyad)}.
#' @param k Number of latent classes (default 2).
#' @param concomitant Optional concomitant formula for class probabilities.
#' @param nrep Number of random restarts (default 3).
#' @param ... Additional arguments passed to \code{\link{remstimate}}.
#' @return A \code{remstimate_mixrem} object.
#' @seealso \code{\link{dlcrem}}, \code{\link{remstimate}}.

#' @references
#' Lakdawala, R., Leenders, R., & Mulder, J. (2026). Not all bonds are created
#' equal: Dyadic latent class models for relational event data.
#' \emph{Social Networks}. \doi{10.1016/j.socnet.2026.06.006}
#'
#' @examples
#' \donttest{
#' library(remstimate)
#' data(tie_data)
#' reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
#' stats <- remstats::remstats(reh = reh,
#'   tie_effects = ~ 1 + remstats::inertia() + remstats::reciprocity())
#'
#' # two-component mixture, dyads as clustering unit
#' fit2 <- remstimate(reh, stats,
#'                    mixture = list(random = ~ (1 + inertia | dyad), k = 2))
#' fit2
#'
#' # coefficients: rows = effects, columns = components
#' fit2$coefficients
#' fit2$prior_probs
#'
#' # which dyads ended up in which component?
#' flexmix::clusters(fit2$backend_fit)
#'
#' # try k = 2:4 and compare via BIC
#' fits <- remstimate(reh, stats,
#'                    mixture = list(random = ~ (1 + inertia | dyad), k = 2:4))
#' bic_table(fits)
#' plot(fits)   # prints BIC table
#'
#' # pick the best k and inspect
#' best <- fits[[ paste0("k", bic_table(fits)$k[1]) ]]
#' best
#'
#' # diagnostics: per-component recall + combined
#' diag2 <- diagnostics(fit2, reh, stats)
#' diag2
#' plot(fit2, reh, stats = stats)
#' }
#'
#' @export
remixture <- function(reh, stats, random, k = 2L, concomitant = NULL,
                      nrep = 3L, ...) {
  remstimate(reh, stats,
             mixture = list(k = k, random = random,
                            concomitant = concomitant, nrep = nrep), ...)
}

#' Dyadic Latent Class REM
#'
#' Fits a mixture REM where dyads are assigned to \code{k} latent classes,
#' each with class-specific coefficients. Special case of
#' \code{\link{remixture}} with a dyad-level random intercept
#' (\code{~ (1 | dyad)}).
#'
#' @param reh A \code{remify} object.
#' @param stats A \code{remstats} object.
#' @param k Number of latent classes (default 2).
#' @param nrep Number of random restarts (default 3).
#' @param ... Additional arguments passed to \code{\link{remstimate}}.
#' @return A \code{remstimate_mixrem} object.
#' @references
#' Lakdawala, R., Leenders, R., & Mulder, J. (2026). Not all bonds are created
#' equal: Dyadic latent class models for relational event data.
#' \emph{Social Networks}. \doi{10.1016/j.socnet.2026.06.006}
#' @seealso \code{\link{remixture}} for the general mixture interface.
#'
#' @examples
#' data(tie_data)
#' reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
#' stats <- remstats::remstats(reh = reh,
#'   tie_effects = ~ 1 + remstats::inertia() + remstats::reciprocity())
#' fit_dlcrem <- dlcrem(reh, stats)
#'
#' @export
dlcrem <- function(reh, stats, k = 2L, nrep = 3L, ...) {
  remixture(reh, stats, random = ~ (1 | dyad), k = k, nrep = nrep, ...)
}


#' Frailty REM (actor / dyad random intercepts)
#'
#' Convenience wrapper around \code{\link{remstimate}} that fits a GLMM with the
#' default \emph{frailty} structure: random intercepts capturing unobserved
#' heterogeneity in actor activity and attractiveness. \code{remfrailty} only
#' builds this default structure; for custom random-effects (random slopes,
#' dyad-level intercepts, etc.) call \code{\link{remstimate}(..., random = ...)}
#' directly.
#'
#' For \strong{directed tie-oriented models}, each actor receives a random
#' intercept for sending (\code{actor1}) and receiving (\code{actor2}), i.e.
#' \code{(1 | actor1) + (1 | actor2)}, capturing sender activity and receiver
#' attractiveness. Actor-level frailty is \strong{not} identified for undirected
#' tie models (\code{actor1}/\code{actor2} are an arbitrary dyad ordering, not
#' sender/receiver), so for undirected data \code{remfrailty()} falls back to
#' symmetric dyad-level frailty \code{(1 | dyad)} and emits a message. (An
#' explicit \code{remstimate(..., random = ~ (1 | actor1) + (1 | actor2))} on
#' undirected data still errors.)
#'
#' For \strong{actor-oriented models}, the sender rate model receives
#' \code{(1 | actor1)} and the receiver choice model receives
#' \code{(1 | actor2)}.
#'
#' For \strong{duration models} (\code{remify_durem}), random intercepts are
#' crossed with the \code{process} indicator (start vs end). When
#' \code{extend_riskset_by_type = TRUE} (tie-oriented only), they are further
#' crossed with the event \code{type}.
#'
#' For interval timing, estimation uses \code{lme4} or \code{glmmTMB}
#' (Poisson GLMM). For ordinal timing or actor-oriented receiver choice
#' models, \code{coxme} is used automatically.
#'
#' @param reh A \code{remify} or \code{remify_durem} object.
#' @param stats A \code{remstats} object (\code{tomstats}, \code{aomstats},
#'   or \code{remstats_durem}).
#' @param approach Either \code{"frequentist"} or \code{"Bayesian"}. Bayesian
#'   frailty is not offered in this version.
#' @param engine For \code{frequentist}, the interval-model backend:
#'   \code{"glmmTMB"} or \code{"lme4"}. The default \code{"auto"} uses
#'   \code{glmmTMB} when installed, otherwise \code{lme4}. Ignored for ordinal
#'   models, which always use \code{coxme}.
#' @param ... Additional arguments passed to \code{\link{remstimate}}.
#'
#' @return A \code{remstimate_glmm} object. See \code{\link{remstimate}}
#'   for details on the return structure.
#'
#' @seealso \code{\link{remstimate}} for the general estimation interface (and
#'   for custom random-effects structures via \code{random = ...}),
#'   \code{\link{remixture}} / \code{\link{dlcrem}} for mixture models.
#'
#' @references
#' Juozaitiene, R., & Wit, E. C. (2024). Nodal heterogeneity can induce ghost
#' triadic effects in relational event models. \emph{Psychometrika}, 89(1),
#' 151-171. \doi{10.1007/s11336-024-09952-x}
#'
#' Mulder, J., & Hoff, P. D. (2024). A latent variable approach for modeling
#' relational data with multiple receivers. \emph{Annals of Applied
#' Statistics}. \doi{10.1214/24-AOAS1885}
#'
#' @examples
#' # GLMM fits stack the full case-control design and call lme4, which is slow
#' # on realistic data, so the examples are shown (via \dontrun) not executed.
#' \dontrun{
#' # Tie-oriented frailty
#' \donttest{
#' library(remstimate)
#'
#' data(tie_data)
#' reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
#' stats <- remstats::remstats(reh = reh,
#'   tie_effects = ~ 1 + remstats::inertia() + remstats::reciprocity())
#'
#' # sender frailty
#' fit_s <- remstimate(reh, stats, random = ~ (1 | actor1))
#' fit_s
#'
#' # sender + receiver frailty
#' fit_sr <- remstimate(reh, stats,
#'                      random = ~ (1 | actor1) + (1 | actor2))
#' fit_sr
#'
#' # inspect random effects (accessor matches the engine used, glmmTMB by default)
#' glmmTMB::ranef(fit_sr$backend_fit)
#' glmmTMB::VarCorr(fit_sr$backend_fit)
#'
#' # diagnostics and recall
#' # recall uses the random effects (BLUPs); Schoenfeld residuals and the
#' # waiting-time Q-Q are on the fixed effects, as for MLE.
#' diag_sr <- diagnostics(fit_sr, reh, stats)
#' diag_sr
#' plot(diag_sr, which = 1:3)
#'
#' # plot 7: random-effects normality Q-Q
#' plot(fit_sr, reh, stats = stats, which = 7)
#'
#' # glmmTMB engine (faster for large data, supports more complex structures)
#' fit_tmb <- remstimate(reh, stats,
#'                       random = ~ (1 | actor1), engine = "glmmTMB")
#' }
#'
#' # Actor-oriented frailty
#' reh_ao <- remify::remify(tie_data$edgelist, model = "actor")
#' stats_ao <- remstats::remstats(reh_ao,
#'   sender_effects = ~ 1 + indegreeSender(),
#'   receiver_effects = ~ inertia())
#' fit_ao <- remfrailty(reh_ao, stats_ao)
#' summary(fit_ao)
#' }
#'
#' @export
remfrailty <- function(reh, stats, approach = c("frequentist", "Bayesian"),
                       engine = "auto", ...) {
  approach    <- match.arg(approach)
  is_durem    <- inherits(reh, "remify_durem")
  is_actor    <- inherits(stats, "aomstats")
  ext_by_type <- isTRUE(reh$meta$with_type_riskset)

  # Actor-level (sender/receiver) frailty is only identified for directed tie
  # models; for undirected data actor1/actor2 are an arbitrary dyad ordering, so
  # remfrailty falls back to symmetric dyad-level frailty (constructed below).
  directed <- if (!is.null(reh$meta)) isTRUE(reh$meta$directed)
              else isTRUE(attr(reh, "directed"))

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
    if (directed) {
      random <- stats::as.formula(sprintf(
        "~ (1 | actor1%s%s) + (1 | actor2%s%s)",
        proc_tag, type_tag, proc_tag, type_tag
      ))
    } else {
      # Undirected: actor-level (sender/receiver) frailty is not identified;
      # default to symmetric dyad-level frailty.
      random <- stats::as.formula(sprintf("~ (1 | dyad%s%s)", proc_tag, type_tag))
      message("Undirected tie model: actor-level frailty is not identified; ",
              "using dyad-level frailty '~ (1 | dyad", proc_tag, type_tag, ")'.")
    }
  }

  remstimate(reh, stats, approach = approach, random = random,
             engine = engine, ...)
}

#' @rdname remfrailty
#' @export
frailty_rem <- function(reh, stats, ...) {
  .Deprecated("remfrailty")
  remfrailty(reh, stats, ...)
}


#' Penalized REM (lasso / ridge / elastic net; horseshoe when Bayesian)
#'
#' Convenience wrapper around \code{\link{remstimate}} for regularised
#' estimation. Frequentist fits use \pkg{glmnet}; Bayesian fits use
#' \pkg{shrinkem} shrinkage priors. The baseline / intercept is never
#' penalised (glmnet keeps it as its unpenalised intercept; shrinkem flags it
#' via \code{unpenalized}).
#'
#' @param reh A \code{remify} or \code{remify_durem} object.
#' @param stats A \code{remstats} object.
#' @param approach \code{"frequentist"} (glmnet) or \code{"Bayesian"}
#'   (shrinkem).
#' @param alpha Elastic-net mixing for the frequentist fit: \code{1} = lasso
#'   (default), \code{0} = ridge.
#' @param prior Shrinkage prior for the Bayesian fit (default
#'   \code{"horseshoe"}).
#' @param nfolds Cross-validation folds (frequentist). Default \code{10}.
#' @param lambda_select Which lambda (frequentist): \code{"1se"} (default) or
#'   \code{"min"}.
#' @param ... Additional arguments passed to \code{\link{remstimate}}.
#'
#' @return A \code{remstimate_glmnet} (frequentist) or
#'   \code{remstimate_shrinkem} (Bayesian) object.
#'
#' @seealso \code{\link{remstimate}}.
#'
#' @references
#' Karimova, D., Leenders, R., Meijerink-Bosman, M., & Mulder, J. (2023).
#' Separating the wheat from the chaff: Bayesian regularization in dynamic
#' social networks. \emph{Social Networks}, 74, 139-155.
#' \doi{10.1016/j.socnet.2023.02.006}
#'
#' Karimova, D., van Erp, S., Leenders, R., & Mulder, J. (2025). Honey, I
#' shrunk the irrelevant effects! Simple and flexible approximate Bayesian
#' regularization. \emph{Journal of Mathematical Psychology}, 126, 102925.
#' \doi{10.1016/j.jmp.2025.102925}
#'
#' Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
#' \emph{Journal of the Royal Statistical Society: Series B (Methodological)},
#' 58(1), 267-288. \doi{10.1111/j.2517-6161.1996.tb02080.x}
#'
#' @examples
#' # ---- MLE for basic tie model ----
#' data(tie_data)
#' reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
#' stats <- remstats::remstats(reh = reh,
#'   tie_effects = ~ 1 + remstats::inertia() + remstats::reciprocity())
#'
#' # ---- Penalised (lasso) ----
#' fit_lasso <- remstimate(reh, stats, penalty = list(alpha = 1))
#'
#' # ---- Penalised using Bayesian horseshoe (shrinkem) ----
#' fit_horseshoe <- remstimate(reh, stats, penalty = list(prior = "horseshoe"))
#'
#' @export
rempenalty <- function(reh, stats, approach = c("frequentist", "Bayesian"),
                       alpha = 1, prior = "horseshoe",
                       nfolds = 10L, lambda_select = c("1se", "min"), ...) {
  approach      <- match.arg(approach)
  lambda_select <- match.arg(lambda_select)

  penalty <- if (approach == "Bayesian")
    list(prior = prior)
  else
    list(alpha = alpha, nfolds = nfolds, lambda_select = lambda_select)

  remstimate(reh, stats, approach = approach, penalty = penalty, ...)
}





# ══════════════════════════════════════════════════════════════════════════
# Moving-window REM estimation
# ══════════════════════════════════════════════════════════════════════════

.subset_remstats_tie <- function(stats, a0, a1) {
  orig_subset <- attr(stats, "subset")
  offset <- orig_subset$start - 1L
  stopifnot(orig_subset$stop - orig_subset$start + 1L == dim(stats)[1])

  attrs <- attributes(stats)
  arr   <- stats[a0:a1, , , drop = FALSE]

  attrs$dim      <- dim(arr)
  attrs$dimnames <- dimnames(arr)
  attrs$subset   <- data.frame(start = a0 + offset, stop = a1 + offset)

  if (inherits(stats, "tomstats_sampled")) {
    idx <- a0:a1
    attrs$sample_map    <- attr(stats, "sample_map")[idx, , drop = FALSE]
    attrs$samp_prob     <- attr(stats, "samp_prob")[idx, , drop = FALSE]
    attrs$log_samp_prob <- attr(stats, "log_samp_prob")[idx, , drop = FALSE]
    attrs$case_pos      <- attr(stats, "case_pos")[idx]
    attrs$events_mult   <- attr(stats, "events_mult")[idx]
  }

  attributes(arr) <- attrs
  arr
}

.subset_remstats_actor <- function(stats, a0, a1) {
  orig_subset <- attr(stats, "subset")
  offset <- orig_subset$start - 1L
  ref_dim <- if (!is.null(stats$sender_stats)) dim(stats$sender_stats)[1] else dim(stats$receiver_stats)[1]
  stopifnot(orig_subset$stop - orig_subset$start + 1L == ref_dim)

  out <- stats
  if (!is.null(out$sender_stats))
    out$sender_stats   <- out$sender_stats[a0:a1, , , drop = FALSE]
  if (!is.null(out$receiver_stats))
    out$receiver_stats <- out$receiver_stats[a0:a1, , , drop = FALSE]
  attr(out, "subset") <- data.frame(start = a0 + offset, stop = a1 + offset)
  out
}

#' Moving-window REM estimation
#'
#' Fits a relational event model repeatedly over slices of events, reusing a
#' single \code{remify} object and full \code{remstats} object throughout.
#' Each window's design is a row-slice of the stats array/list carrying an
#' updated \code{subset} attribute; \code{remstimate} uses that attribute to
#' pull the matching interevent times / dyad or actor IDs from the (unsliced)
#' \code{reh}.
#'
#' In \strong{auto mode} (\code{window.width = NULL}, the default),
#' \code{n.windows} equal contiguous slices spanning the full event range are
#' used, reduced automatically (with a message) if that would put any window
#' below the required minimum width. In \strong{manual mode}
#' (\code{window.width} supplied), windows have a fixed size and advance by
#' \code{step.size.window} (default equal to \code{window.width}, i.e.
#' non-overlapping); the final window absorbs any remainder so no events are
#' dropped.
#'
#' @param reh A \code{remify} object.
#' @param stats A \code{tomstats}, \code{tomstats_sampled}, or \code{aomstats}
#'   object built on the full \code{reh}.
#' @param n.windows Target number of windows in auto mode. Default 5.
#' @param window.width Fixed number of events per window (manual mode).
#' @param step.size.window Events to advance between windows (manual mode
#'   only). Default equals \code{window.width}.
#' @param start.point First array row to start windowing from. Default 1L.
#' @param min.events Floor on window width (auto and manual).
#' @param approach \code{"frequentist"} (default) or \code{"Bayesian"},
#'   forwarded to \code{remstimate}.
#' @param parallel Fit windows in parallel? Default \code{FALSE}.
#' @param ncores.window Parallel workers across windows (Unix/macOS only;
#'   falls back to sequential with a message on Windows).
#' @param ... Further arguments passed to \code{remstimate} for every window.
#'
#' @return A \code{remstimate_window} object: \code{list(fits, windows, call,
#'   type, mode, n.windows)}.
#'
#' @references
#' Meijerink-Bosman, M., Leenders, R., & Mulder, J. (2022). Dynamic relational
#' event modeling: Testing, exploring, and applying. PLoS One, 17(8).
#' \doi{10.1371/journal.pone.0272309}
#'
#' Mulder, J., & Leenders, R. T. A. (2019). Modeling the evolution of interaction
#' behavior in social networks: A dynamic relational event approach for real-time
#' analysis. Chaos, Solitons & Fractals, 119, 73-85.
#' \doi{10.1016/j.chaos.2018.11.027}
#'
#' @examples
#' # ---- MLE for basic tie model ----
#' data(tie_data)
#' reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
#' stats <- remstats::remstats(reh = reh,
#'   tie_effects = ~ 1 + remstats::inertia() + remstats::reciprocity())
#' fit <- remstimate(reh, stats)
#' summary(fit)
#'
#' \donttest{
#' # ---- MLE for a model of events with a duration (tie-oriented only) ----
#' # (the baboons dataset is provided by the 'remdata' package)
#' if (requireNamespace("remdata", quietly = TRUE)) {
#'   data(baboons_obs, package = "remdata")
#'   reh_dur <- remify::remify(baboons_obs$edgelist[1:1000,], model = "tie",
#'     directed = FALSE, duration = TRUE)
#'   remstats_dur <- remstats::remstats(reh_dur,
#'     start_effects = ~ inertia(scaling = "std") +
#'       activeDegreeDyad(scaling = "std"),
#'     end_effects = ~ totaldegreeDyad(scaling = "std"),
#'     first = 50)
#'   remstimate_dur <- remstimate(reh_dur, remstats_dur)
#'   summary(remstimate_dur)
#'   diagnos_dur <- diagnostics(remstimate_dur, reh_dur, remstats_dur)
#'   plot(diagnos_dur)
#' }
#'
#' # ---- Random effects (frequentist) ----
#' fit_glmm <- remstimate(reh, stats,
#'   random = ~ (1 | actor1) + (1 | actor2))
#'
#' # ---- Penalised (lasso) ----
#' fit_lasso <- remstimate(reh, stats, penalty = list(alpha = 1))
#'
#' # ---- Penalised using Bayesian horseshoe (shrinkem) ----
#' fit_horseshoe <- remstimate(reh, stats, penalty = list(prior = "horseshoe"))
#'
#' # ---- Mixture ----
#' fit_mix <- remstimate(reh, stats,
#'   mixture = list(k = 2, random = ~ (1 + inertia | dyad)))
#' }
#'

#' @export
remwindow <- function(reh, stats,
                      n.windows        = 5L,
                      window.width     = NULL,
                      step.size.window = NULL,
                      start.point      = 1L,
                      min.events       = 50L,
                      approach         = c("frequentist", "Bayesian"),
                      parallel         = FALSE,
                      ncores.window    = 1L,
                      ...) {

  approach <- match.arg(approach)
  is_tie   <- inherits(stats, c("tomstats", "tomstats_sampled"))
  is_actor <- inherits(stats, "aomstats")
  if (!is_tie && !is_actor)
    stop("remwindow() currently supports 'tomstats'/'tomstats_sampled' and ",
         "'aomstats' objects only.", call. = FALSE)

  auto_mode <- is.null(window.width)

  if (auto_mode && !is.null(step.size.window))
    stop("'step.size.window' requires 'window.width' to also be specified. ",
         "With window.width = NULL (auto mode), window sizing is controlled ",
         "by 'n.windows' instead.", call. = FALSE)
  if (!auto_mode && !missing(n.windows))
    warning("'window.width' was specified, so 'n.windows' is ignored.", call. = FALSE)

  M <- if (is_tie) dim(stats)[1] else {
    if (!is.null(stats$sender_stats)) dim(stats$sender_stats)[1]
    else if (!is.null(stats$receiver_stats)) dim(stats$receiver_stats)[1]
    else stop("aomstats object has neither 'sender_stats' nor 'receiver_stats'.", call. = FALSE)
  }
  P_floor <- if (is_tie) 20L * dim(stats)[3] else
    20L * max(dim(stats$sender_stats)[3] %||% 0L, dim(stats$receiver_stats)[3] %||% 0L)
  min_width <- max(min.events, P_floor)
  total     <- M - start.point + 1L

  stats_method <- attr(stats, "method") %||% "pt"
  orig_subset  <- attr(stats, "subset")
  offset       <- if (!is.null(orig_subset)) orig_subset$start - 1L else 0L
  time_axis    <- if (stats_method == "pt") unique(reh$edgelist$time) else reh$edgelist$time

  if (auto_mode) {
    nw <- max(1L, n.windows)
    while (nw > 1L && floor(total / nw) < min_width) nw <- nw - 1L
    if (floor(total / nw) < min_width)
      stop("Not enough events (", total, ") for even one window of the required ",
           "minimum width (", min_width, ").", call. = FALSE)
    if (nw < n.windows && !missing(n.windows))
      message("Reduced n.windows from ", n.windows, " to ", nw,
              " to keep each window >= ", min_width, " events.")
    edges  <- floor(seq(start.point - 1L, M, length.out = nw + 1L))
    starts <- edges[-length(edges)] + 1L
    ends   <- edges[-1L]
  } else {
    if (window.width < min_width)
      warning("'window.width' (", window.width, ") is below the recommended ",
              "EPV floor (", min_width, "); estimates may be unstable.", call. = FALSE)
    if (is.null(step.size.window)) step.size.window <- window.width
    starts <- seq(start.point, M, by = step.size.window)
    starts <- starts[starts <= M]
    ends   <- pmin(starts + window.width - 1L, M)
    ends[length(ends)] <- M
  }

  dots <- list(...)
  use_parallel <- parallel && ncores.window > 1L
  if (use_parallel && .Platform$OS.type == "windows") {
    message("Parallel windowing uses mclapply (fork-based); not available on ",
            "Windows. Falling back to sequential.")
    use_parallel <- FALSE
  }

  fit_one <- function(w) {
    s0 <- starts[w]; s1 <- ends[w]
    stats_w <- if (is_tie) .subset_remstats_tie(stats, s0, s1)
    else        .subset_remstats_actor(stats, s0, s1)
    fit <- tryCatch(remstimate(reh, stats_w, approach = approach, ...),
                    error = function(e) e)
    conv <- if (inherits(fit, "remstimate")) {
      if (is_tie) isTRUE(fit$converged)
      else all(vapply(fit, function(x) isTRUE(x$converged), logical(1)))
    } else NA
    list(fit = fit, info = data.frame(
      start_event = s0, end_event = s1,
      start_time  = time_axis[s0 + offset], end_time = time_axis[s1 + offset],
      n_events    = s1 - s0 + 1L, converged = conv))
  }

  results <- if (use_parallel) {
    parallel::mclapply(seq_along(starts), fit_one, mc.cores = ncores.window)
  } else lapply(seq_along(starts), fit_one)

  fits    <- lapply(results, `[[`, "fit")
  windows <- do.call(rbind, lapply(results, `[[`, "info"))
  windows <- cbind(window = seq_along(starts), windows)

  structure(
    list(fits = fits, windows = windows, call = match.call(),
         type = if (is_tie) "tie" else "actor",
         mode = if (auto_mode) "auto" else "manual",
         n.windows = length(starts)),
    class = c("remstimate_window", "list")
  )
}

#' @export
print.remstimate_window <- function(x, digits = 3, ...) {
  cf <- coef.remstimate_window(x)
  w  <- x$windows

  .print_block <- function(block, title) {
    cat(title, "\n")
    header <- rbind(
      time_range = paste0(round(w$start_time, 2), "-", round(w$end_time, 2)),
      n_events   = as.character(w$n_events),
      converged  = as.character(w$converged)
    )
    coefs <- t(round(block$coefficients, digits))
    out <- rbind(header, coefs)
    colnames(out) <- paste0("W", w$window)
    print(out, quote = FALSE, right = TRUE)
    cat("\n")
  }

  cat("Moving-window relational event model (", x$type, ")\n", sep = "")
  cat("windows =", nrow(w), " (mode:", x$mode, ")\n\n")

  if (x$type == "tie") {
    .print_block(cf, "Coefficients:")
  } else {
    if (!is.null(cf$sender))   .print_block(cf$sender,   "Sender (rate) model:")
    if (!is.null(cf$receiver)) .print_block(cf$receiver, "Receiver (choice) model:")
  }
  invisible(x)
}

#' @export
print.diagnostics.remstimate_window <- function(x, digits = 3, ...) {
  .print_recall <- function(recall, title) {
    if (is.null(recall)) return(invisible(NULL))
    rs <- recall$summary
    cat(title, "\n")
    cat(sprintf("  mean rank = %.3f  |  median rank = %.3f  |  prob ratio = %.2f  |  top %g%% = %.1f%%\n",
                rs$mean_rel_rank, rs$median_rel_rank, rs$mean_prob_ratio,
                rs$top_pct * 100, rs$top_pct_prop * 100))
  }

  cat("Moving-window diagnostics (", x$type, ", interpolated coefficients)\n", sep = "")
  cat("windows =", nrow(x$windows), "\n\n")

  if (x$type == "tie") {
    .print_recall(x$recall, "Recall:")
  } else {
    .print_recall(x$sender$recall,   "Recall - sender (rate) model:")
    .print_recall(x$receiver$recall, "Recall - receiver (choice) model:")
  }
  invisible(x)
}

.robust_ylim <- function(y, k = 4) {
  y <- y[is.finite(y)]
  if (length(y) < 2) return(range(y, na.rm = TRUE) + c(-0.1, 0.1))
  med <- stats::median(y); iqr <- stats::IQR(y)
  if (iqr == 0) {
    pad <- if (med == 0) 0.5 else abs(med) * 0.5
    return(c(med - pad, med + pad))
  }
  c(med - k * iqr, med + k * iqr)
}

#' @export
plot.remstimate_window <- function(x, ci = 0.95, max.panels = 6,
                                   robust.ylim = TRUE, k = 4, ...) {
  cf  <- coef.remstimate_window(x, ci = ci)
  mid <- rowMeans(x$windows[, c("start_time", "end_time")])

  .plot_trajectory_block <- function(block, title_prefix) {
    var_names <- colnames(block$coefficients)
    P <- length(var_names)
    pages <- split(seq_len(P), ceiling(seq_len(P) / max.panels))

    for (pg in seq_along(pages)) {
      idx <- pages[[pg]]
      n   <- length(idx)
      old_par <- par(mfrow = c(min(n, 2), ceiling(n / min(n, 2))), mar = c(4, 4, 2, 1))
      for (p in idx) {
        y <- block$coefficients[, p]; lo <- block$lower[, p]; hi <- block$upper[, p]
        ok <- !is.na(y)
        ylims <- if (robust.ylim) .robust_ylim(y[ok], k = k)
        else range(c(lo[ok], hi[ok]), na.rm = TRUE)
        plot(mid, y, type = "n", ylim = ylims,
             xlab = "time (window midpoint)", ylab = var_names[p],
             main = paste0(title_prefix, var_names[p]))
        polygon(c(mid[ok], rev(mid[ok])), c(lo[ok], rev(hi[ok])),
                col = grDevices::adjustcolor("grey70", 0.4), border = NA)
        lines(mid[ok], y[ok], col = "steelblue", lwd = 2)
        points(mid[ok], y[ok], pch = 16, col = "steelblue")
        abline(h = 0, lty = 2, col = "grey40")
      }
      par(old_par)
      if (pg < length(pages) && interactive())
        invisible(readline(sprintf("Press <Enter> for next page (%d/%d)...", pg + 1, length(pages))))
    }
  }

  if (x$type == "tie") {
    .plot_trajectory_block(cf, "")
  } else {
    if (!is.null(cf$sender)) .plot_trajectory_block(cf$sender, "sender: ")
    if (!is.null(cf$sender) && !is.null(cf$receiver) && interactive())
      invisible(readline("Press <Enter> for receiver model plots..."))
    if (!is.null(cf$receiver)) .plot_trajectory_block(cf$receiver, "receiver: ")
  }
  invisible(x)
}

.stability_table <- function(block) {
  cf  <- block$coefficients; se <- block$se; sep <- block$separated
  eff <- colnames(cf)
  cf_clean <- cf; cf_clean[sep] <- NA
  se_clean <- se; se_clean[sep] <- NA

  excluded <- vapply(eff, function(e) {
    w <- which(sep[, e])
    if (length(w) == 0L) "-" else paste0("W", w, collapse = ", ")
  }, character(1))

  out <- data.frame(
    effect   = eff,
    n        = colSums(!is.na(cf_clean)),
    mean     = apply(cf_clean, 2, mean, na.rm = TRUE),
    sd       = apply(cf_clean, 2, sd,   na.rm = TRUE),
    min      = apply(cf_clean, 2, min,  na.rm = TRUE),
    max      = apply(cf_clean, 2, max,  na.rm = TRUE),
    mean_se  = apply(se_clean, 2, mean, na.rm = TRUE)
  )
  out$sd_over_se <- out$sd / out$mean_se
  out$excluded   <- excluded
  rownames(out) <- NULL
  out
}

.flag_separated <- function(se_mat, k = 10) {
  apply(se_mat, 2, function(col) {
    med <- stats::median(col, na.rm = TRUE)
    if (!is.finite(med) || med == 0) return(rep(FALSE, length(col)))
    col > k * med
  })
}

#' @export
coef.remstimate_window <- function(object, ci = 0.95, k = 10, ...) {
  z  <- stats::qnorm(1 - (1 - ci) / 2)
  nW <- length(object$fits)

  .one_matrix <- function(fits, get_coef, get_se, get_converged) {
    var_names <- NULL
    for (f in fits) {
      if (!inherits(f, "remstimate")) next
      cf <- tryCatch(get_coef(f), error = function(e) NULL)
      if (!is.null(cf)) { var_names <- names(cf); break }
    }
    if (is.null(var_names))
      stop("No successful fits in this remstimate_window object.", call. = FALSE)
    P <- length(var_names)

    coef_mat <- matrix(NA_real_, nW, P, dimnames = list(NULL, var_names))
    se_mat   <- matrix(NA_real_, nW, P, dimnames = list(NULL, var_names))

    for (w in seq_len(nW)) {
      f <- fits[[w]]
      if (!inherits(f, "remstimate")) next
      if (!isTRUE(tryCatch(get_converged(f), error = function(e) FALSE))) next
      cf <- tryCatch(get_coef(f), error = function(e) NULL)
      se <- tryCatch(get_se(f),   error = function(e) NULL)
      if (is.null(cf)) next
      coef_mat[w, names(cf)] <- cf
      se_mat[w,   names(se)] <- se
    }

    sep_mat <- .flag_separated(se_mat, k = k)
    dimnames(sep_mat) <- dimnames(se_mat)

    list(coefficients = coef_mat, se = se_mat, separated = sep_mat,
         lower = coef_mat - z * se_mat, upper = coef_mat + z * se_mat)
  }

  if (object$type == "tie") {
    out <- .one_matrix(object$fits,
                       get_coef      = function(f) f$coefficients,
                       get_se        = function(f) f$se %||% sqrt(diag(f$vcov)),
                       get_converged = function(f) f$converged)
    out$windows <- object$windows
    return(structure(out, class = "coef.remstimate_window"))
  }

  has_sender   <- any(vapply(object$fits, function(f) inherits(f, "remstimate") && !is.null(f$sender_model),   logical(1)))
  has_receiver <- any(vapply(object$fits, function(f) inherits(f, "remstimate") && !is.null(f$receiver_model), logical(1)))

  sender <- if (has_sender)
    .one_matrix(object$fits,
                get_coef      = function(f) f$sender_model$coefficients,
                get_se        = function(f) f$sender_model$se %||% sqrt(diag(f$sender_model$vcov)),
                get_converged = function(f) f$sender_model$converged)
  else NULL

  receiver <- if (has_receiver)
    .one_matrix(object$fits,
                get_coef      = function(f) f$receiver_model$coefficients,
                get_se        = function(f) f$receiver_model$se %||% sqrt(diag(f$receiver_model$vcov)),
                get_converged = function(f) f$receiver_model$converged)
  else NULL

  structure(list(sender = sender, receiver = receiver, windows = object$windows),
            class = "coef.remstimate_window")
}

#' @export
summary.remstimate_window <- function(object, k = 10, ...) {
  is_tie <- object$type == "tie"
  nW <- length(object$fits)
  w  <- object$windows

  loglik <- AIC <- BIC <- rep(NA_real_, nW)
  for (i in seq_len(nW)) {
    f <- object$fits[[i]]
    if (!inherits(f, "remstimate")) next
    if (is_tie) {
      if (!isTRUE(f$converged)) next
      loglik[i] <- f$loglik; AIC[i] <- f$AIC; BIC[i] <- f$BIC
    } else {
      sm <- f$sender_model; rm <- f$receiver_model
      sm_ok <- is.null(sm) || isTRUE(sm$converged)
      rm_ok <- is.null(rm) || isTRUE(rm$converged)
      if (!sm_ok || !rm_ok) next
      loglik[i] <- sum(c(sm$loglik, rm$loglik))
      AIC[i]    <- sum(c(sm$AIC,    rm$AIC))
      BIC[i]    <- sum(c(sm$BIC,    rm$BIC))
    }
  }
  fit_table <- cbind(w, loglik = loglik, AIC = AIC, BIC = BIC)

  .stability_table <- function(block) {
    cf <- block$coefficients; se <- block$se; sep <- block$separated
    eff <- colnames(cf)
    cf_clean <- cf; cf_clean[sep] <- NA
    se_clean <- se; se_clean[sep] <- NA

    excluded <- vapply(eff, function(e) {
      wi <- which(sep[, e])
      if (length(wi) == 0L) "-" else paste0("W", wi, collapse = ", ")
    }, character(1))

    out <- data.frame(
      effect  = eff,
      n       = colSums(!is.na(cf_clean)),
      mean    = apply(cf_clean, 2, mean, na.rm = TRUE),
      sd      = apply(cf_clean, 2, sd,   na.rm = TRUE),
      min     = apply(cf_clean, 2, min,  na.rm = TRUE),
      max     = apply(cf_clean, 2, max,  na.rm = TRUE),
      mean_se = apply(se_clean, 2, mean, na.rm = TRUE)
    )
    out$sd_over_se <- out$sd / out$mean_se
    out$excluded   <- excluded
    rownames(out) <- NULL
    out
  }

  cf <- coef.remstimate_window(object, k = k)

  if (is_tie) {
    structure(list(fit = fit_table, stability = .stability_table(cf), type = "tie"),
              class = "summary.remstimate_window")
  } else {
    structure(list(fit = fit_table,
                   stability_sender   = if (!is.null(cf$sender))   .stability_table(cf$sender)   else NULL,
                   stability_receiver = if (!is.null(cf$receiver)) .stability_table(cf$receiver) else NULL,
                   type = "actor"),
              class = "summary.remstimate_window")
  }
}

#' @export
print.summary.remstimate_window <- function(x, digits = 3, ...) {
  cat("Moving-window relational event model (", x$type, ")\n", sep = "")
  cat("windows =", nrow(x$fit), "\n\n")

  cat("Fit:\n")
  ft <- x$fit
  ft$time_range <- paste0(round(ft$start_time, 2), "-", round(ft$end_time, 2))
  ft <- ft[, c("window", "time_range", "n_events", "converged", "loglik", "AIC", "BIC")]
  ft[c("loglik", "AIC", "BIC")] <- lapply(ft[c("loglik", "AIC", "BIC")], round, digits)
  print(ft, row.names = FALSE)
  cat("\n")

  .print_stability <- function(tab, title) {
    cat(title, "\n")
    num_cols <- setdiff(names(tab), c("effect", "n", "excluded"))
    tab[num_cols] <- lapply(tab[num_cols], round, digits)
    print(tab, row.names = FALSE)
    cat("\n")
  }

  if (x$type == "tie") {
    .print_stability(x$stability, "Stability (across converged, non-separated windows):")
  } else {
    if (!is.null(x$stability_sender))   .print_stability(x$stability_sender,   "Stability - sender (rate) model:")
    if (!is.null(x$stability_receiver)) .print_stability(x$stability_receiver, "Stability - receiver (choice) model:")
  }
  invisible(x)
}



#' Compare BIC across a list of MIXREM fits
#'
#' @usage bic_table(x, ...)
#' @param x   A \code{remstimate_mixrem_list} returned when \code{k} is a
#'   vector in \code{remstimate(..., method = "MIXREM")}.
#' @param ... Unused.
#'
#' @return A data frame with columns \code{k}, \code{BIC}, and
#'   \code{delta_BIC}, sorted by ascending BIC.
#'
#' @examples
#' \donttest{
#' library(remstimate)
#' data(history, package = "remstats")
#' data(info,    package = "remstats")
#' colnames(history)[colnames(history) == "setting"] <- "type"
#'
#' reh   <- remify::remify(edgelist = history[1:100, ], model = "tie",
#'                          riskset = "active")
#' stats <- remstats::tomstats(~ inertia(consider_type = FALSE), reh = reh,
#'                              attr_actors = info)
#'
#' fits <- remstimate(reh, stats,
#'                    mixture = list(random = ~ (1 + inertia | dyad), k = 2:5))
#' bic_table(fits)
#' }
#'
#' @name bic_table
#' @export
bic_table

