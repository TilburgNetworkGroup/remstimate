#' @title remstimate
#'
#' @description Fit a relational event model, using a \code{reh} object of class
#' \code{remify} and a \code{stats} object of class \code{remstats}
#'
#'   The model structure is determined by the arguments:
#'   \itemize{
#'     \item No structure arguments: basic REM (tie or actor, depending on
#'       \code{reh}/\code{stats}).
#'     \item \code{random}: mixed-effects model (GLMM).
#'     \item \code{penalty}: penalized relational event modeling (GLMNET or Bayesian via shrinkem).
#'     \item \code{mixture}: finite mixture model (MIXREM / flexmix).
#'   }
#'
#'   The estimation approach is set via \code{approach}:
#'   \itemize{
#'     \item \code{"frequentist"} (default): MLE, lme4/glmmTMB, glmnet, or
#'       flexmix depending on the structure.
#'     \item \code{"Bayesian"}: HMC via C++ (basic models) or shrinkem
#'       (penalized models).
#'   }
#'
#' @param reh       A \code{remify} or \code{remify_durem} object.
#' @param stats     A \code{tomstats}, \code{tomstats_sampled},
#'   \code{aomstats}, or \code{remstats_durem} object.
#' @param approach  \code{"frequentist"} (default) or \code{"Bayesian"}.
#' @param random    One-sided formula for random effects, e.g.
#'   \code{~ (1 | actor1) + (1 | actor2)}. Fits a GLMM (lme4/glmmTMB, or
#'   coxme for ordinal models).
#' @param penalty   List of penalisation settings. Frequentist uses glmnet;
#'   Bayesian uses shrinkem. Recognised elements:
#'   \code{alpha} (elastic-net mixing, \code{1} = lasso (default),
#'   \code{0} = ridge), \code{nfolds} (CV folds, default \code{10}),
#'   \code{lambda_select} (\code{"1se"} (default) or \code{"min"}),
#'   \code{unpenalized} / \code{penalized} (see below), and \code{prior}
#'   (shrinkem prior, default \code{"horseshoe"}).
#'   By default the intercept / baseline structure is left unpenalised: any
#'   statistic that is an indicator (all values in \{0,1\}) is exempt, which
#'   covers the overall \code{baseline}, the duration start/end process
#'   intercepts (\code{baseline.start} / \code{baseline.end}) and fixed-effect
#'   type dummies (e.g. \code{FEtype_*}); count and continuous effect statistics
#'   are penalised. Adjust this default with two additive controls (both
#'   character vectors of statistic names): \code{unpenalized} \emph{adds}
#'   names to the exemption, and \code{penalized} \emph{removes} names from it,
#'   i.e. forces them back into the penalty even though they are 0/1 indicators
#'   (e.g. a p-shift dummy such as \code{psABAB.end} that is a genuine effect,
#'   not an intercept). \code{penalized} takes precedence when a name is given
#'   in both. Names must match the model statistics exactly (as printed by the
#'   remstats object, e.g. \code{"psABAB.end"}, not \code{"psABAB"}); a name that
#'   matches nothing is ignored with a warning.
#' @param mixture   List of finite-mixture settings (flexmix). Recognised
#'   elements: \code{k} (components, default \code{2}), \code{random}
#'   (clustering formula, e.g. \code{~ (1 | dyad)}), \code{concomitant}
#'   (optional concomitant formula), \code{nrep} (random restarts, default
#'   \code{3}). E.g. \code{list(k = 2, random = ~ (1 | dyad))} fits the dyadic
#'   latent class model (Lakdawala et al., 2026).
#' @param engine    GLMM backend: \code{"glmmTMB"} or \code{"lme4"}; ordinal
#'   models automatically use \code{coxme}. The default \code{"auto"} selects
#'   \code{"glmmTMB"} when installed (more robust on the stacked tie-oriented
#'   design), otherwise falls back to \code{"lme4"}.
#' @param bayes     List of Bayesian (C++ HMC) controls for basic tie/actor
#'   models. Recognised elements: \code{nsim} (post-burnin iterations, default
#'   \code{2000}), \code{nchains} (default \code{2}), \code{burnin} (default
#'   \code{1000}), \code{thin} (default \code{1}), \code{init} (initial values,
#'   default MLE estimates), \code{L} (leapfrog steps, default \code{50}),
#'   \code{epsilon} (leapfrog step size, default \code{0.002}), \code{prior}
#'   (list with \code{mean} and \code{vcov}), and \code{nsimWAIC} (WAIC draws,
#'   default \code{100}).
#' @param seed      Random seed.
#' @param ncores    Number of threads (C++ backends). Default \code{1L}.
#' @param WAIC      Compute WAIC? Default \code{FALSE}.
#' @param method    [\emph{Deprecated}] Legacy argument for backward
#'   compatibility. Use \code{approach} instead. Accepts \code{"MLE"} or
#'   \code{"HMC"} for basic models.
#' @param ...       Further arguments passed to the backend. The former
#'   top-level knobs (\code{alpha}, \code{nfolds}, \code{lambda_select},
#'   \code{k}, \code{concomitant}, \code{nrep}, \code{nsim}, \code{nchains},
#'   \code{burnin}, \code{thin}, \code{init}, \code{L}, \code{epsilon},
#'   \code{prior}, \code{nsimWAIC}) are \emph{deprecated} but still accepted
#'   here and routed into \code{penalty} / \code{mixture} / \code{bayes}.
#'
#' @return A \code{remstimate} S3 object.
#'
#' @references
#' Butts, C. T. (2008). A relational event framework for social action.
#' \emph{Sociological Methodology}, 38(1), 155-200.
#' \doi{10.1111/j.1467-9531.2008.00203.x}
#'
#' Stadtfeld, C., & Block, P. (2017). Interactions, actors, and time: Dynamic
#' network actor models for relational events. \emph{Sociological Science}, 4,
#' 318-352. \doi{10.15195/v4.a14}
#'
#' Juozaitiene, R., & Wit, E. C. (2024). Nodal heterogeneity can induce ghost
#' triadic effects in relational event models. \emph{Psychometrika}, 89(1),
#' 151-171. \doi{10.1007/s11336-024-09952-x}
#'
#' Mulder, J., & Hoff, P. D. (2024). A latent variable approach for modeling
#' relational data with multiple receivers. \emph{Annals of Applied
#' Statistics}. \doi{10.1214/24-AOAS1885}
#'
#' Lakdawala, R., Leenders, R., & Mulder, J. (2026). Not all bonds are created
#' equal: Dyadic latent class models for relational event data.
#' \emph{Social Networks}. \doi{10.1016/j.socnet.2026.06.006}
#'
#' Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
#' \emph{Journal of the Royal Statistical Society: Series B (Methodological)},
#' 58(1), 267-288. \doi{10.1111/j.2517-6161.1996.tb02080.x}
#'
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
#' Lakdawala, R., Leenders, R., Ejbye-Ernst, P., & Mulder, J. (2026).
#' Modelling interaction duration in relational event models.
#' \emph{arXiv preprint}. \doi{10.48550/arXiv.2602.21000}
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
remstimate <- function(reh,
                       stats,
                       approach = c("frequentist", "Bayesian"),
                       # Model structure (each list owns its own knobs)
                       random   = NULL,   # GLMM  : one-sided formula
                       penalty  = NULL,   # glmnet/shrinkem: list(alpha, nfolds, lambda_select, prior)
                       mixture  = NULL,   # flexmix: list(k, random, concomitant, nrep)
                       engine   = "auto",
                       # Bayesian (C++ HMC) controls, e.g.
                       #   list(nsim, nchains, burnin, thin, init, L, epsilon, prior, nsimWAIC)
                       bayes    = list(),
                       # Shared
                       seed     = NULL,
                       ncores   = 1L,
                       WAIC     = FALSE,
                       method   = NULL,   # deprecated: "MLE" / "HMC"
                       ...) {

  # ── Deprecated top-level knobs (nsim, alpha, k, nfolds, …) are still accepted
  #    via ... and routed into the penalty/mixture lists or the 'bayes' bundle,
  #    so pre-existing code keeps working. Consumed names are stripped before
  #    the remaining ... are forwarded to the backends. ─────────────────────────
  dots <- list(...)
  .consumed <- c("alpha", "nfolds", "lambda_select", "unpenalized", "penalized",
                 "k", "concomitant", "nrep",
                 "nsim", "nchains", "burnin", "thin", "init", "L", "epsilon",
                 "prior", "nsimWAIC")
  extra <- dots[setdiff(names(dots), .consumed)]
  dep   <- function(nm) dots[[nm]]

  # ── Backward compatibility: method = "MLE" / "HMC" ─────────────────────────
  if (!is.null(method)) {
    method <- toupper(method)
    if (method == "MLE")      approach <- "frequentist"
    else if (method == "HMC") approach <- "Bayesian"
    else stop("'method' is deprecated. ",
              "Use 'approach' with 'random', 'penalty', or 'mixture'.",
              call. = FALSE)
  }

  approach <- match.arg(approach, choices = c("frequentist", "Bayesian"))
  is_durem <- inherits(reh, "remify_durem") && inherits(stats, "remstats_durem")
  is_sampled <- inherits(stats, "tomstats_sampled")

  # ── Validate structure arguments ───────────────────────────────────────────
  has_random  <- !is.null(random)
  has_penalty <- !is.null(penalty)
  has_mixture <- !is.null(mixture)

  if (has_mixture && (has_random || has_penalty))
    stop("'mixture' cannot be combined with 'random' or 'penalty'.", call. = FALSE)
  if (has_random && has_penalty)
    stop("Combining 'random' and 'penalty' is not supported.", call. = FALSE)

  # ── WAIC / nsimWAIC / ncores validation ────────────────────────────────────
  if (!is.logical(WAIC) || length(WAIC) != 1)
    stop("'WAIC' must be a single logical value (TRUE or FALSE).")
  nsimWAIC <- bayes$nsimWAIC %||% dep("nsimWAIC") %||% 100L
  if (!is.numeric(nsimWAIC) || length(nsimWAIC) != 1 || nsimWAIC < 1) {
    warning("'nsimWAIC' must be a positive integer. Using default value of 100.")
    nsimWAIC <- 100L
  }
  nsimWAIC <- as.integer(nsimWAIC)

  if (!is.numeric(ncores) || ncores < 1L) {
    warning("'ncores' must be a positive integer. Using default value of 1.")
    ncores <- 1L
  }
  ncores <- as.integer(ncores)

  # ── Mixture dispatch ───────────────────────────────────────────────────────
  if (has_mixture) {
    if(is_sampled){
      stop("sampled tomstats are not supported for mixture models")
    }
    return(do.call(.remstimate_mixrem, c(list(
      reh, stats,
      random      = mixture$random,
      k           = mixture$k %||% dep("k") %||% 2L,
      concomitant = mixture$concomitant %||% dep("concomitant"),
      nrep        = mixture$nrep %||% dep("nrep") %||% 3L), extra)))
  }

  # ── Penalised dispatch ─────────────────────────────────────────────────────
  #   Bayesian  -> shrinkem shrinkage prior (penalty$prior)
  #   Frequentist -> glmnet (penalty$alpha / nfolds / lambda_select)
  if (has_penalty) {
    if (approach == "Bayesian") {
      return(do.call(.remstimate_shrinkem, c(list(
        reh, stats,
        type        = penalty$prior %||% "horseshoe",
        unpenalized = penalty$unpenalized %||% dep("unpenalized"),
        penalized   = penalty$penalized   %||% dep("penalized"),
        ncores = ncores, seed = seed), extra)))
    }
    return(do.call(.remstimate_glmnet, c(list(
      reh, stats,
      alpha         = penalty$alpha %||% dep("alpha") %||% 1,
      nfolds        = penalty$nfolds %||% dep("nfolds") %||% 10L,
      lambda_select = penalty$lambda_select %||% dep("lambda_select") %||% "1se",
      unpenalized   = penalty$unpenalized %||% dep("unpenalized"),
      penalized     = penalty$penalized   %||% dep("penalized")),
      extra)))
  }

  # ── Random effects dispatch (GLMM) ─────────────────────────────────────────
  if (has_random) {
    if(is_sampled){
      stop("sampled tomstats are not supported for random effects models")
    }
    if (approach == "Bayesian")
      stop("Bayesian random-effects models are not available in this version. ",
           "Use approach = 'frequentist' (GLMM via lme4/glmmTMB).", call. = FALSE)
    eng <- if (identical(engine, "auto")) {
      if (requireNamespace("glmmTMB", quietly = TRUE)) "glmmTMB" else "lme4"
    } else engine
    return(do.call(.remstimate_glmm, c(list(
      reh, stats, random = random, engine = eng), extra)))
  }

  # ── Basic model: duration / pre-stacked → GLM pipeline ─────────────────────
  if (is_durem || inherits(stats, "remstats_stacked")) {
    if (approach == "Bayesian")
      stop("Bayesian estimation is not available for duration or pre-stacked ",
           "models. Use approach = 'frequentist'.", call. = FALSE)
    stacked <- if (inherits(stats, "remstats_stacked")) stats
    else stack_stats(stats, reh, add_actors = FALSE)
    if (isTRUE(stacked$model == "actor"))
      return(.remstimate_stacked_glm_actor(stacked))
    return(.remstimate_glm(stacked))   # tie / durem
  }

  # ── Basic tie/actor → C++ backends ─────────────────────────────────────────
  if (approach == "Bayesian") {
    return(do.call(.remstimate_hmc_cpp, c(list(
      reh = reh, stats = stats, ncores = ncores,
      nsim    = bayes$nsim    %||% dep("nsim")    %||% 2000L,
      nchains = bayes$nchains %||% dep("nchains") %||% 2L,
      burnin  = bayes$burnin  %||% dep("burnin")  %||% 1000L,
      thin    = bayes$thin    %||% dep("thin")    %||% 1L,
      init    = bayes$init    %||% dep("init"),
      L       = bayes$L       %||% dep("L")       %||% 50L,
      epsilon = bayes$epsilon %||% dep("epsilon") %||% 0.002,
      prior   = bayes$prior   %||% dep("prior"),
      seed    = seed, WAIC = WAIC, nsimWAIC = nsimWAIC), extra)))
  }
  return(do.call(.remstimate_mle_cpp, c(list(
    reh = reh, stats = stats, ncores = ncores,
    WAIC = WAIC, nsimWAIC = nsimWAIC), extra)))
}

# ══════════════════════════════════════════════════════════════════════════════
# C++ backend: shared setup
# ══════════════════════════════════════════════════════════════════════════════

.remstimate_cpp_setup <- function(reh, stats, ncores) {
  if (!inherits(reh, "remify")) stop("'reh' must be a 'remify' object.")

  sampled <- inherits(stats, "tomstats_sampled")
  if (!inherits(stats, c("tomstats", "tomstats_sampled", "aomstats")))
    stop("'stats' must be a 'tomstats', 'tomstats_sampled', or 'aomstats' object.")

  # Extract reh fields
  if (!is.null(reh$meta)) {
    model    <- reh$meta$model
    ordinal  <- isTRUE(reh$meta$ordinal)
    riskset  <- reh$meta$riskset_source %||% reh$meta$riskset
    directed <- isTRUE(reh$meta$directed)
    ncores_reh <- reh$meta$ncores %||% 1L
  } else {
    model    <- attr(reh, "model") %||% "tie"
    ordinal  <- isTRUE(attr(reh, "ordinal"))
    riskset  <- attr(reh, "riskset")
    directed <- isTRUE(attr(reh, "directed"))
    ncores_reh <- attr(reh, "ncores") %||% 1L
  }

  if (inherits(stats, "aomstats") && model == "tie")
    stop("'stats' is an 'aomstats' object but 'reh' is tie-oriented")
  if (inherits(stats, c("tomstats", "tomstats_sampled")) && model == "actor")
    stop("'stats' is a 'tomstats' object but 'reh' is actor-oriented.")
  if (!model %in% c("tie", "actor"))
    stop("C++ backends currently support 'tie' and 'actor' models only.")
  if (model == "actor" && !directed)
    stop("actor-oriented modeling can't operate on undirected networks")

  if (is.null(ncores) || ncores < 1L) ncores <- ncores_reh
  max_cores <- parallel::detectCores()
  if (max_cores <= 2L) ncores <- 1L else ncores <- min(ncores, max_cores - 2L)

  # Whether the riskset is a REDUCED one (stats built on the active/manual
  # dyad set). `riskset` may be a *source* label such as "active_saturated"
  # or "active_dynamic", so match any "active*" variant, not just "active".
  reduced_riskset <- grepl("^active", riskset) || identical(riskset, "manual")

  # Dyad IDs and interevent times
  if (!is.null(reh$ids)) {
    if (reduced_riskset)
      dyadID <- as.vector(reh$ids$dyad_active)
    else
      dyadID <- as.vector(reh$ids$dyad)
    intereventTime   <- reh$intereventTime
    simultaneous_idx <- reh$indices_simultaneous_events
    evenly_spaced    <- if (!is.null(reh$simultaneous)) reh$simultaneous$interevent_evenly_spaced else NULL
  } else {
    if (reduced_riskset)
      dyadID <- attr(reh, "dyadIDactive") %||% attr(reh, "dyadID")
    else
      dyadID <- attr(reh, "dyadID")
    intereventTime   <- reh$intereventTime
    simultaneous_idx <- attr(reh, "indices_simultaneous_events")
    evenly_spaced    <- attr(reh, "evenly_spaced_interevent_time")
  }

  D <- if (reduced_riskset) (reh$activeD %||% reh$D) else reh$D
  M <- reh$M

  # Guard: the dyad dimension of a tie-oriented 'stats' array must match the
  # riskset size that 'reh' will index into. When they disagree (e.g. remify
  # and remstats built different risksets), the C++ backend would otherwise
  # fail deep inside Armadillo with a cryptic 'Mat::elem(): index out of
  # bounds'. Fail early with an actionable message instead.
  if (model == "tie" && !sampled && !inherits(stats, "aomstats") && length(dim(stats)) == 3L) {
    stats_D <- dim(stats)[2]
    if (!is.na(stats_D) && stats_D != D) {
      stop(sprintf(
        paste0("'stats' has %d dyad columns but the riskset in 'reh' has %d ",
               "(%s). The remify and remstats objects describe different ",
               "risksets - rebuild both from the same 'reh'."),
        stats_D, D, if (reduced_riskset) "activeD" else "D"),
        call. = FALSE)
    }
  }

  stats_method <- attr(stats, "method")
  if (is.null(stats_method)) stop("'stats' has no 'method' attribute.")

  if (stats_method == "pe") {
    if (!is.null(evenly_spaced)) intereventTime <- as.vector(evenly_spaced)
    if (!is.null(reh$E)) M <- reh$E
  }

  list(
    model = model, ordinal = ordinal, riskset = riskset, directed = directed,
    dyadID = dyadID, intereventTime = intereventTime, D = D, M = M,
    ncores = ncores, sampled = sampled, stats_method = stats_method,
    simultaneous_idx = simultaneous_idx
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# C++ backend: MLE
# ══════════════════════════════════════════════════════════════════════════════

.remstimate_mle_cpp <- function(reh, stats, ncores = 1L,
                                WAIC = FALSE, nsimWAIC = 100L, ...) {
  env <- .remstimate_cpp_setup(reh, stats, ncores)

  if (env$model == "actor") {
    return(.remstimate_mle_cpp_actor(reh, stats, env, WAIC, nsimWAIC))
  } else {
    return(.remstimate_mle_cpp_tie(reh, stats, env, WAIC, nsimWAIC))
  }
}

.remstimate_mle_cpp_tie <- function(reh, stats, env, WAIC, nsimWAIC) {
  subset_attr <- attr(stats, "subset")
  start_stop  <- if (!is.null(subset_attr)) as.integer(unlist(subset_attr)) else c(1L, env$M)

  dyadID         <- env$dyadID
  intereventTime <- env$intereventTime

  if (env$stats_method == "pt") {
    dyadID         <- if (is.list(dyadID)) dyadID[start_stop[1]:start_stop[2]] else dyadID[start_stop[1]:start_stop[2]]
    intereventTime <- intereventTime[start_stop[1]:start_stop[2]]
  } else {
    dyadID         <- if (is.list(dyadID)) unlist(dyadID)[start_stop[1]:start_stop[2]] else dyadID[start_stop[1]:start_stop[2]]
    intereventTime <- intereventTime[start_stop[1]:start_stop[2]]
  }

  M_sub   <- diff(start_stop) + 1L
  ordinal <- env$ordinal
  if (ordinal || is.null(intereventTime)) intereventTime <- rep(1.0, M_sub)

  omit_dyad       <- list()
  variables_names <- dimnames(stats)[[3]]
  if (is.null(variables_names)) stop("'stats' has no dimname labels on third dimension.")
  model_formula   <- attr(stats, "formula") %||%
    stats::as.formula(paste("~", paste(variables_names, collapse = " + ")))
  where_is_baseline <- if (any(tolower(variables_names) == "baseline"))
    which(variables_names == "baseline") else NULL
  P <- length(variables_names)

  dyad_field <- if (is.list(dyadID)) lapply(dyadID, as.integer) else lapply(as.integer(dyadID), identity)

  if (!env$sampled) {
    stats_arr <- aperm(stats, perm = c(2, 3, 1))

    optimum_obj <- trust::trust(
      objfun   = remDerivatives,
      parinit  = rep(0.0, P),
      rinit    = 1, rmax = 100,
      stats    = stats_arr, actor1 = list(), actor2 = list(),
      dyad     = dyad_field, omit_dyad = omit_dyad,
      interevent_time = intereventTime,
      model = "tie", ordinal = ordinal, ncores = env$ncores
    )

    remstimateList <- .build_mle_output(
      optimum_obj, variables_names, P, M_sub, env$D, reh, ordinal,
      stats_arr, dyad_field, omit_dyad, intereventTime, env$ncores,
      WAIC, nsimWAIC, "tie"
    )
  } else {
    stats_arr  <- aperm(stats, perm = c(2, 3, 1))
    sample_map <- attr(stats, "sample_map")
    samp_prob  <- attr(stats, "samp_prob")
    case_pos   <- attr(stats, "case_pos")
    if (is.null(sample_map) || is.null(samp_prob) || is.null(case_pos))
      stop("'stats' is missing required sampling attributes.")
    case_pos_0 <- lapply(case_pos, function(x) as.integer(x) - 1L)

    optimum_obj <- trust::trust(
      objfun    = remDerivativesSampled,
      parinit   = rep(0.0, P),
      rinit     = 1, rmax = 100,
      stats     = stats_arr, case_pos = case_pos_0, samp_prob = samp_prob,
      interevent_time = intereventTime, ordinal = ordinal, ncores = env$ncores
    )

    S <- dim(stats)[2]
    null_obj <- trust::trust(
      objfun    = remDerivativesSampled,
      parinit   = 0.0,
      rinit     = 1, rmax = 100,
      stats     = array(1.0, dim = c(S, 1L, M_sub)),
      case_pos  = case_pos_0, samp_prob = samp_prob,
      interevent_time = intereventTime, ordinal = ordinal, ncores = env$ncores
    )

    remstimateList <- list()
    remstimateList$coefficients     <- as.vector(optimum_obj$argument)
    names(remstimateList$coefficients) <- variables_names
    remstimateList$loglik            <- -optimum_obj$value
    remstimateList$gradient          <- -optimum_obj$gradient
    remstimateList$hessian           <- optimum_obj$hessian
    remstimateList$vcov              <- tryCatch(qr.solve(optimum_obj$hessian),
                                                  error = function(e) matrix(NA, P, P))
    remstimateList$se                <- sqrt(diag(remstimateList$vcov))
    names(remstimateList$se) <- rownames(remstimateList$vcov) <-
      colnames(remstimateList$vcov) <- variables_names
    remstimateList$residual.deviance <- -2 * remstimateList$loglik
    remstimateList$null.deviance     <- 2 * null_obj$value
    remstimateList$model.deviance    <- remstimateList$null.deviance - remstimateList$residual.deviance
    remstimateList$df.null           <- M_sub
    remstimateList$df.model          <- P
    remstimateList$df.residual       <- M_sub - P
    remstimateList$AIC               <- 2 * P - 2 * remstimateList$loglik
    remstimateList$AICC              <- remstimateList$AIC + 2 * P * (P + 1) / max(M_sub - P - 1, 1)
    remstimateList$BIC               <- P * log(M_sub) - 2 * remstimateList$loglik
    remstimateList$converged         <- optimum_obj$converged
    remstimateList$iterations        <- optimum_obj$iterations
    remstimateList$sampled           <- TRUE
    remstimateList$samp_num          <- attr(stats, "samp_num")
    remstimateList$sampling_scheme   <- attr(stats, "sampling_scheme")
  }

  str_out <- structure(remstimateList, class = "remstimate")
  attr(str_out, "formula")           <- model_formula
  attr(str_out, "model")             <- "tie"
  attr(str_out, "ordinal")           <- ordinal
  attr(str_out, "method")            <- "MLE"
  attr(str_out, "approach")          <- "Frequentist"
  attr(str_out, "statistics")        <- variables_names
  attr(str_out, "where_is_baseline") <- where_is_baseline
  attr(str_out, "ncores")            <- env$ncores
  attr(str_out, "sampled")           <- env$sampled
  str_out
}

.remstimate_mle_cpp_actor <- function(reh, stats, env, WAIC, nsimWAIC) {
  start_stop <- as.numeric(unlist(attr(stats, "subset")))

  actor1_ids <- reh$ids$actor1[start_stop[1]:start_stop[2]]
  actor2_ids <- reh$ids$actor2[start_stop[1]:start_stop[2]]
  N <- reh$N

  M_sub_sender   <- length(actor1_ids)
  M_sub_receiver <- length(actor2_ids)

  intereventTime_sender   <- if (env$ordinal) rep(1.0, M_sub_sender) else reh$intereventTime[start_stop[1]:start_stop[2]]
  intereventTime_receiver <- rep(1.0, M_sub_receiver)

  sender_rs   <- if (!is.null(reh$sender_riskset))   as.integer(reh$sender_riskset) - 1L   else integer(0)
  receiver_rs <- if (!is.null(reh$receiver_riskset)) lapply(reh$receiver_riskset, function(x) as.integer(x) - 1L) else list()

  omit_dyad_sender   <- reh$omit_dyad %||% list()
  omit_dyad_receiver <- list()

  which_stats  <- c("sender_stats",   "receiver_stats")
  which_model  <- c("sender_model",   "receiver_model")
  sender_rate  <- c(TRUE, FALSE)

  remstimateList <- list()

  for (i in 1:2) {
    s <- stats[[ which_stats[i] ]]
    if (is.null(s)) next

    s_arr       <- aperm(s, perm = c(2, 3, 1))
    P_i         <- dim(s_arr)[2]
    var_names_i <- dimnames(s)[[3]]

    if (sender_rate[i]) {
      intereventTime <- intereventTime_sender
      M_sub          <- M_sub_sender
      actor1_field   <- lapply(actor1_ids, as.integer)
      actor2_field   <- lapply(actor2_ids, as.integer)
    } else {
      intereventTime <- intereventTime_receiver
      M_sub          <- M_sub_receiver
      actor1_field   <- lapply(actor1_ids, function(x) as.integer(unlist(x)))
      actor2_field   <- lapply(actor2_ids, function(x) as.integer(unlist(x)))
    }

    omit_dyad_i <- if (sender_rate[i]) omit_dyad_sender else omit_dyad_receiver

    opt <- trust::trust(
      objfun          = remDerivatives,
      parinit         = rep(0.0, P_i),
      rinit           = 1, rmax = 100,
      stats           = s_arr,
      actor1          = actor1_field,
      actor2          = actor2_field,
      dyad            = list(),
      omit_dyad       = omit_dyad_i,
      interevent_time = intereventTime,
      model           = "actor",
      ordinal         = env$ordinal,
      ncores          = env$ncores,
      senderRate      = sender_rate[i],
      N               = N,
      sender_riskset   = if (i == 1) sender_rs else integer(0),
      receiver_riskset = if (i == 2) receiver_rs else list()
    )

    res <- list()
    res$coefficients      <- as.vector(opt$argument)
    names(res$coefficients) <- var_names_i
    res$loglik            <- -opt$value
    res$gradient          <- -opt$gradient
    res$hessian           <- opt$hessian
    res$vcov              <- tryCatch(qr.solve(opt$hessian), error = function(e) matrix(NA, P_i, P_i))
    res$se                <- sqrt(diag(res$vcov))
    names(res$se) <- rownames(res$vcov) <- colnames(res$vcov) <- var_names_i
    res$residual.deviance <- -2 * res$loglik
    res$AIC               <- 2 * P_i - 2 * res$loglik
    res$AICC              <- res$AIC + 2 * P_i * (P_i + 1) / max(M_sub - P_i - 1, 1)
    res$BIC               <- P_i * log(M_sub) - 2 * res$loglik
    res$converged         <- opt$converged
    res$iterations        <- opt$iterations
    res$df.null           <- M_sub
    res$df.model          <- P_i
    res$df.residual       <- M_sub - P_i

    null_dim <- dim(s_arr); null_dim[2] <- 1L
    res$null.deviance <- if (env$ordinal) {
      risk_set_size <- if (sender_rate[i]) N else N - 1L
      2 * M_sub * log(risk_set_size)
    } else {
      2 * trust::trust(
        objfun = remDerivatives, parinit = 0.0, rinit = 1, rmax = 100,
        stats = array(1.0, dim = null_dim),
        actor1 = actor1_field, actor2 = actor2_field,
        dyad = list(), omit_dyad = omit_dyad_i,
        interevent_time = intereventTime,
        model = "actor", ordinal = env$ordinal, ncores = env$ncores,
        senderRate = sender_rate[i], N = N
      )$value
    }
    res$model.deviance <- res$null.deviance - res$residual.deviance

    if (WAIC) {
      res$WAIC <- getWAIC(
        mu = res$coefficients, vcov = res$vcov,
        pars = matrix(0, 2, 2), stats = s_arr,
        actor1 = actor1_field, actor2 = actor2_field,
        dyad = list(), omit_dyad = omit_dyad_i,
        interevent_time = intereventTime,
        model = "actor", approach = "Frequentist",
        nsim = nsimWAIC, ordinal = env$ordinal, ncores = env$ncores,
        senderRate = sender_rate[i],
        sender_riskset = if (i == 1) sender_rs else integer(0),
        receiver_riskset = if (i == 2) receiver_rs else list()
      )
    }

    remstimateList[[ which_model[i] ]] <- res
  }

  str_out <- structure(remstimateList, class = "remstimate")
  attr(str_out, "model")      <- "actor"
  attr(str_out, "ordinal")    <- env$ordinal
  attr(str_out, "method")     <- "MLE"
  attr(str_out, "approach")   <- "Frequentist"
  attr(str_out, "formula")    <- list(
    rate_model_formula   = attr(stats, "formula")$rate,
    choice_model_formula = attr(stats, "formula")$choice
  )
  attr(str_out, "statistics") <- list(
    sender_model   = dimnames(stats$sender_stats)[[3]],
    receiver_model = dimnames(stats$receiver_stats)[[3]]
  )
  attr(str_out, "ncores") <- env$ncores
  str_out
}

# ══════════════════════════════════════════════════════════════════════════════
# C++ backend: HMC
# ══════════════════════════════════════════════════════════════════════════════

.remstimate_hmc_cpp <- function(reh, stats, ncores = 1L,
                                nsim = 500L, nchains = 1L, burnin = 500L,
                                thin = 10L, init = NULL, L = 50L,
                                epsilon = 0.002, prior = NULL, seed = NULL,
                                WAIC = FALSE, nsimWAIC = 100L, ...) {
  env <- .remstimate_cpp_setup(reh, stats, ncores)

  # ── HMC parameter validation ─────────────────────────────────────────────
  if (is.null(seed)) seed <- .Random.seed
  else if (length(seed) > 1) { seed <- seed[1]; warning("`seed` length > 1; using first element.") }
  if (is.null(nchains)) nchains <- 1L
  if (is.null(nsim)) nsim <- 1000L
  if (is.null(burnin)) burnin <- 500L
  if (is.null(thin)) {
    thin <- if (nsim >= 100) 10L else 1L
  } else if (thin <= 0) {
    stop("`thin` must be positive. Set thin = 1 for no thinning.")
  }
  if (is.null(L) || !is.numeric(L) || L <= 1) L <- 50L
  if (is.null(epsilon) || !is.numeric(epsilon) || epsilon <= 0) epsilon <- 0.1 / L

  if (env$model == "actor") {
    return(.remstimate_hmc_cpp_actor(reh, stats, env,
                                     nsim, nchains, burnin, thin,
                                     init, L, epsilon, prior, seed))
  } else {
    return(.remstimate_hmc_cpp_tie(reh, stats, env,
                                    nsim, nchains, burnin, thin,
                                    init, L, epsilon, prior, seed))
  }
}

.remstimate_hmc_cpp_tie <- function(reh, stats, env,
                                     nsim, nchains, burnin, thin,
                                     init, L, epsilon, prior, seed) {
  subset_attr <- attr(stats, "subset")
  start_stop  <- if (!is.null(subset_attr)) as.integer(unlist(subset_attr)) else c(1L, env$M)

  dyadID         <- env$dyadID
  intereventTime <- env$intereventTime

  if (env$stats_method == "pt") {
    dyadID         <- if (is.list(dyadID)) dyadID[start_stop[1]:start_stop[2]] else dyadID[start_stop[1]:start_stop[2]]
    intereventTime <- intereventTime[start_stop[1]:start_stop[2]]
  } else {
    dyadID         <- if (is.list(dyadID)) unlist(dyadID)[start_stop[1]:start_stop[2]] else dyadID[start_stop[1]:start_stop[2]]
    intereventTime <- intereventTime[start_stop[1]:start_stop[2]]
  }

  M_sub   <- diff(start_stop) + 1L
  ordinal <- env$ordinal
  if (ordinal || is.null(intereventTime)) intereventTime <- rep(1.0, M_sub)

  omit_dyad       <- list()
  variables_names <- dimnames(stats)[[3]]
  where_is_baseline <- if (any(tolower(variables_names) == "baseline"))
    which(variables_names == "baseline") else NULL
  P <- length(variables_names)

  dyad_field <- if (is.list(dyadID)) lapply(dyadID, as.integer) else lapply(as.integer(dyadID), identity)

  if (!env$sampled) {
    stats_arr <- aperm(stats, perm = c(2, 3, 1))

    if (is.null(init)) {
      mle_obj <- trust::trust(
        objfun = remDerivatives, parinit = rep(0.0, P), rinit = 1, rmax = 100,
        stats = stats_arr, actor1 = list(), actor2 = list(),
        dyad = dyad_field, omit_dyad = omit_dyad,
        interevent_time = intereventTime,
        model = "tie", ordinal = ordinal, ncores = env$ncores
      )
      init <- mle_obj$argument
    }

    if (!is.null(seed)) set.seed(seed)
    prior_mean <- if (!is.null(prior$mean)) prior$mean else rep(0.0, P)
    prior_vcov <- if (!is.null(prior$vcov)) prior$vcov else diag(100, P)

    hmc_out <- HMC(
      pars_init       = matrix(init, nrow = P, ncol = nchains),
      nsim            = as.integer(nsim + burnin),
      nchains         = as.integer(nchains),
      burnin          = as.integer(burnin),
      thin            = as.integer(thin),
      L               = as.integer(L),
      epsilon         = epsilon,
      meanPrior       = prior_mean,
      sigmaPrior      = prior_vcov,
      stats           = stats_arr,
      actor1          = list(),
      actor2          = list(),
      dyad            = dyad_field,
      omit_dyad       = omit_dyad,
      interevent_time = intereventTime,
      model           = "tie",
      ordinal         = ordinal,
      ncores          = env$ncores
    )

    remstimateList <- .build_hmc_output(hmc_out, variables_names, P)
    remstimateList$df.null <- M_sub

  } else {
    # Sampled HMC
    stats_arr  <- aperm(stats, perm = c(2, 3, 1))
    samp_prob  <- attr(stats, "samp_prob")
    case_pos   <- attr(stats, "case_pos")
    case_pos_0 <- lapply(case_pos, function(x) as.integer(x) - 1L)

    if (is.null(init)) {
      mle_obj <- trust::trust(
        objfun = remDerivativesSampled, parinit = rep(0.0, P),
        rinit = 1, rmax = 100, stats = stats_arr,
        case_pos = case_pos_0, samp_prob = samp_prob,
        interevent_time = intereventTime, ordinal = ordinal, ncores = env$ncores
      )
      init <- mle_obj$argument
    }

    if (!is.null(seed)) set.seed(seed)
    prior_mean <- if (!is.null(prior$mean)) prior$mean else rep(0.0, P)
    prior_vcov <- if (!is.null(prior$vcov)) prior$vcov else diag(100, P)

    chain_draws    <- matrix(NA_real_, nrow = nsim, ncol = P)
    chain_logpost  <- numeric(nsim)
    beta_cur       <- init

    logpost_fn <- function(b) {
      ll <- -remDerivativesSampled(b, stats_arr, case_pos_0, samp_prob,
                                    intereventTime, ordinal,
                                    gradient = FALSE, hessian = FALSE,
                                    ncores = env$ncores)$value
      lp <- sum(mvnfast::dmvn(matrix(b, nrow = 1), prior_mean, prior_vcov, log = TRUE))
      ll + lp
    }
    grad_fn <- function(b) {
      g <- -remDerivativesSampled(b, stats_arr, case_pos_0, samp_prob,
                                   intereventTime, ordinal,
                                   gradient = TRUE, hessian = FALSE,
                                   ncores = env$ncores)$gradient
      gp <- solve(prior_vcov, b - prior_mean)
      g - gp
    }

    logpost_cur <- logpost_fn(beta_cur)
    for (iter in seq_len(nsim + burnin)) {
      r <- stats::rnorm(P)
      beta_prop <- beta_cur
      r_prop <- r + 0.5 * epsilon * grad_fn(beta_prop)
      for (step in seq_len(L)) {
        beta_prop <- beta_prop + epsilon * r_prop
        if (step < L) r_prop <- r_prop + epsilon * grad_fn(beta_prop)
      }
      r_prop <- r_prop + 0.5 * epsilon * grad_fn(beta_prop)
      r_prop <- -r_prop

      logpost_prop <- logpost_fn(beta_prop)
      log_accept   <- logpost_prop - logpost_cur + 0.5 * sum(r^2) - 0.5 * sum(r_prop^2)
      if (log(stats::runif(1)) < log_accept) {
        beta_cur    <- beta_prop
        logpost_cur <- logpost_prop
      }
      if (iter > burnin) {
        chain_draws[iter - burnin, ]  <- beta_cur
        chain_logpost[iter - burnin]  <- logpost_cur
      }
    }

    keep <- seq(thin, nsim, by = thin)
    draws <- chain_draws[keep, , drop = FALSE]
    colnames(draws) <- variables_names

    remstimateList <- list()
    remstimateList$draws         <- draws
    remstimateList$log_posterior <- -chain_logpost[keep]
    remstimateList$coefficients  <- draws[which.min(remstimateList$log_posterior), ]
    remstimateList$post.mean     <- colMeans(draws)
    remstimateList$vcov          <- stats::cov(draws)
    remstimateList$sd            <- sqrt(diag(remstimateList$vcov))
    remstimateList$loglik        <- max(-remstimateList$log_posterior)
    remstimateList$df.null       <- M_sub
    names(remstimateList$coefficients) <- names(remstimateList$post.mean) <-
      names(remstimateList$sd) <- rownames(remstimateList$vcov) <-
      colnames(remstimateList$vcov) <- variables_names
    remstimateList$sampled <- TRUE
  }

  model_formula <- attr(stats, "formula") %||%
    stats::as.formula(paste("~", paste(variables_names, collapse = " + ")))

  str_out <- structure(remstimateList, class = "remstimate")
  attr(str_out, "formula")           <- model_formula
  attr(str_out, "model")             <- "tie"
  attr(str_out, "ordinal")           <- ordinal
  attr(str_out, "method")            <- "HMC"
  attr(str_out, "approach")          <- "Bayesian"
  attr(str_out, "statistics")        <- variables_names
  attr(str_out, "where_is_baseline") <- where_is_baseline
  attr(str_out, "ncores")            <- env$ncores
  attr(str_out, "sampled")           <- env$sampled
  attr(str_out, "nsim")    <- nsim
  attr(str_out, "nchains") <- nchains
  attr(str_out, "burnin")  <- burnin
  attr(str_out, "thin")    <- thin
  attr(str_out, "L")       <- L
  attr(str_out, "epsilon") <- epsilon
  attr(str_out, "seed")    <- seed
  attr(str_out, "prior")   <- prior
  str_out
}

.remstimate_hmc_cpp_actor <- function(reh, stats, env,
                                       nsim, nchains, burnin, thin,
                                       init, L, epsilon, prior, seed) {
  start_stop <- as.numeric(unlist(attr(stats, "subset")))

  actor1_ids <- reh$ids$actor1[start_stop[1]:start_stop[2]]
  actor2_ids <- reh$ids$actor2[start_stop[1]:start_stop[2]]
  N <- reh$N

  M_sub_sender   <- length(actor1_ids)
  M_sub_receiver <- length(actor2_ids)

  intereventTime_sender   <- if (env$ordinal) rep(1.0, M_sub_sender) else reh$intereventTime[start_stop[1]:start_stop[2]]
  intereventTime_receiver <- rep(1.0, M_sub_receiver)

  sender_rs   <- if (!is.null(reh$sender_riskset))   as.integer(reh$sender_riskset) - 1L   else integer(0)
  receiver_rs <- if (!is.null(reh$receiver_riskset)) lapply(reh$receiver_riskset, function(x) as.integer(x) - 1L) else list()

  omit_dyad_sender   <- reh$omit_dyad %||% list()
  omit_dyad_receiver <- list()

  which_stats  <- c("sender_stats",   "receiver_stats")
  which_model  <- c("sender_model",   "receiver_model")
  sender_rate  <- c(TRUE, FALSE)

  remstimateList <- list()

  for (i in 1:2) {
    s <- stats[[ which_stats[i] ]]
    if (is.null(s)) next

    s_arr       <- aperm(s, perm = c(2, 3, 1))
    P_i         <- dim(s_arr)[2]
    var_names_i <- dimnames(s)[[3]]

    if (sender_rate[i]) {
      intereventTime <- intereventTime_sender
      M_sub          <- M_sub_sender
      actor1_field   <- lapply(actor1_ids, as.integer)
      actor2_field   <- lapply(actor2_ids, as.integer)
    } else {
      intereventTime <- intereventTime_receiver
      M_sub          <- M_sub_receiver
      actor1_field   <- lapply(actor1_ids, function(x) as.integer(unlist(x)))
      actor2_field   <- lapply(actor2_ids, function(x) as.integer(unlist(x)))
    }

    omit_dyad_i <- if (sender_rate[i]) omit_dyad_sender else omit_dyad_receiver

    if (is.null(init)) init <- list()
    if (is.null(init[[ which_model[i] ]])) {
      init[[ which_model[i] ]] <- matrix(
        stats::runif(P_i * nchains, -0.1, 0.1),
        nrow = P_i, ncol = nchains
      )
    }
    if (!is.null(seed)) set.seed(seed)

    hmc_out <- HMC(
      pars_init       = init[[ which_model[i] ]],
      nsim            = as.integer(nsim + burnin),
      nchains         = as.integer(nchains),
      burnin          = as.integer(burnin),
      thin            = as.integer(thin),
      L               = as.integer(L),
      epsilon         = epsilon,
      meanPrior       = rep(0.0, P_i),
      sigmaPrior      = diag(100, P_i),
      stats           = s_arr,
      actor1          = actor1_field,
      actor2          = actor2_field,
      dyad            = list(),
      omit_dyad       = omit_dyad_i,
      interevent_time = intereventTime,
      model           = "actor",
      ordinal         = env$ordinal,
      ncores          = env$ncores,
      senderRate      = sender_rate[i],
      N               = N,
      sender_riskset  = if (i == 1) sender_rs else integer(0),
      receiver_riskset = if (i == 2) receiver_rs else list()
    )

    res <- list()
    res$coefficients <- hmc_out$draws[which.min(hmc_out$log_posterior), ]
    res$post.mean    <- colMeans(hmc_out$draws)
    res$vcov         <- stats::cov(hmc_out$draws)
    res$sd           <- sqrt(diag(res$vcov))
    res$loglik       <- -min(hmc_out$log_posterior)
    res$draws        <- hmc_out$draws
    res$df.null      <- M_sub
    res$df.model     <- P_i
    res$df.residual  <- M_sub - P_i
    names(res$coefficients) <- names(res$post.mean) <-
      names(res$sd) <- rownames(res$vcov) <- colnames(res$vcov) <- var_names_i

    remstimateList[[ which_model[i] ]] <- res
  }

  str_out <- structure(remstimateList, class = "remstimate")
  attr(str_out, "model")      <- "actor"
  attr(str_out, "ordinal")    <- env$ordinal
  attr(str_out, "method")     <- "HMC"
  attr(str_out, "approach")   <- "Bayesian"
  attr(str_out, "formula")    <- list(
    rate_model_formula   = attr(stats, "formula")$rate,
    choice_model_formula = attr(stats, "formula")$choice
  )
  attr(str_out, "statistics") <- list(
    sender_model   = dimnames(stats$sender_stats)[[3]],
    receiver_model = dimnames(stats$receiver_stats)[[3]]
  )
  attr(str_out, "ncores")  <- env$ncores
  attr(str_out, "nsim")    <- nsim
  attr(str_out, "nchains") <- nchains
  attr(str_out, "burnin")  <- burnin
  attr(str_out, "thin")    <- thin
  attr(str_out, "L")       <- L
  attr(str_out, "epsilon") <- epsilon
  attr(str_out, "seed")    <- seed
  str_out
}

# ══════════════════════════════════════════════════════════════════════════════
# Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

# NULL-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Build MLE output list for full stats path
.build_mle_output <- function(optimum_obj, variables_names, P, M, D,
                                reh, ordinal, stats_arr, dyad_field,
                                omit_dyad, intereventTime, ncores, WAIC,
                                nsimWAIC, model) {
  out <- list()
  out$coefficients     <- as.vector(optimum_obj$argument)
  names(out$coefficients) <- variables_names
  out$loglik           <- -optimum_obj$value
  out$gradient         <- -optimum_obj$gradient
  out$hessian          <- optimum_obj$hessian
  out$vcov             <- tryCatch(qr.solve(optimum_obj$hessian),
                                    error = function(e) matrix(NA, P, P))
  out$se               <- sqrt(diag(out$vcov))
  names(out$se) <- rownames(out$vcov) <- colnames(out$vcov) <- variables_names
  out$residual.deviance <- -2 * out$loglik

  null_stats <- array(1.0, dim = c(D, 1L, M))
  null_obj <- trust::trust(
    objfun   = remDerivatives,
    parinit  = 0.0,
    rinit    = 1, rmax = 100,
    stats    = null_stats, actor1 = list(), actor2 = list(),
    dyad     = dyad_field, omit_dyad = omit_dyad,
    interevent_time = intereventTime,
    model = "tie", ordinal = ordinal, ncores = ncores
  )
  out$null.deviance    <- 2 * null_obj$value
  out$model.deviance   <- out$null.deviance - out$residual.deviance
  out$df.null          <- M
  out$df.model         <- P
  out$df.residual      <- M - P
  out$AIC              <- 2 * P - 2 * out$loglik
  out$AICC             <- out$AIC + 2 * P * (P + 1) / max(M - P - 1, 1)
  out$BIC              <- P * log(M) - 2 * out$loglik

  if (WAIC) {
    out$WAIC <- getWAIC(
      mu = out$coefficients, vcov = out$vcov,
      pars = matrix(0, 2, 2), stats = stats_arr,
      actor1 = list(), actor2 = list(),
      dyad = dyad_field, omit_dyad = omit_dyad,
      interevent_time = intereventTime,
      model = model, approach = "Frequentist",
      nsim = nsimWAIC, ordinal = ordinal, ncores = ncores
    )
  }

  out$converged  <- optimum_obj$converged
  out$iterations <- optimum_obj$iterations
  out$sampled    <- FALSE
  out
}

# Build HMC output list
.build_hmc_output <- function(hmc_out, variables_names, P) {
  out <- list()
  out$draws        <- hmc_out$draws
  out$log_posterior <- hmc_out$log_posterior
  colnames(out$draws) <- variables_names
  out$coefficients <- out$draws[which.min(out$log_posterior), ]
  out$post.mean    <- colMeans(out$draws)
  out$vcov         <- stats::cov(out$draws)
  out$sd           <- sqrt(diag(out$vcov))
  out$loglik       <- -min(out$log_posterior)
  names(out$coefficients) <- names(out$post.mean) <- names(out$sd) <-
    rownames(out$vcov) <- colnames(out$vcov) <- variables_names
  out$sampled <- FALSE
  out
}


#' @export
#' @method logLik remstimate
logLik.remstimate <- function(object, ...) {
  structure(object$loglik, class = "logLik",
            df = length(object$coefficients),
            nobs = object$df.null)
}

