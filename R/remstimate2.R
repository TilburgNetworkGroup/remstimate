#' @title remstimate2 - MLE and HMC estimation for tie-oriented REM
#'
#' @description Estimates tie-oriented relational event model parameters using
#' Maximum Likelihood Estimation (\code{"MLE"}) or Hamiltonian Monte Carlo
#' (\code{"HMC"}). Supports both full statistics (\code{tomstats}) and
#' case-control sampled statistics (\code{tomstats_sampled}).
#'
#' @param reh a \code{remify} object (output of \code{remify::remify2()}).
#' @param stats a \code{tomstats} or \code{tomstats_sampled} object (output of
#'   \code{remstats::tomstats()} or \code{remstats::tomstats2(sampling=TRUE)}).
#' @param method optimization method: \code{"MLE"} (default) or \code{"HMC"}.
#' @param ncores number of threads for parallelization (default 1).
#' @param nsim [\emph{HMC only}] number of MCMC iterations per chain. Default 1000.
#' @param nchains [\emph{HMC only}] number of chains. Default 1.
#' @param burnin [\emph{HMC only}] burn-in iterations. Default 500.
#' @param thin [\emph{HMC only}] thinning interval. Default 10.
#' @param init [\emph{HMC only}] vector of initial parameter values. If NULL,
#'   MLE estimates are used as starting values.
#' @param L [\emph{HMC only}] number of leapfrog steps. Default 50.
#' @param epsilon [\emph{HMC only}] leapfrog step size. Default 0.002.
#' @param prior [\emph{HMC only}] list with elements \code{mean} (prior mean
#'   vector) and \code{vcov} (prior covariance matrix). If NULL, a standard
#'   normal prior is used.
#' @param seed [\emph{HMC only}] random seed for reproducibility.
#' @param WAIC logical. Compute WAIC? Default FALSE. (MLE only for now.)
#'
#' @return a \code{remstimate} S3 object.
#'
#' @export
remstimate2 <- function(reh,
                        stats,
                        method  = c("MLE", "HMC"),
                        ncores  = 1L,
                        nsim    = 1000L,
                        nchains = 1L,
                        burnin  = 500L,
                        thin    = 10L,
                        init    = NULL,
                        L       = 50L,
                        epsilon = 0.002,
                        prior   = NULL,
                        seed    = NULL,
                        WAIC    = FALSE) {

  # ── Input checks ────────────────────────────────────────────────────────────
  if (!inherits(reh, "remify")) stop("'reh' must be a 'remify' object.")

  sampled <- inherits(stats, "tomstats_sampled")
  if (!inherits(stats, "remstats")) {
    stop("'stats' must be a 'tomstats' or 'tomstats_sampled' object.")
  }
  if (!inherits(stats, c("tomstats", "tomstats_sampled"))) {
    stop("'stats' must be a tie-oriented remstats object (tomstats or tomstats_sampled).")
  }

  method <- match.arg(method)

  # ── Extract reh fields (new remify2 structure) ───────────────────────────────
  # Support both old (attr-based) and new ($meta/$ids) remify structure
  if (!is.null(reh$meta)) {
    model   <- reh$meta$model
    ordinal <- isTRUE(reh$meta$ordinal)
    riskset <- reh$meta$riskset_source %||% reh$meta$riskset
    directed <- isTRUE(reh$meta$directed)
    ncores_reh <- reh$meta$ncores %||% 1L
  } else {
    model   <- attr(reh, "model")  %||% "tie"
    ordinal <- isTRUE(attr(reh, "ordinal"))
    riskset <- attr(reh, "riskset")
    directed <- isTRUE(attr(reh, "directed"))
    ncores_reh <- attr(reh, "ncores") %||% 1L
  }

  if (model != "tie") stop("remstimate2 currently supports tie-oriented models only.")
  if (is.null(ncores) || ncores < 1L) ncores <- ncores_reh

  # ── Clamp ncores ────────────────────────────────────────────────────────────
  max_cores <- parallel::detectCores()
  if (max_cores <= 2L) ncores <- 1L else ncores <- min(ncores, max_cores - 2L)

  # ── Extract dyad IDs and interevent times ────────────────────────────────────
  # New remify2 format: ids$dyad_active (active riskset) or ids$dyad
  # Old remify format:  attr(reh, "dyadID") / attr(reh, "dyadIDactive")
  if (!is.null(reh$ids)) {
    # New format — dyad_active is a [N x 1] matrix, flatten to vector
    if (riskset == "active") {
      dyadID <- as.vector(reh$ids$dyad_active)
    } else {
      dyadID <- as.vector(reh$ids$dyad)
    }
    intereventTime <- reh$intereventTime
    simultaneous_idx <- reh$indices_simultaneous_events
    evenly_spaced <- if (!is.null(reh$simultaneous)) reh$simultaneous$interevent_evenly_spaced else NULL
  } else {
    # Old format
    if (riskset == "active") {
      dyadID <- attr(reh, "dyadIDactive") %||% attr(reh, "dyadID")
    } else {
      dyadID <- attr(reh, "dyadID")
    }
    intereventTime <- reh$intereventTime
    simultaneous_idx <- attr(reh, "indices_simultaneous_events")
    evenly_spaced <- attr(reh, "evenly_spaced_interevent_time")
  }

  # ── D: riskset size ──────────────────────────────────────────────────────────
  D <- if (riskset == "active") (reh$activeD %||% reh$D) else reh$D
  M <- reh$M

  # ── stats method (pt / pe) ───────────────────────────────────────────────────
  stats_method <- attr(stats, "method")
  if (is.null(stats_method)) stop("'stats' has no 'method' attribute.")

  if (stats_method == "pe") {
    if (!is.null(evenly_spaced)) intereventTime <- as.vector(evenly_spaced)
    if (!is.null(reh$E)) M <- reh$E
  }

  # ── start / stop subsetting ──────────────────────────────────────────────────
  subset_attr <- attr(stats, "subset")
  start_stop <- if (!is.null(subset_attr)) as.integer(unlist(subset_attr)) else c(1L, M)

  # Subset dyadID and intereventTime
  if (stats_method == "pt") {
    if (is.list(dyadID)) {
      dyadID <- dyadID[start_stop[1]:start_stop[2]]
    } else {
      dyadID <- dyadID[start_stop[1]:start_stop[2]]
    }
    if (!is.null(simultaneous_idx)) {
      # For pt: remove simultaneous event duplicate rows from intereventTime
      intereventTime <- intereventTime[start_stop[1]:start_stop[2]]
    } else {
      intereventTime <- intereventTime[start_stop[1]:start_stop[2]]
    }
  } else { # pe
    if (is.list(dyadID)) dyadID <- unlist(dyadID)[start_stop[1]:start_stop[2]]
    else dyadID <- dyadID[start_stop[1]:start_stop[2]]
    intereventTime <- intereventTime[start_stop[1]:start_stop[2]]
  }
  M_sub <- diff(start_stop) + 1L

  if (ordinal) intereventTime <- rep(1.0, M_sub)

  # ── omit_dyad ────────────────────────────────────────────────────────────────
  omit_dyad <- reh$omit_dyad %||% list()
  if (riskset == "active") omit_dyad <- list() # active riskset: no omit_dyad
  if (length(omit_dyad) > 0 && !is.null(omit_dyad$time)) {
    if (stats_method == "pt" && !is.null(simultaneous_idx)) {
      omit_dyad$time <- omit_dyad$time[-simultaneous_idx]
    }
    omit_dyad$time <- omit_dyad$time[start_stop[1]:start_stop[2]]
  }

  # ── Variable names and formula ───────────────────────────────────────────────
  variables_names <- dimnames(stats)[[3]]
  if (is.null(variables_names)) stop("'stats' has no dimname labels on third dimension.")

  model_formula <- attr(stats, "formula") %||%
    stats::as.formula(paste("~", paste(variables_names, collapse = " + ")))

  where_is_baseline <- if (any(tolower(variables_names) == "baseline"))
    which(variables_names == "baseline") else NULL

  P <- length(variables_names)

  # ── Reshape stats and build inputs ──────────────────────────────────────────
  remstimateList <- list()

  if (!sampled) {
    # ── Full stats path ─────────────────────────────────────────────────────
    # Reshape [M x D x P] -> [D x P x M] for remDerivatives
    stats_arr <- aperm(stats, perm = c(2, 3, 1))

    # Convert dyadID to field<uvec> for C++ (1-based — C++ subtracts 1 internally)
    if (is.list(dyadID)) {
      dyad_field <- lapply(dyadID, function(x) as.integer(x))
    } else {
      dyad_field <- lapply(as.integer(dyadID), function(x) x)
    }

    if (method == "MLE") {
      optimum_obj <- trust::trust(
        objfun   = remDerivatives,
        parinit  = rep(0.0, P),
        rinit    = 1,
        rmax     = 100,
        stats    = stats_arr,
        actor1   = list(),
        actor2   = list(),
        dyad     = dyad_field,
        omit_dyad = omit_dyad,
        interevent_time = intereventTime,
        model    = "tie",
        ordinal  = ordinal,
        ncores   = ncores
      )

      remstimateList <- .build_mle_output(
        optimum_obj, variables_names, P, M_sub, D, reh, ordinal,
        stats_arr, dyad_field, omit_dyad, intereventTime, ncores, WAIC
      )

    } else { # HMC
      if (is.null(init)) {
        # Use MLE as starting point
        mle_obj <- trust::trust(
          objfun   = remDerivatives,
          parinit  = rep(0.0, P),
          rinit    = 1, rmax = 100,
          stats    = stats_arr, actor1 = list(), actor2 = list(),
          dyad     = dyad_field, omit_dyad = omit_dyad,
          interevent_time = intereventTime,
          model = "tie", ordinal = ordinal, ncores = ncores
        )
        init <- mle_obj$argument
      }

      if (!is.null(seed)) set.seed(seed)

      prior_mean <- if (!is.null(prior$mean)) prior$mean else rep(0.0, P)
      prior_vcov <- if (!is.null(prior$vcov)) prior$vcov else diag(100, P)

      hmc_out <- HMC(
        pars_init       = matrix(init, nrow = P, ncol = nchains),
        nsim            = as.integer(nsim),
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
        ncores          = ncores
      )

      remstimateList <- .build_hmc_output(hmc_out, variables_names, P)
    }

  } else {
    # ── Sampled stats path ──────────────────────────────────────────────────
    # stats is [M x S x P], reshape to [S x P x M] for remDerivativesSampled
    stats_arr <- aperm(stats, perm = c(2, 3, 1))

    sample_map  <- attr(stats, "sample_map")   # [M x S], 1-based
    samp_prob   <- attr(stats, "samp_prob")    # [M x S]
    case_pos    <- attr(stats, "case_pos")     # list of length M, 1-based

    if (is.null(sample_map) || is.null(samp_prob) || is.null(case_pos)) {
      stop("'stats' is missing required sampling attributes (sample_map, samp_prob, case_pos).")
    }

    # Convert case_pos to 0-based
    case_pos_0 <- lapply(case_pos, function(x) as.integer(x) - 1L)

    if (method == "MLE") {
      optimum_obj <- trust::trust(
        objfun    = remDerivativesSampled,
        parinit   = rep(0.0, P),
        rinit     = 1,
        rmax      = 100,
        stats     = stats_arr,
        case_pos  = case_pos_0,
        samp_prob = samp_prob,
        interevent_time = intereventTime,
        ordinal   = ordinal,
        ncores    = ncores
      )

      # Null deviance: intercept-only sampled model
      S <- dim(stats)[2]
      null_stats <- array(1.0, dim = c(S, 1L, M_sub))
      null_samp_prob <- samp_prob
      null_obj <- trust::trust(
        objfun    = remDerivativesSampled,
        parinit   = 0.0,
        rinit     = 1, rmax = 100,
        stats     = null_stats,
        case_pos  = case_pos_0,
        samp_prob = null_samp_prob,
        interevent_time = intereventTime,
        ordinal   = ordinal,
        ncores    = ncores
      )

      remstimateList$coefficients    <- as.vector(optimum_obj$argument)
      names(remstimateList$coefficients) <- variables_names
      remstimateList$loglik          <- -optimum_obj$value
      remstimateList$gradient        <- -optimum_obj$gradient
      remstimateList$hessian         <- optimum_obj$hessian
      remstimateList$vcov            <- tryCatch(qr.solve(optimum_obj$hessian),
                                                  error = function(e) matrix(NA, P, P))
      remstimateList$se              <- sqrt(diag(remstimateList$vcov))
      names(remstimateList$se) <- rownames(remstimateList$vcov) <-
        colnames(remstimateList$vcov) <- variables_names
      remstimateList$residual.deviance <- -2 * remstimateList$loglik
      remstimateList$null.deviance    <- 2 * null_obj$value
      remstimateList$model.deviance   <- remstimateList$null.deviance - remstimateList$residual.deviance
      remstimateList$df.null          <- M_sub
      remstimateList$df.model         <- P
      remstimateList$df.residual      <- M_sub - P
      remstimateList$AIC              <- 2 * P - 2 * remstimateList$loglik
      remstimateList$AICC             <- remstimateList$AIC +
        2 * P * (P + 1) / max(M_sub - P - 1, 1)
      remstimateList$BIC              <- P * log(M_sub) - 2 * remstimateList$loglik
      remstimateList$converged        <- optimum_obj$converged
      remstimateList$iterations       <- optimum_obj$iterations
      remstimateList$sampled          <- TRUE
      remstimateList$samp_num         <- attr(stats, "samp_num")
      remstimateList$sampling_scheme  <- attr(stats, "sampling_scheme")

    } else { # HMC + sampled
      if (is.null(init)) {
        mle_obj <- trust::trust(
          objfun    = remDerivativesSampled,
          parinit   = rep(0.0, P),
          rinit     = 1, rmax = 100,
          stats     = stats_arr,
          case_pos  = case_pos_0,
          samp_prob = samp_prob,
          interevent_time = intereventTime,
          ordinal   = ordinal,
          ncores    = ncores
        )
        init <- mle_obj$argument
      }

      if (!is.null(seed)) set.seed(seed)

      prior_mean <- if (!is.null(prior$mean)) prior$mean else rep(0.0, P)
      prior_vcov <- if (!is.null(prior$vcov)) prior$vcov else diag(100, P)

      chain_draws <- matrix(NA_real_, nrow = nsim, ncol = P)
      chain_logpost <- numeric(nsim)

      beta_cur  <- init
      logpost_fn <- function(b) {
        ll <- -remDerivativesSampled(b, stats_arr, case_pos_0, samp_prob,
                                      intereventTime, ordinal,
                                      gradient = FALSE, hessian = FALSE,
                                      ncores = ncores)$value
        lp <- sum(mvnfast::dmvn(matrix(b, nrow = 1), prior_mean,
                                 prior_vcov, log = TRUE))
        ll + lp
      }
      grad_fn <- function(b) {
        g <- -remDerivativesSampled(b, stats_arr, case_pos_0, samp_prob,
                                     intereventTime, ordinal,
                                     gradient = TRUE, hessian = FALSE,
                                     ncores = ncores)$gradient
        gp <- solve(prior_vcov, b - prior_mean)
        g - gp
      }

      logpost_cur <- logpost_fn(beta_cur)

      for (iter in seq_len(nsim + burnin)) {
        # leapfrog
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
        K_cur  <- 0.5 * sum(r^2)
        K_prop <- 0.5 * sum(r_prop^2)

        log_accept <- logpost_prop - logpost_cur + K_cur - K_prop
        if (log(stats::runif(1)) < log_accept) {
          beta_cur   <- beta_prop
          logpost_cur <- logpost_prop
        }
        if (iter > burnin) {
          chain_draws[iter - burnin, ] <- beta_cur
          chain_logpost[iter - burnin]  <- logpost_cur
        }
      }

      # thin
      keep <- seq(thin, nsim, by = thin)
      draws <- chain_draws[keep, , drop = FALSE]
      colnames(draws) <- variables_names

      remstimateList$draws       <- draws
      remstimateList$log_posterior <- -chain_logpost[keep]
      remstimateList$coefficients  <- draws[which.min(remstimateList$log_posterior), ]
      remstimateList$post.mean     <- colMeans(draws)
      remstimateList$vcov          <- stats::cov(draws)
      remstimateList$sd            <- sqrt(diag(remstimateList$vcov))
      remstimateList$loglik        <- max(-remstimateList$log_posterior)
      names(remstimateList$coefficients) <- names(remstimateList$post.mean) <-
        names(remstimateList$sd) <- rownames(remstimateList$vcov) <-
        colnames(remstimateList$vcov) <- variables_names
      remstimateList$sampled      <- TRUE
    }
  }

  # ── Build output object ──────────────────────────────────────────────────────
  str_out <- structure(remstimateList, class = "remstimate")
  attr(str_out, "formula")          <- model_formula
  attr(str_out, "model")            <- model
  attr(str_out, "ordinal")          <- ordinal
  attr(str_out, "method")           <- method
  attr(str_out, "approach")         <- if (method == "MLE") "Frequentist" else "Bayesian"
  attr(str_out, "statistics")       <- variables_names
  attr(str_out, "where_is_baseline") <- where_is_baseline
  attr(str_out, "ncores")           <- ncores
  attr(str_out, "sampled")          <- sampled
  if (method == "HMC") {
    attr(str_out, "nsim")    <- nsim
    attr(str_out, "nchains") <- nchains
    attr(str_out, "burnin")  <- burnin
    attr(str_out, "thin")    <- thin
    attr(str_out, "L")       <- L
    attr(str_out, "epsilon") <- epsilon
    attr(str_out, "seed")    <- seed
    attr(str_out, "prior")   <- prior
  }

  return(str_out)
}

# ── Internal helpers ─────────────────────────────────────────────────────────

# NULL-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Build MLE output list for full stats path
.build_mle_output <- function(optimum_obj, variables_names, P, M, D,
                                reh, ordinal, stats_arr, dyad_field,
                                omit_dyad, intereventTime, ncores, WAIC) {
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

  # Null deviance: intercept-only model
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
  out$converged        <- optimum_obj$converged
  out$iterations       <- optimum_obj$iterations
  out$sampled          <- FALSE
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
  out$sampled      <- FALSE
  out
}
