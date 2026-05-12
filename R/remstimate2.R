#' @title remstimate
#'
#' @description Estimates relational event model parameters. Supports
#'   tie-oriented and actor-oriented models via Maximum Likelihood
#'   (\code{"MLE"}), Hamiltonian Monte Carlo (\code{"HMC"}), mixed-effects
#'   Poisson regression (\code{"GLMM"}), penalised regression (\code{"GLMNET"}),
#'   and finite mixture models (\code{"MIXREM"}).
#'
#' @param reh     A \code{remify} object.
#' @param stats   A \code{tomstats}, \code{tomstats_sampled}, or
#'   \code{aomstats} object.
#' @param method  Estimation method. One of \code{"MLE"} (default),
#'   \code{"GLM"}, \code{"HMC"}, \code{"GLMM"}, \code{"GLMNET"}, or
#'   \code{"MIXREM"}. \code{"GLM"} fits via \code{glm()} — equivalent to
#'   MLE for interval timing, useful for comparison or for duration models.
#' @param ncores  [\emph{MLE/HMC}] Number of threads. Default \code{1L}.
#' @param nsim    [\emph{HMC}] MCMC iterations per chain after burnin. Default \code{500L}.
#' @param nchains [\emph{HMC}] Number of chains. Default \code{1L}.
#' @param burnin  [\emph{HMC}] Burn-in iterations. Default \code{500L}.
#' @param thin    [\emph{HMC}] Thinning interval. Default \code{10L}.
#' @param init    [\emph{HMC}] Initial parameter values; defaults to MLE estimates.
#' @param L       [\emph{HMC}] Leapfrog steps. Default \code{50L}.
#' @param epsilon [\emph{HMC}] Leapfrog step size. Default \code{0.002}.
#' @param prior   [\emph{HMC}] List with \code{mean} and \code{vcov}; defaults
#'   to standard normal.
#' @param seed    [\emph{HMC}] Random seed.
#' @param WAIC    Compute WAIC? Default \code{FALSE}.
#' @param nsimWAIC Number of posterior draws for WAIC. Default \code{100L}.
#' @param random  [\emph{GLMM / MIXREM}] One-sided formula for random
#'   effects or mixture clustering, e.g. \code{~ (1 | actor1)} or
#'   \code{~ (1 + inertia | dyad)}.  Grouping variables must be columns in
#'   the stacked data: \code{actor1} (sender), \code{actor2} (receiver),
#'   \code{dyad} (dyad index).
#' @param engine  [\emph{GLMM}] \code{"lme4"} (default) or \code{"glmmTMB"}.
#' @param alpha   [\emph{GLMNET}] Elastic-net mixing: \code{1} = lasso
#'   (default), \code{0} = ridge.
#' @param nfolds  [\emph{GLMNET}] CV folds. Default \code{10L}.
#' @param lambda_select [\emph{GLMNET}] Which lambda to use: \code{"1se"}
#'   (default) or \code{"min"}.
#' @param k       [\emph{MIXREM}] Number of mixture components (integer or
#'   vector). Default \code{2L}. When a vector, returns a
#'   \code{remstimate_mixrem_list}; use \code{\link{bic_table}} to compare.
#' @param concomitant [\emph{MIXREM}] Optional one-sided formula for the
#'   concomitant model, e.g. \code{~ actor1_degree}.
#' @param nrep    [\emph{MIXREM}] Random restarts for flexmix. Default \code{3L}.
#' @param ...     Further arguments passed to the backend.
#'
#' @return A \code{remstimate} S3 object. The subclass depends on the method:
#'   \code{remstimate} (MLE/HMC), \code{remstimate_glmm},
#'   \code{remstimate_glmnet}, or \code{remstimate_mixrem}.
#'   All objects have a \code{$coefficients} slot and support
#'   \code{summary()}, \code{coef()}, \code{diagnostics()}, and \code{plot()}.
#'
#' @examples
#' # ---- MLE, tie-oriented ----
#' data(tie_data)
#' tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
#' tie_model <- ~ 1 + remstats::indegreeSender() +
#'                    remstats::inertia() + remstats::reciprocity()
#' tie_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)
#' fit_mle <- remstimate(reh = tie_reh, stats = tie_stats, method = "MLE")
#' summary(fit_mle)
#'
#' # ---- MLE, actor-oriented ----
#' data(ao_data)
#' ao_reh <- remify::remify(edgelist = ao_data$edgelist, model = "actor")
#' ao_stats <- remstats::remstats(reh = ao_reh,
#'   sender_effects   = ~ 1 + remstats::indegreeSender(),
#'   receiver_effects = ~ remstats::inertia() + remstats::reciprocity())
#' fit_ao <- remstimate(reh = ao_reh, stats = ao_stats, method = "MLE")
#' summary(fit_ao)
#'
#' \donttest{
#' # ---- GLMM: sender frailty ----
#' data(history, package = "remstats")
#' data(info,    package = "remstats")
#' colnames(history)[colnames(history) == "setting"] <- "type"
#' reh <- remify::remify(edgelist = history[1:100, ], model = "tie",
#'                        riskset = "active")
#' effects <- ~ inertia(consider_type = FALSE) +
#'               indegreeSender(consider_type = FALSE)
#' stats <- remstats::tomstats(effects, reh = reh, attr_actors = info,
#'                              memory = "decay", memory_value = 1000)
#'
#' fit_glmm <- remstimate(reh, stats, method = "GLMM",
#'                         random = ~ (1 | actor1) + (1 | actor2))
#' fit_glmm
#' lme4::ranef(fit_glmm$backend_fit)
#'
#' # ---- GLMNET: lasso ----
#' fit_lasso <- remstimate(reh, stats, method = "GLMNET", alpha = 1)
#' fit_lasso
#' plot(fit_lasso, reh, stats = stats, which = 6)   # regularisation path
#'
#' # ---- MIXREM: two-component mixture ----
#' fit_mix <- remstimate(reh, stats, method = "MIXREM",
#'                        random = ~ (1 + inertia | dyad), k = 2)
#' fit_mix
#'
#' # compare k = 2:4 via BIC
#' fits <- remstimate(reh, stats, method = "MIXREM",
#'                    random = ~ (1 + inertia | dyad), k = 2:4)
#' bic_table(fits)
#' }
#'
#' @export
remstimate <- function(reh,
                       stats,
                       method    = c("MLE", "GLM", "HMC", "GLMM", "GLMNET", "MIXREM"),
                       ncores    = 1L,
                       nsim      = 500L,
                       nchains   = 1L,
                       burnin    = 500L,
                       thin      = 10L,
                       init      = NULL,
                       L         = 50L,
                       epsilon   = 0.002,
                       prior     = NULL,
                       seed      = NULL,
                       WAIC      = FALSE,
                       nsimWAIC  = 100L,
                       # GLMM
                       random    = NULL,
                       engine    = c("lme4", "glmmTMB"),
                       # GLMNET
                       alpha         = 1,
                       nfolds        = 10L,
                       lambda_select = c("1se", "min"),
                       # MIXREM
                       k           = 2L,
                       concomitant = NULL,
                       nrep        = 3L,
                       ...) {

  method <- toupper(match.arg(method))

  # ── Duration REM: MLE/HMC need special backends ───────────────────────────
  if (inherits(reh, "remify_durem") & inherits(stats, "remstats_durem")) {
    if (method %in% c("MLE", "GLM")) {
      stacked <- stack_stats(stats, reh, add_actors = FALSE)
      return(.remstimate_durem_glm(stacked))
    }
    if (method == "HMC") {
      stacked <- stack_stats(stats, reh, add_actors = FALSE)
      return(.remstimate_durem_brms(stacked, prior = prior, nsim = nsim,
                                    nchains = nchains, burnin = burnin,
                                    thin = thin, seed = seed, ...))
    }
    # GLMM / GLMNET / MIXREM fall through — .remstimate_make_stack()
    # already handles remstats_durem objects via stack_stats dispatch
  }

  # ── R backends (any model type including durem) ────────────────────────────
  if (method == "GLMM")
    return(.remstimate_glmm(reh, stats, random = random,
                            engine = match.arg(engine), ...))
  if (method == "GLMNET")
    return(.remstimate_glmnet(reh, stats, alpha = alpha,
                              nfolds = nfolds,
                              lambda_select = match.arg(lambda_select), ...))
  if (method == "MIXREM")
    return(.remstimate_mixrem(reh, stats, random = random, k = k,
                              concomitant = concomitant, nrep = nrep, ...))

  # fit basic REM with tailored fitting algorithms when either method = "MLE" or "HMC"
  if(!is.numeric(nsimWAIC) || length(nsimWAIC) != 1 || nsimWAIC < 1){
    warning("'nsimWAIC' must be a positive integer. Using default value of 100.")
    nsimWAIC <- 100L
  }
  nsimWAIC <- as.integer(nsimWAIC)

  if(!is.logical(WAIC) || length(WAIC) != 1){
    stop("'WAIC' must be a single logical value (TRUE or FALSE).")
  }

  if(!is.numeric(ncores) || ncores < 1L){
    warning("'ncores' must be a positive integer. Using default value of 1.")
    ncores <- 1L
  }
  ncores <- as.integer(ncores)

  # ── Input checks ────────────────────────────────────────────────────────────
  if (!inherits(reh, "remify")) stop("'reh' must be a 'remify' object.")

  sampled <- inherits(stats, "tomstats_sampled")
  if (!inherits(stats, c("tomstats", "tomstats_sampled", "aomstats"))) {
    stop("'stats' must be a 'tomstats', 'tomstats_sampled', or 'aomstats' object.")
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

  if(inherits(stats, "aomstats") && model == "tie"){
    stop("'stats' is an 'aomstats' object but 'reh' is a tie-oriented model.")
  }
  if(inherits(stats, c("tomstats","tomstats_sampled")) && model == "actor"){
    stop("'stats' is a 'tomstats' object but 'reh' is an actor-oriented model.")
  }

  if (!model %in% c("tie", "actor")) {
    stop("remstimate2 currently supports 'tie' and 'actor' models only.")
  }
  if (is.null(ncores) || ncores < 1L) ncores <- ncores_reh

  # ... directed / undirected network
  if((model == "actor") & !directed){stop("actor-oriented modeling can't operate on undirected networks")}

  # ── Clamp ncores ────────────────────────────────────────────────────────────
  max_cores <- parallel::detectCores()
  if (max_cores <= 2L) ncores <- 1L else ncores <- min(ncores, max_cores - 2L)

  # ── Extract dyad IDs and interevent times ────────────────────────────────────
  # New remify2 format: ids$dyad_active (active riskset) or ids$dyad
  # Old remify format:  attr(reh, "dyadID") / attr(reh, "dyadIDactive")
  if (!is.null(reh$ids)) {
    # New format — dyad_active is a [N x 1] matrix, flatten to vector
    if (riskset %in% c("active", "manual")) {
      dyadID <- as.vector(reh$ids$dyad_active)
    } else {
      dyadID <- as.vector(reh$ids$dyad)
    }
    intereventTime <- reh$intereventTime
    simultaneous_idx <- reh$indices_simultaneous_events
    evenly_spaced <- if (!is.null(reh$simultaneous)) reh$simultaneous$interevent_evenly_spaced else NULL
  } else {
    if (riskset %in% c("active", "manual")) {
      dyadID <- attr(reh, "dyadIDactive") %||% attr(reh, "dyadID")
    } else {
      dyadID <- attr(reh, "dyadID")
    }
    intereventTime <- reh$intereventTime
    simultaneous_idx <- attr(reh, "indices_simultaneous_events")
    evenly_spaced <- attr(reh, "evenly_spaced_interevent_time")
  }

  if(method == "HMC"){
    # ... seed
    if(is.null(seed)){
      seed <- .Random.seed
    }
    else if(length(seed) > 1){
      seed <- seed[1]
      warning("`seed` length is greater than 1. Considering only the first element")
    }

    # ... nchains
    if(is.null(nchains)){
      nchains <- 1L # setting default value of 'nchains'
    }

    # ... nsim
    if(is.null(nsim)){
      nsim <- 1e03L # setting default value of 'nsim'
    }

    # ... burnin
    if(is.null(burnin)){
      burnin <- 500L # setting default value of 'burnin'
    }

    # ... thin
    if(is.null(thin)){
      if(nsim >= 100){
        thin <- 10L
        warning("'thin' parameter undefined. Using thin = 10")
      }
      else{
        thin <- 1L # no thinning is applied
        warning("'nsim' is less than 100. No thinning is applied")
      }
    }
    else if(is.numeric(thin)){
      if(thin <= 0){
        stop("`thin` value must be positive. Set thin = 1 if no thinning of chains is required")
      }
      if(thin == 1){
        warning("`thin` value is 1. No thinning applied to chains")
      }
    }

    # ... L (number of leap-frog steps)
    if(is.null(L)){
      L <- 50L # setting default value of 'L'
    }
    else if(length(L)>1){
      L <- 50L
      warning("input 'L' must be a positive (integer) number . 'L' is set to its default value: 50")
    }
    else if((L<=1) | !is.numeric(L)){
      L <- 50L
      warning("input 'L' must be a positive (integer) number . 'L' is set to its default value: 50")
    }
    # ... epsilon (sie of time step in the leap-frog algorithm)
    if(is.null(epsilon)){
      epsilon <- 0.1/L # 0.1 is the full time, reached by steps of 0.1/L
    }
    else if(length(epsilon)>1){
      epsilon <- 0.1/L
      warning("input 'epsilon' must be a positive number. 'epsilon' is set to its default value: 0.1/L")
    }
    else if((epsilon <= 0) | !is.numeric(epsilon)){
      epsilon <- 0.1/L
      warning("input 'epsilon' must be a positive number. 'epsilon' is set to its default value: 0.1/L")
    }
    # 0.1 is a [temporary parameter] and it will change in the future releases of 'remstimate'
  }

  # ── D: riskset size ──────────────────────────────────────────────────────────
  D <- if (riskset %in% c("active", "manual")) (reh$activeD %||% reh$D) else reh$D
  M <- reh$M

  # ── stats method (pt / pe) ───────────────────────────────────────────────────
  stats_method <- attr(stats, "method")
  if (is.null(stats_method)) stop("'stats' has no 'method' attribute.")

  if (stats_method == "pe") {
    if (!is.null(evenly_spaced)) intereventTime <- as.vector(evenly_spaced)
    if (!is.null(reh$E)) M <- reh$E
  }

  # if (ordinal) intereventTime <- rep(1.0, M_sub)

  # ── AOM branch -------------------───────────────────────────────────────────
  if (model == "actor") {

    # Replace the current single start_stop extraction with:
    start_stop <- as.numeric(unlist(attr(stats,"subset")))
    # start_stop_sender   <- as.integer(unlist(subset_attr$sender))
    # start_stop_receiver <- as.integer(unlist(subset_attr$receiver))

    actor1_ids <- reh$ids$actor1[start_stop[1]:start_stop[2]]
    actor2_ids <- reh$ids$actor2[start_stop[1]:start_stop[2]]

    M_sub_sender <- length(actor1_ids)
    M_sub_receiver <- length(actor2_ids)

    intereventTime_sender   <- if (ordinal) rep(1.0, M_sub_sender)   else reh$intereventTime[start_stop[1]:start_stop[2]]
    intereventTime_receiver <- rep(1.0, M_sub_receiver)

    # # Extract actor IDs from new remify2 structure
    # actor1_ids <- reh$ids$actor1[start_stop[1]:start_stop[2]]  # 1-based, keep as list-like for field<uvec>
    # actor2_ids <- reh$ids$actor2[start_stop[1]:start_stop[2]]  # unlisted for receiver model
    #
    N  <- reh$N
    # M_sub <- diff(start_stop) + 1L
    #
    # if (!ordinal) {
    #   intereventTime <- reh$intereventTime[start_stop[1]:start_stop[2]]
    # } else {
    #   intereventTime <- rep(1.0, M_sub)
    # }

    sender_rs   <- if (!is.null(reh$sender_riskset))   as.integer(reh$sender_riskset) - 1L   else integer(0)
    receiver_rs <- if (!is.null(reh$receiver_riskset)) lapply(reh$receiver_riskset, function(x) as.integer(x) - 1L) else list()

    # ── omit_dyad not used because vectorized baniary version could be huge for large N
    omit_dyad_sender   <- reh$omit_dyad %||% list()
    omit_dyad_receiver <- list()

    which_stats  <- c("sender_stats",   "receiver_stats")
    which_model  <- c("sender_model",   "receiver_model")
    sender_rate  <- c(TRUE,             FALSE)

    remstimateList <- list()

    for (i in 1:2) {
      s <- stats[[ which_stats[i] ]]
      if (is.null(s)) next

      # Reshape [M x N x P] -> [N x P x M]
      s_arr <- aperm(s, perm = c(2, 3, 1))

      P_i         <- dim(s_arr)[2]
      var_names_i <- dimnames(s)[[3]]

      if (sender_rate[i]) {
        intereventTime  <- intereventTime_sender
        M_sub           <- M_sub_sender
        actor1_field    <- lapply(actor1_ids, as.integer)
        actor2_field    <- lapply(actor2_ids, as.integer)
      } else {
        intereventTime  <- intereventTime_receiver
        M_sub           <- M_sub_receiver
        actor1_field    <- lapply(actor1_ids, function(x) as.integer(unlist(x)))
        actor2_field    <- lapply(actor2_ids, function(x) as.integer(unlist(x)))
      }

      omit_dyad_i <- if (sender_rate[i]) omit_dyad_sender else omit_dyad_receiver

      if (method == "MLE") {
        opt <- trust::trust(
          objfun          = remDerivatives,
          parinit         = rep(0.0, P_i),
          rinit           = 1,
          rmax            = 100,
          stats           = s_arr,
          actor1          = actor1_field,
          actor2          = actor2_field,
          dyad            = list(),
          omit_dyad       = omit_dyad_i,
          interevent_time = intereventTime,
          model           = "actor",
          ordinal         = ordinal,
          ncores          = ncores,
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

        # null deviance
        null_dim <- dim(s_arr)
        null_dim[2] <- 1L
        res$null.deviance <- if (ordinal) {
          risk_set_size <- if (sender_rate[i]) N else N - 1L
          2 * M_sub * log(risk_set_size)
        } else {
          2 * trust::trust(
            objfun          = remDerivatives,
            parinit         = 0.0,
            rinit           = 1, rmax = 100,
            stats           = array(1.0, dim = null_dim),
            actor1          = actor1_field,
            actor2          = actor2_field,
            dyad            = list(),
            omit_dyad       = omit_dyad_i,
            interevent_time = intereventTime,
            model           = "actor",
            ordinal         = ordinal,
            ncores          = ncores,
            senderRate      = sender_rate[i],
            N               = N
          )$value
        }
        res$model.deviance <- res$null.deviance - res$residual.deviance

        if(WAIC){
          res$WAIC <- getWAIC(
            mu              = res$coefficients,
            vcov            = res$vcov,
            pars            = matrix(0, 2, 2),
            stats           = s_arr,
            actor1          = actor1_field,
            actor2          = actor2_field,
            dyad            = list(),
            omit_dyad       = omit_dyad_i,
            interevent_time = intereventTime,
            model           = "actor",
            approach        = "Frequentist",
            nsim            = nsimWAIC,
            ordinal         = ordinal,
            ncores          = ncores,
            senderRate      = sender_rate[i],
            sender_riskset   = if (i == 1) sender_rs else integer(0),
            receiver_riskset = if (i == 2) receiver_rs else list()
          )
        }

        remstimateList[[ which_model[i] ]] <- res

      } else { # HMC
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
          ordinal         = ordinal,
          ncores          = ncores,
          senderRate      = sender_rate[i],
          N               = N,
          sender_riskset = if(i==1) sender_rs else integer(0),
          receiver_riskset = if(i==2) receiver_rs else list()
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
    }

    # build output
    str_out <- structure(remstimateList, class = "remstimate")
    attr(str_out, "model")   <- "actor"
    attr(str_out, "ordinal") <- ordinal
    attr(str_out, "method")  <- method
    attr(str_out, "approach") <- if (method == "MLE") "Frequentist" else "Bayesian"
    attr(str_out, "formula") <- list(
      rate_model_formula   = attr(stats, "formula")$rate,
      choice_model_formula = attr(stats, "formula")$choice
    )
    attr(str_out, "statistics") <- list(
      sender_model   = dimnames(stats$sender_stats)[[3]],
      receiver_model = dimnames(stats$receiver_stats)[[3]]
    )
    attr(str_out, "ncores") <- ncores
    if (method == "HMC") {
      attr(str_out, "nsim")    <- nsim
      attr(str_out, "nchains") <- nchains
      attr(str_out, "burnin")  <- burnin
      attr(str_out, "thin")    <- thin
      attr(str_out, "L")       <- L
      attr(str_out, "epsilon") <- epsilon
      attr(str_out, "seed")    <- seed
    }
    return(str_out)
  }

  if (model == "tie"){

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
    if (ordinal || is.null(intereventTime)) intereventTime <- rep(1.0, M_sub)

    # ── omit_dyad ────────────────────────────────────────────────────────────────
    omit_dyad <- list() #omit_dyad not used in remstimate2

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
          stats_arr, dyad_field, omit_dyad, intereventTime, ncores, WAIC,
          nsimWAIC, model
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
          nsim            = as.integer(nsim+burnin),
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
        remstimateList$df.null <- M_sub
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
        remstimateList$df.null       <- M_sub
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
  }

  return(str_out)
}

# ── Internal helpers ─────────────────────────────────────────────────────────

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

  # WAIC -- this block was missing
  if(WAIC){
    out$WAIC <- getWAIC(
      mu              = out$coefficients,
      vcov            = out$vcov,
      pars            = matrix(0, 2, 2),
      stats           = stats_arr,
      actor1          = list(),
      actor2          = list(),
      dyad            = dyad_field,
      omit_dyad       = omit_dyad,
      interevent_time = intereventTime,
      model           = model,
      approach        = "Frequentist",
      nsim            = nsimWAIC,
      ordinal         = ordinal,
      ncores          = ncores
    )
  }

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
