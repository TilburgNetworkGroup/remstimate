

# ══════════════════════════════════════════════════════════════════════════════
# Bayesian random effects: brms
# ══════════════════════════════════════════════════════════════════════════════

.remstimate_brms <- function(reh, stats, random,
                             prior = NULL, nsim = 2000L,
                             nchains = 4L, burnin = 1000L,
                             thin = 1L, seed = NULL, ...) {
  if (!requireNamespace("brms", quietly = TRUE))
    stop("install.packages('brms')")

  s <- .remstimate_make_stack(reh, stats, add_actors = TRUE)
  df <- s$df
  stat_names <- s$stat_names

  fixed_part <- paste(stat_names, collapse = " + ")
  rand_part  <- deparse(random[[2]])

  if (!s$ordinal && "log_interevent" %in% names(df)) {
    fml <- stats::as.formula(paste0(
      "obs ~ ", fixed_part, " + ", rand_part, " + offset(log_interevent)"
    ))
    fam <- poisson()
  } else if (s$ordinal) {
    fml <- stats::as.formula(paste0("obs ~ ", fixed_part, " + ", rand_part))
    fam <- brms::bernoulli()
  } else {
    fml <- stats::as.formula(paste0("obs ~ ", fixed_part, " + ", rand_part))
    fam <- poisson()
  }

  if (is.null(seed)) seed <- NA
  fit <- brms::brm(
    formula = fml, family = fam, data = df,
    iter = nsim + burnin, warmup = burnin, chains = nchains,
    thin = thin, seed = seed, prior = prior, ...
  )

  coefs <- brms::fixef(fit)[, "Estimate"]

  .remstimate_wrap(
    coefficients = coefs,
    stat_names   = stat_names,
    loglik       = tryCatch(as.numeric(brms::log_lik(fit) |> rowSums() |> mean()),
                            error = function(e) NA_real_),
    stacked_data = df,
    backend_fit  = fit,
    model        = s$model,
    method       = "GLMM",
    engine       = "brms",
    ordinal      = s$ordinal
  )
}



# ── brms backend (HMC) ──────────────────────────────────────────────────────

#' @keywords internal
.remstimate_durem_brms <- function(stacked,
                                   prior   = NULL,
                                   nsim    = 500L,
                                   nchains = 1L,
                                   burnin  = 500L,
                                   thin    = 10L,
                                   seed    = NULL,
                                   ...) {

  if (!requireNamespace("brms", quietly = TRUE))
    stop("Package 'brms' is required for HMC estimation of duration models. ",
         "Install it with install.packages('brms').", call. = FALSE)

  df         <- stacked$remstats_stack
  stat_names <- stacked$stat_names

  if (isTRUE(stacked$ordinal))
    stop("Ordinal + HMC is not yet supported for duration models. ",
         "Use method = 'MLE' for ordinal durem (clogit backend).",
         call. = FALSE)

  if (length(stat_names) == 0L)
    stop("No statistics found — check start_effects / end_effects.",
         call. = FALSE)

  formula_obj <- stats::as.formula(paste0(
    "obs ~ -1 + offset(log_interevent) + ",
    paste(stat_names, collapse = " + ")
  ))

  if (is.null(prior))
    prior <- brms::set_prior("normal(0, 10)", class = "b")

  fit <- brms::brm(
    formula  = formula_obj,
    family   = stats::poisson(),
    data     = df,
    prior    = prior,
    iter     = nsim + burnin,
    warmup   = burnin,
    chains   = nchains,
    thin     = thin,
    seed     = seed,
    ...
  )

  draws <- brms::as_draws_matrix(fit)
  b_cols <- grep("^b_", colnames(draws))
  coef_draws <- as.matrix(draws[, b_cols])
  colnames(coef_draws) <- stat_names

  coefs    <- colMeans(coef_draws)
  vcov_mat <- stats::cov(coef_draws)
  sd_vec   <- sqrt(diag(vcov_mat))
  P        <- length(coefs)
  M        <- stacked$E

  names(coefs) <- names(sd_vec) <-
    rownames(vcov_mat) <- colnames(vcov_mat) <- stat_names

  res <- list(
    coefficients = coefs,
    post.mean    = coefs,
    vcov         = vcov_mat,
    sd           = sd_vec,
    loglik       = NULL,
    draws        = coef_draws,
    df.null      = M,
    df.model     = P,
    df.residual  = M - P,
    stacked_data = stacked,
    backend_fit  = fit
  )

  structure(
    res,
    class      = c("remstimate_durem", "remstimate"),
    formula    = formula_obj,
    model      = "tie",
    approach   = "Bayesian",
    method     = "HMC",
    engine     = "brms",
    ordinal    = FALSE,
    statistics = stat_names,
    ncores     = 1L
  )
}

