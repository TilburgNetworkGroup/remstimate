# shared helpers for external remstimate backends (GLMM, GLMNET, MIXREM)

.remstimate_make_stack <- function(reh, stats, add_actors = TRUE) {

  model     <- if (inherits(stats, "aomstats")) "actor" else "tie"
  gestapeld <- remstats::stack_stats(stats, reh, add_actors = add_actors)

  if (model == "tie") {
    df <- gestapeld$remstats_stack
    # sampling correction: offset = -log(π_d) per row
    # for uniform case-control sampling this upweights each control by (D-1)/K
    if (inherits(stats, "tomstats_sampled")) {
      sp_mat         <- attr(stats, "samp_prob")   # [M × samp_num], case col = 1
      df$samp_offset <- as.vector(t(-log(sp_mat))) # event-major flatten
    } else {
      df$samp_offset <- rep(0, nrow(df))
    }
    list(df = df, stat_names = dimnames(stats)[[3]], ordinal = gestapeld$ordinal, model = model)
  } else {
    df_s <- gestapeld$sender_stack
    df_r <- gestapeld$receiver_stack
    if (!is.null(df_s)) df_s$samp_offset <- rep(0, nrow(df_s))
    if (!is.null(df_r)) df_r$samp_offset <- rep(0, nrow(df_r))
    list(
      df = list(sender = df_s, receiver = df_r),
      stat_names = list(
        sender_model   = if (!is.null(df_s)) dimnames(stats$sender_stats)[[3]] else character(0),
        receiver_model = if (!is.null(df_r)) dimnames(stats$receiver_stats)[[3]] else character(0)
      ),
      ordinal = gestapeld$ordinal,
      model   = model
    )
  }
}

# rhs string for glmer / glmmTMB / flexmix formula
# interval timing: Poisson with log_interevent + samp_offset
# ordinal timing:  binomial, no time offset but sampling correction still applies
.remstimate_fixed_rhs <- function(stat_names, ordinal) {
  rhs <- paste(c("-1", stat_names), collapse = " + ")
  if (!ordinal)
    paste(rhs, "+ offset(log_interevent + samp_offset)")
  else
    paste(rhs, "+ offset(samp_offset)")
}

# model matrix for glmnet
.remstimate_model_matrix <- function(df, stat_names, ordinal) {
  list(
    X       = as.matrix(df[, stat_names, drop = FALSE]),
    y       = df$obs,
    offset  = if (!ordinal) df$log_interevent + df$samp_offset else df$samp_offset,
    weights = if ("weight" %in% names(df)) df$weight else NULL
  )
}

.remstimate_find_baseline <- function(stat_names) {
  idx <- which(tolower(stat_names) == "baseline")
  if (length(idx)) idx[[1L]] else NULL
}

# build the return object; actor model mirrors MLE structure so
# diagnostics.remstimate works via inheritance without changes
.remstimate_wrap <- function(coefficients,
                              stat_names,
                              loglik       = NULL,
                              formula      = NULL,
                              stacked_data = NULL,
                              backend_fit  = NULL,
                              model        = "tie",
                              method       = "GLMM",
                              engine       = NA_character_,
                              ordinal      = FALSE,
                              extra        = list()) {

  if (model == "actor") {
    kern <- list(
      sender_model = list(
        coefficients = coefficients$sender_model,
        backend_fit  = if (is.list(backend_fit)) backend_fit$sender_model else backend_fit
      ),
      receiver_model = list(
        coefficients = coefficients$receiver_model,
        backend_fit  = if (is.list(backend_fit)) backend_fit$receiver_model else NULL
      ),
      loglik       = loglik,
      formula      = formula,
      stacked_data = stacked_data
    )
    waar_is_baseline <- .remstimate_find_baseline(stat_names$sender_model)
  } else {
    kern <- list(
      coefficients = coefficients,
      loglik       = loglik,
      formula      = formula,
      stacked_data = stacked_data,
      backend_fit  = backend_fit
    )
    waar_is_baseline <- .remstimate_find_baseline(stat_names)
  }

  structure(
    c(kern, extra),
    class             = c(paste0("remstimate_", tolower(method)), "remstimate"),
    model             = model,
    approach          = "Frequentist",
    method            = method,
    engine            = engine,
    ordinal           = ordinal,
    statistics        = stat_names,
    where_is_baseline = waar_is_baseline,
    ncores            = 1L
  )
}

# extract grouping vars from e.g. ~ (1|actor1) + (1|actor2)
.parse_random_groups <- function(random) {
  if (is.null(random)) return(character(0))
  rhs  <- deparse(random[[2]])
  hits <- gregexpr("\\|\\s*([A-Za-z_.][A-Za-z_.0-9]*)", rhs, perl = TRUE)
  trimws(sub("^\\|\\s*", "", regmatches(rhs, hits)[[1]]))
}
