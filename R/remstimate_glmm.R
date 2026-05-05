# GLMM backend — lme4::glmer or glmmTMB::glmmTMB
#
# remstimate(reh, stats, method = "GLMM", random = ~ (1|actor1), engine = "lme4")

.remstimate_glmm <- function(reh, stats,
                              random  = NULL,
                              engine  = c("lme4", "glmmTMB"),
                              verbose = FALSE,
                              ...) {

  engine <- match.arg(engine)
  if (!requireNamespace(engine, quietly = TRUE))
    stop("install.packages('", engine, "')")
  if (is.null(random))
    stop("provide a random formula, e.g. ~ (1 | actor1)")

  s <- .remstimate_make_stack(reh, stats, add_actors = TRUE)

  if (s$model == "tie") {

    fit   <- .glmm_fit_one(s$df, s$stat_names, s$ordinal, random, engine, verbose, ...)
    coefs <- .glmm_fixef(fit, engine)

    .remstimate_wrap(
      coefficients = coefs,
      stat_names   = s$stat_names,
      loglik       = .glmm_loglik(fit),
      stacked_data = s$df,
      backend_fit  = fit,
      model        = "tie",
      method       = "GLMM",
      engine       = engine,
      ordinal      = s$ordinal
    )

  } else {

    fit_s <- fit_r <- NULL
    if (!is.null(s$df$sender))
      fit_s <- .glmm_fit_one(s$df$sender, s$stat_names$sender_model,
                              s$ordinal, random, engine, verbose, ...)
    if (!is.null(s$df$receiver))
      fit_r <- .glmm_fit_one(s$df$receiver, s$stat_names$receiver_model,
                              s$ordinal, random, engine, verbose, ...)

    .remstimate_wrap(
      coefficients = list(sender_model   = .glmm_fixef(fit_s, engine),
                          receiver_model = .glmm_fixef(fit_r, engine)),
      stat_names   = s$stat_names,
      loglik       = c(sender = .glmm_loglik(fit_s), receiver = .glmm_loglik(fit_r)),
      stacked_data = s$df,
      backend_fit  = list(sender_model = fit_s, receiver_model = fit_r),
      model        = "actor",
      method       = "GLMM",
      engine       = engine,
      ordinal      = s$ordinal
    )
  }
}

.glmm_fit_one <- function(df, stat_names, ordinal, random, engine, verbose, ...) {
  vaste_kant <- .remstimate_fixed_rhs(stat_names, ordinal)
  familie    <- if (ordinal) {
    df$event  <- factor(df$event)
    vaste_kant <- sub("^-1 \\+ ", "-1 + event + ", vaste_kant)
    message("ordinal: binomial GLMM with event fixed effects (approx. clogit)")
    binomial()
  } else poisson()

  fml <- as.formula(paste("obs ~", vaste_kant, "+", deparse(random[[2]])))
  if (verbose) message("fit: ", deparse(fml))

  if (engine == "lme4")
    lme4::glmer(fml, family = familie, data = df, ...)
  else
    glmmTMB::glmmTMB(fml, family = familie, data = df, ...)
}

.glmm_fixef <- function(fit, engine) {
  if (is.null(fit)) return(NULL)
  if (engine == "lme4")    return(lme4::fixef(fit))
  if (engine == "glmmTMB") return(glmmTMB::fixef(fit)$cond)
  coef(fit)
}

.glmm_loglik <- function(fit) {
  if (is.null(fit)) return(NULL)
  tryCatch(as.numeric(logLik(fit)), error = function(e) NULL)
}

#' @export
print.remstimate_glmm <- function(x, ...) {
  engine <- attr(x, "engine")
  cat("REM —", attr(x, "model"), "model — GLMM [", engine, "]\n\n")

  if (attr(x, "model") == "tie") {
    cat("Fixed effects:\n"); print(x$coefficients)
    cat("\nRandom effects:\n")
    tryCatch({
      vc <- if (engine == "lme4") lme4::VarCorr(x$backend_fit) else
            glmmTMB::VarCorr(x$backend_fit)
      print(vc, comp = "Variance")
    }, error = function(e) NULL)
    if (!is.null(x$loglik)) cat("\nlogLik:", round(x$loglik, 4), "\n")
    if (engine == "lme4" && isTRUE(tryCatch(lme4::isSingular(x$backend_fit), error = function(e) FALSE)))
      message("singular fit — consider simplifying the random structure")
  } else {
    for (deel in c("sender_model", "receiver_model")) {
      if (is.null(x[[deel]])) next
      cat(deel, "— fixed effects:\n"); print(x[[deel]]$coefficients); cat("\n")
    }
  }
  invisible(x)
}

#' @export
summary.remstimate_glmm <- function(object, ...) {
  if (attr(object, "model") == "tie") summary(object$backend_fit, ...)
  else lapply(list(object$sender_model$backend_fit,
                   object$receiver_model$backend_fit), summary)
}

#' @export
coef.remstimate_glmm   <- function(object, ...) object$coefficients
#' @export
logLik.remstimate_glmm <- function(object, ...) object$loglik
