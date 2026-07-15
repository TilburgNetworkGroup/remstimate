# GLMM backend - lme4::glmer or glmmTMB::glmmTMB
#
# remstimate(reh, stats, method = "GLMM", random = ~ (1|actor1), engine = "lme4")

.remstimate_glmm <- function(reh, stats,
                             random  = NULL,
                             engine  = c("glmmTMB", "lme4"),
                             verbose = FALSE,
                             ...) {
  engine <- match.arg(engine)

  s <- .remstimate_make_stack(reh, stats, add_actors = TRUE)

  # Ordinal models need coxme (conditional logit + random effects)
  if (s$ordinal && engine != "coxme") {
    message(sprintf(
      "Ordinal model: '%s' cannot fit conditional logit with random effects; using 'coxme' instead.",
      engine))
    engine <- "coxme"
  }

  if (!requireNamespace(engine, quietly = TRUE))
    stop("install.packages('", engine, "')")
  if (is.null(random))
    stop("provide a random formula, e.g. ~ (1 | actor1)")

  # Undirected tie models: actor-level frailty is not identified. The actor1 /
  # actor2 columns are an arbitrary canonical dyad ordering, not sender /
  # receiver, so a single actor's sociality gets split across two position
  # labels (two variance components, and actors at the ordering extremes appear
  # in only one column). Dyad-level frailty '(1 | dyad)' is symmetric by
  # construction and remains allowed.
  directed <- if (!is.null(reh$meta)) isTRUE(reh$meta$directed)
              else isTRUE(attr(reh, "directed"))
  if (s$model == "tie" && !directed) {
    actor_grp <- .glmm_actor_grouping(random)
    if (length(actor_grp)) {
      stop("Actor-level random effects (", paste(unique(actor_grp), collapse = ", "),
           ") are not identified for an undirected tie-oriented model: 'actor1' / ",
           "'actor2' are an arbitrary dyad ordering, not sender / receiver, so a ",
           "single actor effect is split across two position labels. Use dyad-level ",
           "frailty 'random = ~ (1 | dyad)' instead, which is symmetric by ",
           "construction.", call. = FALSE)
    }
  }

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

    # The receiver sub-model is a conditional choice: exactly one receiver is
    # picked per event from the sender's candidate set (the receiver stack is a
    # choice-set layout, stratified by time_index, with no timing offset). Its
    # likelihood is a conditional logit, NOT a Poisson rate, so it is ALWAYS fit
    # with coxme + strata(time_index) regardless of interval/ordinal timing.
    # Only the sender rate model follows the interval (Poisson) / ordinal (coxme)
    # split. An unstratified Poisson here would estimate marginal selection
    # rates, not conditional choice probabilities.
    engine_s <- engine
    engine_r <- "coxme"

    fit_s <- fit_r <- NULL
    if (!is.null(s$df$sender))
      fit_s <- .glmm_fit_one(s$df$sender, s$stat_names$sender_model,
                             s$ordinal, random$sender %||% random, engine_s, verbose, ...)
    if (!is.null(s$df$receiver)) {
      if (!requireNamespace("coxme", quietly = TRUE))
        stop("the receiver choice model is a conditional logit and requires ",
             "coxme: install.packages('coxme')", call. = FALSE)
      fit_r <- .glmm_fit_one(s$df$receiver, s$stat_names$receiver_model,
                             s$ordinal, random$receiver %||% random, engine_r, verbose, ...)
    }

    .remstimate_wrap(
      coefficients = list(sender_model   = .glmm_fixef(fit_s, engine_s),
                          receiver_model = .glmm_fixef(fit_r, engine_r)),
      stat_names   = s$stat_names,
      loglik       = c(sender = .glmm_loglik(fit_s), receiver = .glmm_loglik(fit_r)),
      stacked_data = s$df,
      backend_fit  = list(sender_model = fit_s, receiver_model = fit_r),
      model        = "actor",
      method       = "GLMM",
      engine       = engine_s,
      ordinal      = s$ordinal,
      extra        = list(engine_receiver = engine_r)
    )
  }
}

.glmm_fit_one <- function(df, stat_names, ordinal, random, engine, verbose, ...) {
  if (engine == "coxme") {
    if (!requireNamespace("coxme", quietly = TRUE))
      stop("install.packages('coxme')")
    df$.surv_time <- rep(1, nrow(df))
    # coxme cannot parse `a:b` interaction grouping; build real factor columns
    rterms <- paste(deparse(random[[2]]), collapse = " ")
    grp <- trimws(sub("\\|", "", regmatches(rterms, gregexpr("\\|[^)]+", rterms))[[1]]))
    for (g in grp) if (grepl(":", g, fixed = TRUE)) {
      parts   <- trimws(strsplit(g, ":", fixed = TRUE)[[1]])
      newname <- paste(parts, collapse = "_")
      df[[newname]] <- interaction(df[parts], drop = TRUE)
      rterms <- gsub(g, newname, rterms, fixed = TRUE)
    }
    formula_obj <- stats::as.formula(paste0(
      "survival::Surv(.surv_time, obs) ~ ",
      paste(stat_names, collapse = " + "),
      " + ", rterms,
      " + strata(time_index)"
    ))
    # coxme (the only engine able to fit conditional-logit random effects for
    # ordinal timing) is fragile on the many-strata REM design and can return a
    # non-converged fit with runaway coefficients. Capture the convergence
    # warning and flag it so it is surfaced rather than silently trusted.
    conv_msg <- NULL
    fit <- withCallingHandlers(
      coxme::coxme(formula_obj, data = df),
      warning = function(w) {
        if (grepl("did not converge|Ran out of iterations|iteration",
                  conditionMessage(w), ignore.case = TRUE)) {
          conv_msg <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        }
      }
    )
    attr(fit, "remstimate_converged") <- is.null(conv_msg)
    if (!is.null(conv_msg)) attr(fit, "remstimate_conv_msg") <- conv_msg
    fit
  } else {
    # Build formula: obs ~ -1 + stat1 + stat2 + ... + (1 | actor1) + offset(log_interevent)
    # Suppress R's implicit intercept: the REM parameterisation carries its own
    # 'baseline' column (all ones) when an intercept is requested, so an implicit
    # intercept would duplicate it (rank-deficient -> baseline dropped to NA in the
    # sender model) and would invent a spurious intercept in the receiver choice
    # model, which has none. With -1, 'baseline' is the intercept where present.
    fixed_part <- paste(c("-1", stat_names), collapse = " + ")
    rand_part  <- deparse(random[[2]])
    if (!ordinal && "log_interevent" %in% names(df)) {
      fml <- stats::as.formula(paste0(
        "obs ~ ", fixed_part, " + ", rand_part, " + offset(log_interevent)"
      ))
    } else {
      fml <- stats::as.formula(paste0(
        "obs ~ ", fixed_part, " + ", rand_part
      ))
    }
    familie <- if (ordinal) stats::binomial() else stats::poisson()

    if (engine == "lme4") {
      # lme4 fits are S4 (no custom attributes); a boundary/singular fit is
      # reported as an informational message in print() via isSingular().
      lme4::glmer(fml, family = familie, data = df, ...)
    } else {
      # Muffle nlminb's intermediate "NA/NaN function evaluation" chatter: it
      # fires on recoverable optimiser steps and is uninformative on its own
      # (the same message accompanies both harmless hiccups and failed fits).
      # Whether the fit is usable is decided from the FINAL status below.
      fit <- withCallingHandlers(
        glmmTMB::glmmTMB(fml, family = familie, data = df, ...),
        warning = function(w) {
          if (grepl("NA/NaN function evaluation", conditionMessage(w), fixed = TRUE))
            invokeRestart("muffleWarning")
        }
      )
      h <- .glmm_health_glmmTMB(fit)
      attr(fit, "remstimate_converged") <- h$ok
      if (!h$ok) attr(fit, "remstimate_conv_msg") <- h$msg
      fit
    }
  }
}

# glmmTMB post-fit health: usable only if the optimiser returned convergence
# code 0 AND the Hessian is positive definite (so standard errors are
# trustworthy). This is what distinguishes a muffled "NA/NaN function
# evaluation" that was a harmless step from one that signalled a failed fit.
.glmm_health_glmmTMB <- function(fit) {
  code <- tryCatch(fit$fit$convergence, error = function(e) NA_integer_)
  pdh  <- tryCatch(isTRUE(fit$sdr$pdHess), error = function(e) FALSE)
  if (isTRUE(code == 0L) && pdh) return(list(ok = TRUE, msg = NULL))
  parts <- c(
    if (!isTRUE(code == 0L)) "optimiser did not converge",
    if (!pdh) "non-positive-definite Hessian (standard errors unreliable)"
  )
  list(ok = FALSE, msg = paste(parts, collapse = "; "))
}

# Return the actor-level grouping factors (actor / actor1 / actor2, possibly
# crossed as actor1:process:type) referenced in a random spec. `random` may be
# a one-sided formula (tie model) or a list(sender=, receiver=) (actor model).
.glmm_actor_grouping <- function(random) {
  fl <- if (inherits(random, "formula")) list(random) else random
  hits <- character(0)
  for (f in fl) {
    if (is.null(f) || !inherits(f, "formula")) next
    rhs <- paste(deparse(f[[length(f)]]), collapse = " ")
    grp <- trimws(sub("\\|", "", regmatches(rhs, gregexpr("\\|[^)]+", rhs))[[1]]))
    for (g in grp) {
      comps <- trimws(strsplit(g, ":", fixed = TRUE)[[1]])
      if (any(comps %in% c("actor", "actor1", "actor2"))) hits <- c(hits, g)
    }
  }
  unique(hits)
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

# Surface a non-converged backend fit (coxme or glmmTMB; flag set in
# .glmm_fit_one): print an inline marker and raise a warning, so runaway
# coefficients or untrustworthy standard errors are never mistaken for valid
# estimates. Returns TRUE (invisibly) when the fit did not converge.
.glmm_warn_unconverged <- function(fit, label = NULL) {
  if (is.null(fit) || !isFALSE(attr(fit, "remstimate_converged")))
    return(invisible(FALSE))
  msg <- attr(fit, "remstimate_conv_msg") %||% "optimiser did not converge"
  tag <- if (!is.null(label)) paste0(label, ": ") else ""
  cat("\n** ", tag, "fit did NOT converge - estimates unreliable **\n", sep = "")
  warning(tag, "GLMM fit did not converge (", msg, "); fixed effects and ",
          "random-effect variances are unreliable and should not be interpreted.",
          call. = FALSE)
  invisible(TRUE)
}

#' @export
print.remstimate_glmm <- function(x, ...) {
  engine <- attr(x, "engine")
  cat("REM -", attr(x, "model"), "model - GLMM [", engine, "]\n\n")

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
      message("singular fit - consider simplifying the random structure")
    .glmm_warn_unconverged(x$backend_fit)
  } else {
    eng_r <- x$engine_receiver %||% "coxme"
    for (deel in c("sender_model", "receiver_model")) {
      if (is.null(x[[deel]])) next
      eng_lbl <- if (deel == "receiver_model")
        paste0(eng_r, ", conditional logit") else engine
      cat(deel, " [", eng_lbl, "] - fixed effects:\n", sep = "")
      print(x[[deel]]$coefficients); cat("\n")
      .glmm_warn_unconverged(x[[deel]]$backend_fit, deel)
    }
  }
  invisible(x)
}

#' @export
summary.remstimate_glmm <- function(object, ...) {
  if (attr(object, "model") == "tie") {
    .glmm_warn_unconverged(object$backend_fit)
    summary(object$backend_fit, ...)
  } else {
    .glmm_warn_unconverged(object$sender_model$backend_fit,   "sender_model")
    .glmm_warn_unconverged(object$receiver_model$backend_fit, "receiver_model")
    lapply(list(object$sender_model$backend_fit,
                object$receiver_model$backend_fit), summary)
  }
}

#' @export
coef.remstimate_glmm   <- function(object, ...) object$coefficients
#' @export
logLik.remstimate_glmm <- function(object, ...) object$loglik

#' @export
#' @method diagnostics remstimate_glmm
diagnostics.remstimate_glmm <- function(object, reh, stats,
                                        top_pct   = 0.05,
                                        use_ranef = TRUE, ...) {
  model  <- attr(object, "model")
  engine <- attr(object, "engine") %||% "lme4"

  # Duration models are stacked as a tie-like design but need the start/end
  # recall structure, not the Schoenfeld/waiting-time tie diagnostics. Route
  # them to the durem machinery with a random-effect-aware linear predictor.
  if (inherits(reh, "remify_durem") || inherits(stats, "remstats_durem")) {
    return(.glmm_durem_diagnostics(
      object, reh, stats, top_pct = top_pct,
      surprise_threshold = list(...)$surprise_threshold %||% 0.2,
      use_ranef = use_ranef, engine = engine))
  }

  # (a) Fixed-effects diagnostics. Schoenfeld residuals, the waiting-time rates
  #     and the fixed-effects recall all come from the shared MLE machinery by
  #     dropping the GLMM subclass and dispatching diagnostics.remstimate. The
  #     BLUPs deliberately do NOT enter here: Schoenfeld residuals and the
  #     waiting-time Q-Q test the fixed-effect structure, exactly as for MLE.
  obj_fixed <- object
  class(obj_fixed) <- "remstimate"
  diag_obj <- diagnostics(obj_fixed, reh, stats, top_pct = top_pct, ...)

  # (c) BLUPs, kept for the random-effect normality Q-Q (plot.diagnostics which = 7)
  diag_obj$ranef     <- .glmm_ranef(object, engine)
  diag_obj$.engine   <- engine
  diag_obj$use_ranef <- FALSE

  # (b) Random-effect-aware recall: the linear predictor includes the BLUPs, so
  #     each observed event is ranked within its risk set under the fitted mixed
  #     model. Stored alongside (not over) the fixed-effects recall so that the
  #     Schoenfeld residuals, surprises and offender tables in diag_obj stay
  #     internally consistent (they are fixed-effect quantities). plot.diagnostics
  #     prefers $recall_ranef for the recall/prob-ratio panels when present.
  if (isTRUE(use_ranef)) {
    rr <- switch(model,
                 tie = .glmm_ranef_recall_tie(object, top_pct),
                 NULL)  # actor: added in a later stage
    if (!is.null(rr)) {
      diag_obj$recall_ranef <- rr$recall
      # per-type recall for typed events (ext = TRUE); NULL for untyped models
      if (!is.null(rr$recall_by_type))
        diag_obj$recall_by_type <- rr$recall_by_type
      diag_obj$use_ranef <- TRUE
    }
  }
  diag_obj
}

# BLUPs per grouping factor, engine-agnostic (lme4 / glmmTMB / coxme). Returns a
# named list keyed by sub-model ('fit' for tie; 'sender_model'/'receiver_model'
# for actor), each element being the engine's native ranef structure.
.glmm_ranef <- function(object, engine) {
  fits <- if (attr(object, "model") == "tie")
    list(fit = object$backend_fit)
  else
    list(sender_model   = object$sender_model$backend_fit,
         receiver_model = object$receiver_model$backend_fit)

  out <- list()
  for (nm in names(fits)) {
    f <- fits[[nm]]
    if (is.null(f)) next
    re <- tryCatch(
      if (inherits(f, "coxme"))        coxme::ranef(f)
      else if (inherits(f, "glmmTMB")) glmmTMB::ranef(f)$cond
      else                             lme4::ranef(f),
      error = function(e) NULL)
    if (!is.null(re)) out[[nm]] <- re
  }
  if (length(out)) out else NULL
}

# Random-effect-aware recall for the tie model. predict() returns the linear
# predictor including the BLUPs; the timing offset is constant within an event's
# risk set and therefore cancels in the within-event softmax used by
# .recall_block, so it does not need to be stripped.
.glmm_ranef_recall_tie <- function(object, top_pct) {
  df  <- object$stacked_data
  fit <- object$backend_fit
  if (is.null(df) || is.null(fit)) return(NULL)
  lp <- tryCatch(stats::predict(fit, type = "link"),
                 error = function(e)
                   tryCatch(stats::predict(fit, type = "lp"),
                            error = function(e) NULL))
  if (is.null(lp)) return(NULL)
  # .diagnostics_recall returns $recall plus $recall_by_type when df carries a
  # 'type' column with >1 level (typed events, ext = TRUE).
  dr <- .diagnostics_recall(lp, df, top_pct)
  if (is.null(dr$recall)) return(NULL)
  # plot.diagnostics keys the recall x-axis on 'event'; add a dense index over
  # the time groups so each recall table plots like every other one.
  .add_event <- function(rc) {
    if (!is.null(rc) && !is.null(rc$per_event))
      rc$per_event$event <- match(rc$per_event$time, sort(unique(rc$per_event$time)))
    rc
  }
  dr$recall <- .add_event(dr$recall)
  if (!is.null(dr$recall_by_type))
    dr$recall_by_type <- lapply(dr$recall_by_type, .add_event)
  dr
}

# Duration-model GLMM diagnostics. Builds the durem start/end recall tables
# (via the shared .durem_recall_tables) from a random-effect-aware linear
# predictor: predict() on the coxme / glmmTMB fit includes the BLUPs; the timing
# offset is constant within an event's risk set and cancels in the recall
# softmax. Falls back to the fixed-effect linear predictor if predict() fails.
# Returns a 'diagnostics_durem' object so plot()/print() dispatch as for MLE.
.glmm_durem_diagnostics <- function(object, reh, stats, top_pct,
                                    surprise_threshold, use_ranef, engine) {
  df  <- object$stacked_data
  fit <- object$backend_fit
  if (is.null(df) || is.null(fit))
    stop("No stacked data / backend fit stored in the durem GLMM object.", call. = FALSE)

  sn_start <- stats$stacked$stat_names_start
  sn_end   <- stats$stacked$stat_names_end

  lp <- NULL
  if (isTRUE(use_ranef))
    lp <- tryCatch(stats::predict(fit, type = "link"),
                   error = function(e)
                     tryCatch(stats::predict(fit, type = "lp"),
                              error = function(e) NULL))
  ranef_used <- !is.null(lp)
  if (is.null(lp)) {                       # fixed-effects fallback
    sn    <- attr(object, "statistics")
    coefs <- object$coefficients
    X     <- as.matrix(df[, sn, drop = FALSE])
    lp    <- as.numeric(X %*% coefs[sn])
  }

  out <- .durem_recall_tables(lp, df, reh, sn_start, sn_end,
                              top_pct = top_pct,
                              surprise_threshold = surprise_threshold)

  dev <- tryCatch(stats::residuals(fit, type = "deviance"), error = function(e) NULL)
  if (!is.null(dev)) out$deviance_residuals <- dev

  out$ranef          <- .glmm_ranef(object, engine)
  out$.engine        <- engine
  out$use_ranef      <- ranef_used
  out$.reh.processed <- reh
  class(out) <- c("diagnostics_durem", "diagnostics", "remstimate")
  out
}
