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

#' Dyadic Latent Class REM (Lakdawala et al.)
#'
#' Fits a mixture REM where dyads are assigned to K latent classes,
#' each with class-specific regression coefficients.
#'
#' @param reh A \code{remify} object.
#' @param stats A \code{remstats} object.
#' @param k Number of latent classes (default 2).
#' @param nrep Number of random restarts (default 3).
#' @param ... Additional arguments passed to \code{remstimate}.
#' @return A \code{remstimate_mixrem} object.
#' @references Lakdawala, Leenders, and Mulder (2025). Not all bonds are created
#' equal: Dyadic latent class models for relational event data.
#' @export
dlcrem <- function(reh, stats, k = 2L, nrep = 3L, ...) {
  sn <- .get_durem_or_standard_stat_names(stats)
  random <- stats::as.formula(
    paste0("~ (1 + ", paste(sn, collapse = " + "), " | dyad)")
  )
  remstimate(reh, stats, method = "MIXREM", random = random,
             k = k, nrep = nrep, ...)
}


frailty_rem <- function(reh, stats, engine = "lme4", ...) {
  is_durem    <- inherits(reh, "remify_durem")
  is_actor    <- inherits(stats, "aomstats")
  ext_by_type <- isTRUE(reh$meta$with_type_riskset)

  if (is_actor) {
    random <- list(
      sender   = ~ (1 | actor1),
      receiver = ~ (1 | actor2)
    )
  } else if (is_durem && ext_by_type) {
    random <- ~ (1 | actor1:process:type) + (1 | actor2:process:type)
  } else if (is_durem) {
    random <- ~ (1 | actor1:process) + (1 | actor2:process)
  } else if (ext_by_type) {
    random <- ~ (1 | actor1:type) + (1 | actor2:type)
  } else {
    random <- ~ (1 | actor1) + (1 | actor2)
  }

  remstimate(reh, stats, method = "GLMM", random = random,
             engine = engine, ...)
}

