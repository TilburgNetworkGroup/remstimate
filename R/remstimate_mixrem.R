# MIXREM backend — eindige mengsels via flexmix::flexmix
#
# remstimate(reh, stats, method = "MIXREM",
#            random = ~ (1 + inertia | dyad), k = 2)
#
# random specificeert welke effecten per component variëren en de groeperingsvariabele
# k als vector geeft meerdere fits; gebruik bic_table() om te vergelijken

.remstimate_mixrem <- function(reh, stats,
                                random      = NULL,
                                k           = 2L,
                                concomitant = NULL,
                                nrep        = 3L,
                                ...) {

  if (!requireNamespace("flexmix", quietly = TRUE))
    stop("install.packages('flexmix')")
  if (is.null(random))
    stop("specify clustering, e.g. random = ~ (1 + inertia | dyad)")

  s       <- .remstimate_make_stack(reh, stats, add_actors = TRUE)
  geparsed <- .parse_mixrem_formula(random)

  if (length(k) == 1L) {
    if (s$model == "tie")
      return(.mixrem_fit_one(s$df, s$stat_names, s$ordinal,
                              geparsed, k, concomitant, nrep, model = "tie", ...))

    fit_s <- fit_r <- NULL
    if (!is.null(s$df$sender))
      fit_s <- .mixrem_fit_one(s$df$sender, s$stat_names$sender_model, s$ordinal,
                                geparsed, k, concomitant, nrep, model = "actor_sender", ...)
    if (!is.null(s$df$receiver))
      fit_r <- .mixrem_fit_one(s$df$receiver, s$stat_names$receiver_model, s$ordinal,
                                geparsed, k, concomitant, nrep, model = "actor_receiver", ...)
    return(list(sender_model = fit_s, receiver_model = fit_r))
  }

  # meerdere k-waarden
  fits <- lapply(k, function(ki) {
    if (s$model == "tie")
      .mixrem_fit_one(s$df, s$stat_names, s$ordinal,
                      geparsed, ki, concomitant, nrep, model = "tie", ...)
    else
      list(
        sender_model   = if (!is.null(s$df$sender))
          .mixrem_fit_one(s$df$sender, s$stat_names$sender_model, s$ordinal,
                          geparsed, ki, concomitant, nrep, model = "actor_sender", ...) else NULL,
        receiver_model = if (!is.null(s$df$receiver))
          .mixrem_fit_one(s$df$receiver, s$stat_names$receiver_model, s$ordinal,
                          geparsed, ki, concomitant, nrep, model = "actor_receiver", ...) else NULL
      )
  })
  names(fits) <- paste0("k", k)
  structure(fits, class = c("remstimate_mixrem_list", "list"), k = k)
}

.mixrem_fit_one <- function(df, stat_names, ordinal, geparsed, k,
                              concomitant, nrep, model = "tie", ...) {

  vaste_kant <- .remstimate_fixed_rhs(stat_names, ordinal)
  if (ordinal) {
    df$event   <- factor(df$event)
    vaste_kant <- sub("^-1 \\+ ", "-1 + event + ", vaste_kant)
    familie    <- flexmix::FLXMRglm(family = "binomial")
  } else {
    familie <- flexmix::FLXMRglm(family = "poisson")
  }

  fml  <- as.formula(paste("obs ~", vaste_kant, "|", geparsed$group))
  conc <- if (!is.null(concomitant)) flexmix::FLXPmultinom(concomitant) else flexmix::FLXPconstant()

  fit <- flexmix::flexmix(fml, data = df, k = k,
                           model       = familie,
                           concomitant = conc,
                           control     = list(nrep = nrep), ...)

  coef_mat <- flexmix::parameters(fit)
  if (is.null(dim(coef_mat)))
    coef_mat <- matrix(coef_mat, ncol = 1, dimnames = list(names(coef_mat), "Comp.1"))

  kansen <- flexmix::prior(fit)
  volgorde <- order(kansen, decreasing = TRUE)
  coef_mat <- coef_mat[, volgorde, drop = FALSE]
  kansen   <- kansen[volgorde]
  colnames(coef_mat) <- paste0("Component.", seq_len(k))

  .remstimate_wrap(
    coefficients = coef_mat,
    stat_names   = stat_names,
    loglik       = flexmix::logLik(fit),
    formula      = fml,
    stacked_data = df,
    backend_fit  = fit,
    model        = sub("_.*", "", model),
    method       = "MIXREM",
    engine       = "flexmix",
    ordinal      = ordinal,
    extra        = list(k = k, prior_probs = kansen, bic = flexmix::BIC(fit))
  )
}

.parse_mixrem_formula <- function(random) {
  rhs   <- deparse(random[[2]])
  groep_match <- regmatches(rhs,
    regexpr("\\|\\s*([A-Za-z_.][A-Za-z_.0-9]*)", rhs, perl = TRUE))
  if (!length(groep_match))
    stop("random must contain a grouping variable after |, e.g. ~ (1 + inertia | dyad)")

  groep  <- trimws(sub("^\\|\\s*", "", groep_match))
  binnen <- sub("\\|.*", "", sub("^\\(", "", rhs))
  termen <- trimws(strsplit(binnen, "\\+")[[1]])
  termen <- termen[termen != "1" & nchar(termen) > 0]

  list(effects = if (length(termen)) termen else NULL, group = groep)
}

#' @export
bic_table <- function(x, ...) UseMethod("bic_table")

#' @export
bic_table.remstimate_mixrem_list <- function(x, ...) {
  k_vals <- attr(x, "k")
  bics   <- sapply(x, function(fit)
    if (inherits(fit, "remstimate_mixrem")) fit$bic else NA_real_)
  out <- data.frame(k = k_vals, BIC = bics)
  out$delta_BIC <- out$BIC - min(out$BIC, na.rm = TRUE)
  out[order(out$BIC), ]
}

#' @export
print.remstimate_mixrem <- function(x, ...) {
  cat("REM —", attr(x, "model"), "model — MIXREM [flexmix, k =", x$k, "]\n\n")
  cat("Mixing proportions:\n"); print(round(x$prior_probs, 4))
  cat("\nCoefficients per component:\n"); print(round(x$coefficients, 4))
  cat("\nlogLik:", round(as.numeric(x$loglik), 4), " BIC:", round(x$bic, 2), "\n")
  invisible(x)
}

#' @export
print.remstimate_mixrem_list <- function(x, ...) {
  cat("MIXREM fits voor k =", paste(attr(x, "k"), collapse = ", "), "\n")
  print(bic_table(x))
  invisible(x)
}

#' @export
summary.remstimate_mixrem  <- function(object, ...) summary(object$backend_fit, ...)
#' @export
coef.remstimate_mixrem     <- function(object, ...) object$coefficients
#' @export
logLik.remstimate_mixrem   <- function(object, ...) object$loglik
