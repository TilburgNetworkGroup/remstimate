# GLMNET backend — ridge / lasso / elastic net via glmnet::cv.glmnet
#
# remstimate(reh, stats, method = "GLMNET", alpha = 1, lambda_select = "1se")

.remstimate_glmnet <- function(reh, stats,
                                alpha         = 1,
                                nfolds        = 10L,
                                lambda_select = c("1se", "min"),
                                ...) {

  lambda_select <- match.arg(lambda_select)
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("install.packages('glmnet')")

  s <- .remstimate_make_stack(reh, stats, add_actors = FALSE)

  if (s$model == "tie") {

    res <- .glmnet_fit_one(df = s$df, stat_names = s$stat_names, ordinal = s$ordinal,
                           alpha = alpha, nfolds = nfolds, lambda_select = lambda_select, ...)

    .remstimate_wrap(
      coefficients = res$coefs,
      stat_names   = s$stat_names,
      stacked_data = s$df,
      backend_fit  = res$cv_fit,
      model        = "tie",
      method       = "GLMNET",
      engine       = "glmnet",
      ordinal      = s$ordinal,
      extra = list(
        lambda_min = res$lambda_min,
        lambda_1se = res$lambda_1se,
        lambda_sel = res$lambda_sel,
        alpha      = alpha
      )
    )

  } else {

    res_s <- res_r <- NULL
    if (!is.null(s$df$sender))
      res_s <- .glmnet_fit_one(s$df$sender, s$stat_names$sender_model,
                                s$ordinal, alpha, nfolds, lambda_select, ...)
    if (!is.null(s$df$receiver))
      res_r <- .glmnet_fit_one(s$df$receiver, s$stat_names$receiver_model,
                                s$ordinal, alpha, nfolds, lambda_select, ...)

    .remstimate_wrap(
      coefficients = list(sender_model   = if (!is.null(res_s)) res_s$coefs else NULL,
                          receiver_model = if (!is.null(res_r)) res_r$coefs else NULL),
      stat_names   = s$stat_names,
      stacked_data = s$df,
      backend_fit  = list(sender_model   = if (!is.null(res_s)) res_s$cv_fit else NULL,
                          receiver_model = if (!is.null(res_r)) res_r$cv_fit else NULL),
      model   = "actor",
      method  = "GLMNET",
      engine  = "glmnet",
      ordinal = s$ordinal,
      extra   = list(alpha = alpha, lambda_select = lambda_select)
    )
  }
}

.glmnet_fit_one <- function(df, stat_names, ordinal, alpha, nfolds, lambda_select, ...) {
  mm     <- .remstimate_model_matrix(df, stat_names, ordinal)
  cv_fit <- glmnet::cv.glmnet(mm$X, mm$y,
                               family  = if (!ordinal) "poisson" else "binomial",
                               alpha   = alpha,
                               nfolds  = nfolds,
                               offset  = mm$offset,
                               weights = mm$weights,
                               ...)

  lambda_sel <- if (lambda_select == "min") cv_fit$lambda.min else cv_fit$lambda.1se
  ruw        <- glmnet::coef.glmnet(cv_fit, s = lambda_sel)
  coefs      <- as.vector(ruw)
  names(coefs) <- rownames(ruw)
  coefs      <- coefs[names(coefs) != "(Intercept)"]  # baseline dient als intercept

  list(coefs      = coefs,
       cv_fit     = cv_fit,
       lambda_min = cv_fit$lambda.min,
       lambda_1se = cv_fit$lambda.1se,
       lambda_sel = lambda_sel)
}

#' @export
print.remstimate_glmnet <- function(x, ...) {
  methode <- if (x$alpha == 1) "lasso" else if (x$alpha == 0) "ridge" else
    paste0("elastic net (alpha=", x$alpha, ")")
  cat("REM —", attr(x, "model"), "model — GLMNET [", methode, "]\n\n")
  if (attr(x, "model") == "tie") {
    cat("Coefficients (lambda.1se):\n")
    print(x$coefficients)
    cat("\nlambda.min =", round(x$lambda_min, 6),
        " | lambda.1se =", round(x$lambda_1se, 6), "\n")
  } else {
    for (deel in c("sender_model", "receiver_model")) {
      if (is.null(x[[deel]])) next
      cat(deel, ":\n"); print(x[[deel]]$coefficients); cat("\n")
    }
  }
  invisible(x)
}

#' @export
summary.remstimate_glmnet <- function(object, ...) {
  if (attr(object, "model") == "tie") print(object$backend_fit)
  else lapply(object$backend_fit, print)
}

#' @export
coef.remstimate_glmnet <- function(object, ...) object$coefficients

#' @export
#' @method diagnostics remstimate_glmnet
diagnostics.remstimate_glmnet <- function(object, reh = NULL, stats = NULL,
                                           top_pct = 0.05, ...) {
  df <- object$stacked_data
  if (is.null(df)) stop("No stacked data stored in fit object.", call. = FALSE)

  stat_names <- attr(object, "statistics")
  coefs      <- object$coefficients

  X  <- as.matrix(df[, stat_names, drop = FALSE])
  lp <- as.numeric(X %*% coefs[stat_names])

  out <- .diagnostics_recall(lp, df, top_pct)
  class(out) <- c("diagnostics_remstimate", "diagnostics", "remstimate")
  out
}

#' @export
#' @method print diagnostics_remstimate
print.diagnostics_remstimate <- function(x, ...) {
  cat("Diagnostics for Relational Event Model\n")
  cat(strrep("-", 50), "\n")

  if (!is.null(x$recall)) {
    cat("\nRecall:\n")
    .print_recall_summary(x$recall, "Joint")
  }
  if (!is.null(x$recall_by_type)) {
    cat("\nRecall by type:\n")
    for (tp in names(x$recall_by_type))
      .print_recall_summary(x$recall_by_type[[tp]], paste0("  ", tp))
  }
  invisible(x)
}

#' @export
#' @method plot remstimate_glmnet
plot.remstimate_glmnet <- function(x, reh = NULL, stats = NULL,
                                    diagnostics_object = NULL,
                                    which = 1L, ...) {
  if (is.null(diagnostics_object))
    diagnostics_object <- diagnostics(x, reh, stats)

  if (which == 1L) {
    .plot_recall_scatter(diagnostics_object$recall, "Recall: GLMNET", ...)
  } else if (which == 2L && !is.null(diagnostics_object$recall_by_type)) {
    rbt <- diagnostics_object$recall_by_type
    old_par <- graphics::par(mfrow = c(1, length(rbt)))
    on.exit(graphics::par(old_par))
    for (tp in names(rbt))
      .plot_recall_scatter(rbt[[tp]], paste("Type:", tp), ...)
  }
  invisible(x)
}
