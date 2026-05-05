#' Regularised REM via glmnet (lasso / ridge / elastic net)
#'
#' Extends \code{remstimate()} with \code{method = "GLMNET"} to fit a
#' penalised Poisson model on the stacked statistics.  Cross-validation via
#' \code{cv.glmnet} is used to select \eqn{\lambda}; both
#' \code{lambda.min} and \code{lambda.1se} are stored so the user can
#' reselect after fitting.
#'
#' @param reh           A \code{remify} object.
#' @param stats         A \code{tomstats} or \code{aomstats} object.
#' @param method        Must be \code{"GLMNET"}.
#' @param alpha         Elastic-net mixing parameter: \code{1} = lasso
#'   (default), \code{0} = ridge, values in between give elastic net.
#' @param nfolds        Number of cross-validation folds. Default \code{10}.
#' @param lambda_select Which lambda to use for the returned coefficients:
#'   \code{"1se"} (default, more regularised) or \code{"min"} (lowest CV
#'   error).  The full path is always stored in \code{$backend_fit}.
#' @param ...           Passed to \code{cv.glmnet()}.
#'
#' @return A \code{remstimate_glmnet} object with slots:
#'   \describe{
#'     \item{\code{coefficients}}{Coefficients at the selected \eqn{\lambda}.}
#'     \item{\code{backend_fit}}{The \code{cv.glmnet} object; call
#'       \code{plot(fit$backend_fit)} for the full regularisation path.}
#'     \item{\code{lambda_min}, \code{lambda_1se}}{Cross-validated lambdas.}
#'     \item{\code{alpha}}{The mixing parameter used.}
#'   }
#'
#' @seealso \code{\link{remstimate}}, \code{\link[glmnet]{cv.glmnet}}
#'
#' @examples
#' \donttest{
#' library(remstimate)
#' data(history, package = "remstats")
#' data(info,    package = "remstats")
#' colnames(history)[colnames(history) == "setting"] <- "type"
#'
#' reh <- remify::remify(edgelist = history[1:100, ], model = "tie",
#'                        riskset = "active")
#' effects <- ~ inertia(consider_type = FALSE) +
#'               indegreeSender(consider_type = FALSE) +
#'               outdegreeSender(consider_type = FALSE) +
#'               reciprocity(consider_type = FALSE)
#' stats <- remstats::tomstats(effects, reh = reh, attr_actors = info,
#'                              memory = "decay", memory_value = 1000)
#'
#' # lasso (alpha = 1): shrinks small effects to zero
#' fit_lasso <- remstimate(reh, stats, method = "GLMNET", alpha = 1)
#' fit_lasso
#'
#' # ridge (alpha = 0): shrinks all effects towards zero but keeps them in
#' fit_ridge <- remstimate(reh, stats, method = "GLMNET", alpha = 0)
#'
#' # elastic net
#' fit_en <- remstimate(reh, stats, method = "GLMNET", alpha = 0.5)
#'
#' # use lambda.min instead of lambda.1se
#' fit_min <- remstimate(reh, stats, method = "GLMNET",
#'                        alpha = 1, lambda_select = "min")
#'
#' # regularisation path (plot 6)
#' plot(fit_lasso, reh, stats = stats, which = 6)
#'
#' # re-extract coefficients at a different lambda after fitting
#' glmnet::coef.glmnet(fit_lasso$backend_fit, s = fit_lasso$lambda_min)
#'
#' # diagnostics work the same as for MLE
#' diag_lasso <- diagnostics(fit_lasso, reh, stats)
#' plot(diag_lasso)
#' }
#'
#' @name remstimate_glmnet
NULL
