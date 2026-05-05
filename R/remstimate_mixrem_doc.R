#' Finite mixture REM via flexmix
#'
#' Extends \code{remstimate()} with \code{method = "MIXREM"} to fit a finite
#' mixture of Poisson REM components.  Each component gets its own coefficient
#' vector; the clustering unit (dyad, actor, …) is specified via the
#' \code{random} argument using flexmix's \code{|} syntax.
#'
#' When \code{k} is a vector, all values are fitted and a
#' \code{remstimate_mixrem_list} is returned.  Use \code{bic_table()} to
#' compare fits and select \code{k}.
#'
#' An optional \code{concomitant} formula adds a multinomial model for
#' component membership probabilities as a function of actor- or dyad-level
#' covariates.
#'
#' @param reh         A \code{remify} object.
#' @param stats       A \code{tomstats} or \code{aomstats} object.
#' @param method      Must be \code{"MIXREM"}.
#' @param random      One-sided formula specifying which effects vary across
#'   components and the clustering unit, e.g.
#'   \code{~ (1 + inertia | dyad)}.  The grouping variable after \code{|}
#'   must be a column in the stacked data (\code{dyad}, \code{actor1},
#'   \code{actor2}).
#' @param k           Integer or integer vector of component counts to fit.
#'   If a vector, returns a \code{remstimate_mixrem_list}.
#' @param concomitant Optional one-sided formula for the concomitant model,
#'   e.g. \code{~ actor1_degree}.
#' @param nrep        Number of random restarts passed to \code{flexmix}
#'   (default \code{3}).
#' @param ...         Passed to \code{flexmix::flexmix()}.
#'
#' @return For a single \code{k}: a \code{remstimate_mixrem} object with:
#'   \describe{
#'     \item{\code{coefficients}}{Matrix \eqn{[p \times k]} of per-component
#'       coefficients, sorted by decreasing mixing proportion.}
#'     \item{\code{prior_probs}}{Mixing proportions (sorted).}
#'     \item{\code{bic}}{BIC of the fitted model.}
#'     \item{\code{backend_fit}}{The raw \code{flexmix} object.}
#'   }
#'   For a vector \code{k}: a \code{remstimate_mixrem_list}; use
#'   \code{bic_table()} to compare.
#'
#' @seealso \code{\link{bic_table}}, \code{\link{remstimate}},
#'   \code{\link[flexmix]{flexmix}}
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
#'               outdegreeSender(consider_type = FALSE)
#' stats <- remstats::tomstats(effects, reh = reh, attr_actors = info,
#'                              memory = "decay", memory_value = 1000)
#'
#' # two-component mixture, dyads as clustering unit
#' fit2 <- remstimate(reh, stats, method = "MIXREM",
#'                    random = ~ (1 + inertia | dyad), k = 2)
#' fit2
#'
#' # coefficients: rows = effects, columns = components
#' fit2$coefficients
#' fit2$prior_probs
#'
#' # which dyads ended up in which component?
#' flexmix::clusters(fit2$backend_fit)
#'
#' # try k = 2:4 and compare via BIC
#' fits <- remstimate(reh, stats, method = "MIXREM",
#'                    random = ~ (1 + inertia | dyad), k = 2:4)
#' bic_table(fits)
#' plot(fits)   # prints BIC table
#'
#' # pick the best k and inspect
#' best <- fits[[ paste0("k", bic_table(fits)$k[1]) ]]
#' best
#'
#' # diagnostics: per-component recall + combined
#' diag2 <- diagnostics(fit2, reh, stats)
#' diag2
#' plot(fit2, reh, stats = stats)
#'
#' # concomitant model: component membership depends on sender out-degree
#' fit_conc <- remstimate(reh, stats, method = "MIXREM",
#'                         random      = ~ (1 + inertia | dyad),
#'                         k           = 2,
#'                         concomitant = ~ outdegreeSender)
#' }
#'
#' @name remstimate_mixrem
NULL


#' Compare BIC across a list of MIXREM fits
#'
#' @param x   A \code{remstimate_mixrem_list} returned when \code{k} is a
#'   vector in \code{remstimate(..., method = "MIXREM")}.
#' @param ... Unused.
#'
#' @return A data frame with columns \code{k}, \code{BIC}, and
#'   \code{delta_BIC}, sorted by ascending BIC.
#'
#' @examples
#' \donttest{
#' library(remstimate)
#' data(history, package = "remstats")
#' data(info,    package = "remstats")
#' colnames(history)[colnames(history) == "setting"] <- "type"
#'
#' reh   <- remify::remify(edgelist = history[1:100, ], model = "tie",
#'                          riskset = "active")
#' stats <- remstats::tomstats(~ inertia(consider_type = FALSE), reh = reh,
#'                              attr_actors = info)
#'
#' fits <- remstimate(reh, stats, method = "MIXREM",
#'                    random = ~ (1 + inertia | dyad), k = 2:5)
#' bic_table(fits)
#' }
#'
#' @name bic_table
#' @export
bic_table
