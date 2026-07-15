# ── strata() shim ────────────────────────────────────────────────────────────
#
# 'survival' is only a Suggested dependency, so it is never attached to the
# search path (and cannot be @importFrom'd). When survival::clogit(),
# survival::coxph() or coxme::coxme() build their model frames, the strata()
# term in a formula is resolved as a *function* using the formula's environment
# and its parents. For formulas constructed inside this package that environment
# is the package namespace, where 'strata' is not defined -> the lookup falls
# through to the (unattached) search path and fails with:
#
#   Error in strata(time_index) : could not find function "strata"
#
# Whether it fails is platform-dependent: it only works when some *other*
# package (coxme, lme4, glmmTMB, ...) happens to have attached 'survival'
# first. That is why R-hub Windows passed but macOS (and clean CRAN checks)
# fail.
#
# Defining an internal 'strata' in the package namespace makes a bare
# strata(...) in any package-built formula resolve here, on every platform,
# regardless of load order. It simply forwards to survival::strata(), so the
# object returned is identical and each backend's "specials" detection - which
# keys on the symbol name "strata", not on what it evaluates to - keeps working.
#
#' @keywords internal
#' @noRd
strata <- function(...) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("Package 'survival' is required for stratified (clogit/coxph/coxme) ",
         "duration models. Install it with install.packages('survival').",
         call. = FALSE)
  survival::strata(...)
}
