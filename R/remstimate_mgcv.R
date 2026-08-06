# GAM backend - smooth effects via mgcv
#
# remstimate(reh, stats, gam = list(smooths = "inertia", bs = "tp", k = 10))
# remstimate(reh, stats, gam = list(tve = "inertia", baseline = TRUE))

# Smoothing bases accepted by mgcv::s(). Kept as a constant so the default
# argument, the validation and the per-smooth matching all agree.
# Restricted to the bases that can actually be built with univariate statistics
# ("sc" is the shape-constrained B-spline basis, see .gam_constraint_bs()).
.gam_bs_choices <- c("tp", "ts", "ds", "cr", "cs", "cc", "bs", "ps", "cp",
                     "re", "gp", "sc")
.gam_bs_not_smooth <- "re"

# Shape constraints implemented by mgcv's "sc" basis
.gam_constraint_codes <- c("m+", "m-", "c+", "c-", "+")

# Constraints that conflict with each other and cannot be applied to the same smooth.
.gam_constraint_conflicts <- list(c("m+", "m-"), c("c+", "c-"))

# Bases a time-varying effect can be built on. 're' would put one coefficient on
# every distinct event time, and 'sc' should not be implementable on a time-varying effect because it is a B-spline basis with no anchor point to constrain, so both are excluded.
.gam_tve_bs_choices <- setdiff(.gam_bs_choices, c("re", "sc"))

# Column the time axis is written to. Exposed in the mgcv labels ('s(t)',
# 's(t):inertia')
.gam_tve_time_col <- "t"


.remstimate_gam <- function(reh, stats,
                            smooths = NULL, # vector of variable names to be smoothed (s() in mgcv)
                            tve = NULL,     # vector of variable names with a time-varying coefficient
                            tve_baseline = TRUE,
                            bs = .gam_bs_choices,
                            k = NULL,
                            pc = 0,
                            pc_supplied = FALSE,
                            constraints = NULL,
                            select = FALSE,
                            method = NULL,
                            engine = c("auto", "gam", "bam", "scasm"),
                            ...) {

  if (!requireNamespace("mgcv", quietly = TRUE))
    stop("install.packages('mgcv')")

  if (is.null(bs) || !length(bs)) bs <- "tp"
  if (identical(bs, .gam_bs_choices)) bs <- "tp"
  bs <- vapply(bs, function(b) match.arg(b, .gam_bs_choices),
               character(1), USE.NAMES = FALSE)

  engine <- match.arg(engine)

  if (inherits(reh, "remify_durem") || inherits(stats, "remstats_durem"))
    stop("The GAM backend supports standard (tie-oriented) relational event ",
         "models only, not duration models.", call. = FALSE)

  s <- .remstimate_make_stack(reh, stats, add_actors = FALSE)

  # A duration design handed over PRE-STACKED ('remstats_stacked') carries no
  # 'remstats_durem' class, so the check above lets it through; its start/end
  # layout is identified by the 'process' column that stack_stats() adds.
  if (is.data.frame(s$df) && !is.null(s$df$process))
    stop("The GAM backend supports standard (tie-oriented) relational event ",
         "models only; this stacked design is a duration model (it carries a ",
         "'process' column).", call. = FALSE)

  .check_smooths_names(smooths,
    valid = if (s$model == "tie") s$stat_names
            else union(s$stat_names$sender_model, s$stat_names$receiver_model))

  # check if request is valid
  .gam_tve_check(smooths, tve, tve_baseline, bs, constraints, pc_supplied, s)

  # 'bs' / 'k' follow the same positional contract in both modes; 'pc' belongs
  # to the nonlinear mode only.
  if (length(tve)) .remstimate_gam_spec_check(tve, bs, k, NULL, arg = "tve")
  else             .remstimate_gam_spec_check(smooths, bs, k, pc)

  con    <- .gam_parse_constraints(constraints, smooths)
  engine <- .gam_constraint_engine(engine, con, method)

  if (s$model == "tie") {

    res <- .mgcv_fit(s, reh = reh,
                     smooths = smooths, tve = tve, tve_baseline = tve_baseline,
                     bs = bs, k = k, pc = pc, con = con,
                     select = select, method = method, engine = engine, ...)

    return(res)

  } else {
    stop("GAM backend is not yet implemented for actor models (sender/receiver).",
         call. = FALSE)
  }
}

# 'bs', 'k' and 'pc' must be length 1 (one value shared by every smooth) or
# length(smooths) (one value per smooth, in order). Any other length is a user
# error that must be caught here rather than left to the rep_len() broadcast in
# .gam_smooth_terms(): a length that is neither 1 nor length(smooths) means the
# caller's per-smooth intent cannot be honoured, and silently recycling it would
# attach 'k' or 'bs' to the wrong statistic.
#
# A repeated name in 'smooths' is the same failure by another route: it emits
# two s() terms for one variable, which mgcv only warns about, so it is an error
# here too. It cannot be silently de-duplicated either, because that would
# misalign a per-smooth 'bs' / 'k' / 'pc' with the smooths that remain.
.remstimate_gam_spec_check <- function(smooths, bs, k, pc = NULL,
                                       arg = "smooths") {
  n <- length(smooths)
  if (n == 0L) return(invisible(NULL))
  dup <- unique(smooths[duplicated(smooths)])
  if (length(dup))
    stop("'", arg, "' names the same statistic more than once: ",
         paste(dup, collapse = ", "),
         ".\n  One s() term per statistic; give a single entry per name.",
         call. = FALSE)
  for (nm in c("bs", "k", "pc")) {
    v <- switch(nm, bs = bs, k = k, pc = pc)
    if (is.null(v)) next
    if (!length(v) %in% c(1L, n))
      stop("'", nm, "' must be length 1 (shared by all smooths) or length ", n,
           " (one per smooth, matching '", arg, "'); got length ", length(v),
           ".", call. = FALSE)
  }
  invisible(NULL)
}

# Normalise 'constraints' - positional (as 'bs' / 'k' / 'pc') or named after the
# smooths - into one character vector of mgcv codes per smooth.
.gam_parse_constraints <- function(constraints, smooths) {
  n <- length(smooths)
  if (is.null(constraints) || !length(constraints))
    return(stats::setNames(rep(list(character(0)), n), smooths))
  if (n == 0L)
    stop("'constraints' was given but 'smooths' is empty: a shape constraint ",
         "applies to a smooth term, so there is nothing to constrain.",
         call. = FALSE)

  # collapse the list form, list(inertia = c("m+", "c-")), onto the vector form
  if (is.list(constraints))
    constraints <- vapply(constraints, function(z) paste(z, collapse = " "),
                          character(1))
  if (!is.character(constraints))
    stop("'constraints' must be a character vector, or a list of character ",
         "vectors, of mgcv shape-constraint codes (",
         paste0('"', .gam_constraint_codes, '"', collapse = ", "), ").",
         call. = FALSE)

  nms <- names(constraints)
  if (!is.null(nms) && any(nzchar(nms))) {
    if (!all(nzchar(nms)))
      stop("'constraints' is only partly named. Name every element (each name ",
           "one of 'smooths'), or none of them and give one entry per smooth ",
           "in order.", call. = FALSE)
    bad <- setdiff(nms, smooths)
    if (length(bad))
      stop("'constraints' names statistic(s) that are not among 'smooths': ",
           paste(bad, collapse = ", "),
           ".\n  A shape constraint applies to a smooth term; smooths are: ",
           paste(smooths, collapse = ", "), ".", call. = FALSE)
    if (anyDuplicated(nms))
      stop("'constraints' names the same smooth more than once: ",
           paste(unique(nms[duplicated(nms)]), collapse = ", "), ".",
           call. = FALSE)
    v <- rep(NA_character_, n)
    v[match(nms, smooths)] <- constraints
  } else {
    if (!length(constraints) %in% c(1L, n))
      stop("'constraints' must be length 1 (shared by all smooths), length ", n,
           " (one per smooth, matching 'smooths'), or named after the ",
           "smooths it applies to; got length ", length(constraints), ".",
           call. = FALSE)
    v <- rep_len(constraints, n)
  }

  stats::setNames(lapply(seq_len(n), function(i) .gam_codes(v[i], smooths[i])),
                  smooths)
}

# One entry of 'constraints' -> the mgcv codes it asks for.
.gam_codes <- function(x, nm) {
  if (is.na(x)) return(character(0))
  parts <- strsplit(x, "[,[:space:]]+")[[1L]]
  parts <- parts[nzchar(parts)]
  parts <- parts[!parts %in% c("none", "NA")]
  if (!length(parts)) return(character(0))
  bad <- setdiff(parts, .gam_constraint_codes)
  if (length(bad))
    stop("unknown shape constraint for '", nm, "': ",
         paste0('"', bad, '"', collapse = ", "),
         ".\n  mgcv's codes are ",
         paste0('"', .gam_constraint_codes, '"', collapse = ", "),
         ' (monotone increasing / decreasing, convex, concave, positive).',
         '\n  Combine them for one smooth with a space, e.g. "m+ c-"; ',
         '"none" leaves the smooth unconstrained.', call. = FALSE)
  codes <- unique(parts)
  for (p in .gam_constraint_conflicts)
    if (all(p %in% codes))
      stop("contradictory shape constraints for '", nm, "': ",
           paste0('"', p, '"', collapse = " and "),
           " cannot both hold. mgcv would silently apply only one of them.",
           call. = FALSE)
  codes
}

# Only scasm() imposes shape constraints; gam()/bam() build the 'sc' basis but
# ignore its constraint matrix, so they would fit silently unconstrained.
.gam_constraint_engine <- function(engine, con, method) {
  constrained <- any(lengths(con) > 0L)
  if (constrained && engine %in% c("gam", "bam"))
    stop("engine = '", engine, "' cannot impose shape constraints: ", engine,
         "() constructs the 'sc' basis but ignores its constraint matrix, so ",
         "the fit would come back unconstrained. Use engine = 'scasm', which ",
         "engine = 'auto' (the default) selects whenever 'constraints' is ",
         "given.", call. = FALSE)
  if (!constrained && !identical(engine, "scasm")) return(engine)
  .gam_scasm_check(method)
  "scasm"
}

# scasm() is a recent addition to mgcv and has no 'method' argument, so one
# passed here would be swallowed by '...' and silently ignored.
.gam_scasm_check <- function(method) {
  if (!exists("scasm", envir = asNamespace("mgcv"), inherits = FALSE))
    stop("shape constraints are fitted with mgcv::scasm(), which the ",
         "installed mgcv (", getNamespaceVersion("mgcv"), ") does not ",
         "provide. Update mgcv.", call. = FALSE)
  if (!is.null(method))
    stop("'method' does not apply to engine = 'scasm': scasm() selects the ",
         "smoothing parameters by the extended Fellner-Schall method and ",
         "takes no smoothness criterion. Drop 'method'.", call. = FALSE)
  invisible(NULL)
}

# Only the 'sc' basis carries a constraint, so a constrained smooth is built on
# it. "tp" is the package default and counts as "not chosen"; any other basis is
# an error rather than a silently dropped constraint.
.gam_constraint_bs <- function(smooths, bs, con) {
  bs_v <- rep_len(bs, length(smooths))
  has  <- lengths(con) > 0L
  bad  <- which(has & !bs_v %in% c("tp", "sc"))
  if (length(bad))
    stop("shape-constrained smooth(s) asked for a basis that cannot carry the ",
         "constraint: ",
         paste0(smooths[bad], ' (bs = "', bs_v[bad], '")', collapse = "; "),
         '.\n  Only bs = "sc" holds a shape constraint. Drop \'bs\' for these ',
         'smooths (it then defaults to "sc"), or set it to "sc".',
         call. = FALSE)
  bs_v[has] <- "sc"
  bs_v
}

# Build the 's(...)' term strings. bs / k / pc have already been validated as
# length 1 or length(smooths), so rep_len() only broadcasts a shared scalar;
# building each term explicitly (rather than one vectorised sprintf over all
# three) also means a zero-length argument can never recycle 'smooths' and
# silently emit duplicate or missing smooth terms.
# k = NULL omits the k argument, leaving mgcv's per-basis default.
# pc = NULL omits the point constraint, leaving mgcv's default sum-to-zero
# centring; it is also dropped for bases with no anchor point to constrain.
.gam_smooth_terms <- function(smooths, bs, k, pc, con = NULL) {
  n <- length(smooths)
  if (n == 0L) return(character(0))
  if (is.null(con)) con <- rep(list(character(0)), n)
  bs_v <- .gam_constraint_bs(smooths, bs, con)
  k_v  <- if (is.null(k))  rep(NA_real_, n) else rep_len(as.numeric(k),  n)
  pc_v <- if (is.null(pc)) rep(NA_real_, n) else rep_len(as.numeric(pc), n)
  pc_v[bs_v %in% .gam_bs_not_smooth] <- NA_real_
  # "+" already pins the level, and f(pc) = 0 on the boundary of f >= 0 makes
  # mgcv fail with a singular backsolve.
  pc_v[vapply(con, function(z) "+" %in% z, logical(1))] <- NA_real_
  vapply(seq_len(n), function(i) {
    a <- sprintf('%s, bs = "%s"', smooths[i], bs_v[i])
    if (length(con[[i]]))
      a <- paste0(a, ", xt = c(",
                  paste0('"', con[[i]], '"', collapse = ", "), ")")
    if (!is.na(k_v[i]))  a <- paste0(a, ", k = ",  k_v[i])
    if (!is.na(pc_v[i])) a <- paste0(a, ", pc = ", pc_v[i])
    paste0("s(", a, ")")
  }, character(1))
}

# ── Time-varying effects ─────────────────────────────────────────────────────
# A TVE replaces the constant coefficient of a statistic by a smooth function of
# time, s(t, by = STAT), so the linear predictor carries beta(t) * STAT instead
# of beta * STAT. Two properties of mgcv drive every rule below.
#
# (1) A smooth with a NUMERIC, varying 'by' variable is NOT centred: smoothCon()
#     replaces the sum-to-zero constraint by a zero-row constraint matrix
#     (mgcv:::smoothCon, "if (sd(by) > mean(by) * .Machine$double.eps * 1000)").
#     The smooth therefore carries the WHOLE coefficient function, its mean
#     included, which is exactly what makes f(t) readable as beta(t).
# (2) Because of (1) the statistic must NOT also enter linearly. Its column is
#     then exactly the constant component of the by-smooth's span, and mgcv
#     resolves the rank deficiency silently: the linear coefficient comes back
#     as 0 with standard error 0 and a NaN p-value, with no warning. See
#     .mgcv_fit(), where 'tve' statistics are removed from 'linear_terms'.

# 'baseline' in the gam list: TRUE / FALSE / a basis dimension. Returns NULL
# when the time-varying baseline is off, otherwise the k to give s(t) (NA for
# mgcv's default).
.gam_tve_baseline_spec <- function(x) {
  if (is.null(x) || isFALSE(x)) return(NULL)
  if (isTRUE(x)) return(NA_real_)
  if (is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 3)
    return(as.numeric(x))
  stop("'baseline' in the 'gam' list must be TRUE (a time-varying baseline ",
       "with mgcv's default basis dimension), FALSE (a constant baseline), or ",
       "a single basis dimension of at least 3; got an object of class '",
       paste(class(x), collapse = "/"), "' of length ", length(x), ".",
       call. = FALSE)
}

# Check validity of a time-varying effect specification
.gam_tve_check <- function(smooths, tve, tve_baseline, bs, constraints,
                           pc_supplied, stacked) {

  if (!length(tve)) {
    # 'baseline' only means anything alongside 'tve'.
    if (!is.null(tve_baseline) && !isTRUE(tve_baseline))
      stop("'baseline' was given but 'tve' is empty: 'baseline' switches the ",
           "time-varying baseline b0 + s(t) of a time-varying-effects model on ",
           "or off, and there is no such model to configure.", call. = FALSE)
    return(invisible(NULL))
  }

  if (length(smooths))
    stop("'smooths' and 'tve' cannot both be given: a fit is either a ",
         "nonlinear-effects model (s(statistic), relaxing the functional form) ",
         "or a time-varying-effects model (s(t, by = statistic), relaxing ",
         "constancy over time), not both.\n  smooths: ",
         paste(smooths, collapse = ", "), "\n  tve: ",
         paste(tve, collapse = ", "), call. = FALSE)

  if (!identical(stacked$model, "tie"))
    stop("'tve' is implemented for tie-oriented models only.", call. = FALSE)

  if (isTRUE(stacked$ordinal))
    stop("'tve' is not available for ordinal models. A time-varying effect is ",
         "a smooth of the event time, and remify() replaces the event times of ",
         "an ordinal model by their ranks, so there is no time axis to smooth ",
         "over.\n  Refit with an interval-timing remify object, or use ",
         "remwindow() to let the effects vary over windows of events.",
         call. = FALSE)

  .check_smooths_names(tve, valid = stacked$stat_names, arg = "tve")

  dup <- unique(tve[duplicated(tve)])
  if (length(dup))
    stop("'tve' names the same statistic more than once: ",
         paste(dup, collapse = ", "),
         ".\n  One s(t, by = ...) term per statistic; give a single entry per ",
         "name.", call. = FALSE)

  if ("baseline" %in% tolower(tve))
    stop("'tve' cannot name 'baseline': s(t, by = baseline) with an all-ones ",
         "column IS the time-varying baseline, which mgcv would then centre ",
         "and leave without an intercept to sit on.\n  Ask for it with ",
         "baseline = TRUE (the default) instead.", call. = FALSE)

  if (!is.null(constraints) && length(constraints))
    stop("'constraints' does not apply to 'tve': shape constraints are fitted ",
         "on the 'sc' basis through mgcv::scasm() and are a nonlinear-effects ",
         "feature, so they constrain the shape of s(statistic), not of a ",
         "coefficient function over time. Drop 'constraints'.", call. = FALSE)

  if (isTRUE(pc_supplied))
    stop("'pc' does not apply to 'tve'. A point constraint pins f(pc) = 0, and ",
         "on a time-varying effect that would force the coefficient of the ",
         "statistic to be exactly zero at t = pc rather than anchoring a level. ",
         "Drop 'pc'.", call. = FALSE)

  bad_bs <- setdiff(unique(bs), .gam_tve_bs_choices)
  if (length(bad_bs))
    stop("basis/bases not available for a time-varying effect: ",
         paste0('"', bad_bs, '"', collapse = ", "), ".\n  ",
         if ("re" %in% bad_bs)
           '"re" would place one coefficient on every distinct event time. ',
         if ("sc" %in% bad_bs)
           '"sc" only carries a shape constraint, which is a nonlinear-effects feature. ',
         "Available: ",
         paste0('"', .gam_tve_bs_choices, '"', collapse = ", "), ".",
         call. = FALSE)

  .gam_tve_baseline_spec(tve_baseline)

  df <- stacked$df
  tcol <- .gam_tve_time_col
  if (tcol %in% stacked$stat_names || tcol %in% colnames(df))
    stop("a statistic or design column is named '", tcol, "', which is the ",
         "column the time axis of a time-varying effect is written to. Rename ",
         "it in the remstats formula.", call. = FALSE)

  # mgcv only drops the centring constraint when the 'by' column actually
  # varies; reproduce its exact test so a near-constant statistic fails here
  # with the reason rather than as an unidentified level later.
  centred <- vapply(tve, function(v) {
    x <- df[[v]]; x <- x[!is.na(x)]
    !(length(x) > 1L && stats::sd(x) > mean(x) * .Machine$double.eps * 1000)
  }, logical(1))
  if (any(centred))
    stop("statistic(s) with no usable variation for a time-varying effect: ",
         paste(tve[centred], collapse = ", "),
         ".\n  mgcv centres a by-smooth whose 'by' column is constant, which ",
         "would leave the level of the coefficient function unidentified.",
         call. = FALSE)

  invisible(NULL)
}

# The time axis. reh$edgelist$time is the event time (remify converts Date /
# POSIXct input to numeric time from the origin, so this is always numeric and
# starts near 0), and df$time_index indexes reh$edgelist - so the lookup is
# row-exact for the plain, sampled and pre-stacked layouts alike, without
# assuming any particular row ordering or risk-set size.
.gam_time_column <- function(df, reh) {
  if (is.null(reh) || is.null(reh$edgelist) || is.null(reh$edgelist$time))
    stop("time-varying effects need the event times, which are taken from ",
         "reh$edgelist$time; this 'reh' does not carry them.", call. = FALSE)
  ti <- df$time_index
  if (is.null(ti))
    stop("the stacked design has no 'time_index' column, so the event times ",
         "cannot be matched to its rows.", call. = FALSE)
  tt <- as.numeric(reh$edgelist$time)
  if (max(ti, na.rm = TRUE) > length(tt) || min(ti, na.rm = TRUE) < 1L)
    stop("'time_index' in the stacked design runs outside the ", length(tt),
         " events of 'reh', so the two do not describe the same sequence.",
         call. = FALSE)
  tt[ti]
}

# Build the 's(t, ...)' term strings. Kept separate from .gam_smooth_terms() so
# the pc / xt / 'sc' logic of the nonlinear path stays untouched: a TVE term
# never carries a point constraint or a shape constraint.
# The baseline smooth has no 'by' and is therefore centred by mgcv, which is
# what keeps 'baseline' interpretable as the row-average log-rate. It uses the
# first entry of 'bs'; its k comes from 'baseline', not from 'k'.
.gam_tve_terms <- function(tve, bs, k, tve_baseline,
                           tcol = .gam_tve_time_col) {
  n <- length(tve)
  if (n == 0L) return(character(0))
  bs_v <- rep_len(bs, n)
  k_v  <- if (is.null(k)) rep(NA_real_, n) else rep_len(as.numeric(k), n)

  terms <- vapply(seq_len(n), function(i) {
    a <- sprintf('%s, bs = "%s", by = %s', tcol, bs_v[i], tve[i])
    if (!is.na(k_v[i])) a <- paste0(a, ", k = ", k_v[i])
    paste0("s(", a, ")")
  }, character(1))

  kb <- .gam_tve_baseline_spec(tve_baseline)
  if (!is.null(kb)) {
    a <- sprintf('%s, bs = "%s"', tcol, bs_v[1L])
    if (!is.na(kb)) a <- paste0(a, ", k = ", kb)
    terms <- c(paste0("s(", a, ")"), terms)
  }
  terms
}

# Basis-dimension check for the TVE path. The smoothed covariate is TIME, so the
# ceiling is the number of distinct event times - NOT the number of distinct
# values of the 'by' statistic, which is typically far smaller (a small count
# such as inertia) and would reject the default k = 10 out of hand.
.remstimate_gam_tve_k_check <- function(n_unique_t, tve, k, tve_baseline) {
  n <- length(tve)
  if (!n) return(invisible(NULL))
  k_eff <- if (is.null(k)) rep(10, n) else rep_len(as.numeric(k), n)
  kb    <- .gam_tve_baseline_spec(tve_baseline)

  # k = NULL leaves mgcv's own default of 10 for these bases, so the ceiling
  # still applies and is better reported here than from inside smoothCon().
  bad <- which(k_eff > n_unique_t)
  if (length(bad))
    stop("time-varying effect(s) with a basis dimension above the number of ",
         "distinct event times (", n_unique_t, "): ",
         paste0(tve[bad], " (k = ", k_eff[bad], ")", collapse = "; "),
         ".\n  Reduce 'k'.", call. = FALSE)

  if (is.null(kb)) return(invisible(NULL))
  kb_eff <- if (is.na(kb)) 10 else kb
  if (kb_eff > n_unique_t)
    stop("the time-varying baseline needs ", kb_eff, " basis functions ",
         if (is.na(kb)) "(mgcv's default) " else "",
         "but there are only ", n_unique_t, " distinct event times.\n  Set ",
         "baseline = FALSE or a smaller basis dimension.", call. = FALSE)
  if (kb_eff > n_unique_t / 2)
    warning("the time-varying baseline uses ", kb_eff, " basis functions for ",
            n_unique_t, " risk sets, which approaches a separate baseline per ",
            "risk set. The smoothing penalty still applies, but consider a ",
            "smaller 'baseline'.", call. = FALSE)
  invisible(NULL)
}


# mgcv::bam cannot fit general families ("general families not supported by
# bam"), and cox.ph is one - so an ordinal GAM is gam(), or scasm() when shape
# constraints were asked for (scasm() handles cox.ph, see ?scasm).
.gam_ordinal_engine <- function(engine) {
  if (identical(engine, "bam"))
    stop("engine = 'bam' is not available for ordinal models: an ordinal GAM ",
         "is fit with mgcv's 'cox.ph' family, and bam() supports no general ",
         "families. Use engine = 'gam' (the default).", call. = FALSE)
  if (identical(engine, "scasm")) return("scasm")
  "gam"
}

# General families (cox.ph) and scasm() leave gam()$converged empty: scasm
# reports a logical in outer.info$converged, gam() a message in outer.info$conv.
# With no smooth terms there is no outer problem, so only finiteness is left.
.gam_converged <- function(fit) {
  if (!is.null(fit$converged)) return(isTRUE(fit$converged))
  if (is.logical(fit$outer.info$converged))
    return(isTRUE(fit$outer.info$converged))
  conv <- fit$outer.info$conv
  if (!is.null(conv)) return(isTRUE(grepl("full convergence", conv, fixed = TRUE)))
  if (!length(fit$sp)) return(all(is.finite(stats::coef(fit))))
  NA
}

# scasm() puts its Fellner-Schall iteration count in outer.info, not $iter.
.gam_iterations <- function(fit) fit$iter %||% fit$outer.info$niter

# A constraint the data run against pins every coefficient of the smooth to an
# active constraint, and mgcv 1.9.4 then dies in its covariance projection with
# a bare "subscript out of bounds". Translate that into the cause.
.gam_scasm_try <- function(fit_call, con, smooths) {
  tryCatch(fit_call(), error = function(e) {
    msg <- conditionMessage(e)
    if (!grepl("subscript out of bounds|not positive definite|singular", msg))
      stop(e)
    lab <- smooths[lengths(con) > 0L]
    stop("mgcv::scasm() failed with \"", msg, "\" while fitting the shape-",
         "constrained smooth(s): ", paste(lab, collapse = ", "),
         ".\n  That is what mgcv reports when a constraint binds over the ",
         "whole range of the statistic: the data run against it, so the fit ",
         "is held on the edge of the feasible set with nothing left to ",
         "estimate.",
         "\n  Look at the unconstrained shape first (the same model without ",
         "'constraints'), then reverse or relax the constraint, or lower 'k'.",
         call. = FALSE)
  })
}

# Under ordinal timing each risk set has its own baseline hazard, which is
# profiled out of the partial likelihood. A statistic that never varies WITHIN a
# risk set is therefore absorbed by it and carries no information. Left in, such a
# column makes mgcv fail with a bare "non-conformable arguments". The shared
# .drop_constant_stats() only removes globally constant columns (and exempts
# 'baseline' by design, since the interval backends need it), so the
# within-stratum case is caught here.
# An aliased statistic the caller never named is dropped with a warning. One
# named in 'smooths' is an error instead: silently discarding a requested smooth
# would answer a different question than the one asked, and 'bs' / 'k' / 'pc'
# are positional, so dropping an entry from 'smooths' would shift every
# remaining setting onto the wrong statistic.
.gam_drop_aliased_ordinal <- function(df, stat_names, smooths) {
  present <- intersect(stat_names, colnames(df))
  aliased <- present[vapply(present, function(v) {
    x <- df[[v]]
    if (all(is.na(x))) return(TRUE)
    # one distinct value per risk set. Deliberately type-generic: a pre-stacked
    # design can carry a categorical statistic (e.g. an event-type factor).
    n_uniq <- tapply(x, df$time_index,
                     function(z) length(unique(z[!is.na(z)])))
    all(n_uniq <= 1L, na.rm = TRUE)
  }, logical(1))]
  named <- intersect(aliased, smooths)
  if (length(named))
    stop("'smooths' names statistic(s) that are constant within every risk ",
         "set: ", paste(named, collapse = ", "),
         ".\n  Under ordinal timing the risk-set baseline hazard absorbs them, ",
         "so no shape is identifiable. Drop them from 'smooths'.",
         call. = FALSE)
  if (length(aliased))
    warning("Statistic(s) constant within every risk set: ",
            paste(aliased, collapse = ", "),
            ". Absorbed by the risk-set baseline hazard under ordinal timing ",
            "(not identifiable); dropped from the fit.", call. = FALSE)
  setdiff(stat_names, aliased)
}

# cox.ph handles ties by Peto's approximation, which is clogit(method =
# "approximate") - exact only with a single event per risk set. A stacked REM
# normally has exactly that, but 'obs' comes from tabulate() and counts
# simultaneous events on one dyad, so more than one event per time point is
# possible and the reported partial likelihood is then approximate.
.gam_ordinal_tie_check <- function(df) {
  n_ev <- tapply(df$obs, df$time_index, sum)
  if (any(n_ev > 1, na.rm = TRUE))
    warning(sum(n_ev > 1, na.rm = TRUE), " of ", length(n_ev), " risk sets ",
            "contain more than one event. mgcv's cox.ph handles ties with ",
            "Peto's approximation, so the partial likelihood (and the ",
            "information criteria built on it) are approximate, as in ",
            "survival::clogit(method = 'approximate').", call. = FALSE)
  invisible(NULL)
}

# A point constraint constrains the smooth at f(pc) = 0, which is what makes
# 'baseline' interpretable as the log-rate when the statistic equals pc. That
# only holds if pc is inside the observed range: outside it, mgcv still fits
# but the anchor is extrapolated and 'baseline' becomes unstable (moving from
# -4.95 to -9.27 on the test design). Warn rather than stop - an out-of-range
# anchor is unusual but not intrinsically invalid.
# Under ordinal timing there is no 'baseline' coefficient to anchor - the
# risk-set baseline hazards absorb any constant - so pc only sets the level at
# which the smooth is drawn, and the warning says so.
.remstimate_gam_pc_check <- function(stacked, smooths, bs, pc, ordinal = FALSE,
                                     con = NULL) {
  if (is.null(pc) || !length(smooths)) return(invisible(NULL))
  if (is.null(con)) con <- rep(list(character(0)), length(smooths))
  pc_v <- rep_len(as.numeric(pc), length(smooths))
  bs_v <- rep_len(bs, length(smooths))
  bad  <- character(0)
  for (i in seq_along(smooths)) {
    if (bs_v[i] %in% .gam_bs_not_smooth) next   # no pc emitted for these
    if ("+" %in% con[[i]]) next                 # nor for a positive smooth
    x <- stacked$df[[smooths[i]]]
    x <- x[!is.na(x)]
    if (!length(x)) next
    if (pc_v[i] < min(x) || pc_v[i] > max(x))
      bad <- c(bad, sprintf("%s (pc = %g, observed range %g to %g)",
                            smooths[i], pc_v[i], min(x), max(x)))
  }
  if (length(bad))
    warning("point constraint outside the observed range of: ",
            paste(bad, collapse = "; "),
            if (ordinal)
              paste0(".\n  The smooth is anchored by extrapolation, so the ",
                     "level it is drawn at is not the effect at an observed ",
                     "value. Consider pc = NULL.")
            else
              paste0(".\n  The smooth is anchored by extrapolation, so ",
                     "'baseline' is not the log-rate at an observed value. ",
                     "Consider pc = NULL."),
            call. = FALSE)
  invisible(NULL)
}

# A basis dimension above the number of distinct values a statistic takes in the
# stacked design cannot be built: mgcv aborts with "A term has fewer unique
# covariate combinations than specified maximum degrees of freedom". Many REM
# statistics are small counts (inertia, reciprocity on a short history), so this
# is a routine user error rather than an exotic one.
.remstimate_gam_k_check <- function(stacked, smooths, bs, k) {
  if (!length(smooths)) return(invisible(NULL))

  n_uniq <- vapply(smooths, function(v)
    length(unique(stacked$df[[v]][!is.na(stacked$df[[v]])])), integer(1))

  # k = NULL means mgcv picks the basis dimension, which is 10 for the common
  # 1-d spline bases - so the check still applies, and fails here with an
  # actionable message rather than later inside smoothCon(). 're' sizes itself
  # from the data (one coefficient per level) and carries no such ceiling, so it
  # is exempt.
  # bs / k are already length 1 or length(smooths) (.remstimate_gam_spec_check),
  # so rep_len() only broadcasts a shared scalar - it can never mask a
  # mismatch here.
  bs_v    <- rep_len(bs, length(smooths))
  k_eff   <- if (is.null(k)) rep(10L, length(smooths))
             else rep_len(k, length(smooths))
  checked <- !bs_v %in% .gam_bs_not_smooth

  error_index <- which(checked & n_uniq < k_eff)

  if (length(error_index) > 0) {
    stop("smooth term(s) with fewer unique values than the basis dimension k: ",
         paste0(smooths[error_index], " (", n_uniq[error_index],
                " unique, k = ", k_eff[error_index], ")", collapse = "; "),
         ".\n  Reduce 'k' for these terms, or drop them from 'smooths'.",
         call. = FALSE)
  }
  invisible(NULL)
}

# 'smooths' and 'tve' are the two (mutually exclusive) term specifications, both
# optional: with neither, this fits the plain stacked Poisson / cox.ph GLM.
# 'reh' is only read when 'tve' is given, to build the time axis.
.mgcv_fit <- function(stacked, reh = NULL, smooths = NULL, tve = NULL,
                      tve_baseline = TRUE, bs, k, pc = 0, con = NULL, select,
                      method = NULL, engine = "auto", ...) {

  if (is.null(con))
    con <- stats::setNames(rep(list(character(0)), length(smooths)), smooths)

  df <- stacked$df
  stat_names <- stacked$stat_names
  ordinal <- stacked$ordinal
  tve_mode <- length(tve) > 0L

  # Case-control sampling correction: offset += -log(pi_d) per row.
  # .remstimate_make_stack always supplies 'samp_offset' for a tie design - all
  # zeros when the stats are unsampled, and -log(pi_d) taken from
  # attr(stats, "samp_prob") when they are (cases 0, controls log((D-1)/K)).
  # Use that canonical column rather than recomputing log(weight): the two
  # agree for freshly stacked sampled stats, but 'weight' is absent entirely
  # for unsampled stats, and a pre-stacked 'remstats_stacked' object can carry
  # a samp_offset that was not derived from 'weight'.
  if (is.null(df$samp_offset)) df$samp_offset <- 0

  if (length(stat_names) == 0L)
    stop("No statistics found in the stacked design.", call. = FALSE)

  # An ordinal design carries no 'baseline' / 'log_interevent' column, and any
  # statistic that is constant within a risk set is absorbed by that risk set's
  # baseline hazard - drop those before the k / pc checks inspect them.
  if (ordinal) {
    stat_names <- .gam_drop_aliased_ordinal(df, stat_names, smooths)
    if (!length(stat_names))
      stop("No identifiable statistics left in the ordinal design: every ",
           "statistic is constant within the risk sets, so the conditional ",
           "likelihood carries no information about it.", call. = FALSE)
  }

  if (tve_mode) {

    # .gam_tve_check() has already rejected an ordinal design; assert it here
    # too, since the branch below would otherwise build an unidentifiable model.
    if (ordinal)
      stop("time-varying effects are not available for ordinal models.",
           call. = FALSE)

    df[[.gam_tve_time_col]] <- .gam_time_column(df, reh)
    .remstimate_gam_tve_k_check(length(unique(df[[.gam_tve_time_col]])),
                                tve, k, tve_baseline)

    smooth_terms <- .gam_tve_terms(tve, bs, k, tve_baseline)
    # A statistic with a time-varying coefficient must NOT also enter linearly:
    # s(t, by = STAT) is uncentred, so its span already contains the constant
    # effect, and mgcv resolves the resulting rank deficiency silently (the
    # linear coefficient comes back as 0 with standard error 0).
    linear_terms <- setdiff(stat_names, tve)
    constrained  <- character(0)

  } else {

    .remstimate_gam_k_check(stacked, smooths, bs, k)
    .remstimate_gam_pc_check(stacked, smooths, bs, pc, ordinal, con)

    smooth_terms <- .gam_smooth_terms(smooths, bs, k, pc, con)
    linear_terms <- setdiff(stat_names, smooths)
    constrained  <- smooths[lengths(con) > 0L]
  }

  if (ordinal) {

    # ── Ordinal: stratified Cox partial likelihood (= conditional logit) ─────
    # mgcv's cox.ph family takes a two-column response (time, stratum index) and
    # the event indicator through 'weights'. With a constant dummy time every
    # risk set is a single stratum whose members are all at risk at that one
    # time, so the partial likelihood is exactly the conditional-logit
    # likelihood the C++ MLE and the clogit/coxme backends maximise. The
    # stratum-specific baseline hazards are profiled out, so - unlike the
    # interval branch - there is no intercept to fit and none is needed.
    mgcv_ver <- package_version(getNamespaceVersion("mgcv"))
    if (mgcv_ver < "1.9.0")
      stop("Ordinal models need the stratified 'cox.ph' family from mgcv ",
           ">= 1.9.0; the installed version is ", mgcv_ver, ".", call. = FALSE)

    engine <- .gam_ordinal_engine(engine)
    if (engine != "scasm") {
      if (is.null(method)) method <- "REML"
      if (identical(method, "fREML"))
        stop("method = 'fREML' is a bam() criterion, and an ordinal GAM is fit ",
             "with gam(). Use method = 'REML'.", call. = FALSE)
    }

    .gam_ordinal_tie_check(df)

    df$.surv_time <- 1
    formula_obj <- stats::as.formula(paste0(
      "cbind(.surv_time, time_index) ~ offset(samp_offset) + ",
      paste(c(smooth_terms, linear_terms), collapse = " + ")
    ))

    fit <- if (engine == "scasm")
      .gam_scasm_try(function()
        mgcv::scasm(formula_obj, family = mgcv::cox.ph(), data = df,
                    weights = df$obs, select = select, ...),
        con, smooths)
    else
      mgcv::gam(formula_obj, family = mgcv::cox.ph(), data = df,
                weights = df$obs, select = select, method = method, ...)

  } else {

    # ── Interval: Poisson GAM ────────────────────────────────────────────────
    # Require an all-ones column (baseline) among the terms that enter linearly
    if (!any(vapply(linear_terms, function(v) {
          x <- df[[v]][!is.na(df[[v]])]
          length(x) > 0L && all(x == 1)
        }, logical(1))))
      stop("The GAM backend requires an intercept column ('baseline') among the ",
           "statistics that enter linearly: mgcv constrains the ",
           if (tve_mode) "time-varying baseline" else "smooth terms",
           ", so nothing else can carry the overall event rate.\n  Rebuild the ",
           "statistics with the baseline included (include '~ 1 + ...' in ",
           "remstats), and keep it out of '", if (tve_mode) "tve" else "smooths",
           "'.", call. = FALSE)

    formula_obj <- stats::as.formula(paste0(
      "obs ~ -1 + offset(log_interevent + samp_offset) + ",
      paste(c(smooth_terms, linear_terms), collapse = " + ")
    ))

    # "auto": gam() below the threshold, bam() above it. Shape constraints have
    # already resolved the engine to "scasm" upstream, so the size rule never
    # overrides them.
    if (engine == "auto")
      engine <- if (nrow(df) > 1e5L) "bam" else "gam"

    # method = NULL: each engine's own name for restricted maximum likelihood.
    if (is.null(method) && engine != "scasm")
      method <- if (engine == "bam") "fREML" else "REML"

    fit <- if (engine == "scasm")
      .gam_scasm_try(function()
        mgcv::scasm(formula_obj, family = stats::poisson(), data = df,
                    select = select, ...),
        con, smooths)
    else if (engine == "bam")
      mgcv::bam(formula_obj, family = stats::poisson(), data = df,
                select = select, method = method, ...)
    else
      mgcv::gam(formula_obj, family = stats::poisson(), data = df,
                select = select, method = method, ...)
  }

  coefs <- stats::coef(fit)
  P     <- length(coefs)

  # mgcv drops incomplete rows, so off_obs, M and the null model must be taken
  # from the rows it actually used, otherwise loglik is off by their log(dt).
  used <- seq_len(nrow(df))
  if (!is.null(fit$na.action)) used <- used[-as.integer(fit$na.action)]
  df_fit <- df
  if (length(used) < nrow(df)) {
    warning(nrow(df) - length(used), " of ", nrow(df), " rows of the stacked ",
            "design were dropped because a statistic was NA. The fit, the ",
            "log-likelihood and the information criteria refer to the ",
            "remaining rows only.", call. = FALSE)
    df_fit <- df[used, , drop = FALSE]
  }

  M <- length(unique(df_fit$time_index))

  # Ordinal: the cox.ph partial loglik IS the conditional-logit loglik (it
  # matches clogit's loglik[2] exactly), so it needs no offset correction.
  # Interval: Poisson loglik -> REM loglik by dropping the offset constant
  # Σ log(dt) over the events.
  ll_obj     <- stats::logLik(fit)
  off_obs    <- if (ordinal) 0 else sum(df_fit$log_interevent[df_fit$obs == 1L])
  loglik_val <- as.numeric(ll_obj) - off_obs
  # mgcv carries two effective-df vectors: fit$edf (trace of the influence
  # matrix - the per-smooth edf that summary.gam() tabulates) and fit$edf2
  # (Wood, Pya & Safken 2016: edf inflated to account for smoothing-parameter
  # uncertainty)
  edf_total <- sum(fit$edf)
  if (!length(fit$edf) || !is.finite(edf_total)) edf_total <- P
  edf_corrected <- if (is.null(fit$edf2)) NA_real_ else sum(fit$edf2)
  edf_parametric <- if (is.null(fit$nsdf) || fit$nsdf < 1L) 0
                    else sum(fit$edf[seq_len(fit$nsdf)])

  vcov_mat <- tryCatch(stats::vcov(fit),
                       error = function(e) matrix(NA_real_, P, P))
  se <- sqrt(diag(vcov_mat)); names(se) <- names(coefs)

  if (ordinal) {
    # beta = 0 conditional logit: -Σ_m log Σ_{d in R_m} exp(samp_offset_d),
    # which is -Σ_m log|R_m| for unsampled statistics.
    # Only risk sets that actually contain an event contribute: a stratum with
    # no case drops out of a partial likelihood entirely
    ev   <- tapply(df_fit$obs, df_fit$time_index, sum)
    wsum <- tapply(exp(df_fit$samp_offset), df_fit$time_index, sum)
    null_loglik <- -sum(log(wsum[ev > 0]))
  } else {
    # baseline-intercept null (matches C++), same offset correction, same rows
    null_fit    <- stats::glm(obs ~ offset(log_interevent + samp_offset),
                              family = stats::poisson(), data = df_fit)
    null_loglik <- as.numeric(stats::logLik(null_fit)) - off_obs
  }

  aic_val  <- 2 * edf_total - 2 * loglik_val
  bic_val  <- edf_total * log(M) - 2 * loglik_val    # n = E, not nrow
  aicc_val <- aic_val +
    2 * edf_total * (edf_total + 1) / max(M - edf_total - 1, 1)

  where_bl <- which(tolower(names(coefs)) == "baseline")
  where_bl <- if (length(where_bl)) where_bl[1L] else NULL

  # map of smooth terms to the statistics they were built from
  smooth_map     <- .gam_smooth_map(fit)
  par_idx        <- if (is.null(fit$nsdf) || fit$nsdf < 1L) integer(0)
                    else seq_len(fit$nsdf)
  smoothed_stats <- if (is.null(smooth_map)) character(0)
                    else intersect(stat_names, smooth_map$statistic)

  res <- list(
    coefficients      = coefs,
    tve               = if (tve_mode) intersect(tve, smoothed_stats)
                        else character(0),
    loglik            = loglik_val,
    gradient          = NULL,
    hessian           = NULL,
    vcov              = vcov_mat,
    se                = se,
    residual.deviance = -2 * loglik_val,
    null.deviance     = -2 * null_loglik,
    model.deviance    = -2 * null_loglik - (-2 * loglik_val),
    AIC               = aic_val,
    AICC              = aicc_val,
    BIC               = bic_val,
    df.null           = M,
    df.model          = edf_total,
    df.model.corrected = edf_corrected,
    df.parametric     = edf_parametric,
    df.residual       = M - edf_total,
    converged         = .gam_converged(fit),
    iterations        = .gam_iterations(fit),
    smooth_map        = smooth_map,
    smoothed          = smoothed_stats,
    linear            = setdiff(stat_names, smoothed_stats),
    parametric_index  = par_idx,
    parametric        = coefs[par_idx],
    parametric_se     = se[par_idx],
    stacked_data      = df,
    backend_fit       = fit
  )

  structure(
    res,
    class             = c("remstimate_gam", "remstimate"),
    formula           = formula_obj,
    model             = "tie",
    approach          = "Frequentist",
    method            = "GAM",
    engine            = engine,
    family            = if (ordinal) "cox.ph" else "poisson",
    ordinal           = ordinal,
    statistics        = stat_names,
    gam_type          = if (tve_mode) "tve" else "nle",
    smooths           = if (tve_mode) character(0) else smoothed_stats,
    tve               = if (tve_mode) intersect(tve, smoothed_stats)
                        else character(0),
    time_var          = if (tve_mode) .gam_tve_time_col else NULL,
    time_range        = if (tve_mode) range(df[[.gam_tve_time_col]],
                                           na.rm = TRUE) else NULL,
    constraints       = con[intersect(names(con), smoothed_stats)],
    constrained       = intersect(constrained, smoothed_stats),
    coef_contract_1to1 = is.null(fit$smooth) || !length(fit$smooth),
    where_is_baseline = where_bl,
    ncores            = 1L
  )
}

# Check whether the fitted GAM's coefficient vector holds one entry per statistic
.gam_coef_contract <- function(object) {
  isTRUE(attr(object, "coef_contract_1to1"))
}

# Time-varying-effects fit (s(t, by = ...)) rather than a nonlinear-effects one.
.gam_is_tve <- function(object) identical(attr(object, "gam_type"), "tve")


# ── S3 methods ───────────────────────────────────────────────────────────────
# NB: these are deliberately NOT copies of the GLMNET methods. A GAM's
# coefficient vector holds one entry per BASIS function of every smooth
# (e.g. 's(inertia).1' ... 's(inertia).9'), not one per statistic, so the
# GLMNET/SHRINKEM pattern of dropping the subclass and dispatching the shared
# remstimate methods does not apply here - see diagnostics below.

#' @export
print.remstimate_gam <- function(x, ...) {
  fit <- x$backend_fit
  cat("REM - ", attr(x, "model"), " model - GAM [mgcv::",
      attr(x, "engine") %||% "gam", ", ",
      if (isTRUE(attr(x, "ordinal")))
        "cox.ph: stratified partial likelihood (conditional logit)"
      else "poisson: interval likelihood",
      "]\n\n", sep = "")

  is_tve <- .gam_is_tve(x)
  sm <- .gam_smooth_table(fit)
  if (!is.null(sm)) {
    cat(if (is_tve) "Time-varying effects (smooths of time):\n"
        else "Smooth terms (nonlinear effects):\n")
    print(sm, row.names = FALSE)
    if (!is.null(sm$shape))
      cat("  shape: m+/m- monotone increasing/decreasing,",
          "c+/c- convex/concave, + positive\n")
    if (is_tve) {
      tr <- attr(x, "time_range")
      if (!is.null(tr))
        cat("  s(", .gam_tve_time_col, ") is the time-varying baseline; ",
            "s(", .gam_tve_time_col, "):X is the coefficient of X over time ",
            sprintf("(t in [%g, %g])", tr[1L], tr[2L]), "\n", sep = "")
    }
    cat("\n")
  }

  par_coefs <- .gam_parametric_coefs(x)
  if (length(par_coefs)) {
    cat("Parametric (linear) coefficients:\n")
    print(par_coefs)
  }

  cat("\nNull deviance:", x$null.deviance,
      "\nResidual deviance:", x$residual.deviance, "\n")
  cat("AIC:", x$AIC, "AICC:", x$AICC, "BIC:", x$BIC, "\n")
  cat("Effective df (sum of edf):", round(x$df.model, 3), "\n")
  if (!.gam_coef_contract(x))
    cat("\n",
        if (is_tve)
          "Time-varying statistic(s) have no single coefficient; their whole effect\nis the curve beta(t), drawn by plot(which = 9) and tested in summary().\n"
        else
          "Smoothed statistic(s) have no single coefficient; their shape is in\nplot(which = 9) and their whole-term test in summary().\n",
        sep = "")
  cat("\n")
  invisible(x)
}

#' @export
summary.remstimate_gam <- function(object, ...) {
  summary.remstimate(object, ...)
}

# coef() has to choose between two vectors that are both defensible, so it makes
# the choice explicit instead of guessing:
#   type = "all"        every fitted coefficient, i.e. the parametric block
#                       followed by one entry per basis function of each smooth
#                       ('s(inertia).1' ... 's(inertia).4'). Identical to
#                       object$coefficients and to coef() on the mgcv fit, so
#                       coef() never disagrees with the stored vector.
#   type = "parametric" only the statistics that entered linearly, which is the
#                       "one number per statistic" vector the rest of the
#                       package expects - but note it OMITS the smoothed
#                       statistics entirely, since a smooth has no single
#                       coefficient. Its effect is summarised by edf and the
#                       whole-term test in summary(), and drawn by plot(which=9).
# The default is "all": it is what is stored, and a silent switch to a shorter
# vector would be indexed against attr(object, "statistics") by callers that
# assume the usual 1:1 contract.
#' @export
coef.remstimate_gam <- function(object, type = c("all", "parametric"), ...) {
  type <- match.arg(type)
  if (type == "parametric") return(.gam_parametric_coefs(object))
  object$coefficients
}

#' @export
logLik.remstimate_gam <- function(object, ...) object$loglik

# Statistic <-> coefficient mapping for a fitted GAM. This is the object that
# makes the broken 1:1 contract explicit: attr(object, "statistics") has one
# entry per statistic, while object$coefficients has one entry per BASIS
# function of every smooth, so downstream consumers cannot index one by the
# other. One row per smooth term:
#   statistic  the model statistic that was smoothed (as named in 'statistics')
#   term       the mgcv label, e.g. 's(inertia)'
#   bs         smoothing basis actually constructed
#   k          basis dimension
#   edf        effective degrees of freedom of the whole smooth
#   shape      mgcv shape-constraint codes imposed on the smooth, "" if none
#   first/last index block this smooth occupies in object$coefficients
# Returns NULL when the model has no smooth terms (a GAM with an empty
# 'smooths' is an ordinary Poisson GLM and its contract IS 1:1).
.gam_smooth_map <- function(fit) {
  if (is.null(fit) || !length(fit$smooth)) return(NULL)
  do.call(rbind, lapply(fit$smooth, function(sm) {
    idx <- sm$first.para:sm$last.para
    # A time-varying effect is s(t, by = STAT): the smoothed covariate is time
    # and the STATISTIC it belongs to is the 'by' variable. mgcv writes the
    # string "NA" (not a missing value) when there is no by variable.
    by_var <- if (!is.null(sm$by) && !identical(sm$by, "NA")) sm$by
              else NA_character_
    data.frame(
      statistic = if (is.na(by_var)) paste(sm$term, collapse = ", ") else by_var,
      by        = by_var,
      term      = sm$label,
      bs        = sub("\\.smooth$", "", class(sm)[1L]),
      k         = sm$bs.dim,
      edf       = round(sum(fit$edf[idx]), 3),
      shape     = if (is.character(sm$xt)) paste(sm$xt, collapse = " ") else "",
      first     = sm$first.para,
      last      = sm$last.para,
      stringsAsFactors = FALSE
    )
  }))
}

# Display subset of the mapping, used by print.remstimate_gam. The 'shape'
# column only earns its space when something was actually constrained.
.gam_smooth_table <- function(fit) {
  map <- .gam_smooth_map(fit)
  if (is.null(map)) return(NULL)
  cols <- c("term", "bs", "k", "edf")
  if (any(nzchar(map$shape))) cols <- c(cols, "shape")
  map[, cols, drop = FALSE]
}

# The parametric (non-smooth) block of the coefficient vector, i.e. the
# statistics that entered the formula linearly. fit$nsdf is mgcv's count of
# parametric terms and they always occupy the leading coefficients.
.gam_parametric_coefs <- function(object) {
  fit <- object$backend_fit
  cf  <- object$coefficients
  if (is.null(fit) || is.null(fit$nsdf) || fit$nsdf < 1L) return(cf[0])
  cf[seq_len(fit$nsdf)]
}


# Fitted linear predictor WITHOUT the offset, row-aligned to the design.
# predict.gam(type = "link") returns eta + offset; both parts must come off
# before ranking. log_interevent is constant within a risk set and would cancel,
# but samp_offset is not (0 for the case, log((D-1)/K) for each control), so
# leaving it in would shift the case relative to its own controls.
# Also copies 'time_index' into 'time': .diagnostics_recall() keys risk sets on
# df$time, which stack_stats() supplies for the durem/sampled layouts but not
# for the plain tie stack.
.gam_linear_predictor <- function(fit, df) {
  lp <- tryCatch(as.numeric(stats::predict(fit, type = "link")),
                 error = function(e) NULL)
  if (is.null(lp)) return(NULL)
  keep <- seq_len(nrow(df))          # rows mgcv dropped (NA) are absent from lp
  if (!is.null(fit$na.action)) keep <- keep[-as.integer(fit$na.action)]
  if (length(keep) != length(lp)) return(NULL)
  df  <- df[keep, , drop = FALSE]
  off <- rep(0, length(lp))
  if (!is.null(df$log_interevent)) off <- off + df$log_interevent
  if (!is.null(df$samp_offset))    off <- off + df$samp_offset
  if (is.null(df$time)) df$time <- df$time_index
  if (is.null(df$time)) return(NULL)
  list(lp = lp - off, df = df)
}

# GAM diagnostics. The GLMNET route - drop the subclass and dispatch
# diagnostics.remstimate - cannot work here: that route hands
# object$coefficients to computeDiagnostics(), which forms the linear predictor
# as stats %*% pars and so needs one coefficient per statistic. A GAM has one
# per basis function, and stats %*% pars is not this model's linear predictor
# anyway. So diagnostics follow the GLMM pattern: take the fitted linear
# predictor from mgcv (smooth contributions included by construction) and feed
# it to the shared .diagnostics_recall(), giving recall, relative ranks,
# probability ratios, log-loss and the per-type breakdown as for every other
# backend. The Schoenfeld residuals (which = 2) and the C++ rates behind the
# waiting-time Q-Q (which = 1) are per-coefficient quantities with no GAM
# analogue; they are simply absent and plot.diagnostics skips them. The honest
# substitutes are plot(which = 9) and plot(which = 10).
#' @export
#' @method diagnostics remstimate_gam
diagnostics.remstimate_gam <- function(object, reh, stats = NULL,
                                       top_pct = 0.05,
                                       surprise_threshold = 0.2, ...) {
  if (!identical(attr(object, "model"), "tie"))
    stop("GAM diagnostics are implemented for tie-oriented models only.",
         call. = FALSE)
  if (inherits(reh, "remify_durem") || inherits(stats, "remstats_durem"))
    stop("The GAM backend does not support duration models, so no durem ",
         "diagnostics are available.", call. = FALSE)

  df  <- object$stacked_data
  fit <- object$backend_fit
  if (is.null(df) || is.null(fit))
    stop("No stacked data / backend fit stored in the GAM object.", call. = FALSE)

  pred <- .gam_linear_predictor(fit, df)
  if (is.null(pred))
    stop("Could not obtain a linear predictor from the mgcv fit, so recall ",
         "cannot be computed.", call. = FALSE)

  out <- .diagnostics_recall(pred$lp, pred$df, top_pct)
  if (is.null(out$recall))
    stop("Recall could not be computed from the stacked design.", call. = FALSE)

  # .plot_recall keys the x-axis on 'event': add a dense index over the time
  # groups so each recall table plots like every other one.
  .add_event <- function(rc) {
    if (!is.null(rc) && !is.null(rc$per_event))
      rc$per_event$event <-
        match(rc$per_event$time, sort(unique(rc$per_event$time)))
    rc
  }
  out$recall <- .add_event(out$recall)
  if (!is.null(out$recall_by_type))
    out$recall_by_type <- lapply(out$recall_by_type, .add_event)

  out$surprises          <- .surprises_from_recall(out$recall, surprise_threshold)
  out$surprise_threshold <- surprise_threshold
  out$smooth_map         <- object$smooth_map
  # With several smooths of the same covariate - a time-varying baseline s(t)
  # alongside one or more s(t, by = X) - the terms can be close to mutually
  # predictable, which leaves the individual curves poorly determined even
  # though the fit as a whole is fine. mgcv measures exactly that.
  if (length(fit$smooth) > 1L)
    out$concurvity <- tryCatch(mgcv::concurvity(fit), error = function(e) NULL)
  out$.reh.processed     <- denormalize_reh(reh)
  class(out) <- c("diagnostics_gam", "diagnostics", "remstimate")
  out
}

#' @export
#' @method print diagnostics_gam
print.diagnostics_gam <- function(x, ...) {
  reh <- x$.reh.processed
  cat("Diagnostics - Relational Event Model (GAM [mgcv])\n")
  if (!is.null(reh)) cat(sprintf("Actors: %d  Events: %d\n", reh$N, reh$M))

  if (!is.null(x$smooth_map)) {
    cols <- c("term", "bs", "k", "edf")
    if (any(nzchar(x$smooth_map$shape))) cols <- c(cols, "shape")
    cat("\nSmooth terms:\n")
    print(x$smooth_map[, cols, drop = FALSE], row.names = FALSE)
  }
  if (!is.null(x$concurvity) && "worst" %in% rownames(x$concurvity)) {
    cat("\nConcurvity (worst; 0 = independent, 1 = fully confounded):\n")
    print(round(x$concurvity["worst", ], 3))
  }
  if (!is.null(x$recall)) {
    cat("\nRecall:\n")
    .print_recall_summary(x$recall, "Tie model",
                          x$surprises, x$surprise_threshold)
  }
  if (!is.null(x$recall_by_type)) {
    cat("\nRecall by type:\n")
    for (tp in names(x$recall_by_type))
      .print_recall_summary(x$recall_by_type[[tp]], tp)
  }
  cat("\nSchoenfeld residuals and the fitted-rate Q-Q have no GAM analogue; ",
      "see\nplot(which = 9) for the smooth shapes and plot(which = 10) for ",
      "mgcv's\nresidual / basis-dimension check.\n", sep = "")
  invisible(x)
}

# which 1-8: shared plot.diagnostics interface (recall = 3, prob-ratio = 6,
#            per-type = 8). Panels 1-2 are absent from a GAM diagnostics object
#            and are skipped, so asking for them warns rather than drawing
#            nothing.
# which 9:   partial-effect curves - the reason to fit a GAM at all.
# which 10:  mgcv's residual / basis-dimension check.
#' @export
#' @method plot remstimate_gam
plot.remstimate_gam <- function(x, reh, diagnostics = NULL,
                                which = c(3, 9), effects = NULL,
                                sender_effects = NULL, receiver_effects = NULL,
                                ...) {
  fit <- x$backend_fit

  # which = 10 sets a 2x2 layout for gam.check(); restore the caller's par.
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  unavailable <- intersect(which, c(1L, 2L))
  if (length(unavailable))
    warning("plot which = ", paste(unavailable, collapse = ", "),
            " is not available for a GAM: the fitted rates and the scaled ",
            "Schoenfeld residuals are per-coefficient quantities that a ",
            "penalised smooth has no analogue for. Use which = 9 (partial ",
            "effects) or which = 10 (mgcv residual / k check) instead.",
            call. = FALSE)

  # Compute the diagnostics here rather than in plot.remstimate: that insists on
  # a 'stats' argument when 'diagnostics' is NULL, while diagnostics.remstimate_gam
  # works entirely off the stored stack.
  # 1 and 2 are dropped, not forwarded: they draw nothing and would still cost
  # the diagnostics computation.
  base_which <- setdiff(which[which <= 8L], c(1L, 2L))
  if (length(base_which)) {
    d <- diagnostics %||% diagnostics(x, reh)
    plot.remstimate(x, reh, diagnostics = d, which = base_which,
                    effects = effects, sender_effects = sender_effects,
                    receiver_effects = receiver_effects, ...)
  }

  # (9) one panel per smooth: the estimated f(statistic) on the log-rate scale.
  # seWithMean = TRUE adds the intercept uncertainty to the band, giving roughly
  # correct across-the-function coverage rather than the pointwise band of the
  # centred smooth alone.
  # A shape-constrained fit is the exception: its covariance is conditional on
  # the active constraints and mgcv declines seWithMean there ("seWithMean
  # unavailable"), so it is only asked for when it can be delivered.
  # A time-varying fit mixes two kinds of panel - the centred baseline s(t) is a
  # deviation on the log-rate, s(t):X is the coefficient of X - so the panels are
  # drawn one at a time and labelled individually rather than sharing one ylab.
  if (9L %in% which) {
    if (is.null(fit) || !length(fit$smooth)) {
      warning("no smooth terms in this fit, so there are no partial-effect ",
              "curves to plot (which = 9).", call. = FALSE)
    } else {
      se_with_mean <- !identical(attr(x, "engine"), "scasm")
      if (!.gam_is_tve(x)) {
        mgcv::plot.gam(fit, pages = 1, shade = TRUE,
                       seWithMean = se_with_mean,
                       rug = FALSE, ylab = "Partial effect on log-rate")
      } else {
        n_sm <- length(fit$smooth)
        nc   <- ceiling(sqrt(n_sm))
        graphics::par(mfrow = c(ceiling(n_sm / nc), nc))
        for (i in seq_len(n_sm)) {
          smi <- fit$smooth[[i]]
          yl  <- if (is.null(smi$by) || identical(smi$by, "NA"))
                   "Time-varying baseline (log-rate)"
                 else paste0("Coefficient of ", smi$by, " over time")
          mgcv::plot.gam(fit, select = i, shade = TRUE,
                         seWithMean = se_with_mean, rug = FALSE, ylab = yl)
        }
      }
    }
  }

  # (10) gam.check draws mgcv's residual panels and prints the k.check table it
  # runs internally: a k-index below 1 with a small p-value means the residuals
  # still hold pattern with respect to the covariate, i.e. k was set too low.
  # Its defaults cap the cost (the check subsamples at 5000 rows), so they are
  # left alone here.
  if (10L %in% which) {
    if (is.null(fit)) {
      warning("no backend fit stored, so which = 10 cannot be drawn.",
              call. = FALSE)
    } else if (isTRUE(attr(x, "ordinal"))) {
      # gam.check()'s residual panels need one residual per fitted value, which
      # a cox.ph fit does not provide (its response is a two-column
      # time/stratum matrix), so it fails there with "'x' and 'y' lengths
      # differ". The basis-dimension table is the part that transfers, and
      # k.check() produces it without the plots.
      message("Ordinal GAM: mgcv's residual panels are not defined for the ",
              "'cox.ph' family; showing the basis-dimension check only.")
      print(mgcv::k.check(fit))
    } else {
      graphics::par(mfrow = c(2, 2))
      mgcv::gam.check(fit)
    }
  }

  invisible(x)
}
