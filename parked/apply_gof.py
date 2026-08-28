import io, os, sys, shutil

ROOT = os.path.expanduser("~/mnt/remstimate")

def read(p):
    with io.open(os.path.join(ROOT, p), encoding="utf-8") as f:
        return f.read()

def write(p, s):
    # write to a temp file and rename, so a hardlinked original is not
    # truncated through its other link
    full = os.path.join(ROOT, p)
    tmp  = full + ".tmp_gof"
    with io.open(tmp, "w", encoding="utf-8") as f:
        f.write(s)
    shutil.copymode(full, tmp)
    os.rename(tmp, full)

def sub1(s, old, new, tag):
    n = s.count(old)
    if n != 1:
        sys.exit("ABORT [%s]: anchor found %d times, expected 1" % (tag, n))
    return s.replace(old, new)

# ─────────────────────────────────────────────── R/diagnostics.R
d = read("R/diagnostics.R")

# 1. signature + roxygen
d = sub1(d,
"#' @param ... further arguments to be passed to the 'diagnostics' method.",
"""#' @param gof_R integer; number of replicates for the simulation-based
#'   goodness-of-fit check (default \\code{0L}, which skips it). When
#'   \\code{gof_R > 0}, replicate outcomes are drawn from the fitted per-event
#'   probabilities -- holding the endogenous statistics fixed at their observed
#'   values -- and aggregated into reference distributions for actor and dyad
#'   frequencies, returned in \\code{$gof} and shown by plots 10-12 of
#'   \\code{plot.diagnostics}. See \\code{.gof_compute} in \\code{R/gof.R} for
#'   what the check does and does not test.
#' @param gof_top_dyads integer; how many dyads to retain in the GOF dyad
#'   table, ranked by observed plus expected count (default \\code{30L}).
#' @param ... further arguments to be passed to the 'diagnostics' method.""",
"roxygen params")

d = sub1(d,
"diagnostics.remstimate <- function(object, reh, stats, top_pct = 0.05, surprise_threshold = 0.2, ...) {",
"diagnostics.remstimate <- function(object, reh, stats, top_pct = 0.05, surprise_threshold = 0.2,\n"
"                                   gof_R = 0L, gof_top_dyads = 30L, ...) {",
"signature")

# 2. tie branch
d = sub1(d,
"""    diagnostics$recall <- .recall_block_3d(
      pars     = as.vector(object$coefficients)[select_vars],
      baseline = baseline_value,
      stats_3d = stats,
      obs_ids  = obs_dyad_ids,
      top_pct  = top_pct
    )""",
"""    diagnostics$recall <- .recall_block_3d(
      pars     = as.vector(object$coefficients)[select_vars],
      baseline = baseline_value,
      stats_3d = stats,
      obs_ids  = obs_dyad_ids,
      top_pct  = top_pct
    )
    # simulation-based GOF (opt-in): reuses the same pars/stats/riskset
    if (gof_R > 0L) {
      diagnostics$gof <- .gof_tie(
        pars         = as.vector(object$coefficients)[select_vars],
        baseline     = baseline_value,
        stats_3d     = stats,
        obs_dyad_ids = obs_dyad_ids,
        reh          = reh,
        R            = gof_R,
        top_dyads    = gof_top_dyads
      )
    }""",
"tie gof")

# 3. actor branch
d = sub1(d,
"""      diagnostics[[which_model[i]]]$recall <- .recall_block_3d(
        pars      = as.vector(object[[which_model[i]]]$coefficients)[select_vars],
        baseline  = baseline_value,
        stats_3d  = stats[[which_stats[i]]],
        obs_ids   = obs_ids_rc,
        valid_ids = valid_ids_rc,
        top_pct   = top_pct
      )""",
"""      diagnostics[[which_model[i]]]$recall <- .recall_block_3d(
        pars      = as.vector(object[[which_model[i]]]$coefficients)[select_vars],
        baseline  = baseline_value,
        stats_3d  = stats[[which_stats[i]]],
        obs_ids   = obs_ids_rc,
        valid_ids = valid_ids_rc,
        top_pct   = top_pct
      )
      # simulation-based GOF (opt-in). The sender is drawn from the rate model;
      # the receiver is drawn conditional on the OBSERVED sender, because
      # receiver_stats at event m are built for the sender that actually acted.
      if (gof_R > 0L) {
        if (is.null(diagnostics$gof)) diagnostics$gof <- list(R = gof_R)
        if (senderRate[i]) {
          diagnostics$gof$sender <- .gof_actor_sender(
            pars      = as.vector(object[[which_model[i]]]$coefficients)[select_vars],
            baseline  = baseline_value,
            stats_3d  = stats[[which_stats[i]]],
            obs_ids   = obs_ids_rc,
            valid_ids = valid_ids_rc,
            reh       = reh,
            R         = gof_R
          )
        } else {
          .rec <- .gof_actor_receiver(
            pars      = as.vector(object[[which_model[i]]]$coefficients)[select_vars],
            stats_3d  = stats[[which_stats[i]]],
            obs_ids   = obs_ids_rc,
            valid_ids = valid_ids_rc,
            reh       = reh,
            R         = gof_R,
            top_dyads = gof_top_dyads
          )
          diagnostics$gof$receiver <- .rec$receiver
          diagnostics$gof$dyad     <- .rec$dyad
        }
      }""",
"actor gof")

write("R/diagnostics.R", d)

# ─────────────────────────────────────────────── R/plot.R
p = read("R/plot.R")

p = sub1(p,
"  which    <- rep(FALSE, max(9L, max(selected)))",
"  which    <- rep(FALSE, max(12L, max(selected)))",
"which width")

p = sub1(p,
"""  # (7) random-effect normality Q-Q (GLMM only; model-agnostic).""",
"""  # (10-12) simulation-based GOF. Stored at the top level of the diagnostics
  # object for both models, so it is plotted once, outside the model branches.
  if (any(which[10:12])) {
    if (is.null(x$gof)) {
      warning("no GOF output on this 'diagnostics' object; re-run diagnostics() ",
              "with gof_R > 0 to enable plots 10-12.", call. = FALSE)
    } else {
      if (which[10L] && !is.null(x$gof$sender)) {
        par(mfrow = c(1, 1))
        .gof_panel(x$gof$sender, "GOF: sender frequencies")
      }
      if (which[11L] && !is.null(x$gof$receiver)) {
        par(mfrow = c(1, 1))
        .gof_panel(x$gof$receiver, "GOF: receiver frequencies")
      }
      if (which[12L] && !is.null(x$gof$dyad)) {
        par(mfrow = c(1, 1))
        .gof_panel(x$gof$dyad, "GOF: dyad frequencies")
      }
    }
  }

  # (7) random-effect normality Q-Q (GLMM only; model-agnostic).""",
"gof plots")

# document the new plot numbers
p = sub1(p,
"#' @param which integer vector of plots to produce. Default \\code{1:3} covers\n"
"#'   waiting times, Schoenfeld residuals, and recall. Add \\code{4} or \\code{5}\n"
"#'   for HMC posterior density and trace plots.",
"#' @param which integer vector of plots to produce. Default \\code{1:3} covers\n"
"#'   waiting times, Schoenfeld residuals, and recall. Add \\code{4} or \\code{5}\n"
"#'   for HMC posterior density and trace plots. Plots \\code{10}, \\code{11} and\n"
"#'   \\code{12} show the simulation-based goodness-of-fit check for actor\n"
"#'   frequencies as sender, as receiver, and for dyads; they require a\n"
"#'   \\code{diagnostics} object built with \\code{gof_R > 0}.",
"which roxygen")

write("R/plot.R", p)
print("OK: R/diagnostics.R and R/plot.R patched")
