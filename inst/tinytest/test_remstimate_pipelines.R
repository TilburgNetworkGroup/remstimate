# ══════════════════════════════════════════════════════════════════════════════
# test_remstimate2.R — comprehensive tests for the refactored remstimate API
# ══════════════════════════════════════════════════════════════════════════════
#
# Covers:
#   - New dispatch: approach + random/penalty/mixture
#   - Backward compat: method = "MLE" / "HMC"
#   - Error handling for invalid inputs
#   - Tie model: directed/undirected, ordinal/interval, full/active riskset
#   - Actor model: directed, ordinal/interval
#   - Duration model: directed/undirected end, ordinal/interval, ext=TRUE/FALSE
#   - GLMM: lme4, glmmTMB, coxme (ordinal)
#   - GLMNET: lasso, ridge, elasticnet
#   - SHRINKEM (Bayesian penalty): tie + actor
#   - MIXREM / dlcrem
#   - remwindow (moving-window), remtribute (attribute model), actor HMC
#   - Wrappers: remfrailty / remixture / rempenalty (+ deprecated frailty_rem)
#   - Argument API: penalty/mixture lists own their knobs; bayes=list() controls
#
# Uses tinytest. Run: tinytest::run_test_file("test_remstimate2.R")
# ══════════════════════════════════════════════════════════════════════════════

library(tinytest)

# ── Test data: small edgelists ────────────────────────────────────────────────

# Directed tie edgelist (20 actors, ~200 events for estimability) with frailty
set.seed(42)
n_actors <- 4
actors   <- paste0("A", seq_len(n_actors))
n_events <- 200
el_tie <- data.frame(
  time   = sort(cumsum(rexp(n_events, 0.5))),
  actor1 = sample(actors, n_events, replace = TRUE,prob=seq_len(n_actors)),
  actor2 = sample(actors, n_events, replace = TRUE,prob=seq_len(n_actors)),
  stringsAsFactors = FALSE
)
el_tie <- el_tie[el_tie$actor1 != el_tie$actor2, ]

# Undirected edgelist
el_undir <- el_tie

# Typed edgelist (2 types)
el_typed <- el_tie
el_typed$setting <- sample(c("X", "Y"), nrow(el_typed), replace = TRUE)

# Duration edgelist with types
el_dur_typed <- el_typed
el_dur_typed$duration <- rexp(nrow(el_dur_typed),1)
el_dur <- el_dur_typed

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Dispatch and error handling
# ══════════════════════════════════════════════════════════════════════════════

# ── 1.1 Invalid structure combinations ────────────────────────────────────────
reh_t <- remify(el_tie, model = "tie", actors = actors)
stats_t <- remstats(reh_t, tie_effects = ~ inertia() + reciprocity())

expect_error(
  remstimate(reh_t, stats_t, random = ~ (1 | actor1), mixture = list(k = 2)),
  pattern = "cannot be combined",
  info = "'mixture' cannot be combined with 'random' or 'penalty'."
)

# ── 1.3 Backward compat: method = "MLE" still works ──────────────────────────
fit_compat <- remstimate(reh_t, stats_t, method = "MLE")
expect_inherits(fit_compat, "remstimate",
                info = "method='MLE' backward compat works")
expect_true(all(is.finite(coef(fit_compat))),
            info = "method='MLE' produces finite coefficients")

# ── 1.4 Default approach = frequentist ────────────────────────────────────────
fit_default <- remstimate(reh_t, stats_t)
expect_inherits(fit_default, "remstimate",
                info = "default call (no approach) works")
expect_equal(attr(fit_default, "approach"), "Frequentist",
             info = "default approach is Frequentist")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Basic MLE — tie model variations
# ══════════════════════════════════════════════════════════════════════════════

# ── 2.1 Directed, interval, full riskset ──────────────────────────────────────
fit_dir <- remstimate(reh_t, stats_t)
expect_true(length(coef(fit_dir)) == 3L,
            info = "directed tie MLE: 3 coefficients")
expect_true(!is.null(fit_dir$loglik),
            info = "directed tie MLE: loglik present")

# ── 2.2 Directed, ordinal, full riskset ──────────────────────────────────────
reh_ord <- remify(el_tie, model = "tie", actors = actors, ordinal = TRUE)
stats_ord <- remstats(reh_ord, tie_effects = ~ inertia() + reciprocity())
fit_ord <- remstimate(reh_ord, stats_ord)
expect_inherits(fit_ord, "remstimate",
                info = "ordinal tie MLE returns remstimate")

# ── 2.3 Undirected, interval ─────────────────────────────────────────────────
reh_und <- remify(el_undir, model = "tie", actors = actors, directed = FALSE)
stats_und <- remstats(reh_und, tie_effects = ~ inertia())
fit_und <- remstimate(reh_und, stats_und)
expect_true(all(is.finite(coef(fit_und))),
            info = "undirected tie MLE: finite coefficients")

# ── 2.4 Active riskset ───────────────────────────────────────────────────────
reh_act <- remify(el_typed, model = "tie", actors = actors, riskset = "active",
                  event_type = "setting")
stats_act <- remstats(reh_act, tie_effects = ~ inertia(consider_type = FALSE))
fit_act <- remstimate(reh_act, stats_act)
expect_inherits(fit_act, "remstimate",
                info = "active riskset tie MLE works")

# ── 2.5 ext = TRUE (typed riskset) ──────────────────────────────────────────
reh_ext <- remify(el_typed, model = "tie", actors = actors,
                  riskset = "active", extend_riskset_by_type = TRUE,
                  event_type = "setting")
stats_ext <- remstats(reh_ext, tie_effects = ~ inertia(consider_type = "interact"))
fit_ext <- remstimate(reh_ext, stats_ext)
expect_inherits(fit_ext, "remstimate",
                info = "ext=TRUE tie MLE works")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Basic MLE — actor model
# ══════════════════════════════════════════════════════════════════════════════

reh_ao <- remify(el_tie, model = "actor", actors = actors)
stats_ao <- remstats(reh_ao,
                     sender_effects   = ~ 1 + indegreeSender() + outdegreeSender(),
                     receiver_effects = ~ inertia() + reciprocity())

# ── 3.1 Actor MLE ────────────────────────────────────────────────────────────
fit_ao <- remstimate(reh_ao, stats_ao)
expect_inherits(fit_ao, "remstimate",
                info = "actor MLE returns remstimate")
expect_true(!is.null(fit_ao$sender_model),
            info = "actor MLE has sender_model")
expect_true(!is.null(fit_ao$receiver_model),
            info = "actor MLE has receiver_model")

# ── 3.2 Actor ordinal MLE ────────────────────────────────────────────────────
reh_ao_ord <- remify(el_tie, model = "actor", actors = actors, ordinal = TRUE,
                     riskset = "active")
stats_ao_ord <- remstats(reh_ao_ord,
                         sender_effects   = ~ 1 + indegreeSender(),
                         receiver_effects = ~ inertia())
fit_ao_ord <- remstimate(reh_ao_ord, stats_ao_ord)
expect_inherits(fit_ao_ord, "remstimate",
                info = "actor ordinal MLE works")

# ── 3.3 Cross-model mismatch errors ──────────────────────────────────────────
expect_error(
  remstimate(reh_t, stats_ao),
  info = "tie reh + actor stats errors"
)
expect_error(
  remstimate(reh_ao, stats_t),
  info = "actor reh + tie stats errors"
)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Basic HMC (backward compat)
# ══════════════════════════════════════════════════════════════════════════════

# ── 4.1 Tie HMC via method= ──────────────────────────────────────────────────
fit_hmc <- remstimate(reh_t, stats_t, method = "HMC",
                      nsim = 20, burnin = 5, thin = 1, nchains = 1, L = 5)
expect_inherits(fit_hmc, "remstimate",
                info = "HMC backward compat works")
expect_equal(attr(fit_hmc, "approach"), "Bayesian",
             info = "HMC sets approach to Bayesian")
expect_true(!is.null(fit_hmc$draws),
            info = "HMC produces draws")

# ── 4.2 Tie HMC via approach= ────────────────────────────────────────────────
fit_hmc2 <- remstimate(reh_t, stats_t, approach = "Bayesian",
                       nsim = 20, burnin = 5, thin = 1, nchains = 1, L = 5)
expect_inherits(fit_hmc2, "remstimate",
                info = "approach='bayesian' C++ HMC works")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: Duration model (durem)
# ══════════════════════════════════════════════════════════════════════════════

# ── 5.1 Basic durem MLE (start only) ─────────────────────────────────────────
reh_dur <- remify(el_dur_typed, duration = TRUE, model = "tie", ordinal = TRUE)
stats_dur <- remstats(reh_dur, start_effects = ~ inertia(), psi_start = 1)
fit_dur <- remstimate(reh_dur, stats_dur)
expect_inherits(fit_dur, "remstimate_durem",
                info = "durem MLE returns remstimate_durem")
expect_true(all(is.finite(coef(fit_dur))),
            info = "durem MLE: finite coefficients")

# ── 5.2 Durem with start + end effects ──────────────────────────────────────
stats_dur_both <- remstats(reh_dur,
                           start_effects = ~ itp(),
                           end_effects   = ~ inertia(),
                           psi_start = 1, psi_end = 1)
fit_dur_both <- remstimate(reh_dur, stats_dur_both)
expect_inherits(fit_dur_both, "remstimate_durem",
                info = "durem start+end MLE works")

# ── 5.3 Durem ordinal ────────────────────────────────────────────────────────
reh_dur_ord <- remify(el_dur, duration = TRUE, model = "tie", ordinal = TRUE)
stats_dur_ord <- remstats(reh_dur_ord, start_effects = ~ inertia(), psi_start = 1)
fit_dur_ord <- remstimate(reh_dur_ord, stats_dur_ord)
expect_equal(attr(fit_dur_ord, "engine"), "clogit",
             info = "ordinal durem uses clogit engine")

# ── 5.4 Durem with dur_directed_end = TRUE ──────────────────────────────────────
reh_dur_de <- remify(el_dur, duration = TRUE, model = "tie", dur_directed_end = TRUE)
stats_dur_de <- remstats(reh_dur_de, start_effects = ~ inertia(), psi_start = 1)
fit_dur_de <- remstimate(reh_dur_de, stats_dur_de)
expect_inherits(fit_dur_de, "remstimate_durem",
                info = "dur_directed_end=TRUE durem works")

# # ── 5.5 Durem Bayesian (brms) ────────────────────────────────────────────────
# if (requireNamespace("brms", quietly = TRUE)) {
#   fit_dur_bayes <- remstimate(reh_dur, stats_dur, approach = "Bayesian",
#                               nsim = 50, burnin = 10, nchains = 1)
#   expect_inherits(fit_dur_bayes, "remstimate_durem",
#                   info = "durem bayesian (brms) works")
# }

# ── 5.6 Durem: start_effects/end_effects rejected for non-durem ──────────────
expect_error(
  remstats(reh_t, start_effects = ~ inertia()),
  info = "start_effects rejected for non-durem reh"
)

# ── 5.7 Durem with ext=TRUE ──────────────────────────────────────────────────
reh_dur_ext <- remify(el_dur_typed, duration = TRUE, model = "tie",
                      riskset = "active", extend_riskset_by_type = TRUE, event_type = "setting")
stats_dur_ext <- remstats(reh_dur_ext,
                          start_effects = ~ inertia(consider_type = "separate"),
                          end_effects = ~ inertia(),
                          psi_start = 1, first = 1)
fit_dur_ext <- remstimate(reh_dur_ext, stats_dur_ext)
expect_inherits(fit_dur_ext, "remstimate_durem",
                info = "durem ext=TRUE works")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: GLMM
# ══════════════════════════════════════════════════════════════════════════════

if (requireNamespace("lme4", quietly = TRUE)) {

  # ── 6.1 Tie model, lme4, interval ──────────────────────────────────────────
  fit_glmm <- remstimate(reh_t, stats_t,
                         random = ~ (1 | actor1) + (1 | actor2))
  expect_true(!is.null(fit_glmm$backend_fit),
              info = "GLMM lme4 tie: backend_fit present")

  # ── 6.2 Tie model, random intercept + slope ────────────────────────────────
  fit_glmm_slope <- remstimate(reh_t, stats_t,
                               random = ~ (1 + inertia | actor1))
  expect_true(!is.null(coef(fit_glmm_slope)),
              info = "GLMM with random slope works")

  # ── 6.3 Actor model GLMM (receiver choice model always uses coxme) ─────────
  if (requireNamespace("coxme", quietly = TRUE)) {
    fit_glmm_ao <- remstimate(reh_ao, stats_ao,
                              random = ~ (1 | actor_label))
    expect_true(!is.null(fit_glmm_ao),
                info = "GLMM actor model works")
  }

  # ── 6.4 Durem GLMM ────────────────────────────────────────────────────────
  fit_glmm_dur <- remstimate(reh_dur, stats_dur,
                             random = ~ (1 | actor1) + (1 | actor2))
  expect_true(!is.null(fit_glmm_dur),
              info = "GLMM durem works")
}

if (requireNamespace("glmmTMB", quietly = TRUE)) {
  # ── 6.5 glmmTMB engine ─────────────────────────────────────────────────────
  fit_tmb <- remstimate(reh_t, stats_t,
                        random = ~ (1 | actor1),
                        engine = "glmmTMB")
  expect_true(!is.null(fit_tmb),
              info = "GLMM glmmTMB engine works")
}

if (requireNamespace("coxme", quietly = TRUE)) {
  # ── 6.6 Ordinal GLMM (auto-selects coxme) ──────────────────────────────────
  fit_coxme <- remstimate(reh_ord, stats_ord,
                          random = ~ (1 | actor1))
  expect_true(!is.null(fit_coxme),
              info = "ordinal GLMM auto-selects coxme")
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6b: revised GLMM behaviour
#   default engine (glmmTMB), baseline intercept, undirected-frailty guard +
#   dyad fallback, receiver conditional logit, convergence flag
# ══════════════════════════════════════════════════════════════════════════════

if (requireNamespace("glmmTMB", quietly = TRUE)) {

  # ── 6b.1 Default engine resolves to glmmTMB ────────────────────────────────
  fit_def_eng <- remstimate(reh_t, stats_t, random = ~ (1 | actor1))
  expect_equal(attr(fit_def_eng, "engine"), "glmmTMB",
               info = "default GLMM engine is glmmTMB")

  # ── 6b.2 baseline is the intercept: kept (finite), no phantom (Intercept) ───
  cf_def <- coef(fit_def_eng)
  expect_true("baseline" %in% names(cf_def) && all(is.finite(cf_def)),
              info = "baseline kept as a finite intercept (not dropped to NA)")
  expect_false("(Intercept)" %in% names(cf_def),
               info = "no duplicate implicit intercept in the GLMM formula")

  # ── 6b.3 Healthy glmmTMB fit is flagged as converged ───────────────────────
  expect_true(isTRUE(attr(fit_def_eng$backend_fit, "remstimate_converged")),
              info = "healthy glmmTMB fit flagged as converged")

  # ── 6b.4 Undirected: actor-level frailty is rejected ───────────────────────
  expect_error(
    remstimate(reh_und, stats_und, random = ~ (1 | actor1) + (1 | actor2)),
    pattern = "not identified",
    info = "undirected actor-level frailty errors"
  )

  # ── 6b.5 Undirected: dyad-level frailty is allowed ─────────────────────────
  fit_und_dyad <- remstimate(reh_und, stats_und, random = ~ (1 | dyad))
  expect_true(!is.null(fit_und_dyad$backend_fit),
              info = "undirected dyad frailty works")

  # ── 6b.6 remfrailty on undirected falls back to dyad-level (with message) ───
  fb_msgs <- character(0)
  withCallingHandlers(
    remfrailty(reh_und, stats_und),
    message = function(m) {
      fb_msgs <<- c(fb_msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("dyad", fb_msgs)),
              info = "remfrailty undirected -> dyad-level fallback (message)")
}

if (requireNamespace("glmmTMB", quietly = TRUE) &&
    requireNamespace("coxme",   quietly = TRUE)) {

  # ── 6b.7 Actor model: receiver choice model is fitted with coxme ───────────
  fit_ao_re <- remstimate(reh_ao, stats_ao,
                          random = list(sender   = ~ (1 | actor),
                                        receiver = ~ (1 | actor)))
  expect_equal(fit_ao_re$engine_receiver, "coxme",
               info = "actor receiver choice model uses coxme")
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: GLMNET (penalised)
# ══════════════════════════════════════════════════════════════════════════════

if (requireNamespace("glmnet", quietly = TRUE)) {

  # ── 7.0 Lasso (alpha = 1) ──────────────────────────────────────────────────
  reh_act_ext <- remify(el_typed, model = "tie", actors = actors, riskset = "active",
                    event_type = "setting", extend_riskset_by_type = TRUE)
  stats_act_ext <- remstats(reh_act_ext, tie_effects = ~ inertia(consider_type = FALSE)
                        + FEtype())
  fit_act_ext_MLE <- remstimate(reh_act_ext, stats_act_ext)
  fit_lasso <- remstimate(reh_act_ext, stats_act_ext, penalty = list(alpha = 1))
  expect_true(!is.null(coef(fit_lasso)),
              info = "GLMNET lasso works")

  # ── 7.1 Lasso (alpha = 1) ──────────────────────────────────────────────────
  fit_lasso <- remstimate(reh_t, stats_t, penalty = list(alpha = 1))
  expect_true(!is.null(coef(fit_lasso)),
              info = "GLMNET lasso works")

  # ── 7.2 Ridge (alpha = 0) ──────────────────────────────────────────────────
  fit_ridge <- remstimate(reh_t, stats_t, penalty = list(alpha = 0))
  expect_true(!is.null(coef(fit_ridge)),
              info = "GLMNET ridge works")

  # ── 7.3 Elastic net (alpha = 0.5) ──────────────────────────────────────────
  fit_enet <- remstimate(reh_t, stats_t, penalty = list(alpha = 0.5))
  expect_true(!is.null(coef(fit_enet)),
              info = "GLMNET elastic net works")

  # ── 7.4 GLMNET with actor model ────────────────────────────────────────────
  fit_lasso_ao <- remstimate(reh_ao, stats_ao, penalty = list(alpha = 1))

  # ── 7.5 GLMNET with durem ──────────────────────────────────────────────────
  expect_error(remstimate(reh_dur, stats_dur, penalty = list(alpha = 1)), info="number of effects too low")

  # path that should work: >= 2 statistics so glmnet runs
  stats_dur <- remstats(reh_dur, start_effects = ~ inertia() + reciprocity(scaling = "std"), psi_start = 1)
  fit_lasso_dur <- remstimate(reh_dur, stats_dur, penalty = list(alpha = 1))
  expect_true(!is.null(fit_lasso_dur), info = "GLMNET durem (>=2 stats) fits")
  expect_equal(length(coef(fit_lasso_dur)), 2L,
               info = "two penalized coefficients returned")

  # path that should work: >= 2 statistics so glmnet runs
  stats_dur2 <- remstats(reh_dur_de, start_effects = ~ 1 + inertia() + reciprocity(scaling = "std"),
                        end_effects = ~ 1 + inertia() + psABAB())
  fit_lasso_dur2 <- remstimate(reh_dur_de, stats_dur2, penalty = list(alpha = 1, penalized = "psABAB.end"))
  expect_true(!is.null(fit_lasso_dur2), info = "GLMNET durem (>=2 stats) fits")
  expect_equal(length(coef(fit_lasso_dur)), 2L,
               info = "two penalized coefficients returned")

  # ── 7.6 lambda_select = "min" ──────────────────────────────────────────────
  fit_lmin <- remstimate(reh_t, stats_t,
                         penalty = list(alpha = 1),
                         lambda_select = "min")
  expect_true(!is.null(coef(fit_lmin)),
              info = "GLMNET lambda_select='min' works")
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8: MIXREM / dlcrem
# ══════════════════════════════════════════════════════════════════════════════

if (requireNamespace("flexmix", quietly = TRUE)) {

  # ── 8.1 Basic MIXREM ───────────────────────────────────────────────────────
  fit_mix <- remstimate(reh_t, stats_t,
                        mixture = list(k = 2, random = ~ (1 | dyad)))
  expect_true(!is.null(fit_mix),
              info = "MIXREM k=2 works")

  # ── 8.2 dlcrem wrapper ─────────────────────────────────────────────────────
  fit_dlc <- dlcrem(reh_t, stats_t, k = 2)
  expect_true(!is.null(fit_dlc),
              info = "dlcrem wrapper works")

  # ── 8.3 MIXREM with multiple k ─────────────────────────────────────────────
  fit_mix_multi <- remstimate(reh_t, stats_t,
                              mixture = list(k = 2:3, random = ~ (1 | dyad)))
  expect_true(length(fit_mix_multi) == 2L,
              info = "MIXREM with k=2:3 returns 2 fits")

  # ── 8.4 MIXREM with durem ──────────────────────────────────────────────────
  set.seed(123)
  stats_dur_mixrem <- stats_dur
  stats_dur_mixrem$stacked$remstats_stack <- stats_dur_mixrem$stacked$remstats_stack[1:1000,]
  fit_mix_dur <- remstimate(reh_dur, stats_dur_mixrem,
                            mixture = list(k = 2, random = ~ (1 | dyad)), nrep = 1)
  expect_true(!is.null(fit_mix_dur),
              info = "MIXREM durem works")
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9: frailty_rem wrapper
# ══════════════════════════════════════════════════════════════════════════════

if (requireNamespace("lme4", quietly = TRUE)) {

  # ── 9.1 Tie model frailty ──────────────────────────────────────────────────
  fit_fr <- remfrailty(reh_t, stats_t)
  expect_true(!is.null(fit_fr),
              info = "frailty_rem tie model works")

  # ── 9.2 Actor model frailty ────────────────────────────────────────────────
  fit_fr_ao <- remfrailty(reh_ao, stats_ao)
  expect_true(!is.null(fit_fr_ao),
              info = "frailty_rem actor model works")

  # ── 9.3 Durem frailty ──────────────────────────────────────────────────────
  fit_fr_dur <- remfrailty(reh_dur, stats_dur)
  expect_true(!is.null(fit_fr_dur),
              info = "frailty_rem durem works")

  # ── 9.4 Typed tie frailty (ext=TRUE) ───────────────────────────────────────
  fit_fr_ext <- remfrailty(reh_ext, stats_ext)
  expect_true(!is.null(fit_fr_ext),
              info = "frailty_rem ext=TRUE works")

  # ── 9.5 Typed durem frailty ────────────────────────────────────────────────
  fit_fr_dur_ext <- remfrailty(reh_dur_ext, stats_dur_ext)
  expect_true(!is.null(fit_fr_dur_ext),
              info = "frailty_rem durem ext=TRUE works")
}

# # ── 9.6 Bayesian frailty ─────────────────────────────────────────────────────
# if (requireNamespace("brms", quietly = TRUE)) {
#   fit_fr_bayes <- remfrailty(reh_t, stats_t, approach = "bayesian",
#                               nsim = 50, burnin = 10, nchains = 1)
#   expect_true(!is.null(fit_fr_bayes),
#               info = "frailty_rem bayesian (brms) works")
# }

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 10: Bayesian random effects (brms)
# ══════════════════════════════════════════════════════════════════════════════

# if (requireNamespace("brms", quietly = TRUE)) {
#
#   # ── 10.1 Tie model, random intercept ────────────────────────────────────────
#   fit_brms <- remstimate(reh_t, stats_t, approach = "bayesian",
#                          random = ~ (1 | actor1),
#                          nsim = 50, burnin = 10, nchains = 1)
#   expect_true(!is.null(fit_brms),
#               info = "brms random intercept works")
#
#   # ── 10.2 Tie model, random slope ───────────────────────────────────────────
#   fit_brms_slope <- remstimate(reh_t, stats_t, approach = "bayesian",
#                                random = ~ (1 + inertia | actor1),
#                                nsim = 50, burnin = 10, nchains = 1)
#   expect_true(!is.null(fit_brms_slope),
#               info = "brms random slope works")
#
#   # ── 10.3 Durem brms ────────────────────────────────────────────────────────
#   fit_brms_dur <- remstimate(reh_dur, stats_dur, approach = "bayesian",
#                              random = ~ (1 | actor1),
#                              nsim = 50, burnin = 10, nchains = 1)
#   expect_true(!is.null(fit_brms_dur),
#               info = "brms durem works")
# }

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 11: Diagnostics and logLik
# ══════════════════════════════════════════════════════════════════════════════

# ── 11.1 Tie MLE diagnostics ─────────────────────────────────────────────────
diag_tie <- diagnostics(fit_default, reh_t, stats_t)
expect_true(!is.null(diag_tie$residuals),
            info = "tie MLE diagnostics: residuals present")
expect_true(!is.null(diag_tie$recall),
            info = "tie MLE diagnostics: recall present")

# ── 11.2 Actor MLE diagnostics ───────────────────────────────────────────────
diag_ao <- diagnostics(fit_ao, reh_ao, stats_ao)
expect_true(!is.null(diag_ao$sender_model),
            info = "actor MLE diagnostics: sender_model present")
expect_true(!is.null(diag_ao$receiver_model),
            info = "actor MLE diagnostics: receiver_model present")

# ── 11.3 Durem diagnostics ───────────────────────────────────────────────────
diag_dur <- diagnostics(fit_dur, reh_dur, stats_dur)
expect_inherits(diag_dur, "diagnostics_durem",
                info = "durem diagnostics class correct")
expect_true(!is.null(diag_dur$recall_joint),
            info = "durem diagnostics: recall_joint present")
expect_true(!is.null(diag_dur$recall_start),
            info = "durem diagnostics: recall_start present")

# ── 11.4 logLik ──────────────────────────────────────────────────────────────
ll_tie <- logLik(fit_default)
expect_true(is.finite(ll_tie), info = "tie MLE logLik is finite")

ll_dur <- logLik(fit_dur)
expect_inherits(ll_dur, "logLik", info = "durem logLik has correct class")
expect_true(is.finite(ll_dur), info = "durem logLik is finite")

# ── 11.5 Durem ext=TRUE diagnostics (per-type recall) ────────────────────────
diag_dur_ext <- diagnostics(fit_dur_ext, reh_dur_ext, stats_dur_ext)
expect_true(!is.null(diag_dur_ext$recall_by_type),
            info = "durem ext=TRUE: recall_by_type present")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 12: Stack stats validation
# ══════════════════════════════════════════════════════════════════════════════

# ── 12.1 Durem stacked output structure ──────────────────────────────────────
stacked <- stack_stats(stats_dur, reh_dur)
expect_true("obs" %in% names(stacked$remstats_stack),
            info = "stacked durem has 'obs' column")
expect_true("time_index" %in% names(stacked$remstats_stack),
            info = "stacked durem has 'time_index' column")
expect_true("dyad" %in% names(stacked$remstats_stack),
            info = "stacked durem has 'dyad' column")
expect_true("process" %in% names(stacked$remstats_stack),
            info = "stacked durem has 'process' column")
expect_true(all(stacked$remstats_stack$process %in% c("start", "end")),
            info = "process column has only 'start'/'end' values")

# ── 12.2 Ordinal durem has no log_interevent ─────────────────────────────────
stacked_ord <- stack_stats(stats_dur_ord, reh_dur_ord)
expect_true(!"log_interevent" %in% names(stacked_ord$remstats_stack),
            info = "ordinal durem stack has no log_interevent")

# ── 12.3 Type column present when ext=TRUE ───────────────────────────────────
stacked_ext <- stack_stats(stats_dur_ext, reh_dur_ext, add_actors = FALSE)
expect_true("type" %in% names(stacked_ext$remstats_stack),
            info = "ext=TRUE durem stack has 'type' even with add_actors=FALSE")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 13: Model selection helpers
# ══════════════════════════════════════════════════════════════════════════════

# ── 13.1 AIC / BIC for tie MLE ───────────────────────────────────────────────
expect_true(is.finite(fit_default$AIC), info = "tie MLE AIC is finite")
expect_true(is.finite(fit_default$BIC), info = "tie MLE BIC is finite")

# ── 13.2 AIC / BIC for durem ─────────────────────────────────────────────────
expect_true(is.finite(fit_dur$AIC), info = "durem AIC is finite")
expect_true(is.finite(fit_dur$BIC), info = "durem BIC is finite")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 14: Edge cases
# ══════════════════════════════════════════════════════════════════════════════

# ── 14.1 WAIC validation ─────────────────────────────────────────────────────
expect_error(
  remstimate(reh_t, stats_t, WAIC = "yes"),
  pattern = "logical",
  info = "WAIC='yes' errors"
)

# ── 14.2 ncores validation ───────────────────────────────────────────────────
expect_warning(
  remstimate(reh_t, stats_t, ncores = -1),
  pattern = "positive",
  info = "negative ncores warns"
)

# ── 14.3 nsimWAIC validation ─────────────────────────────────────────────────
expect_warning(
  remstimate(reh_t, stats_t, nsimWAIC = -5),
  pattern = "positive",
  info = "negative nsimWAIC warns"
)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 15: Diagnostics for all pipelines
# ══════════════════════════════════════════════════════════════════════════════

# ── 15.1 Tie MLE diagnostics (directed, interval) ────────────────────────────
diag_tie_mle <- diagnostics(fit_default, reh_t, stats_t)
expect_inherits(diag_tie_mle, "diagnostics",
                info = "tie MLE diagnostics class correct")
expect_true(!is.null(diag_tie_mle$residuals),
            info = "tie MLE diagnostics: residuals present")
expect_true(!is.null(diag_tie_mle$recall),
            info = "tie MLE diagnostics: recall present")
expect_true(!is.null(diag_tie_mle$rates),
            info = "tie MLE diagnostics: rates present")

# ── 15.2 Tie MLE diagnostics (ordinal) ───────────────────────────────────────
diag_tie_ord <- diagnostics(fit_ord, reh_ord, stats_ord)
expect_inherits(diag_tie_ord, "diagnostics",
                info = "ordinal tie MLE diagnostics works")
expect_true(!is.null(diag_tie_ord$recall),
            info = "ordinal tie diagnostics: recall present")

# ── 15.3 Tie HMC diagnostics ─────────────────────────────────────────────────
diag_hmc <- diagnostics(fit_hmc, reh_t, stats_t)
expect_inherits(diag_hmc, "diagnostics",
                info = "HMC diagnostics class correct")
expect_true(!is.null(diag_hmc$recall),
            info = "HMC diagnostics: recall present")

# ── 15.4 Actor MLE diagnostics ───────────────────────────────────────────────
diag_ao_mle <- diagnostics(fit_ao, reh_ao, stats_ao)
expect_inherits(diag_ao_mle, "diagnostics",
                info = "actor MLE diagnostics class correct")
expect_true(!is.null(diag_ao_mle$sender_model$residuals),
            info = "actor MLE diagnostics: sender residuals present")
expect_true(!is.null(diag_ao_mle$sender_model$recall),
            info = "actor MLE diagnostics: sender recall present")
expect_true(!is.null(diag_ao_mle$receiver_model$recall),
            info = "actor MLE diagnostics: receiver recall present")

# ── 15.5 Actor ordinal diagnostics ───────────────────────────────────────────
diag_ao_ord <- diagnostics(fit_ao_ord, reh_ao_ord, stats_ao_ord)
expect_inherits(diag_ao_ord, "diagnostics",
                info = "actor ordinal diagnostics works")

# ── 15.6 Durem MLE diagnostics (start only) ──────────────────────────────────
diag_dur_mle <- diagnostics(fit_dur, reh_dur, stats_dur)
expect_inherits(diag_dur_mle, "diagnostics_durem",
                info = "durem MLE diagnostics class correct")
expect_true(!is.null(diag_dur_mle$recall_joint),
            info = "durem diagnostics: recall_joint present")
expect_true(!is.null(diag_dur_mle$recall_start),
            info = "durem diagnostics: recall_start present")
expect_true(!is.null(diag_dur_mle$deviance_residuals),
            info = "durem diagnostics: deviance_residuals present")

# ── 15.7 Durem with start + end diagnostics ──────────────────────────────────
diag_dur_both <- diagnostics(fit_dur_both, reh_dur, stats_dur_both)
expect_true(!is.null(diag_dur_both$recall_joint),
            info = "durem start+end diagnostics: recall_joint present")
expect_true(!is.null(diag_dur_both$recall_start),
            info = "durem start+end diagnostics: recall_start present")
expect_true(!is.null(diag_dur_both$recall_end),
            info = "durem start+end diagnostics: recall_end present")

# ── 15.8 Durem ordinal diagnostics ───────────────────────────────────────────
diag_dur_ord <- diagnostics(fit_dur_ord, reh_dur_ord, stats_dur_ord)
expect_inherits(diag_dur_ord, "diagnostics_durem",
                info = "durem ordinal diagnostics works")

# ── 15.9 Durem ext=TRUE diagnostics (per-type recall) ────────────────────────
diag_dur_ext <- diagnostics(fit_dur_ext, reh_dur_ext, stats_dur_ext)
expect_true(!is.null(diag_dur_ext$recall_by_type),
            info = "durem ext=TRUE: recall_by_type present")

# ── 15.10 GLMM diagnostics ───────────────────────────────────────────────────
if (requireNamespace("lme4", quietly = TRUE)) {
  diag_glmm_tie <- diagnostics(fit_glmm, reh_t, stats_t)
  expect_true(!is.null(diag_glmm_tie),
              info = "GLMM tie diagnostics runs")
  expect_true(!is.null(diag_glmm_tie$recall),
              info = "GLMM tie diagnostics: recall present")

  # GLMM actor diagnostics
  diag_glmm_ao <- diagnostics(fit_glmm_ao, reh_ao, stats_ao)
  # expect_true(!is.null(diag_glmm_ao),
  #             info = "GLMM actor diagnostics runs")

  # GLMM durem diagnostics
  diag_glmm_dur <- diagnostics(fit_glmm_dur, reh_dur, stats_dur)
  expect_true(!is.null(diag_glmm_dur),
              info = "GLMM durem diagnostics runs")
}

# ── 15.11 GLMNET diagnostics ─────────────────────────────────────────────────
if (requireNamespace("glmnet", quietly = TRUE)) {
  diag_lasso_tie <- diagnostics(fit_lasso, reh_t, stats_t)
  expect_true(!is.null(diag_lasso_tie),
              info = "GLMNET tie diagnostics runs")
  expect_true(!is.null(diag_lasso_tie$recall),
              info = "GLMNET tie diagnostics: recall present")

  # GLMNET actor diagnostics
  diag_lasso_ao <- diagnostics(fit_lasso_ao, reh_ao, stats_ao)
  expect_true(!is.null(diag_lasso_ao),
              info = "GLMNET actor diagnostics runs")

  # GLMNET durem diagnostics
  diag_lasso_dur <- diagnostics(fit_lasso_dur, reh_dur, stats_dur)
  expect_true(!is.null(diag_lasso_dur),
              info = "GLMNET durem diagnostics runs")
}

# ── 15.12 MIXREM diagnostics ─────────────────────────────────────────────────
if (requireNamespace("flexmix", quietly = TRUE)) {
  diag_mix_tie <- diagnostics(fit_mix, reh_t, stats_t)
  expect_inherits(diag_mix_tie, "diagnostics_mixrem",
                  info = "MIXREM diagnostics class correct")
  expect_true(!is.null(diag_mix_tie$recall_by_component),
              info = "MIXREM diagnostics: per_component present")

  # dlcrem diagnostics
  diag_dlc <- diagnostics(fit_dlc, reh_t, stats_t)
  expect_true(!is.null(diag_dlc),
              info = "dlcrem diagnostics runs")

  # MIXREM durem diagnostics
  diag_mix_dur <- diagnostics(fit_mix_dur, reh_dur, stats_dur_mixrem)
  expect_true(!is.null(diag_mix_dur),
              info = "MIXREM durem diagnostics runs")
}

# ── 15.13 frailty_rem diagnostics ────────────────────────────────────────────
if (requireNamespace("lme4", quietly = TRUE)) {
  diag_fr_tie <- diagnostics(fit_fr, reh_t, stats_t)
  expect_true(!is.null(diag_fr_tie),
              info = "frailty_rem tie diagnostics runs")

  diag_fr_ao <- diagnostics(fit_fr_ao, reh_ao, stats_ao)
  # expect_true(!is.null(diag_fr_ao),
  #             info = "frailty_rem actor diagnostics runs")

  diag_fr_dur <- diagnostics(fit_fr_dur, reh_dur, stats_dur)
  expect_true(!is.null(diag_fr_dur),
              info = "frailty_rem durem diagnostics runs")
}

# # ── 15.14 Bayesian (brms) diagnostics ────────────────────────────────────────
# if (requireNamespace("brms", quietly = TRUE)) {
#   diag_brms <- diagnostics(fit_brms, reh_t, stats_t)
#   expect_true(!is.null(diag_brms),
#               info = "brms diagnostics runs")
# }

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 16: Plot methods for all pipelines
# ══════════════════════════════════════════════════════════════════════════════
#
# All plot tests use pdf(NULL) to suppress graphical output, then dev.off().
# We just verify they don't error.

.test_plot_no_error <- function(expr, info) {
  pdf(NULL)
  err <- tryCatch({ expr; NULL }, error = function(e) e)
  dev.off()
  expect_true(is.null(err), info = info)
}

# ── 16.1 Tie MLE plots ───────────────────────────────────────────────────────
# which=1: waiting-time Q-Q, which=2: Schoenfeld, which=3: recall
.test_plot_no_error(
  plot(fit_default, reh_t, stats = stats_t, which = 1),
  info = "tie MLE plot which=1 (QQ) runs"
)
.test_plot_no_error(
  plot(fit_default, reh_t, stats = stats_t, which = 2),
  info = "tie MLE plot which=2 (Schoenfeld) runs"
)
.test_plot_no_error(
  plot(fit_default, reh_t, stats = stats_t, which = 3),
  info = "tie MLE plot which=3 (recall) runs"
)

# ── 16.2 Tie MLE ordinal plot (no QQ) ────────────────────────────────────────
.test_plot_no_error(
  plot(fit_ord, reh_ord, stats = stats_ord, which = 3),
  info = "ordinal tie MLE plot which=3 (recall) runs"
)

# ── 16.3 Tie HMC plots (posterior + trace) ───────────────────────────────────
.test_plot_no_error(
  plot(fit_hmc, reh_t, stats = stats_t, which = 3),
  info = "HMC plot which=3 (recall) runs"
)
.test_plot_no_error(
  plot(fit_hmc, reh_t, stats = stats_t, which = 4),
  info = "HMC plot which=4 (posterior density) runs"
)
.test_plot_no_error(
  plot(fit_hmc, reh_t, stats = stats_t, which = 5),
  info = "HMC plot which=5 (trace) runs"
)

# ── 16.4 Actor MLE plots ─────────────────────────────────────────────────────
.test_plot_no_error(
  plot(fit_ao, reh_ao, stats = stats_ao, which = 1),
  info = "actor MLE plot which=1 (QQ) runs"
)
.test_plot_no_error(
  plot(fit_ao, reh_ao, stats = stats_ao, which = 2),
  info = "actor MLE plot which=2 (Schoenfeld) runs"
)
.test_plot_no_error(
  plot(fit_ao, reh_ao, stats = stats_ao, which = 3),
  info = "actor MLE plot which=3 (recall) runs"
)

# ── 16.5 Durem plots ─────────────────────────────────────────────────────────
# which=1: joint recall, which=2: deviance residuals,
# which=3: start recall, which=5: all recalls side by side
.test_plot_no_error(
  plot(fit_dur, reh = reh_dur, stats = stats_dur, which = 1L),
  info = "durem plot which=1 (joint recall) runs"
)
.test_plot_no_error(
  plot(fit_dur, reh = reh_dur, stats = stats_dur, which = 2L),
  info = "durem plot which=2 (deviance residuals) runs"
)
.test_plot_no_error(
  plot(fit_dur, reh = reh_dur, stats = stats_dur, which = 3L),
  info = "durem plot which=3 (start recall) runs"
)
.test_plot_no_error(
  plot(fit_dur, reh = reh_dur, stats = stats_dur, which = 5L),
  info = "durem plot which=5 (all recalls) runs"
)

# Durem with start + end: also test end recall (which=4)
.test_plot_no_error(
  plot(fit_dur_both, reh = reh_dur, stats = stats_dur_both, which = 4L),
  info = "durem start+end plot which=4 (end recall) runs"
)

# Durem coefficient CI fallback (no diagnostics)
.test_plot_no_error(
  plot(fit_dur, which = 0L),
  info = "durem coefficient CI fallback plot runs"
)

# Durem prob ratio plots
.test_plot_no_error(
  plot(fit_dur, reh = reh_dur, stats = stats_dur, which = 9L),
  info = "durem plot which=9 (prob ratio joint) runs"
)
.test_plot_no_error(
  plot(fit_dur, reh = reh_dur, stats = stats_dur, which = 10L),
  info = "durem plot which=10 (prob ratio start+end) runs"
)

# Durem ext=TRUE: per-type recall (which=6)
.test_plot_no_error(
  plot(fit_dur_ext, reh = reh_dur_ext, stats = stats_dur_ext, which = 6L),
  info = "durem ext=TRUE plot which=6 (per-type recall) runs"
)

# ── 16.6 GLMM plots ──────────────────────────────────────────────────────────
if (requireNamespace("lme4", quietly = TRUE)) {
  # Tie GLMM: standard diagnostic plots
  .test_plot_no_error(
    plot(fit_glmm, reh_t, stats = stats_t, which = 3),
    info = "GLMM tie plot which=3 (recall) runs"
  )
  # GLMM: random effects Q-Q (which=6)
  .test_plot_no_error(
    plot(fit_glmm, reh_t, stats = stats_t, which = 6),
    info = "GLMM tie plot which=6 (ranef QQ) runs"
  )
}

# ── 16.7 GLMNET plots ────────────────────────────────────────────────────────
if (requireNamespace("glmnet", quietly = TRUE)) {
  # Tie GLMNET: recall
  .test_plot_no_error(
    plot(fit_lasso, reh_t, stats = stats_t, which = 3),
    info = "GLMNET tie plot which=3 (recall) runs"
  )
  # Regularization path (which=6)
  .test_plot_no_error(
    plot(fit_lasso, reh_t, stats = stats_t, which = 6),
    info = "GLMNET tie plot which=6 (regularization path) runs"
  )
}

# ── 16.8 MIXREM plots ────────────────────────────────────────────────────────
if (requireNamespace("flexmix", quietly = TRUE)) {
  # MIXREM recall: per component + combined
  .test_plot_no_error(
    plot(fit_mix, reh_t, stats = stats_t, which = 3L),
    info = "MIXREM tie plot which=3 (recall) runs"
  )

  # dlcrem plot
  .test_plot_no_error(
    plot(fit_dlc, reh_t, stats = stats_t, which = 3L),
    info = "dlcrem plot which=3 (recall) runs"
  )
}

# ── 16.9 frailty_rem plots ───────────────────────────────────────────────────
if (requireNamespace("lme4", quietly = TRUE)) {
  .test_plot_no_error(
    plot(fit_fr, reh_t, stats = stats_t, which = 3),
    info = "frailty_rem tie plot which=3 (recall) runs"
  )
  .test_plot_no_error(
    plot(fit_fr, reh_t, stats = stats_t, which = 6),
    info = "frailty_rem tie plot which=6 (ranef QQ) runs"
  )
}

# ── 16.10 Durem ordinal plot ─────────────────────────────────────────────────
.test_plot_no_error(
  plot(fit_dur_ord, reh = reh_dur_ord, stats = stats_dur_ord, which = 1L),
  info = "durem ordinal plot which=1 (joint recall) runs"
)

# ── 16.11 Undirected tie MLE plot ─────────────────────────────────────────────
.test_plot_no_error(
  plot(fit_und, reh_und, stats = stats_und, which = 3),
  info = "undirected tie MLE plot which=3 (recall) runs"
)

# ── 16.12 Active riskset plot ─────────────────────────────────────────────────
.test_plot_no_error(
  plot(fit_act, reh_act, stats = stats_act, which = 3),
  info = "active riskset tie plot which=3 (recall) runs"
)

# # ── 16.13 Bayesian (brms) plot ────────────────────────────────────────────────
# if (requireNamespace("brms", quietly = TRUE)) {
#   .test_plot_no_error(
#     plot(fit_brms, reh_t, stats = stats_t, which = 3),
#     info = "brms plot which=3 (recall) runs"
#   )
# }

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 17: SHRINKEM (Bayesian penalty) — basic "does it run" checks
# ══════════════════════════════════════════════════════════════════════════════

if (requireNamespace("shrinkem", quietly = TRUE)) {

  # ── 17.1 Tie model, horseshoe ──────────────────────────────────────────────
  fit_hs <- remstimate(reh_t, stats_t,
                       approach = "Bayesian",
                       penalty  = list(prior = "horseshoe"))
  expect_inherits(fit_hs, "remstimate_shrinkem",
                  info = "shrinkem tie returns remstimate_shrinkem")
  expect_true(all(is.finite(coef(fit_hs))),
              info = "shrinkem tie: finite coefficients")

  # ── 17.2 Tie model, lasso prior ────────────────────────────────────────────
  fit_hs_lasso <- remstimate(reh_t, stats_t,
                             approach = "Bayesian",
                             penalty  = list(prior = "lasso"))
  expect_inherits(fit_hs_lasso, "remstimate_shrinkem",
                  info = "shrinkem tie (lasso prior) works")

  # ── 17.3 Actor model, horseshoe ────────────────────────────────────────────
  fit_hs_ao <- remstimate(reh_ao, stats_ao,
                          approach = "Bayesian",
                          penalty  = list(prior = "horseshoe"))
  expect_inherits(fit_hs_ao, "remstimate_shrinkem",
                  info = "shrinkem actor returns remstimate_shrinkem")
  expect_true(!is.null(fit_hs_ao$sender_model) &&
              !is.null(fit_hs_ao$receiver_model),
              info = "shrinkem actor has sender + receiver models")

  # ── 17.4 shrinkem diagnostics (tie + actor) ─────────────────────────────────
  diag_hs <- diagnostics(fit_hs, reh_t, stats_t)
  expect_true(!is.null(diag_hs),
              info = "shrinkem tie diagnostics runs")
  diag_hs_ao <- diagnostics(fit_hs_ao, reh_ao, stats_ao)
  expect_true(!is.null(diag_hs_ao),
              info = "shrinkem actor diagnostics runs")

  # ── 17.5 shrinkem plots (recall + shrinkage scatter) ────────────────────────
  .test_plot_no_error(
    plot(fit_hs, reh_t, stats = stats_t, which = 3),
    info = "shrinkem tie plot which=3 (recall) runs"
  )
  .test_plot_no_error(
    plot(fit_hs, reh_t, which = 6),
    info = "shrinkem tie plot which=6 (shrinkage scatter) runs"
  )
  .test_plot_no_error(
    plot(fit_hs_ao, reh_ao, which = 6),
    info = "shrinkem actor plot which=6 (shrinkage scatter) runs"
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 18: remwindow (moving-window, basic REM) — basic checks
# ══════════════════════════════════════════════════════════════════════════════

# ── 18.1 Tie model, auto mode ────────────────────────────────────────────────
fit_win <- remwindow(reh_t, stats_t, n.windows = 2L)
expect_inherits(fit_win, "remstimate_window",
                info = "remwindow tie returns remstimate_window")
expect_true(length(fit_win$fits) >= 1L,
            info = "remwindow tie: at least one window fit")
expect_true(!is.null(coef(fit_win)),
            info = "remwindow tie: coef() runs")

# ── 18.2 Tie model, manual mode (fixed window width) ─────────────────────────
fit_win_manual <- remwindow(reh_t, stats_t, window.width = 70L)
expect_inherits(fit_win_manual, "remstimate_window",
                info = "remwindow tie manual mode works")

# ── 18.3 Actor model, auto mode ──────────────────────────────────────────────
fit_win_ao <- remwindow(reh_ao, stats_ao, n.windows = 2L)
expect_inherits(fit_win_ao, "remstimate_window",
                info = "remwindow actor returns remstimate_window")
expect_equal(fit_win_ao$type, "actor",
             info = "remwindow actor: type is 'actor'")

# ── 18.4 remwindow diagnostics + plots (tie) ─────────────────────────────────
diag_win <- diagnostics(fit_win, reh_t, stats_t)
expect_inherits(diag_win, "diagnostics.remstimate_window",
                info = "remwindow diagnostics class correct")
.test_plot_no_error(
  plot(fit_win),
  info = "remwindow plot (coefficient trajectories) runs"
)
.test_plot_no_error(
  plot(diag_win),
  info = "remwindow plot(diagnostics) runs"
)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 19: remtribute (event-attribute model) — basic checks
# ══════════════════════════════════════════════════════════════════════════════

# Tie edgelist carrying a nominal event attribute 'type'
el_attr <- el_tie
el_attr$type <- sample(c("X", "Y"), nrow(el_attr), replace = TRUE)

# ── 19.1 Tie model, pre-computed stats route ─────────────────────────────────
reh_attr <- remify(el_attr, model = "tie", actors = actors,
                   event_type = "type", event_attributes = "type")
stats_attr <- remstats(reh_attr, tie_effects = ~ inertia() + reciprocity())
fit_attr <- remtribute(reh_attr, stats = stats_attr,
                       attribute = "type", attribute_type = "nominal")
expect_inherits(fit_attr, "remtribute",
                info = "remtribute tie (stats route) returns remtribute")

# ── 19.2 Actor model, effects-formula route ──────────────────────────────────
reh_attr_ao <- remify(el_attr, model = "actor", actors = actors,
                      directed = TRUE, event_attributes = "type")
fit_attr_ao <- remtribute(reh_attr_ao,
                          effects = ~ inertia() + reciprocity(),
                          attribute = "type", attribute_type = "nominal")
expect_inherits(fit_attr_ao, "remtribute",
                info = "remtribute actor (effects route) returns remtribute")

# ── 19.3 remtribute has no diagnostics/plot methods (attribute regression);
#         smoke-test the inspection methods it does provide (coef + summary) ────
expect_true(!is.null(coef(fit_attr)),
            info = "remtribute coef() runs")
expect_false(inherits(tryCatch(summary(fit_attr), error = function(e) e), "error"),
             info = "remtribute summary() runs")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 20: Actor Bayesian HMC — basic check
# ══════════════════════════════════════════════════════════════════════════════

# ── 20.1 Actor HMC ───────────────────────────────────────────────────────────
# Minimal HMC settings: this is a "does it run / right class" check only, not a
# convergence check, so keep iterations and leapfrog steps (L) tiny for speed.
fit_hmc_ao <- remstimate(reh_ao, stats_ao, approach = "Bayesian",
                         nsim = 20, burnin = 5, thin = 1, nchains = 1, L = 5)
expect_inherits(fit_hmc_ao, "remstimate",
                info = "actor HMC returns remstimate")
expect_equal(attr(fit_hmc_ao, "approach"), "Bayesian",
             info = "actor HMC sets approach to Bayesian")
expect_true(!is.null(fit_hmc_ao$sender_model) &&
            !is.null(fit_hmc_ao$receiver_model),
            info = "actor HMC has sender + receiver models")
expect_true(!is.null(fit_hmc_ao$sender_model$draws),
            info = "actor HMC sender model has draws")

# ── 20.2 actor HMC diagnostics + plot ────────────────────────────────────────
diag_hmc_ao <- diagnostics(fit_hmc_ao, reh_ao, stats_ao)
expect_true(!is.null(diag_hmc_ao$sender_model) &&
            !is.null(diag_hmc_ao$receiver_model),
            info = "actor HMC diagnostics: sender + receiver present")
.test_plot_no_error(
  plot(fit_hmc_ao, reh_ao, stats = stats_ao, which = 3),
  info = "actor HMC plot which=3 (recall) runs"
)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 21: Non-default list-argument passthrough (refactor guards)
# ══════════════════════════════════════════════════════════════════════════════
# Exercise sub-arguments at NON-default values. After the argument-collapse
# refactor, nfolds / lambda_select / concomitant / nrep live inside the
# penalty / mixture lists; the old TOP-LEVEL forms used here remain accepted
# (deprecated) and are routed into those lists via '...'. Section 22 exercises
# the canonical list form.

if (requireNamespace("glmnet", quietly = TRUE)) {
  # penalty: nfolds (never varied elsewhere) + lambda_select, both non-default
  fit_pen_nd <- remstimate(reh_t, stats_t,
                           penalty = list(alpha = 1),
                           nfolds = 5, lambda_select = "min")
  expect_true(!is.null(coef(fit_pen_nd)),
              info = "penalty non-defaults (nfolds=5, lambda_select='min') pass through")
}

if (requireNamespace("flexmix", quietly = TRUE)) {
  # mixture: concomitant (never tested elsewhere) + non-default nrep
  fit_mix_nd <- remstimate(reh_t, stats_t,
                           mixture = list(k = 2, random = ~ (1 | dyad)),
                           concomitant = ~ 1, nrep = 1)
  expect_true(!is.null(fit_mix_nd),
              info = "mixture with non-default concomitant + nrep passes through")
}

if (requireNamespace("shrinkem", quietly = TRUE)) {
  # shrinkem actor with non-default prior (lasso rather than horseshoe)
  fit_hs_ao_lasso <- remstimate(reh_ao, stats_ao, approach = "Bayesian",
                                penalty = list(prior = "lasso"))
  expect_inherits(fit_hs_ao_lasso, "remstimate_shrinkem",
                  info = "shrinkem actor non-default prior (lasso) passes through")
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 22: New wrapper API (remfrailty / remixture / rempenalty + bayes list)
# ══════════════════════════════════════════════════════════════════════════════

# ── 22.1 remfrailty (renamed frailty_rem; default intercepts only) ───────────
if (requireNamespace("lme4", quietly = TRUE)) {
  fit_rf <- remfrailty(reh_t, stats_t)
  expect_true(!is.null(fit_rf$backend_fit),
              info = "remfrailty tie works")
}

# ── 22.2 remixture (general mixture) + dlcrem (special case) ──────────────────
if (requireNamespace("flexmix", quietly = TRUE)) {
  fit_rmix <- remixture(reh_t, stats_t, random = ~ (1 | dyad), k = 2)
  expect_true(!is.null(fit_rmix),
              info = "remixture works")
  fit_dlc2 <- dlcrem(reh_t, stats_t, k = 2)
  expect_true(!is.null(fit_dlc2),
              info = "dlcrem (special case of remixture) works")
}

# ── 22.3 rempenalty: frequentist (lasso) + Bayesian (horseshoe) ──────────────
if (requireNamespace("glmnet", quietly = TRUE)) {
  fit_rpen <- rempenalty(reh_t, stats_t)                  # lasso default
  expect_true(!is.null(coef(fit_rpen)),
              info = "rempenalty frequentist (lasso) works")
}
if (requireNamespace("shrinkem", quietly = TRUE)) {
  fit_rpen_b <- rempenalty(reh_t, stats_t, approach = "Bayesian")  # horseshoe default
  expect_inherits(fit_rpen_b, "remstimate_shrinkem",
                  info = "rempenalty Bayesian (horseshoe) works")
}

# ── 22.4 bayes = list(...) HMC control bundle (canonical form) ────────────────
fit_bayes_list <- remstimate(reh_t, stats_t, approach = "Bayesian",
                             bayes = list(nsim = 20, burnin = 5,
                                          nchains = 1, thin = 1, L = 5))
expect_inherits(fit_bayes_list, "remstimate",
                info = "bayes=list(...) HMC control bundle works")
expect_true(!is.null(fit_bayes_list$draws),
            info = "bayes=list(...) produces draws")

# ── 22.5 penalty / mixture lists own ALL their knobs (post-refactor form) ─────
if (requireNamespace("glmnet", quietly = TRUE)) {
  fit_pen_list <- remstimate(reh_t, stats_t,
                             penalty = list(alpha = 1, nfolds = 5,
                                            lambda_select = "min"))
  expect_true(!is.null(coef(fit_pen_list)),
              info = "penalty list owns alpha/nfolds/lambda_select")
}
if (requireNamespace("flexmix", quietly = TRUE)) {
  fit_mix_list <- remstimate(reh_t, stats_t,
                             mixture = list(k = 2, random = ~ (1 | dyad),
                                            concomitant = ~ 1, nrep = 1))
  expect_true(!is.null(fit_mix_list),
              info = "mixture list owns k/random/concomitant/nrep")
}

# ── 22.6 deprecated aliases still function ───────────────────────────────────
if (requireNamespace("lme4", quietly = TRUE)) {
  fit_dep <- suppressWarnings(frailty_rem(reh_t, stats_t))
  expect_true(!is.null(fit_dep),
              info = "deprecated frailty_rem alias still works")
}

cat("\n═══ All tests completed ═══\n")
