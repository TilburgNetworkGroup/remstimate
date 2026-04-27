library(tinytest)
library(remify)
library(remstats)
library(remstimate)

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────

make_edgelist <- function(typed = FALSE) {
  el <- data.frame(
    time   = c(1, 2, 3, 3, 4, 5, 6, 7, 8, 9),
    actor1 = c(1, 2, 1, 3, 2, 1, 3, 2, 4, 1),
    actor2 = c(2, 3, 4, 5, 1, 3, 2, 4, 1, 5)
  )
  if (typed) el$type <- rep(c("a", "b"), length.out = nrow(el))
  el
}

make_attr_actors <- function() {
  data.frame(
    name = rep(1:5, each = 2),
    time = rep(c(0, 5), times = 5),
    x    = rnorm(10)
  )
}

el <- make_edgelist()

# manual riskset: all observed dyads + their mirrors (as matrix to avoid
# column-name issues with rbind on data.frames)
manual_rs <- as.data.frame(rbind(as.matrix(el[, c("actor1","actor2")]),
                                  as.matrix(el[, c("actor2","actor1")])))

# ─────────────────────────────────────────────────────────────────────────────
# 1.  TIE MODEL
# ─────────────────────────────────────────────────────────────────────────────

# tie | ordinal=FALSE | full riskset
reh <- remify(el, model = "tie", ordinal = FALSE)
st  <- remstats(reh, tie_effects = ~ inertia() + reciprocity())
fit <- remstimate(reh, stats = st, method = "MLE")

expect_inherits(fit, "remstimate",                        info = "tie/ordinal=F/full: class")
expect_equal(attr(fit, "model"), "tie",                   info = "tie/ordinal=F/full: model attr")
expect_equal(attr(fit, "ordinal"), FALSE,                 info = "tie/ordinal=F/full: ordinal attr")
expect_true(fit$converged,                                info = "tie/ordinal=F/full: converged")
expect_equal(names(fit$coefficients),
             c("baseline", "inertia", "reciprocity"),     info = "tie/ordinal=F/full: coef names")

# tie | ordinal=TRUE | full riskset (no baseline)
reh_ord <- remify(el, model = "tie", ordinal = TRUE)
st_ord  <- remstats(reh_ord, tie_effects = ~ inertia() + reciprocity())
fit_ord <- remstimate(reh_ord, stats = st_ord, method = "MLE")

expect_inherits(fit_ord, "remstimate",                    info = "tie/ordinal=T/full: class")
expect_equal(attr(fit_ord, "ordinal"), TRUE,              info = "tie/ordinal=T/full: ordinal attr")
expect_equal(names(fit_ord$coefficients),
             c("inertia", "reciprocity"),                 info = "tie/ordinal=T/full: coef names (no baseline)")

# tie | active riskset
reh_act <- remify(el, model = "tie", ordinal = FALSE, riskset = "active")
st_act  <- remstats(reh_act, tie_effects = ~ inertia())
fit_act <- remstimate(reh_act, stats = st_act, method = "MLE")

expect_inherits(fit_act, "remstimate",                    info = "tie/active riskset: class")
expect_true("inertia" %in% names(fit_act$coefficients),  info = "tie/active riskset: inertia in coefs")

# tie | manual riskset
reh_man <- remify(el, model = "tie", ordinal = FALSE,
                   riskset = "manual", manual.riskset = manual_rs)
st_man  <- remstats(reh_man, tie_effects = ~ inertia())
fit_man <- remstimate(reh = reh_man, stats = st_man, method = "MLE")

expect_inherits(fit_man, "remstimate",                    info = "tie/manual riskset: class")
expect_true("inertia" %in% names(fit_man$coefficients),  info = "tie/manual riskset: inertia in coefs")
expect_true(fit_man$converged,                            info = "tie/manual riskset: converged")

# tie | typed events | consider_type = "ignore" (default)
el_t  <- make_edgelist(typed = TRUE)
reh_t <- remify(el_t, model = "tie", ordinal = FALSE)
st_t  <- remstats(reh_t, tie_effects = ~ inertia())
fit_t <- remstimate(reh_t, stats = st_t, method = "MLE")

expect_inherits(fit_t, "remstimate",                      info = "tie/typed/ignore: class")
expect_true("inertia" %in% names(fit_t$coefficients),    info = "tie/typed/ignore: inertia in coefs")

# tie | typed events | consider_type = "separate"
st_ts  <- remstats(reh_t, tie_effects = ~ inertia(consider_type = "separate"))
fit_ts <- remstimate(reh_t, stats = st_ts, method = "MLE")

expect_inherits(fit_ts, "remstimate",                     info = "tie/typed/separate: class")
expect_true(length(fit_ts$coefficients) > 2,             info = "tie/typed/separate: more coefs than untyped")

# tie | actor attributes
attr  <- make_attr_actors()
st_a  <- remstats(reh, tie_effects = ~ inertia() + send(variable = "x"),
                  attr_actors = attr)
fit_a <- remstimate(reh, stats = st_a, method = "MLE")

expect_inherits(fit_a, "remstimate",                      info = "tie/attr: class")
expect_true("inertia" %in% names(fit_a$coefficients),    info = "tie/attr: inertia in coefs")
expect_true("send_x"  %in% names(fit_a$coefficients),    info = "tie/attr: send.x in coefs")

# ─────────────────────────────────────────────────────────────────────────────
# 2.  ACTOR MODEL
# ─────────────────────────────────────────────────────────────────────────────

# actor | ordinal=FALSE | sender only
reh_a <- remify(el, model = "actor", ordinal = FALSE)
st_s  <- remstats(reh_a, sender_effects  = ~ outdegreeSender(),
                         receiver_effects = NULL)
fit_s <- remstimate(reh_a, stats = st_s, method = "MLE")

expect_inherits(fit_s, "remstimate",                      info = "actor/sender-only: class")
expect_equal(attr(fit_s, "model"), "actor",               info = "actor/sender-only: model attr")
expect_false(is.null(fit_s$sender_model),                 info = "actor/sender-only: sender_model present")
expect_true(is.null(fit_s$receiver_model),                info = "actor/sender-only: receiver_model absent")
expect_true(fit_s$sender_model$converged,                 info = "actor/sender-only: converged")
expect_true("baseline" %in% names(fit_s$sender_model$coefficients),
                                                          info = "actor/sender-only/ordinal=F: baseline present")

# actor | ordinal=TRUE | sender only (no baseline)
reh_ao <- remify(el, model = "actor", ordinal = TRUE)
st_so  <- remstats(reh_ao, sender_effects  = ~ outdegreeSender(),
                           receiver_effects = NULL)
fit_so <- remstimate(reh = reh_ao, stats = st_so, method = "MLE")

expect_false("baseline" %in% names(fit_so$sender_model$coefficients),
                                                          info = "actor/sender-only/ordinal=T: no baseline")

# actor | ordinal=FALSE | receiver only (no baseline in receiver)
st_r  <- remstats(reh_a, sender_effects  = NULL,
                         receiver_effects = ~ inertia() + reciprocity())
fit_r <- remstimate(reh_a, stats = st_r, method = "MLE")

expect_inherits(fit_r, "remstimate",                      info = "actor/receiver-only: class")
expect_true(is.null(fit_r$sender_model),                  info = "actor/receiver-only: sender_model absent")
expect_false(is.null(fit_r$receiver_model),               info = "actor/receiver-only: receiver_model present")
expect_false("baseline" %in% names(fit_r$receiver_model$coefficients),
                                                          info = "actor/receiver-only: no baseline in receiver")
expect_true(fit_r$receiver_model$converged,               info = "actor/receiver-only: converged")

# actor | ordinal=FALSE | sender + receiver
st_sr  <- remstats(reh_a, sender_effects   = ~ outdegreeSender(),
                          receiver_effects  = ~ inertia() + reciprocity())
fit_sr <- remstimate(reh_a, stats = st_sr, method = "MLE")

expect_inherits(fit_sr, "remstimate",                     info = "actor/both/ordinal=F: class")
expect_true(fit_sr$sender_model$converged,                info = "actor/both/ordinal=F: sender converged")
expect_true(fit_sr$receiver_model$converged,              info = "actor/both/ordinal=F: receiver converged")
expect_equal(names(fit_sr$sender_model$coefficients),
             c("baseline", "outdegreeSender"),            info = "actor/both/ordinal=F: sender coef names")
expect_equal(names(fit_sr$receiver_model$coefficients),
             c("inertia", "reciprocity"),                 info = "actor/both/ordinal=F: receiver coef names")

# actor | ordinal=TRUE | sender + receiver (no baseline)
fit_sr_ord <- remstimate(reh_a, stats = st_sr, method = "MLE")

expect_inherits(fit_sr_ord, "remstimate",                 info = "actor/both/ordinal=T: class")
expect_true("baseline" %in% names(fit_sr_ord$sender_model$coefficients),
                                                          info = "actor/both/ordinal=T: no baseline in sender")

# actor | tied events: stat array dimensions and no crash
n_unique_times <- length(unique(el$time))   # 9
n_events       <- nrow(el)                  # 10

expect_equal(dim(st_sr$sender_stats)[1],   n_unique_times-1, # minus 1 for new default start=2
             info = "actor/tied: sender stats rows == unique times")
expect_equal(dim(st_sr$receiver_stats)[1], n_unique_times - 1,
             info = "actor/tied: receiver stats rows == n events")
expect_true(fit_sr$receiver_model$converged,
             info = "actor/tied: receiver model converges despite ties")

# actor | typed events
el_t2  <- make_edgelist(typed = TRUE)
reh_at <- remify(el_t2, model = "actor", ordinal = FALSE)
st_at  <- remstats(reh_at, sender_effects   = ~ outdegreeSender(),
                           receiver_effects  = ~ inertia())
fit_at <- remstimate(reh_at, stats = st_at, method = "MLE")

expect_inherits(fit_at, "remstimate",                     info = "actor/typed: class")
expect_false(is.null(fit_at$sender_model),                info = "actor/typed: sender_model present")
expect_false(is.null(fit_at$receiver_model),              info = "actor/typed: receiver_model present")

# actor | actor attributes
attr_df <- make_attr_actors()
st_aa  <- remstats(reh_a,
                   sender_effects   = ~ outdegreeSender() + send(variable = "x"),
                   receiver_effects = ~ inertia() + receive(variable = "x"),
                   attr_actors = attr_df)
fit_aa <- remstimate(reh_a, stats = st_aa, method = "MLE")

expect_inherits(fit_aa, "remstimate",                     info = "actor/attr: class")
expect_true("send"    %in% names(fit_aa$sender_model$coefficients),
                                                          info = "actor/attr: send.x in sender coefs")
expect_true("receive" %in% names(fit_aa$receiver_model$coefficients),
                                                          info = "actor/attr: receive.x in receiver coefs")

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Output structure sanity
# ─────────────────────────────────────────────────────────────────────────────

# tie model: all expected fields present (flat structure)
for (field in c("coefficients", "loglik", "AIC", "AICC", "BIC",
                "vcov", "se", "null.deviance", "residual.deviance",
                "model.deviance", "converged", "iterations",
                "df.null", "df.model", "df.residual")) {
  expect_false(is.null(fit_a[[field]]),
               info = paste("tie output field present:", field))
}

# actor model: all expected fields in sender and receiver submodels
for (nm in c("sender_model", "receiver_model")) {
  m <- fit_sr[[nm]]
  for (field in c("coefficients", "loglik", "AIC", "AICC", "BIC",
                  "vcov", "se", "null.deviance", "residual.deviance",
                  "model.deviance", "converged", "iterations",
                  "df.null", "df.model", "df.residual")) {
    expect_false(is.null(m[[field]]),
                 info = paste("actor", nm, "field present:", field))
  }
}

