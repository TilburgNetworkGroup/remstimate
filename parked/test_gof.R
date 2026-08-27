# Synthetic test of the GOF internals (no remstimate needed).
# Builds a small tie-model event sequence from a known model, then checks:
#  (a) under the TRUE parameters the observed counts sit inside the reference
#      band about (1-alpha) of the time  -> calibrated
#  (b) under badly WRONG parameters many actors fall outside -> has power

source("/home/claude/gof_prototype.R")

set.seed(42)
N <- 12L
# full directed riskset
dm <- expand.grid(a2 = 1:N, a1 = 1:N)
dm <- dm[dm$a1 != dm$a2, c("a1", "a2")]
dm$dyadID <- seq_len(nrow(dm))
D <- nrow(dm)
M <- 400L
P <- 2L

# stats array: D x P x M, with actor-level heterogeneity so counts are uneven
sender_pull <- rnorm(N, 0, 0.8)
recv_pull   <- rnorm(N, 0, 0.8)
stats_3d <- array(0, dim = c(D, P, M))
for (m in seq_len(M)) {
  stats_3d[, 1L, m] <- sender_pull[dm$a1] + rnorm(D, 0, 0.2)   # "activity"
  stats_3d[, 2L, m] <- recv_pull[dm$a2]   + rnorm(D, 0, 0.2)   # "popularity"
}
beta_true <- c(1.0, 0.7)

# generate the observed sequence from the true model
obs <- integer(M)
for (m in seq_len(M)) {
  p <- .gof_softmax(beta_true, 0, stats_3d, m, seq_len(D))
  obs[m] <- sample.int(D, 1L, prob = p)
}
obs_ids <- as.list(obs)

a1_of <- dm$a1[order(dm$dyadID)]
a2_of <- dm$a2[order(dm$dyadID)]

run <- function(beta, R = 1000L) {
  blk <- .gof_sample_block(beta, 0, stats_3d, obs_ids, R = R)
  to_role <- function(mat, lk) {
    m <- matrix(lk[mat], nrow = nrow(mat), ncol = ncol(mat)); storage.mode(m) <- "integer"; m
  }
  list(
    sender = .gof_tally(list(sim = to_role(blk$sim, a1_of), obs = a1_of[blk$obs]),
                        nbins = N),
    receiver = .gof_tally(list(sim = to_role(blk$sim, a2_of), obs = a2_of[blk$obs]),
                          nbins = N),
    dyad = .gof_tally(blk, nbins = D)
  )
}

cat("=== (a) TRUE parameters ===\n")
set.seed(1); g_true <- run(beta_true)
for (nm in c("sender", "receiver", "dyad")) {
  tb <- g_true[[nm]]$table
  act <- tb[tb$observed > 0 | tb$sim_mean > 0, ]
  cat(sprintf("%-9s units=%3d  inside 95%% band: %5.1f%%  max|z| = %.2f\n",
              nm, nrow(act), 100 * mean(!act$outside), max(abs(act$z), na.rm = TRUE)))
}
cat("\nsender table (true params):\n")
print(format(g_true$sender$table[, c("id","observed","sim_mean","sim_lo","sim_hi","z","p_two")],
             digits = 3), row.names = FALSE)

cat("\n=== (b) WRONG parameters (signs flipped) ===\n")
set.seed(1); g_bad <- run(c(-1.0, -0.7))
for (nm in c("sender", "receiver", "dyad")) {
  tb <- g_bad[[nm]]$table
  act <- tb[tb$observed > 0 | tb$sim_mean > 0, ]
  cat(sprintf("%-9s units=%3d  inside 95%% band: %5.1f%%  max|z| = %.2f\n",
              nm, nrow(act), 100 * mean(!act$outside), max(abs(act$z), na.rm = TRUE)))
}

cat("\n=== (c) totals conserved? ===\n")
blk <- .gof_sample_block(beta_true, 0, stats_3d, obs_ids, R = 50L)
cat("observed events:", length(blk$obs), " simulated per replicate:",
    unique(colSums(!is.na(blk$sim))), "\n")

cat("\n=== (d) simultaneous events handled? ===\n")
obs_sim <- obs_ids; obs_sim[[3]] <- c(obs[3], obs[4])   # two events at time 3
b2 <- .gof_sample_block(beta_true, 0, stats_3d, obs_sim, R = 20L)
cat("rows:", nrow(b2$sim), "(expected", M + 1L, ")\n")

cat("\n=== (e) baseline is irrelevant (softmax invariance)? ===\n")
set.seed(7); s1 <- .gof_sample_block(beta_true, 0,    stats_3d, obs_ids, R = 5L)
set.seed(7); s2 <- .gof_sample_block(beta_true, 99.5, stats_3d, obs_ids, R = 5L)
cat("identical draws:", identical(s1$sim, s2$sim), "\n")

cat("\n=== (f) overflow guard vs the current recall code ===\n")
big <- array(0, dim = c(5, 1, 1)); big[, 1, 1] <- c(0, 200, 400, 800, 1000)
cat("with max-subtraction: ", paste(round(.gof_softmax(1, 0, big, 1, 1:5), 4), collapse = " "), "\n")
naive <- exp(as.numeric(0 + matrix(big[1:5,,1], ncol = 1) %*% 1)); naive <- naive / sum(naive)
cat("without (as in .recall_block_3d):", paste(naive, collapse = " "), "\n")

cat("\n=== (g) plot renders ===\n")
png("/home/claude/gof_test_plot.png", width = 1500, height = 1100, res = 130)
par(mfrow = c(2, 2))
.gof_panel(g_true$sender,  "Sender - true params")
.gof_panel(g_true$receiver,"Receiver - true params")
.gof_panel(g_bad$sender,   "Sender - wrong params")
.gof_panel(g_true$dyad,    "Dyads - true params", max_units = 20)
dev.off()
cat("written\n")
