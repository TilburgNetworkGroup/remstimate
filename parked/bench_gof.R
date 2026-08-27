source("/home/claude/gof_prototype.R")

mk <- function(N, M, P = 3L) {
  D <- N * (N - 1L)
  list(N = N, M = M, D = D,
       stats_3d = array(rnorm(D * P * M, 0, 0.5), dim = c(D, P, M)),
       obs = as.list(sample.int(D, M, replace = TRUE)),
       beta = rnorm(P, 0, 0.5))
}

bench <- function(N, M, R) {
  s <- mk(N, M)
  gc()
  t_draw <- system.time(
    blk <- .gof_sample_block(s$beta, 0, s$stats_3d, s$obs, R = R))[["elapsed"]]
  t_tal <- system.time(.gof_tally(blk, nbins = s$D))[["elapsed"]]
  mem   <- as.numeric(object.size(blk$sim)) / 1024^2
  data.frame(N = N, M = M, D = s$D, R = R,
             draw_s = round(t_draw, 2), tally_s = round(t_tal, 2),
             total_s = round(t_draw + t_tal, 2), sim_MB = round(mem, 1))
}

cat("### scaling ###\n")
res <- do.call(rbind, list(
  bench(  20,  500, 1000L),
  bench(  50,  500, 1000L),
  bench( 100,  500, 1000L),
  bench(  50, 2000, 1000L),
  bench(  50, 5000, 1000L),
  bench( 100, 5000, 1000L),
  bench( 200, 2000, 1000L),
  bench(  50, 5000,  200L),
  bench(  50, 5000,   50L)
))
print(res, row.names = FALSE)

# ---- how much of that is the softmax diagnostics() already computes? ----
cat("\n### softmax-only cost (already paid by .recall_block_3d) ###\n")
s <- mk(100, 5000)
t_sm <- system.time({
  for (m in seq_len(s$M)) .gof_softmax(s$beta, 0, s$stats_3d, m, seq_len(s$D))
})[["elapsed"]]
cat(sprintf("N=100 M=5000 D=%d : softmax pass = %.2fs\n", s$D, t_sm))

# ---- streaming variant: accumulate counts, never hold the M x R matrix ----
cat("\n### streaming variant (counts accumulated per event) ###\n")
.gof_stream <- function(pars, baseline, stats_3d, obs_ids, lookup, nbins,
                        valid_ids = NULL, R = 1000L) {
  M <- dim(stats_3d)[3L]; Dall <- dim(stats_3d)[1L]
  acc <- matrix(0L, nrow = nbins, ncol = R)
  for (m in seq_len(M)) {
    n_m <- length(obs_ids[[m]]); if (!n_m) next
    valid <- if (is.null(valid_ids)) seq_len(Dall) else valid_ids[[m]]
    if (!length(valid)) next
    p <- .gof_softmax(pars, baseline, stats_3d, m, valid)
    d <- lookup[valid[sample.int(length(valid), n_m * R, replace = TRUE, prob = p)]]
    # d has n_m*R entries laid out replicate-major; bin into (unit, replicate)
    rep_id <- rep(seq_len(R), each = n_m)
    idx <- (rep_id - 1L) * nbins + as.integer(d)
    acc[] <- acc + matrix(tabulate(idx, nbins = nbins * R), nrow = nbins)
  }
  acc
}
s <- mk(100, 5000)
lk <- rep_len(seq_len(s$N), s$D)   # stand-in dyad -> sender map
gc()
t_str <- system.time(a <- .gof_stream(s$beta, 0, s$stats_3d, s$obs, lk,
                                      nbins = s$N, R = 1000L))[["elapsed"]]
cat(sprintf("N=100 M=5000 R=1000 : %.2fs, acc matrix %.1f MB (vs %.1f MB for M x R)\n",
            t_str, as.numeric(object.size(a))/1024^2, 5000*1000*4/1024^2))
