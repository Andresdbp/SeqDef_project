# Recover the Grace runtime data from the SLURM .out log (the job OOM-killed at
# n=150k before write.csv ran, but every completed n was printed to stdout).
lines <- readLines("results/grace_runtime_18719664.out")
dl <- grep("^n=", lines, value = TRUE)
rx <- "n=\\s*([0-9]+) reps=\\s*([0-9]+) time=\\s*([0-9.]+)s \\[([0-9.]+), ([0-9.]+)\\] peak=\\s*([0-9]+) MB"
m  <- regmatches(dl, regexec(rx, dl))
df <- do.call(rbind, lapply(m, function(x) data.frame(
  n = as.numeric(x[2]), reps = as.numeric(x[3]), time_s = as.numeric(x[4]),
  time_lo = as.numeric(x[5]), time_hi = as.numeric(x[6]), peak_mb = as.numeric(x[7]))))
df$status <- "ok"
df$theory_gb <- 8 * df$n^2 / 1e9
write.csv(df, "results/runtime_benchmark_grace.csv", row.names = FALSE)

cat("parsed", nrow(df), "rows; n", min(df$n), "to", max(df$n), "\n")
cb <- df[df$n >= 1000, ]
cat("time slope (n>=1000):", round(coef(lm(log(time_s) ~ log(n), cb))[2], 3),
    "| memory slope:", round(coef(lm(log(peak_mb) ~ log(n), cb))[2], 3), "\n")
cat(sprintf("n=100000: %.0f s, %.0f GB peak\n", df$time_s[df$n==100000], df$peak_mb[df$n==100000]/1024))
