#!/usr/bin/env Rscript
# hello_world.R — a tiny "does my Talapas R job work?" demo.
#
# Run it on an interactive compute node:
#     module load R
#     Rscript hello_world.R
# or submit it as a batch job via the wrapper:
#     sbatch --job-name=hello run_rscript.sbatch hello_world.R
#
# It uses only base R (no extra packages), prints where it is running, and
# writes two output files you can open directly in VS Code:
#     hello_world_plot.pdf   — a simple plot
#     hello_world_summary.csv — a one-row summary table

cat("Hello from", Sys.info()[["nodename"]], "at", format(Sys.time()), "\n")
cat("R version:", R.version.string, "\n")

set.seed(2026)
x <- 1:100
y <- cumsum(rnorm(100))

# 1) a PDF plot
pdf("hello_world_plot.pdf", width = 7, height = 4.5)
plot(x, y, type = "l", lwd = 2, col = "#1c7293",
     main = "hello_world.R — random walk",
     xlab = "step", ylab = "value")
dev.off()

# 2) a one-row summary CSV
summary_df <- data.frame(
  node     = Sys.info()[["nodename"]],
  n_points = length(y),
  mean_y   = round(mean(y), 3),
  sd_y     = round(sd(y), 3),
  min_y    = round(min(y), 3),
  max_y    = round(max(y), 3),
  when     = format(Sys.time())
)
write.csv(summary_df, "hello_world_summary.csv", row.names = FALSE)

cat("Wrote hello_world_plot.pdf and hello_world_summary.csv\n")
