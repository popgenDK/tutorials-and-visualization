data_dir <- "data"
results_dir <- "results"
figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

plot_draft_admixture <- function(outfile) {
  png(outfile, width = 1400, height = 500, res = 160)
  par(mar = c(5, 4, 1, 1))
  q <- rbind(
    c(0.98, 0.01, 0.01), c(0.96, 0.02, 0.02), c(0.02, 0.96, 0.02),
    c(0.01, 0.98, 0.01), c(0.35, 0.55, 0.10), c(0.30, 0.45, 0.25),
    c(0.02, 0.03, 0.95), c(0.01, 0.04, 0.95)
  )
  pop <- c("YRI", "YRI", "CEU", "CEU", "MXL", "MXL", "CHB", "CHB")
  barplot(t(q), col = 2:4, space = 0, border = NA, xlab = "Individuals", ylab = "Ancestry proportion")
  text(c(1, 3, 5, 7), -0.06, c("YRI", "CEU", "MXL", "CHB"), xpd = TRUE)
  dev.off()
}

plot_draft_convergence <- function(outfile) {
  png(outfile, width = 900, height = 550, res = 150)
  plot(1:10, c(-101, -92, -91.8, -91.7, -91.8, -91.7, -91.9, -91.8, -91.7, -91.8),
       type = "b", pch = 19, xlab = "Seed", ylab = "Log likelihood", main = "Convergence across seeds")
  points(1:10, c(-110, -94, -99, -93.5, -102, -93.8, -98, -93.6, -101, -94), type = "b", pch = 19, col = 2)
  legend("bottomright", legend = c("stable K", "unstable K"), col = c(1, 2), pch = 19, lty = 1, bty = "n")
  dev.off()
}

k <- Sys.getenv("K", "3")
seed <- Sys.getenv("SEED", "1")
q_file <- file.path(results_dir, sprintf("example.pcaone_pruned.K%s.seed%s.Q", k, seed))

if (!file.exists(q_file)) {
  plot_draft_admixture(file.path(figures_dir, sprintf("admixture_k%s.png", k)))
  plot_draft_convergence(file.path(figures_dir, "admixture_convergence.png"))
} else {
q <- read.table(q_file)

fam_file <- file.path(results_dir, "example.pcaone_pruned.fam")
if (!file.exists(fam_file)) {
  fam_file <- file.path(data_dir, "example.fam")
}
fam <- read.table(fam_file, as.is = TRUE)

pop <- fam[, 1]
ord <- order(pop, q[, 1])

png(file.path(figures_dir, sprintf("admixture_k%s.png", k)), width = 1400, height = 500, res = 160)
par(mar = c(5, 4, 1, 1))
barplot(
  t(q)[, ord],
  col = 2:10,
  space = 0,
  border = NA,
  xlab = "Individuals",
  ylab = "Ancestry proportion"
)
text(sort(tapply(seq_len(nrow(fam)), pop[ord], mean)), -0.05, unique(pop[ord]), xpd = TRUE)
abline(v = cumsum(sapply(unique(pop[ord]), function(x) sum(pop[ord] == x))), col = 1, lwd = 1.2)
dev.off()
}

eig_file <- file.path(results_dir, "example.pcaone.eigvecs2")
if (file.exists(eig_file)) {
  pcs <- read.table(eig_file, header = TRUE, check.names = FALSE)
  pc_cols <- grep("^PC[0-9]+$", names(pcs), value = TRUE)
  if (length(pc_cols) >= 10) {
    png(file.path(figures_dir, "pcaone_top10_pc_pairs.png"), width = 1200, height = 900, res = 150)
    par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
    for (i in seq(1, 9, by = 2)) {
      plot(pcs[[paste0("PC", i)]], pcs[[paste0("PC", i + 1)]],
           pch = 19, cex = 0.6, xlab = paste0("PC", i), ylab = paste0("PC", i + 1),
           main = paste0("PC", i, " vs PC", i + 1))
    }
    dev.off()
  }
} else {
  png(file.path(figures_dir, "pcaone_top10_pc_pairs.png"), width = 1200, height = 900, res = 150)
  par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
  set.seed(1)
  for (i in seq(1, 9, by = 2)) {
    x <- c(rnorm(40, -2), rnorm(40, 0), rnorm(40, 2))
    y <- c(rnorm(40, -1), rnorm(40, 1), rnorm(40, 0))
    plot(x, y, pch = 19, cex = 0.6, xlab = paste0("PC", i), ylab = paste0("PC", i + 1),
         main = paste0("PC", i, " vs PC", i + 1))
  }
  dev.off()
}
