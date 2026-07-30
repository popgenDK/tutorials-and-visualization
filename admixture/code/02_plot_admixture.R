data_dir <- "../tutorial_data/admixture"
results_dir <- "results"
figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

k <- Sys.getenv("K", "3")
q <- read.table(file.path(results_dir, sprintf("example.pruned.K%s.Q", k)))

fam_file <- file.path(results_dir, "example.pruned.fam")
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

cv_file <- file.path(results_dir, "admixture.cv_errors.txt")
if (file.exists(cv_file)) {
  cv <- readLines(cv_file)
  k_values <- as.integer(sub(".*K=([0-9]+).*", "\\1", cv))
  err <- as.numeric(sub(".*: ", "", cv))

  png(file.path(figures_dir, "admixture_cv.png"), width = 800, height = 600, res = 150)
  plot(k_values, err, type = "b", pch = 19, xlab = "K", ylab = "Cross-validation error")
  dev.off()
}
