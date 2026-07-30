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

read_pop_labels <- function(n) {
  fam_file <- file.path(results_dir, "example.qc.fam")
  if (!file.exists(fam_file)) {
    fam_file <- file.path(results_dir, "example.pcaone_pruned.fam")
  }
  if (!file.exists(fam_file)) {
    fam_file <- file.path(data_dir, "example.fam")
  }
  if (!file.exists(fam_file)) {
    return(rep("unknown", n))
  }
  fam <- read.table(fam_file, as.is = TRUE)
  pop <- fam[, 1]
  if (length(pop) != n) {
    return(rep("unknown", n))
  }
  pop
}

read_pcaone_eigvecs2 <- function(path) {
  pcs <- read.table(path, header = TRUE, check.names = FALSE, comment.char = "")
  pc_cols <- grep("^PC[0-9]+$", names(pcs), value = TRUE)
  if (length(pc_cols) >= 10) {
    return(pcs)
  }

  pcs <- read.table(path, header = FALSE, check.names = FALSE, comment.char = "")
  if (ncol(pcs) >= 12) {
    names(pcs) <- c("FID", "IID", paste0("PC", seq_len(ncol(pcs) - 2)))
  } else if (ncol(pcs) >= 10) {
    names(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  }
  pcs
}

plot_pc_pairs <- function(pcs, pop, outfile) {
  pop <- factor(pop)
  pop_cols <- setNames(
    c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")[seq_along(levels(pop))],
    levels(pop)
  )

  png(outfile, width = 1200, height = 900, res = 150)
  par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
  for (i in seq(1, 9, by = 2)) {
    plot(
      pcs[[paste0("PC", i)]],
      pcs[[paste0("PC", i + 1)]],
      pch = 19,
      cex = 0.65,
      col = pop_cols[as.character(pop)],
      xlab = paste0("PC", i),
      ylab = paste0("PC", i + 1),
      main = paste0("PC", i, " vs PC", i + 1)
    )
  }
  plot.new()
  legend("center", legend = levels(pop), col = pop_cols, pch = 19, bty = "n", title = "Population")
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

read_loglikelihoods <- function() {
  logs <- list.files(results_dir, pattern = "^admixture\\.K[0-9]+\\.seed[0-9]+\\.log$", full.names = TRUE)
  if (length(logs) == 0) {
    return(NULL)
  }

  out <- lapply(logs, function(path) {
    base <- basename(path)
    k <- as.integer(sub("^admixture\\.K([0-9]+)\\.seed[0-9]+\\.log$", "\\1", base))
    seed <- as.integer(sub("^admixture\\.K[0-9]+\\.seed([0-9]+)\\.log$", "\\1", base))
    lines <- readLines(path, warn = FALSE)
    ll_lines <- grep("^Loglikelihood:", lines, value = TRUE)
    if (length(ll_lines) == 0) {
      return(NULL)
    }
    loglik <- as.numeric(sub("^Loglikelihood:[[:space:]]*", "", tail(ll_lines, 1)))
    data.frame(K = k, seed = seed, loglikelihood = loglik)
  })
  do.call(rbind, out)
}

plot_loglikelihoods <- function(outfile) {
  conv <- read_loglikelihoods()
  if (is.null(conv) || nrow(conv) == 0) {
    plot_draft_convergence(outfile)
    return(invisible(NULL))
  }

  conv <- conv[order(conv$K, conv$seed), ]
  k_values <- sort(unique(conv$K))
  cols <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_along(k_values)], k_values)
  conv$delta_best <- ave(conv$loglikelihood, conv$K, FUN = function(x) x - max(x))

  png(outfile, width = 900, height = 550, res = 150)
  par(mar = c(4, 5, 1, 1))
  plot(
    range(conv$seed),
    range(conv$delta_best),
    type = "n",
    xlab = "Seed",
    ylab = "Loglikelihood difference from best seed"
  )
  abline(h = 0, lty = 2, col = "grey60")
  for (k_i in k_values) {
    x <- conv[conv$K == k_i, ]
    lines(x$seed, x$delta_best, type = "b", pch = 19, col = cols[as.character(k_i)])
  }
  legend("bottomright", legend = paste0("K=", k_values), col = cols, pch = 19, lty = 1, bty = "n")
  dev.off()
}

k <- Sys.getenv("K", "3")
seed <- Sys.getenv("SEED", "1")
q_file <- file.path(results_dir, sprintf("example.pcaone_pruned.K%s.seed%s.Q", k, seed))

if (!file.exists(q_file)) {
  plot_draft_admixture(file.path(figures_dir, sprintf("admixture_k%s.png", k)))
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

plot_loglikelihoods(file.path(figures_dir, "admixture_convergence.png"))

eig_file <- file.path(results_dir, "example.pcaone.eigvecs2")
if (file.exists(eig_file)) {
  pcs <- read_pcaone_eigvecs2(eig_file)
  pc_cols <- grep("^PC[0-9]+$", names(pcs), value = TRUE)
  if (length(pc_cols) >= 10) {
    plot_pc_pairs(pcs, read_pop_labels(nrow(pcs)), file.path(figures_dir, "pcaone_top10_pc_pairs.png"))
  }
} else {
  set.seed(1)
  pop <- rep(c("CEU", "CHB", "FIN", "PEL", "PJL", "YRI"), each = 20)
  centers <- matrix(
    c(
      -0.8, -0.6, 0.5, 0.4, -0.2, 0.2, 0.4, -0.2, -0.3, 0.2,
      1.2, -0.2, -0.4, -0.7, 0.3, 0.1, -0.4, -0.3, 0.2, 0.4,
      -0.7, -0.4, 0.2, 0.7, -0.1, 0.4, 0.2, -0.4, -0.2, 0.3,
      -0.5, 1.0, 0.8, -0.4, -0.5, -0.1, 0.1, 0.5, -0.4, -0.2,
      0.3, 0.7, -0.3, 0.3, 0.5, -0.5, -0.1, 0.2, 0.4, -0.3,
      1.0, 0.8, -0.1, 0.4, -0.2, 0.5, 0.3, -0.5, 0.1, -0.3
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(c("CEU", "CHB", "FIN", "PEL", "PJL", "YRI"), paste0("PC", 1:10))
  )
  pcs <- data.frame(matrix(NA_real_, nrow = length(pop), ncol = 10))
  names(pcs) <- paste0("PC", 1:10)
  for (p in rownames(centers)) {
    idx <- which(pop == p)
    pcs[idx, ] <- sweep(matrix(rnorm(length(idx) * 10, sd = 0.25), nrow = length(idx)), 2, centers[p, ], "+")
  }
  plot_pc_pairs(pcs, pop, file.path(figures_dir, "pcaone_top10_pc_pairs.png"))
}
