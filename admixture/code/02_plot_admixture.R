data_dir <- "data"
results_dir <- "results"
figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

read_pop_labels <- function(pcs) {
  fam_file <- file.path(results_dir, "example.qc.fam")
  if (!file.exists(fam_file)) {
    fam_file <- file.path(results_dir, "example.pcaone_pruned.fam")
  }
  if (!file.exists(fam_file)) {
    fam_file <- file.path(data_dir, "example.fam")
  }
  if (!file.exists(fam_file)) {
    return(rep("unknown", nrow(pcs)))
  }
  fam <- read.table(fam_file, as.is = TRUE)
  names(fam)[1:2] <- c("#FID", "IID")
  if (all(c("#FID", "IID") %in% names(pcs))) {
    pop <- fam[match(pcs$IID, fam$IID), "#FID"]
  } else {
    pop <- fam[, "#FID"]
  }
  if (length(pop) != nrow(pcs) || any(is.na(pop))) {
    return(rep("unknown", nrow(pcs)))
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

read_loglikelihoods <- function() {
  likes <- list.files(results_dir, pattern = "^admixture\\.K[0-9]+\\.likes$", full.names = TRUE)
  if (length(likes) > 0) {
    out <- lapply(likes, function(path) {
      k <- as.integer(sub("^admixture\\.K([0-9]+)\\.likes$", "\\1", basename(path)))
      x <- read.table(path, col.names = c("seed", "loglikelihood"))
      if (nrow(x) == 0) {
        return(NULL)
      }
      data.frame(K = k, seed = x$seed, loglikelihood = x$loglikelihood)
    })
    out <- Filter(Negate(is.null), out)
    if (length(out) > 0) {
      return(do.call(rbind, out))
    }
  }

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

plot_evaladmix <- function(k, seed = 1, suffix = paste0("seed", seed)) {
  cor_file <- file.path(results_dir, sprintf("evaladmix.K%s.%s.corres", k, suffix))
  if (!file.exists(cor_file) && suffix == paste0("seed", seed)) {
    cor_file <- file.path(results_dir, sprintf("evaladmix.K%s.seed%s.corres", k, seed))
  }
  if (!file.exists(cor_file)) {
    return(FALSE)
  }

  r <- as.matrix(read.table(cor_file))
  fam_file <- file.path(results_dir, "example.pcaone_pruned.fam")
  if (!file.exists(fam_file)) {
    fam_file <- file.path(data_dir, "example.fam")
  }
  pop <- read.table(fam_file, as.is = TRUE)[, 1]
  ord <- order(pop)
  r <- r[ord, ord]
  pop <- pop[ord]

  pal <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
  lim <- max(abs(r), na.rm = TRUE)
  png(file.path(figures_dir, sprintf("evaladmix_k%s_%s.png", k, suffix)), width = 900, height = 800, res = 150)
  par(mar = c(4, 4, 2, 5))
  image(
    seq_len(nrow(r)),
    seq_len(ncol(r)),
    r[nrow(r):1, ],
    col = pal,
    zlim = c(-lim, lim),
    axes = FALSE,
    xlab = "Individuals",
    ylab = "Individuals",
    main = sprintf("evalAdmix residual correlations, K=%s", k)
  )
  axis(1, at = tapply(seq_along(pop), pop, mean), labels = names(table(pop)), las = 2, cex.axis = 0.75)
  axis(2, at = nrow(r) - tapply(seq_along(pop), pop, mean) + 1, labels = names(table(pop)), las = 2, cex.axis = 0.75)
  abline(v = cumsum(table(pop)) + 0.5, col = "grey30", lwd = 0.8)
  abline(h = nrow(r) - cumsum(table(pop)) + 0.5, col = "grey30", lwd = 0.8)
  dev.off()
  TRUE
}

plot_admixture_q <- function(k, seed = 1) {
  q_file <- file.path(results_dir, sprintf("example.pcaone_pruned.K%s.seed%s.Q", k, seed))
  if (!file.exists(q_file)) {
    return(FALSE)
  }

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
    col = seq_len(ncol(q)) + 1,
    space = 0,
    border = NA,
    xlab = "Individuals",
    ylab = "Ancestry proportion"
  )
  text(sort(tapply(seq_len(nrow(fam)), pop[ord], mean)), -0.05, unique(pop[ord]), xpd = TRUE)
  abline(v = cumsum(sapply(unique(pop[ord]), function(x) sum(pop[ord] == x))), col = 1, lwd = 1.2)
  dev.off()
  TRUE
}

plot_admixture_q(as.integer(Sys.getenv("K", "3")), as.integer(Sys.getenv("SEED", "1")))
invisible(lapply(2:5, plot_admixture_q, seed = 1))
invisible(lapply(2:5, plot_evaladmix, seed = 1))
invisible(lapply(2:5, function(k) plot_evaladmix(k, seed = 1, suffix = "best")))

plot_loglikelihoods(file.path(figures_dir, "admixture_convergence.png"))

eig_file <- file.path(results_dir, "example.pcaone_plot.eigvecs2")
if (!file.exists(eig_file)) {
  eig_file <- file.path(results_dir, "example.pcaone.eigvecs2")
}
if (file.exists(eig_file)) {
  pcs <- read_pcaone_eigvecs2(eig_file)
  pc_cols <- grep("^PC[0-9]+$", names(pcs), value = TRUE)
  if (length(pc_cols) >= 10) {
    plot_pc_pairs(pcs, read_pop_labels(pcs), file.path(figures_dir, "pcaone_top10_pc_pairs.png"))
  }
}
