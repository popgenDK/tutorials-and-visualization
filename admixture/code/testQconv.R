get_fast <- function(q, q_best) {
  n_pop <- nrow(q_best)
  res <- matrix(NA_real_, nrow = nrow(q), ncol = 2)
  for (g in seq_len(nrow(q))) {
    dist <- rowSums((matrix(q[g, ], nrow = n_pop, ncol = ncol(q), byrow = TRUE) - q_best)^2)
    res[g, ] <- c(which.min(dist), min(dist))
  }

  duplicated_match <- duplicated(res[, 1])
  duplicated_ids <- res[duplicated_match, 1]
  for (id in unique(duplicated_ids)) {
    idx <- which(res[, 1] == id)
    if (length(idx) > 1) {
      res[idx[which.max(res[idx, 2])], 1] <- n_pop + 1
    }
  }
  q[order(res[, 1]), ]
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript code/testQconv.R <loglikelihood-file> <Q-list-file> <threshold>")
}

ll <- read.table(args[1], header = FALSE)[, 1]
q_list <- read.table(args[2], header = FALSE, as.is = TRUE)[, 1]
threshold <- as.numeric(args[3])

best <- which.max(ll)
q_best <- t(as.matrix(read.table(q_list[best])))

converged <- logical(length(ll))
for (i in seq_along(ll)) {
  q_current <- t(as.matrix(read.table(q_list[i])))
  q_current <- get_fast(q_current, q_best)
  converged[i] <- max(abs(q_current - q_best)) < threshold
}

cat(sum(converged))
