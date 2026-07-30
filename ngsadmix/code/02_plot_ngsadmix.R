source("https://raw.githubusercontent.com/GenisGE/evalAdmix/master/visFuns.R")

data_dir <- "../tutorial_data/ngsadmix"
results_dir <- "results"
figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

k <- Sys.getenv("K", "3")
seed <- Sys.getenv("SEED", "3")

pop <- read.table(file.path(data_dir, "1000G5pops.pop.info"), as.is = TRUE)
q <- read.table(file.path(results_dir, sprintf("1000G5pops.ngsadmix.K%s.seed%s.qopt", k, seed)))

ord <- orderInds(pop = pop[, 1], q = q)

png(file.path(figures_dir, sprintf("ngsadmix_k%s.png", k)), width = 1400, height = 500, res = 160)
par(mar = c(5, 4, 1, 1))
barplot(
  t(q)[, ord],
  col = 2:10,
  space = 0,
  border = NA,
  xlab = "Individuals",
  ylab = "Admixture proportions"
)
text(sort(tapply(seq_len(nrow(pop)), pop[ord, 1], mean)), -0.05, unique(pop[ord, 1]), xpd = TRUE)
abline(v = cumsum(sapply(unique(pop[ord, 1]), function(x) sum(pop[ord, 1] == x))), col = 1, lwd = 1.2)
dev.off()
