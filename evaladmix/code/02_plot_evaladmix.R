source("https://raw.githubusercontent.com/GenisGE/evalAdmix/master/visFuns.R")

data_dir <- "../tutorial_data/ngsadmix"
ngsadmix_results_dir <- "../ngsadmix/results"
results_dir <- "results"
figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

k <- Sys.getenv("K", "3")
seed <- Sys.getenv("SEED", "3")

pop <- read.table(file.path(data_dir, "1000G5pops.pop.info"), as.is = TRUE)
q <- read.table(file.path(ngsadmix_results_dir, sprintf("1000G5pops.ngsadmix.K%s.seed%s.qopt", k, seed)))
r <- as.matrix(read.table(file.path(results_dir, sprintf("1000G5pops.K%s.seed%s.corres", k, seed))))

ord <- orderInds(pop = pop[, 1], q = q)

png(file.path(figures_dir, sprintf("evaladmix_k%s_residuals.png", k)), width = 900, height = 800, res = 150)
plotCorRes(r, pop = pop[, 1], ord = ord, max_z = 0.2)
dev.off()
