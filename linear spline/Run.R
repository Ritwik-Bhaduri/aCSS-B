remove(list = ls())
library(Rcpp)
seed = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

### load parameters and generate data from aCSS code
example = "linear_spline"

dir <- ifelse(requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable(),
              dirname(rstudioapi::getActiveDocumentContext()$path),
              dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1],fixed=TRUE))))
setwd(dir) # if this fails, change this to wherever this file is stored.

source(paste0(example,'_source.R'))
source(paste0("BaCSS_",example,'_source.R'))

parameters = generate_parameters()
nsignal = length(parameters$signal)

## set parameters for aCSS-Bayes
## Prior on beta: Beta_j iid ~ N(bayes_acss_parameters$mu_1, bayes_acss_parameters$tau_1^2)
## Prior on knot: t ~  N(bayes_acss_parameters$mu_2, bayes_acss_parameters$tau_2^2)

if(parameters$k != 2)
  stop("Current implementation only works for one knot (k=2)")

bayes_acss_parameters <- list(
  mu_1 = 0, mu_2 = 0, tau_1 = 1, tau_2 = 1, # prior parameters
  burnin_posterior = 100,
  M = 300, B = 25, n = parameters$n, L = 1, var_min = 1e-8
)

start_time <- Sys.time()
cat(paste0("\nseed = ", seed, "\n"))
pval = run_one_trial(seed = seed, example = example, print_progress = FALSE, bayes_acss_parameters = bayes_acss_parameters)
end_time <- Sys.time()
cat(paste0("\nseed = ", seed, ", time=", end_time - start_time, "\n"))

result = list("seed" = seed, "pval" = pval)

saveRDS(object = result, file = paste0("result files/Result-", "seed=", seed, ".rds"))

