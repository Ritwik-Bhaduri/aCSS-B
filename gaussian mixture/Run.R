remove(list = ls())
library(Rcpp)
seed = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

### load parameters and generate data from aCSS code
example = "gaussian_mixture"
dir <- ifelse(requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable(),
              dirname(rstudioapi::getActiveDocumentContext()$path),
              dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1],fixed=TRUE))))
setwd(dir) # if this fails, change this to wherever this file is stored.

source(paste0(example,'_source.R'))
source(paste0("BaCSS_",example,'_source.R'))

parameters = generate_parameters()
nsignal = length(parameters$signal)

## Set parameters for aCSS-Bayes.
## This code is written for the following model: X ~ w_1 N(mu_1, sigma_1^2) + (1-w_1) N(mu_2, sigma_2^2) (number of mixture components is fixed at 2)
## Prior on mean:     mu_i ~ N(bayes_acss_parameters$prior_params["mu"]["mean"], sigma_i^2 / bayes_acss_parameters$prior_params["mu"]["kappa"]) i = 1,2
## Prior on variance: sigma_1^2, sigma_2^2 iid ~ Inv-Gamma(bayes_acss_parameters$prior_params["sigma2"]["alpha"], bayes_acss_parameters$prior_params["sigma2"]["beta"])
## Prior on weights:  w ~ beta(bayes_acss_parameters$prior_params["w"]["pi"])

bayes_acss_parameters = list(prior_params = list("mu" = list("mean" = 0, "kappa" = 1),
                                                 "sigma2" = list("alpha" = 1, "beta" = 0.5),
                                                 "w" = list("pi" = rep(2, parameters$J))),
                             burnin = 500, thinning = 10, B = 25, 
                             M = 300, L = 1, n = parameters$example$n)

start_time <- Sys.time()
cat(paste0("\nseed = ", seed, "\n"))
pval = run_one_trial(seed = seed, example = example, print_progress = FALSE, bayes_acss_parameters = bayes_acss_parameters)
end_time <- Sys.time()
cat(paste0("\nseed = ", seed, ", time=", end_time - start_time, "\n"))

result = list("seed" = seed, "pval" = pval)

saveRDS(object = result, file = paste0("result files/Result-", "seed=", seed, ".rds"))

