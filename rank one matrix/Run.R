remove(list = ls())
seed = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

### load parameters and generate data from aCSS code
example = "rank_one_matrix"
dir <- ifelse(requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable(),
              dirname(rstudioapi::getActiveDocumentContext()$path),
              dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1],fixed=TRUE))))
setwd(dir) # if this fails, change this to wherever this file is stored.

source(paste0(example,'_source.R'))
source(paste0("BaCSS_",example, '_source.R'))

parameters = generate_parameters()
if(parameters$sigma != 0.5)
  stop("the sampling distribution for \\tilde X is written for sigma = 0.5. If sigma is changed one would need to change functions for log_posterior and its (double) derivatives")


## Set parameters for aCSS-Bayes
## Prior on U: U_i iid ~ N(bayes_acss_parameters$prior_mean[1], bayes_acss_parameters$prior_sd[1]^2) 
## Prior on V: V_j iid ~ N(bayes_acss_parameters$prior_mean[2], bayes_acss_parameters$prior_sd[2]^2) (indep of U_i's)

bayes_acss_parameters = list(prior_mean = c(0, 0), prior_sd = c(1, 1), 
                             burnin_posterior = 500, thinning_posterior = 10, L = 1,
                             M=300, B = 25, n = parameters$example$n)

if(any(bayes_acss_parameters$prior_mean != 0))
  stop("Both prior means must be 0")

start_time <- Sys.time()
cat(paste0("\nseed = ", seed, "\n"))
pval = run_one_trial(seed = seed, example = example, print_progress = FALSE, bayes_acss_parameters = bayes_acss_parameters)
end_time <- Sys.time()
cat(paste0("\nseed = ", seed, ", time=", end_time - start_time, "\n"))

result = list("seed" = seed, "pval" = pval)

saveRDS(object = result, file = paste0("result files/Result-", "seed=", seed, ".rds"))

