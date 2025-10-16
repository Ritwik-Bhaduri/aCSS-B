remove(list = ls())
seed = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

### load parameters and generate data from aCSS code
example = "group_sparsity"
dir <- ifelse(requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable(),
              dirname(rstudioapi::getActiveDocumentContext()$path),
              dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1],fixed=TRUE))))
setwd(dir) # if this fails, change this to wherever this file is stored.

source(paste0(example,'_source.R'))
source(paste0("BaCSS_",example, '_source.R'))

## set parameters for aCSS-Bayes
## Prior: beta_{I_{g^\star}} ~ N(rep(bayes_acss_parameters$mu, n_{g^\star}), bayes_acss_parameters$tau^2 * diag(n_{g^\star})), n_{g^\star} = dimension of beta_{I_{g^\star}}

parameters = generate_parameters()
bayes_acss_parameters = list(mu = 0, tau = 1, # prior parameters
                             M=300, B = 25, n = parameters$n, L = 1)

start_time <- Sys.time()
cat(paste0("\nseed = ", seed, "\n"))
pval = run_one_trial(seed = seed, example = example, print_progress = FALSE, bayes_acss_parameters = bayes_acss_parameters)
end_time <- Sys.time()
cat(paste0("\nseed = ", seed, ", time=", end_time - start_time, "\n"))

result = list("seed" = seed, "pval" = pval)

saveRDS(object = result, file = paste0("result files/Result-", "seed=", seed, ".rds"))

