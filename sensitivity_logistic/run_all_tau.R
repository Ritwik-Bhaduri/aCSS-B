remove(list = ls())
library('sigmoid')
library(ggplot2)
library(dplyr)
library(tidyr)
#library(parallel)
library('pbapply')
#library('ggpubr')

dir <- ifelse(requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable(),
              dirname(rstudioapi::getActiveDocumentContext()$path),
              dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1],fixed=TRUE))))
setwd(dir) # if this fails, change this to wherever this file is stored.

### load parameters and generate data from aCSS code
example = "logistic"
source(paste0("../aCSS code/", example,'_source.R'))
source("BaCSS_logistic_source.R")
parameters = generate_parameters()

## Set parameters for aCSS-Bayes.
## This code is written for prior: beta ~ N(0, tau^2)
bayes_acss_parameters = list(tau = rep(1, parameters$d), burnin = 500, thinning = 10,
                             B = 25, M = 300, L = 1, n = parameters$example$n)

##################### Parallel computing ###############################################################
start_time <- Sys.time()

parameter_grid <- expand.grid(k = c(1:500),tau = 10^((-2:2)))
chunk_size = 50
n_workers = ceiling(nrow(parameter_grid)/chunk_size)
print(n_workers)
task_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
current_idx = (task_id - 1)*chunk_size + 1:chunk_size
current_idx = current_idx[current_idx<=nrow(parameter_grid)]
df_list = list()
for (i in current_idx){
  cat(sprintf("Running iteration i: %d, seed: %d\n", i, parameter_grid[i,"k"]))
  bayes_acss_parameters_tmp <- bayes_acss_parameters
  bayes_acss_parameters_tmp$tau <- rep(parameter_grid[i,"tau"],parameters$d)
  df = run_one_trial(M=bayes_acss_parameters$M, seed= parameter_grid[i,"k"], example = 'logistic', print_progress = FALSE, parameters = parameters, bayes_acss_parameters = bayes_acss_parameters_tmp)
  df = cbind(parameter_grid[i,"k"],parameter_grid[i,"tau"],df)
  colnames(df)[1:2] = c("k","tau")
  df_list[[i]] = df
}
df = do.call(rbind,df_list)
end_time <- Sys.time()
print(end_time-start_time)
saveRDS(df, file = paste0("result_files/pvalues_aCSS_B_",task_id,".rds"))
