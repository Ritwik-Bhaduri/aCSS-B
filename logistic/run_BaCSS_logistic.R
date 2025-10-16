remove(list = ls())
library('sigmoid')
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)
library('pbapply')
library('ggpubr')

dir <- ifelse(requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable(),
              dirname(rstudioapi::getActiveDocumentContext()$path),
              dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1],fixed=TRUE))))
setwd(dir) # if this fails, change this to wherever this file is stored.

### load parameters and generate data from aCSS code
example = "logistic"
source(paste0("../aCSS code/", example,'_source.R'))
parameters = generate_parameters()
parameters$example$n = 100
parameters$theta0=rep(0.2, 5)
nsignal = length(parameters$signal)
isignal = 6 # this specifies signal=0.5 (signal is basically c which is multiplied to beta_0 and beta_1 in end of page 24 in aCSS paper)
experiment = generate_experiment(isignal,parameters)

## Set parameters for aCSS-Bayes.
## This code is written for prior: beta ~ N(0, tau^2)
bayes_acss_parameters = list(tau = rep(1, parameters$d), burnin = 500, thinning = 10,
                             B = 25, M = 300, L = 1, n = parameters$example$n)

##################### Parallel computing ###############################################################
lapply(c("pbapply", "parallelly", "foreach", "doParallel"), require, character.only = TRUE)
n.cores <- parallel::detectCores()
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "FORK",
  outfile=paste0("parallel_output.txt")
)
doParallel::registerDoParallel(cl = my.cluster)
output_conn <- file("parallel_output.txt", open = "wt")

start_time <- Sys.time()
pval_matrix <- foreach(k = 1:500, .combine = rbind, .inorder = FALSE) %dopar%{
  cat(sprintf("Running iteration k: %d, seed: %d\n", k, k))
  run_one_trial(M=bayes_acss_parameters$M, seed= k, example = 'logistic', print_progress = FALSE, parameters = parameters, bayes_acss_parameters = bayes_acss_parameters)
}
parallel::stopCluster(my.cluster)
close(output_conn)
end_time <- Sys.time()
print(end_time-start_time)

pval_list = lapply(1:ncol(pval_matrix), function(i) matrix(pval_matrix[,i], nrow = length(parameters$signal), byrow = FALSE))
names(pval_list) = colnames(pval_matrix)
saveRDS(pval_list,file="files/pval_BACSS.rds")

saveRDS(end_time-start_time,file="files/time_BACSS.rds")
