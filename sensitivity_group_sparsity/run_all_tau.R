remove(list = ls())
task_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

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
bayes_acss_parameters = list(mu = 5, tau = 1, # prior parameters # mu changed from 0 to 5
                             M=300, B = 25, n = parameters$n, L = 1)

## sensitivity analysis
parameter_grid <- expand.grid(k = c(1:500),tau = 10^(-2:2))
chunk_size = 10
n_workers = ceiling(nrow(parameter_grid)/chunk_size)
cat(paste0("number of workers needed: ", n_workers, "\n"))

current_idx = (task_id - 1)*chunk_size + 1:chunk_size
current_idx = current_idx[current_idx<=nrow(parameter_grid)]
df_list = list()

start_time <- Sys.time()
cat(paste0("Running indices ", current_idx, "\n"))
for(i in current_idx){
  k = parameter_grid[i,"k"] # seed
  tau = parameter_grid[i, "tau"]
  cat(paste0("Running iteration i: ", i, ", seed: ", k, ", tau: ", tau, "\n"))
  bayes_acss_parameters_temp = bayes_acss_parameters
  bayes_acss_parameters_temp$tau <- tau
  df = run_one_trial(seed = k, example = example, print_progress = FALSE, bayes_acss_parameters = bayes_acss_parameters_temp)
  df = cbind(parameter_grid[i,"k"], parameter_grid[i,"tau"], df)
  colnames(df)[1:2] = c("k","tau")
  df_list[[i]] = df
}
df = do.call(rbind,df_list)
end_time <- Sys.time()
print(end_time-start_time)

saveRDS(object = df, file = paste0("result files/pvalues_aCSS_B_", task_id, ".rds"))
