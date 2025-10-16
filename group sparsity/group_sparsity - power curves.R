remove(list = ls())
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
example = "group_lasso"
source(paste0(example,'_source.R'))
source(paste0("BaCSS_",example,'_source.R'))

parameters = generate_parameters()

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

## set parameters for aCSS-Bayes
## Prior: beta_{I_{g^\star}} ~ N(rep(bayes_acss_parameters$mu, n_{g^\star}), bayes_acss_parameters$tau^2 * diag(n_{g^\star})), n_{g^\star} = dimension of beta_{I_{g^\star}}

bayes_acss_parameters = list(mu = 0, tau = 1, # prior parameters
                             M=300, B = 25, n = parameters$n, L = 1)

start_time <- Sys.time()
pval_matrix <- foreach(k = 1:500, .combine = rbind, .inorder = FALSE) %dopar%{
  cat(sprintf("Running iteration k: %d, seed: %d\n", k, k))
  run_one_trial(seed= k, example = 'group_lasso', print_progress = FALSE, bayes_acss_parameters = bayes_acss_parameters)
}
parallel::stopCluster(my.cluster)
close(output_conn)
end_time <- Sys.time()
print(end_time-start_time)

pval_list = lapply(1:ncol(pval_matrix), function(i) matrix(pval_matrix[,i], nrow = length(parameters$signal), byrow = FALSE))
names(pval_list) = colnames(pval_matrix)

saveRDS(pval_list,file="files/pval_oracle_BACSS.rds")

pval_list = readRDS("files/pval_oracle_BACSS.rds")

power_list = lapply(pval_list, function(pval) rowMeans(pval < 0.05))
stderr_list = lapply(1:length(pval_list), function(l) sqrt(power_list[[l]] * (1 - power_list[[l]]) / ncol(pval_list[[l]])))

power_list <- power_list[!sapply(power_list, function(x) all(is.na(x)))]
stderr_list <- stderr_list[!sapply(stderr_list, function(x) all(is.na(x)))]
methods = names(power_list)
methods[methods == "Bayesian_aCSS"] = "aCSS-B"

plot_data <- data.frame(
  Signal = rep(parameters$signal, length(methods)),
  Power = do.call(c, lapply(power_list, function(arr) arr[1:length(parameters$signal)])),
  StdErr = do.call(c, lapply(stderr_list, function(arr) arr[1:length(parameters$signal)])),
  Method = factor(rep(methods, each = length(parameters$signal)), levels = c("oracle", "aCSS-B"))
)

colors = c("#0943A8", "#E74C3C")
p <- ggplot(plot_data, aes(x = Signal, y = Power, label = Method)) +
  geom_line(aes(linetype = Method, color = Method), linewidth = 1) +  # Differentiate lines by method
  geom_errorbar(aes(ymin = Power - StdErr, ymax = Power + StdErr, color = Method), width = 0.01, show.legend = FALSE) +
  geom_point(aes(color = Method), show.legend = FALSE) +
  labs(title = "Power Comparison", x = "\n Signal strength (c)", y = "Power \n") +
  theme_light() +
  scale_color_manual(values = colors)+
  ylim(c(0,1)) +
  scale_linetype_manual(values = c("dashed", "solid", "dotdash"))+
  geom_hline(yintercept = 0.05, col = 'black', linetype = "dotted") + 
  theme(
    legend.position = "right",  # Keep legend on the right
    legend.title = element_text(size = 22),  # Larger legend title
    legend.text = element_text(size = 22),  # Larger legend text
    axis.text = element_text(size = 22),  # Increase axis text size
    axis.title = element_text(size = 22),  # Bold axis titles
    plot.title = element_blank()#element_text(size = 18, face = "bold", hjust = 0.5)  # Center title
  )

p 
ggsave("plots/Power group sparsity.pdf", p, width = 9, height = 6, dpi = 600)



