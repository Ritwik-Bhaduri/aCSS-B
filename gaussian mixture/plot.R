remove(list = ls())
example = "gaussian_mixture"
dir <- ifelse(requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable(),
              dirname(rstudioapi::getActiveDocumentContext()$path),
              dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1],fixed=TRUE))))
setwd(dir) # if this fails, change this to wherever this file is stored.
source(paste0(example,'_source.R'))

parameters = generate_parameters()
library(ggplot2)

# load regularized aCSS data
load("files/pval_mixturegaussian_sig8.RData")
pval_ACSS = pval_mixturegaussian_sig8
pval_oracle = do.call(cbind, lapply(pval_ACSS, function(mat) mat[,"oracle"]))
pval_ACSS = do.call(cbind, lapply(pval_ACSS, function(mat) mat[,"aCSS"]))

# load BACSS data
file_list = list.files(path="./result files/", pattern=NULL, all.files=FALSE,full.names=FALSE)
file_list <- file_list[grep("^Result-", file_list)]

result_list = list()
for(k in 1:length(file_list)){
  file = file_list[k]
  result = readRDS(paste0("result files/", file))
  result_list[[k]] = result$pval
}
pval_matrix = do.call(rbind,result_list)
signal = as.numeric(unique(rownames(pval_matrix)))
pval_BACSS = matrix(pval_matrix[,"Bayesian_aCSS"],nrow = length(signal))

pval_list = list(pval_oracle, pval_ACSS, pval_BACSS)

power_list = lapply(pval_list, function(pval) rowMeans(pval < 0.05))
stderr_list = lapply(1:length(pval_list), function(l) sqrt(power_list[[l]] * (1 - power_list[[l]]) / ncol(pval_list[[l]])))

power_list <- power_list[!sapply(power_list, function(x) all(is.na(x)))]
stderr_list <- stderr_list[!sapply(stderr_list, function(x) all(is.na(x)))]
methods = c("oracle", "reg-aCSS", "aCSS-B")

plot_data <- data.frame(
  Signal = rep(parameters$signal, times = length(methods)),
  Power = do.call(c, power_list),
  StdErr = do.call(c, stderr_list),
  Method = factor(rep(methods, each = length(parameters$signal)), 
                  levels = c("reg-aCSS", "aCSS-B", "oracle"))
)

colors = c("#2CA02C", "#E74C3C", "#0943A8")
p <- ggplot(plot_data, aes(x = Signal, y = Power)) +
  geom_line(aes(linetype = Method, color = Method), linewidth = 1) +
  geom_errorbar(aes(ymin = Power - StdErr, ymax = Power + StdErr, color = Method), width = 0.005, show.legend = FALSE) +
  geom_point(aes(color = Method), show.legend = FALSE) +
  labs(title = "Power Comparison", x = expression(atop("", pi[0])), y = "Power\n") +
  theme_light() +
  scale_color_manual(values = colors)+
  scale_linetype_manual(values = c("dotdash", "solid", "dashed"))+
  geom_hline(yintercept = 0.05, col = 'black', linetype = "dotted") + 
  theme(legend.position = "right", legend.title = element_text(size = 22), legend.text = element_text(size = 22), 
        axis.text = element_text(size = 22), axis.title = element_text(size = 22), plot.title = element_blank())

p; ggsave("plots/Power mixture of gaussians.pdf", p, width = 9, height = 6, dpi = 600)

