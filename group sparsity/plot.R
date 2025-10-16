remove(list = ls())

library(ggplot2)
library(ggpubr)

dir <- ifelse(requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable(),
              dirname(rstudioapi::getActiveDocumentContext()$path),
              dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1],fixed=TRUE))))
setwd(dir) # if this fails, change this to wherever this file is stored.

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
pval_oracle = matrix(pval_matrix[,"oracle"],nrow = length(signal))
pval_bacss = matrix(pval_matrix[,"Bayesian_aCSS"],nrow = length(signal))

power_bacss <- rowMeans(pval_bacss < 0.05, na.rm = TRUE)
stderr_bacss <- sqrt(power_bacss * (1 - power_bacss) / ncol(pval_bacss))
power_oracle <- rowMeans(pval_oracle < 0.05, na.rm = TRUE)
stderr_oracle <- sqrt(power_oracle * (1 - power_oracle) / ncol(pval_oracle))

plot_data <- data.frame(
  Signal = rep(as.numeric(unique(rownames(pval_matrix))), 2),
  Power = c(power_bacss, power_oracle),
  StdErr = c(stderr_bacss, stderr_oracle),
  Method = factor(rep(c("aCSS-B", "oracle"), each = length(signal)), levels = c("oracle","aCSS-B"))
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



