remove(list = ls())

library(ggplot2)
#library(ggpubr)
library(dplyr)

dir <- ifelse(requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable(),
              dirname(rstudioapi::getActiveDocumentContext()$path),
              dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1],fixed=TRUE))))
setwd(dir) # if this fails, change this to wherever this file is stored.
example = gsub(" ", "_", tail(strsplit(dir, "- ")[[1]], 1))

file_list = list.files(path="./result_files/", pattern=NULL, all.files=FALSE,full.names=FALSE)

result_list = list()
for(k in 1:length(file_list)){
  file = file_list[k]
  result_list[[k]] = readRDS(paste0("result_files/", file))
}
pval_mat = do.call(rbind,result_list)
pval_df = data.frame(pval_mat)
pval_df$Signal = as.numeric(rownames(pval_mat))
rownames(pval_df) = NULL

pval_df_oracle = pval_df[pval_df[,"tau"] == pval_df[1,"tau"], ] %>% 
  group_by(Signal) %>% 
  summarize(
    n = dplyr::n(),
    Power = mean(oracle < 0.05, na.rm = TRUE),
    StdErr = sqrt(Power * (1 - Power) / n)
  )
pval_df_oracle$Method = "oracle"

pval_df_acss_b = pval_df %>% 
  group_by(tau, Signal) %>% 
  summarize(
    n = dplyr::n(),
    Power = mean(Bayesian_aCSS < 0.05, na.rm = TRUE),
    StdErr = sqrt(Power * (1 - Power) / n),
    .groups = "drop"
  )
pval_df_acss_b$Method = sapply(pval_df_acss_b$tau, function(tau) paste0("aCSS-B (tau = ", tau, ")"))

plot_data = rbind(pval_df_acss_b[,c("Signal","Power","StdErr","Method")], pval_df_oracle[,c("Signal","Power","StdErr","Method")])
plot_data$Method = factor(plot_data$Method, 
                          levels = c("oracle", sapply(sort(unique(pval_df$tau)), function(tau) paste0("aCSS-B (tau = ", tau, ")"))))

levs <- levels(plot_data$Method)

# create plotmath labels: Oracle stays plain; aCSS-B shows tau as Greek
legend_labs <- c(
  "oracle" = "oracle",
  setNames(
    lapply(sort(unique(pval_df_acss_b$tau)), function(tau) bquote("aCSS-B ("*tau*" = "*.(tau)*")")),
    sapply(sort(unique(pval_df_acss_b$tau)), function(tau) paste0("aCSS-B (tau = ", tau, ")"))
  )
)

oracle_color <- "#000000"
acss_colors <- RColorBrewer::brewer.pal(5, "Set1")[c(2,3,1,4,5)]   # or "Dark2", "Set2"
colors <- c(oracle_color, acss_colors)


p <- ggplot(plot_data, aes(x = Signal, y = Power)) +
  geom_line(aes(linetype = Method, color = Method), linewidth = 1) +
  geom_errorbar(aes(ymin = Power - StdErr, ymax = Power + StdErr, color = Method),
                width = 0.01, show.legend = FALSE) +
  geom_point(aes(color = Method), show.legend = FALSE) +
  labs(x = "\n Signal strength (c)", y = "Power \n") +
  theme_light() +
  scale_color_manual(values = colors, labels = legend_labs) +
  scale_linetype_manual(values = c("dashed", rep("solid", 5)), labels = legend_labs) +
  ylim(c(0,1)) +
  geom_hline(yintercept = 0.05, col = "black", linetype = "dotted") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 22),
    plot.title = element_blank()
  )

p

ggsave(paste0("../../plot_logistic_B_25.pdf"), p, width = 9, height = 6, dpi = 600)
