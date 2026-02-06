remove(list = ls())

library(ggplot2)
#library(ggpubr)
library(dplyr)

dir <- ifelse(requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable(),
              dirname(rstudioapi::getActiveDocumentContext()$path),
              dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1],fixed=TRUE))))
setwd(dir) # if this fails, change this to wherever this file is stored.
example = gsub(" ", "_", tail(strsplit(dir, "- ")[[1]], 1))

file_list = list.files(path="./result files_varying_B/", pattern=NULL, all.files=FALSE,full.names=FALSE)

result_list = list()
for(k in 1:length(file_list)){
  file = file_list[k]
  result_list[[k]] = readRDS(paste0("result files_varying_B/", file))
}
pval_mat = do.call(rbind,result_list)
pval_df = data.frame(pval_mat)
pval_df$Signal = as.numeric(rownames(pval_mat))
rownames(pval_df) = NULL

library(dplyr)
library(ggplot2)
library(RColorBrewer)

alpha <- 0.05

# -------------------------
# Oracle: does not depend on tau or B
# Fix any tau and any B, then compute power vs Signal
# -------------------------
tau0 <- pval_df$tau[1]
B0   <- pval_df$B[1]

pval_df_oracle <- pval_df %>%
  filter(tau == tau0, B == B0) %>%
  group_by(Signal) %>%
  summarise(
    n = dplyr::n(),
    Power = mean(oracle < alpha, na.rm = TRUE),
    StdErr = sqrt(Power * (1 - Power) / n),
    .groups = "drop"
  ) %>%
  mutate(Method = "oracle")

# We'll replicate oracle curve across taus so it appears in every facet
taus <- sort(unique(pval_df$tau))
pval_df_oracle_rep <- bind_rows(lapply(taus, function(t) mutate(pval_df_oracle, tau = t)))

# -------------------------
# Bayesian aCSS: depends on tau and B
# Compute power vs Signal for each (tau, B)
# -------------------------
pval_df_bacss <- pval_df %>%
  group_by(tau, B, Signal) %>%
  summarise(
    n = dplyr::n(),
    Power = mean(Bayesian_aCSS < alpha, na.rm = TRUE),
    StdErr = sqrt(Power * (1 - Power) / n),
    .groups = "drop"
  ) %>%
  mutate(Method = paste0("aCSS-B (B = ", B, ")"))

# -------------------------
# Combine for plotting
# -------------------------
plot_data <- bind_rows(
  pval_df_bacss %>% select(tau, Signal, Power, StdErr, Method),
  pval_df_oracle_rep %>% select(tau, Signal, Power, StdErr, Method)
)

# Method ordering: oracle first, then B levels (sorted)
B_levels <- sort(unique(pval_df$B))
method_levels <- c("oracle", paste0("aCSS-B (B = ", B_levels, ")"))
plot_data$Method <- factor(plot_data$Method, levels = method_levels)

# Legend labels (math-style B; oracle plain)
legend_labels <- setNames(
  lapply(method_levels, function(m) {
    if (m == "oracle") {
      expression(oracle)
    } else {
      B_val <- as.numeric(sub(".*B = ([^)]*).*", "\\1", m))
      bquote("aCSS-B ("*B*" = "*.(B_val)*")")
    }
  }),
  method_levels
)

# Colors: oracle black, then one per B
oracle_color <- "#000000"
acss_colors <- brewer.pal(max(3, min(8, length(B_levels))), "Set1")
if (length(B_levels) > length(acss_colors)) {
  acss_colors <- colorRampPalette(acss_colors)(length(B_levels))
} else {
  acss_colors <- acss_colors[seq_len(length(B_levels))]
}
colors <- c(oracle_color, acss_colors)

# Linetypes: oracle dashed, bacss solid
linetypes <- c("dashed", rep("solid", length(B_levels)))

# -------------------------
# Plot (similar aesthetics, facet over tau)
# -------------------------
p <- ggplot(plot_data, aes(x = Signal, y = Power)) +
  geom_line(aes(linetype = Method, color = Method), linewidth = 1) +
  geom_errorbar(
    aes(ymin = Power - StdErr, ymax = Power + StdErr, color = Method),
    width = 0.005, show.legend = FALSE
  ) +
  geom_point(aes(color = Method), show.legend = FALSE) +
  facet_wrap(~ tau, labeller = label_bquote(tau == .(tau))) +
  labs(x = "\n Signal strength (c)", y = "Power \n") +
  theme_light() +
  scale_color_manual(values = colors, labels = legend_labels) +
  scale_linetype_manual(values = linetypes, labels = legend_labels) +
  ylim(c(0, 1)) +
  geom_hline(yintercept = 0.05, col = "black", linetype = "dotted") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 22),
    plot.title = element_blank(),
    strip.text = element_text(size = 22)
  )

p


ggsave(paste0("plots/Sensitivity_group_sparsity_vs_B.pdf"), p, width = 15, height = 6, dpi = 600)
