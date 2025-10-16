example = "logistic"
dir <- ifelse(requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable(),
              dirname(rstudioapi::getActiveDocumentContext()$path),
              dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1],fixed=TRUE))))
setwd(dir) # if this fails, change this to wherever this file is stored.

source(paste0("../aCSS code/", example,'_source.R'))
parameters = generate_parameters()

pval_ACSS = readRDS("files/pval_ACSS.rds")
pval_BACSS = readRDS("files/pval_BACSS.rds")
pval_ACSS = pval_ACSS["aCSS"]

pval_list = c(pval_ACSS, pval_BACSS)

power_list = lapply(pval_list, function(pval) rowMeans(pval < 0.05))
stderr_list = lapply(1:length(pval_list), function(l) sqrt(power_list[[l]] * (1 - power_list[[l]]) / ncol(pval_list[[l]])))

power_list <- power_list[!sapply(power_list, function(x) all(is.na(x)))]
stderr_list <- stderr_list[!sapply(stderr_list, function(x) all(is.na(x)))]
methods = names(power_list)
methods[methods == "Bayesian_aCSS"] = "aCSS-B"

plot_data <- data.frame(
  Signal = rep(parameters$signal, length(methods)),
  Power = do.call(c, power_list),
  StdErr = do.call(c, stderr_list),
  Method = rep(methods, each = length(parameters$signal))
)
plot_data <- data.frame(
  Signal = rep(parameters$signal, length(methods)),
  Power = do.call(c, power_list),
  StdErr = do.call(c, stderr_list),
  Method = rep(methods, each = length(parameters$signal))
)
plot_data$Method <- factor(plot_data$Method,levels = c("oracle", "aCSS-B", "aCSS"))
my_cols <- c(oracle="#0943A8",`aCSS-B`="#E74C3C", aCSS="#2CA02C")
my_ltys <- c(oracle   = "dashed",`aCSS-B` = "solid",aCSS="dotdash")

p <- ggplot(plot_data, aes(x = Signal, y = Power, colour = Method, linetype = Method)) +
  geom_hline(yintercept = 0.05, col = "black", linetype = "dotted") +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = Power - StdErr, ymax = Power + StdErr), width = 0.01, show.legend = FALSE) +
  geom_point(show.legend = FALSE) +
  scale_colour_manual(values = my_cols) +
  scale_linetype_manual(values = my_ltys) +
  labs(title = "Power Comparison", x = "\nSignal strength (c)", y = "Power\n") + 
  ylim(0,1) +
  theme_light() +
  theme(legend.position = "right", legend.title = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22), axis.title = element_text(size = 22), plot.title = element_blank())

p; ggsave("plots/Power logistic.pdf", p, width = 9, height = 6, dpi = 600)

