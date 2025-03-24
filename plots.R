#------------------------------------------------------------------------------#
#                     Code to plot simulation results
#------------------------------------------------------------------------------#

library(tidyverse)
library(ggh4x)

#------------------------#
# Load data
#------------------------#

df <- readRDS(file.path("files", "results_all.RDS"))

#------------------------#
# Model recovery plot
#------------------------#

df_recov <- df[,c("n", "k", "s", "mods", "mcc", "Sensitivity", "Specificity", "label")]

df_recov_long <- pivot_longer(df_recov,
                              cols = c("mcc", "Sensitivity", "Specificity"),
                              names_to = "metrics",
                              values_to = "recovery")

df_recov_long$n <- as.factor(df_recov_long$n)
old_names <- c(50, 100)
new_names <- c("Time = 50", "Time = 100")
df_recov_long$n <- plyr::mapvalues(
  df_recov_long$n, 
  from = old_names, 
  to = new_names
)
df_recov_long$n <- factor(df_recov_long$n, levels = c("Time = 50", "Time = 100"))

df_recov_long$k <- as.factor(df_recov_long$k)
old_names <- c(50, 100)
new_names <- c("Individuals = 50", "Individuals = 100")
df_recov_long$k <- plyr::mapvalues(
  df_recov_long$k, 
  from = old_names, 
  to = new_names
)
df_recov_long$k <- factor(df_recov_long$k, levels = c("Individuals = 50", "Individuals = 100"))

df_recov_long$s <- as.factor(df_recov_long$s)
old_names <- c(2, 3)
new_names <- c("Subgroups = 2", "Subgroups = 3")
df_recov_long$s <- plyr::mapvalues(
  df_recov_long$s, 
  from = old_names, 
  to = new_names
)
df_recov_long$s <- factor(df_recov_long$s, levels = c("Subgroups = 2", "Subgroups = 3"))

df_recov_long$label <- as.factor(df_recov_long$label)
old_names <- c("bal", "unbal")
new_names <- c("Balanced Subgroups", "Unbalanced Subgroups")
df_recov_long$label <- plyr::mapvalues(
  df_recov_long$label, 
  from = old_names, 
  to = new_names
)
df_recov_long$label <- factor(df_recov_long$label, levels = c("Balanced Subgroups", "Unbalanced Subgroups"))

df_recov_long$mods <- as.factor(df_recov_long$mods)
old_names <- c("conf", "nosub", "sgimme", "walkunique", "sgcvar")
new_names <- c("Confirmatory Multi-VAR", "Standard Multi-VAR", "S-GIMME-VAR", "Subgrouping Multi-VAR", "scGVAR")
df_recov_long$mods <- plyr::mapvalues(
  df_recov_long$mods, 
  from = old_names, 
  to = new_names
)
df_recov_long$mods <- factor(df_recov_long$mods, 
                             levels = c("Subgrouping Multi-VAR", "Confirmatory Multi-VAR", "Standard Multi-VAR", "S-GIMME-VAR", "scGVAR"))

df_recov_long$metrics <- as.factor(df_recov_long$metrics)
old_names <- c("mcc", "Sensitivity", "Specificity")
new_names <- c("MCC", "Sensitivity", "Specificity")
df_recov_long$metrics <- plyr::mapvalues(
  df_recov_long$metrics, 
  from = old_names, 
  to = new_names
)
df_recov_long$metrics <- factor(df_recov_long$metrics, levels = c("MCC", "Sensitivity", "Specificity"))


ggplot(df_recov_long, aes(x = metrics, y = recovery, fill = mods)) +
  geom_boxplot(varwidth = F, lwd=.1, outlier.size = .1, outlier.shape=NA, fatten = 5) +
  facet_nested(label + s ~ k + n) + 
  theme_bw() +
  xlab("Metrics") +
  ylab("Model Recovery") +
  theme(legend.position = "top") +
  scale_fill_grey(labels = c("Subgrouping Multi-VAR", "Confirmatory Multi-VAR", "Standard Multi-VAR", "S-GIMME-VAR", "scGVAR"),  
                  start = 0.5, end = 0.9) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.title = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size = 8)) +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(color = "black", face = "bold")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", linewidth = .2) +
  geom_hline(yintercept = .25, linetype = "dashed", color = "gray", linewidth = .2) + 
  geom_hline(yintercept = .50, linetype = "dashed", color = "gray", linewidth = .2) + 
  geom_hline(yintercept = .75, linetype = "dashed", color = "gray", linewidth = .2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray", linewidth = .2) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 12, angle = 90, hjust = .5, vjust = .5) ) +
  theme(axis.title.x = element_text(color = "black", size = 12, angle = 0, hjust = .5, vjust = 0) ) +
  theme(axis.text.y = element_text(color = "black", size = 8, angle = 0, hjust = 1, vjust = 0) ) +
  theme(axis.text.x = element_text(color = "black", size = 8, angle = 90, hjust = .5, vjust = .5) ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

ggsave(
  "model_recovery.pdf",
  plot = last_plot(),
  width = 7,
  height =9,
  units = c("in"),
  dpi = 300,
  path = "plots"
)

#---------------------------#
# Parameter estimate plot
#---------------------------#

df_params <- df[,c("n", "k", "s", "mods", "biasabs", "rmse", "label")]

df_params_long <- pivot_longer(df_params,
                               cols = c("biasabs", "rmse"),
                               names_to = "metrics",
                               values_to = "param_values")

df_params_long$n <- as.factor(df_params_long$n)
old_names <- c(50, 100)
new_names <- c("Time = 50", "Time = 100")
df_params_long$n <- plyr::mapvalues(
  df_params_long$n, 
  from = old_names, 
  to = new_names
)
df_params_long$n <- factor(df_params_long$n, levels = c("Time = 50", "Time = 100"))

df_params_long$k <- as.factor(df_params_long$k)
old_names <- c(50, 100)
new_names <- c("Individuals = 50", "Individuals = 100")
df_params_long$k <- plyr::mapvalues(
  df_params_long$k, 
  from = old_names, 
  to = new_names
)
df_params_long$k <- factor(df_params_long$k, levels = c("Individuals = 50", "Individuals = 100"))

df_params_long$s <- as.factor(df_params_long$s)
old_names <- c(2, 3)
new_names <- c("Subgroups = 2", "Subgroups = 3")
df_params_long$s <- plyr::mapvalues(
  df_params_long$s, 
  from = old_names, 
  to = new_names
)
df_params_long$s <- factor(df_params_long$s, levels = c("Subgroups = 2", "Subgroups = 3"))

df_params_long$label <- as.factor(df_params_long$label)
old_names <- c("bal", "unbal")
new_names <- c("Balanced Subgroups", "Unbalanced Subgroups")
df_params_long$label <- plyr::mapvalues(
  df_params_long$label, 
  from = old_names, 
  to = new_names
)
df_params_long$label <- factor(df_params_long$label, levels = c("Balanced Subgroups", "Unbalanced Subgroups"))

df_params_long$mods <- as.factor(df_params_long$mods)
old_names <- c("conf", "nosub", "sgimme", "walkunique", "sgcvar")
new_names <- c("Confirmatory Multi-VAR", "Standard Multi-VAR", "S-GIMME-VAR", "Subgrouping Multi-VAR", "scGVAR")
df_params_long$mods <- plyr::mapvalues(
  df_params_long$mods, 
  from = old_names, 
  to = new_names
)
df_params_long$mods <- factor(df_params_long$mods, 
                              levels = c("Subgrouping Multi-VAR", "Confirmatory Multi-VAR", "Standard Multi-VAR", "S-GIMME-VAR", "scGVAR"))

df_params_long$metrics <- as.factor(df_params_long$metrics)
old_names <- c("biasabs", "rmse")
new_names <- c("Absolute Bias", "RMSE")
df_params_long$metrics <- plyr::mapvalues(
  df_params_long$metrics, 
  from = old_names, 
  to = new_names
)
df_params_long$metrics <- factor(df_params_long$metrics, levels = c("Absolute Bias", "RMSE"))

ggplot(df_params_long, aes(x = metrics, y = param_values, fill = mods)) +
  geom_boxplot(varwidth = F, lwd = .1, outlier.size = .1, outlier.shape = NA, fatten = 5) +
  facet_nested(label + s ~ k + n) + 
  theme_bw() +
  xlab("Metrics") +
  ylab("Values") +
  theme(legend.position="top") +
  scale_fill_grey(labels = c("Subgrouping Multi-VAR", "Confirmatory Multi-VAR", "Standard Multi-VAR", "S-GIMME-VAR","scGVAR"),  
                  start = 0.5, end = 0.9) +
  coord_cartesian(ylim = c(0, 0.2)) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.title = element_text(size = 12,face = "bold")) +
  theme(legend.text = element_text(size = 8)) +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(color = "black", face = "bold")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", linewidth = .2) +
  geom_hline(yintercept = .05, linetype = "dashed", color = "gray", linewidth = .2) + 
  geom_hline(yintercept = .10, linetype = "dashed", color = "gray", linewidth = .2) + 
  geom_hline(yintercept = .15, linetype = "dashed", color = "gray", linewidth = .2) +
  geom_hline(yintercept = .2, linetype = "dashed", color = "gray", linewidth = .2) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 12, angle = 90, hjust = .5, vjust = .5) ) +
  theme(axis.title.x = element_text(color = "black", size = 12, angle = 0, hjust = .5, vjust = 0) ) +
  theme(axis.text.y = element_text(color = "black", size = 8, angle = 0, hjust = 1, vjust = 0) ) +
  theme(axis.text.x = element_text(color = "black", size = 8, angle = 90, hjust = .5, vjust = .5) ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

ggsave(
  "param_quality.pdf",
  plot = last_plot(),
  width = 7,
  height =9,
  units = c("in"),
  dpi = 300,
  path = "plots"
)

