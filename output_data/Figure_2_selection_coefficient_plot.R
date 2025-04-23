library(cancereffectsizeR)
library(tidyverse)
library(data.table)

# Load in cesa if necessary
# cesa <- load_cesa("analysis/eso_cesa_after_analysis.rds") 

selection_results <- rbind(cesa$selection$TP53,
                           cesa$selection$NOTCH1,
                           cesa$selection$NOTCH2,
                           cesa$selection$ERBB4,
                           cesa$selection$NFE2L2,
                           cesa$selection$PIK3CA,
                           cesa$selection$CDKN2A,
                           cesa$selection$ARID1A,
                           cesa$selection$FAT1,
                           cesa$selection$EGFR,
                           cesa$selection$ERBB2,
                           cesa$selection$FBXW7,
                           cesa$selection$FGFR3,
                           cesa$selection$RB1,
                           cesa$selection$SMAD4,
                           cesa$selection$SOX2)

selection_results <- selection_results %>%
  mutate(progression = cesa$groups) #take selection results and add progression (Pre/Pri) column

selection_results <- selection_results %>%
  mutate(gene = gsub("\\.1.*","",variant_name)) #extract gene name from variant_name


# Compute 95% CIs ----
ci_results <- list()

for (gene in genes) {
  fit_list <- attr(cesa@selection_results[[gene]], "fit")
  
  fit <- fit_list[[1]]
  
  coef_vals <- coef(fit)
  ci_vals <- tryCatch(
    confint(fit, method = "quad"),  # Wald CI
    error = function(e) matrix(rep(NA, 4), ncol = 2, byrow = TRUE,
                               dimnames = list(c("si_Pre", "si_Pri"), c("2.5 %", "97.5 %")))
  )
  
  ci_results[[gene]] <- data.frame(
    gene = gene,
    si_Pre = coef_vals["si_Pre"],
    si_Pri = coef_vals["si_Pri"],
    si_Pre_low = ci_vals["si_Pre", 1],
    si_Pre_high = ci_vals["si_Pre", 2],
    si_Pri_low = ci_vals["si_Pri", 1],
    si_Pri_high = ci_vals["si_Pri", 2]
  )
}

ci_df <- do.call(rbind, ci_results)

ci_df <- ci_df %>%
  mutate(si_Pre_low = ifelse(si_Pre_low < 0.001, 0.001, si_Pre_low)) %>% # set 0.001 as lower bound for selection
  mutate(si_Pri_low = ifelse(si_Pri_low < 0.001, 0.001, si_Pri_low))
ci_df$gene <- factor(ci_df$gene, levels = unique(ci_df$gene))

# Plot step-specific selection results ----
pre_ci_df <- ci_df %>%
  select(gene, selection_intensity = si_Pre, ci_low = si_Pre_low, ci_high = si_Pre_high) %>%
  mutate(group = "Pre")
pri_ci_df <- ci_df %>%
  select(gene, selection_intensity = si_Pri, ci_low = si_Pri_low, ci_high = si_Pri_high) %>%
  mutate(group = "Pri")

stage_data <- rbind(pre_ci_df, pri_ci_df)

stage_data <- stage_data %>%
  left_join(loglik_df, by = "gene") # add significance results from LRT
stage_data$gene <- reorder(stage_data$gene, stage_data$p_value) # order by significance

# define significance levels
stage_data <- stage_data %>%
  mutate(significance = case_when(
    p_value <= 0.001 ~ "***",
    p_value <= 0.01  ~ "**",
    p_value <= 0.05  ~ "*",
    p_value > 0.05   ~ "ns"
  ))

label_data <- stage_data %>%
  group_by(gene) %>%
  summarize(x = 1.5,
            y = max(ci_high) * 1.05,
            label = significance[1],  # same across Pre and Pri
            .groups = "drop")

selection_plots <- ggplot(stage_data, aes(x = group, y = selection_intensity, color = group)) +
  geom_point(size=2.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, width = 0.25)) +
  labs(x="Evolutionary trajectory", y="Scaled selection coefficient", color = "Tissue type") +
  scale_color_manual(labels = c("Normal", "Tumor"), values = c("#F8766D", "#00BFC4")) +
  theme_bw() +
  facet_wrap(~gene, ncol=4, scales = "free_y") +
  expand_limits(y = 0) +
  geom_hline(yintercept = 0, lwd = 0.5, color = "lightgrey", linetype = "dotted") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.ticks.x = element_blank(),
        axis.title=element_text(size = 22),
        legend.position = "none",
        text = element_text(size = 22)) +
  scale_x_discrete(labels= expression(O %->% N, N %->% T)) +
  geom_text(data = label_data,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 6) + 
  geom_text(data = label_data,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 6)

ggsave("output_data/fig_2_selection_plots.png", selection_plots, width=12.7, height=6.38)
ggsave("output_data/fig_2_selection_plots.pdf", selection_plots, width=12.7, height=6.38)
ggsave("output_data/fig_2_selection_plots.jpg", selection_plots, width=12.7, height=6.38)
