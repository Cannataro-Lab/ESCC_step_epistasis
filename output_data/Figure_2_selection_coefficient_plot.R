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

pre <- selection_results[, .(variant = variant_name, ci_low_95_si_Pre, ci_high_95_si_Pre, variant_type, gene, si_Pre, progression = "Pre")] #filter normal samples
pri <- selection_results[, .(variant = variant_name, ci_low_95_si_Pri, ci_high_95_si_Pri, variant_type, gene, si_Pri, progression = "Pri")] #filter tumor samples

# Rename columns and combine data frames
setnames(pre, c("si_Pre", "ci_low_95_si_Pre", "ci_high_95_si_Pre"),
         c("selection_intensity", "ci_low_95", "ci_high_95"))
setnames(pri, c("si_Pri", "ci_low_95_si_Pri", "ci_high_95_si_Pri"),
         c("selection_intensity", "ci_low_95", "ci_high_95"))

twostage_results <- rbind(pre, pri) 


# Restrict variants recurrent within each progression
twostage_results_normal <- twostage_results[progression == "Pre"]
twostage_results_primary <- twostage_results[progression == "Pri"]


stage_data <- data.frame(variant = twostage_results$variant,
                         gene = twostage_results$gene,
                         selection_intensity = twostage_results$selection_intensity,
                         group = twostage_results$progression,
                         ci_low = twostage_results$ci_low_95,
                         ci_high = twostage_results$ci_high_95)

stage_data <- stage_data %>%
  mutate(ci_low_95 = if_else(is.na(ci_low), 0, ci_low)) #set lower bound for ci_low at 0 because it can't be negative

stage_data <- stage_data %>%
  filter(gene %in% c("TP53", 
                     "NOTCH1", 
                     "NOTCH2", 
                     "ERBB4", 
                     "NFE2L2", 
                     "PIK3CA", 
                     "CDKN2A.p16INK4a",
                     "CDKN2A.p14arf",
                     "CDKN2A",
                     "ARID1A", 
                     "FAT1", 
                     "EGFR",
                     "ERBB2",
                     "FBXW7",
                     "FGFR3",
                     "RB1",
                     "SMAD4",
                     "SOX2")) #select genes of interest (genes with clear trends in this case)

stage_data$gene <- factor(stage_data$gene, levels = unique(stage_data$gene))
stage_data$gene <- reorder(stage_data$gene, -stage_data$selection_intensity)

selection_plots <- ggplot(stage_data, aes(x = group, y = selection_intensity, color = group)) +
  geom_point(size=2.5) +
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high, width = 0.25)) +
  labs(x="Evolutionary trajectory", y="Scaled selection coefficient", color = "Tissue type") +
  scale_color_manual(labels = c("Normal", "Tumor"), values = c("#F8766D", "#00BFC4")) +
  theme_bw() +
  facet_wrap(~gene, ncol=4, scales = "free_y") +
  expand_limits(y = 0) +
  geom_vline(xintercept = 1.5, lwd = 0.5, color = "lightgrey") +
  geom_hline(yintercept = 0, lwd = 0.5, color = "lightgrey", linetype = "dotted") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.ticks.x = element_blank(),
        axis.title=element_text(size = 22),
        legend.position = "none",
        text = element_text(size = 22)) +
  scale_x_discrete(labels= expression(O %->% N, N %->% T))
selection_plots

ggsave("output_data/fig_2_selection_plots.png", selection_plots, width=12.7, height=6.38)
ggsave("output_data/fig_2_selection_plots.pdf", selection_plots, width=12.7, height=6.38)
ggsave("output_data/fig_2_selection_plots.jpg", selection_plots, width=12.7, height=6.38)
