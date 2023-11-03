library(cancereffectsizeR)
library(tidyverse)
library(patchwork)
library(cowplot)

# Load in  cesa if necessary
# cesa <- load_cesa("analysis/eso_cesa_after_analysis.rds")

gene_ep_results <- cesa$epistasis$epistasis_compound_variants_all_samples

plot_epistasis_results_with_ci = function(ep_results, current_gene, upper_limit1, upper_limit2) {
  data_for_compare <- ep_results %>%
    mutate(variant1 = gsub("\\.1.*", "", variant_A)) %>%
    mutate(variant2 = gsub("\\.1.*", "", variant_B)) %>%
    filter(variant1 == current_gene | variant2 == current_gene) #filter for GOI results
  
  data_for_compare <- data_for_compare %>%
    mutate(gene_of_interest = current_gene) %>%
    mutate(other_gene = case_when(
      variant1 == current_gene ~ variant2,
      variant2 == current_gene ~ variant1)) %>%
    mutate(ces_GOI = case_when(
      variant1 == current_gene ~ ces_A0,
      variant2 == current_gene ~ ces_B0)) %>%
    mutate(ces_OTHER = case_when(
      variant2 == current_gene ~ ces_A0,
      variant1 == current_gene ~ ces_B0)) %>%
    mutate(ces_GOI_after_OTHER = case_when(
      variant1 == current_gene ~ ces_A_on_B,
      variant2 == current_gene ~ ces_B_on_A)) %>%
    mutate(ces_OTHER_after_GOI = case_when(
      variant2 == current_gene ~ ces_A_on_B,
      variant1 == current_gene ~ ces_B_on_A)) %>%
    mutate(joint_cov_samples_just_GOI = case_when(
      variant2 == current_gene ~ nB0,
      variant1 == current_gene ~ nA0)) %>%
    mutate(joint_cov_samples_just_OTHER = case_when(
      variant1 == current_gene ~ nB0,
      variant2 == current_gene ~ nA0)) %>%
    mutate(ci_low_95_ces_GOI = case_when(
      variant1 == current_gene ~ ci_low_95_ces_A0,
      variant2 == current_gene ~ ci_low_95_ces_B0)) %>%
    mutate(ci_high_95_ces_GOI = case_when(
      variant1 == current_gene ~ ci_high_95_ces_A0,
      variant2 == current_gene ~ ci_high_95_ces_B0)) %>%
    mutate(ci_low_95_ces_OTHER = case_when(
      variant1 == current_gene ~ ci_low_95_ces_B0,
      variant2 == current_gene ~ ci_low_95_ces_A0)) %>%
    mutate(ci_high_95_ces_OTHER = case_when(
      variant1 == current_gene ~ ci_high_95_ces_B0,
      variant2 == current_gene ~ ci_high_95_ces_A0)) %>%
    mutate(ci_low_95_ces_GOI_after_OTHER = case_when(
      variant1 == current_gene ~ ci_low_95_ces_A_on_B,
      variant2 == current_gene ~ ci_low_95_ces_B_on_A)) %>%
    mutate(ci_high_95_ces_GOI_after_OTHER = case_when(
      variant1 == current_gene ~ ci_high_95_ces_A_on_B,
      variant2 == current_gene ~ ci_high_95_ces_B_on_A)) %>%
    mutate(ci_low_95_ces_OTHER_after_GOI = case_when(
      variant1 == current_gene ~ ci_low_95_ces_B_on_A,
      variant2 == current_gene ~ ci_low_95_ces_A_on_B)) %>%
    mutate(ci_high_95_ces_OTHER_after_GOI = case_when(
      variant1 == current_gene ~ ci_high_95_ces_B_on_A,
      variant2 == current_gene ~ ci_high_95_ces_A_on_B)) %>%
    # select(gene_of_interest, other_gene, ends_with("OTHER"),ends_with("GOI"), joint_cov_samples_just_GOI, joint_cov_samples_just_OTHER, 
    #        nAB, n00, ci_low_95_ces_GOI, ci_high_95_ces_GOI, ci_low_95_ces_OTHER, ci_high_95_ces_OTHER, 
    #        ci_low_95_ces_A_on_B, ci_high_95_ces_A_on_B, ci_low_95_ces_B_on_A, ci_high_95_ces_B_on_A, p_epistasis) %>%
    mutate(ci_low_95_ces_GOI = if_else(is.na(ci_low_95_ces_GOI), 0, ci_low_95_ces_GOI)) %>%
    mutate(ci_low_95_ces_OTHER = if_else(is.na(ci_low_95_ces_OTHER), 0, ci_low_95_ces_OTHER)) %>%
    mutate(ci_low_95_ces_GOI_after_OTHER = if_else(is.na(ci_low_95_ces_GOI_after_OTHER), 0, ci_low_95_ces_GOI_after_OTHER)) %>%
    mutate(ci_low_95_ces_OTHER_after_GOI = if_else(is.na(ci_low_95_ces_OTHER_after_GOI), 0, ci_low_95_ces_OTHER_after_GOI)) %>%
    mutate(ci_high_95_ces_GOI = if_else(is.na(ci_high_95_ces_GOI), 50000, ci_high_95_ces_GOI)) %>%
    mutate(ci_high_95_ces_OTHER = if_else(is.na(ci_high_95_ces_OTHER), 50000, ci_high_95_ces_OTHER)) %>%
    mutate(ci_high_95_ces_GOI_after_OTHER = if_else(is.na(ci_high_95_ces_GOI_after_OTHER), 50000, ci_high_95_ces_GOI_after_OTHER)) %>%
    mutate(ci_high_95_ces_OTHER_after_GOI = if_else(is.na(ci_high_95_ces_OTHER_after_GOI), 50000, ci_high_95_ces_OTHER_after_GOI)) %>%
    mutate(p_epistasis_scientific = formatC(p_epistasis, format = "e", digits = 2)) %>%
    mutate(significance = case_when(p_epistasis <= 0.001 ~ "***", 
                                    p_epistasis <= 0.01 ~ "**",
                                    p_epistasis <= 0.05 ~ "*", 
                                    p_epistasis > 0.05 ~ "ns")) %>%
    mutate(change_GOI = abs(ces_GOI - ces_GOI_after_OTHER)) %>%
    mutate(change_other = abs(ces_OTHER - ces_OTHER_after_GOI))
  
  
  data_for_compare$other_gene <- factor(data_for_compare$other_gene, levels = unique(data_for_compare$other_gene))
  data_for_compare$other_gene <- reorder(data_for_compare$other_gene, -data_for_compare$change_other)
  
  
  text_size <- 18
  geom_text_size <- text_size * (5/14)
  
  GOI_plot <- ggplot(data = data_for_compare) +
    geom_segment(aes(x=as.numeric(other_gene)-0.1,xend=as.numeric(other_gene)+0.1, y=ces_GOI, yend=ces_GOI_after_OTHER), 
                 arrow = arrow(length = unit(0.1, "inches"), type ="closed")) +
    geom_point(aes(x=as.numeric(other_gene)-0.1, y=ces_GOI), size=2) +
    geom_errorbar(aes(x=as.numeric(other_gene)-0.1, 
                      ymin= ci_low_95_ces_GOI, ymax= ci_high_95_ces_GOI),width = 0,alpha=0.2) +
    geom_errorbar(aes(x=as.numeric(other_gene)+0.1, 
                      ymin= ci_low_95_ces_GOI_after_OTHER, ymax= ci_high_95_ces_GOI_after_OTHER),width = 0,alpha=0.2) +
    scale_x_discrete(limits = levels(data_for_compare$other_gene)) +
    theme_classic() +
    # coord_cartesian(ylim=c(0,3500)) +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), lwd = 0.5, color = "lightgrey") +
    geom_hline(yintercept = 0, lwd = 0.5, color = "lightgrey", linetype = "dotted") +
    labs(x = paste0("Gene (paired with ", current_gene, ")"), y = paste0("\nScaled selection \ncoefficients of \n", current_gene, "\nwhen paired \ngene is \nwildtype (\u25cf) \nand mutated (\u25b8)")) +
    #scale_color_discrete(name="p-value", breaks=c("yes", "no"), labels=c("p < 0.05", "p ≥ 0.05")) +
    theme(axis.text = element_text(size = text_size), 
          axis.title = element_text(size = text_size), 
          axis.title.x = element_text(vjust = -5),
          axis.title.y = element_text(angle = 0, margin = margin(r = 30)),
          plot.margin = margin(0, 1, 2, 1, "cm")) +
    annotate("text", x = as.numeric(data_for_compare$other_gene), y = -(upper_limit1/2.8), label = data_for_compare$significance, size = geom_text_size) +
    coord_cartesian(ylim = c(0, upper_limit1), clip = "off")
  # guides(color = guide_legend(override.aes = aes(label = "|")))
  
  OTHER_plot <- ggplot(data = data_for_compare) +
    geom_segment(aes(x=as.numeric(other_gene)-0.1, xend=as.numeric(other_gene)+0.1, y=ces_OTHER, yend=ces_OTHER_after_GOI), 
                 arrow = arrow(length = unit(0.1, "inches"), type ="closed")) +
    geom_point(aes(x=as.numeric(other_gene)-0.1, y=ces_OTHER), size=2) +
    geom_errorbar(aes(x=as.numeric(other_gene)-0.1, 
                      ymin= ci_low_95_ces_OTHER, ymax= ci_high_95_ces_OTHER),width = 0,alpha=0.2) +
    geom_errorbar(aes(x=as.numeric(other_gene)+0.1, 
                      ymin= ci_low_95_ces_OTHER_after_GOI, 
                      ymax= ci_high_95_ces_OTHER_after_GOI),
                  width = 0,alpha=0.2) +
    scale_x_discrete(limits = levels(data_for_compare$other_gene)) +
    theme_classic() +
    theme(axis.text = element_text(size = text_size), 
          axis.title = element_text(size = text_size), 
          axis.title.x = element_text(vjust = -5),
          axis.title.y = element_text(angle=0, margin = margin(r = 30)),
          plot.margin = margin(1, 1, 2, 1, "cm")) +
    # coord_cartesian(ylim=c(0,3500)) +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), lwd = 0.5, color = "lightgrey") +
    geom_hline(yintercept = 0, lwd = 0.5, color = "lightgrey", linetype = "dotted") +
    labs(x = paste0("Gene (paired with ",current_gene, ")\n"), y = paste0("\nScaled selection \ncoefficients of \npaired gene when\n", current_gene, "\nis wildtype (\u25cf) \nand mutated (\u25b8)"))+
    scale_color_discrete(name="p-value", breaks=c("yes", "no"), labels=c("p < 0.05", "p ≥ 0.05")) +
    annotate("text", x = as.numeric(data_for_compare$other_gene), y = -(upper_limit2/2.8), label = data_for_compare$significance, size = geom_text_size) +
    coord_cartesian(ylim = c(0, upper_limit2), clip = "off")

  #guides(color = guide_legend(override.aes = aes(label = "|")))
  
  both_plots <- (OTHER_plot/GOI_plot)
  
  
  return(both_plots)
}

# Define genes for plots
# genes = c("NOTCH1", "TP53", "FAT1", "NOTCH2")
# 
# for (gene in genes) {
#   epistasis_with_ci = plot_epistasis_results_with_ci(gene_ep_results, gene)
#   ggsave(paste0("figures/epistasis_", gene, ".png"), epistasis_with_ci, width = 15, dpi = 300, height = 8)
# }

notch1_ep_plot <- plot_epistasis_results_with_ci(gene_ep_results, "NOTCH1", 2000, 4500)
tp53_ep_plot <- plot_epistasis_results_with_ci(gene_ep_results, "TP53", 9000, 3500)

ep_plots_combined <- plot_grid(notch1_ep_plot, tp53_ep_plot, labels = c("A", "B"), ncol = 1, label_size = 20, scale = 0.95)
ggsave("output_data/fig_3_combined_ep_plots_errorbars.png", ep_plots_combined, width=18, height=14)
ggsave("output_data/fig_3_combined_ep_plots_errorbars.pdf", ep_plots_combined, width=18, height=14, device=cairo_pdf)
ggsave("output_data/fig_3_combined_ep_plots_errorbars.jpg", ep_plots_combined, width=18, height=14)
